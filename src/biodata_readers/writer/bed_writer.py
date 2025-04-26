import typing as tp
import os
from pathlib import Path
from fuc.api.pybed import BedFrame
from biodata_readers.writer.abstraction import Writer
import pandas as pd
import pandas.api.typing as pdtp
import json
from loguru import logger
from tqdm.auto import tqdm


class BedWriter(Writer):  # type: ignore

    SEQUENCE_COL: tp.ClassVar[str] = "Chromosome"
    START_COL: tp.ClassVar[str] = "Start"
    END_COL: tp.ClassVar[str] = "End"
    NAME_COL: tp.ClassVar[str] = "Name"
    SCORE_COL: tp.ClassVar[str] = "Score"
    STRAND_COL: tp.ClassVar[str] = "Strand"

    MINIMAL_REQUIRED_COLUMNS: tp.List[str] = [
        SEQUENCE_COL,
        START_COL,
        END_COL,
        NAME_COL
    ]

    ADDITIONAL_COLUMN_ORDER: tp.List[str] = [
        SCORE_COL,
        STRAND_COL
    ]

    UNIQUE_GROUP_COLUMNS: tp.List[str] = [
        SEQUENCE_COL,
        START_COL,
        END_COL
    ]

    CANONICAL_CHROMOSOMES: tp.List[str] = list(
        map(lambda x: f"chr{x}", list(range(1, 23)) + ["X", "Y"])
    )

    def __init__(
        self,
        data: pd.DataFrame,
        output_file: tp.Union[str, Path],
        only_required: bool = True,
        only_canonical: bool = True,
        only_unique: bool = True,
        additional_column: tp.Optional[tp.List[str]] = None,
    ) -> None:
        super().__init__(output_file=output_file)
        self._data: pd.DataFrame = data
        self.only_required: bool = only_required
        self.only_canonical: bool = only_canonical
        self.only_unique: bool = only_unique
        self.additional_column: tp.Optional[tp.List[str]] = additional_column

    def _getColumnOrder(self, data: tp.Optional[pd.DataFrame] = None) -> tp.List[str]:
        
        if data is None:
            data = self._data
        
        column_order: tp.List[str] = list(data.columns)
        for rq in self.MINIMAL_REQUIRED_COLUMNS:
            column_order.remove(rq)

        additional_columns: tp.List[str] = []
        for c in self.ADDITIONAL_COLUMN_ORDER:
            if c in column_order:
                additional_columns.append(c)

        column_order = self.MINIMAL_REQUIRED_COLUMNS + additional_columns
        return column_order

    def _getOnlyUnique(self, data: pd.DataFrame) -> pd.DataFrame:

        group_by_columns: tp.List[str] = self.UNIQUE_GROUP_COLUMNS.copy()

        if self.STRAND_COL in data.columns:
            logger.debug(f"{self.STRAND_COL} found in data columns. Adding it to group by columns list.")
            group_by_columns.append(self.STRAND_COL)

        logger.info("Merging non unique columns.")
        group: pdtp.DataFrameGroupBy = data.groupby(by=group_by_columns)
        merge_keys: tp.List[str] = list(set(data.columns) - set(group_by_columns))
        sampled_df: pd.DataFrame
        key_holder: tp.List[str]
        data_holder: tp.List[tp.List[tp.Any]] = []

        logger.debug(f"Keys selected for merging: {merge_keys}")

        gr: tp.Any
        for gr in tqdm(group.groups.keys(), total=len(group), desc="Merging non unique"):
            sampled_df = group.get_group(gr)[merge_keys]
            key_holder = []
            for k in merge_keys:
                key_holder.append(
                    ",".join(set(sampled_df[k]))
                )
            gr = list(gr)
            gr.extend(key_holder)
            data_holder.append(gr)

        column_header: tp.List[str] = group_by_columns + merge_keys
        return pd.DataFrame(data_holder, columns=column_header)

    def _write(self, data: pd.DataFrame, output_file: str) -> None:

        meta: tp.List[str] = []
        req_columns: tp.List[str] = self.MINIMAL_REQUIRED_COLUMNS.copy()
        if not self.only_required:
            if self.additional_column is not None:
                req_columns.extend(self.additional_column)

        data = data[req_columns]

        logger.debug(f"Columns in selected data frame: {data.columns}")
        logger.debug(f"Total number of samples in dataframe: {data.shape}")

        column_order: tp.List[str]
        if self.only_unique:
            data = self._getOnlyUnique(data)
            column_order = self._getColumnOrder(data)
            data = data[column_order]

        logger.debug(f"Total number of samples in dataframe (AFTER selecting uniques): {data.shape}")
        _bedDataFrame: BedFrame = BedFrame.from_frame(meta=meta, data=data)
        _bedDataFrame.to_file(output_file)
        logger.success(f"BED file written at location: {output_file}")

    def write(self, chromosome_wise: bool) -> None:

        _ldata: pd.DataFrame = self._data
        if self.only_canonical:
            _ldata = _ldata[_ldata[self.SEQUENCE_COL].isin(self.CANONICAL_CHROMOSOMES)]

        if chromosome_wise:
            _sample_location: str = os.path.splitext(self.output_file)[0]
            os.makedirs(_sample_location, exist_ok=True)
            unique_chromosomes: tp.List[str] = list(_ldata[self.SEQUENCE_COL].unique())
            for _uq in unique_chromosomes:
                _sample_df: pd.DataFrame = _ldata[_ldata[self.SEQUENCE_COL] == _uq]
                _sample_path: str = os.path.join(_sample_location, f"{_uq}.bed")
                self._write(_sample_df, _sample_path)
        else:
            self._write(_ldata, self.output_file)

    @property
    def data(self) -> pd.DataFrame:
        return self._data

    @data.setter
    def data(self, val: pd.DataFrame) -> None:
        self._data = val


class BedConverter(BedWriter):

    # from: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
    # chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671). Many assemblies also support several different chromosome aliases (e.g. '1' or 'NC_000001.11' in place of 'chr1').
    # chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
    # chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature, however, the number in position format will be represented. For example, the first 100 bases of chromosome 1 are defined as chrom=1, chromStart=0, chromEnd=100, and span the bases numbered 0-99 in our software (not 0-100), but will represent the position notation chr1:1-100. Read more here.
    # chromStart and chromEnd can be identical, creating a feature of length 0, commonly used for insertions. For example, use chromStart=0, chromEnd=0 to represent an insertion before the first nucleotide of a chromosome.

    def __init__(
        self,
        data: pd.DataFrame,
        column_map_file: tp.Union[str, Path, tp.Dict[str, str]],
        output_file: tp.Union[str, Path],
        only_required: bool = True,
        only_canonical: bool = True,
        only_unique: bool = True,
        additional_column: tp.Optional[tp.List[str]] = None,
    ) -> None:
        super().__init__(
            data=data,
            output_file=output_file,
            only_required=only_required,
            only_canonical=only_canonical,
            additional_column=additional_column,
            only_unique=only_unique
        )
        self.column_map: tp.Union[str, Path, tp.Dict[str, str]] = column_map_file

        if isinstance(self.column_map, str):
            with open(self.column_map, "r") as fid:
                self.column_map = json.load(fid)

        if (
            isinstance(self.column_map, dict)
            and len(set(self.MINIMAL_REQUIRED_COLUMNS) - set(self.column_map.values()))
            != 0
        ):
            raise ValueError(
                f"Please provide the minimal columns {self.MINIMAL_REQUIRED_COLUMNS}"
            )

        self._data = self._data.rename(self.column_map, axis=1)

    def convert(
        self, is_one_indexed: bool = True, is_end_included: bool = True
    ) -> None:

        if is_one_indexed:
            self._data[self.START_COL] = self._data[self.START_COL] - 1
            if not is_end_included:
                self._data[self.END_COL] = self._data[self.END_COL] - 1
        else:
            if is_end_included:
                self._data[self.END_COL] = self._data[self.END_COL] + 1

        column_order: tp.List[str] = self._getColumnOrder()
        self._data = self._data[column_order]
