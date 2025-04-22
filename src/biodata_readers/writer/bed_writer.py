import typing as tp
import os
from pathlib import Path
from fuc.api.pybed import BedFrame
from biodata_readers.writer.abstraction import Writer
import pandas as pd
import json
from loguru import logger


class BedWriter(Writer):  # type: ignore

    SEQUENCE_COL: tp.ClassVar[str] = "Chromosome"
    START_COL: tp.ClassVar[str] = "Start"
    END_COL: tp.ClassVar[str] = "End"
    NAME_COL: tp.ClassVar[str] = "Type"

    MINIMAL_REQUIRED_COLUMNS: tp.List[str] = [
        SEQUENCE_COL,
        START_COL,
        END_COL,
        NAME_COL,
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
    ) -> None:
        super().__init__(output_file=output_file)
        self._data: pd.DataFrame = data
        self.only_required: bool = only_required
        self.only_canonical: bool = only_canonical

    def _write(self, data: pd.DataFrame, output_file: str) -> None:

        meta: tp.List[str] = []
        if self.only_required:
            data = data[self.MINIMAL_REQUIRED_COLUMNS]

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
    ) -> None:
        super().__init__(
            data=data,
            output_file=output_file,
            only_required=only_required,
            only_canonical=only_canonical,
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

        column_order: tp.List[str] = list(self._data.columns)
        for rq in self.MINIMAL_REQUIRED_COLUMNS:
            column_order.remove(rq)

        column_order = self.MINIMAL_REQUIRED_COLUMNS + column_order
        self._data = self._data[column_order]
