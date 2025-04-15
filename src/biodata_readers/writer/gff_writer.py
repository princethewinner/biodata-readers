import typing as tp
import os
from pathlib import Path
from biodata_readers.writer.abstraction import Writer
import pandas as pd
import json
from loguru import logger
from gffpandas.gffpandas import Gff3DataFrame


class GFFWriter(Writer):  # type: ignore

    # seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    # source - name of the program that generated this feature, or the data source (database or project name)
    # feature - feature type name, e.g. Gene, Variation, Similarity
    # start - Start position* of the feature, with sequence numbering starting at 1.
    # end - End position* of the feature, with sequence numbering starting at 1.
    # score - A floating point value.
    # strand - defined as + (forward) or - (reverse).
    # frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    # attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

    SEQUENCE_COL: tp.ClassVar[str] = "seq_id"
    SOURCE_COL: tp.ClassVar[str] = "source"
    FEATURE_COL: tp.ClassVar[str] = "type"
    START_COL: tp.ClassVar[str] = "start"
    END_COL: tp.ClassVar[str] = "end"
    SCORE_COL: tp.ClassVar[str] = "score"
    STRAND_COL: tp.ClassVar[str] = "strand"
    FRAME_COL: tp.ClassVar[str] = "phase"
    ATTR_COL: tp.ClassVar[str] = "attributes"

    REQUIRED_COLUMN_ORDER: tp.List[str] = [
        SEQUENCE_COL,
        SOURCE_COL,
        FEATURE_COL,
        START_COL,
        END_COL,
        SCORE_COL,
        STRAND_COL,
        FRAME_COL,
        ATTR_COL,
    ]

    CANONICAL_CHROMOSOMES: tp.List[str] = list(
        map(lambda x: f"chr{x}", list(range(1, 23)) + ["X", "Y"])
    )

    def __init__(
        self,
        data: pd.DataFrame,
        output_file: tp.Union[str, Path],
        only_required: bool = True,  # Unused in this, keeping this here for consistency across APIs.
        only_canonical: bool = True,
    ) -> None:
        super().__init__(output_file=output_file)
        self._data: pd.DataFrame = data
        self.only_required: bool = only_required
        self.only_canonical: bool = only_canonical

    def _write(self, data: pd.DataFrame, output_file: str) -> None:

        _gff3_data_frame = Gff3DataFrame(input_df=data, input_header="")
        _gff3_data_frame.to_gff3(gff_file=output_file)
        logger.success(f"GFF file written at location: {output_file}")

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


class GFFConverter(GFFWriter):

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
            and len(set(self.REQUIRED_COLUMN_ORDER) - set(self.column_map.values()))
            != 0
        ):
            raise ValueError(
                f"Please provide the minimal columns {self.REQUIRED_COLUMN_ORDER}"
            )

        self._data = self._data.rename(self.column_map, axis=1)

    def convert(
        self, is_one_indexed: bool = True, is_end_included: bool = True
    ) -> None:

        if not is_one_indexed:
            self._data[self.START_COL] = self._data[self.START_COL] + 1
            if is_end_included:
                self._data[self.END_COL] = self._data[self.END_COL] + 1

        self._data = self._data[self.REQUIRED_COLUMN_ORDER]
