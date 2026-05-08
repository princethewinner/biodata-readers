import typing as tp
from pathlib import Path
from biodata_readers.writer.abstraction import Writer


class VCFWriter(Writer):  # type: ignore

    def __init__(self, output_file: tp.Union[str, Path]) -> None:
        super().__init__(output_file=output_file)

    def write(self) -> None:
        pass


class VCFConverter(VCFWriter):

    def __init__(
        self, column_map_file: tp.Union[str, Path], output_file: tp.Union[str, Path]
    ) -> None:
        super().__init__(output_file=output_file)
