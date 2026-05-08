import typing as tp
from pathlib import Path
import pandas as pd
from abc import ABC, abstractmethod


class Writer(ABC):

    def __init__(self, output_file: tp.Union[str, Path]) -> None:
        self.output_file: tp.Union[str, Path] = output_file

    @abstractmethod
    def write(self, chromosome_wise: bool) -> None:
        pass

    @property
    @abstractmethod
    def data(self) -> pd.DataFrame:
        raise NotImplementedError

    @data.setter
    @abstractmethod
    def data(self, val: pd.DataFrame) -> None:
        raise NotADirectoryError
