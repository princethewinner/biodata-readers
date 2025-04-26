from collections import defaultdict
import pandas as pd
import typing as tp
from loguru import logger
from tqdm.auto import tqdm

_T_Str_Generic_Dict: tp.TypeAlias = tp.Dict[str, tp.Any]
_T_Str_Str_Dict: tp.TypeAlias = tp.Dict[str, str]
_T_Str_Int_Dict: tp.TypeAlias = tp.Dict[str, int]
_T_Set_Str: tp.TypeAlias = tp.Set[str]
_T_List_Str: tp.TypeAlias = tp.List[str]


ATTR_DELIM: str = " "
TAG_ATTRIBUTE_COLUMN: str = "tag"
SPECIAL_HANDLING_BINARY: str = "binary"
SPECIAL_HANDLING_DURING_EXPANSION: _T_List_Str = [TAG_ATTRIBUTE_COLUMN]
SPECIAL_HANDLING_MAP: _T_Str_Str_Dict = {TAG_ATTRIBUTE_COLUMN: SPECIAL_HANDLING_BINARY}
SPECIAL_HANDLING_SEP: str = "__"


class AttributeHolder:

    def __init__(self) -> None:
        self.attributes: _T_Str_Generic_Dict = defaultdict(list)
        self.attributes_length: _T_Str_Int_Dict = defaultdict(int)

    def _insertRow(self, dict_values: _T_Str_Generic_Dict) -> None:

        for key, value in dict_values.items():
            self.attributes[key].append(value)
            self.attributes_length[key] += 1

        missing_keys: _T_Set_Str = self.available_attributes - set(dict_values.keys())
        for key in missing_keys:
            self.attributes[key].append(None)
            self.attributes_length[key] += 1

    def _reconcileNewKeys(self, new_keys: _T_Set_Str) -> None:
        """
        Assumption: All list are of same length
        """

        random_key: str = next(iter(self.attributes.keys()))
        list_length: int = self.attributes_length[random_key]
        for key in new_keys:
            self.attributes[key] = [None] * list_length
            self.attributes_length[key] = list_length

    def insertRow(self, dict_values: _T_Str_Generic_Dict) -> None:

        if len(self.available_attributes) > 0:
            new_keys: _T_Set_Str = set(dict_values.keys()) - self.available_attributes
            if len(new_keys) > 0:
                self._reconcileNewKeys(new_keys)

        self._insertRow(dict_values)

    @property
    def as_df(self) -> pd.DataFrame:
        return pd.DataFrame.from_dict(dict(self.attributes))

    @property
    def available_attributes(self) -> _T_Set_Str:
        return set(self.attributes.keys())


def fixSpecialAttributes(attribute_df: pd.DataFrame) -> pd.DataFrame:

    attr_name: str
    attr_split_holder: _T_List_Str
    for columns in attribute_df.columns:
        attr_split_holder = columns.split(SPECIAL_HANDLING_SEP)
        if len(attr_split_holder) > 1:
            attr_name = attr_split_holder[0]
            if attr_name in SPECIAL_HANDLING_DURING_EXPANSION:
                if SPECIAL_HANDLING_MAP[attr_name] == SPECIAL_HANDLING_BINARY:
                    attribute_df[columns] = attribute_df[columns].fillna(0).astype(int)
    return attribute_df


def expandAttributes(gtf_data: pd.DataFrame) -> pd.DataFrame:  # noqa: C901

    attribute_holder: AttributeHolder = AttributeHolder()
    list_attribute: _T_List_Str
    attr_name: str
    attr_val: tp.Any
    split_holder: _T_List_Str
    handling_type: str

    attribute_dict: defaultdict[str, str]
    for attribute in tqdm(gtf_data.attribute, total=gtf_data.shape[0]):
        attribute_dict = defaultdict(str)
        list_attribute = attribute.split(";")
        for _la in list_attribute:
            if len(_la) > 0:
                # attr_name, attr_val = _la.strip(" ").split(self.ATTR_DELIM)
                split_holder = _la.strip(" ").split(ATTR_DELIM)
                attr_name = split_holder[0]
                attr_val = ATTR_DELIM.join(split_holder[1:])
                attr_name = attr_name.strip('"')
                attr_val = attr_val.strip('"')
                if attr_name in SPECIAL_HANDLING_DURING_EXPANSION:
                    handling_type = SPECIAL_HANDLING_MAP[attr_name]
                    if handling_type == SPECIAL_HANDLING_BINARY:
                        attr_name = SPECIAL_HANDLING_SEP.join(
                            [attr_name, str(attr_val)]
                        )
                        attr_val = 1
                attribute_dict[attr_name] = attr_val
        attribute_holder.insertRow(attribute_dict)

    logger.debug("Creating final dataframe.")
    attribute_df: pd.DataFrame = attribute_holder.as_df
    logger.debug("Fixing special attributes")
    attribute_df = fixSpecialAttributes(attribute_df)
    gtf_data = pd.concat([gtf_data, attribute_df], axis=1)
    gtf_data.drop("attribute", inplace=True, axis=1)
    return gtf_data
