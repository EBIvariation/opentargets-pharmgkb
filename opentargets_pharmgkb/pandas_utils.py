import numpy as np
import pandas as pd


def read_tsv_to_df(path):
    return pd.read_csv(path, sep='\t', dtype=str, na_values=[''], keep_default_na=False)


def none_to_nan(x):
    return np.nan if x is None else x


def split_and_explode_column(df, source_col, target_col, sep=';', split_only=False):
    """
    Splits a string-valued column in dataframe and explodes on the values, storing them in the specified target column.
    Any white space around the separator will be stripped.

    :param df: Pandas dataframe
    :param source_col: name of column in df to split
    :param target_col: destination column name for exploded values
    :param sep: string separator to split source_col by (default ';')
    :param split_only: if True will only split on separator, leaving target_col as a list (default False)
    :return: dataframe with target_col added
    """
    split_cols = df.assign(**{target_col: df[source_col].str.split(pat=f'\s*{sep}\s*')})
    if not split_only:
        split_cols = split_cols.explode(target_col)
    return split_cols.reset_index(drop=True)
