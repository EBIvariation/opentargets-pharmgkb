import numpy as np
import pandas as pd


def read_tsv_to_df(path):
    return pd.read_csv(path, sep='\t', dtype=str)


def none_to_nan(x):
    return np.nan if x is None else x


def explode_column(df, source_col, target_col, sep=';'):
    """
    Splits a string-valued column in dataframe and explodes on the values, storing them in the specified target column.
    Any white space around the separator will be stripped.

    :param df: Pandas dataframe
    :param source_col: name of column in df to split
    :param target_col: destination column name for exploded values
    :param sep: string separator to split source_col by (default ';')
    :return: dataframe with target_col added
    """
    split_cols = df.assign(**{target_col: df[source_col].str.split(sep)}).explode(target_col).reset_index(drop=True)
    split_cols[target_col] = split_cols[target_col].map(lambda x: str(x).strip() if pd.notna(x) else np.nan)
    return split_cols
