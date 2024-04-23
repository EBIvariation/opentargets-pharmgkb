import numpy as np
import pandas as pd

from opentargets_pharmgkb.pandas_utils import explode_column


def test_explode_column():
    df = pd.DataFrame([
        [1, 'apple; pear; banana'],
        [2, 'cat;frog'],
        [3, 'something'],
        [4, np.nan]
    ], columns=['A', 'B'])

    expected = pd.DataFrame([
        [1, 'apple; pear; banana', 'apple'],
        [1, 'apple; pear; banana', 'pear'],
        [1, 'apple; pear; banana', 'banana'],
        [2, 'cat;frog', 'cat'],
        [2, 'cat;frog', 'frog'],
        [3, 'something', 'something'],
        [4, np.nan, np.nan]
    ], columns=['A', 'B', 'C'])
    result = explode_column(df, 'B', 'C')

    assert result.equals(expected)


def test_explode_column_split_only():
    df = pd.DataFrame([
        [1, 'apple1;apple2 / pear / banana'],
        [2, 'cat/frog'],
        [3, 'something'],
        [4, np.nan]
    ], columns=['A', 'B'])

    expected = pd.DataFrame([
        [1, 'apple1;apple2 / pear / banana', ['apple1;apple2', 'pear', 'banana']],
        [2, 'cat/frog', ['cat', 'frog']],
        [3, 'something', ['something']],
        [4, np.nan, np.nan]
    ], columns=['A', 'B', 'C'])
    result = explode_column(df, 'B', 'C', sep='/', split_only=True)

    assert result.equals(expected)
