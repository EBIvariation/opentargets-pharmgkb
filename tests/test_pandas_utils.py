import pandas as pd

from opentargets_pharmgkb.pandas_utils import explode_column


def test_explode_column():
    df = pd.DataFrame([
        [1, 'apple; pear; banana'],
        [2, 'cat;frog'],
        [3, 'something']
    ], columns=['A', 'B'])

    expected = pd.DataFrame([
        [1, 'apple; pear; banana', 'apple'],
        [1, 'apple; pear; banana', 'pear'],
        [1, 'apple; pear; banana', 'banana'],
        [2, 'cat;frog', 'cat'],
        [2, 'cat;frog', 'frog'],
        [3, 'something', 'something']
    ], columns=['A', 'B', 'C'])
    result = explode_column(df, 'B', 'C')

    assert result.equals(expected)
