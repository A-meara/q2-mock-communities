import pandas as pd

from .plugin_setup import plugin
from ._formats import CommunityMapFormat


@plugin.register_transformer
def _1(data: pd.DataFrame) -> CommunityMapFormat:
    ff = CommunityMapFormat()
    data.to_csv(str(ff), sep='\t', index=True, index_label='sample-id')
    return ff


@plugin.register_transformer
def _2(ff: CommunityMapFormat) -> pd.DataFrame:
    df = pd.read_csv(str(ff), sep='\t', index_col=0)
    df.index.name = 'sample-id'
    return df
