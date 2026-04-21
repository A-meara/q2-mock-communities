import biom
import numpy as np
import pandas as pd

from ._core import make_communities


def _resolve_list_param(value: list, n_communities: int, name: str):
    """Collapse a 1-element list to scalar, validate length otherwise."""
    if len(value) == 1:
        return value[0]
    if len(value) != n_communities:
        raise ValueError(
            f"--p-{name} must have length 1 (uniform) or {n_communities} "
            f"(one per community), got {len(value)}"
        )
    return value


def generate(
    n_communities: int = 3,
    taxa_per_community: list = None,
    overlap: float = 0.0,
    core_overlap: float = 0.0,
    overlap_mode: str = 'chain',
    alpha: float = 0.5,
    n_samples: list = None,
    library_size: int = 10_000,
    shuffle_taxa: bool = False,
    seed: int = None,
) -> (biom.Table, biom.Table, pd.DataFrame):
    """Generate mock microbial communities."""
    if taxa_per_community is None:
        taxa_per_community = [10]
    if n_samples is None:
        n_samples = [1]

    taxa_per_community = _resolve_list_param(taxa_per_community, n_communities, 'taxa-per-community')
    n_samples = _resolve_list_param(n_samples, n_communities, 'n-samples')

    # treat whole-number floats >= 1 as exact taxon counts (e.g. 3.0 → 3)
    if core_overlap >= 1.0 and core_overlap == int(core_overlap):
        core_overlap = int(core_overlap)

    result = make_communities(
        n_communities=n_communities,
        taxa_per_community=taxa_per_community,
        overlap=overlap,
        core_overlap=core_overlap,
        overlap_mode=overlap_mode,
        alpha=alpha,
        n_samples=n_samples,
        library_size=library_size,
        shuffle_taxa=shuffle_taxa,
        seed=seed,
    )

    def df_to_biom(df, dtype):
        return biom.Table(
            df.T.values.astype(dtype),
            observation_ids=list(df.columns),
            sample_ids=list(df.index),
        )

    counts_biom = df_to_biom(result['counts'], int)
    props_biom = df_to_biom(result['proportions'], float)

    community_df = result['metadata'].copy()
    community_df['community'] = community_df['community'].astype(str)
    community_df.index.name = 'sample-id'

    return counts_biom, props_biom, community_df
