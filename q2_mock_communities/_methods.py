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


def _parse_groups(groups: list | None, n_communities: int) -> list[list[int]] | None:
    """Parse List[Str] groups param (e.g. ['0,1', '2,3']) to list[list[int]]."""
    if not groups:
        return None
    try:
        parsed = [[int(c) for c in g.split(',')] for g in groups]
    except ValueError:
        raise ValueError(
            "--p-groups must be comma-separated integers per group, "
            "e.g. --p-groups '0,1' '2,3'"
        )
    all_indices = sorted(c for g in parsed for c in g)
    if all_indices != list(range(n_communities)):
        raise ValueError(
            f"--p-groups must partition community indices 0..{n_communities - 1} exactly once"
        )
    return parsed


def generate(
    n_communities: int = 3,
    taxa_per_community: list = None,
    overlap: float = 0.0,
    core_overlap: float = 0.0,
    overlap_mode: str = 'chain',
    groups: list = None,
    group_core_overlap: float = 0.0,
    alpha: float = 0.5,
    core_signature_strength: float = 1.0,
    core_weight: float = None,
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
    if group_core_overlap >= 1.0 and group_core_overlap == int(group_core_overlap):
        group_core_overlap = int(group_core_overlap)

    parsed_groups = _parse_groups(groups, n_communities)

    result = make_communities(
        n_communities=n_communities,
        taxa_per_community=taxa_per_community,
        overlap=overlap,
        core_overlap=core_overlap,
        overlap_mode=overlap_mode,
        groups=parsed_groups,
        group_core_overlap=group_core_overlap,
        alpha=alpha,
        core_signature_strength=core_signature_strength,
        core_weight=core_weight,
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
