from qiime2.plugin import (
    Plugin, Int, Float, Bool, Str, List, Range, Choices,
)
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency

from . import __version__
from ._types import CommunityMap
from ._formats import CommunityMapFormat, CommunityMapDirFmt
from ._methods import generate
from ._visualizers import summarize

plugin = Plugin(
    name='mock-communities',
    version=__version__,
    website='https://github.com/example/q2-mock-communities',
    package='q2_mock_communities',
    description=(
        'Generate structured mock microbial communities with controllable '
        'overlap, sparsity, and library depth for benchmarking purposes.'
    ),
    short_description='Mock microbial community simulator.',
)

plugin.register_semantic_types(CommunityMap)
plugin.register_formats(CommunityMapFormat, CommunityMapDirFmt)
plugin.register_artifact_class(
    CommunityMap,
    directory_format=CommunityMapDirFmt,
)

plugin.methods.register_function(
    function=generate,
    inputs={},
    parameters={
        'n_communities': Int % Range(1, None),
        'taxa_per_community': List[Int % Range(1, None)],
        'overlap': Float % Range(0, 1, inclusive_end=True),
        'core_overlap': Float % Range(0, None),
        'overlap_mode': Str % Choices('chain', 'core', 'both'),
        'alpha': Float % Range(0, None),
        'n_samples': List[Int % Range(1, None)],
        'library_size': Int % Range(1, None),
        'shuffle_taxa': Bool,
        'seed': Int,
    },
    outputs=[
        ('counts', FeatureTable[Frequency]),
        ('proportions', FeatureTable[RelativeFrequency]),
        ('community_map', CommunityMap),
    ],
    input_descriptions={},
    parameter_descriptions={
        'n_communities': 'Number of distinct communities to generate.',
        'taxa_per_community': (
            'Number of taxa per community. A single value applies uniformly; '
            'a list of n_communities values sets each community independently '
            '(e.g. --p-taxa-per-community 5 20 8).'
        ),
        'overlap': (
            'Fraction (0–1) of taxa shared between adjacent communities (chain overlap). '
            '0 = fully disjoint, 1 = fully overlapping.'
        ),
        'core_overlap': (
            'Taxa shared by ALL communities. Values in (0, 1) are treated as a fraction '
            'of the smallest community; whole numbers >= 1 (e.g. 3.0) are treated as an '
            'exact count. Used when overlap_mode is core or both.'
        ),
        'overlap_mode': (
            'chain = adjacent communities share taxa; '
            'core = all communities share a common block; '
            'both = core shared by all plus chain overlap among unique portions.'
        ),
        'alpha': (
            'Dirichlet concentration parameter controlling within-community '
            'evenness. <1 → sparse/uneven, =1 → uniform, >1 → even.'
        ),
        'n_samples': (
            'Samples per community. A single value applies uniformly; '
            'a list of n_communities values sets each community independently '
            '(e.g. --p-n-samples 2 5 3).'
        ),
        'library_size': 'Total sequencing depth per sample (multinomial draws).',
        'shuffle_taxa': 'Randomly permute taxon column order.',
        'seed': 'Random seed for reproducibility.',
    },
    output_descriptions={
        'counts': 'Integer count table (samples × taxa) in BIOM format.',
        'proportions': 'True relative abundance table (samples × taxa).',
        'community_map': (
            'TSV mapping each sample-id to its community label. '
            'Pass to the summarize visualizer for community-colored plots.'
        ),
    },
    name='Generate mock communities',
    description=(
        'Generate structured mock microbial communities with controllable '
        'overlap, sparsity, sample count, and library depth. Returns integer '
        'counts, true proportions, and a sample-to-community mapping.'
    ),
)

plugin.visualizers.register_function(
    function=summarize,
    inputs={
        'proportions': FeatureTable[RelativeFrequency],
        'community_map': CommunityMap,
    },
    parameters={
        'log_scale': Bool,
        'show_boundaries': Bool,
    },
    input_descriptions={
        'proportions': 'Proportions table from the generate action.',
        'community_map': 'Community mapping from the generate action.',
    },
    parameter_descriptions={
        'log_scale': 'If true, plot log₁₀(proportion + 1e-6).',
        'show_boundaries': (
            'Draw dashed lines between community sample groups '
            '(only meaningful when shuffle_taxa=False).'
        ),
    },
    name='Summarize mock communities',
    description=(
        'Generate an interactive heatmap of community proportions with '
        'samples colored by community assignment.'
    ),
)

# Register transformers — imported last to avoid circular imports
from . import _transformers  # noqa: F401, E402
