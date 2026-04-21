# q2-mock-communities

A [QIIME 2](https://qiime2.org) plugin for generating structured mock microbial communities with controllable overlap, sparsity, and library depth. Useful for benchmarking microbiome analysis pipelines.

## Actions

| Action | Type | Description |
|---|---|---|
| `generate` | Method | Generate mock communities with configurable taxa overlap, Dirichlet evenness, and sequencing depth |
| `summarize` | Visualizer | Community-colored heatmap of sample proportions |

## Outputs

`generate` returns three artifacts:
- `FeatureTable[Frequency]` — integer count table
- `FeatureTable[RelativeFrequency]` — true proportions
- `CommunityMap` — TSV mapping each sample to its community label

## Installation

```bash
pip install -e . --no-deps
```

Install into your QIIME 2 conda environment. Verify with `qiime info`.

## Example

```bash
qiime mock-communities generate \
  --p-n-communities 3 \
  --p-taxa-per-community 20 \
  --p-overlap 0.2 \
  --p-n-samples 5 \
  --p-library-size 10000 \
  --p-seed 42 \
  --o-counts counts.qza \
  --o-proportions proportions.qza \
  --o-community-map community_map.qza
```
