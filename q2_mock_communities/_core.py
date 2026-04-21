"""
Core community simulator — copied from mock_communities.py (Part 1).
Edit mock_communities.py (for notebooks) and this file independently.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors  # noqa: F401
from matplotlib.patches import Patch


def make_communities(
    n_communities: int = 3,
    taxa_per_community: int | list[int] = 10,
    overlap: float = 0.0,
    core_overlap: float | int = 0.0,
    overlap_mode: str = 'chain',
    alpha: float = 0.5,
    n_samples: int | list[int] = 1,
    library_size: int = 10_000,
    shuffle_taxa: bool = False,
    seed: int | None = None,
) -> dict:
    """Generate mock microbial communities with controlled structure."""
    if overlap_mode not in ('chain', 'core', 'both'):
        raise ValueError(f"overlap_mode must be 'chain', 'core', or 'both', got {overlap_mode!r}")

    rng = np.random.default_rng(seed)

    if isinstance(n_samples, int):
        samples_per = [n_samples] * n_communities
    else:
        if len(n_samples) != n_communities:
            raise ValueError(
                f"n_samples list length ({len(n_samples)}) != "
                f"n_communities ({n_communities})"
            )
        samples_per = list(n_samples)

    if isinstance(taxa_per_community, int):
        taxa_per = [taxa_per_community] * n_communities
    else:
        if len(taxa_per_community) != n_communities:
            raise ValueError(
                f"taxa_per_community list length ({len(taxa_per_community)}) != "
                f"n_communities ({n_communities})"
            )
        taxa_per = list(taxa_per_community)

    # --- taxa layout ---
    if overlap_mode not in ('core', 'both'):
        n_core = 0
    elif isinstance(core_overlap, int):
        n_core = core_overlap
    else:
        n_core = round(core_overlap * min(taxa_per))

    if overlap_mode == 'chain':
        overlap_counts = [
            round(overlap * min(taxa_per[i], taxa_per[i + 1]))
            for i in range(n_communities - 1)
        ]
        for i, ov in enumerate(overlap_counts):
            if taxa_per[i] - ov < 1:
                raise ValueError(
                    f"overlap={overlap} between communities {i} and {i+1} leaves stride<1"
                )
        starts = [0]
        for i in range(n_communities - 1):
            starts.append(starts[-1] + taxa_per[i] - overlap_counts[i])
        total_taxa = starts[-1] + taxa_per[-1]
        all_taxa = [f"taxon_{i + 1:03d}" for i in range(total_taxa)]
        taxon_membership: list[list[int]] = [[] for _ in range(total_taxa)]
        community_indices = []
        for c in range(n_communities):
            idx = list(range(starts[c], starts[c] + taxa_per[c]))
            community_indices.append(np.array(idx))
            for t in idx:
                taxon_membership[t].append(c)
        boundary_lines = starts[1:]

    else:  # 'core' or 'both'
        for c, tp in enumerate(taxa_per):
            if n_core > tp:
                raise ValueError(
                    f"core_overlap={core_overlap} gives n_core={n_core} which exceeds "
                    f"community {c} (taxa_per={tp})"
                )
        unique_per = [tp - n_core for tp in taxa_per]

        if overlap_mode == 'core':
            overlap_counts = []
            unique_starts = [n_core]
            for c in range(n_communities - 1):
                unique_starts.append(unique_starts[-1] + unique_per[c])
        else:  # 'both'
            overlap_counts = [
                round(overlap * min(unique_per[i], unique_per[i + 1]))
                for i in range(n_communities - 1)
            ]
            for i, ov in enumerate(overlap_counts):
                if unique_per[i] - ov < 1:
                    raise ValueError(
                        f"overlap={overlap} on unique portion of communities {i} and {i+1} "
                        f"leaves stride<1; reduce overlap or increase taxa_per_community"
                    )
            unique_starts = [n_core]
            for i in range(n_communities - 1):
                unique_starts.append(unique_starts[-1] + unique_per[i] - overlap_counts[i])

        total_taxa = unique_starts[-1] + unique_per[-1]
        all_taxa = [f"taxon_{i + 1:03d}" for i in range(total_taxa)]
        taxon_membership = [[] for _ in range(total_taxa)]

        core_idx = list(range(n_core))
        for t in core_idx:
            taxon_membership[t] = list(range(n_communities))

        community_indices = []
        for c in range(n_communities):
            unique_idx = list(range(unique_starts[c], unique_starts[c] + unique_per[c]))
            community_indices.append(np.array(core_idx + unique_idx))
            for t in unique_idx:
                if c not in taxon_membership[t]:
                    taxon_membership[t].append(c)

        boundary_lines = unique_starts

    # --- generate mean proportions per community ---
    community_means = np.zeros((n_communities, total_taxa))
    for c in range(n_communities):
        idx = community_indices[c]
        community_means[c, idx] = rng.dirichlet(np.full(len(idx), alpha))

    all_proportions = []
    all_counts = []
    meta_rows = []
    sample_idx = 0

    for c in range(n_communities):
        mean_p = community_means[c]
        for _ in range(samples_per[c]):
            concentration = alpha * taxa_per[c] * mean_p
            nonzero = mean_p > 0
            p = np.zeros(total_taxa)
            p[nonzero] = rng.dirichlet(concentration[nonzero])
            counts = rng.multinomial(library_size, p)
            sid = f"sample_{sample_idx:03d}"
            all_proportions.append(p)
            all_counts.append(counts)
            meta_rows.append({"sample_id": sid, "community": c})
            sample_idx += 1

    proportions_arr = np.array(all_proportions)
    counts_arr = np.array(all_counts)
    sample_ids = [r["sample_id"] for r in meta_rows]

    col_order = np.arange(total_taxa)
    if shuffle_taxa:
        col_order = rng.permutation(total_taxa)
    ordered_taxa = [all_taxa[i] for i in col_order]

    counts_df = pd.DataFrame(counts_arr[:, col_order], index=sample_ids, columns=ordered_taxa)
    proportions_df = pd.DataFrame(proportions_arr[:, col_order], index=sample_ids, columns=ordered_taxa)
    metadata = pd.DataFrame(meta_rows).set_index("sample_id")

    taxa_table = pd.DataFrame({
        "taxon": [all_taxa[i] for i in col_order],
        "communities": [taxon_membership[i] for i in col_order],
    })

    return {
        "counts": counts_df,
        "proportions": proportions_df,
        "metadata": metadata,
        "taxa_table": taxa_table,
        "community_means": pd.DataFrame(
            community_means[:, col_order],
            index=[f"community_{c}" for c in range(n_communities)],
            columns=ordered_taxa,
        ),
        "params": {
            "n_communities": n_communities,
            "taxa_per_community": taxa_per,
            "overlap": overlap,
            "core_overlap": core_overlap,
            "overlap_mode": overlap_mode,
            "n_core": n_core,
            "overlap_counts": overlap_counts,
            "boundary_lines": boundary_lines,
            "alpha": alpha,
            "n_samples": samples_per,
            "library_size": library_size,
            "shuffle_taxa": shuffle_taxa,
            "seed": seed,
            "total_taxa": total_taxa,
        },
    }


def plot_community_heatmap(
    result: dict,
    log_scale: bool = True,
    show_boundaries: bool = True,
    figsize: tuple | None = None,
    ax=None,
):
    """Heatmap of community proportions with community annotations."""
    props = result["proportions"]
    metadata = result["metadata"]
    params = result["params"]
    n_samples_total, n_taxa = props.shape

    if figsize is None:
        w = min(max(6, n_taxa * 0.25), 20)
        h = min(max(4, n_samples_total * 0.35), 20)
        figsize = (w, h)

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure

    data = props.values.copy()
    if log_scale:
        data = np.log10(data + 1e-6)
        label = "log₁₀(proportion + 1e-6)"
    else:
        label = "proportion"

    im = ax.imshow(data, aspect="auto", cmap="viridis", interpolation="nearest")
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label(label)

    comm_colors_list = plt.cm.Set2.colors
    communities = metadata["community"].values
    if n_samples_total <= 50:
        ax.set_yticks(range(n_samples_total))
        ax.set_yticklabels(props.index, fontsize=7)
        for i, c in enumerate(communities):
            ax.get_yticklabels()[i].set_color(comm_colors_list[c % len(comm_colors_list)])
    else:
        cum = 0
        positions, labels = [], []
        for c in range(params["n_communities"]):
            n = params["n_samples"][c]
            positions.append(cum + (n - 1) / 2.0)
            labels.append(f"Community {c}  (n={n})")
            cum += n
        ax.set_yticks(positions)
        ax.set_yticklabels(labels, fontsize=8)
        for i, c in enumerate(range(params["n_communities"])):
            ax.get_yticklabels()[i].set_color(comm_colors_list[c % len(comm_colors_list)])

    if n_taxa <= 50:
        ax.set_xticks(range(n_taxa))
        ax.set_xticklabels(props.columns, fontsize=6, rotation=90)
        ax.set_xlabel("Taxa")
    else:
        ax.set_xticks([])
        ax.set_xlabel(f"Taxa (n={n_taxa})")

    if show_boundaries and not params["shuffle_taxa"] and n_taxa <= 200 and n_samples_total <= 100:
        for x in params["boundary_lines"]:
            ax.axvline(x - 0.5, color="white", linewidth=1.5, linestyle="--", alpha=0.7)
        cum = 0
        for c in range(params["n_communities"] - 1):
            cum += params["n_samples"][c]
            ax.axhline(cum - 0.5, color="white", linewidth=1.5, linestyle="--", alpha=0.7)

    handles = [
        Patch(facecolor=comm_colors_list[c % len(comm_colors_list)], label=f"Community {c}")
        for c in range(params["n_communities"])
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.15, 1.0), fontsize=7, frameon=False)

    ax.set_ylabel("Samples")
    ax.set_title("Mock Community Heatmap")
    return fig, ax
