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
    groups: list[list[int]] | None = None,
    group_core_overlap: float | int = 0.0,
    alpha: float = 0.5,
    core_signature_strength: float = 1.0,
    core_weight: float | None = None,
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

    if groups is None:
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

        n_group_core_list: list[int] = []

    else:
        # ---- groups layout ----
        all_in_groups = sorted(c for g in groups for c in g)
        if all_in_groups != list(range(n_communities)):
            raise ValueError(
                "groups must partition all community indices 0..n_communities-1 exactly once"
            )

        community_to_group = {c: g_idx for g_idx, g in enumerate(groups) for c in g}

        n_group_core_list = []
        for g in groups:
            if len(g) <= 1:
                n_group_core_list.append(0)
            elif isinstance(group_core_overlap, int):
                n_group_core_list.append(group_core_overlap)
            else:
                n_group_core_list.append(
                    round(group_core_overlap * min(taxa_per[c] for c in g))
                )

        for c in range(n_communities):
            g_idx = community_to_group[c]
            n_unique_c = taxa_per[c] - n_core - n_group_core_list[g_idx]
            if n_unique_c < 1:
                raise ValueError(
                    f"community {c}: taxa_per={taxa_per[c]} - n_core={n_core} "
                    f"- n_group_core={n_group_core_list[g_idx]} = {n_unique_c} < 1"
                )

        taxa_offset = n_core
        boundary_lines_list: list[int] = []
        if n_core > 0:
            boundary_lines_list.append(n_core)

        group_core_idx_map: dict[int, list[int]] = {}
        community_unique_idx_map: dict[int, list[int]] = {}
        overlap_counts: list[int] = []

        for g_idx, g in enumerate(groups):
            gc = n_group_core_list[g_idx]
            group_core_idx_map[g_idx] = list(range(taxa_offset, taxa_offset + gc))
            taxa_offset += gc
            if gc > 0:
                boundary_lines_list.append(taxa_offset)

            unique_sizes = [taxa_per[c] - n_core - gc for c in g]

            if overlap_mode in ('chain', 'both') and len(g) > 1:
                ov_g = [
                    round(overlap * min(unique_sizes[i], unique_sizes[i + 1]))
                    for i in range(len(g) - 1)
                ]
                for i, ov in enumerate(ov_g):
                    if unique_sizes[i] - ov < 1:
                        raise ValueError(
                            f"overlap={overlap} between communities {g[i]} and {g[i+1]} "
                            f"in group {g_idx} leaves stride<1"
                        )
                overlap_counts.extend(ov_g)
                u_starts = [taxa_offset]
                for i in range(len(g) - 1):
                    u_starts.append(u_starts[-1] + unique_sizes[i] - ov_g[i])
                taxa_offset = u_starts[-1] + unique_sizes[-1]
                boundary_lines_list.extend(u_starts[1:])
            else:
                u_starts = [taxa_offset + sum(unique_sizes[:i]) for i in range(len(g))]
                taxa_offset += sum(unique_sizes)
                if len(g) > 1:
                    boundary_lines_list.extend(u_starts[1:])

            for i, c in enumerate(g):
                community_unique_idx_map[c] = list(
                    range(u_starts[i], u_starts[i] + unique_sizes[i])
                )

            if g_idx < len(groups) - 1:
                boundary_lines_list.append(taxa_offset)

        total_taxa = taxa_offset
        all_taxa = [f"taxon_{i + 1:03d}" for i in range(total_taxa)]
        taxon_membership: list[list[int]] = [[] for _ in range(total_taxa)]
        community_indices: list = [None] * n_communities

        for t in range(n_core):
            taxon_membership[t] = list(range(n_communities))

        for g_idx, g in enumerate(groups):
            gc_idx = group_core_idx_map[g_idx]
            for t in gc_idx:
                taxon_membership[t] = list(g)
            for c in g:
                u_idx = community_unique_idx_map[c]
                community_indices[c] = np.array(list(range(n_core)) + gc_idx + u_idx)
                for t in u_idx:
                    taxon_membership[t].append(c)

        boundary_lines = sorted(set(boundary_lines_list))

    # --- core signature matrix ---
    use_signatures = core_signature_strength != 1.0 and n_core > 0
    if use_signatures:
        sig = np.ones((n_communities, n_core))
        for c, chunk in enumerate(np.array_split(range(n_core), n_communities)):
            sig[c, chunk] = core_signature_strength

    # --- generate mean proportions per community ---
    community_means = np.zeros((n_communities, total_taxa))
    for c in range(n_communities):
        idx = community_indices[c]
        if use_signatures:
            core_pos = idx[:n_core]
            other_pos = idx[n_core:]
            w = core_weight if core_weight is not None else n_core / len(idx)
            community_means[c, core_pos] = rng.dirichlet(sig[c]) * w
            if len(other_pos) > 0:
                community_means[c, other_pos] = rng.dirichlet(np.full(len(other_pos), alpha)) * (1 - w)
        else:
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
            "groups": groups,
            "group_core_overlap": group_core_overlap,
            "n_group_core": n_group_core_list,
            "overlap_counts": overlap_counts,
            "boundary_lines": boundary_lines,
            "alpha": alpha,
            "core_signature_strength": core_signature_strength,
            "core_weight": core_weight,
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
