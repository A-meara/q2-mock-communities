import os

import biom
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def summarize(output_dir: str, proportions: biom.Table, community_map: pd.DataFrame,
              log_scale: bool = True, show_boundaries: bool = True) -> None:
    """Heatmap of community proportions with sample community labels."""
    # BIOM → DataFrame (samples × taxa)
    props_df = pd.DataFrame(
        proportions.matrix_data.toarray().T,
        index=proportions.ids('sample'),
        columns=proportions.ids('observation'),
    )

    n_samples, n_taxa = props_df.shape
    w = min(max(6, n_taxa * 0.25), 20)
    h = min(max(4, n_samples * 0.35), 20)
    figsize = (w, h)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    data = props_df.values.copy()
    if log_scale:
        data = np.log10(data + 1e-6)
        cbar_label = "log₁₀(proportion + 1e-6)"
    else:
        cbar_label = "proportion"

    im = ax.imshow(data, aspect="auto", cmap="viridis", interpolation="nearest")
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label(cbar_label)

    # Community color palette
    comm_colors = plt.cm.Set2.colors
    communities = community_map.reindex(props_df.index)['community'].values

    unique_comms = sorted(set(communities))
    comm_to_idx = {c: i for i, c in enumerate(unique_comms)}

    if n_samples <= 50:
        ax.set_yticks(range(n_samples))
        ax.set_yticklabels(props_df.index, fontsize=7)
        for i, c in enumerate(communities):
            ax.get_yticklabels()[i].set_color(
                comm_colors[comm_to_idx[c] % len(comm_colors)]
            )
    else:
        # One label per community at the midpoint of its block
        positions, labels, label_comms = [], [], []
        cum = 0
        prev = communities[0]
        block_start = 0
        for i in range(1, len(communities)):
            if communities[i] != prev:
                mid = block_start + (i - block_start - 1) / 2.0
                positions.append(mid)
                labels.append(f"Community {prev}  (n={i - block_start})")
                label_comms.append(prev)
                block_start = i
                prev = communities[i]
        # last block
        mid = block_start + (n_samples - block_start - 1) / 2.0
        positions.append(mid)
        labels.append(f"Community {prev}  (n={n_samples - block_start})")
        label_comms.append(prev)
        ax.set_yticks(positions)
        ax.set_yticklabels(labels, fontsize=8)
        for i, c in enumerate(label_comms):
            ax.get_yticklabels()[i].set_color(
                comm_colors[comm_to_idx[c] % len(comm_colors)]
            )

    if n_taxa <= 50:
        ax.set_xticks(range(n_taxa))
        ax.set_xticklabels(props_df.columns, fontsize=6, rotation=90)
    else:
        ax.set_xticks([])

    # Horizontal boundaries between community blocks (suppress at large scale)
    if show_boundaries and n_taxa <= 200 and n_samples <= 100:
        prev = communities[0]
        for i in range(1, len(communities)):
            if communities[i] != prev:
                ax.axhline(i - 0.5, color='white', linewidth=1.5,
                           linestyle='--', alpha=0.7)
            prev = communities[i]

    handles = [
        Patch(facecolor=comm_colors[comm_to_idx[c] % len(comm_colors)],
              label=f"Community {c}")
        for c in unique_comms
    ]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.15, 1.0),
              fontsize=7, frameon=False)

    ax.set_xlabel("Taxa")
    ax.set_ylabel("Samples")
    ax.set_title("Mock Community Heatmap")
    fig.tight_layout()

    fig.savefig(os.path.join(output_dir, 'heatmap.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)

    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write(
            '<!DOCTYPE html><html><head><title>Mock Community Summary</title>'
            '<style>body{font-family:sans-serif;padding:20px}'
            'img{max-width:100%;height:auto}</style></head>'
            '<body><h2>Mock Community Heatmap</h2>'
            '<img src="heatmap.png" alt="Community heatmap"/>'
            '</body></html>'
        )
