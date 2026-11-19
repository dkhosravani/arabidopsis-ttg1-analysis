"""Create publication-style plots directly from the reproducible result tables."""

from __future__ import annotations

import csv
import json
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .align_analysis import OUTPUT

PLOTS = OUTPUT / "figures"
GREEN = "#126f65"
TEAL = "#3fa594"
GOLD = "#da944c"
INK = "#1a2c30"


def load_csv(name: str):
    with (OUTPUT / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def export(fig, stem: str):
    fig.savefig(PLOTS / f"{stem}.png", dpi=210, bbox_inches="tight", facecolor="white")
    fig.savefig(PLOTS / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_gc():
    data = defaultdict(list)
    for row in load_csv("sliding_gc.csv"):
        data[row["accession"]].append((float(row["position_nt"]), float(row["gc_percent"])))
    fig, ax = plt.subplots(figsize=(9.1, 4.5))
    colors = [GREEN, GOLD, "#5586bd", "#a77eb4"]
    for (accession, points), color in zip(data.items(), colors):
        x, y = zip(*points)
        ax.plot(x, y, linewidth=1.9, label=accession, color=color)
    ax.set(title="TTG1 mRNA: sliding-window GC content", xlabel="mRNA position (nt, window midpoint)",
           ylabel="GC content (%)")
    ax.set_ylim(0, 70)
    ax.grid(alpha=0.15)
    ax.legend(title="NCBI RefSeq transcript", fontsize=8, frameon=False, loc="best")
    fig.text(0.11, -0.04, "100-nt windows, 10-nt step; transcript (not genomic DNA)",
             fontsize=9, color="#52636a")
    fig.tight_layout()
    export(fig, "gc_profile")


def plot_identity():
    rows = load_csv("pairwise_identity_matrix.csv")
    labels = list(rows[0])[1:]
    values = np.array([[float(row[label]) for label in labels] for row in rows])
    display = [label.replace("_", " ") for label in labels]
    fig, ax = plt.subplots(figsize=(8.2, 6.6))
    image = ax.imshow(values, cmap="YlGn", vmin=20, vmax=100)
    ax.set_xticks(range(len(labels)), display, rotation=35, ha="right")
    ax.set_yticks(range(len(labels)), display)
    ax.tick_params(axis="both", length=0)
    for i in range(len(labels)):
        for j in range(len(labels)):
            score = values[i, j]
            ax.text(j, i, f"{score:.1f}%", ha="center", va="center",
                    color="white" if score > 78 else INK, fontsize=10, fontweight="semibold")
    ax.set_title("TTG1-family proteins: pairwise sequence identity", pad=19)
    fig.colorbar(image, ax=ax, fraction=0.05, pad=0.03, label="Exact identity (%)")
    fig.text(0.12, 0.01, "Global BLOSUM62 alignment; identity = identical residues / all columns (gaps included)",
             fontsize=8, color="#52636a")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    export(fig, "identity_heatmap")


def plot_wd40():
    repeats = json.loads((OUTPUT / "wd40_annotations.json").read_text(encoding="utf-8"))
    info = json.loads((OUTPUT / "summary.json").read_text(encoding="utf-8"))
    length = info["reference_protein_length_aa"]
    fig, ax = plt.subplots(figsize=(10, 3.3))
    ax.plot([1, length], [0, 0], color="#a0b2ac", linewidth=12, solid_capstyle="round")
    for index, repeat in enumerate(repeats):
        left, right = repeat["start_aa"], repeat["end_aa"]
        ax.plot([left, right], [0, 0], linewidth=13,
                color=[GREEN, TEAL, GOLD, "#6085bb"][index % 4], solid_capstyle="butt")
        y = 0.24 if index % 2 == 0 else -0.29
        ax.annotate(f"{repeat['label']}\n{left}–{right}", xy=((left + right) / 2, 0),
                    xytext=((left + right) / 2, y), ha="center", va="center", fontsize=9,
                    arrowprops={"arrowstyle": "-", "color": "#82948e"})
    ax.set(xlim=(-10, length + 12), ylim=(-0.62, 0.60),
           xlabel="Amino acid coordinate (UniProt Q9XGN1)", title="Annotated WD repeats in Arabidopsis TTG1")
    ax.set_xticks([1, 50, 100, 150, 200, 250, 300, length])
    ax.set_yticks([])
    for edge in ["left", "top", "right"]:
        ax.spines[edge].set_visible(False)
    fig.tight_layout()
    export(fig, "wd40_map")


def plot_orfs():
    orfs = [row for row in load_csv("orf_candidates.csv") if row["accession"] == "NM_122360.2"]
    fig, ax = plt.subplots(figsize=(9.1, 3.4))
    for item in orfs:
        start, end, frame = int(item["start_nt"]), int(item["end_nt"]), int(item["frame"])
        ax.broken_barh([(start, end - start + 1)], (frame - 0.24, 0.48), facecolors=GREEN, alpha=0.80)
    ax.set(xlabel="mRNA position (nt)", ylabel="Reading frame (+ strand)",
           yticks=[1, 2, 3], title="TTG1 candidate complete ORFs (≥30 aa)")
    ax.grid(axis="x", alpha=0.15)
    ax.invert_yaxis()
    fig.text(0.12, 0.01, "Candidate ORFs are not a substitute for annotated CDS coordinates.",
             fontsize=9, color="#52636a")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    export(fig, "orf_candidates")


def main():
    PLOTS.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.titlesize": 13, "axes.titleweight": "bold",
                         "axes.labelcolor": INK, "text.color": INK})
    plot_gc()
    plot_identity()
    plot_wd40()
    plot_orfs()
    print(f"Created four figures (PNG and SVG) in {PLOTS.relative_to(PLOTS.parents[1])}")


if __name__ == "__main__":
    main()

