#!/usr/bin/env python3
# [S6.6] Figures of the trios
# Supplementary 6 of the manuscript. See README.md for the full map.
"""Evidencia dos trios recomendados: fracao de amostras com todos os membros, por estudo."""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from make_figures import BLUE, ORANGE, LGRAY, INK
from figures_extra import read, save

ROT = {"usual": "Arthrobacter + Sphingomicrobium + Streptomyces",
       "novo": "Devosia_A + Phenylobacterium + Ramlibacter"}


def painel(ax, letra):
    ax.set_title(f"{letra})", loc="left", fontsize=9, fontweight="bold", color=INK, pad=6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", default="Figuras")
    a = ap.parse_args()
    T = read(os.path.join(a.root, "merged", "meta", "consorcios", "trio_por_estudo.tsv"))
    if T is None:
        raise SystemExit("falta trio_por_estudo.tsv")
    fig = plt.figure(figsize=(7.2, 4.6))

    axA = fig.add_axes([0.10, 0.58, 0.85, 0.34])
    for nome, cor in (("usual", BLUE), ("novo", ORANGE)):
        v = T[(T.trio == nome)].frac.dropna().sort_values(ascending=False).values
        axA.plot(np.arange(1, len(v) + 1), v, color=cor, lw=1.6)
    axA.axhline(0.5, color=LGRAY, lw=.9, ls="--", zorder=1)
    axA.set_xlabel("Studies, ranked")
    axA.set_ylabel("Samples with all three genera")
    axA.set_ylim(0, 1); axA.grid(color="#e8e6e1", lw=.6); axA.set_axisbelow(True)
    axA.legend(handles=[Line2D([], [], color=BLUE, label=ROT["usual"]),
                        Line2D([], [], color=ORANGE, label=ROT["novo"])],
               frameon=False, fontsize=6.6, loc="lower left")
    painel(axA, "a")

    axB = fig.add_axes([0.10, 0.10, 0.85, 0.30])
    ordem = sorted(T.cultura.dropna().unique())
    pos = {c: i for i, c in enumerate(ordem)}
    for nome, cor, dx in (("usual", BLUE, -0.16), ("novo", ORANGE, 0.16)):
        d = T[(T.trio == nome)].dropna(subset=["frac"])
        x = [pos[c] + dx + np.random.RandomState(1).uniform(-.05, .05) for c in d.cultura]
        axB.scatter(x, d.frac, color=cor, s=12, alpha=.6, lw=0)
        med = d.groupby("cultura").frac.median()
        for c, m in med.items():
            axB.plot([pos[c] + dx - .12, pos[c] + dx + .12], [m, m], color=cor, lw=2)
    axB.set_xticks(range(len(ordem)), ordem, fontsize=7)
    axB.set_ylabel("Samples with all three genera")
    axB.set_ylim(0, 1); axB.grid(axis="y", color="#e8e6e1", lw=.6); axB.set_axisbelow(True)
    painel(axB, "b")

    save(fig, a.out, "Fig3_trios", T, a.out)
    print(f"Fig3_trios: {T.study.nunique()} estudos")


if __name__ == "__main__":
    main()
