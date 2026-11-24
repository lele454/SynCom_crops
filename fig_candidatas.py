#!/usr/bin/env python3
# [S6.5] Figures of the candidate combinations
# Supplementary 6 of the manuscript. See README.md for the full map.
"""Fig: combinacoes testaveis in planta (usuais x nao usuais) e parceiros de inoculantes."""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from make_figures import BLUE, ORANGE, LGRAY, INK
from figures_extra import read, save

USUAL = ["Arthrobacter", "Sphingomicrobium", "Streptomyces"]
NOVO = ["Devosia_A", "Phenylobacterium", "Ramlibacter"]


def painel(ax, letra):
    ax.set_title(f"{letra})", loc="left", fontsize=9, fontweight="bold", color=INK, pad=6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", default="Figuras")
    ap.add_argument("--n-pares", type=int, default=12)
    a = ap.parse_args()
    MM = os.path.join(a.root, "merged", "meta")
    S = read(os.path.join(MM, "consorcios", "candidatas_testaveis.tsv"))
    Q = read(os.path.join(MM, "consorcios", "parceiros_elegiveis.tsv"))
    if S is None:
        raise SystemExit("falta candidatas_testaveis.tsv")
    S["chave"] = S.membros.map(lambda x: " + ".join(sorted(str(x).split(" + "))))
    rec = {" + ".join(sorted(USUAL)), " + ".join(sorted(NOVO))}

    fig = plt.figure(figsize=(7.2, 7.2))

    axA = fig.add_axes([0.11, 0.60, 0.84, 0.34])
    for cond, cor, z in ((S.n_pgp > 0, LGRAY, 2), (S.n_pgp == 0, ORANGE, 3)):
        d = S[cond]
        axA.scatter(d.copresenca, d.r_min, s=6 + d.estudos_copresente * 0.55,
                    color=cor, alpha=.55, lw=0, zorder=z)
    d = S[S.chave.isin(rec)]
    axA.scatter(d.copresenca, d.r_min, s=6 + d.estudos_copresente * 0.55,
                facecolor="none", edgecolor=BLUE, lw=1.5, zorder=5)
    axA.set_xlabel("Co-presence (fraction of samples)")
    axA.set_ylabel("Weakest internal association (r)")
    axA.grid(color="#e8e6e1", lw=.6); axA.set_axisbelow(True)
    axA.legend(handles=[Line2D([], [], marker="o", ls="", color=LGRAY,
                               label="With a reported growth promoter"),
                        Line2D([], [], marker="o", ls="", color=ORANGE,
                               label="No reported growth promoter"),
                        Line2D([], [], marker="o", ls="", mfc="none", mec=BLUE,
                               label="Recommended for testing")],
               frameon=False, fontsize=6.6, loc="upper left")
    painel(axA, "a")

    axB = fig.add_axes([0.34, 0.07, 0.61, 0.40])
    if Q is not None and len(Q):
        q = Q.sort_values("copresenca_leg").tail(a.n_pares)
        y = np.arange(len(q))
        cor = [ORANGE if g == "rizobios" else BLUE for g in q.grupo]
        axB.hlines(y, 0, q.copresenca_leg, color=cor, lw=1.8, alpha=.8)
        axB.scatter(q.copresenca_leg, y, color=cor, s=30, zorder=4)
        axB.set_yticks(y, [f"{c} \u2013 {p}" for c, p in zip(q.rizobio, q.parceiro)],
                       fontsize=7, style="italic")
        axB.set_ylim(-0.8, len(q) - 0.2)
        axB.set_xlabel("Co-presence within the host group")
        axB.set_xlim(0, max(q.copresenca_leg) * 1.08)
        axB.grid(axis="x", color="#e8e6e1", lw=.6); axB.set_axisbelow(True)
        for s in ("top", "right", "left"):
            axB.spines[s].set_visible(False)
        axB.tick_params(axis="y", length=0)
    painel(axB, "b")

    save(fig, a.out, "Fig2_candidatas", S, a.out)
    print(f"Fig2_candidatas: {len(S)} candidatas | pares: {0 if Q is None else len(Q)}")


if __name__ == "__main__":
    main()
