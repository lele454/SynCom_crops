#!/usr/bin/env python3
# [S6.4] Figures: data loading and saving helpers
# Supplementary 6 of the manuscript. See README.md for the full map.
"""
Figuras adicionais do manuscrito (mesma convencao de make_figures.py: sem titulos nem anotacoes
dentro dos paineis; so eixos, rotulos de categoria e legendas de cor; tudo explicado na legenda
do manuscrito). Cada figura sai em PNG (300 dpi) e PDF, e os dados vao para <saida>/dados/.

  Fig4_consorcios    saturacao por tamanho: n de consorcios, elo mais fraco e copresenca
  Fig5_rizobios      parceiros dos rizobios: efeito em leguminosas x gramineas
  Fig6_funcao        tracos PGP e clusters biossinteticos por genero dos consorcios
  FigS3_sequeiro_arroz   pares do nucleo: sequeiro x arroz
  FigS4_vies_filos   razao observado/esperado de associacoes por filo
  FigS5_moderadores  efeito sob manejo (N mineral x controle; seca x irrigado)

Uso:
  LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 python3 figures_extra.py \
      --root ~/syncom --out Figuras
"""
import argparse, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_figures import (BLUE, ORANGE, AQUA, VIOLET, RED, YELLOW, MAGENTA, GREEN, GRAY, LGRAY,
                          INK, MUTED, CROP_EN, log)

PHY_EN = {"Bacillota": "Bacillota", "Actinomycetota": "Actinomycetota"}


def save(fig, out, name, data=None, dd=None):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"))
    plt.close(fig)
    if data is not None and dd:
        data.to_csv(os.path.join(dd, f"{name}.tsv"), sep="\t", index=False, float_format="%.5g")
    log(f"  {name}")


def read(p, **kw):
    return pd.read_csv(p, sep="\t", **kw) if os.path.exists(p) else None


# ------------------------------------------------------------------ Fig4: saturacao
def fig_consorcios(mm, out, dd):
    S = read(os.path.join(mm, "consorcios", "estrelas_saturacao.tsv"))
    E = read(os.path.join(mm, "consorcios", "estrelas.tsv"))
    if S is None or E is None:
        return log("  Fig4: faltam estrelas_saturacao.tsv/estrelas.tsv")
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.3))
    x = S.tamanho.values
    ax[0].bar(x, S.n.values, color=BLUE, width=.62)
    ax[0].set_yscale("log"); ax[0].set_xlabel("Consortium size (genera)")
    ax[0].set_ylabel("Candidate consortia")
    ax[1].plot(x, S.rcond_min_melhor, "o-", color=ORANGE, lw=1.2, ms=4)
    ax[1].set_xlabel("Consortium size (genera)"); ax[1].set_ylabel("Best weakest link (r)")
    ax[1].set_ylim(0, None)
    bp = [E.copresenca[E.tamanho == t].values for t in x]
    b = ax[2].boxplot(bp, positions=x, widths=.55, patch_artist=True, showfliers=False,
                      medianprops=dict(color=INK, lw=1))
    for patch in b["boxes"]:
        patch.set(facecolor=LGRAY, edgecolor=MUTED, lw=.6)
    ax[2].plot(x, S.copres_melhor, "o", color=RED, ms=4)
    ax[2].set_xlabel("Consortium size (genera)"); ax[2].set_ylabel("Co-presence (fraction of samples)")
    ax[2].legend(handles=[Line2D([], [], marker="o", ls="", color=RED, label="Best consortium"),
                          Patch(facecolor=LGRAY, edgecolor=MUTED, label="All consortia")],
                 frameon=False, fontsize=6.5, loc="upper right")
    for a in ax:
        a.set_xticks(x)
    fig.tight_layout()
    save(fig, out, "Fig4_consorcios", S, dd)


# ------------------------------------------------------------------ Fig5: rizobios
def fig_rizobios(mm, out, dd, n=26):
    P = read(os.path.join(mm, "rizobios", "parceiros.tsv"))
    if P is None:
        return log("  Fig5: falta rizobios/parceiros.tsv")
    d = P[P.ok_leg & P.r_leg.notna()].copy()
    if d.empty:
        return log("  Fig5: nenhum parceiro confirmado")
    d["rotulo"] = d.rizobio + " – " + d.parceiro
    d = d.nlargest(min(n, len(d)), "r_leg").iloc[::-1]
    hubs = sorted(d.rizobio.unique())
    cor = dict(zip(hubs, [BLUE, ORANGE, AQUA, VIOLET, GREEN, YELLOW, MAGENTA, GRAY, RED][:len(hubs)]))
    fig, ax = plt.subplots(1, 2, figsize=(7.2, max(0.19 * len(d) + 1.1, 3.2)),
                           gridspec_kw={"width_ratios": [2.1, 1]}, sharey=True)
    y = np.arange(len(d))
    ax[0].barh(y, d.r_leg, color=[cor[h] for h in d.rizobio], height=.66)
    m = d.r_gram.notna().values
    ax[0].plot(d.r_gram.values[m], y[m], "o", ms=3.4, mfc="white", mec=INK, mew=.8)
    ax[0].set_yticks(y); ax[0].set_yticklabels(d.rotulo, fontsize=6.4, style="italic")
    ax[0].set_xlabel("Effect size in legumes (r)"); ax[0].axvline(0, color=MUTED, lw=.6)
    ax[0].legend(handles=[Patch(facecolor=cor[h], label=h) for h in hubs] +
                         [Line2D([], [], marker="o", ls="", mfc="white", mec=INK, label="Grasses")],
                 frameon=False, fontsize=6.2, loc="lower right")
    ax[1].barh(y, d.copresenca_leg, color=LGRAY, edgecolor=MUTED, lw=.5, height=.66)
    ax[1].set_xlabel("Co-presence in legumes")
    fig.tight_layout()
    save(fig, out, "Fig5_rizobios", d, dd)


# ------------------------------------------------------------------ Fig6: funcao
def fig_funcao(root, mm, out, dd):
    T = read(os.path.join(root, "genomas", "pgp", "tracos_por_genero.tsv"), index_col=0)
    B = read(os.path.join(root, "genomas", "bgc", "bgc_por_genero.tsv"), index_col=0)
    R = read(os.path.join(mm, "consorcios", "ranking_aptos.tsv"), index_col=0)
    if T is None or R is None:
        return log("  Fig6: faltam tracos_por_genero.tsv/ranking_aptos.tsv")
    gen = []
    for m in R.membros.head(12):
        gen += [x.strip() for x in str(m).split("+")]
    gen = [g for g in dict.fromkeys(gen) if g in T.index]
    if not gen:
        return log("  Fig6: nenhum genero com tracos")
    tr = [c for c in T.columns if c != "n_genomas"]
    M = T.loc[gen, tr]
    cols = list(M.columns)
    if B is not None:
        tipos = [c for c in B.columns if c.startswith("BGC_") and c not in ("BGC_total", "BGC_Unknown")]
        M = M.join(B.loc[B.index.intersection(gen), tipos], how="left")
        cols = cols + tipos
    M = M.reindex(columns=cols).fillna(0)
    nice = {c: c.replace("_", " ").replace("BGC ", "") for c in M.columns}
    cmap = LinearSegmentedColormap.from_list("seq", ["#f4f3f0", BLUE])
    fig, ax = plt.subplots(figsize=(max(0.34 * len(M.columns) + 2.2, 5.6), 0.26 * len(M) + 1.5))
    im = ax.imshow(M.values, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(M.columns)))
    ax.set_xticklabels([nice[c] for c in M.columns], rotation=55, ha="right", fontsize=6.4)
    ax.set_yticks(range(len(M))); ax.set_yticklabels(M.index, fontsize=6.6, style="italic")
    ax.set_xticks(np.arange(-.5, len(M.columns), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(M), 1), minor=True)
    ax.grid(which="minor", color="white", lw=.8); ax.tick_params(which="minor", length=0)
    cb = fig.colorbar(im, ax=ax, fraction=.022, pad=.02)
    cb.set_label("Fraction of genomes", fontsize=6.6); cb.ax.tick_params(labelsize=6)
    fig.tight_layout()
    save(fig, out, "Fig6_funcao", M.reset_index(), dd)


# ------------------------------------------------- FigS3: sequeiro x arroz
def fig_sequeiro_arroz(mm, out, dd):
    A = read(os.path.join(mm, "arroz", "meta_pairs_arroz.tsv.gz"))
    M = read(os.path.join(mm, "meta_pairs.tsv.gz"))
    if A is None or M is None:
        return log("  FigS3: faltam meta_pairs_arroz.tsv.gz/meta_pairs.tsv.gz")
    key = lambda d: d.genus_a + "|" + d.genus_b
    A.index = key(A); M.index = key(M)
    both = A.index.intersection(M.index)
    d = pd.DataFrame({"r_total": M.loc[both, "r_mu"], "r_arroz": A.loc[both, "r_mu"],
                      "sig_total": M.loc[both, "fdr"] < 0.05, "sig_arroz": A.loc[both, "fdr"] < 0.05})
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw={"width_ratios": [1.15, 1]})
    ax[0].axhline(0, color=LGRAY, lw=.6); ax[0].axvline(0, color=LGRAY, lw=.6)
    ax[0].plot([-.6, .8], [-.6, .8], color=GRAY, lw=.6, ls="--")
    ax[0].scatter(d.r_total[~(d.sig_total & d.sig_arroz)], d.r_arroz[~(d.sig_total & d.sig_arroz)],
                  s=3, color=LGRAY, lw=0)
    s = d.sig_total & d.sig_arroz
    ax[0].scatter(d.r_total[s], d.r_arroz[s], s=5, color=BLUE, lw=0)
    ax[0].set_xlabel("Effect size, all crops (r)"); ax[0].set_ylabel("Effect size, paddy rice (r)")
    ax[0].legend(handles=[Line2D([], [], marker="o", ls="", color=BLUE, ms=3.6,
                                 label="Significant in both"),
                          Line2D([], [], marker="o", ls="", color=LGRAY, ms=3.6, label="Other pairs")],
                 frameon=False, fontsize=6.5, loc="upper left")
    cat = ["Both", "All crops only", "Rice only"]
    v = [int((d.sig_total & d.sig_arroz).sum()), int((d.sig_total & ~d.sig_arroz).sum()),
         int((~d.sig_total & d.sig_arroz).sum())]
    ax[1].bar(cat, v, color=[BLUE, GRAY, AQUA], width=.6)
    ax[1].set_ylabel("Consistent pairs")
    ax[1].tick_params(axis="x", labelsize=6.8)
    fig.tight_layout()
    save(fig, out, "FigS3_sequeiro_arroz", d.reset_index(names="pair"), dd)


# ------------------------------------------------- FigS4: vies por filo
def fig_vies_filos(mm, out, dd):
    P = read(os.path.join(mm, "vies_filos", "por_filo.tsv"))
    if P is None:
        return log("  FigS4: falta vies_filos/por_filo.tsv")
    d = P.sort_values("razao_obs_esp")
    fig, ax = plt.subplots(figsize=(4.6, 0.23 * len(d) + 1.2))
    y = np.arange(len(d))
    cor = [BLUE if v >= 1 else RED for v in d.razao_obs_esp]
    ax.barh(y, d.razao_obs_esp, color=cor, height=.66)
    ax.axvline(1, color=MUTED, lw=.7)
    ax.set_yticks(y); ax.set_yticklabels(d.filo, fontsize=6.8, style="italic")
    ax.set_xlabel("Observed / expected associations")
    if "ssu_mediano" in d:
        ax2 = ax.twiny()
        ax2.plot(d.ssu_mediano, y, "o", ms=3.4, mfc="white", mec=INK, mew=.8)
        ax2.set_xlabel("16S copies per genome (median)", fontsize=7)
        ax2.spines["top"].set_visible(True)
    fig.tight_layout()
    save(fig, out, "FigS4_vies_filos", d, dd)


# ------------------------------------------------- FigS5: moderadores
def fig_moderadores(mm, out, dd):
    P = read(os.path.join(mm, "moderadores", "por_par.tsv.gz"))
    if P is None:
        return log("  FigS5: falta moderadores/por_par.tsv.gz")
    cs = [c for c in P.contraste.unique()]
    fig, ax = plt.subplots(1, len(cs), figsize=(3.4 * len(cs), 3.0), squeeze=False)
    for a, c in zip(ax[0], cs):
        d = P[P.contraste == c]
        n1, n2 = c.split("_vs_")
        c1, c2 = f"r_{n1}", f"r_{n2}"
        if c1 not in d or c2 not in d:
            continue
        lim = np.nanmax(np.abs(np.r_[d[c1].values, d[c2].values])) * 1.1
        a.plot([-lim, lim], [-lim, lim], color=GRAY, lw=.6, ls="--")
        a.axhline(0, color=LGRAY, lw=.6); a.axvline(0, color=LGRAY, lw=.6)
        sig = d.fdr < 0.05
        a.scatter(d[c2][~sig], d[c1][~sig], s=4, color=LGRAY, lw=0)
        a.scatter(d[c2][sig], d[c1][sig], s=7, color=RED, lw=0)
        a.set_xlabel(f"Effect size, {n2.replace('_', ' ')} (r)")
        a.set_ylabel(f"Effect size, {n1.replace('_', ' ')} (r)")
        a.set_xlim(-lim, lim); a.set_ylim(-lim, lim); a.set_aspect("equal")
        if sig.any():
            a.legend(handles=[Line2D([], [], marker="o", ls="", color=RED, ms=3.6,
                                     label="Significant moderation")], frameon=False,
                     fontsize=6.5, loc="upper left")
    fig.tight_layout()
    save(fig, out, "FigS5_moderadores", P, dd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", default="Figuras")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    mm = os.path.join(root, "merged", "meta")
    out = os.path.abspath(os.path.expanduser(a.out)); dd = os.path.join(out, "dados")
    os.makedirs(dd, exist_ok=True)
    log("Figuras adicionais:")
    fig_consorcios(mm, out, dd)
    fig_rizobios(mm, out, dd)
    fig_funcao(root, mm, out, dd)
    fig_sequeiro_arroz(mm, out, dd)
    fig_vies_filos(mm, out, dd)
    fig_moderadores(mm, out, dd)
    log(f"Figuras em {out}")


if __name__ == "__main__":
    main()
