#!/usr/bin/env python3
# [S6.9] Final body and supplementary figures
# Supplementary 6 of the manuscript. See README.md for the full map.
"""
Conjunto final de figuras, depois da reorganizacao dos Resultados.

  CORPO
  Fig3_corpo  a) isolados x clinicos  b) classe por grupo
              c) copresenca dos tres conjuntos  d) matriz do conjunto executavel
  Fig4_corpo  a) trios por estudo  b) trios por cultura
              c) copresenca por tamanho  d) compromisso copresenca x complementaridade

  SUPLEMENTAR
  FigS13_empates     frequencia de cada genero entre os cliques maximos
  FigS14_remocao     a) remocao in silico  b) estruturacao da comunidade
  FigS15_parceiros   pares inoculante-parceiro por copresenca

Saem do conjunto: a rede de consenso (hairball, redundante com S6 e S7), o painel de
dispersao das candidatas e o funil de quatro degraus, que contradiz o funil de sete
degraus descrito nos Resultados.

Uso: LD_PRELOAD=... python3 figs_final.py --root ~/syncom --out Figuras
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
from make_figures import (BLUE, ORANGE, AQUA, VIOLET, RED, YELLOW, GREEN,
                          GRAY, LGRAY, INK, MUTED, log, grid)
from figures_extra import read, save
from figures_new import painel, CULT_EN

ROT = {"usual": "Arthrobacter + Sphingomicrobium + Streptomyces",
       "novo": "Devosia_A + Phenylobacterium + Ramlibacter"}
# gramineas primeiro, leguminosas depois, como nas colunas da Fig. 2
ORDEM_CULT = ["milho", "arroz", "sorgo", "cana", "trigo", "soja", "feijao/leguminosas"]


def grupos_genero(root, mm):
    """Tabela por genero com o grupo (parceiro de rizobio, de PGP, ou outro)."""
    G = read(os.path.join(root, "genomas", "isolados", "por_genero.tsv"), index_col=0)
    if G is None:
        return None
    grupo = pd.Series("Other candidate genera", index=G.index)
    for sub, nome in (("rizobios", "Rhizobial partners"), ("pgp", "PGP partners")):
        P = read(os.path.join(mm, sub, "parceiros.tsv"))
        if P is None:
            continue
        p = P[P.ok_leg] if "ok_leg" in P else P
        grupo.loc[grupo.index.intersection(set(p.parceiro))] = nome
    G = G.assign(grupo=grupo)
    iso = G.genomas_isolados.astype(float)
    G["frac_clinica"] = np.where(iso > 0, G.isolados_clinicos.astype(float) / iso.clip(lower=1), np.nan)
    return G


# --------------------------------------------- Fig3: linhagens e nucleo recorrente
def fig3_corpo(root, mm, out, dd):
    G = grupos_genero(root, mm)
    md = os.path.join(mm, "modulo")
    P = read(os.path.join(md, "copresenca.tsv"))
    Q = read(os.path.join(md, "pares.tsv"))
    if G is None or P is None or Q is None:
        return log("  Fig3: faltam por_genero.tsv ou merged/meta/modulo/")

    fig = plt.figure(figsize=(7.2, 6.0))
    COR = {"Rhizobial partners": ORANGE, "PGP partners": BLUE, "Other candidate genera": GRAY}

    axA = fig.add_axes([0.09, 0.60, 0.37, 0.33])
    for g, c in COR.items():
        d = G[G.grupo == g]
        axA.scatter(d.frac_isolados, d.frac_clinica,
                    s=6 + 14 * np.log10(d.genomas.clip(lower=1) + 1),
                    color=c, alpha=.75, lw=.3, edgecolor="white", zorder=3)
    axA.axvline(.5, color=MUTED, lw=.6, ls="--"); axA.axhline(.25, color=MUTED, lw=.6, ls="--")
    axA.set_xlabel("Genomes from isolates (fraction)")
    axA.set_ylabel("Isolates of clinical origin (fraction)")
    axA.set_xlim(-0.03, 1.03); axA.set_ylim(-0.03, None)
    axA.legend(handles=[Line2D([], [], marker="o", ls="", ms=5, mfc=c, mec="white", label=g)
                        for g, c in COR.items()], frameon=False, fontsize=6.3,
               loc="upper left", bbox_to_anchor=(0, 1.0))
    grid(axA, "both"); painel(axA, "a")

    ordem = ["acessivel", "dificil", "cultivado_sem_genoma", "conferir_fatia", "sem_cultivado"]
    rot = {"acessivel": "Accessible", "dificil": "Few isolates",
           "cultivado_sem_genoma": "Cultured, no isolate genome",
           "conferir_fatia": "Check circumscription", "sem_cultivado": "No cultured representative"}
    cores = [BLUE, YELLOW, AQUA, VIOLET, RED]
    grupos = ["PGP partners", "Rhizobial partners", "Other candidate genera"]
    T = (G.groupby(["grupo", "classe"]).size().unstack(fill_value=0)
         .reindex(index=grupos, columns=ordem, fill_value=0))
    frac = T.div(T.sum(axis=1).clip(lower=1), axis=0)
    axB = fig.add_axes([0.70, 0.60, 0.26, 0.33])
    y = np.arange(len(grupos))[::-1]
    left = np.zeros(len(grupos))
    for cl, c in zip(ordem, cores):
        axB.barh(y, frac[cl].values, left=left, color=c, height=.6, edgecolor="white", lw=.5)
        left += frac[cl].values
    axB.set_yticks(y, [f"{g}\n(n = {int(T.loc[g].sum())})" for g in grupos], fontsize=6.4)
    axB.set_xlabel("Genera (fraction)"); axB.set_xlim(0, 1)
    axB.tick_params(axis="y", length=0); axB.spines["left"].set_visible(False)
    axB.legend(handles=[Patch(facecolor=c, label=rot[cl]) for cl, c in zip(ordem, cores)],
               frameon=False, fontsize=6.2, ncol=2, loc="upper center", bbox_to_anchor=(.5, -.18))
    painel(axB, "b")

    axC = fig.add_axes([0.09, 0.08, 0.37, 0.33])
    estilo = {"dez": (GRAY, "o--", "Densest set"),
              "replicado": (ORANGE, "o-", "Best replicated"),
              "acessivel": (BLUE, "o-", "Replicated and accessible")}
    hs = []
    for nome, (cor, st, rotulo) in estilo.items():
        d = P[P.conjunto == nome]
        if d.empty:
            continue
        axC.plot(d.tamanho, d.copres_media, st, color=cor, lw=1.3, ms=4.2)
        hs.append(Line2D([], [], color=cor, marker="o", ms=4,
                         ls="--" if "--" in st else "-", label=rotulo))
    axC.set_xlabel("Set size (genera)")
    axC.set_ylabel("Co-presence (fraction of samples)")
    axC.set_ylim(0, None)
    axC.legend(handles=hs, frameon=False, fontsize=6.4, loc="upper right")
    grid(axC, "y"); painel(axC, "c")

    A = Q[Q.conjunto == "acessivel"]
    g = sorted(set(A.genus_a) | set(A.genus_b))
    n = len(g)
    Z = np.full((n, n), np.nan)
    for _, r in A.iterrows():
        i1, j1 = g.index(r.genus_a), g.index(r.genus_b)
        Z[i1, j1] = Z[j1, i1] = r.r_mu
    cmap = LinearSegmentedColormap.from_list("seq", ["#f4f3f0", BLUE])
    axD = fig.add_axes([0.61, 0.10, 0.27, 0.31])
    im = axD.imshow(np.ma.masked_invalid(Z), cmap=cmap, vmin=0,
                    vmax=max(0.35, np.nanmax(Z) if n else 0.35))
    axD.set_xticks(range(n), g, rotation=55, ha="right", fontsize=6.3, style="italic")
    axD.set_yticks(range(n), g, fontsize=6.3, style="italic")
    axD.set_xticks(np.arange(-.5, n, 1), minor=True)
    axD.set_yticks(np.arange(-.5, n, 1), minor=True)
    axD.grid(which="minor", color="white", lw=.8); axD.tick_params(which="minor", length=0)
    cb = fig.colorbar(im, ax=axD, fraction=.045, pad=.04)
    cb.set_label("Pooled effect (r)", fontsize=6.6); cb.ax.tick_params(labelsize=6)
    painel(axD, "d")

    save(fig, out, "Fig3_corpo", G.reset_index(), dd)
    log(f"  Fig3_corpo: conjunto executavel = {' + '.join(g)}")


# ------------------------------------------- Fig4: combinacoes e limite de tamanho
def fig4_corpo(mm, out, dd):
    T = read(os.path.join(mm, "consorcios", "trio_por_estudo.tsv"))
    S = read(os.path.join(mm, "consorcios", "estrelas_saturacao.tsv"))
    E = read(os.path.join(mm, "consorcios", "estrelas.tsv"))
    C = read(os.path.join(mm, "consorcios", "complementaridade.tsv"))
    if T is None or S is None or E is None:
        return log("  Fig4: faltam trio_por_estudo.tsv ou estrelas*.tsv")

    fig = plt.figure(figsize=(7.2, 6.0))

    axA = fig.add_axes([0.09, 0.73, 0.87, 0.22])
    for nome, cor in (("usual", BLUE), ("novo", ORANGE)):
        v = T[T.trio == nome].frac.dropna().sort_values(ascending=False).values
        axA.plot(np.arange(1, len(v) + 1), v, color=cor, lw=1.6)
    axA.axhline(0.5, color=LGRAY, lw=.9, ls="--", zorder=1)
    axA.set_xlabel("Studies, ranked")
    axA.set_ylabel("Samples with all three genera")
    axA.set_ylim(0, 1)
    axA.legend(handles=[Line2D([], [], color=BLUE, label=ROT["usual"]),
                        Line2D([], [], color=ORANGE, label=ROT["novo"])],
               frameon=False, fontsize=6.6, loc="lower left")
    grid(axA, "y"); painel(axA, "a")

    axB = fig.add_axes([0.09, 0.42, 0.87, 0.20])
    presentes = set(T.cultura.dropna().unique())
    ordem = [c for c in ORDEM_CULT if c in presentes]
    ordem += sorted(presentes - set(ordem))
    pos = {c: i for i, c in enumerate(ordem)}
    rng = np.random.RandomState(1)
    for nome, cor, dx in (("usual", BLUE, -0.16), ("novo", ORANGE, 0.16)):
        d = T[T.trio == nome].dropna(subset=["frac"])
        x = [pos[c] + dx + rng.uniform(-.05, .05) for c in d.cultura]
        axB.scatter(x, d.frac, color=cor, s=12, alpha=.6, lw=0)
        for c, m in d.groupby("cultura").frac.median().items():
            axB.plot([pos[c] + dx - .12, pos[c] + dx + .12], [m, m], color=cor, lw=2)
    axB.set_xticks(range(len(ordem)), [CULT_EN.get(c, c) for c in ordem], fontsize=7)
    axB.set_ylabel("Samples with all three genera")
    axB.set_ylim(0, 1)
    grid(axB, "y"); painel(axB, "b")

    axC = fig.add_axes([0.09, 0.07, 0.37, 0.25])
    x = S.tamanho.values
    b = axC.boxplot([E.copresenca[E.tamanho == t].values for t in x], positions=x, widths=.55,
                    patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1))
    for p in b["boxes"]:
        p.set(facecolor=LGRAY, edgecolor=MUTED, lw=.6)
    axC.plot(x, S.copres_melhor, "o", color=RED, ms=4, zorder=4)
    axC.set_xticks(x); axC.set_xlabel("Combination size (genera)")
    axC.set_ylabel("Co-presence (fraction of samples)")
    axC.legend(handles=[Line2D([], [], marker="o", ls="", color=RED, label="Best combination"),
                        Patch(facecolor=LGRAY, edgecolor=MUTED, label="All combinations")],
               frameon=False, fontsize=6.5, loc="upper right")
    grid(axC, "y"); painel(axC, "c")

    axD = fig.add_axes([0.60, 0.07, 0.33, 0.25])
    K = C[C.n_com_tracos == C.tamanho] if (C is not None and "n_com_tracos" in C) else E
    col = "tracos_compl" if "tracos_compl" in K else None
    g = K.groupby("tamanho")
    cop = g.copresenca.mean()
    axD.plot(cop.index, cop.values, "o-", color=BLUE, lw=1.3, ms=4.5, zorder=4)
    axD.set_xlabel("Combination size (genera)")
    axD.set_ylabel("Co-presence (fraction of samples)", color=BLUE)
    axD.tick_params(axis="y", colors=BLUE)
    axD.set_ylim(0, None); axD.set_xticks(sorted(K.tamanho.unique()))
    if col:
        tr = g[col].mean()
        a2 = axD.twinx()
        a2.plot(tr.index, tr.values, "s--", color=ORANGE, lw=1.3, ms=4, zorder=4)
        a2.set_ylabel("Complementary traits (mean)", color=ORANGE)
        a2.tick_params(axis="y", colors=ORANGE)
        a2.spines["right"].set_visible(True); a2.spines["top"].set_visible(False)
        a2.set_ylim(0, None)
    grid(axD, "y"); painel(axD, "d")

    save(fig, out, "Fig4_corpo", T, dd)
    log(f"  Fig4_corpo: {T.study.nunique()} estudos nos trios")


# ------------------------------------------------------------- suplementares
def figS13(mm, out, dd):
    F = read(os.path.join(mm, "modulo", "frequencia.tsv"))
    if F is None:
        return log("  FigS13: falta modulo/frequencia.tsv")
    fig, ax = plt.subplots(figsize=(7.2, 2.8))
    H = F.head(16)
    cores = [BLUE if bool(v) else LGRAY for v in H.acessivel]
    ax.bar(range(len(H)), H.frac, color=cores, width=.74)
    ax.set_xticks(range(len(H)), H.genus, rotation=55, ha="right", fontsize=6.3, style="italic")
    ax.set_ylabel("Fraction of maximal sets")
    ax.set_ylim(0, 1)
    ax.legend(handles=[Patch(facecolor=BLUE, label="Accessible isolate"),
                       Patch(facecolor=LGRAY, label="No accessible isolate")],
              frameon=False, fontsize=6.4, loc="upper right")
    grid(ax, "y")
    fig.tight_layout()
    save(fig, out, "FigS13_empates", F, dd)


def figS14(mm, out, dd):
    R = read(os.path.join(mm, "remocao", "remocao.tsv"))
    K = read(os.path.join(mm, "consorcios", "ranking_aptos.tsv"))
    if R is None:
        return log("  FigS14: falta remocao/remocao.tsv")
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
    met = [(c, n) for c, n in (("frac_lcc", "Largest\ncomponent"), ("eficiencia", "Global\nefficiency"),
                               ("Q", "Modularity")) if f"{c}_delta" in R]
    ax[0].axhline(0, color=MUTED, lw=.6, zorder=2)
    rng = np.random.default_rng(2)
    for i, (c, _) in enumerate(met):
        v = (100 * R[f"{c}_delta"] / R[f"{c}_nulo"].replace(0, np.nan)).dropna().values
        ax[0].scatter(rng.normal(i + 1, .07, len(v)), v, s=5, color=BLUE, alpha=.55, lw=0, zorder=3)
    ax[0].set_xticks(range(1, len(met) + 1), [n for _, n in met], fontsize=6.5)
    ax[0].set_ylabel("Change after removal (%)")
    grid(ax[0], "y"); painel(ax[0], "a")

    if K is not None and "estrut_R2_Z" in K:
        ax[1].axhspan(-1.96, 1.96, color=LGRAY, alpha=.55, zorder=0)
        ax[1].axvspan(-1.96, 1.96, color=LGRAY, alpha=.55, zorder=0)
        ax[1].axhline(0, color=MUTED, lw=.6); ax[1].axvline(0, color=MUTED, lw=.6)
        ax[1].scatter(K.estrut_R2_Z, K.estrut_alfa_Z, s=7, color=VIOLET, alpha=.7, lw=0, zorder=3)
        ax[1].set_xlabel("Community composition (z vs. null)")
        ax[1].set_ylabel("Alpha diversity (z vs. null)")
    grid(ax[1], "both"); painel(ax[1], "b")
    fig.tight_layout()
    save(fig, out, "FigS14_remocao", R, dd)


def figS15(mm, out, dd, n_pares=12):
    Q = read(os.path.join(mm, "consorcios", "parceiros_elegiveis.tsv"))
    if Q is None or not len(Q):
        return log("  FigS15: falta parceiros_elegiveis.tsv")
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    q = Q.sort_values("copresenca_leg").tail(n_pares)
    y = np.arange(len(q))
    cor = [ORANGE if g == "rizobios" else BLUE for g in q.grupo]
    ax.hlines(y, 0, q.copresenca_leg, color=cor, lw=1.8, alpha=.8)
    ax.scatter(q.copresenca_leg, y, color=cor, s=30, zorder=4)
    ax.set_yticks(y, [f"{c} – {p}" for c, p in zip(q.rizobio, q.parceiro)],
                  fontsize=7, style="italic")
    ax.set_ylim(-0.8, len(q) - 0.2)
    ax.set_xlabel("Co-presence within the host group")
    ax.set_xlim(0, max(q.copresenca_leg) * 1.08)
    ax.grid(axis="x", color="#e8e6e1", lw=.6); ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(handles=[Line2D([], [], color=ORANGE, marker="o", ls="", label="Rhizobial centre"),
                       Line2D([], [], color=BLUE, marker="o", ls="", label="Other inoculant centre")],
              frameon=False, fontsize=6.5, ncol=2, loc="upper center",
              bbox_to_anchor=(0.5, -0.16))
    fig.tight_layout()
    save(fig, out, "FigS15_parceiros", q, dd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", default="Figuras")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    mm = os.path.join(root, "merged", "meta")
    out = os.path.abspath(os.path.expanduser(a.out)); dd = os.path.join(out, "dados")
    os.makedirs(dd, exist_ok=True)
    import traceback

    def roda(f, *args):
        try:
            f(*args)
        except Exception:
            log(f"  !! {f.__name__} falhou:"); traceback.print_exc()

    log("Corpo:")
    roda(fig3_corpo, root, mm, out, dd)
    roda(fig4_corpo, mm, out, dd)
    log("Suplementares:")
    roda(figS13, mm, out, dd)
    roda(figS14, mm, out, dd)
    roda(figS15, mm, out, dd)
    log(f"Figuras em {out}")


if __name__ == "__main__":
    main()
