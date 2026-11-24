# [S6.3] Figures: panel layout and crop labels
# Supplementary 6 of the manuscript. See README.md for the full map.
"""
Figuras do manuscrito, numeracao definitiva e sem colisao de nomes.

  CORPO
  Fig3_rede            rede de consenso; cor = papel do genero, nao modulo algoritmico
  Fig4_limites         a) remocao in silico  b) estruturacao  c) saturacao da copresenca
  Fig5_cultivabilidade a) isolados x clinicos  b) classe por grupo  c) funil de criterios
  Fig6_modulo          a) matriz de r  b) copresenca por tamanho  c) tracos funcionais

  SUPLEMENTAR (geradas aqui)
  FigS3_sequeiro_arroz, FigS4_vies_filos, FigS5_moderadores,
  FigS9_consorcios, FigS10_rizobios, FigS11_funcao

Fig1_dataset, Fig2_reprodutibilidade, FigS1, FigS2, FigS6, FigS7 e FigS8 vem de make_figures.py.
Todo painel de figura com mais de um painel recebe a), b), c) acima do eixo, a esquerda.

Uso: LD_PRELOAD=... python3 figures_new.py --root ~/syncom --out Figuras
"""
import argparse, itertools, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_figures import (BLUE, ORANGE, AQUA, VIOLET, RED, YELLOW, MAGENTA, GREEN,
                          GRAY, LGRAY, INK, MUTED, PGPR, log, grid)
from figures_extra import read, save, fig_vies_filos

MODULO = ["Arthrobacter", "Nocardioides", "Agromyces", "Streptomyces", "Kribbella", "Glycomyces"]
TRACO_EN = {"fixacao_N": "N fixation", "ACC_desaminase": "ACC deaminase",
            "PQQ_solub_fosfato": "P solubilization (PQQ)", "glicose_desidrogenase_gcd": "Glucose dehydrogenase",
            "fosfatases_fitase": "Phosphatase / phytase", "AIA_via_IPyA_ou_IAM": "IAA (IPyA or IAM)",
            "quitinase": "Chitinase", "HCN": "HCN", "acetoina_butanodiol": "Acetoin / butanediol",
            "sideroforo_NIS": "NIS siderophore"}
NIVEL_EN = {"seca": "Drought", "hidrico_ok": "Well watered", "N_mineral": "Mineral N",
            "controle": "Control", "organica": "Organic amendment"}
CULT_EN = {"milho": "Maize", "trigo": "Wheat", "arroz": "Rice", "sorgo": "Sorghum",
           "cana": "Sugarcane", "soja": "Soybean", "feijao/leguminosas": "Other legumes"}


def painel(ax, letra):
    ax.set_title(f"{letra})", loc="left", fontsize=9, fontweight="bold", color=INK, pad=6)


def en(c):
    return TRACO_EN.get(c, c.replace("BGC_", "").replace("_", " "))


def nivel(s):
    return NIVEL_EN.get(s, s.replace("_", " "))


# ------------------------------------------------------------- Fig3: a rede
def fig_rede(mm, root, out, dd, generos, min_abs_r, rotular=()):
    """Sub-rede dos generos que compoem os consorcios candidatos, colorida por
    disponibilidade de linhagem isolada."""
    cons = read(os.path.join(mm, "rede_entrada.tsv"))
    if cons is None:
        cons = read(os.path.join(mm, "consistentes.tsv"))
    occ = read(os.path.join(root, "merged", "genus_occupancy.tsv"))
    G0 = read(os.path.join(root, "genomas", "isolados", "por_genero.tsv"), index_col=0)
    if cons is None or occ is None or G0 is None:
        return log("  Fig3: faltam rede_entrada/genus_occupancy/por_genero")
    alvo = set(G0.index)
    e = cons[cons.genus_a.isin(alvo) & cons.genus_b.isin(alvo)]
    if e.empty:
        return log("  Fig3: nenhuma associacao entre generos candidatos")
    G = nx.Graph()
    for _, r in e.iterrows():
        G.add_edge(r.genus_a, r.genus_b, w=r.r_mu)
    sozinhos = len(alvo) - len(G)
    prev = occ.set_index("genus").ocupacao.to_dict()
    cls = G0.classe.to_dict()
    L = read(os.path.join(mm, "consorcios", "lista_pgp.tsv"))
    lista = dict(zip(L.genus, L.lista)) if L is not None and "lista" in (L.columns if L is not None else []) else {}
    def grupo(n):
        v = lista.get(n)
        if v == "classica":
            return "classica"
        if v or PGPR.match(n):
            return "relatada"
        return "outros"
    cat = {n: grupo(n) for n in G}
    COR = {"classica": ORANGE, "relatada": BLUE, "outros": LGRAY}
    ROT = {"classica": "Established plant-growth promoter",
           "relatada": "Reported plant-growth promoter",
           "outros": "No growth-promotion record"}
    ordem = ["outros", "relatada", "classica"]
    pos = nx.kamada_kawai_layout(G)
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.set_position([0.0, 0.0, 0.78, 1.0])
    for u, v, d in sorted(G.edges(data=True), key=lambda x: abs(x[2]["w"])):
        ax.plot(*zip(pos[u], pos[v]), color=BLUE if d["w"] > 0 else RED,
                lw=0.3 + 3.2 * abs(d["w"]), alpha=0.40, solid_capstyle="round", zorder=1)
    for c in ordem:
        for sem in (False, True):
            ns = [n for n in G if cat[n] == c and (cls.get(n) == "sem_cultivado") == sem]
            if not ns:
                continue
            ax.scatter([pos[n][0] for n in ns], [pos[n][1] for n in ns],
                       s=[10 + 110 * prev.get(n, 0.2) for n in ns], c=COR[c],
                       edgecolor=INK if sem else "white", lw=1.0 if sem else 0.6,
                       zorder=3 + ordem.index(c) + (1 if sem else 0))
    for n in rotular:
        if n in pos:
            ax.annotate(n, pos[n], fontsize=6.4, style="italic", color=INK,
                        xytext=(5, 4), textcoords="offset points", zorder=6)
    ax.axis("off"); ax.set_aspect("equal")
    h = [Line2D([], [], color=BLUE, lw=2, label="Positive association"),
         Line2D([], [], color=RED, lw=2, label="Negative association")]
    h += [Line2D([], [], marker="o", ls="", ms=7, mfc=COR[c], mec="white", label=ROT[c])
          for c in ordem[::-1] if any(cat[n] == c for n in G)]
    if any(cls.get(n) == "sem_cultivado" for n in G):
        h += [Line2D([], [], marker="o", ls="", ms=7, mfc=LGRAY, mec=INK, mew=1.0,
                     label="No cultured representative")]
    h += [Line2D([], [], marker="o", ls="", ms=np.sqrt(10 + 110 * v), mfc=LGRAY, mec="white",
                 label=f"Present in {v:.0%} of studies") for v in (0.25, 0.5, 1.0)]
    fig.legend(handles=h, loc="center left", bbox_to_anchor=(0.79, 0.5), frameon=False,
               fontsize=7, labelspacing=1.0)
    save(fig, out, "Fig3_rede")
    pd.DataFrame([{"genus": n, "pgp_group": cat[n], "class": cls.get(n), "degree": G.degree(n),
                   "degree_pos": sum(1 for _, _, d in G.edges(n, data=True) if d["w"] > 0),
                   "occupancy": prev.get(n, np.nan), "x": pos[n][0], "y": pos[n][1]} for n in G]) \
        .sort_values("degree", ascending=False) \
        .to_csv(os.path.join(dd, "Fig3_rede_nodes.tsv"), sep="\t", index=False)
    nc = sum(1 for n in G if cls.get(n) == "sem_cultivado")
    log(f"    grupos PGP: " + ", ".join(f"{c} {sum(1 for n in G if cat[n] == c)}" for c in ordem))
    log(f"    Fig3_rede: {len(G)} generos candidatos, {G.number_of_edges()} associacoes | "
        f"sem isolado: {nc} | sem nenhuma associacao interna (fora da figura): {sozinhos}")


# --------------------------------------------------------- Fig4: limites
def fig_limites(mm, out, dd):
    R = read(os.path.join(mm, "remocao", "remocao.tsv"))
    K = read(os.path.join(mm, "consorcios", "ranking_aptos.tsv"))
    S = read(os.path.join(mm, "consorcios", "estrelas_saturacao.tsv"))
    E = read(os.path.join(mm, "consorcios", "estrelas.tsv"))
    if R is None or S is None or E is None:
        return log("  Fig4: faltam remocao.tsv/estrelas*.tsv")
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.0))
    ax = ax.ravel()
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

    x = S.tamanho.values
    b = ax[2].boxplot([E.copresenca[E.tamanho == t].values for t in x], positions=x, widths=.55,
                      patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1))
    for p in b["boxes"]:
        p.set(facecolor=LGRAY, edgecolor=MUTED, lw=.6)
    ax[2].plot(x, S.copres_melhor, "o", color=RED, ms=4, zorder=4)
    ax[2].set_xticks(x); ax[2].set_xlabel("Consortium size (genera)")
    ax[2].set_ylabel("Co-presence (fraction of samples)")
    ax[2].legend(handles=[Line2D([], [], marker="o", ls="", color=RED, label="Best consortium"),
                          Patch(facecolor=LGRAY, edgecolor=MUTED, label="All consortia")],
                 frameon=False, fontsize=6.5, loc="upper right")
    grid(ax[2], "y"); painel(ax[2], "c")

    # painel d: compromisso entre viabilidade observacional e ganho funcional esperado
    C = read(os.path.join(mm, "consorcios", "complementaridade.tsv"))
    if C is not None and "n_com_tracos" in C:
        C = C[C.n_com_tracos == C.tamanho]
    K = C if C is not None and len(C) else K
    if K is not None and {"tamanho", "copresenca"} <= set(K.columns):
        col = "tracos_compl" if "tracos_compl" in K else ("tracos_PGP" if "tracos_PGP" in K else None)
        g = K.groupby("tamanho")
        cop = g.copresenca.mean()
        ax[3].plot(cop.index, cop.values, "o-", color=BLUE, lw=1.3, ms=4.5, zorder=4)
        ax[3].set_xlabel("Consortium size (genera)")
        ax[3].set_ylabel("Co-presence (fraction of samples)", color=BLUE)
        ax[3].tick_params(axis="y", colors=BLUE)
        ax[3].set_ylim(0, None); ax[3].set_xticks(sorted(K.tamanho.unique()))
        if col:
            tr = g[col].mean()
            a2 = ax[3].twinx()
            a2.plot(tr.index, tr.values, "s--", color=ORANGE, lw=1.3, ms=4, zorder=4)
            a2.set_ylabel("Complementary traits (mean)", color=ORANGE)
            a2.tick_params(axis="y", colors=ORANGE)
            a2.spines["right"].set_visible(True); a2.spines["top"].set_visible(False)
            a2.set_ylim(0, None)
        grid(ax[3], "y")
    painel(ax[3], "d")
    fig.tight_layout()
    save(fig, out, "Fig4_limites", R, dd)


# ------------------------------------------------- Fig5: cultivabilidade
def fig_cultivabilidade(root, mm, out, dd):
    G = read(os.path.join(root, "genomas", "isolados", "por_genero.tsv"), index_col=0)
    if G is None:
        return log("  Fig5: falta genomas/isolados/por_genero.tsv")
    grupo = pd.Series("Other candidate genera", index=G.index)
    npar = {}
    for sub, nome in (("rizobios", "Rhizobial partners"), ("pgp", "PGP partners")):
        P = read(os.path.join(mm, sub, "parceiros.tsv"))
        if P is None:
            continue
        p = P[P.ok_leg] if "ok_leg" in P else P
        npar[nome] = int(len(p))
        grupo.loc[grupo.index.intersection(set(p.parceiro))] = nome
    COR = {"Rhizobial partners": ORANGE, "PGP partners": BLUE, "Other candidate genera": GRAY}
    G = G.assign(grupo=grupo)
    iso = G.genomas_isolados.astype(float)
    G["frac_clinica"] = np.where(iso > 0, G.isolados_clinicos.astype(float) / iso.clip(lower=1), np.nan)

    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1.25, 1, .85]})
    for g, c in COR.items():
        d = G[G.grupo == g]
        ax[0].scatter(d.frac_isolados, d.frac_clinica,
                      s=6 + 14 * np.log10(d.genomas.clip(lower=1) + 1),
                      color=c, alpha=.75, lw=.3, edgecolor="white", zorder=3)
    ax[0].axvline(.5, color=MUTED, lw=.6, ls="--"); ax[0].axhline(.25, color=MUTED, lw=.6, ls="--")
    ax[0].set_xlabel("Genomes from isolates (fraction)")
    ax[0].set_ylabel("Isolates of clinical origin (fraction)")
    ax[0].set_xlim(-0.03, 1.03); ax[0].set_ylim(-0.03, None)
    ax[0].legend(handles=[Line2D([], [], marker="o", ls="", ms=5, mfc=c, mec="white", label=g)
                          for g, c in COR.items()], frameon=False, fontsize=6.3,
                 loc="upper center", bbox_to_anchor=(.5, -.26), ncol=1)
    grid(ax[0], "both"); painel(ax[0], "a")

    ordem = ["acessivel", "dificil", "cultivado_sem_genoma", "conferir_fatia", "sem_cultivado"]
    rot = {"acessivel": "Accessible", "dificil": "Few isolates",
           "cultivado_sem_genoma": "Cultured, no isolate genome",
           "conferir_fatia": "Check circumscription", "sem_cultivado": "No cultured representative"}
    cores = [BLUE, YELLOW, AQUA, VIOLET, RED]
    grupos = ["PGP partners", "Rhizobial partners", "Other candidate genera"]
    T = (G.groupby(["grupo", "classe"]).size().unstack(fill_value=0)
         .reindex(index=grupos, columns=ordem, fill_value=0))
    frac = T.div(T.sum(axis=1).clip(lower=1), axis=0)
    y = np.arange(len(grupos))[::-1]
    left = np.zeros(len(grupos))
    for cl, c in zip(ordem, cores):
        ax[1].barh(y, frac[cl].values, left=left, color=c, height=.6, edgecolor="white", lw=.5)
        left += frac[cl].values
    ax[1].set_yticks(y, [f"{g}\n(n = {int(T.loc[g].sum())})" for g in grupos], fontsize=6.4)
    ax[1].set_xlabel("Genera (fraction)"); ax[1].set_xlim(0, 1)
    ax[1].tick_params(axis="y", length=0); ax[1].spines["left"].set_visible(False)
    ax[1].legend(handles=[Patch(facecolor=c, label=rot[cl]) for cl, c in zip(ordem, cores)],
                 frameon=False, fontsize=6.2, ncol=2, loc="upper center", bbox_to_anchor=(.5, -.26))
    painel(ax[1], "b")

    # funil: do par positivo ao par testavel
    P = read(os.path.join(mm, "pgp", "parceiros.tsv"))
    etapas, vals = [], []
    if P is not None:
        ok = P[P.ok_leg] if "ok_leg" in P else P
        cult = G.classe.isin(["acessivel", "cultivado_sem_genoma"])
        acess = set(G.index[cult])
        limpo = set(G.index[cult & (G.frac_clinica.fillna(0) < .25)])
        etapas = ["All", "Confirmed", "Culturable", "Low clinical"]
        vals = [int(len(P)), int(len(ok)), int(ok.parceiro.isin(acess).sum()),
                int(ok.parceiro.isin(limpo).sum())]
        ax[2].bar(range(len(vals)), vals, color=[GRAY, GRAY, BLUE, GREEN], width=.62)
        for i, v in enumerate(vals):
            ax[2].text(i, v, str(v), ha="center", va="bottom", fontsize=6.4, color=INK)
        ax[2].set_xticks(range(len(vals)), etapas, fontsize=6.0, rotation=30, ha="right")
        ax[2].set_ylabel("Centre-partner pairs")
        grid(ax[2], "y")
    painel(ax[2], "c")
    fig.tight_layout()
    save(fig, out, "Fig5_cultivabilidade", G.reset_index(), dd)


# ------------------------------------------------------------- Fig6: modulo
def copresenca_modulo(root, generos):
    try:
        from meta_pairs import load
        from consortia import load_presence, copresence
    except Exception as e:
        log(f"  (copresenca do modulo nao calculada: {e})"); return None
    md = os.path.join(root, "merged")
    d, qc = load(root, False)
    pres = load_presence(md, set(generos), qc)
    rows = []
    for n in range(2, len(generos) + 1):
        melhor = None
        for c in itertools.combinations(generos, n):
            cp, nst, _ = copresence(list(c), pres)
            if melhor is None or cp > melhor[0]:
                melhor = (cp, nst, " + ".join(c))
        rows.append({"tamanho": n, "copresenca": melhor[0], "estudos": melhor[1], "membros": melhor[2]})
    return pd.DataFrame(rows)


def fig_modulo(root, mm, out, dd, generos):
    """Fig6: substituibilidade do conjunto denso, copresenca dos tres conjuntos e
    matriz do conjunto executavel (linhagens acessiveis)."""
    md = os.path.join(mm, "modulo")
    F = read(os.path.join(md, "frequencia.tsv"))
    P = read(os.path.join(md, "copresenca.tsv"))
    Q = read(os.path.join(md, "pares.tsv"))
    if F is None or P is None or Q is None:
        return log("  Fig6: falta merged/meta/modulo/ (rode modulo_fig.py)")
    cmap = LinearSegmentedColormap.from_list("seq", ["#f4f3f0", BLUE])
    fig = plt.figure(figsize=(7.2, 6.2))

    # a) frequencia de cada genero entre os cliques maximos
    axA = fig.add_axes([0.11, 0.70, 0.82, 0.25])
    H = F.head(16)
    cores = [BLUE if bool(v) else LGRAY for v in H.acessivel]
    axA.bar(range(len(H)), H.frac, color=cores, width=.74)
    axA.set_xticks(range(len(H)), H.genus, rotation=55, ha="right", fontsize=6.3, style="italic")
    axA.set_ylabel("Fraction of maximal cliques")
    axA.set_ylim(0, 1)
    axA.legend(handles=[Patch(facecolor=BLUE, label="Accessible isolate"),
                        Patch(facecolor=LGRAY, label="No accessible isolate")],
               frameon=False, fontsize=6.4, loc="upper right")
    axA.grid(axis="y", color="#e8e6e1", lw=.6); axA.set_axisbelow(True)
    painel(axA, "a")

    # b) copresenca por tamanho, tres conjuntos
    axB = fig.add_axes([0.11, 0.07, 0.36, 0.35])
    estilo = {"dez": (GRAY, "o--", "Densest clique"),
              "replicado": (ORANGE, "o-", "Best replicated"),
              "acessivel": (BLUE, "o-", "Replicated and accessible")}
    hs = []
    for nome, (cor, st, rot) in estilo.items():
        d = P[P.conjunto == nome]
        if d.empty:
            continue
        axB.plot(d.tamanho, d.copres_media, st, color=cor, lw=1.3, ms=4.2)
        hs.append(Line2D([], [], color=cor, marker="o", ms=4,
                         ls="--" if "--" in st else "-", label=rot))
    axB.set_xlabel("Consortium size (genera)")
    axB.set_ylabel("Co-presence (fraction of samples)")
    axB.set_ylim(0, None)
    axB.legend(handles=hs, frameon=False, fontsize=6.4, loc="upper right")
    axB.grid(axis="y", color="#e8e6e1", lw=.6); axB.set_axisbelow(True)
    painel(axB, "b")

    # c) matriz do conjunto executavel
    A = Q[Q.conjunto == "acessivel"]
    g = sorted(set(A.genus_a) | set(A.genus_b))
    n = len(g)
    Z = np.full((n, n), np.nan)
    for _, r in A.iterrows():
        i1, j1 = g.index(r.genus_a), g.index(r.genus_b)
        Z[i1, j1] = Z[j1, i1] = r.r_mu
    axC = fig.add_axes([0.60, 0.07, 0.27, 0.35])
    im = axC.imshow(np.ma.masked_invalid(Z), cmap=cmap, vmin=0,
                    vmax=max(0.35, np.nanmax(Z) if n else 0.35))
    axC.set_xticks(range(n), g, rotation=55, ha="right", fontsize=6.3, style="italic")
    axC.set_yticks(range(n), g, fontsize=6.3, style="italic")
    axC.set_xticks(np.arange(-.5, n, 1), minor=True)
    axC.set_yticks(np.arange(-.5, n, 1), minor=True)
    axC.grid(which="minor", color="white", lw=.8); axC.tick_params(which="minor", length=0)
    cb = fig.colorbar(im, ax=axC, fraction=.045, pad=.04)
    cb.set_label("Pooled effect (r)", fontsize=6.6); cb.ax.tick_params(labelsize=6)
    painel(axC, "c")

    save(fig, out, "Fig6_modulo", P, dd)
    log(f"  Fig6_modulo: {len(F)} generos entre os empates | conjunto executavel: "
        f"{' + '.join(g)}")

def figS_consorcios(mm, out, dd):
    S = read(os.path.join(mm, "consorcios", "estrelas_saturacao.tsv"))
    E = read(os.path.join(mm, "consorcios", "estrelas.tsv"))
    if S is None or E is None:
        return log("  FigS9: faltam estrelas*.tsv")
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4))
    x = S.tamanho.values
    ax[0].bar(x, S.n.values, color=BLUE, width=.62)
    ax[0].set_yscale("log"); ax[0].set_xlabel("Consortium size (genera)")
    ax[0].set_ylabel("Candidate consortia"); painel(ax[0], "a")
    ax[1].plot(x, S.rcond_min_melhor, "o-", color=ORANGE, lw=1.2, ms=4)
    ax[1].set_xlabel("Consortium size (genera)"); ax[1].set_ylabel("Best weakest link (r)")
    ax[1].set_ylim(0, None); painel(ax[1], "b")
    b = ax[2].boxplot([E.copresenca[E.tamanho == t].values for t in x], positions=x, widths=.55,
                      patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1))
    for p in b["boxes"]:
        p.set(facecolor=LGRAY, edgecolor=MUTED, lw=.6)
    ax[2].plot(x, S.copres_melhor, "o", color=RED, ms=4)
    ax[2].set_xlabel("Consortium size (genera)")
    ax[2].set_ylabel("Co-presence (fraction of samples)"); painel(ax[2], "c")
    for a in ax:
        a.set_xticks(x); grid(a, "y")
    fig.tight_layout()
    save(fig, out, "FigS9_consorcios", S, dd)


def figS_rizobios(mm, out, dd, n=26):
    P = read(os.path.join(mm, "rizobios", "parceiros.tsv"))
    if P is None:
        return log("  FigS10: falta rizobios/parceiros.tsv")
    d = P[P.ok_leg & P.r_leg.notna()].copy()
    if d.empty:
        return log("  FigS10: nenhum parceiro confirmado")
    d["rotulo"] = d.rizobio + " – " + d.parceiro
    d = d.nlargest(min(n, len(d)), "r_leg").iloc[::-1]
    hubs = sorted(d.rizobio.unique())
    cor = dict(zip(hubs, [BLUE, ORANGE, AQUA, VIOLET, GREEN, YELLOW, MAGENTA, GRAY, RED][:len(hubs)]))
    fig, ax = plt.subplots(1, 2, figsize=(7.2, max(0.19 * len(d) + 1.3, 3.2)),
                           gridspec_kw={"width_ratios": [2.1, 1]}, sharey=True)
    y = np.arange(len(d))
    ax[0].barh(y, d.r_leg, color=[cor[h] for h in d.rizobio], height=.66)
    m = d.r_gram.notna().values
    ax[0].plot(d.r_gram.values[m], y[m], "o", ms=3.4, mfc="white", mec=INK, mew=.8)
    ax[0].set_yticks(y); ax[0].set_yticklabels(d.rotulo, fontsize=6.4, style="italic")
    ax[0].set_xlabel("Effect size in legumes (r)"); ax[0].axvline(0, color=MUTED, lw=.6)
    ax[0].legend(handles=[Patch(facecolor=cor[h], label=h) for h in hubs] +
                         [Line2D([], [], marker="o", ls="", mfc="white", mec=INK,
                                 label="Same pair in grasses")],
                 frameon=False, fontsize=6.2, loc="lower right")
    painel(ax[0], "a")
    ax[1].barh(y, d.copresenca_leg, color=LGRAY, edgecolor=MUTED, lw=.5, height=.66)
    ax[1].set_xlabel("Co-presence in legumes (fraction of samples)"); painel(ax[1], "b")
    fig.tight_layout()
    save(fig, out, "FigS10_rizobios", d, dd)


def figS_funcao(root, mm, out, dd):
    T = read(os.path.join(root, "genomas", "pgp", "tracos_por_genero.tsv"), index_col=0)
    B = read(os.path.join(root, "genomas", "bgc", "bgc_por_genero.tsv"), index_col=0)
    R = read(os.path.join(mm, "consorcios", "ranking_aptos.tsv"), index_col=0)
    if T is None or R is None:
        return log("  FigS11: faltam tracos/ranking")
    gen = []
    for m in R.membros.head(12):
        gen += [x.strip() for x in str(m).split("+")]
    gen = [g for g in dict.fromkeys(gen) if g in T.index]
    if not gen:
        return log("  FigS11: nenhum genero com tracos")
    M = T.loc[gen, [c for c in T.columns if c != "n_genomas"]]
    if B is not None:
        tipos = [c for c in B.columns if c.startswith("BGC_") and c not in ("BGC_total", "BGC_Unknown")]
        M = M.join(B.loc[B.index.intersection(gen), tipos], how="left")
    M = M.fillna(0)
    cmap = LinearSegmentedColormap.from_list("seq", ["#f4f3f0", BLUE])
    fig, ax = plt.subplots(figsize=(max(0.34 * len(M.columns) + 2.2, 5.6), 0.26 * len(M) + 1.5))
    im = ax.imshow(M.values, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(M.columns)), [en(c) for c in M.columns], rotation=55, ha="right", fontsize=6.4)
    ax.set_yticks(range(len(M)), list(M.index), fontsize=6.6, style="italic")
    ax.set_xticks(np.arange(-.5, len(M.columns), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(M), 1), minor=True)
    ax.grid(which="minor", color="white", lw=.8); ax.tick_params(which="minor", length=0)
    cb = fig.colorbar(im, ax=ax, fraction=.022, pad=.02)
    cb.set_label("Fraction of genomes", fontsize=6.6); cb.ax.tick_params(labelsize=6)
    fig.tight_layout()
    save(fig, out, "FigS11_funcao", M.reset_index(), dd)


def figS_sequeiro_arroz(mm, out, dd):
    A = read(os.path.join(mm, "arroz", "meta_pairs_arroz.tsv.gz"))
    M = read(os.path.join(mm, "meta_pairs.tsv.gz"))
    if A is None or M is None:
        return log("  FigS3: faltam meta_pairs_arroz/meta_pairs")
    key = lambda d: d.genus_a + "|" + d.genus_b
    A.index = key(A); M.index = key(M)
    both = A.index.intersection(M.index)
    d = pd.DataFrame({"r_total": M.loc[both, "r_mu"], "r_arroz": A.loc[both, "r_mu"],
                      "sig_total": M.loc[both, "fdr"] < 0.05, "sig_arroz": A.loc[both, "fdr"] < 0.05})
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1.15, 1]})
    ax[0].axhline(0, color=LGRAY, lw=.6); ax[0].axvline(0, color=LGRAY, lw=.6)
    ax[0].plot([-.6, .8], [-.6, .8], color=GRAY, lw=.6, ls="--")
    s = d.sig_total & d.sig_arroz
    ax[0].scatter(d.r_total[~s], d.r_arroz[~s], s=3, color=LGRAY, lw=0)
    ax[0].scatter(d.r_total[s], d.r_arroz[s], s=5, color=BLUE, lw=0)
    ax[0].set_xlabel("Effect size, all crops (r)"); ax[0].set_ylabel("Effect size, paddy rice (r)")
    ax[0].legend(handles=[Line2D([], [], marker="o", ls="", color=BLUE, ms=3.6, label="Significant in both"),
                          Line2D([], [], marker="o", ls="", color=LGRAY, ms=3.6, label="Other pairs")],
                 frameon=False, fontsize=6.5, loc="upper left")
    painel(ax[0], "a")
    v = [int(s.sum()), int((d.sig_total & ~d.sig_arroz).sum()), int((~d.sig_total & d.sig_arroz).sum())]
    ax[1].bar(["Both", "All crops only", "Rice only"], v, color=[BLUE, GRAY, AQUA], width=.6)
    ax[1].set_ylabel("Significant pairs"); ax[1].tick_params(axis="x", labelsize=6.8)
    grid(ax[1], "y"); painel(ax[1], "b")
    fig.tight_layout()
    save(fig, out, "FigS3_sequeiro_arroz", d.reset_index(names="pair"), dd)


def figS_moderadores(mm, out, dd):
    P = read(os.path.join(mm, "moderadores", "por_par.tsv.gz"))
    if P is None:
        return log("  FigS5: falta moderadores/por_par.tsv.gz")
    cs = list(P.contraste.unique())
    fig, ax = plt.subplots(1, len(cs), figsize=(3.4 * len(cs), 3.1), squeeze=False)
    for i, (a, c) in enumerate(zip(ax[0], cs)):
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
        a.set_xlabel(f"Effect size, {nivel(n2)} (r)"); a.set_ylabel(f"Effect size, {nivel(n1)} (r)")
        a.set_xlim(-lim, lim); a.set_ylim(-lim, lim); a.set_aspect("equal")
        if sig.any():
            a.legend(handles=[Line2D([], [], marker="o", ls="", color=RED, ms=3.6,
                                     label="Significant moderation")], frameon=False,
                     fontsize=6.5, loc="upper left")
        painel(a, "abcdefg"[i])
    fig.tight_layout()
    save(fig, out, "FigS5_moderadores", P, dd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", default="Figuras")
    ap.add_argument("--modulo", default=",".join(MODULO))
    ap.add_argument("--min-abs-r", type=float, default=0.15)
    ap.add_argument("--rotular", default="", help="generos a rotular na Fig3")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    mm = os.path.join(root, "merged", "meta")
    out = os.path.abspath(os.path.expanduser(a.out)); dd = os.path.join(out, "dados")
    os.makedirs(dd, exist_ok=True)
    gen = [x.strip() for x in a.modulo.split(",") if x.strip()]
    import traceback

    def roda(f, *args):
        try:
            f(*args)
        except Exception:
            log(f"  !! {f.__name__} falhou:"); traceback.print_exc()

    log("Corpo:")
    roda(fig_rede, mm, root, out, dd, gen, a.min_abs_r, tuple(x.strip() for x in a.rotular.split(',') if x.strip()))
    roda(fig_limites, mm, out, dd)
    roda(fig_cultivabilidade, root, mm, out, dd)
    roda(fig_modulo, root, mm, out, dd, gen)
    log("Suplementares:")
    roda(figS_sequeiro_arroz, mm, out, dd)
    roda(fig_vies_filos, mm, out, dd)
    roda(figS_moderadores, mm, out, dd)
    roda(figS_consorcios, mm, out, dd)
    roda(figS_rizobios, mm, out, dd)
    roda(figS_funcao, root, mm, out, dd)
    log(f"Figuras em {out}")


if __name__ == "__main__":
    main()
