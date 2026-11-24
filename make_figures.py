#!/usr/bin/env python3
# [S6.2] Figures: shared palette and plotting helpers
# Supplementary 6 of the manuscript. See README.md for the full map.
"""
Figuras do manuscrito a partir dos resultados reais.
Convencao: sem titulos, legendas ou anotacoes dentro dos paineis (so eixos e rotulos de
categoria); tudo e explicado na legenda do manuscrito. Cada figura sai em PNG (300 dpi) e
PDF (vetorial), e os dados de cada figura vao para <saida>/dados/.

  Fig1_acervo        mapa dos estudos + composicao (cultura, tratamento, regiao/modo)
  Fig2_rede          rede de consenso (pares consistentes), modulos por cor
  Fig3_heatmap       efeito por estudo dos pares de destaque + efeito combinado
  FigS1_nucleo       ocupacao x abundancia dos generos
  FigS2_nulo         dados reais x modelo nulo (efeitos combinados)

Uso:
  LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 python3 make_figures.py \
      --root ~/syncom --studies subset_studies.tsv --world world_110m.geojson --out Figuras
"""
import argparse, glob, gzip, json, os, re, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, TwoSlopeNorm
import networkx as nx

log = lambda *x: print(*x, file=sys.stderr, flush=True)
BLUE, ORANGE, AQUA, VIOLET, RED = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e34948"
YELLOW, MAGENTA, GREEN = "#eda100", "#e87ba4", "#008300"
GRAY, LGRAY, INK, MUTED = "#9a9994", "#d9d8d3", "#0b0b0b", "#52514e"
GRUPO_COR = {"gramínea": BLUE, "leguminosa": ORANGE}
MOD_COR = [AQUA, VIOLET, YELLOW, MAGENTA, GREEN]          # modulos (azul/vermelho = arestas)
mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "savefig.dpi": 300, "savefig.bbox": "tight",
    "pdf.fonttype": 42})
PGPR = re.compile(r"^(Pseudomonas|Bacillus|Paenibacillus|Streptomyces|Bradyrhizobium|Rhizobium|"
                  r"Mesorhizobium|Sinorhizobium|Ensifer|Azospirillum|Herbaspirillum|Burkholderia|"
                  r"Paraburkholderia|Enterobacter|Pantoea|Serratia|Klebsiella|Azotobacter|"
                  r"Gluconacetobacter|Arthrobacter|Microbacterium|Variovorax|Brucella|Stenotrophomonas|"
                  r"Lysobacter|Flavobacterium|Achromobacter|Delftia|Acinetobacter|Kosakonia|"
                  r"Priestia|Peribacillus|Lysinibacillus|Rhodococcus|Agrobacterium)(_[A-Z]+)?$")
CROP_EN = {"milho": "Maize", "trigo": "Wheat", "arroz": "Rice", "sorgo": "Sorghum", "cana": "Sugarcane",
           "soja": "Soybean", "feijao/leguminosas": "Other legumes"}
GROUP_EN = {"gramínea": "Grasses", "leguminosa": "Legumes"}
TRT_EN = {"nenhum": "None detected", "aduba_mineral": "Mineral fertilization", "seca_agua": "Water/drought",
          "genotipo_planta": "Plant genotype", "aduba_organica": "Organic amendment",
          "patogeno": "Pathogen", "rotacao_cultivo": "Crop rotation", "contaminacao": "Contamination",
          "manejo_org_conv": "Organic vs. conventional", "salinidade": "Salinity",
          "temperatura_co2": "Temperature/CO₂", "preparo_solo": "Tillage"}
CROP_ORDER = ["milho", "trigo", "arroz", "sorgo", "cana", "soja", "feijao/leguminosas"]
# culturas: sem azul/laranja (reservados a gramineas/leguminosas) e sem vermelho (negativo)
CROP_COR = [AQUA, YELLOW, VIOLET, GREEN, MUTED, MAGENTA, "#b5b3ad"]
COUNTRY_ALIAS = {"USA": "United States of America", "United States": "United States of America",
                 "Korea": "South Korea", "Republic of Korea": "South Korea",
                 "Russian Federation": "Russia", "Czech Republic": "Czechia",
                 "Viet Nam": "Vietnam", "UK": "United Kingdom"}


LINEAGE = {}


def carrega_linhagem(md):
    f = os.path.join(md, "genus_lineage.tsv")
    if not os.path.exists(f):
        log("  (sem genus_lineage.tsv: rotulos ficam so com o codigo do genero)"); return
    d = pd.read_csv(f, sep="\t").fillna("")
    for _, r in d.iterrows():
        LINEAGE[r.genus] = (r.family, r["order"], r["class"], r.phylum)


def formal(n):
    """Nome taxonomico publicado: sem digitos, sem hifen, em caixa de nome proprio."""
    n = str(n)
    return bool(n) and not any(c.isdigit() for c in n) and "-" not in n \
        and n[:1].isupper() and n[1:2].islower()


def nome_bonito(g):
    """Genero com designacao provisoria do GTDB recebe como prefixo o menor nivel
    taxonomico com nome formal; se nenhum nivel tiver, fica so o codigo."""
    if formal(g):
        return g
    for t in LINEAGE.get(g, ()):
        if formal(t) and t != g:
            return f"{t} {g}"
    return g


def grid(ax, axis="y"):
    ax.grid(axis=axis, color=LGRAY, lw=0.5, zorder=0); ax.set_axisbelow(True)


def save(fig, out, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"))
    plt.close(fig)
    log(f"  {name}")


def draw_world(ax, world):
    for f in world["features"]:
        g = f["geometry"]
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        for p in polys:
            xy = np.array(p[0])
            ax.fill(xy[:, 0], xy[:, 1], color="#ecebe7", ec="white", lw=0.3, zorder=1)


# --------------------------------------------------------------------- Fig. 1
def fig1(qc, studies, world, out, dd):
    cent = {f["properties"]["name"]: (f["properties"]["lon"], f["properties"]["lat"])
            for f in world["features"]}
    rng = np.random.default_rng(1)
    q = qc[qc.status == "incluido"].copy()
    q["pais"] = q.country.str.split(":").str[0].str.strip().replace(COUNTRY_ALIAS)
    q["lon"] = q.pais.map(lambda c: cent.get(c, (np.nan, np.nan))[0])
    q["lat"] = q.pais.map(lambda c: cent.get(c, (np.nan, np.nan))[1])
    q["lon"] += rng.normal(0, 2.5, len(q)); q["lat"] += rng.normal(0, 2.0, len(q))
    q["amostras_ok"] = q.amostras_ok.astype(int)
    sem = q[q.lon.isna()].pais.unique()
    if len(sem):
        log(f"  (paises sem coordenada no mapa: {', '.join(map(str, sem))})")
    fig = plt.figure(figsize=(7.2, 5.4))
    ax = fig.add_axes([0, 0.42, 1, 0.58]); draw_world(ax, world)
    for g, c in GRUPO_COR.items():
        m = q.grupo == g
        ax.scatter(q.lon[m], q.lat[m], s=q.amostras_ok[m] * 1.1, c=c, alpha=0.75,
                   edgecolor="white", lw=0.5, zorder=3)
    ax.set_xlim(-170, 180); ax.set_ylim(-58, 80); ax.axis("off")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", ms=7, mfc=c, mec="white", label=GROUP_EN[g])
         for g, c in GRUPO_COR.items()]
    h += [Line2D([], [], marker="o", ls="", ms=np.sqrt(v * 1.1), mfc=GRAY, mec="white",
                 label=f"{v} samples") for v in (25, 50, 100)]
    ax.legend(handles=h, loc="lower left", frameon=False, fontsize=7, labelspacing=0.9,
              bbox_to_anchor=(0.0, 0.02))
    # composicao
    trat = studies.set_index("study_accession").reindex(q.study).tipos_tratamento.fillna("")
    tc = pd.Series([t for x in trat for t in (x.split(";") if x else ["nenhum"])]).value_counts()
    tc = tc.rename(index=TRT_EN)
    cult = q.cultura.map(CROP_EN).value_counts()
    modo = (q.region.str.replace("V4|V3-V4?", "V4", regex=False).str.replace("V3-V4", "V3–V4")
            + q.modo.map({"PE": "", "R1": " (R1 only)"})).value_counts()
    blocos = [(cult, [ORANGE if c in ("Soybean", "Other legumes") else BLUE for c in cult.index]),
              (tc, [GRAY if i == "None detected" else BLUE for i in tc.index]),
              (modo, [BLUE] * len(modo))]
    for i, (s, cols) in enumerate(blocos):
        ax = fig.add_axes([[0.13, 0.52, 0.85][i], 0.04, 0.15, 0.32])
        y = np.arange(len(s))[::-1]
        ax.barh(y, s.values, color=cols, height=0.7, zorder=3)
        ax.set_yticks(y, list(s.index)); ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False); grid(ax, "x"); ax.set_xlabel("Studies")
    save(fig, out, "Fig1_dataset")
    q[["study", "pais", "cultura", "grupo", "region", "modo", "amostras_ok"]].to_csv(
        os.path.join(dd, "Fig1_estudos.tsv"), sep="\t", index=False)


# --------------------------------------------------------------------- Fig. 2
def fig2(cons, occ, out, dd, min_abs_r, name="Fig2_network", only=None):
    from matplotlib.lines import Line2D
    e = cons[np.abs(cons.r_mu) >= min_abs_r].copy()
    if only is not None:                     # subrede: arestas com pelo menos um genero da lista
        e = e[e.genus_a.isin(only) | e.genus_b.isin(only)]
    if e.empty:
        log(f"  ({name} pulada: sem arestas)"); return
    G = nx.Graph()
    for _, r in e.iterrows():
        G.add_edge(r.genus_a, r.genus_b, w=r.r_mu, k=r.k)
    Gp = nx.Graph([(u, v, d) for u, v, d in G.edges(data=True) if d["w"] > 0])
    Gp.add_nodes_from(G)
    mods = sorted(nx.community.greedy_modularity_communities(Gp, weight="w"), key=len, reverse=True)
    mod_of = {n: i for i, m in enumerate(mods) for n in m}
    prev = occ.set_index("genus").ocupacao.to_dict()
    pos = nx.kamada_kawai_layout(G) if len(G) <= 400 else nx.spring_layout(G, seed=3, k=0.15)
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    ax.set_position([0.0, 0.0, 0.78, 1.0])
    for u, v, d in sorted(G.edges(data=True), key=lambda x: abs(x[2]["w"])):
        ax.plot(*zip(pos[u], pos[v]), color=BLUE if d["w"] > 0 else RED,
                lw=0.3 + 3.5 * abs(d["w"]), alpha=0.45, solid_capstyle="round", zorder=1)
    nmods = [i for i, m in enumerate(mods) if len(m) >= 5][:len(MOD_COR)]
    for n in G:
        m = mod_of.get(n, 99)
        c = MOD_COR[nmods.index(m)] if m in nmods else GRAY
        ring = only is not None and n in only
        ax.scatter(*pos[n], s=8 + 90 * prev.get(n, 0.2), c=c, edgecolor=INK if ring else "white",
                   lw=1.2 if ring else 0.6, zorder=3)
    ax.axis("off"); ax.set_aspect("equal")
    h = [Line2D([], [], color=BLUE, lw=2, label="Positive association"),
         Line2D([], [], color=RED, lw=2, label="Negative association")]
    h += [Line2D([], [], marker="o", ls="", ms=7, mfc=MOD_COR[j], mec="white",
                 label=f"Module {j + 1} ({len(mods[m])} genera)") for j, m in enumerate(nmods)]
    h += [Line2D([], [], marker="o", ls="", ms=7, mfc=GRAY, mec="white", label="Other genera")]
    if only is not None:
        h += [Line2D([], [], marker="o", ls="", ms=7, mfc="white", mec=INK, mew=1.2, label="PGP genus")]
    h += [Line2D([], [], marker="o", ls="", ms=np.sqrt(8 + 90 * v), mfc=LGRAY, mec="white",
                 label=f"Occupancy {v:.0%}") for v in (0.25, 0.5, 1.0)]
    fig.legend(handles=h, loc="center left", bbox_to_anchor=(0.79, 0.5), frameon=False,
               fontsize=7, labelspacing=1.0)
    save(fig, out, name)
    pd.DataFrame([{"genus": n, "module": (nmods.index(mod_of[n]) + 1) if mod_of.get(n) in nmods else 0,
                   "degree": G.degree(n),
                   "degree_pos": sum(1 for _, _, d in G.edges(n, data=True) if d["w"] > 0),
                   "occupancy": prev.get(n, np.nan), "pgp_genus": bool(PGPR.match(n)),
                   "x": pos[n][0], "y": pos[n][1]} for n in G]) \
        .sort_values(["module", "degree"], ascending=[True, False]) \
        .to_csv(os.path.join(dd, f"{name}_nodes.tsv"), sep="\t", index=False)
    e.to_csv(os.path.join(dd, f"{name}_edges.tsv"), sep="\t", index=False)
    log(f"    {name}: {len(G)} genera, {G.number_of_edges()} edges (|r| >= {min_abs_r}), "
        f"{len(nmods)} modules with >= 5 genera")


# --------------------------------------------------------------------- Fig. 3
def fig3(cons, qc, md, out, dd, n_pos, n_neg, name="Fig3_heatmap", only=None, pares=None):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    if pares is not None:
        k = cons.genus_a + "|" + cons.genus_b
        alvo = ["|".join(x) for x in pares] + ["|".join(x[::-1]) for x in pares]
        sel = cons[k.isin(alvo)].sort_values("mu", ascending=False)
    else:
        c2 = cons if only is None else cons[cons.genus_a.isin(only) | cons.genus_b.isin(only)]
        sel = pd.concat([c2[c2.mu > 0].nlargest(n_pos, "mu"), c2[c2.mu < 0].nsmallest(n_neg, "mu")])
    if sel.empty:
        log(f"  ({name} pulada: sem pares)"); return
    keys = set(sel.genus_a + "|" + sel.genus_b)
    q = qc[qc.status == "incluido"].copy()
    q["oc"] = q.cultura.map({c: i for i, c in enumerate(CROP_ORDER)}).fillna(7)
    q["og"] = q.grupo.map({"gramínea": 0, "leguminosa": 1}).fillna(2)
    q = q.sort_values(["og", "oc", "study"])
    studies = list(q.study)
    Z = pd.DataFrame(np.nan, index=list(keys), columns=studies)
    for s in studies:
        f = os.path.join(md, "effects", f"pairs_{s}.tsv.gz")
        d = pd.read_csv(f, sep="\t", usecols=["genus_a", "genus_b", "z"])
        d = d[(d.genus_a + "|" + d.genus_b).isin(keys)]
        for _, r in d.iterrows():
            Z.loc[r.genus_a + "|" + r.genus_b, s] = r.z
    order = list(sel.genus_a + "|" + sel.genus_b)
    Z = Z.loc[order]
    P, S = Z.shape
    cmap = LinearSegmentedColormap.from_list("div", [RED, "#f0efec", BLUE]); cmap.set_bad("white")
    fig = plt.figure(figsize=(7.2, max(0.18 * P + 1.6, 4.2)))
    h = 0.80
    axg = fig.add_axes([0.30, 0.08 + h + 0.035, 0.52, 0.02])
    axc = fig.add_axes([0.30, 0.08 + h + 0.01, 0.52, 0.02])
    axg.imshow(q.og.values[None, :], aspect="auto", cmap=ListedColormap([BLUE, ORANGE, GRAY]),
               vmin=0, vmax=2, interpolation="nearest")
    ccols = CROP_COR + [GRAY]
    axc.imshow(q.oc.clip(upper=7).values[None, :], aspect="auto", cmap=ListedColormap(ccols),
               vmin=0, vmax=7, interpolation="nearest")
    for a in (axg, axc):
        a.set_xticks([]); a.set_yticks([]); [sp.set_visible(False) for sp in a.spines.values()]
    ax = fig.add_axes([0.30, 0.08, 0.52, h])
    im = ax.imshow(np.ma.masked_invalid(Z.values), aspect="auto", cmap=cmap,
                   norm=TwoSlopeNorm(0, -0.8, 0.8), interpolation="nearest")
    ax.set_xticks(np.arange(-.5, S, 1), minor=True); ax.set_yticks(np.arange(-.5, P, 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.4); ax.tick_params(which="minor", length=0)
    nleg = int((q.og == 0).sum())
    ax.axvline(nleg - 0.5, color=INK, lw=0.8)
    ax.set_xticks([]); ax.set_xlabel("Studies")
    labs = [" – ".join(nome_bonito(g) for g in o.split("|")) for o in order]
    ax.set_yticks(range(P), labs, fontsize=6)
    for t, o in zip(ax.get_yticklabels(), order):          # pares com genero PGP em negrito
        if only is None and any(PGPR.match(g) for g in o.split("|")):
            t.set_fontweight("bold")
    ax.tick_params(axis="y", length=0); [sp.set_visible(False) for sp in ax.spines.values()]
    axf = fig.add_axes([0.835, 0.08, 0.15, h], sharey=ax)
    for i, (_, r) in enumerate(sel.iterrows()):
        col = BLUE if r.mu > 0 else RED
        axf.plot([r.r_pi_lo, r.r_pi_hi], [i, i], color=col, lw=0.9, alpha=0.5)
        axf.plot([r.r_ci_lo, r.r_ci_hi], [i, i], color=col, lw=2.6, solid_capstyle="round")
        axf.scatter(r.r_mu, i, s=12, marker="D", c=col, edgecolor="white", lw=0.5, zorder=3)
    axf.axvline(0, color=MUTED, lw=0.6)
    axf.tick_params(axis="y", which="both", left=False, length=0, labelleft=False)
    axf.spines["left"].set_visible(False); grid(axf, "x"); axf.set_ylim(P - 0.5, -0.5)
    axf.set_xlabel("Pooled effect (r)")
    # legendas abaixo da figura
    cax = fig.add_axes([0.30, -0.06, 0.2, 0.018])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal"); cb.outline.set_visible(False)
    cb.set_label("Partial correlation in each study (Fisher z)", fontsize=7); cb.ax.tick_params(labelsize=6)
    hs = [Patch(fc=BLUE, label="Grasses"), Patch(fc=ORANGE, label="Legumes")]
    hc = [Patch(fc=c, label=CROP_EN[k]) for k, c in zip(CROP_ORDER, CROP_COR)]
    hf = [Line2D([], [], marker="D", ls="", ms=5, mfc=BLUE, mec="white", label="Pooled effect"),
          Line2D([], [], color=BLUE, lw=2.6, label="95% confidence interval"),
          Line2D([], [], color=BLUE, lw=0.9, alpha=0.5, label="95% prediction interval"),
          Patch(fc="white", ec=LGRAY, label="Pair absent from study")]
    fig.legend(handles=hs + hc, loc="upper left", bbox_to_anchor=(0.55, -0.01), ncol=3,
               frameon=False, fontsize=6.5, columnspacing=1.0)
    fig.legend(handles=hf, loc="upper left", bbox_to_anchor=(0.02, -0.01), ncol=1,
               frameon=False, fontsize=6.5)
    save(fig, out, name)
    Z.to_csv(os.path.join(dd, f"{name}_z_by_study.tsv"), sep="\t")
    sel.to_csv(os.path.join(dd, f"{name}_pairs.tsv"), sep="\t", index=False)
    q[["study", "grupo", "cultura"]].to_csv(os.path.join(dd, f"{name}_study_order.tsv"), sep="\t", index=False)


# --------------------------------------------------------------- Figs. suplementares
def figS1(occ, out, dd):
    from matplotlib.lines import Line2D
    o = occ.copy(); o["abund_media"] = o.abund_media.astype(float)
    m = o.genus.map(lambda g: bool(PGPR.match(g)))
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    ax.scatter(o.ocupacao[~m], o.abund_media[~m], s=7, c=GRAY, alpha=0.45, lw=0, zorder=3)
    ax.scatter(o.ocupacao[m], o.abund_media[m], s=22, c=AQUA, edgecolor="white", lw=0.7, zorder=4)
    ax.set_yscale("log"); ax.set_xlim(0, 1.02); grid(ax, "both")
    ax.set_xlabel("Occupancy (fraction of studies)"); ax.set_ylabel("Mean relative abundance")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=6, mfc=AQUA, mec="white", label="PGP genera"),
                       Line2D([], [], marker="o", ls="", ms=5, mfc=GRAY, mec="white", label="Other genera")],
              loc="upper left", frameon=False, fontsize=7)
    save(fig, out, "FigS1_core")
    o.assign(pgpr=m).to_csv(os.path.join(dd, "FigS1_ocupacao.tsv"), sep="\t", index=False)


def figS2(md, out, dd):
    real = os.path.join(md, "meta", "meta_pairs.tsv.gz")
    nulos = sorted(glob.glob(os.path.join(md, "meta", "meta_pairs_nulo*.tsv.gz")))
    if not nulos:
        log("  (FigS2 pulada: sem resultado do modelo nulo)"); return
    a = pd.read_csv(real, sep="\t", usecols=["r_mu"]).r_mu
    b = pd.read_csv(nulos[0], sep="\t", usecols=["r_mu"]).r_mu
    bins = np.linspace(-0.5, 0.5, 101)
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    ax.hist(b, bins=bins, color=GRAY, alpha=0.8, zorder=2)
    ax.hist(a, bins=bins, histtype="step", color=BLUE, lw=1.2, zorder=3)
    ax.set_yscale("log"); grid(ax, "y")
    ax.set_xlabel("Pooled effect (r)"); ax.set_ylabel("Genus pairs")
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color=BLUE, lw=1.2, label="Observed data"),
                       Patch(fc=GRAY, alpha=0.8, label="Permutation null")],
              loc="upper left", frameon=False, fontsize=7)
    save(fig, out, "FigS2_null")
    pd.DataFrame({"r_real_quantis": a.quantile([.001, .01, .5, .99, .999]).values,
                  "r_nulo_quantis": b.quantile([.001, .01, .5, .99, .999]).values},
                 index=[.001, .01, .5, .99, .999]).to_csv(os.path.join(dd, "FigS2_quantis.tsv"), sep="\t")


def pgpr_summary(md, dd):
    m = pd.read_csv(os.path.join(md, "meta", "meta_pairs.tsv.gz"), sep="\t",
                    usecols=["genus_a", "genus_b", "k", "I2", "classe", "r_mu"])
    rows = []
    for col in ("genus_a", "genus_b"):
        d = m[m[col].map(lambda g: bool(PGPR.match(g)))].copy()
        d["pgp"] = d[col]; rows.append(d)
    d = pd.concat(rows)
    out = d.groupby("pgp").agg(pairs_tested=("classe", "size"),
                               consistent_pos=("classe", lambda x: (x == "consistente_pos").sum()),
                               consistent_neg=("classe", lambda x: (x == "consistente_neg").sum()),
                               associated=("classe", lambda x: x.str.startswith("associado").sum()),
                               median_I2=("I2", "median"), median_k=("k", "median")) \
        .sort_values("consistent_pos", ascending=False)
    out.to_csv(os.path.join(dd, "PGP_genera_summary.tsv"), sep="\t")
    log(f"    PGP genera tested: {len(out)} | with >= 1 consistent positive partner: "
        f"{(out.consistent_pos > 0).sum()}")
    return set(out.index)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--studies", required=True, help="subset_studies.tsv")
    ap.add_argument("--world", default="world_110m.geojson")
    ap.add_argument("--out", default="Figuras")
    ap.add_argument("--pares", default="pares_fig2.tsv")
    ap.add_argument("--min-abs-r", type=float, default=0.15, help="Fig2: |r| minimo das arestas")
    ap.add_argument("--n-pos", type=int, default=20, help="Fig3: pares positivos")
    ap.add_argument("--n-neg", type=int, default=10, help="Fig3: pares negativos")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root)); md = os.path.join(root, "merged")
    out = os.path.abspath(a.out); dd = os.path.join(out, "dados"); os.makedirs(dd, exist_ok=True)
    qc = pd.read_csv(os.path.join(md, "study_qc.tsv"), sep="\t", dtype=str).fillna("")
    # estudos sem efeitos (ex.: < 20 amostras apos excluir amostras fora da rizosfera) saem das figuras
    sem = qc.study.map(lambda s: not os.path.exists(os.path.join(md, "effects", f"pairs_{s}.tsv.gz")))
    qc.loc[(qc.status == "incluido") & sem, "status"] = "excluido"
    log(f"Estudos nas figuras: {(qc.status == 'incluido').sum()}")
    studies = pd.read_csv(a.studies, sep="\t", dtype=str).fillna("")
    fr = os.path.join(md, "meta", "rede_entrada.tsv")
    cons = pd.read_csv(fr if os.path.exists(fr) else os.path.join(md, "meta", "consistentes.tsv"), sep="\t")
    occ = pd.read_csv(os.path.join(md, "genus_occupancy.tsv"), sep="\t")
    carrega_linhagem(md)
    world = json.load(open(a.world))
    log(f"Gerando figuras em {out}")
    fig1(qc, studies, world, out, dd)
    pg = pgpr_summary(md, dd)
    fig2(cons, occ, out, dd, a.min_abs_r, name="FigS6_rede")
    fig2(cons, occ, out, dd, a.min_abs_r, name="FigS7_rede_PGP", only=pg)
    P = pd.read_csv(a.pares, sep="\t") if os.path.exists(a.pares) else None
    fig3(cons, qc, md, out, dd, a.n_pos, a.n_neg, name="Fig2_reprodutibilidade",
         pares=None if P is None else list(zip(P.genus_a, P.genus_b)))
    fig3(cons, qc, md, out, dd, a.n_pos, a.n_neg, name="FigS12_maiores_efeitos")
    fig3(cons, qc, md, out, dd, a.n_pos, a.n_neg, name="FigS8_heatmap_PGP", only=pg)
    figS1(occ, out, dd)
    figS2(md, out, dd)
    log("Pronto.")


if __name__ == "__main__":
    main()
