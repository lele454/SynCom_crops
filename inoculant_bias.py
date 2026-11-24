#!/usr/bin/env python3
# [S3.17] Representation of commercial inoculant genera
# Supplementary 3 of the manuscript. See README.md for the full map.
"""
Os generos usados como inoculante comercial estao sub-representados na rede de associacoes?

Pergunta: generos como Pseudomonas e Bacillus formam menos associacoes consistentes do que
generos de abundancia e prevalencia equivalentes? Ou a ausencia deles no nucleo e so reflexo
de abundancia/poder estatistico?

Por genero (so os testados na meta-analise):
  grau_pos     parceiros com associacao positiva significativa (FDR < 0,05)
  grau_t1      parceiros na rede do nucleo tier 1, sinal positivo
  n_testados   pares em que o genero foi testavel (poder/oportunidade)
  covariaveis  log10 abundancia media, ocupacao (fracao de estudos), log n_testados
  grupos       inoculante comercial / PGP descrito (sem produto) / outros

Analises:
  1. Poisson (dispersao quasi-Poisson) com offset log(n_testados):
       grau_pos ~ inoculante + pgp + log_abund + ocupacao       (idem grau_t1)
     O coeficiente de "inoculante" e a razao de taxas de associacao, ja descontados
     abundancia, prevalencia e numero de pares testados.
  2. Pareamento por estratos (quintis de abundancia x tercis de ocupacao):
       razao de Mantel-Haenszel e diferenca media da taxa de associacao
  3. Tabela genero a genero: taxa observada e percentil dentro do proprio estrato
Nota: o desfecho binario "esta no nucleo (sim/nao)" foi retirado porque satura - generos testados
em centenas de pares quase sempre tem ao menos uma associacao. O informativo e a densidade
(associacoes por par testado).

Saidas: <root>/merged/meta/vies_inoculantes/ (por_genero.tsv, resumo.txt)
Uso: LD_PRELOAD=... python3 inoculant_bias.py --root ~/syncom
"""
import argparse, os, re, sys
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from consortia import PGP_CLASSIC, PGP_EXTRA

# Generos com inoculante/biofertilizante ou biocontrole bacteriano comercializado
# (conferir referencias e registros antes da submissao; a lista e o "tratamento" do teste)
INOCULANTE = ["Bradyrhizobium", "Rhizobium", "Sinorhizobium", "Ensifer", "Mesorhizobium",
              "Azospirillum", "Azotobacter", "Herbaspirillum", "Gluconacetobacter", "Nitrosomonas",
              "Bacillus", "Priestia", "Peribacillus", "Lysinibacillus", "Paenibacillus", "Brevibacillus",
              "Pseudomonas", "Burkholderia", "Paraburkholderia", "Pantoea", "Serratia", "Streptomyces",
              "Delftia", "Variovorax", "Enterobacter", "Klebsiella", "Kosakonia", "Rhodococcus",
              "Arthrobacter", "Micrococcus", "Lactobacillus", "Frateuria"]


def base(g):
    return re.sub(r"_[A-Z]+$", "", str(g))


def irls(y, X, offset, poisson=True, names=None):
    """Poisson (log, com offset) ou logistica, por IRLS; EP quasi para Poisson."""
    X = np.asarray(X, float); names = list(names or [])
    keep = [i for i in range(X.shape[1]) if np.nanstd(X[:, i]) > 0]      # tira colunas constantes
    X = X[:, keep]; names = [names[i] for i in keep]
    X = np.column_stack([np.ones(len(y)), X]); names = ["intercepto"] + names
    b = np.zeros(X.shape[1]); lo = np.log(offset) if poisson else 0
    for _ in range(100):
        eta = X @ b + lo
        if poisson:
            mu = np.exp(np.clip(eta, -30, 30)); W = mu; z = eta - lo + (y - mu) / np.maximum(mu, 1e-9)
        else:
            mu = 1 / (1 + np.exp(-eta)); W = np.maximum(mu * (1 - mu), 1e-9)
            z = eta + (y - mu) / W
        H = X.T @ (X * W[:, None])
        bn = np.linalg.lstsq(H, X.T @ (W * z), rcond=None)[0]
        if np.abs(bn - b).max() < 1e-9:
            b = bn; break
        b = bn
    V = np.linalg.pinv(H)
    phi = 1.0
    if poisson:
        phi = max(((y - mu) ** 2 / np.maximum(mu, 1e-9)).sum() / (len(y) - X.shape[1]), 1.0)
    se = np.sqrt(np.diag(V) * phi); z_ = b / se
    return pd.DataFrame({"coef": b, "razao": np.exp(b), "EP": se, "z": z_,
                         "p": 2 * stats.norm.sf(np.abs(z_))}, index=names), phi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, "vies_inoculantes"); os.makedirs(out, exist_ok=True)

    m = pd.read_csv(os.path.join(mm, "meta_pairs.tsv.gz"), sep="\t",
                    usecols=["genus_a", "genus_b", "mu", "fdr"])
    pos = m[(m.mu > 0) & (m.fdr < 0.05)]
    core = pd.read_csv(os.path.join(mm, "nucleo_robusto.tsv"), sep="\t")
    t1 = core[(core.tier == 1) & (core.sinal == "positivo")]
    cnt = lambda d: pd.concat([d.genus_a, d.genus_b]).value_counts()
    G = pd.DataFrame({"n_testados": cnt(m), "grau_pos": cnt(pos), "grau_t1": cnt(t1)}).fillna(0)
    G = G[G.n_testados >= 10]
    occ = pd.read_csv(os.path.join(md, "genus_occupancy.tsv"), sep="\t").set_index("genus")
    G = G.join(occ[["ocupacao", "abund_media"]])
    G = G.dropna(subset=["ocupacao", "abund_media"])
    G["log_abund"] = np.log10(G.abund_media.clip(lower=1e-8))
    ino = {base(x) for x in INOCULANTE}
    pgp = {base(x) for x in PGP_CLASSIC + PGP_EXTRA} - ino
    G["grupo"] = np.where(G.index.map(lambda g: base(g) in ino), "inoculante",
                          np.where(G.index.map(lambda g: base(g) in pgp), "PGP sem produto", "outros"))
    G["inoculante"] = (G.grupo == "inoculante").astype(int)
    G["pgp"] = (G.grupo == "PGP sem produto").astype(int)
    G["taxa_pos"] = G.grau_pos / G.n_testados
    G["no_nucleo"] = (G.grau_t1 > 0).astype(int)

    L = [f"Generos testados (>= 10 pares): {len(G)} | " +
         " | ".join(f"{k}: {v}" for k, v in G.grupo.value_counts().items()),
         f"Taxa de associacao positiva (grau/pares testados) - mediana por grupo: " +
         ", ".join(f"{k} {v:.4f}" for k, v in G.groupby('grupo').taxa_pos.median().items()), ""]

    cols = ["inoculante", "pgp", "log_abund", "ocupacao"]
    for y, lab in (("grau_pos", "grau positivo (FDR<0,05)"), ("grau_t1", "grau no nucleo tier 1")):
        if G[y].nunique() < 2:
            L += [f"1. {lab}: sem variacao entre generos, analise nao aplicavel", ""]; continue
        tab, phi = irls(G[y].values.astype(float), G[cols].values, G.n_testados.values, True, cols)
        L.append(f"1. Poisson com offset log(pares testados) - {lab} (dispersao quasi = {phi:.1f}):")
        L.append(tab.to_string(float_format=lambda v: f"{v:.3g}"))
        L.append("")
    G["estrato"] = (pd.qcut(G.log_abund.rank(method="first"), 5, labels=False).astype(str) + "_" +
                    pd.qcut(G.ocupacao.rank(method="first"), 3, labels=False).astype(str))
    sub = G[G.grupo != "PGP sem produto"]
    L.append("2. Pareado por estratos (abundancia x ocupacao), inoculante vs outros:")
    e = sub.groupby("estrato").taxa_pos.mean()
    sub = sub.assign(esperado=sub.estrato.map(e))
    agg = sub.groupby("grupo").apply(lambda d: pd.Series(
        {"n": len(d), "taxa_obs": d.taxa_pos.mean(), "taxa_esperada_estrato": d.esperado.mean()}),
        include_groups=False)
    agg["razao_obs_esp"] = agg.taxa_obs / agg.taxa_esperada_estrato
    L.append("   taxa de associacao observada vs esperada pelo estrato:")
    L.append(agg.to_string(float_format=lambda v: f"{v:.4g}"))

    G["pct_no_estrato"] = G.groupby("estrato").taxa_pos.rank(pct=True)
    show = G[G.inoculante == 1].sort_values("taxa_pos", ascending=False)
    L += ["", "3. Generos de inoculante, um a um (pct_no_estrato: percentil da taxa entre generos de "
               "abundancia e ocupacao parecidas):",
          show[["n_testados", "grau_pos", "grau_t1", "taxa_pos", "pct_no_estrato", "ocupacao",
                "abund_media"]].to_string(float_format=lambda v: f"{v:.4g}")]
    pr = G.loc[G.inoculante == 1, "pct_no_estrato"]
    L.append(f"   percentil mediano = {pr.median():.2f} (0,5 = igual aos pares; < 0,5 = sub-representados) | "
             f"Wilcoxon vs 0,5: p = {stats.wilcoxon(pr - 0.5).pvalue:.3g}")
    G.to_csv(os.path.join(out, "por_genero.tsv"), sep="\t", float_format="%.5g")
    txt = "\n".join(L)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
