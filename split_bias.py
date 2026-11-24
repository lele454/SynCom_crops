#!/usr/bin/env python3
# [S3.16] GTDB taxonomic splitting
# Supplementary 3 of the manuscript. See README.md for the full map.
"""
Ultimo teste do mecanismo (etapa 1): fatiamento taxonomico do GTDB.

Hipotese: o GTDB dividiu os generos de Bacillota mais que os de outros filos (Bacillus,
Bacillus_A, Bacillus_BD, Bacillus_G, Bacillus_H, Priestia, Peribacillus...). Se linhagens
ecologicamente parecidas foram repartidas, cada fatia aparece em um subconjunto das amostras,
com mais zeros, e o sinal de associacao se dilui entre as fatias.

CRITERIO DEFINIDO ANTES DE VER O RESULTADO (evita garimpagem de hipoteses):
  o mecanismo e considerado sustentado se, e somente se,
    (i)  o coeficiente de log(fatias) for negativo com p < 0,05, E
    (ii) a atenuacao do coeficiente de Bacillota for >= 50%.
  Caso contrario, a linha mecanistica se encerra e o fenomeno e reportado como descritivo.

Metricas por genero:
  fatias_gtdb    generos do GTDB que compartilham o mesmo nome-base (sp_clusters)
  fatias_dados   fatias do mesmo nome-base que aparecem na nossa meta-analise (as que de fato
                 competem pelo sinal)
Modelos de Poisson quasi (desfecho = parceiros com associacao positiva; offset = pares testados):
  1. base + Bacillota
  2. base + log(fatias_gtdb) + log(fatias_dados)
  3. 2 + Bacillota   -> atenuacao

Saidas: <root>/merged/meta/vies_filos/ (fatias_por_genero.tsv, fatias_resumo.txt)
Uso: LD_PRELOAD=... python3 split_bias.py --root ~/syncom
"""
import argparse, os, sys
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inoculant_bias import irls, base
from spore_bias import tax_map
from consortia import pk  # noqa: F401  (mantem a importacao do modulo de utilidades)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, "vies_filos"); os.makedirs(out, exist_ok=True)

    import glob, gzip
    f = sorted(glob.glob(os.path.join(root, "genomas", "sp_clusters*")))
    if not f:
        raise SystemExit("sp_clusters do GTDB nao encontrado em ~/syncom/genomas")
    op = gzip.open if f[0].endswith(".gz") else open
    with op(f[0], "rt") as fh:
        sp = pd.read_csv(fh, sep="\t", low_memory=False)
    col = next(c for c in sp.columns if "gtdb taxonomy" in c.lower())
    gtdb = sp[col].str.extract(r"g__([^;]+)")[0].dropna().unique()
    fat_gtdb = pd.Series([base(g) for g in gtdb]).value_counts()

    m = pd.read_csv(os.path.join(mm, "meta_pairs.tsv.gz"), sep="\t",
                    usecols=["genus_a", "genus_b", "mu", "fdr"])
    pos = m[(m.mu > 0) & (m.fdr < 0.05)]
    cnt = lambda d: pd.concat([d.genus_a, d.genus_b]).value_counts()
    G = pd.DataFrame({"n_testados": cnt(m), "grau_pos": cnt(pos)}).fillna(0)
    G = G[G.n_testados >= 10]
    occ = pd.read_csv(os.path.join(md, "genus_occupancy.tsv"), sep="\t").set_index("genus")
    G = G.join(occ[["ocupacao", "abund_media"]]).join(tax_map(md)).dropna(subset=["ocupacao", "abund_media"])
    G["nome_base"] = [base(g) for g in G.index]
    G["fatias_gtdb"] = G.nome_base.map(fat_gtdb).fillna(1)
    G["fatias_dados"] = G.nome_base.map(G.nome_base.value_counts())
    G["log_fatias_gtdb"] = np.log(G.fatias_gtdb)
    G["log_fatias_dados"] = np.log(G.fatias_dados)
    G["log_abund"] = np.log10(G.abund_media.clip(lower=1e-8))
    G["bacillota"] = G.filo.fillna("").str.match(r"^(Bacillota|Firmicutes)").astype(int)
    G["taxa_pos"] = G.grau_pos / G.n_testados
    y = G.grau_pos.values.astype(float); off = G.n_testados.values

    L = ["Critério definido antes do resultado: mecanismo sustentado se o coeficiente de "
         "log(fatias) for negativo com p < 0,05 E a atenuacao de Bacillota for >= 50%.", "",
         f"Generos: {len(G)} | fatias no GTDB (mediana): Bacillota "
         f"{G.loc[G.bacillota == 1, 'fatias_gtdb'].median():.0f} | demais "
         f"{G.loc[G.bacillota == 0, 'fatias_gtdb'].median():.0f}",
         f"fatias presentes nos nossos dados (mediana): Bacillota "
         f"{G.loc[G.bacillota == 1, 'fatias_dados'].median():.0f} | demais "
         f"{G.loc[G.bacillota == 0, 'fatias_dados'].median():.0f}", "",
         "Fatiamento por filo (filos com >= 8 generos):",
         G.dropna(subset=["filo"]).groupby("filo").filter(lambda d: len(d) >= 8)
          .groupby("filo").agg(n_generos=("fatias_gtdb", "size"),
                               fatias_gtdb_mediana=("fatias_gtdb", "median"),
                               fatias_dados_mediana=("fatias_dados", "median"),
                               frac_com_sufixo=("nome_base", lambda s: np.mean(
                                   [b != g for b, g in zip(s, s.index)])))
          .sort_values("fatias_gtdb_mediana", ascending=False)
          .to_string(float_format=lambda v: f"{v:.3g}"), ""]

    bcols = ["log_abund", "ocupacao"]
    t1, p1 = irls(y, G[bcols + ["bacillota"]].values, off, True, bcols + ["bacillota"])
    c2 = bcols + ["log_fatias_gtdb", "log_fatias_dados"]
    t2, p2 = irls(y, G[c2].values, off, True, c2)
    t3, p3 = irls(y, G[c2 + ["bacillota"]].values, off, True, c2 + ["bacillota"])
    L += [f"1. Base + Bacillota (quasi = {p1:.1f}):", t1.to_string(float_format=lambda v: f"{v:.3g}"), "",
          f"2. Base + fatiamento (quasi = {p2:.1f}):", t2.to_string(float_format=lambda v: f"{v:.3g}"), "",
          f"3. 2 + Bacillota (quasi = {p3:.1f}):", t3.to_string(float_format=lambda v: f"{v:.3g}"), ""]
    c1, c3 = t1.loc["bacillota", "coef"], t3.loc["bacillota", "coef"]
    aten = 100 * (1 - c3 / c1) if c1 else np.nan
    neg = [(k, t2.loc[k, "coef"], t2.loc[k, "p"]) for k in ("log_fatias_gtdb", "log_fatias_dados")
           if k in t2.index]
    if not neg:
        L.append("Fatiamento sem variacao entre generos; teste nao aplicavel")
    ok_i = bool(neg) and any(c < 0 and p < 0.05 for _, c, p in neg)
    ok_ii = aten >= 50
    L += [f"Mediacao pelo fatiamento: coeficiente de Bacillota {c1:+.3f} -> {c3:+.3f} "
          f"(atenuacao de {aten:.0f}%)",
          "   " + " | ".join(f"{k}: coef {c:+.3f}, p = {p:.3g}" for k, c, p in neg), "",
          f"CRITERIO (i) coeficiente de fatias negativo e significativo: {'SIM' if ok_i else 'NAO'}",
          f"CRITERIO (ii) atenuacao >= 50%: {'SIM' if ok_ii else 'NAO'}",
          "CONCLUSAO: " + ("mecanismo SUSTENTADO - vale a etapa 2 (reagregar por nome-base/familia "
                           "e refazer a meta-analise)" if (ok_i and ok_ii) else
                           "mecanismo NAO sustentado - encerrar a linha mecanistica e reportar o "
                           "fenomeno como descritivo, com as explicacoes candidatas na discussao")]
    r = stats.spearmanr(G.fatias_dados, G.taxa_pos)
    L.append(f"   (correlacao bruta fatias nos dados x taxa de associacao: rho = {r.statistic:.3f}, "
             f"p = {r.pvalue:.3g})")
    G.to_csv(os.path.join(out, "fatias_por_genero.tsv"), sep="\t", float_format="%.5g")
    txt = "\n".join(L)
    open(os.path.join(out, "fatias_resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
