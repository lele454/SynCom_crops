#!/usr/bin/env python3
# [S3.12] Trait richness as a predictor of the positive-partner rate
# Supplementary 3 of the manuscript. See README.md for the full map.
"""Pergunta 6: traços de promoção de crescimento predizem a taxa de parceiros positivos?"""
import os, numpy as np, pandas as pd, statsmodels.api as sm
R = os.path.expanduser("~/syncom"); MM = f"{R}/merged/meta"

M = pd.read_csv(f"{MM}/meta_pairs.tsv.gz", sep="\t", usecols=["genus_a","genus_b","fdr","r_mu"])
testados = pd.concat([M.genus_a, M.genus_b]).value_counts().rename("pares_testados")
pos = M[(M.fdr < 0.05) & (M.r_mu > 0)]
parc = pd.concat([pos.genus_a, pos.genus_b]).value_counts().rename("parceiros_pos")

occ = pd.read_csv(f"{R}/merged/genus_occupancy.tsv", sep="\t").set_index("genus")
T = pd.read_csv(f"{R}/genomas/pgp/tracos_por_genero.tsv", sep="\t").set_index("genus")
L = pd.read_csv("generos_q6.tsv", sep="\t").set_index("genus")
tcols = [c for c in T.columns if c != "n_genomas"]
T["n_tracos"] = (T[tcols] >= 0.5).sum(axis=1)

df = pd.concat([testados, parc, occ[["ocupacao","abund_media"]],
                T[["n_tracos","n_genomas"]], L["grupo"]], axis=1)
df = df.dropna(subset=["pares_testados","n_tracos","ocupacao"]).fillna({"parceiros_pos":0})
df = df[df.pares_testados >= 10].copy()
df["log_ab"] = np.log10(df.abund_media.clip(lower=1e-9))
df["log_gen"] = np.log10(df.n_genomas.clip(lower=1))
res = sm.OLS(df.n_tracos, sm.add_constant(df.log_gen)).fit()
df["tracos_resid"] = res.resid
print(f"generos no teste: {len(df)} | {df.grupo.value_counts().to_dict()}")
print(f"tracos ~ log(genomas): R2={res.rsquared:.3f} (por isso usamos o residuo)")

def modelo(d, rotulo, preditor, extra=()):
    X = sm.add_constant(d[["log_ab","ocupacao",preditor] + list(extra)])
    m = sm.GLM(d.parceiros_pos, X, family=sm.families.Poisson(),
               offset=np.log(d.pares_testados)).fit(scale="X2")
    print(f"{rotulo:26s} n={len(d):4d}  beta={m.params[preditor]:+.4f}  p={m.pvalues[preditor]:.4g}  "
          f"IRR={np.exp(m.params[preditor]):.3f}  disp={m.scale:.2f}")

F = pd.read_csv(f"{MM}/vies_filos/por_genero.tsv", sep="\t", index_col=0)
if "filo" in F.columns:
    df["filo"] = F.filo.reindex(df.index)
    df["bacillota"] = df.filo.fillna("").str.match(r"^(Bacillota|Firmicutes)").astype(int)
    print(f"com filo atribuido: {df.filo.notna().sum()} | Bacillota: {int(df.bacillota.sum())}")
    print(f"tracos (mediana): Bacillota {df.loc[df.bacillota==1,'n_tracos'].median():.1f} | "
          f"demais {df.loc[df.bacillota==0,'n_tracos'].median():.1f}")

print("\n--- desfecho: parceiros positivos por par testavel")
modelo(df, "todos (residuo)", "tracos_resid")
modelo(df, "todos (bruto)", "n_tracos")
for g in ("candidato", "contraste"):
    s = df[df.grupo == g]
    if len(s) > 20:
        modelo(s, f"so {g}", "tracos_resid")
print("\nmedianas por grupo:")
print(df.groupby("grupo")[["n_tracos","n_genomas","parceiros_pos","pares_testados"]].median().to_string())

if "bacillota" in df:
    print("\n--- com o filo controlado (indicador de Bacillota no modelo)")
    modelo(df, "todos + Bacillota", "tracos_resid", ["bacillota"])
    s = df[df.bacillota == 0]
    modelo(s, "excluindo Bacillota", "tracos_resid")
    s = df[df.bacillota == 1]
    if len(s) > 20:
        modelo(s, "so Bacillota", "tracos_resid")
