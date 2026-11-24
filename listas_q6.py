#!/usr/bin/env python3
# [S3.13] Contrast group matched on abundance and occupancy
# Supplementary 3 of the manuscript. See README.md for the full map.
"""Lista candidatos (consorcios) e um grupo de contraste pareado por abundancia/ocupacao."""
import os, numpy as np, pandas as pd
R = os.path.expanduser("~/syncom"); MM = f"{R}/merged/meta"
cand = set(pd.read_csv(f"{R}/genomas/isolados/por_genero.tsv", sep="\t", index_col=0).index)
nr = pd.read_csv(f"{MM}/nucleo_robusto.tsv", sep="\t")
rede = set(nr.genus_a) | set(nr.genus_b)
occ = pd.read_csv(f"{R}/merged/genus_occupancy.tsv", sep="\t").set_index("genus")
occ = occ[occ.index.isin(rede)].copy()
occ["la"] = np.log10(occ.abund_media.clip(lower=1e-9))
pool = occ[~occ.index.isin(cand)]
alvo = occ[occ.index.isin(cand)]
print(f"candidatos: {len(cand)} ({len(alvo)} com ocupacao) | rede: {len(rede)} | pool de contraste: {len(pool)}")
usados, ctrl = set(), []
for g, r in alvo.iterrows():
    d = (pool.la - r.la).abs() + (pool.ocupacao - r.ocupacao).abs() * 2
    for c in d.sort_values().index:
        if c not in usados:
            usados.add(c); ctrl.append(c); break
print(f"contraste pareado: {len(ctrl)}")
tudo = sorted(set(alvo.index) | set(ctrl))
open("generos_q6.txt", "w").write(",".join(tudo))
pd.DataFrame({"genus": tudo,
              "grupo": ["candidato" if g in cand else "contraste" for g in tudo]}) \
  .to_csv("generos_q6.tsv", sep="\t", index=False)
print(f"total a anotar: {len(tudo)} generos -> generos_q6.txt")
