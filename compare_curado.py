#!/usr/bin/env python3
# [S1.14] Curated against uncurated analysis
# Supplementary 1 of the manuscript. See README.md for the full map.
"""
Compara a analise principal com a analise curada (fatores curados + amostras fora da rizosfera excluidas).

- correlacao de Spearman dos efeitos (pares testados nas duas)
- nucleo robusto: fracao dos pares de cada tier confirmados na curada (mesmo sinal, FDR < 0,05)
- pares dos consorcios da shortlist: confirmados ou nao

Requer: meta_pairs.tsv.gz, meta_pairs_curado.tsv.gz (meta_pairs.py --effects effects_curado --tag curado),
        nucleo_robusto.tsv, consorcios/shortlist.tsv
Uso: LD_PRELOAD=... python3 compare_curado.py --root ~/syncom
"""
import argparse, itertools, os
import numpy as np, pandas as pd
from scipy import stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    mm = os.path.join(os.path.abspath(os.path.expanduser(a.root)), "merged", "meta")
    key = lambda d: d.genus_a + "|" + d.genus_b
    m0 = pd.read_csv(os.path.join(mm, "meta_pairs.tsv.gz"), sep="\t"); m0.index = key(m0)
    m1 = pd.read_csv(os.path.join(mm, "meta_pairs_curado.tsv.gz"), sep="\t"); m1.index = key(m1)
    both = m0.index.intersection(m1.index)
    rho = stats.spearmanr(m0.loc[both, "mu"], m1.loc[both, "mu"]).correlation
    L = [f"Pares testados: principal {len(m0)} | curada {len(m1)} | em comum {len(both)} | rho = {rho:.3f}",
         "Classes na curada: " + ", ".join(f"{k} {v}" for k, v in m1.classe.value_counts().items()), ""]
    core = pd.read_csv(os.path.join(mm, "nucleo_robusto.tsv"), sep="\t"); core.index = key(core)
    j = core.join(m1[["mu", "fdr", "r_mu", "classe"]].add_suffix("_cur"))
    j["status"] = np.where(j.mu_cur.isna(), "nao_testavel",
                           np.where((np.sign(j.mu_cur) == np.sign(j.r_mu)) & (j.fdr_cur < 0.05), "confirmado",
                                    "nao_confirmado"))
    L.append("Nucleo robusto na analise curada:")
    L.append(pd.crosstab([j.tier, j.sinal], j.status, margins=True, margins_name="total").to_string())
    sl = pd.read_csv(os.path.join(mm, "consorcios", "shortlist.tsv"), sep="\t")
    L += ["", "Pares dos consorcios da shortlist (r principal -> r curado):"]
    for _, r in sl.iterrows():
        mem = [r.centro] + r.parceiros.split(" + ")
        parts = []
        for x, y in itertools.combinations(sorted(mem), 2):
            k = f"{x}|{y}"
            if k in j.index and k in core.index:
                rr = j.loc[k]
                parts.append(f"{x}-{y} {rr.r_mu:.2f}->{'' if pd.isna(rr.r_mu_cur) else f'{rr.r_mu_cur:.2f}'} {rr.status}")
        L.append(f"  {r.id}: " + " | ".join(parts))
    txt = "\n".join(L)
    open(os.path.join(mm, "comparacao_curado.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
