#!/usr/bin/env python3
# [S3.1] Confirmation tiers
# Supplementary 3 of the manuscript. See README.md for the full map.
"""
Nucleo robusto: classifica os pares consistentes da analise principal pela
confirmacao nas analises de sensibilidade (sem China; so estudos com fatores de desenho).

Confirmado numa sensibilidade = par testavel la (k >= 10), mesmo sinal e FDR < 0,05.
  tier 1 (nucleo robusto) : confirmado nas duas
  tier 2                  : confirmado em uma; na outra nao testavel ou nao confirmado
  tier 3                  : consistente so na analise principal

Requer: meta_pairs.tsv.gz (principal) e sensibilidade/meta_sem_China.tsv.gz,
        sensibilidade/meta_com_fatores.tsv.gz (gerados por sensitivity.py)
Saidas: <root>/merged/meta/nucleo_robusto.tsv e resumo no terminal
Uso: LD_PRELOAD=... python3 robust_core.py --root ~/syncom
"""
import argparse, os, re, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meta_pairs import entrada

PGPR = re.compile(r"^(Pseudomonas|Bacillus|Paenibacillus|Streptomyces|Bradyrhizobium|Rhizobium|"
                  r"Mesorhizobium|Sinorhizobium|Ensifer|Azospirillum|Herbaspirillum|Burkholderia|"
                  r"Paraburkholderia|Enterobacter|Pantoea|Serratia|Klebsiella|Azotobacter|"
                  r"Gluconacetobacter|Arthrobacter|Microbacterium|Variovorax|Brucella|Stenotrophomonas|"
                  r"Lysobacter|Flavobacterium|Achromobacter|Delftia|Acinetobacter|Kosakonia|"
                  r"Priestia|Peribacillus|Lysinibacillus|Rhodococcus|Agrobacterium)(_[A-Z]+)?$")


def key(d):
    return d.genus_a + "|" + d.genus_b


def confirm(main, sens, tag):
    s = sens.set_index(key(sens))[["k", "mu", "fdr", "r_mu"]].astype(float)
    s.columns = [f"{c}_{tag}" for c in s.columns]
    m = main.join(s, on="pair")
    testable = m[f"k_{tag}"].notna()
    ok = testable & (np.sign(m[f"mu_{tag}"]) == np.sign(m.mu)) & (m[f"fdr_{tag}"] < 0.05)
    m[f"status_{tag}"] = np.where(~testable, "nao_testavel", np.where(ok, "confirmado", "nao_confirmado"))
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    md = os.path.join(os.path.abspath(os.path.expanduser(a.root)), "merged", "meta")
    m = pd.read_csv(os.path.join(md, "meta_pairs.tsv.gz"), sep="\t")
    m = m[entrada(m)].copy()
    m["pair"] = key(m)
    for tag, f in (("semChina", "meta_sem_China.tsv.gz"), ("comFatores", "meta_com_fatores.tsv.gz")):
        m = confirm(m, pd.read_csv(os.path.join(md, "sensibilidade", f), sep="\t"), tag)
    c1 = m.status_semChina == "confirmado"; c2 = m.status_comFatores == "confirmado"
    m["tier"] = np.where(c1 & c2, 1, np.where(c1 | c2, 2, 3))
    m["sinal"] = np.where(m.mu > 0, "positivo", "negativo")
    m["pgp"] = [bool(PGPR.match(x)) or bool(PGPR.match(y)) for x, y in zip(m.genus_a, m.genus_b)]
    cols = ["tier", "sinal", "genus_a", "genus_b", "pgp", "k", "r_mu", "r_ci_lo", "r_ci_hi", "r_pi_lo",
            "r_pi_hi", "I2", "fdr", "status_semChina", "r_mu_semChina", "k_semChina",
            "status_comFatores", "r_mu_comFatores", "k_comFatores"]
    m = m.sort_values(["tier", "sinal", "r_mu"], ascending=[True, False, False])[cols]
    m.to_csv(os.path.join(md, "nucleo_robusto.tsv"), sep="\t", index=False, float_format="%.4g")

    print(pd.crosstab(m.tier, m.sinal, margins=True, margins_name="total").to_string())
    t1 = m[(m.tier == 1) & (m.sinal == "positivo")]
    gen = pd.unique(t1[["genus_a", "genus_b"]].values.ravel())
    print(f"\nTier 1 positivo: {len(t1)} pares, {len(gen)} generos | com genero PGP: {int(t1.pgp.sum())} pares")
    print("\nTier 1 positivos (top 15 por efeito):")
    for _, r in t1.head(15).iterrows():
        print(f"  {r.genus_a} - {r.genus_b}: r={r.r_mu:.2f} [PI {r.r_pi_lo:.2f}; {r.r_pi_hi:.2f}] "
              f"k={r.k}{'  (PGP)' if r.pgp else ''}")
    tp = t1[t1.pgp]
    if len(tp):
        print("\nTier 1 positivos com genero PGP:")
        for _, r in tp.iterrows():
            print(f"  {r.genus_a} - {r.genus_b}: r={r.r_mu:.2f} k={r.k}")


if __name__ == "__main__":
    main()
