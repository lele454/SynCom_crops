#!/usr/bin/env python3
# [S6.10] Body and supplementary tables
# Supplementary 6 of the manuscript. See README.md for the full map.
"""
Monta as tabelas do manuscrito a partir dos dados-fonte das figuras.

Gera: Table1_proposed_sets, TableS3_management, TableS5_strains,
      TableS7_pgp_genera, TableS8_phylum_mediation, TableS9_sensitivity

Nao gera (fontes ainda nao localizadas): S4 validacao temporal e S6 diagnostico
de zeros par a par.

Uso: python3 build_tables.py --dados Figuras/dados --out Tabelas
"""
import argparse, os
import pandas as pd

TABELA1 = [
    ("Recommended trio", "Arthrobacter + Sphingomicrobium + Streptomyces", 0.74,
     "highest co-presence among testable combinations, in 90 studies and all seven crops; "
     "all member strains accessible"),
    ("Contrast trio", "Devosia_A + Phenylobacterium + Ramlibacter", 0.45,
     "no member with reported growth promotion; 56 studies and all seven crops"),
    ("Executable set", "Kribbella + Mycobacterium + Nocardioides + Pseudonocardia + Streptomyces", 0.42,
     "largest replicated set with accessible strains; strain screening required for Mycobacterium"),
    ("Replicated set", "six genera, five without a cultured representative", 0.33,
     "all internal pairs replicate in the validation studies; cannot be assembled"),
    ("Densest set", "ten genera, in 45 equivalent sets", 0.02,
     "largest mutually associated set; composition is substitutable"),
]

TABELAS9 = [
    ("Studies with the complete V4 region only", "81%",
     "eight forward-read-only studies discarded"),
    ("Two-stage meta-analysis by country", "88%",
     "confidence interval 1.29x wider"),
    ("Exclusion of sugarcane studies", "91%",
     "rho = 0.996 between networks; three studies excluded"),
    ("Stratification by mineral nitrogen fertilisation", "87%",
     "median difference in z of -0.026; no significant moderation"),
    ("Idem, associations involving rhizobia", "86%", "-"),
    ("Minimum depth of 5,000 instead of 2,000 reads", "97%",
     "rho = 0.998 between effect sizes"),
]


TABELAS8 = [
    ("16S gene copies per genome", "16%"),
    ("r versus k life strategy", "5%"),
    ("Genome size", "0%"),
    ("Guanine and cytosine content", "-2%"),
    ("Genus splitting in the GTDB taxonomy", "-3%"),
    ("The three genomic variables together", "26%"),
    ("Endospore formation", "not estimable"),
]


def primeiro(d, *nomes):
    for n in nomes:
        p = os.path.join(d, n)
        if os.path.exists(p):
            return p
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dados", default=os.path.join("Figuras", "dados"))
    ap.add_argument("--out", default="Tabelas")
    a = ap.parse_args()
    d, out = a.dados, a.out
    os.makedirs(out, exist_ok=True)

    def grava(df, nome):
        df.to_csv(os.path.join(out, nome + ".tsv"), sep="\t", index=False)
        print(f"  {nome}.tsv  {len(df)} linhas x {df.shape[1]} colunas")

    grava(pd.DataFrame(TABELA1, columns=["set", "composition", "co_presence", "credential"]),
          "Table1_proposed_sets")

    grava(pd.DataFrame(TABELAS9, columns=["test", "associations_retained", "effect_on_estimate"]),
          "TableS9_sensitivity")

    grava(pd.DataFrame(TABELAS8, columns=["explanation_tested", "attenuation_of_phylum_effect"]),
          "TableS8_phylum_mediation")

    p = primeiro(d, "FigS5_moderadores.tsv")
    if p:
        M = pd.read_csv(p, sep="\t")
        cols = [c for c in ("genus_a", "genus_b", "contraste", "k", "r_mu", "r_ci_lo", "r_ci_hi",
                            "r_global", "r_N_mineral", "r_controle", "fdr") if c in M.columns]
        grava(M[cols].sort_values("r_mu", ascending=False), "TableS3_management")
    else:
        print("  !! FigS5_moderadores.tsv nao encontrado")

    p = primeiro(d, "Fig3_corpo.tsv", "Fig5_cultivabilidade.tsv", "Fig4_cultivabilidade.tsv")
    if p:
        G = pd.read_csv(p, sep="\t")
        cols = [c for c in ("genus", "classe", "genomas", "genomas_isolados", "linhagens_tipo",
                            "especies", "isolados_de_solo", "isolados_clinicos", "frac_isolados",
                            "nome_base") if c in G.columns]
        grava(G[cols].sort_values(["classe", "genus"]), "TableS5_strains")
    else:
        print("  !! tabela de linhagens por genero nao encontrada")

    p = primeiro(d, "PGP_genera_summary.tsv")
    if p:
        P = pd.read_csv(p, sep="\t")
        grava(P.sort_values("associated", ascending=False), "TableS7_pgp_genera")
    else:
        print("  !! PGP_genera_summary.tsv nao encontrado")

    print("\nFaltam S4 (validacao temporal) e S6 (diagnostico de zeros par a par).")
    print(f"Saida em {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
