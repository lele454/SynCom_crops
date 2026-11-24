#!/usr/bin/env python3
# [S4.3] Functional complementarity
# Supplementary 4 of the manuscript. See README.md for the full map.
"""
Complementaridade funcional das combinacoes: acrescenta a estrelas.tsv o numero de
membros anotados, de tracos reunidos e de tracos complementares.

Um traco e considerado presente num genero quando ocorre em pelo menos --min-frac
de seus genomas anotados. Complementares = presentes em parte dos membros e
ausentes em outros.
Saida: <root>/merged/meta/consorcios/complementaridade.tsv
Uso: python3 complementarity.py --root ~/syncom
"""
import argparse, os, shutil, pandas as pd


def membros(s):
    s = str(s)
    for sep in (" + ", ";", ","):
        if sep in s:
            return [x.strip() for x in s.split(sep) if x.strip()]
    return [s.strip()] if s.strip() else []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--min-frac", type=float, default=0.5)
    ap.add_argument("--entrada", default="estrelas.tsv")
    a = ap.parse_args()
    cd = os.path.join(a.root, "merged", "meta", "consorcios")
    E = pd.read_csv(os.path.join(cd, a.entrada), sep="\t")
    T = pd.read_csv(os.path.join(a.root, "genomas", "pgp", "tracos_por_genero.tsv"),
                    sep="\t").set_index("genus")
    tcols = [c for c in T.columns if c != "n_genomas"]
    pres = {g: {c for c in tcols if r[c] >= a.min_frac} for g, r in T.iterrows()}
    print(f"generos anotados: {len(pres)} | tracos avaliados: {len(tcols)} | "
          f"combinacoes: {len(E)}")

    n_ok, compl, total = [], [], []
    for _, r in E.iterrows():
        m = membros(r.membros)
        com = [g for g in m if g in pres]
        sets = [pres[g] for g in com]
        uni = set().union(*sets) if sets else set()
        inter = set.intersection(*sets) if sets else set()
        n_ok.append(len(com)); total.append(len(uni)); compl.append(len(uni - inter))
    E["n_com_tracos"] = n_ok
    E["tracos_compl"] = compl
    E["tracos_total"] = total

    f = os.path.join(cd, "complementaridade.tsv")
    if os.path.exists(f):
        shutil.copy(f, f + ".bak")
    E.to_csv(f, sep="\t", index=False, float_format="%.5g")
    A = E[E.n_com_tracos == E.tamanho]
    print(f"\ncombinacoes totalmente anotadas: {len(A)} de {len(E)}")
    print(A.groupby("tamanho").agg(n=("membros", "size"), copres=("copresenca", "median"),
                                   compl=("tracos_compl", "median"),
                                   total=("tracos_total", "median")).round(3).to_string())
    print(f"\nescrito em {f}")


if __name__ == "__main__":
    main()
