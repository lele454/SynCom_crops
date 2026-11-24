#!/usr/bin/env python3
# [S4.11] Sampling against rarity as the cause of absence in a crop
# Supplementary 4 of the manuscript. See README.md for the full map.
"""
Por que alguns consorcios nao tem dados em certas culturas? Amostragem ou raridade?

Para cada cultura: n de estudos e amostras; para cada genero da shortlist:
  ocupacao (fracao de estudos onde aparece), prevalencia mediana nas amostras,
  fracao de estudos onde e testavel (prevalencia >= --min-prev, criterio de effect_sizes.py),
  abundancia relativa mediana (nos estudos onde aparece)
Para cada consorcio: n de estudos da cultura onde TODOS os membros sao testaveis.

Leitura:
  poucos estudos na cultura (< 5)             -> limitacao de amostragem
  muitos estudos, membros raros/nao testaveis -> raridade biologica na cultura

Saidas (<root>/merged/meta/consorcios/): cobertura_generos.tsv, cobertura_consorcios.tsv,
  cobertura_resumo.txt
Uso: LD_PRELOAD=... python3 crop_coverage.py --root ~/syncom
"""
import argparse, os
import numpy as np, pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--min-prev", type=float, default=0.2)
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); out = os.path.join(md, "meta", "consorcios")

    sl = pd.read_csv(os.path.join(out, "shortlist.tsv"), sep="\t")
    mem = {r.id: [r.centro] + r.parceiros.split(" + ") for _, r in sl.iterrows()}
    gen = sorted({g for v in mem.values() for g in v})

    qc = pd.read_csv(os.path.join(md, "study_qc.tsv"), sep="\t")
    qc = qc[qc.status == "incluido"]
    per = []                                    # estudo x genero
    for s, c in zip(qc.study, qc.cultura.fillna("NA")):
        d = pd.read_csv(os.path.join(md, "genus", f"{s}.tsv.gz"), sep="\t", index_col=0)
        rel = d.div(d.sum(axis=1), axis=0)
        for g in gen:
            if g in d.columns:
                x = d[g] > 0
                per.append((s, c, len(d), g, x.mean(), rel.loc[x, g].median() if x.any() else 0.0))
            else:
                per.append((s, c, len(d), g, 0.0, 0.0))
    P = pd.DataFrame(per, columns=["study", "cultura", "n", "genus", "prev", "abund"])
    P["testavel"] = P.prev >= a.min_prev

    crop_n = P.drop_duplicates("study").groupby("cultura").agg(estudos=("study", "size"), amostras=("n", "sum"))
    G = P.groupby(["cultura", "genus"]).agg(
        ocupacao=("prev", lambda v: (v > 0).mean()), prev_mediana=("prev", "median"),
        frac_testavel=("testavel", "mean"),
        abund_mediana=("abund", lambda v: np.median(v[v > 0]) if (v > 0).any() else 0.0)).reset_index()
    G.to_csv(os.path.join(out, "cobertura_generos.tsv"), sep="\t", index=False, float_format="%.4g")

    T = P.pivot_table(index="study", columns="genus", values="testavel", aggfunc="first")
    cu = P.drop_duplicates("study").set_index("study").cultura
    rows = []
    for i, m in mem.items():
        allt = T[m].all(axis=1)
        for c in crop_n.index:
            idx = cu[cu == c].index
            rows.append({"id": i, "membros": " + ".join(m), "cultura": c,
                         "estudos_cultura": len(idx), "estudos_todos_testaveis": int(allt[idx].sum())})
    C = pd.DataFrame(rows)
    C.to_csv(os.path.join(out, "cobertura_consorcios.tsv"), sep="\t", index=False)

    L = ["Estudos e amostras por cultura:", crop_n.to_string(), ""]
    L.append("Ocupacao (fracao de estudos onde o genero aparece) por cultura:")
    L.append(G.pivot(index="genus", columns="cultura", values="ocupacao").to_string(float_format=lambda v: f"{v:.2f}"))
    L += ["", f"Fracao de estudos onde o genero e testavel (prevalencia >= {a.min_prev:.0%}):"]
    L.append(G.pivot(index="genus", columns="cultura", values="frac_testavel").to_string(float_format=lambda v: f"{v:.2f}"))
    L += ["", "Estudos com TODOS os membros testaveis, por consorcio (total de estudos da cultura no cabecalho):"]
    tab = C.pivot(index="id", columns="cultura", values="estudos_todos_testaveis")
    tab.columns = [f"{c} ({crop_n.estudos[c]})" for c in tab.columns]
    L.append(tab.to_string())
    txt = "\n".join(L)
    open(os.path.join(out, "cobertura_resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
