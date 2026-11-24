#!/usr/bin/env python3
# [S3.5] Read type
# Supplementary 3 of the manuscript. See README.md for the full map.
"""
Sensibilidade ao tipo de leitura: estudos cuja V4 foi reconstituida por completo (modo PE)
contra estudos em que so a leitura direta foi aproveitada (modo R1).

Refaz a meta-analise separadamente nos dois conjuntos e compara, nos pares testaveis em
ambos: correlacao entre os efeitos, fracao de mesmo sinal e fracao significativa nos dois.
Tambem compara cada conjunto com a rede principal, para dizer quanto da rede sobrevive se
os estudos de fragmento forem descartados.

Saida: <root>/merged/meta/frag/ (comparacao.tsv, resumo.txt)
Uso: LD_PRELOAD=... python3 frag_sens.py --root ~/syncom
"""
import argparse, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meta_pairs import load, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--min-k", type=int, default=5)
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    mm = os.path.join(root, "merged", "meta")
    out = os.path.join(mm, "frag"); os.makedirs(out, exist_ok=True)

    d, qc = load(root, False)
    if "modo" not in qc.columns:
        raise SystemExit("study_qc.tsv sem coluna 'modo'")
    modo = dict(zip(qc.study, qc.modo.fillna("")))
    d["modo"] = d.study.map(modo)
    L = [f"Estudos por modo: " + ", ".join(f"{k} {v}"
         for k, v in qc[qc.status == "incluido"].modo.value_counts().items())]

    M = {}
    for g in ("PE", "R1"):
        sub = d[d.modo == g]
        v = meta(sub, a.min_k)
        if len(v):
            v.index = v.genus_a + "|" + v.genus_b
        M[g] = v
        L.append(f"{g}: {sub.study.nunique()} estudos | {len(v)} pares testaveis (k >= {a.min_k}) | "
                 f"significativos: {int((v.fdr < 0.05).sum()) if len(v) else 0}")

    principal = pd.read_csv(os.path.join(mm, "meta_pairs.tsv.gz"), sep="\t")
    principal.index = principal.genus_a + "|" + principal.genus_b
    entrada = principal[principal.fdr < 0.05]
    L.append(f"Rede principal (FDR < 0,05): {len(entrada)} pares")

    if len(M["PE"]) and len(M["R1"]):
        comum = M["PE"].index.intersection(M["R1"].index)
        C = pd.DataFrame({"r_PE": M["PE"].loc[comum, "r_mu"], "r_R1": M["R1"].loc[comum, "r_mu"],
                          "fdr_PE": M["PE"].loc[comum, "fdr"], "fdr_R1": M["R1"].loc[comum, "fdr"]})
        C["mesmo_sinal"] = np.sign(C.r_PE) == np.sign(C.r_R1)
        C["sig_ambos"] = (C.fdr_PE < 0.05) & (C.fdr_R1 < 0.05)
        C.to_csv(os.path.join(out, "comparacao.tsv"), sep="\t", float_format="%.4g")
        rho = C.r_PE.corr(C.r_R1, method="spearman")
        L += ["", f"Pares testaveis nos dois conjuntos: {len(C)}",
              f"  correlacao de Spearman entre os efeitos: {rho:.3f}",
              f"  mesmo sinal: {C.mesmo_sinal.mean():.1%}",
              f"  significativos nos dois: {C.sig_ambos.mean():.1%}"]
        sig = C[C.sig_ambos]
        if len(sig):
            L.append(f"  entre os significativos nos dois, mesmo sinal: {sig.mesmo_sinal.mean():.1%}")

    # quanto da rede principal se mantem usando so os estudos de V4 completa
    v = M["PE"]
    if len(v):
        alvo = entrada.index.intersection(v.index)
        mant = ((v.loc[alvo, "fdr"] < 0.05) &
                (np.sign(v.loc[alvo, "mu"]) == np.sign(entrada.loc[alvo, "mu"]))).mean()
        L += ["", f"Pares da rede principal testaveis so com V4 completa: {len(alvo)} de {len(entrada)}",
              f"  mantidos com mesmo sinal e FDR < 0,05: {mant:.1%}"]
    txt = "\n".join(L)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
