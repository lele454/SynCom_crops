#!/usr/bin/env python3
# [S3.10] Entry criterion tested against replication
# Supplementary 3 of the manuscript. See README.md for the full map.
"""
Teste justo da replicacao temporal: o criterio "consistente" (intervalo de predicao exclui zero)
prediz replicacao melhor que "associado" (so FDR) quando tamanho de efeito e poder sao iguais?

Le validacao_temporal/validacao.tsv (temporal_validation.py) e faz:
  1. Regressao logistica (IRLS): replicado ~ consistente + |z| descoberta + |r| descoberta
       + log k descoberta + log k validacao + sinal negativo   (so pares consistente/associado)
  2. Comparacao pareada: estratos por quintis de |r| descoberta x tercis de k validacao;
       diferenca media de replicacao consistente - associado dentro dos estratos (Mantel-Haenszel)
  3. Heterogeneidade como preditor continuo (todos os pares, inclusive ns):
       replicado ~ |z| + I2 + log k validacao
  4. Maldicao do vencedor: inclinacao de r validacao ~ r descoberta por classe
Desfecho principal: replicado (mesmo sinal e p unilateral < 0,05); secundario: mesmo sinal.

Saida: <root>/merged/meta/validacao_temporal/teste_replicacao.txt
Uso: python3 replication_test.py --root ~/syncom
"""
import argparse, os
import numpy as np, pandas as pd
from scipy import stats


def logit(y, X, names):
    X = np.column_stack([np.ones(len(y)), X]); names = ["intercepto"] + names
    b = np.zeros(X.shape[1])
    for _ in range(50):
        p = 1 / (1 + np.exp(-X @ b)); W = p * (1 - p) + 1e-9
        H = X.T @ (X * W[:, None]); g = X.T @ (y - p)
        step = np.linalg.solve(H, g); b += step
        if np.abs(step).max() < 1e-8:
            break
    se = np.sqrt(np.diag(np.linalg.inv(H)))
    z = b / se
    return pd.DataFrame({"coef": b, "OR": np.exp(b), "EP": se, "z": z, "p": 2 * stats.norm.sf(np.abs(z))}, index=names)


def mh(d, out):
    """Mantel-Haenszel: razao de chances comum consistente vs associado e diferenca ponderada."""
    num = den = 0.0; diffs = []; ws = []
    for _, g in d.groupby("estrato"):
        c = g[g.grupo == "consistente"]; a = g[g.grupo == "associado"]
        if len(c) == 0 or len(a) == 0:
            continue
        n = len(g)
        A = c[out].sum(); B = len(c) - A; C = a[out].sum(); D = len(a) - C
        num += A * D / n; den += B * C / n
        w = len(c) * len(a) / n
        diffs.append(c[out].mean() - a[out].mean()); ws.append(w)
    return (num / den if den > 0 else np.nan), np.average(diffs, weights=ws), len(diffs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    od = os.path.join(os.path.abspath(os.path.expanduser(a.root)), "merged", "meta", "validacao_temporal")
    d = pd.read_csv(os.path.join(od, "validacao.tsv"), sep="\t", index_col=0)
    d["absz"] = (d.mu / d.se_mu).abs(); d["absr"] = d.r_mu.abs() * 10                 # por 0,1 de r
    d["logk"] = np.log(d.k); d["logkv"] = np.log(d.k_val); d["neg"] = (d.mu < 0).astype(float)
    d["repl"] = d.replicado.astype(float); d["sinal"] = d.mesmo_sinal.astype(float)
    L = [f"Pares: {len(d)} (" + ", ".join(f"{k} {v}" for k, v in d.grupo.value_counts().items()) + ")", ""]

    s = d[d.grupo.isin(["consistente", "associado"])].copy()
    s["cons"] = (s.grupo == "consistente").astype(float)
    cols = ["cons", "absz", "absr", "logk", "logkv", "neg"]
    for out, lab in (("repl", "replicado"), ("sinal", "mesmo sinal")):
        L.append(f"1. Logistica ({lab}) ~ consistente + |z| + |r| (por 0,1) + log k + log k_val + negativo:")
        L.append(logit(s[out].values, s[cols].values, cols).to_string(float_format=lambda v: f"{v:.3g}"))
        L.append("")

    s["estrato"] = (pd.qcut(s.absr.rank(method="first"), 5, labels=False).astype(str) + "_" +
                    pd.qcut(s.k_val.rank(method="first"), 3, labels=False).astype(str))
    L.append("2. Comparacao pareada por estratos (|r| descoberta x k validacao):")
    for out, lab in (("repl", "replicado"), ("sinal", "mesmo sinal")):
        orr, dif, ns = mh(s, out)
        L.append(f"   {lab}: OR de Mantel-Haenszel (consistente vs associado) = {orr:.2f} | "
                 f"diferenca media = {dif:+.3f} | estratos = {ns}")
    L.append("")

    cols3 = ["absz", "I2", "logkv"]
    L.append("3. Heterogeneidade como preditor (todos os pares) ~ |z| + I2 + log k_val:")
    for out, lab in (("repl", "replicado"), ("sinal", "mesmo sinal")):
        L.append(f"   {lab}:")
        L.append(logit(d[out].values, d[cols3].values, cols3).to_string(float_format=lambda v: f"{v:.3g}"))
    L.append("")

    L.append("4. Maldicao do vencedor: inclinacao de r_validacao ~ r_descoberta (sem intercepto):")
    for g, x in d.groupby("grupo"):
        sl = (x.r_mu * x.r_mu_val).sum() / (x.r_mu ** 2).sum()
        L.append(f"   {g}: inclinacao = {sl:.2f} (1 = sem atenuacao) | n = {len(x)}")
    txt = "\n".join(L)
    open(os.path.join(od, "teste_replicacao.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
