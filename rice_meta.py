#!/usr/bin/env python3
# [S3.8] Rice-only network
# Supplementary 3 of the manuscript. See README.md for the full map.
"""
Consorcios especificos de arroz (rizosfera alagada).

1. Meta-analise DL so com os estudos de arroz (k >= --min-k)
2. Robustez por exclusao de um estudo por vez (leave-one-out), no lugar dos tiers:
     tier 1 = em todas as exclusoes o par mantem o sinal e p < 0,05; tier 2 = demais consistentes
3. Diagnostico de zeros (correlacao condicional + coocorrencia), como em zero_diagnostic.py
4. Estrelas (star_consortia.py) com copresenca so nos estudos de arroz

Saidas: <root>/merged/meta/arroz/ (meta_pairs_arroz.tsv.gz, diag_pares_arroz.tsv, resumo_arroz.txt)
        e <root>/merged/meta/consorcios/estrelas_arroz*.tsv/txt
Uso: LD_PRELOAD=... python3 rice_meta.py --root ~/syncom --runs runs_subset.tsv --attrs stage2/attrs
"""
import argparse, csv, os, subprocess, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from meta_pairs import load, meta, entrada
from zero_diagnostic import conditional_effects
from effect_sizes import load_curation


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--attrs", default="stage2/attrs")
    ap.add_argument("--design", default=None)
    ap.add_argument("--excluir", default=None)
    ap.add_argument("--cultura", default="arroz")
    ap.add_argument("--min-k", type=int, default=5)
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, a.cultura.replace("/", "_")); os.makedirs(out, exist_ok=True)

    d, qc = load(root, False)
    qc = qc[qc.cultura == a.cultura]
    studies = sorted(qc.study)
    d = d[d.study.isin(studies)]
    L = [f"Cultura: {a.cultura} | estudos: {len(studies)}",
         "Paises: " + ", ".join(f"{k} ({v})" for k, v in qc.country.fillna("NA").value_counts().items())]

    S = meta(d, a.min_k)
    S.to_csv(os.path.join(out, f"meta_pairs_{a.cultura}.tsv.gz"), sep="\t", index=False, float_format="%.5g")
    vc = S.classe.value_counts()
    L.append(f"Pares testados (k >= {a.min_k}): {len(S)} | " +
             " | ".join(f"{c}: {vc.get(c, 0)}" for c in ("consistente_pos", "consistente_neg",
                                                          "associado_pos", "associado_neg")))
    cons = S[entrada(S)].copy()
    cons["pair"] = cons.genus_a + "|" + cons.genus_b
    if cons.empty:
        L.append("Nenhum par consistente."); print("\n".join(L)); return

    # leave-one-out
    dc = d[d.pair.isin(set(cons.pair))]
    ok = pd.Series(True, index=cons.pair.values)
    for s in studies:
        m = meta(dc[dc.study != s], max(a.min_k - 1, 3))
        m.index = m.genus_a + "|" + m.genus_b
        m = m.reindex(ok.index)
        sign = np.sign(cons.set_index("pair").mu.reindex(ok.index))
        ok &= (np.sign(m.mu) == sign) & (m.p < 0.05)
    cons["tier"] = np.where(cons.pair.map(ok).fillna(False), 1, 2)
    cons["sinal"] = np.where(cons.mu > 0, "positivo", "negativo")

    # diagnostico de zeros
    run2bs = {}
    with open(a.runs) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            run2bs[r["run_accession"]] = r["sample_accession"]
    forced, excl = load_curation(a.design, a.excluir)
    dA, dB = conditional_effects(md, studies, set(cons.pair), run2bs, a.attrs, forced=forced, excl=excl)
    mA = meta(dA, max(a.min_k - 1, 3)); mA.index = mA.genus_a + "|" + mA.genus_b
    mB = meta(dB, max(a.min_k - 1, 3)); mB.index = mB.genus_a + "|" + mB.genus_b
    cons = cons.join(mA[["k", "r_mu", "fdr"]].add_suffix("_cond"), on="pair")
    cons = cons.join(mB[["mu", "fdr"]].rename(columns={"mu": "logOR"}).add_suffix("_cooc"), on="pair")
    sg = np.sign(cons.r_mu)
    a_ok = (np.sign(cons.r_mu_cond) == sg) & (cons.fdr_cond < 0.05)
    a_same = np.sign(cons.r_mu_cond) == sg
    b_ok = (np.sign(cons.logOR_cooc) == sg) & (cons.fdr_cooc < 0.05)
    cons["diag"] = np.select([cons.k_cond.isna(), a_ok, a_same & ~b_ok, b_ok],
                             ["inconclusivo", "coabundancia", "coabund_fraca", "so_coocorrencia"],
                             "nao_confirmado")
    cols = ["tier", "sinal", "genus_a", "genus_b", "k", "r_mu", "r_pi_lo", "r_pi_hi", "I2", "fdr",
            "k_cond", "r_mu_cond", "fdr_cond", "logOR_cooc", "fdr_cooc", "diag"]
    cons = cons.sort_values(["tier", "sinal", "r_mu"], ascending=[True, False, False])[cols]
    fdiag = os.path.join(out, f"diag_pares_{a.cultura}.tsv")
    cons.to_csv(fdiag, sep="\t", index=False, float_format="%.4g")

    L += ["", "Consistentes por tier (leave-one-out) e diagnostico:",
          pd.crosstab([cons.tier, cons.sinal], cons["diag"], margins=True, margins_name="total").to_string(), "",
          "Positivos tier 1 (top 20):"]
    for _, r in cons[(cons.tier == 1) & (cons.sinal == "positivo")].head(20).iterrows():
        rc = "-" if pd.isna(r.r_mu_cond) else f"{r.r_mu_cond:.2f}"
        L.append(f"  {r.genus_a} - {r.genus_b}: r={r.r_mu:.2f} k={r.k} r_cond={rc} -> {r['diag']}")
    txt = "\n".join(L)
    open(os.path.join(out, f"resumo_{a.cultura}.txt"), "w").write(txt)
    print(txt)

    print("\n================ Estrelas (" + a.cultura + ") ================", flush=True)
    for extra in ([], ["--incluir-cooc"]):
        tag = "estrelas_" + a.cultura + ("_cooc" if extra else "")
        print(f"\n--- {'coabundancia + coocorrencia' if extra else 'so coabundancia'} ---", flush=True)
        subprocess.run([sys.executable, os.path.join(HERE, "star_consortia.py"), "--root", root,
                        "--diag", fdiag, "--meta", os.path.join(out, f"meta_pairs_{a.cultura}.tsv.gz"),
                        "--cultura", a.cultura, "--prefixo", tag, "--max-size", "4"] + extra)


if __name__ == "__main__":
    main()
