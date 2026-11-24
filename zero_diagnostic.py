#!/usr/bin/env python3
# [S3.3] Zero diagnostic
# Supplementary 3 of the manuscript. See README.md for the full map.
"""
Diagnostico de zeros: a associacao dos pares do nucleo vem de coabundancia ou so de
coocorrencia (presenca/ausencia conjunta)?

Para cada par consistente (nucleo_robusto.tsv, todos os tiers) e cada estudo onde e testavel:
  A. correlacao condicional: mesmo CLR e mesmos fatores de desenho de effect_sizes.py,
     residuos calculados em todas as amostras; correlacao de Pearson so nas amostras em
     que os DOIS generos foram detectados (n_ambos >= --min-n). EP = 1/sqrt(n_ambos-3-k)
  B. coocorrencia: log odds ratio da deteccao conjunta (tabela 2x2, correcao de Haldane
     +0,5, EP de Woolf); pulado se um dos generos esta em todas as amostras
Meta-analise DerSimonian-Laird (funcao de meta_pairs.py) para A e B, FDR entre os pares do nucleo.

Classe do diagnostico:
  coabundancia      : A significativo (FDR < 0,05) com o mesmo sinal da analise principal
  coabund_fraca     : A com o mesmo sinal, nao significativo, e B nao significativo
  so_coocorrencia   : A ns ou de sinal oposto, B significativo com o mesmo sinal
  nao_confirmado    : A e B nao confirmam o sinal
  inconclusivo      : A testavel em < --min-k estudos

Saidas: <root>/merged/meta/diagnostico_zeros/diag_pares.tsv e resumo.txt
Uso: LD_PRELOAD=... python3 zero_diagnostic.py --root ~/syncom --runs runs_subset.tsv --attrs stage2/attrs
"""
import argparse, csv, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from effect_sizes import load_genus, load_attrs, detect_design, design_matrix, load_curation, drop_samples
from meta_pairs import meta

log = lambda *x: print(*x, file=sys.stderr, flush=True)


def residuals(Y, X):
    n = Y.shape[0]
    Z = np.column_stack([np.ones(n), X])
    beta, *_ = np.linalg.lstsq(Z, Y, rcond=None)
    return Y - Z @ beta


def conditional_effects(md, inc, want, run2bs, attrs_dir, min_prev=0.2, max_df_frac=0.25, min_n=10,
                        forced=None, excl=None):
    """Correlacao condicional (amostras com os dois generos) e log OR de coocorrencia por estudo."""
    gen_core = {g for p in want for g in p.split("|")}
    A, B = [], []
    for i, s in enumerate(inc, 1):
        samples, genera, M = load_genus(os.path.join(md, "genus", f"{s}.tsv.gz"))
        samples, M = drop_samples(samples, M, (excl or {}).get(s))
        if len(samples) < 20:
            continue
        n = len(samples)
        named = np.array([not g.startswith("NA") for g in genera])
        prev = (M > 0).mean(0)
        test = named & (prev >= min_prev)
        G = [g for g, t in zip(genera, test) if t]
        Y = np.column_stack([M[:, test], M[:, ~test].sum(1)]) + 0.5
        L = np.log(Y); clr = (L - L.mean(1, keepdims=True))[:, :-1]
        attrs = load_attrs(os.path.join(attrs_dir, f"{s}.tsv"))
        amap = {x: attrs.get(run2bs.get(x, x), {}) for x in samples}
        keys = detect_design(amap, samples, (forced or {}).get(s))
        X, _ = design_matrix(amap, samples, keys, max(0, int(max_df_frac * n)))
        k = X.shape[1]
        R = residuals(clr, X)
        P = M[:, test] > 0
        idx = {g: j for j, g in enumerate(G) if g in gen_core}
        names = sorted(idx)
        for x in range(len(names)):
            for y in range(x + 1, len(names)):
                ga, gb = names[x], names[y]
                pr = f"{ga}|{gb}"
                if pr not in want:
                    continue
                ia, ib = idx[ga], idx[gb]
                pa, pb = P[:, ia], P[:, ib]
                both = pa & pb
                nb = int(both.sum())
                if nb >= min_n and nb - 3 - k >= 3:
                    ra, rb = R[both, ia], R[both, ib]
                    ra = ra - ra.mean(); rb = rb - rb.mean()
                    den = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
                    if den > 0:
                        r = np.clip((ra * rb).sum() / den, -0.999999, 0.999999)
                        A.append((pr, s, np.arctanh(r), 1 / np.sqrt(nb - 3 - k), nb, n))
                if pa.all() or pb.all():
                    continue
                n11 = nb + .5; n10 = int((pa & ~pb).sum()) + .5
                n01 = int((~pa & pb).sum()) + .5; n00 = int((~pa & ~pb).sum()) + .5
                B.append((pr, s, np.log(n11 * n00 / (n10 * n01)),
                          np.sqrt(1 / n11 + 1 / n10 + 1 / n01 + 1 / n00)))
        if i % 10 == 0 or i == len(inc):
            log(f"[{i}/{len(inc)}] estudos lidos")
    return (pd.DataFrame(A, columns=["pair", "study", "z", "se", "n_ambos", "n"]),
            pd.DataFrame(B, columns=["pair", "study", "z", "se"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--attrs", default="stage2/attrs")
    ap.add_argument("--design", default=None, help="design_curado.tsv")
    ap.add_argument("--excluir", default=None, help="excluir_amostras.tsv")
    ap.add_argument("--min-prev", type=float, default=0.2)
    ap.add_argument("--max-df-frac", type=float, default=0.25)
    ap.add_argument("--min-n", type=int, default=10, help="minimo de amostras com os dois generos")
    ap.add_argument("--min-k", type=int, default=5, help="minimo de estudos na meta condicional")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, "diagnostico_zeros"); os.makedirs(out, exist_ok=True)

    core = pd.read_csv(os.path.join(mm, "nucleo_robusto.tsv"), sep="\t")
    core["pair"] = core.genus_a + "|" + core.genus_b
    want = set(core.pair)
    gen_core = set(core.genus_a) | set(core.genus_b)

    run2bs = {}
    with open(a.runs) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            run2bs[r["run_accession"]] = r["sample_accession"]
    qc = pd.read_csv(os.path.join(md, "study_qc.tsv"), sep="\t", dtype=str)
    inc = sorted(qc.study[qc.status == "incluido"])

    forced, excl = load_curation(a.design, a.excluir)
    dA, dB = conditional_effects(md, inc, want, run2bs, a.attrs, a.min_prev, a.max_df_frac, a.min_n,
                                 forced, excl)

    mA = meta(dA, a.min_k)
    mA.index = mA.genus_a + "|" + mA.genus_b
    mB = meta(dB, a.min_k)
    mB.index = mB.genus_a + "|" + mB.genus_b
    nb_med = dA.groupby("pair").n_ambos.median()
    frac_amostras = (dA.n_ambos / dA.n).groupby(dA.pair).median()

    d = core[["tier", "sinal", "genus_a", "genus_b", "pgp", "k", "r_mu", "pair"]].copy()
    d = d.join(mA[["k", "r_mu", "r_ci_lo", "r_ci_hi", "fdr", "I2"]].add_suffix("_cond"), on="pair")
    d = d.join(mB[["k", "mu", "fdr"]].rename(columns={"mu": "logOR"}).add_suffix("_cooc"), on="pair")
    d["n_ambos_mediano"] = d.pair.map(nb_med)
    d["frac_amostras_ambos"] = d.pair.map(frac_amostras)
    sg = np.sign(d.r_mu)
    a_ok = (np.sign(d.r_mu_cond) == sg) & (d.fdr_cond < 0.05)
    a_same = np.sign(d.r_mu_cond) == sg
    b_ok = (np.sign(d.logOR_cooc) == sg) & (d.fdr_cooc < 0.05)
    d["diag"] = np.select(
        [d.k_cond.isna(), a_ok, a_same & ~b_ok, b_ok],
        ["inconclusivo", "coabundancia", "coabund_fraca", "so_coocorrencia"], "nao_confirmado")
    d["atenuacao"] = d.r_mu_cond / d.r_mu
    d = d.drop(columns=["pair"])
    d.to_csv(os.path.join(out, "diag_pares.tsv"), sep="\t", index=False, float_format="%.4g")

    L = ["Diagnostico de zeros nos pares consistentes (todos os tiers)",
         f"Pares com correlacao condicional testavel (k >= {a.min_k}): {int(d.k_cond.notna().sum())}/{len(d)}", ""]
    for t in (1, 2, 3):
        sub = d[d.tier == t]
        if not len(sub):
            continue
        L.append(f"Tier {t}:")
        L.append(pd.crosstab(sub.sinal, sub.diag, margins=True, margins_name="total").to_string())
        ok = sub[sub.k_cond.notna()]
        if len(ok):
            L.append(f"  atenuacao mediana (r_cond / r_principal): {ok.atenuacao.median():.2f} | "
                     f"amostras com os dois detectados (mediana): {ok.frac_amostras_ambos.median():.0%}")
        L.append("")
    t1 = d[(d.tier == 1) & (d.sinal == "positivo")].sort_values("r_mu", ascending=False)
    L.append("Tier 1 positivos:")
    for _, r in t1.iterrows():
        rc = "  -  " if pd.isna(r.r_mu_cond) else f"{r.r_mu_cond:.2f}"
        L.append(f"  {r.genus_a} - {r.genus_b}: r={r.r_mu:.2f} r_cond={rc} k_cond="
                 f"{'' if pd.isna(r.k_cond) else int(r.k_cond)} -> {r['diag']}{'  (PGP)' if r.pgp else ''}")
    txt = "\n".join(L)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
