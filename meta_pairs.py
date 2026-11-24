#!/usr/bin/env python3
# [S2.2] Random-effects pooling, both intervals, FDR and the null model
# Supplementary 2 of the manuscript. See README.md for the full map.
"""
Meta-analise de efeitos aleatorios por par de generos (z de Fisher), vetorizada.

Para cada par presente em >= --min-k estudos:
  - tau^2 por DerSimonian-Laird; efeito combinado, EP, IC 95%, p (teste z)
  - intervalo de predicao 95% (t com k-2 gl): onde se espera o efeito num estudo novo
  - I^2, Q, fracao de estudos com o mesmo sinal do efeito combinado
  - FDR (Benjamini-Hochberg) sobre todos os pares testados
Classe:
  consistente_pos / consistente_neg : FDR < 0,05 e intervalo de predicao exclui zero
  associado_pos / associado_neg     : FDR < 0,05, mas o intervalo de predicao cruza zero
  ns                                : FDR >= 0,05

Saidas (em <root>/merged/meta/):
  meta_pairs.tsv.gz   todos os pares testados
  consistentes.tsv    pares consistentes (base da rede de consenso)
  resumo.txt

Uso: python3 meta_pairs.py --root ~/syncom [--min-k 10] [--excluir-r1]
"""
import argparse, glob, os, sys
import numpy as np, pandas as pd
from scipy import stats

log = lambda *x: print(*x, file=sys.stderr, flush=True)


def entrada(S, criterio=None):
    """Pares que entram na rede. fdr (padrao): FDR < 0,05; consistente: intervalo de predicao exclui zero.
    Definido pela variavel de ambiente SYNCOM_ENTRADA (run_all.sh)."""
    criterio = criterio or os.environ.get("SYNCOM_ENTRADA", "fdr")
    return (S.fdr < 0.05) if criterio == "fdr" else S.classe.str.startswith("consistente")


def load(root, excluir_r1, effects="effects"):
    md = os.path.join(root, "merged")
    qc = pd.read_csv(os.path.join(md, "study_qc.tsv"), sep="\t", dtype=str)
    qc = qc[qc.status == "incluido"]
    if excluir_r1:
        qc = qc[qc.modo != "R1"]
    keep = set(qc.study)
    parts = []
    for f in sorted(glob.glob(os.path.join(md, effects, "pairs_*.tsv.gz"))):
        s = os.path.basename(f)[6:-7]
        if s not in keep:
            continue
        d = pd.read_csv(f, sep="\t", usecols=["genus_a", "genus_b", "z", "se"])
        d["study"] = s
        parts.append(d)
    d = pd.concat(parts, ignore_index=True)
    d["pair"] = d.genus_a + "|" + d.genus_b
    return d, qc


def meta(d, min_k):
    d = d[np.isfinite(d.z) & (d.se > 0)].copy()
    k = d.groupby("pair").z.transform("size")
    d = d[k >= min_k].copy()
    if d.empty:
        return pd.DataFrame(columns=["genus_a", "genus_b", "k", "mu", "se_mu", "ci_lo", "ci_hi", "pi_lo", "pi_hi",
                                     "p", "fdr", "tau2", "I2", "Q", "frac_mesmo_sinal", "r_mu", "r_ci_lo",
                                     "r_ci_hi", "r_pi_lo", "r_pi_hi", "classe"])
    d["v"] = d.se ** 2
    d["w"] = 1 / d.v
    d["wz"] = d.w * d.z
    d["wz2"] = d.w * d.z ** 2
    d["w2"] = d.w ** 2
    g = d.groupby("pair")
    S = g[["w", "wz", "wz2", "w2"]].sum()
    S["k"] = g.size()
    S["Q"] = S.wz2 - S.wz ** 2 / S.w
    C = S.w - S.w2 / S.w
    S["tau2"] = np.maximum(0, (S.Q - (S.k - 1)) / C)
    d = d.join(S["tau2"], on="pair")
    d["ws"] = 1 / (d.v + d.tau2)
    d["wsz"] = d.ws * d.z
    g = d.groupby("pair")
    S2 = g[["ws", "wsz"]].sum()
    S["mu"] = S2.wsz / S2.ws
    S["se_mu"] = np.sqrt(1 / S2.ws)
    d = d.join(S["mu"], on="pair")
    S["frac_mesmo_sinal"] = (np.sign(d.z) == np.sign(d.mu)).groupby(d.pair).mean()
    S["ci_lo"] = S.mu - 1.96 * S.se_mu
    S["ci_hi"] = S.mu + 1.96 * S.se_mu
    tq = stats.t.ppf(0.975, np.maximum(S.k - 2, 1))
    half = tq * np.sqrt(S.tau2 + S.se_mu ** 2)
    S["pi_lo"] = S.mu - half
    S["pi_hi"] = S.mu + half
    S["p"] = 2 * stats.norm.sf(np.abs(S.mu / S.se_mu))
    S["I2"] = np.where(S.Q > 0, np.maximum(0, (S.Q - (S.k - 1)) / S.Q), 0)
    # FDR (Benjamini-Hochberg)
    p = S.p.values; o = np.argsort(p); m = len(p)
    q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    S["fdr"] = np.minimum(q, 1)
    for c in ("mu", "ci_lo", "ci_hi", "pi_lo", "pi_hi"):
        S["r_" + c] = np.tanh(S[c])
    sig = S.fdr < 0.05
    S["classe"] = "ns"
    S.loc[sig & (S.mu > 0), "classe"] = "associado_pos"
    S.loc[sig & (S.mu < 0), "classe"] = "associado_neg"
    S.loc[sig & (S.pi_lo > 0), "classe"] = "consistente_pos"
    S.loc[sig & (S.pi_hi < 0), "classe"] = "consistente_neg"
    S = S.drop(columns=["w", "wz", "wz2", "w2"]).reset_index()
    S[["genus_a", "genus_b"]] = S.pair.str.split("|", n=1, expand=True)
    cols = ["genus_a", "genus_b", "k", "mu", "se_mu", "ci_lo", "ci_hi", "pi_lo", "pi_hi", "p", "fdr",
            "tau2", "I2", "Q", "frac_mesmo_sinal", "r_mu", "r_ci_lo", "r_ci_hi", "r_pi_lo", "r_pi_hi",
            "classe"]
    return S[cols].sort_values(["classe", "mu"], ascending=[True, False])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--min-k", type=int, default=10)
    ap.add_argument("--excluir-r1", action="store_true", help="sensibilidade: sem estudos so-R1")
    ap.add_argument("--tag", default="", help="sufixo dos arquivos de saida (sensibilidades)")
    ap.add_argument("--effects", default="effects", help="pasta de efeitos (ex.: effects_perm1)")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    od = os.path.join(root, "merged", "meta"); os.makedirs(od, exist_ok=True)
    d, qc = load(root, a.excluir_r1, a.effects)
    log(f"Estudos: {d.study.nunique()} | estimativas estudo x par: {len(d)} | pares distintos: {d.pair.nunique()}")
    S = meta(d, a.min_k)
    t = a.tag and f"_{a.tag}"
    S.to_csv(os.path.join(od, f"meta_pairs{t}.tsv.gz"), sep="\t", index=False, float_format="%.5g")
    cons = S[S.classe.str.startswith("consistente")]
    cons.to_csv(os.path.join(od, f"consistentes{t}.tsv"), sep="\t", index=False, float_format="%.4g")
    rede = S[entrada(S)]
    rede.to_csv(os.path.join(od, f"rede_entrada{t}.tsv"), sep="\t", index=False, float_format="%.4g")
    vc = S.classe.value_counts()
    gen = pd.unique(cons[["genus_a", "genus_b"]].values.ravel())
    lines = [f"Pares testados (k >= {a.min_k}): {len(S)}",
             *(f"  {c}: {vc.get(c, 0)}" for c in
               ("consistente_pos", "consistente_neg", "associado_pos", "associado_neg", "ns")),
             f"Generos na rede de consenso: {len(gen)}",
             f"Rede de entrada ({os.environ.get('SYNCOM_ENTRADA', 'fdr')}): {len(rede)} pares "
             f"(pos {(rede.mu > 0).sum()}, neg {(rede.mu < 0).sum()})",
             f"k mediano dos pares testados: {S.k.median():.0f}",
             f"I2 mediano dos consistentes: {cons.I2.median():.2f}" if len(cons) else ""]
    open(os.path.join(od, f"resumo{t}.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if len(cons):
        print("\nTop 10 positivos consistentes:")
        top = cons[cons.classe == "consistente_pos"].nlargest(10, "mu")
        for _, r in top.iterrows():
            print(f"  {r.genus_a} - {r.genus_b}: r={r.r_mu:.2f} [PI {r.r_pi_lo:.2f}; {r.r_pi_hi:.2f}] k={r.k}")


if __name__ == "__main__":
    main()
