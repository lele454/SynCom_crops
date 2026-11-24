#!/usr/bin/env python3
# [S3.2] Sensitivity subsets
# Supplementary 3 of the manuscript. See README.md for the full map.
"""
Analises de sensibilidade da meta-analise e teste PGP x demais generos.

1. Sensibilidade: refaz a meta-analise em subconjuntos de estudos e compara com a
   analise principal (todos os estudos):
     sem_R1          sem estudos processados so com a leitura R1
     sem_China       sem estudos da China (~1/3 do acervo)
     com_fatores     so estudos com fatores de desenho ajustados
     sem_fatores     so estudos sem fatores de desenho
   Para cada um: numero de pares por classe; fracao dos consistentes da analise
   principal que continuam consistentes; correlacao (Spearman) dos efeitos combinados.
2. PGP: por genero, heterogeneidade mediana (I2) e fracao de pares so associados,
   em funcao de ser PGP, controlando k medio, ocupacao e abundancia (regressao linear).

Saidas em <root>/merged/meta/sensibilidade/
Uso: LD_PRELOAD=... python3 sensitivity.py --root ~/syncom
"""
import argparse, glob, os, re, sys
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meta_pairs import meta, entrada

log = lambda *x: print(*x, file=sys.stderr, flush=True)
PGPR = re.compile(r"^(Pseudomonas|Bacillus|Paenibacillus|Streptomyces|Bradyrhizobium|Rhizobium|"
                  r"Mesorhizobium|Sinorhizobium|Ensifer|Azospirillum|Herbaspirillum|Burkholderia|"
                  r"Paraburkholderia|Enterobacter|Pantoea|Serratia|Klebsiella|Azotobacter|"
                  r"Gluconacetobacter|Arthrobacter|Microbacterium|Variovorax|Brucella|Stenotrophomonas|"
                  r"Lysobacter|Flavobacterium|Achromobacter|Delftia|Acinetobacter|Kosakonia|"
                  r"Priestia|Peribacillus|Lysinibacillus|Rhodococcus|Agrobacterium)(_[A-Z]+)?$")


def load_all(md, studies):
    parts = []
    for s in studies:
        f = os.path.join(md, "effects", f"pairs_{s}.tsv.gz")
        if os.path.exists(f):
            d = pd.read_csv(f, sep="\t", usecols=["genus_a", "genus_b", "z", "se"])
            d["study"] = s; parts.append(d)
    d = pd.concat(parts, ignore_index=True)
    d["pair"] = d.genus_a + "|" + d.genus_b
    return d


def compare(base, sub, name, n_st):
    b = base.set_index(base.genus_a + "|" + base.genus_b)
    s = sub.set_index(sub.genus_a + "|" + sub.genus_b)
    common = b.index.intersection(s.index)
    rho = stats.spearmanr(b.loc[common, "mu"], s.loc[common, "mu"]).correlation if len(common) > 2 else np.nan
    row = {"analise": name, "estudos": n_st, "pares_testados": len(s)}
    for c in ("consistente_pos", "consistente_neg"):
        bc = b.index[b.classe == c]
        row[c] = int((s.classe == c).sum())
        tested = bc.intersection(s.index)
        same = (s.loc[tested, "classe"] == c).sum() if len(tested) else 0
        sig_same = (s.loc[tested, "fdr"] < 0.05) & (np.sign(s.loc[tested, "mu"]) == (1 if c.endswith("pos") else -1))
        row[f"{c}_mantidos"] = f"{same}/{len(tested)}"
        row[f"{c}_mesmo_sinal_signif"] = f"{int(sig_same.sum())}/{len(tested)}"
    eb = b.index[entrada(b)]; tested = eb.intersection(s.index)
    ok = (s.loc[tested, "fdr"] < 0.05) & (np.sign(s.loc[tested, "mu"]) == np.sign(b.loc[tested, "mu"]))
    row["rede_entrada"] = int(entrada(s).sum())
    row["rede_mantidos_mesmo_sinal_signif"] = f"{int(ok.sum())}/{len(tested)}"
    row["spearman_mu"] = round(rho, 3)
    return row


def ols(y, X, names):
    X = np.column_stack([np.ones(len(y)), X]); names = ["intercepto"] + names
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta; n, p = X.shape
    s2 = res @ res / (n - p)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    t = beta / se; pv = 2 * stats.t.sf(np.abs(t), n - p)
    r2 = 1 - (res @ res) / ((y - y.mean()) @ (y - y.mean()))
    return pd.DataFrame({"coef": beta, "ep": se, "t": t, "p": pv}, index=names), r2, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--min-k", type=int, default=10)
    ap.add_argument("--min-pairs-genus", type=int, default=50)
    a = ap.parse_args()
    md = os.path.join(os.path.abspath(os.path.expanduser(a.root)), "merged")
    od = os.path.join(md, "meta", "sensibilidade"); os.makedirs(od, exist_ok=True)

    qc = pd.read_csv(os.path.join(md, "study_qc.tsv"), sep="\t", dtype=str).fillna("")
    qc = qc[qc.status == "incluido"]
    des = pd.read_csv(os.path.join(md, "effects", "study_design.tsv"), sep="\t")
    com_f = set(des[des.k > 0].study)
    subsets = {
        "principal": set(qc.study),
        "sem_R1": set(qc[qc.modo != "R1"].study),
        "sem_China": set(qc[~qc.country.str.strip().str.startswith("China")].study),
        "com_fatores": set(qc.study) & com_f,
        "sem_fatores": set(qc.study) - com_f,
    }
    d = load_all(md, sorted(subsets["principal"]))
    log(f"Estimativas carregadas: {len(d)}")
    base = meta(d, a.min_k)
    rows = []
    for name, st in subsets.items():
        sub = base if name == "principal" else meta(d[d.study.isin(st)], a.min_k)
        rows.append(compare(base, sub, name, len(st)))
        log(f"  {name}: {len(st)} estudos | pos {rows[-1]['consistente_pos']} "
            f"(mantidos {rows[-1]['consistente_pos_mantidos']}) | neg {rows[-1]['consistente_neg']} "
            f"| Spearman {rows[-1]['spearman_mu']}")
        if name != "principal":
            sub.to_csv(os.path.join(od, f"meta_{name}.tsv.gz"), sep="\t", index=False, float_format="%.5g")
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(od, "sensibilidade_resumo.tsv"), sep="\t", index=False)

    # ---- PGP x demais, controlando prevalencia/poder
    occ = pd.read_csv(os.path.join(md, "genus_occupancy.tsv"), sep="\t").set_index("genus")
    long = pd.concat([base.rename(columns={"genus_a": "g"})[["g", "I2", "classe", "k"]],
                      base.rename(columns={"genus_b": "g"})[["g", "I2", "classe", "k"]]])
    gs = long.groupby("g").agg(n=("I2", "size"), I2=("I2", "median"), k=("k", "mean"),
                               assoc=("classe", lambda x: x.str.startswith("associado").mean()),
                               cons=("classe", lambda x: (x == "consistente_pos").mean()))
    gs = gs[gs.n >= a.min_pairs_genus].join(occ[["ocupacao", "abund_media"]], how="inner")
    gs["pgp"] = [1.0 if PGPR.match(g) else 0.0 for g in gs.index]
    gs["log_k"] = np.log(gs.k); gs["log_ab"] = np.log10(gs.abund_media.astype(float))
    X = gs[["pgp", "log_k", "ocupacao", "log_ab"]].values
    out = []
    for y in ("I2", "assoc", "cons"):
        tb, r2, n = ols(gs[y].values, X, ["PGP", "log(k medio)", "ocupacao", "log10(abundancia)"])
        out.append(f"== {y} ~ PGP + log(k) + ocupacao + log10(abundancia)   (n = {n} generos, R2 = {r2:.3f})")
        out.append(tb.round(5).to_string())
        out.append("")
    txt = "\n".join(out)
    open(os.path.join(od, "pgp_regressao.txt"), "w").write(txt)
    gs.to_csv(os.path.join(od, "pgp_por_genero.tsv"), sep="\t")
    print("\n" + tab.to_string(index=False))
    print("\n" + txt)


if __name__ == "__main__":
    main()
