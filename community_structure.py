#!/usr/bin/env python3
# [S4.10] Independent test of community structuring
# Supplementary 4 of the manuscript. See README.md for the full map.
"""
Os consorcios estruturam a comunidade? Teste dentro de cada estudo + combinacao entre estudos.

Para cada estudo e cada consorcio da shortlist com todos os membros testaveis (prevalencia >= 20%):
  x  = log da abundancia relativa somada dos membros (pseudocontagem: metade do menor valor > 0)
  composicao: CLR dos generos testaveis (como em effect_sizes.py), membros excluidos da resposta
     R2 = fracao da variancia da composicao (distancia de Aitchison) explicada por x,
          ajustada pelos fatores de desenho (RDA parcial com um preditor)
  alfa: Shannon da comunidade sem os membros da shortlist; correlacao parcial com x
Nulo: --reps consorcios aleatorios do mesmo tamanho, cada membro trocado por um genero da mesma
  faixa de abundancia (decis no estudo), fora da shortlist.
Por estudo: z = (obs - media nula) / dp nula; percentil do observado no nulo.
Entre estudos: Z de Stouffer (R2: unilateral, obs > nulo; alfa: bilateral), FDR entre consorcios,
  fracao de estudos com obs acima do percentil 95 do nulo, e r de alfa combinado (DerSimonian-Laird).

Saidas: <root>/merged/meta/estruturacao/ (por_estudo.tsv, consorcios.tsv, resumo.txt)
Uso: LD_PRELOAD=... python3 community_structure.py --root ~/syncom --runs runs_subset.tsv --attrs stage2/attrs
"""
import argparse, csv, os, sys
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from effect_sizes import load_genus, load_attrs, detect_design, design_matrix, load_curation, drop_samples
from meta_pairs import meta


def resid(Y, X):
    Z = np.column_stack([np.ones(Y.shape[0]), X])
    b, *_ = np.linalg.lstsq(Z, Y, rcond=None)
    return Y - Z @ b


def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p)
    q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--attrs", default="stage2/attrs")
    ap.add_argument("--design", default=None)
    ap.add_argument("--excluir", default=None)
    ap.add_argument("--reps", type=int, default=199)
    ap.add_argument("--min-prev", type=float, default=0.2)
    ap.add_argument("--max-df-frac", type=float, default=0.25)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, "estruturacao"); os.makedirs(out, exist_ok=True)

    sl = pd.read_csv(os.path.join(mm, "consorcios", "shortlist.tsv"), sep="\t")
    cons = {r.id: [r.centro] + r.parceiros.split(" + ") for _, r in sl.iterrows()}
    U = {g for v in cons.values() for g in v}
    run2bs = {}
    with open(a.runs) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            run2bs[r["run_accession"]] = r["sample_accession"]
    qc = pd.read_csv(os.path.join(md, "study_qc.tsv"), sep="\t", dtype=str)
    qc = qc[qc.status == "incluido"]
    crop = dict(zip(qc.study, qc.cultura))

    forced, excl = load_curation(a.design, a.excluir)
    rows = []
    for i, s in enumerate(sorted(qc.study), 1):
        samples, genera, M = load_genus(os.path.join(md, "genus", f"{s}.tsv.gz"))
        samples, M = drop_samples(samples, M, excl.get(s))
        if len(samples) < 20:
            continue
        n = len(samples)
        named = np.array([not g.startswith("NA") for g in genera])
        prev = (M > 0).mean(0)
        test = named & (prev >= a.min_prev)
        G = [g for g, t in zip(genera, test) if t]
        gi = {g: j for j, g in enumerate(G)}
        todo = {c: m for c, m in cons.items() if all(g in gi for g in m)}
        if not todo:
            continue
        Y = np.column_stack([M[:, test], M[:, ~test].sum(1)]) + 0.5
        L = np.log(Y); clr = (L - L.mean(1, keepdims=True))[:, :-1]
        attrs = load_attrs(os.path.join(a.attrs, f"{s}.tsv"))
        amap = {x: attrs.get(run2bs.get(x, x), {}) for x in samples}
        X, _ = design_matrix(amap, samples, detect_design(amap, samples, forced.get(s)), max(0, int(a.max_df_frac * n)))
        k = X.shape[1]
        if n - 3 - k < 5:
            continue
        Yr = resid(clr, X)
        ss = (Yr ** 2).sum(0)
        rel = M / M.sum(1, keepdims=True)
        # alfa sem os membros da shortlist
        keepU = np.array([g not in U for g in genera])
        q = rel[:, keepU]; q = q / q.sum(1, keepdims=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            H = -(np.where(q > 0, q * np.log(q), 0)).sum(1)
        Hr = resid(H[:, None], X)[:, 0]
        # faixas de abundancia para o nulo
        relT = rel[:, test]
        ab = relT.mean(0)
        cand = np.array([g not in U for g in G])
        cidx = np.where(cand)[0]
        if len(cidx) < 20:
            continue
        edges = np.quantile(ab[cand], np.linspace(0, 1, 11))

        def binof(v):
            return int(np.clip(np.searchsorted(edges, v, side="right") - 1, 0, 9))
        bins = {}
        for j in cidx:
            bins.setdefault(binof(ab[j]), []).append(j)

        for c, mem in todo.items():
            idx = [gi[g] for g in mem]
            sets = [idx]
            for _ in range(a.reps):
                pick = []
                for j in idx:
                    pool = [u for u in bins.get(binof(ab[j]), []) if u not in pick]
                    pick.append(pool[rng.integers(len(pool))] if pool else rng.choice(cidx))
                sets.append(pick)
            xs = []
            for st in sets:
                v = relT[:, st].sum(1)
                eps = v[v > 0].min() / 2 if (v > 0).any() else 1e-6
                xs.append(np.log(v + eps))
            Xs = resid(np.column_stack(xs), X)
            nx2 = (Xs ** 2).sum(0)
            P = Yr.T @ Xs                                    # generos x conjuntos
            num = (P ** 2).sum(0) / nx2
            den = np.full(len(sets), ss.sum())
            for t, st in enumerate(sets):
                num[t] -= (P[st, t] ** 2).sum() / nx2[t]
                den[t] -= ss[st].sum()
            R2 = num / den
            ra = (Hr @ Xs) / np.sqrt((Hr ** 2).sum() * nx2)
            nR, nA = R2[1:], ra[1:]
            rows.append({"study": s, "cultura": crop.get(s, ""), "id": c, "n": n, "k": k,
                         "R2_obs": R2[0], "R2_nulo": nR.mean(),
                         "R2_z": (R2[0] - nR.mean()) / nR.std(ddof=1),
                         "R2_pct": (nR < R2[0]).mean(),
                         "alfa_r": ra[0], "alfa_r_nulo": nA.mean(),
                         "alfa_z": (ra[0] - nA.mean()) / nA.std(ddof=1),
                         "alfa_pct": (nA < ra[0]).mean()})
        if i % 10 == 0:
            print(f"[{i}/{len(qc)}] estudos", file=sys.stderr, flush=True)

    D = pd.DataFrame(rows)
    D.to_csv(os.path.join(out, "por_estudo.tsv"), sep="\t", index=False, float_format="%.4g")
    res = []
    for c, d in D.groupby("id"):
        k = len(d)
        Zr = d.R2_z.sum() / np.sqrt(k); Za = d.alfa_z.sum() / np.sqrt(k)
        mz = pd.DataFrame({"pair": f"{c}|alfa", "z": np.arctanh(np.clip(d.alfa_r, -.999, .999)),
                           "se": 1 / np.sqrt(d.n - 3 - d.k)})
        ma = meta(mz, 1).iloc[0] if k >= 2 else None
        res.append({"id": c, "membros": " + ".join(cons[c]), "estudos": k,
                    "R2_obs_mediano": d.R2_obs.median(), "R2_nulo_mediano": d.R2_nulo.median(),
                    "R2_razao": d.R2_obs.median() / d.R2_nulo.median(),
                    "R2_Z": Zr, "R2_p": stats.norm.sf(Zr), "frac_estudos_R2_p95": (d.R2_pct >= 0.95).mean(),
                    "alfa_r_meta": ma.r_mu if ma is not None else d.alfa_r.iloc[0],
                    "alfa_r_nulo_mediano": d.alfa_r_nulo.median(),
                    "alfa_Z": Za, "alfa_p": 2 * stats.norm.sf(abs(Za)),
                    "frac_estudos_alfa_extremo": ((d.alfa_pct >= 0.975) | (d.alfa_pct <= 0.025)).mean()})
    R = pd.DataFrame(res)
    R["R2_fdr"] = bh(R.R2_p); R["alfa_fdr"] = bh(R.alfa_p)
    R = R.sort_values("R2_Z", ascending=False)
    R.to_csv(os.path.join(out, "consorcios.tsv"), sep="\t", index=False, float_format="%.4g")

    L = [f"Estudos x consorcios testados: {len(D)} | nulo: {a.reps} consorcios pareados por abundancia", "",
         "R2 = variancia da composicao (Aitchison) explicada pela abundancia do consorcio, ajustada pelo desenho;",
         "alfa = correlacao parcial com o Shannon da comunidade restante. Z de Stouffer contra o nulo.", "",
         R[["id", "membros", "estudos", "R2_obs_mediano", "R2_nulo_mediano", "R2_razao", "R2_Z", "R2_fdr",
            "frac_estudos_R2_p95", "alfa_r_meta", "alfa_r_nulo_mediano", "alfa_Z", "alfa_fdr"]]
         .to_string(index=False, float_format=lambda v: f"{v:.3f}")]
    txt = "\n".join(L)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
