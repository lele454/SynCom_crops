#!/usr/bin/env python3
# [S2.1] Per-study effect sizes: genus selection, partial correlation, Fisher z
# Supplementary 2 of the manuscript. See README.md for the full map.
"""
Tamanho de efeito por par de generos em cada estudo: correlacao parcial em CLR,
ajustada pelos fatores de desenho do estudo.

Para cada estudo incluido no controle de qualidade:
  1. generos testados = presentes em >= --min-prev das amostras; os demais (e os sem
     genero) sao somados num componente "outros", que entra no CLR mas nao e testado
  2. CLR com pseudocontagem 0,5
  3. fatores de desenho: atributos do BioSample com 2-20 niveis (tratamento, cultivar,
     tempo, local...), detectados automaticamente ou definidos em --design (curadoria)
  4. residuos do CLR contra os fatores (regressao linear) -> correlacao de Pearson entre
     os residuos = correlacao parcial; z de Fisher, EP = 1/sqrt(n - 3 - k)

Saidas (em <root>/merged/effects/):
  pairs_<estudo>.tsv.gz : genus_a, genus_b, r, z, se, n, k
  study_design.tsv      : fatores usados, n, k, generos testados, por estudo

Uso:
  python3 effect_sizes.py --root ~/syncom --runs runs_subset.tsv --attrs stage2/attrs
  (opcional: --design design_curado.tsv  com colunas study, fatores [separados por ;])
"""
import argparse, csv, gzip, os, re, sys
from collections import Counter, defaultdict
import numpy as np

csv.field_size_limit(10 ** 8)
log = lambda *x: print(*x, file=sys.stderr, flush=True)

# prioridade dos tipos de fator quando ha fatores demais para o n de amostras
DESIGN_GROUPS = [
    ("tratamento", r"treat|perturb|chem.?admin|agrochem|watering|fertili[sz]|experimental.?factor|"
                   r"condition|management|regm|irrigat|amend|tillage|rotation|drought|inocul"),
    ("genotipo", r"genotype|cultivar|variety|pedigree|host.?subspecific|line\b"),
    ("tempo", r"collection.?date|^date|time|stage|season|year|month|age\b|days?"),
    ("local", r"geo.?loc|location|site|field|farm|plot|lat.?lon|region|province|country"),
    ("solo", r"soil.?type|soil\b|ph\b|depth"),
]
SKIP = re.compile(r"ena-|ncbi|submission|insdc|replicate|block|^rep\b|sample.?name|title|alias|"
                  r"description|accession|^_|target|primer|seq|library|instrument|platform", re.I)


def load_attrs(path):
    d = defaultdict(dict)
    if os.path.exists(path):
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                d[r["sample"]].setdefault(r["key"], r["value"])
    return d


def detect_design(attrs, samples, forced=None):
    """Retorna lista de (chave, grupo) com 2..20 niveis presentes em >= 80% das amostras."""
    keys = Counter(k for s in samples for k in attrs.get(s, {}))
    out = []
    for k, n in keys.items():
        if n < 0.8 * len(samples) or (SKIP.search(k) and not (forced and k in forced)):
            continue
        if forced is not None and k not in forced:
            continue
        lv = {attrs.get(s, {}).get(k, "") for s in samples}
        if not (2 <= len(lv) <= 20):
            continue
        grp = next((g for g, p in DESIGN_GROUPS if re.search(p, k, re.I)), None)
        if grp or forced is not None:
            out.append((k, grp or "curado"))
    order = {g: i for i, (g, _) in enumerate(DESIGN_GROUPS)}
    return sorted(out, key=lambda x: (order.get(x[1], -1), x[0]))


def design_matrix(attrs, samples, keys, max_df):
    """Dummies (sem a 1a categoria) dos fatores, respeitando o limite de graus de liberdade."""
    cols, used = [], []
    for k, g in keys:
        vals = [attrs.get(s, {}).get(k, "NA") for s in samples]
        levels = sorted(set(vals))[1:]
        if len(cols) + len(levels) > max_df:
            continue
        D = np.array([[1.0 if v == l else 0.0 for l in levels] for v in vals])
        # descarta colunas redundantes com as ja incluidas (fatores aninhados)
        trial = np.column_stack(cols + [D]) if cols else D
        if np.linalg.matrix_rank(np.column_stack([np.ones(len(samples)), trial])) < 1 + trial.shape[1]:
            continue
        cols.append(D); used.append(f"{k}[{g}]:{len(levels) + 1}")
    X = np.column_stack(cols) if cols else np.zeros((len(samples), 0))
    return X, used


def partial_corr(Y, X):
    """Y: amostras x generos (CLR); X: covariaveis. Correlacao entre residuos."""
    n = Y.shape[0]
    Z = np.column_stack([np.ones(n), X])
    beta, *_ = np.linalg.lstsq(Z, Y, rcond=None)
    R = Y - Z @ beta
    R -= R.mean(0)
    sd = np.sqrt((R ** 2).sum(0))
    sd[sd == 0] = np.nan
    C = (R.T @ R) / np.outer(sd, sd)
    return C


def load_curation(design=None, excluir=None):
    """Fatores curados (study -> chaves) e amostras excluidas (study -> runs)."""
    forced, excl = {}, defaultdict(set)
    if design and os.path.exists(design):
        for r in csv.DictReader(open(design), delimiter="\t"):
            forced[r["study"]] = {x.strip() for x in r["fatores"].split(";") if x.strip()}
    if excluir and os.path.exists(excluir):
        for r in csv.DictReader(open(excluir), delimiter="\t"):
            excl[r["study"]].add(r["run"])
    return forced, excl


def drop_samples(samples, M, runs_out):
    if not runs_out:
        return samples, M
    keep = np.array([x not in runs_out for x in samples])
    return [x for x, k in zip(samples, keep) if k], M[keep]


def load_genus(path):
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")[1:]
        samples, rows = [], []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            samples.append(f[0]); rows.append([float(x) for x in f[1:]])
    return samples, hdr, np.array(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--runs", required=True, help="runs_subset.tsv (corrida -> BioSample)")
    ap.add_argument("--attrs", default="stage2/attrs")
    ap.add_argument("--design", default=None, help="TSV curado: study, fatores (sep. ;)")
    ap.add_argument("--min-prev", type=float, default=0.2)
    ap.add_argument("--permute", type=int, default=None,
                    help="modelo nulo: embaralha cada genero entre amostras (semente); "
                         "grava em effects_perm<semente>/")
    ap.add_argument("--excluir", default=None, help="TSV study, run: amostras a excluir")
    ap.add_argument("--sufixo", default="",
                    help="usa merged/genus<sufixo> e study_qc<sufixo>.tsv (ex.: _d5000)")
    ap.add_argument("--outdir", default=None, help="pasta de saida em merged/ (padrao: effects)")
    ap.add_argument("--max-df-frac", type=float, default=0.25,
                    help="no maximo esta fracao de n em graus de liberdade de covariaveis")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged")
    ed = os.path.join(md, a.outdir or ("effects" if a.permute is None else f"effects_perm{a.permute}"))
    os.makedirs(ed, exist_ok=True)
    rng = np.random.default_rng(a.permute) if a.permute is not None else None

    run2bs = {}
    with open(a.runs) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            run2bs[r["run_accession"]] = r["sample_accession"]
    forced = {}
    if a.design:
        for r in csv.DictReader(open(a.design), delimiter="\t"):
            forced[r["study"]] = {x.strip() for x in r["fatores"].split(";") if x.strip()}
    excl = defaultdict(set)
    if a.excluir:
        for r in csv.DictReader(open(a.excluir), delimiter="\t"):
            excl[r["study"]].add(r["run"])
    inc = [r["study"] for r in csv.DictReader(open(os.path.join(md, "study_qc" + a.sufixo + ".tsv")), delimiter="\t")
           if r["status"] == "incluido"]
    log(f"Estudos incluidos: {len(inc)}")

    rows_design = []
    for i, s in enumerate(sorted(inc), 1):
        out = os.path.join(ed, f"pairs_{s}.tsv.gz")
        samples, genera, M = load_genus(os.path.join(md, "genus" + a.sufixo, f"{s}.tsv.gz"))
        if excl.get(s):
            keep = np.array([x not in excl[s] for x in samples])
            samples = [x for x, k_ in zip(samples, keep) if k_]; M = M[keep]
            if len(samples) < 20:
                log(f"[{i}/{len(inc)}] {s}: {len(samples)} amostras apos exclusao (< 20), estudo pulado")
                continue
        n = len(samples)
        if rng is not None:                    # nulo: quebra as associacoes, mantem marginais
            M = np.column_stack([rng.permutation(M[:, j]) for j in range(M.shape[1])])
        named = np.array([not g.startswith("NA") for g in genera])
        prev = (M > 0).mean(0)
        test = named & (prev >= a.min_prev)
        G = [g for g, t in zip(genera, test) if t]
        Y = np.column_stack([M[:, test], M[:, ~test].sum(1)]) + 0.5          # + "outros"
        L = np.log(Y); clr = L - L.mean(1, keepdims=True)
        clr = clr[:, :-1]                                                     # nao testa "outros"
        # fatores de desenho (atributos do BioSample de cada amostra)
        attrs = load_attrs(os.path.join(a.attrs, f"{s}.tsv"))
        bs = [run2bs.get(x, x) for x in samples]
        amap = {x: attrs.get(b, {}) for x, b in zip(samples, bs)}
        keys = detect_design(amap, samples, forced.get(s))
        X, used = design_matrix(amap, samples, keys, max(0, int(a.max_df_frac * n)))
        k = X.shape[1]
        C = partial_corr(clr, X)
        dfree = n - 3 - k
        iu = np.triu_indices(len(G), 1)
        r = np.clip(C[iu], -0.999999, 0.999999)
        z = np.arctanh(r); se = 1 / np.sqrt(dfree)
        ab = (M / M.sum(1, keepdims=True))[:, test].mean(0)
        with gzip.open(out, "wt") as fo:
            fo.write("genus_a\tgenus_b\tr\tz\tse\tn\tk\tprev_a\tprev_b\tab_a\tab_b\n")
            pa, pb = prev[test], prev[test]
            for a_, b_, rr, zz in zip(iu[0], iu[1], r, z):
                if np.isnan(rr):
                    continue
                ga, gb = sorted((G[a_], G[b_]))
                ia, ib = (a_, b_) if G[a_] == ga else (b_, a_)
                fo.write(f"{ga}\t{gb}\t{rr:.4f}\t{zz:.4f}\t{se:.4f}\t{n}\t{k}\t"
                         f"{pa[ia]:.3f}\t{pb[ib]:.3f}\t{ab[ia]:.2e}\t{ab[ib]:.2e}\n")
        rows_design.append({"study": s, "n": n, "k": k, "generos_testados": len(G),
                            "pares": len(r), "fatores": "; ".join(used) or "nenhum",
                            "fonte": "curado" if s in forced else "automatico"})
        log(f"[{i}/{len(inc)}] {s}: n={n} k={k} generos={len(G)} | {'; '.join(used) or 'sem fatores'}")
    with open(os.path.join(ed, "study_design.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows_design[0]), delimiter="\t"); w.writeheader(); w.writerows(rows_design)
    nf = sum(1 for r in rows_design if r["k"] > 0)
    log(f"\nEstudos com fatores de desenho: {nf}/{len(rows_design)} | "
        f"pares (estudo x par): {sum(r['pares'] for r in rows_design)}")


if __name__ == "__main__":
    main()
