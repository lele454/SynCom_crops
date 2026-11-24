#!/usr/bin/env python3
# [S1.6] Stratified subsampling of runs
# Supplementary 1 of the manuscript. See README.md for the full map.
"""
Gera a lista de download do recorte "grandes culturas", com no maximo N amostras por
estudo, escolhidas por amostragem estratificada pelo desenho experimental.

Estratos = combinacao dos valores dos atributos de tratamento/genotipo/tempo de cada
amostra (de stage2/attrs/<estudo>.tsv). Alocacao proporcional, com pelo menos 1
amostra por estrato; dentro do estrato, sorteio com semente fixa (reprodutivel).

Uso:
  python3 make_subset.py --runs runs_download.tsv.gz --studies subset_studies.tsv \
      --attrs stage2/attrs --max-samples 100 --out runs_subset.tsv
"""
import argparse, csv, gzip, os, random, re
from collections import Counter, defaultdict

csv.field_size_limit(10 ** 8)

# chaves de atributo que descrevem o desenho (tratamento, genotipo, tempo)
DESIGN = re.compile(r"treat|perturb|chem.?admin|agrochem|watering|fertiliz|experimental.?factor|"
                    r"condition|management|regm|irrigat|amend|tillage|rotation|genotype|cultivar|"
                    r"variety|collection.?date|time|stage|season|year|crop", re.I)
MAX_RUNS_PER_SAMPLE = 4
SKIP = re.compile(r"ena-|ncbi|submission|insdc|replicate|block|^rep|plot|sample.?name|title|alias", re.I)


def load_attrs(path):
    d = defaultdict(dict)
    if not os.path.exists(path):
        return d
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d[r["sample"]].setdefault(r["key"], r["value"])
    return d


def design_keys(attrs, samples):
    """Chaves de desenho com 2..20 niveis entre as amostras do estudo."""
    keys = Counter()
    for s in samples:
        for k in attrs.get(s, {}):
            if DESIGN.search(k) and not SKIP.search(k):
                keys[k] += 1
    out = []
    for k, n in keys.items():
        if n < 0.8 * len(samples):
            continue
        lv = {attrs[s].get(k, "") for s in samples}
        if 2 <= len(lv) <= 20:
            out.append(k)
    return sorted(out)


def stratified(samples, strata, n, rng):
    groups = defaultdict(list)
    for s in samples:
        groups[strata[s]].append(s)
    for g in groups.values():
        rng.shuffle(g)
    tot = len(samples)
    if len(groups) >= n:                      # mais estratos que vagas: 1 de cada, sorteados
        keys = sorted(groups); rng.shuffle(keys)
        return [groups[k][0] for k in keys[:n]]
    alloc = {k: max(1, round(n * len(v) / tot)) for k, v in groups.items()}
    while sum(alloc.values()) > n:            # ajusta arredondamento
        k = max(alloc, key=lambda x: alloc[x]); alloc[k] -= 1
    while sum(alloc.values()) < n:
        k = max(groups, key=lambda x: len(groups[x]) - alloc[x])
        if len(groups[k]) <= alloc[k]:
            break
        alloc[k] += 1
    return [s for k, v in groups.items() for s in v[: min(alloc[k], len(v))]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True)
    ap.add_argument("--studies", required=True)
    ap.add_argument("--attrs", default="stage2/attrs")
    ap.add_argument("--max-samples", type=int, default=100)
    ap.add_argument("--out", default="runs_subset.tsv")
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()

    with open(a.studies) as fh:
        keep = {r["study_accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
    op = gzip.open if a.runs.endswith(".gz") else open
    with op(a.runs, "rt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        fields = rd.fieldnames
        by = defaultdict(list)
        for r in rd:
            if r["study_accession"] in keep:
                by[r["study_accession"]].append(r)

    rows, log = [], []
    for s in sorted(by):
        rs = by[s]
        # unidade amostral: BioSample, ou a propria corrida quando o submetedor
        # reaproveitou BioSamples (>4 corridas por BioSample)
        per = Counter(r["sample_accession"] for r in rs)
        unit_run = max(per.values()) > MAX_RUNS_PER_SAMPLE
        uid = (lambda r: r["run_accession"]) if unit_run else (lambda r: r["sample_accession"])
        samples = sorted({uid(r) for r in rs})
        if len(samples) <= a.max_samples:
            rows += rs; log.append((s, len(samples), len(samples), "todas" + (" (unidade=corrida)" if unit_run else ""), "")); continue
        attrs = load_attrs(os.path.join(a.attrs, f"{s}.tsv"))
        if unit_run:                       # atributos vem do BioSample de cada corrida
            r2s = {r["run_accession"]: r["sample_accession"] for r in rs}
            attrs = {run: attrs.get(r2s[run], {}) for run in samples}
        keys = design_keys(attrs, samples)
        strata = {x: " | ".join(attrs.get(x, {}).get(k, "") for k in keys) for x in samples}
        rng = random.Random(f"{a.seed}-{s}")
        chosen = set(stratified(samples, strata, a.max_samples, rng))
        rows += [r for r in rs if uid(r) in chosen]
        log.append((s, len(samples), len(chosen), f"{len(set(strata.values()))} estratos" + (" (unidade=corrida)" if unit_run else ""),
                    ";".join(keys)))

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fields, delimiter="\t"); w.writeheader(); w.writerows(rows)
    with open(os.path.splitext(a.out)[0] + "_amostragem.tsv", "w") as fh:
        fh.write("study\tamostras_total\tamostras_usadas\testratificacao\tchaves_desenho\n")
        for x in log:
            fh.write("\t".join(map(str, x)) + "\n")
    sub = [x for x in log if x[3] != "todas"]
    print(f"{len(by)} estudos | {len(rows)} corridas | "
          f"{sum(x[2] for x in log)} amostras -> {a.out}")
    print(f"{len(sub)} estudos subamostrados (> {a.max_samples} amostras); "
          f"{sum(1 for x in sub if x[4])} com estratificacao por desenho")


if __name__ == "__main__":
    main()
