#!/usr/bin/env python3
# [S1.1] ENA survey of candidate studies
# Supplementary 1 of the manuscript. See README.md for the full map.
"""
Levantamento de corridas 16S de rizosfera no ENA (espelha SRA/DDBJ - INSDC).

Etapas:
  1. Consulta o ENA Portal API: tax_tree(939928) [rhizosphere metagenome],
     AMPLICON, ILLUMINA, PAIRED  ->  runs.tsv
  2. Agrupa por estudo e filtra por n minimo de corridas
  3. Infere os primers baixando so as primeiras N leituras (streaming) de
     algumas corridas por estudo  ->  sniff_cache.tsv (retomavel)
  4. Resume por estudo e por regiao  ->  studies.tsv

Sem dependencias externas (Python >= 3.8).

Uso:
  python3 rhizo_survey.py --outdir survey --min-runs 20 --runs-per-study 3 --threads 8
  python3 rhizo_survey.py --count-only      # so conta corridas
"""
import argparse, csv, gzip, os, re, sys, time, urllib.parse, urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

ENA = "https://www.ebi.ac.uk/ena/portal/api"

QUERY = None
DEFAULT_QUERY = ('tax_tree(939928) AND library_strategy="AMPLICON" '
                 'AND instrument_platform="ILLUMINA" AND library_layout="PAIRED"')

WANTED_FIELDS = [
    "run_accession", "study_accession", "secondary_study_accession",
    "sample_accession", "study_title", "sample_title", "experiment_title",
    "library_layout", "instrument_model", "read_count", "base_count",
    "nominal_length", "host", "country", "collection_date", "lat", "lon",
    "environment_biome", "environment_material", "first_public",
    "tax_id", "scientific_name", "fastq_ftp", "fastq_bytes",
]

# --- primers (5'->3', IUPAC). Ordem importa: o 338F contem o 341F. ---------
PRIMERS = [
    ("338F",   "F", "ACTCCTACGGGAGGCAGCAG"),
    ("341F",   "F", "CCTACGGGNGGCWGCAG"),
    ("341F_b", "F", "CCTACGGGNBGCASCAG"),
    ("515F",   "F", "GTGYCAGCMGCCGCGGTAA"),
    ("799F",   "F", "AACMGGATTAGATACCCKG"),
    ("27F",    "F", "AGAGTTTGATYMTGGCTCAG"),
    ("ITS1F",  "F", "CTTGGTCATTTAGAGGAAGTAA"),
    ("gITS7",  "F", "GTGARTCATCGARTCTTTG"),
    ("ITS3",   "F", "GCATCGATGAAGAACGCAGC"),
    ("806R",   "R", "GGACTACHVGGGTWTCTAAT"),
    ("805R",   "R", "GACTACHVGGGTATCTAATCC"),
    ("1193R",  "R", "ACGTCATCCCCACCTTCC"),
    ("ITS4",   "R", "TCCTCCGCTTATTGATATGC"),
    ("ITS2",   "R", "GCTGCGTTCTTCATCGATGC"),
]
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
         "S": "[CG]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
         "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
PRIMER_RE = [(n, d, re.compile("".join(IUPAC[c] for c in s))) for n, d, s in PRIMERS]
WINDOW = 30          # primer pode comecar ate a base 30 (espacadores/phasing)
MIN_FRAC = 0.3       # fracao minima de leituras com o primer


def http_get(url, params=None, retries=4, timeout=120):
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read().decode()
        except Exception as e:
            if i == retries - 1:
                raise
            time.sleep(2 ** i)


def ena_fields():
    txt = http_get(f"{ENA}/returnFields", {"result": "read_run", "format": "tsv"})
    return {l.split("\t")[0] for l in txt.splitlines()[1:] if l}


def ena_count(query):
    return int(http_get(f"{ENA}/count", {"result": "read_run", "query": query}).strip().splitlines()[-1])


LITE_FIELDS = ["run_accession", "study_accession", "secondary_study_accession",
               "sample_accession", "instrument_model", "read_count",
               "base_count", "first_public"]


def study_details(study, fields, outdir):
    """Metadados completos de um estudo (so corridas de rizosfera da query)."""
    f = os.path.join(outdir, f"{study}.tsv")
    if not os.path.exists(f):
        q = f'({QUERY}) AND study_accession="{study}"'
        txt = http_get(f"{ENA}/search", {"result": "read_run", "query": q,
                                         "fields": ",".join(fields), "format": "tsv",
                                         "limit": 0}, timeout=300)
        with open(f + ".part", "w") as fh:
            fh.write(txt)
        os.replace(f + ".part", f)
    with open(f) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ena_search(query, fields, out, meta_name="meta_lite"):
    """Baixa por ano de publicacao (first_public), gravando em disco por streaming.
    Cada ano fica em <outdir>/meta/runs_<ano>.tsv; anos ja baixados sao pulados."""
    import shutil
    meta = os.path.join(os.path.dirname(out) or ".", meta_name)
    os.makedirs(meta, exist_ok=True)

    def fetch(q, part, label, tries):
        """Baixa um pedaco; confere n de linhas contra /count. True se ok."""
        url = f"{ENA}/search?" + urllib.parse.urlencode(
            {"result": "read_run", "query": q, "fields": ",".join(fields),
             "format": "tsv", "limit": 0})
        for i in range(tries):
            tmp = part + ".part"
            try:
                expected = ena_count(q)
                with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as fh:
                    shutil.copyfileobj(r, fh, 1 << 16)
                with open(tmp) as fh:
                    n = max(sum(1 for _ in fh) - 1, 0)
                if n != expected:
                    raise IOError(f"{n} linhas, esperado {expected}")
                os.replace(tmp, part)
                print(f"  {label}: {n} corridas", file=sys.stderr, flush=True)
                return True
            except Exception as e:
                print(f"  {label}: tentativa {i + 1} falhou ({e})", file=sys.stderr, flush=True)
                time.sleep(5 * 2 ** i)
        return False

    import datetime as dt
    D = dt.date
    FLOOR = D(1900, 1, 1)

    def name(a, b):
        # nomes compativeis com as versoes anteriores (ano / mes)
        lo = D(2008, 1, 1) if a == FLOOR else a
        if lo.month == 1 and lo.day == 1 and b == D(lo.year + 1, 1, 1):
            return f"runs_{lo.year}.tsv"
        nxt = D(lo.year + (lo.month == 12), lo.month % 12 + 1, 1)
        if lo.day == 1 and b == nxt:
            return f"runs_{lo.year}_{lo.month:02d}.tsv"
        return f"runs_{a.isoformat()}_{b.isoformat()}.tsv"

    def children(a, b):
        lo = D(2008, 1, 1) if a == FLOOR else a
        if lo.month == 1 and lo.day == 1 and b == D(lo.year + 1, 1, 1):   # ano -> meses
            ms = [D(lo.year, m, 1) for m in range(1, 13)] + [b]
            ms[0] = a
            return list(zip(ms[:-1], ms[1:]))
        mid = a + (b - a) // 2                                           # bissecao
        return [(a, mid), (mid, b)]

    def get(a, b):
        """Baixa [a, b); se falhar, divide o intervalo recursivamente."""
        part = os.path.join(meta, name(a, b))
        if os.path.exists(part):
            return [part]
        subs = children(a, b)
        if any(os.path.exists(os.path.join(meta, name(x, y))) for x, y in subs):
            return [p for x, y in subs for p in get(x, y)]               # retomada
        q = f"({query}) AND first_public>={a.isoformat()} AND first_public<{b.isoformat()}"
        days = (b - a).days
        tries = 2 if days > 1 else 8
        print(f"  {a}..{b}: baixando", file=sys.stderr, flush=True)
        if fetch(q, part, f"{a}..{b}", tries):
            return [part]
        if days <= 1:
            raise RuntimeError(f"Falha persistente em {a}; rode de novo para retomar")
        print(f"  {a}..{b}: dividindo intervalo", file=sys.stderr, flush=True)
        return [p for x, y in subs for p in get(x, y)]

    y_end = time.gmtime().tm_year
    parts = []
    for y in range(2008, y_end + 1):
        a = FLOOR if y == 2008 else D(y, 1, 1)
        parts += get(a, D(y + 1, 1, 1))
    header_done = False
    with open(out, "w") as fo:
        for part in parts:
            with open(part) as fh:
                hdr = fh.readline()
                if not hdr:
                    continue
                if not header_done:
                    fo.write(hdr); header_done = True
                shutil.copyfileobj(fh, fo)


# --- inferencia de primers ---------------------------------------------------
def head_seqs(url, n_reads):
    seqs = []
    with urllib.request.urlopen(url, timeout=90) as r:
        gz = gzip.GzipFile(fileobj=r)
        for i, line in enumerate(gz):
            if i % 4 == 1:
                seqs.append(line.decode(errors="ignore").strip().upper())
                if len(seqs) >= n_reads:
                    break
    return seqs


def primer_hits(seqs):
    hits = Counter()
    for s in seqs:
        w = s[:WINDOW + 25]
        for name, _, rx in PRIMER_RE:          # primeiro match (ordem de prioridade)
            m = rx.search(w)
            if m and m.start() <= WINDOW:
                hits[name] += 1
                break
    return hits


def classify(fwd, rev):
    its = {"ITS1F", "gITS7", "ITS3", "ITS4", "ITS2"}
    if fwd in its or rev in its:
        return "ITS"
    r806 = rev in {"806R", "805R"}
    if fwd == "515F" and (r806 or rev is None):
        return "V4" if r806 else "V4?"          # reverso nao detectado (aparado/R2 ausente)
    if fwd in {"338F", "341F", "341F_b"} and (r806 or rev is None):
        return "V3-V4" if r806 else "V3-V4?"    # recortavel in silico para 515F-806R
    if fwd is None and r806:
        return "V4|V3-V4?"                      # forward aparado; resolver por alinhamento
    if fwd == "799F":
        return "V5-V7"
    if fwd == "27F":
        return "V1-V3"
    if fwd is None and rev is None:
        return "sem_primer"     # primers ja removidos? checar por alinhamento
    return f"{fwd or 'NA'}/{rev or 'NA'}"


def sniff_run(run, ftp_field, n_reads):
    files = [f for f in ftp_field.split(";") if f]
    pair = [f for f in files if re.search(r"_[12]\.f(ast)?q\.gz$", f)] or files
    total, hits = 0, Counter()
    for f in pair[:2]:
        url = f if "://" in f else "https://" + f
        seqs = head_seqs(url, n_reads)
        total += len(seqs)
        hits += primer_hits(seqs)
    # leituras podem ter orientacao mista: procura F e R nos dois arquivos
    nfile = max(1, len(pair[:2]))
    per_file = max(1, total // nfile)
    best = {}
    for d in ("F", "R"):
        cands = [(hits[n], n) for n, dd, _ in PRIMER_RE if dd == d and hits[n]]
        c, n = max(cands) if cands else (0, None)
        best[d] = n if c / per_file >= MIN_FRAC else None
    return {"run_accession": run, "fwd": best["F"] or "", "rev": best["R"] or "",
            "region": classify(best["F"], best["R"]), "n_reads": total,
            "hits": ";".join(f"{k}:{v}" for k, v in hits.most_common())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="survey")
    ap.add_argument("--query", default=DEFAULT_QUERY)
    ap.add_argument("--min-runs", type=int, default=20)
    ap.add_argument("--runs-per-study", type=int, default=3)
    ap.add_argument("--reads", type=int, default=500)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--count-only", action="store_true")
    ap.add_argument("--skip-sniff", action="store_true")
    a = ap.parse_args()

    print(f"Query: {a.query}", file=sys.stderr)
    print(f"Corridas no ENA: {ena_count(a.query)}", file=sys.stderr, flush=True)
    if a.count_only:
        return
    global QUERY
    QUERY = a.query
    os.makedirs(a.outdir, exist_ok=True)
    runs_tsv = os.path.join(a.outdir, "runs_lite.tsv")
    full_tsv = os.path.join(a.outdir, "runs.tsv")
    have_full = os.path.exists(full_tsv) and not os.path.exists(runs_tsv)
    if have_full:                       # metadados completos de versao anterior
        runs_tsv = full_tsv

    # 1. metadados leves de todas as corridas
    avail = ena_fields()
    if not os.path.exists(runs_tsv):
        ena_search(a.query, [f for f in LITE_FIELDS if f in avail], runs_tsv)
    with open(runs_tsv) as fh:
        runs = list(csv.DictReader(fh, delimiter="\t"))
    print(f"Corridas baixadas: {len(runs)}", file=sys.stderr)

    # 2. por estudo
    by_study = defaultdict(list)
    for r in runs:
        by_study[r["study_accession"]].append(r)
    keep_ids = [s for s, rs in by_study.items() if len(rs) >= a.min_runs]
    print(f"Estudos: {len(by_study)} | com >= {a.min_runs} corridas: {len(keep_ids)}",
          file=sys.stderr, flush=True)

    # 2b. metadados completos so dos estudos mantidos
    fields = [f for f in WANTED_FIELDS if f in avail]
    miss = set(WANTED_FIELDS) - set(fields)
    if miss:
        print(f"Campos indisponiveis (ignorados): {sorted(miss)}", file=sys.stderr)
    det_dir = os.path.join(a.outdir, "study_runs")
    os.makedirs(det_dir, exist_ok=True)
    keep = {}
    if have_full:
        keep = {s: by_study[s] for s in keep_ids}
        keep_ids = []
    print(f"Baixando metadados completos de {len(keep_ids)} estudos...", file=sys.stderr, flush=True)
    with ThreadPoolExecutor(min(a.threads, 8)) as ex:
        futs = {ex.submit(study_details, s, fields, det_dir): s for s in keep_ids}
        for i, f in enumerate(as_completed(futs), 1):
            s = futs[f]
            try:
                keep[s] = f.result() or by_study[s]
            except Exception as e:
                print(f"  {s}: falhou ({e}); usando metadados leves", file=sys.stderr)
                keep[s] = by_study[s]
            if i % 50 == 0:
                print(f"  {i}/{len(keep_ids)}", file=sys.stderr, flush=True)
    if not have_full:
      with open(full_tsv, "w", newline="") as fh:
          w = csv.DictWriter(fh, fields, delimiter="\t", extrasaction="ignore")
          w.writeheader()
          for rs in keep.values():
              w.writerows(rs)

    # 3. primers (cache retomavel)
    cache_f = os.path.join(a.outdir, "sniff_cache.tsv")
    cache = {}
    if os.path.exists(cache_f):
        with open(cache_f) as fh:
            cache = {r["run_accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
    todo = []
    for s, rs in keep.items():
        rs = sorted(rs, key=lambda r: -int(r.get("read_count") or 0))
        rs = [r for r in rs if r.get("fastq_ftp")][: a.runs_per_study]
        todo += [r for r in rs if r["run_accession"] not in cache]
    cols = ["run_accession", "fwd", "rev", "region", "n_reads", "hits"]
    if todo and not a.skip_sniff:
        print(f"Inferindo primers em {len(todo)} corridas...", file=sys.stderr)
        new = not os.path.exists(cache_f)
        with open(cache_f, "a", newline="") as fh, ThreadPoolExecutor(a.threads) as ex:
            w = csv.DictWriter(fh, cols, delimiter="\t")
            if new:
                w.writeheader()
            futs = {ex.submit(sniff_run, r["run_accession"], r["fastq_ftp"], a.reads): r
                    for r in todo}
            for i, f in enumerate(as_completed(futs), 1):
                r = futs[f]
                try:
                    res = f.result()
                except Exception as e:
                    res = {"run_accession": r["run_accession"], "fwd": "", "rev": "",
                           "region": f"erro:{type(e).__name__}", "n_reads": 0, "hits": ""}
                cache[res["run_accession"]] = res
                w.writerow(res); fh.flush()
                if i % 10 == 0 or i == len(todo):
                    print(f"  {i}/{len(todo)}", file=sys.stderr, flush=True)

    # 4. resumo por estudo
    out = os.path.join(a.outdir, "studies.tsv")
    reg_count, reg_runs = Counter(), Counter()
    with open(out, "w", newline="") as fh:
        hdr = ["study_accession", "secondary_study_accession", "n_runs", "region",
               "region_votes", "fwd", "rev", "median_reads", "instruments",
               "hosts", "countries", "first_public", "study_title"]
        w = csv.writer(fh, delimiter="\t"); w.writerow(hdr)
        for s, rs in sorted(keep.items(), key=lambda kv: -len(kv[1])):
            sn = [cache[r["run_accession"]] for r in rs if r["run_accession"] in cache]
            votes = Counter(x["region"] for x in sn)
            region = votes.most_common(1)[0][0] if votes else "nao_inferido"
            fwd = Counter(x["fwd"] for x in sn).most_common(1)[0][0] if sn else ""
            rev = Counter(x["rev"] for x in sn).most_common(1)[0][0] if sn else ""
            rc = sorted(int(r.get("read_count") or 0) for r in rs)
            uniq = lambda k, n=5: ";".join(x for x, _ in Counter(r.get(k, "") for r in rs).most_common(n) if x)
            w.writerow([s, rs[0].get("secondary_study_accession", ""), len(rs), region,
                        ";".join(f"{k}:{v}" for k, v in votes.items()), fwd, rev,
                        rc[len(rc) // 2], uniq("instrument_model", 3), uniq("host"),
                        uniq("country"), min(r.get("first_public", "") for r in rs),
                        rs[0].get("study_title", "")])
            reg_count[region] += 1; reg_runs[region] += len(rs)

    print(f"\nResumo (estudos com >= {a.min_runs} corridas) -> {out}")
    print(f"{'regiao':<22}{'estudos':>9}{'corridas':>10}")
    for reg, n in reg_count.most_common():
        print(f"{reg:<22}{n:>9}{reg_runs[reg]:>10}")


if __name__ == "__main__":
    main()
