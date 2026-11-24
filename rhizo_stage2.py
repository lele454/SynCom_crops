#!/usr/bin/env python3
# [S1.2] Amplicon region per run and BioSample attributes
# Supplementary 1 of the manuscript. See README.md for the full map.
"""
Etapa 2 - classificacao fina e curadoria dos estudos candidatos.

Entradas (da etapa 1, rhizo_survey.py):  survey/studies.tsv, survey/runs.tsv

  A. Regiao por corrida: amostra K corridas por estudo; detecta primers OU, se
     ja removidos, o motivo conservado logo apos o sitio do primer.
     Estudos mistos (16S + ITS etc.): corridas nao inspecionadas ficam como
     'a_confirmar' e sao checadas no download (ou --mixed all para classificar ja).
  B. Atributos de cada amostra (EBI BioSamples, cobre SAMN/SAMEA/SAMD).
  C. Curadoria: compartimento, campo vs. controlado, fatores experimentais
     (atributos que variam dentro do estudo).

Saidas em stage2/:
  region_cache.tsv     classificacao das corridas inspecionadas (retomavel)
  attrs/<estudo>.tsv   atributos longos: sample, key, value (retomavel)
  runs_classified.tsv  todas as corridas dos estudos 16S V4/V3-V4
  studies_curation.tsv um estudo por linha, para curadoria manual

Uso:
  python3 rhizo_stage2.py --survey survey --out stage2 --threads 4
  python3 rhizo_stage2.py --only region     # so a etapa A
"""
import argparse, csv, gzip, json, os, re, sys, time, urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

csv.field_size_limit(10 ** 8)

# ---------------------------------------------------------------- regioes ---
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
         "S": "[CG]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
         "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
iu = lambda s: "".join(IUPAC[c] for c in s)

# (nome, classe, regex, janela maxima de inicio)
# classe: F_xxx = evidencia do inicio do amplicon; R_806 = fim em 806; ITS; OUTRO
SIGNALS = [
    # primers presentes
    ("338F",  "F_V3",  iu("ACTCCTACGGGAGGCAGCAG"), 60),
    ("341F",  "F_V3",  iu("CCTACGGGNGGCWGCAG"), 60),
    ("341Fb", "F_V3",  iu("CCTACGGGNBGCASCAG"), 60),
    ("515F",  "F_V4",  iu("GTGYCAGCMGCCGCGGTAA"), 60),
    ("799F",  "F_V5",  iu("AACMGGATTAGATACCCKG"), 60),
    ("27F",   "F_V1",  iu("AGAGTTTGATYMTGGCTCAG"), 60),
    ("806R",  "R_806", iu("GGACTACHVGGGTWTCTAAT"), 60),
    ("805R",  "R_806", iu("GACTACHVGGGTATCTAATCC"), 60),
    ("ITS1F", "ITS",   iu("CTTGGTCATTTAGAGGAAGTAA"), 60),
    ("gITS7", "ITS",   iu("GTGARTCATCGARTCTTTG"), 60),
    ("ITS3",  "ITS",   iu("GCATCGATGAAGAACGCAGC"), 60),
    ("ITS4",  "ITS",   iu("TCCTCCGCTTATTGATATGC"), 60),
    ("ITS2",  "ITS",   iu("GCTGCGTTCTTCATCGATGC"), 60),
    # primers ja removidos: motivo conservado logo apos o sitio do primer
    ("m515F", "F_V4",  r"TAC[AG][GT][AG]GG",        6),   # E. coli 534: TACGGAGGG
    ("m341F", "F_V3",  r"T[ACG]GGGAAT[ACT]TT",      6),   # E. coli 358: TGGGGAATATT
    ("m799F", "F_V5",  r"GT[AG]AACGATG",            6),   # E. coli 820: GTAAACGATG
    ("m27F",  "F_V1",  r"A[AGT][CT]GAACGCTGGCGG",   6),   # E. coli 28:  ATTGAACGCTGGCGG
    ("m806R", "R_806", r"CCTGTT[CT]GCT[ACG]CCC[AC][CT]GC", 6),  # revcomp E. coli 763-786
    ("mITS",  "ITS",   r"GTCGTAACAAGG|AACGCA[CT]ATTGCGC|GCATCGATGAAGAACGC", 40),
]
SIG_RE = [(n, c, re.compile(p), w) for n, c, p, w in SIGNALS]
MIN_FRAC = 0.3


def head_seqs(url, n):
    out = []
    with urllib.request.urlopen(url, timeout=45) as r:
        for i, line in enumerate(gzip.GzipFile(fileobj=r)):
            if i % 4 == 1:
                out.append(line.decode(errors="ignore").strip().upper())
                if len(out) >= n:
                    break
    return out


def scan(seqs):
    """Conta, por classe, as leituras com algum sinal (primeiro sinal que bater)."""
    cls, names = Counter(), Counter()
    for s in seqs:
        for n, c, rx, w in SIG_RE:
            m = rx.search(s[: w + 30])
            if m and m.start() <= w:
                cls[c] += 1; names[n] += 1
                break
    return cls, names


CLASS_OF = {n: c for n, c, _, _ in SIGNALS}


def decide(cls, n):
    """Regras (fracoes sobre leituras por arquivo):
    - ITS com >=10% -> ITS (exclusao conservadora)
    - forward: classe com mais leituras, se >=15%; V3 e V4 juntos e proximos -> conflito
    - so reverso 806/805 (>=30%) -> V4|V3-V4? (contem V4; resolvido no processamento)"""
    fr = lambda c: cls[c] / max(n, 1)
    if fr("ITS") >= 0.10:
        return "ITS"
    fw = sorted(((fr(c), c) for c in ("F_V4", "F_V3", "F_V5", "F_V1")), reverse=True)
    (f1, c1), (f2, c2) = fw[0], fw[1]
    if f1 >= 0.15:
        if {c1, c2} == {"F_V4", "F_V3"} and f2 >= 0.15 and f1 < 2 * f2:
            return "V3-V4|V4"
        return {"F_V4": "V4", "F_V3": "V3-V4", "F_V5": "V5-V7", "F_V1": "V1-V3"}[c1]
    if fr("R_806") >= MIN_FRAC:
        return "V4|V3-V4?"
    return "desconhecido"          # nao 16S (nifH, AMF...) ou regiao nao prevista


def redecide(row):
    """Reaplica as regras a partir da coluna 'signals' do cache (sem baixar nada)."""
    if row["region"].startswith("erro") or not row.get("signals"):
        return row
    cls = Counter()
    for kv in row["signals"].split(";"):
        k, _, v = kv.partition(":")
        if k in CLASS_OF and v.isdigit():
            cls[CLASS_OF[k]] += int(v)
    row["region"] = decide(cls, int(row.get("n_reads") or 0))
    return row


FASTQ_DUMP = None   # caminho do fastq-dump (definido em main)


def fastqdump_seqs(run, n, timeout=90):
    """Primeiras n leituras de cada mate via SRA Toolkit (le so o necessario)."""
    import subprocess, tempfile, glob
    with tempfile.TemporaryDirectory() as td:
        subprocess.run([FASTQ_DUMP, "-X", str(n), "--split-files", "--skip-technical",
                        "-O", td, run], check=True, timeout=timeout,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        files = sorted(glob.glob(os.path.join(td, "*.fastq")))
        if not files:
            raise IOError("fastq-dump sem saida")
        out = []
        for f in files[:2]:
            with open(f) as fh:
                out.append([l.strip().upper() for i, l in enumerate(fh) if i % 4 == 1])
        return out


def classify_run(run, ftp, n_reads, retries=2):
    err = ""
    for i in range(retries):
        try:
            if FASTQ_DUMP:
                try:
                    per_file = fastqdump_seqs(run, n_reads)
                except Exception as e:          # cai para o ENA
                    err = f"fqd:{type(e).__name__}"
                    per_file = None
            else:
                per_file = None
            if per_file is None:
                files = [f for f in ftp.split(";") if f]
                pair = [f for f in files if re.search(r"_[12]\.f(ast)?q\.gz$", f)] or files
                per_file = [head_seqs(f if "://" in f else "https://" + f, n_reads)
                            for f in pair[:2]]
            cls, names, n = Counter(), Counter(), 0
            for seqs in per_file:
                c, nm = scan(seqs)
                cls += c; names += nm; n = max(n, len(seqs))
            # leituras podem vir em orientacao mista: n = leituras por arquivo
            return {"run_accession": run, "region": decide(cls, n), "n_reads": n,
                    "signals": ";".join(f"{k}:{v}" for k, v in names.most_common())}
        except Exception as e:
            err = type(e).__name__
            time.sleep(2)
    return {"run_accession": run, "region": f"erro:{err}", "n_reads": 0, "signals": ""}


# ------------------------------------------------------------- atributos ---
BIOSAMPLES = "https://www.ebi.ac.uk/biosamples/samples/"
SKIP_KEYS = {"organism", "scientific_name", "scientific name", "taxid", "ena-checklist",
             "insdc center name", "insdc first public", "insdc last update",
             "insdc status", "submitter id", "sra accession", "ena first public",
             "ena last update", "external id", "insdc center alias", "common name",
             "title", "description", "sample name", "sample_name", "alias"}


def fetch_sample(acc, retries=4):
    for i in range(retries):
        try:
            req = urllib.request.Request(BIOSAMPLES + acc, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                js = json.load(r)
            rows = [(acc, "_name", js.get("name", ""))]
            for k, vs in (js.get("characteristics") or {}).items():
                for v in vs:
                    rows.append((acc, k, str(v.get("text", ""))))
            return rows
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return [(acc, "_erro", "404")]
            time.sleep(2 * 2 ** i)
        except Exception:
            time.sleep(2 * 2 ** i)
    return [(acc, "_erro", "falha")]


def study_attrs(study, samples, outdir, threads):
    f = os.path.join(outdir, f"{study}.tsv")
    if os.path.exists(f):
        with open(f) as fh:
            return [tuple(r) for r in csv.reader(fh, delimiter="\t")][1:]
    rows = []
    with ThreadPoolExecutor(threads) as ex:
        for res in ex.map(fetch_sample, samples):
            rows += res
    with open(f + ".part", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["sample", "key", "value"]); w.writerows(rows)
    os.replace(f + ".part", f)
    return rows


COMP = [("bulk",        r"\bbulk\b|unplanted|non[- ]?rhizosph|unvegetated|fallow soil|control soil"),
        ("endosfera",   r"endosph|endophyt|root (tissue|interior|endo)|inside root|surface[- ]steril"),
        ("rizoplano",   r"rhizoplane"),
        ("rizosfera",   r"rhizosph|rhizo soil|root[- ]adher"),
        ("raiz",        r"\broots?\b"),
        ("solo",        r"\bsoil\b")]
SETTING = [("controlado", r"greenhouse|glasshouse|growth chamber|phytotron|\bpots?\b|potted|mesocosm|microcosm|climate chamber"),
           ("campo",      r"\bfield\b|\bfarm|\bplots?\b|orchard|plantation|paddy")]
INOC = r"inocul|synthetic communit|syncom|biofertili|bioinoc|biocontrol agent"


def first_match(text, rules):
    for lab, p in rules:
        if re.search(p, text, re.I):
            return lab
    return ""


def sample_features(attrs, extra_text=""):
    """attrs: {key: [values]} de uma amostra."""
    vals = " | ".join(v for k, vs in attrs.items() if k.lower() not in SKIP_KEYS
                      and not k.startswith("_") for v in vs)
    txt = f"{vals} | {extra_text}"
    return {"compartimento": first_match(txt, COMP),
            "ambiente": first_match(txt, SETTING),
            "inoculacao": "sim" if re.search(INOC, txt, re.I) else ""}


def study_factors(per_sample):
    """Atributos que variam entre amostras do estudo = fatores candidatos."""
    keys = defaultdict(dict)
    for s, attrs in per_sample.items():
        for k, vs in attrs.items():
            if k.lower() in SKIP_KEYS or k.startswith("_"):
                continue
            keys[k][s] = ";".join(vs)
    n = len(per_sample)
    out = []
    for k, d in keys.items():
        u = Counter(d.values())
        if len(u) < 2:
            continue
        numeric = all(re.fullmatch(r"-?[\d.]+", v or "x") for v in u)
        if len(u) >= 0.9 * n and not numeric:      # parece identificador
            continue
        ex = ",".join(v[:25] for v, _ in u.most_common(4))
        out.append(f"{k}[{len(u)}]: {ex}")
    return " || ".join(sorted(out))


# ----------------------------------------------------------------- main ---
USE = {"V4", "V3-V4", "V3-V4|V4", "V4|V3-V4?"}   # todos contem a regiao V4


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--survey", default="survey")
    ap.add_argument("--out", default="stage2")
    ap.add_argument("--runs-per-study", type=int, default=5)
    ap.add_argument("--reads", type=int, default=200)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--only", choices=["region", "attrs", "final"])
    ap.add_argument("--mixed", choices=["defer", "all"], default="defer",
                    help="estudos mistos: 'defer' marca corridas nao inspecionadas como "
                         "a_confirmar (checadas no download); 'all' classifica todas agora")
    ap.add_argument("--no-fastq-dump", action="store_true", help="ler FASTQ so pelo ENA")
    a = ap.parse_args()
    global FASTQ_DUMP
    import shutil
    FASTQ_DUMP = None if a.no_fastq_dump else shutil.which("fastq-dump")
    os.makedirs(os.path.join(a.out, "attrs"), exist_ok=True)
    log = lambda *x: print(*x, file=sys.stderr, flush=True)

    with open(os.path.join(a.survey, "studies.tsv")) as fh:
        studies = {r["study_accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
    drop = {"ITS", "V5-V7", "V1-V3"}
    cand = [s for s, r in studies.items() if r["region"] not in drop]
    with open(os.path.join(a.survey, "runs.tsv")) as fh:
        runs = [r for r in csv.DictReader(fh, delimiter="\t") if r["study_accession"] in set(cand)]
    by_study = defaultdict(list)
    for r in runs:
        by_study[r["study_accession"]].append(r)
    for s in by_study:
        by_study[s].sort(key=lambda r: r["run_accession"])
    log(f"Estudos candidatos: {len(cand)} | corridas: {len(runs)}")
    log(f"Leitura de FASTQ: {'fastq-dump (' + FASTQ_DUMP + ')' if FASTQ_DUMP else 'ENA (https)'}")

    # ---- A. regiao
    cache_f = os.path.join(a.out, "region_cache.tsv")
    cache = {}
    if os.path.exists(cache_f):
        with open(cache_f) as fh:
            cache = {r["run_accession"]: redecide(r) for r in csv.DictReader(fh, delimiter="\t")
                     if not r["region"].startswith("erro")}
        with open(cache_f, "w", newline="") as fh:      # regrava sem os erros
            w = csv.DictWriter(fh, ["run_accession", "region", "n_reads", "signals"], delimiter="\t")
            w.writeheader(); w.writerows(cache.values())

    def run_batch(todo, label):
        if not todo:
            return
        log(f"{label}: {len(todo)} corridas")
        with open(cache_f, "a", newline="") as fh, ThreadPoolExecutor(a.threads) as ex:
            w = csv.DictWriter(fh, ["run_accession", "region", "n_reads", "signals"], delimiter="\t")
            if fh.tell() == 0:
                w.writeheader()
            futs = [ex.submit(classify_run, r["run_accession"], r.get("fastq_ftp", ""), a.reads)
                    for r in todo]
            for i, f in enumerate(as_completed(futs), 1):
                res = f.result(); cache[res["run_accession"]] = res
                w.writerow(res); fh.flush()
                if i % 20 == 0 or i == len(todo):
                    log(f"  {i}/{len(todo)}")

    def sample_runs(rs, k):
        rs = [r for r in rs if r.get("fastq_ftp")]
        if len(rs) <= k:
            return rs
        step = len(rs) / k                                  # espalhadas pelo estudo
        return [rs[int(i * step)] for i in range(k)]

    if a.only in (None, "region"):
        todo = [r for s in cand for r in sample_runs(by_study[s], a.runs_per_study)
                if r["run_accession"] not in cache]
        run_batch(todo, "Regiao (amostra por estudo)")

    def study_region(s):
        """Regiao do estudo. Variantes que contem V4 (V4, V3-V4...) contam como um grupo;
        'misto' so quando ha corridas com V4 E corridas sem V4 (ITS, desconhecido...)."""
        v = Counter(cache[r["run_accession"]]["region"] for r in by_study[s]
                    if r["run_accession"] in cache)
        v = Counter({k: n for k, n in v.items() if not k.startswith("erro")})
        if not v:
            return "nao_inferido", v
        groups = {("V4grp" if k in USE else k) for k in v}
        if len(groups) == 1:
            if groups == {"V4grp"}:
                firm = Counter({k: n for k, n in v.items() if k in ("V4", "V3-V4")})
                return (firm or v).most_common(1)[0][0], v
            return next(iter(v)), v
        return "misto", v

    if a.only in (None, "region") and a.mixed == "all":
        mixed = [s for s in cand if study_region(s)[0] == "misto"]
        todo = [r for s in mixed for r in by_study[s]
                if r.get("fastq_ftp") and r["run_accession"] not in cache]
        run_batch(todo, f"Estudos mistos ({len(mixed)}): todas as corridas")

    def run_region(r):
        s = r["study_accession"]
        sr, _ = study_region(s)
        if sr == "misto":
            # corrida ja inspecionada -> sua regiao; senao confirmar no download
            return cache.get(r["run_accession"], {}).get("region", "a_confirmar")
        return sr

    keep_runs = [r for r in runs if run_region(r) in USE | {"a_confirmar"}]
    keep_studies = sorted({r["study_accession"] for r in keep_runs})
    n_conf = sum(run_region(r) == "a_confirmar" for r in keep_runs)
    mixed_n = sum(study_region(s)[0] == "misto" for s in cand)
    log(f"16S V4/V3-V4: {len(keep_studies)} estudos, {len(keep_runs)} corridas "
        f"({n_conf} a confirmar no download, de {mixed_n} estudos mistos)")
    st_reg = Counter(study_region(s)[0] for s in cand)
    log("Regiao por estudo: " + ", ".join(f"{k}:{v}" for k, v in st_reg.most_common()))
    if a.only == "region":
        return

    # ---- B. atributos
    samples_by_study = defaultdict(list)
    for r in keep_runs:
        if r.get("sample_accession"):
            samples_by_study[r["study_accession"]].append(r["sample_accession"])
    attrs = {}
    if a.only in (None, "attrs", "final"):
        pending = [s for s in keep_studies
                   if not os.path.exists(os.path.join(a.out, "attrs", f"{s}.tsv"))]
        if a.only != "final":
            log(f"Atributos de amostras: {len(pending)} estudos a baixar")
        for i, s in enumerate(keep_studies, 1):
            if a.only == "final" and s in pending:
                continue
            rows = study_attrs(s, sorted(set(samples_by_study[s])),
                               os.path.join(a.out, "attrs"), a.threads)
            d = defaultdict(lambda: defaultdict(list))
            for smp, k, v in rows:
                d[smp][k].append(v)
            attrs[s] = d
            if s in pending and (i % 10 == 0 or i == len(keep_studies)):
                log(f"  {i}/{len(keep_studies)}")
    if a.only == "attrs":
        return

    # ---- C. tabelas finais
    out_runs = os.path.join(a.out, "runs_classified.tsv")
    cols = ["run_accession", "study_accession", "sample_accession", "region",
            "compartimento", "ambiente", "inoculacao", "read_count", "host", "country",
            "fastq_ftp"]
    per_study = defaultdict(lambda: defaultdict(Counter))
    with open(out_runs, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore"); w.writeheader()
        for r in keep_runs:
            s = r["study_accession"]
            sa = attrs.get(s, {}).get(r.get("sample_accession", ""), {})
            extra = f'{r.get("sample_title", "")} | {r.get("experiment_title", "")}'
            feat = sample_features(sa, extra)
            row = dict(r, region=run_region(r), **feat)
            w.writerow(row)
            for k in ("region", "compartimento", "ambiente", "inoculacao"):
                per_study[s][k][row[k] or "NA"] += 1

    out_st = os.path.join(a.out, "studies_curation.tsv")
    with open(out_st, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["study_accession", "n_runs_16S", "region", "compartimento", "ambiente",
                    "inoculacao", "fatores_candidatos", "hosts", "countries", "study_title"])
        for s in keep_studies:
            ps = per_study[s]
            fmt = lambda c: ";".join(f"{k}:{v}" for k, v in c.most_common())
            st = studies.get(s, {})
            w.writerow([s, sum(ps["region"].values()), fmt(ps["region"]),
                        fmt(ps["compartimento"]), fmt(ps["ambiente"]), fmt(ps["inoculacao"]),
                        study_factors(attrs.get(s, {})), st.get("hosts", ""),
                        st.get("countries", ""), st.get("study_title", "")])

    # resumo
    comp = Counter(); amb = Counter(); inoc = 0
    for s in keep_studies:
        ps = per_study[s]
        comp[ps["compartimento"].most_common(1)[0][0]] += 1
        amb[ps["ambiente"].most_common(1)[0][0]] += 1
        inoc += ps["inoculacao"].get("sim", 0) > 0
    print(f"\n{len(keep_studies)} estudos 16S V4/V3-V4, {len(keep_runs)} corridas -> {out_st}")
    print("compartimento predominante:", dict(comp.most_common()))
    print("ambiente predominante:", dict(amb.most_common()))
    print("estudos com indicio de inoculacao:", inoc)


if __name__ == "__main__":
    main()
