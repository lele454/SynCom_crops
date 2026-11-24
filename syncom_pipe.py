#!/usr/bin/env python3
# [S1.7] Per-study processing driver
# Supplementary 1 of the manuscript. See README.md for the full map.
"""
Pipeline estudo a estudo: download subamostrado -> regiao -> cutadapt -> DADA2 -> V4.

Dois modos, para rodar em paralelo (cada um retomavel):
  download : baixa os primeiros N pares de cada corrida (fastq-dump -X) para o HD
  process  : processa os estudos ja baixados (um por vez, DADA2 multithread)

Uso:
  python3 syncom_pipe.py download --root /media/HD/syncom --runs runs_download.tsv --threads 8
  python3 syncom_pipe.py process  --root /media/HD/syncom --runs runs_download.tsv --threads 6
  python3 syncom_pipe.py status   --root /media/HD/syncom --runs runs_download.tsv

Requer no PATH: fastq-dump (SRA Toolkit), cutadapt, Rscript (+ pacote dada2).
rhizo_stage2.py e dada2_study.R devem estar na mesma pasta deste script.

Estrutura no HD:
  raw/<estudo>/<run>_{1,2}.fastq.gz   leituras subamostradas (mantidas)
  work/<estudo>/                      temporario (apagado ao fim)
  out/<estudo>/                       asv_long.tsv, track.tsv, params.tsv, v4_counts.tsv.gz,
                                      runs_region.tsv, stats.json
  state/<estudo>.{downloaded,done,failed}
"""
import argparse, csv, gzip, json, os, re, shutil, subprocess, sys, time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rhizo_stage2 as R2          # scan(), decide(), IUPAC

csv.field_size_limit(10 ** 8)
log = lambda *x: print(time.strftime("%H:%M:%S"), *x, file=sys.stderr, flush=True)

PRIMERS = {"515F": "GTGYCAGCMGCCGCGGTAA", "338F": "ACTCCTACGGGAGGCAGCAG",
           "341F": "CCTACGGGNGGCWGCAG", "341Fb": "CCTACGGGNBGCASCAG",
           "806R": "GGACTACHVGGGTWTCTAAT", "805R": "GACTACHVGGGTATCTAATCC"}
V4_GROUP = {"V4", "V3-V4", "V3-V4|V4", "V4|V3-V4?"}
_COMP = str.maketrans("ACGTYRMKSWBDHVN", "TGCARYKMSWVHDBN")
revcomp = lambda x: x[::-1].translate(_COMP)


def load_runs(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    by = defaultdict(list)
    for r in rows:
        by[r["study_accession"]].append(r)
    return by


def study_order(by, studies_file=None):
    """'manter' primeiro (se studies_final.tsv for dado), depois por n de corridas."""
    status = {}
    if studies_file and os.path.exists(studies_file):
        with open(studies_file) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                status[r["study_accession"]] = r.get("status", "")
    rank = lambda s: (0 if status.get(s, "manter") == "manter" else 1, -len(by[s]), s)
    order = sorted(by, key=rank)
    return [s for s in order if not status.get(s, "").startswith("excluir")]


def free_gb(path):
    return shutil.disk_usage(path).free / 1e9


# ------------------------------------------------------------------ download ---
def dl_run(run, dest, n, timeout):
    f1, f2 = (os.path.join(dest, f"{run}_{i}.fastq.gz") for i in (1, 2))
    if os.path.exists(f1):
        return run, "ok"
    tmp = os.path.join(dest, f".tmp_{run}")
    shutil.rmtree(tmp, ignore_errors=True); os.makedirs(tmp)
    for i in range(2):
        try:
            subprocess.run(["fastq-dump", "-X", str(n), "--split-3", "--skip-technical",
                            "--gzip", "-O", tmp, run], check=True, timeout=timeout,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            got = sorted(os.listdir(tmp))
            p1 = [g for g in got if g.endswith("_1.fastq.gz")]
            p2 = [g for g in got if g.endswith("_2.fastq.gz")]
            if not p1 or not p2:
                shutil.rmtree(tmp, ignore_errors=True)
                return run, "sem_par"
            os.replace(os.path.join(tmp, p2[0]), f2)
            os.replace(os.path.join(tmp, p1[0]), f1)       # _1 por ultimo = marcador
            shutil.rmtree(tmp, ignore_errors=True)
            return run, "ok"
        except Exception as e:
            err = type(e).__name__
            time.sleep(5)
    shutil.rmtree(tmp, ignore_errors=True)
    return run, f"erro:{err}"


def pending_count(root):
    sd = os.path.join(root, "state")
    if not os.path.isdir(sd):
        return 0
    names = set(os.listdir(sd))
    return sum(1 for f in names if f.endswith(".downloaded")
               and f[:-11] + ".done" not in names and f[:-11] + ".failed" not in names)


def cmd_download(a):
    by = load_runs(a.runs)
    order = study_order(by, a.studies)
    if a.limit:
        order = order[: a.limit]
    os.makedirs(os.path.join(a.root, "state"), exist_ok=True)
    for i, s in enumerate(order, 1):
        mark = os.path.join(a.root, "state", f"{s}.downloaded")
        if os.path.exists(mark):
            continue
        if free_gb(a.root) < a.min_free_gb:
            log(f"Espaco livre < {a.min_free_gb} GB; parando o download."); return
        # controle de fila: nao deixa acumular estudos baixados e nao processados
        warned = False
        while a.max_pending and pending_count(a.root) >= a.max_pending:
            if not warned:
                log(f"{a.max_pending} estudos aguardando processamento; download em pausa...")
                warned = True
            time.sleep(60)
        dest = os.path.join(a.root, "raw", s); os.makedirs(dest, exist_ok=True)
        runs = [r["run_accession"] for r in by[s]]
        res = Counter()
        with ThreadPoolExecutor(a.threads) as ex:
            for run, st in ex.map(lambda r: dl_run(r, dest, a.reads, a.timeout), runs):
                res[st.split(":")[0]] += 1
        with open(mark, "w") as fh:
            json.dump(dict(res), fh)
        log(f"[{i}/{len(order)}] {s}: {dict(res)} | livre {free_gb(a.root):.0f} GB")


# ------------------------------------------------------------------- process ---
def head_fastq(path, n=200):
    out = []
    with gzip.open(path, "rt") as fh:
        for i, l in enumerate(fh):
            if i % 4 == 1:
                out.append(l.strip().upper())
                if len(out) >= n:
                    break
    return out


def run_signals(raw, run):
    """Regiao e primers da corrida a partir das leituras ja baixadas."""
    cls, names, n = Counter(), Counter(), 0
    per = {}
    for i in (1, 2):
        seqs = head_fastq(os.path.join(raw, f"{run}_{i}.fastq.gz"))
        c, nm = R2.scan(seqs)
        cls += c; names += nm; n = max(n, len(seqs)); per[i] = nm
    return R2.decide(cls, n), names, n, per


def primer_suffix_hits(raw, run, n=200):
    """Fracao de leituras (R1 e R2) que trazem o primer inteiro OU so a sua parte 3'
    (primers depositados parcialmente cortados). Usa os ultimos 10 nt do primer."""
    hits = Counter(); tot = 0
    for k in (1, 2):
        seqs = head_fastq(os.path.join(raw, f"{run}_{k}.fastq.gz"), n)
        tot = max(tot, len(seqs))
        for name, pr in PRIMERS.items():
            rx = re.compile(R2.iu(pr[-10:]))
            lim = len(pr) + 8
            hits[name] += sum(1 for x in seqs if rx.search(x[:lim + 10]) and
                              rx.search(x[:lim + 10]).start() <= lim)
    return hits, tot


def pick_primers_suffix(hits_list):
    tot = Counter(); n = 0
    for h, t in hits_list:
        tot += h; n += t
    frac = lambda p: tot[p] / max(n, 1)
    fw = max(("515F", "338F", "341F", "341Fb"), key=lambda p: (round(frac(p), 1), len(PRIMERS[p])))
    rv = max(("806R", "805R"), key=lambda p: (round(frac(p), 1), len(PRIMERS[p])))
    return (fw if frac(fw) >= 0.3 else None), (rv if frac(rv) >= 0.3 else None)


def pick_primers(names_list, n_total):
    """Primers presentes no estudo (nomes reais, nao motivos)."""
    tot = Counter()
    for nm in names_list:
        tot += nm
    fw = [(tot[p], p) for p in ("515F", "338F", "341F", "341Fb")]
    rv = [(tot[p], p) for p in ("806R", "805R")]
    f = max(fw); r = max(rv)
    frac = lambda c: c / max(n_total, 1)
    return (f[1] if frac(f[0]) >= 0.3 else None), (r[1] if frac(r[0]) >= 0.3 else None)


def fq_records(path):
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            yield h, fh.readline(), fh.readline(), fh.readline()


def spot_id(h):
    # "@SRR123.45 45 length=250" ou "@SRR123.45.1 ..." -> "SRR123.45"
    t = h[1:].split()[0]
    p = t.split(".")
    return ".".join(p[:2]) if len(p) >= 2 else t


def sync_pairs(i1, i2, work, run):
    """Se R1 e R2 tem numeros diferentes de leituras, mantem so os pares completos."""
    id1 = [spot_id(r[0]) for r in fq_records(i1)]
    id2 = [spot_id(r[0]) for r in fq_records(i2)]
    if id1 == id2:
        return i1, i2, 0
    n1, n2 = len(id1), len(id2)
    ks = set(id1) & set(id2)
    o1, o2 = (os.path.join(work, f".sync_{run}.r{k}.fq.gz") for k in (1, 2))  # nao casa *_1.fastq.gz
    for src, dst in ((i1, o1), (i2, o2)):
        seen = set()
        with gzip.open(dst, "wt", compresslevel=1) as fo:
            for r in fq_records(src):
                sid = spot_id(r[0])
                if sid in ks and sid not in seen:
                    seen.add(sid); fo.write("".join(r))
    return o1, o2, max(n1, n2) - len(ks)


def cutadapt_run(raw, work, run, fwd, rev, threads, mixed=False):
    try:
        _cutadapt_run(raw, work, run, fwd, rev, threads, mixed)
    finally:
        for k in (1, 2):
            t = os.path.join(work, f".sync_{run}.r{k}.fq.gz")
            if os.path.exists(t):
                os.remove(t)


def _cutadapt_run(raw, work, run, fwd, rev, threads, mixed=False):
    i1, i2 = (os.path.join(raw, f"{run}_{k}.fastq.gz") for k in (1, 2))
    i1, i2, _ = sync_pairs(i1, i2, work, run)
    o1, o2 = (os.path.join(work, f"{run}_{k}.fastq.gz") for k in (1, 2))
    if not fwd and not rev:                       # primers ja removidos
        subprocess.run(["cutadapt", "-j", str(threads), "-m", "100", "--quiet",
                        "-a", "AGATCGGAAGAGC", "-A", "AGATCGGAAGAGC",
                        "-a", "CTGTCTCTTATACACATCT", "-A", "CTGTCTCTTATACACATCT",
                        "-o", o1, "-p", o2, i1, i2], check=True, stdout=subprocess.DEVNULL)
        return

    def cut(g, G, x1, x2):
        # adaptadores ligados: primer no 5' (obrigatorio) ... rc do outro primer no 3'
        # (opcional, remove read-through em leituras mais longas que o amplicon),
        # mais adaptadores Illumina (TruSeq/Nextera) no 3'
        cmd = ["cutadapt", "-j", str(threads), "-e", "0.15", "--discard-untrimmed",
               "-m", "100", "-o", x1, "-p", x2, "--quiet",
               "-a", "AGATCGGAAGAGC", "-A", "AGATCGGAAGAGC",
               "-a", "CTGTCTCTTATACACATCT", "-A", "CTGTCTCTTATACACATCT"]
        if g and not G:
            cmd += ["--pair-filter", "first"]     # so o R1 precisa ter o primer
        elif G and not g:
            cmd += ["--pair-filter", "both"]      # so o R2 tem primer: nao exige do R1
        if g:
            cmd += ["-g", f"{PRIMERS[g]};required...{revcomp(PRIMERS[G])};optional" if G
                    else PRIMERS[g]]
        if G:
            cmd += ["-G", f"{PRIMERS[G]};required...{revcomp(PRIMERS[g])};optional" if g
                    else PRIMERS[G]]
        subprocess.run(cmd + [i1, i2], check=True, stdout=subprocess.DEVNULL)

    cut(fwd, rev, o1, o2)
    if mixed and fwd and rev:
        # pares invertidos (R1 comeca no primer reverso): apara e troca R1<->R2
        b1, b2 = o1 + ".inv1.gz", o2 + ".inv2.gz"
        cut(rev, fwd, b1, b2)
        for o, b in ((o1, b2), (o2, b1)):
            with open(o, "ab") as fo, open(b, "rb") as fb:   # gzip concatenado e valido
                shutil.copyfileobj(fb, fo)
            os.remove(b)


def orientation_mixed(raw, run, fwd, rev):
    """>=20% das leituras do R1 comecam pelo primer reverso?"""
    if not (fwd and rev):
        return False
    rx = re.compile(R2.iu(PRIMERS[rev]))
    seqs = head_fastq(os.path.join(raw, f"{run}_1.fastq.gz"))
    hit = sum(1 for x in seqs if (m := rx.search(x[:80])) and m.start() <= 60)
    return hit / max(len(seqs), 1) >= 0.2


V4_515 = re.compile(R2.iu("GTGYCAGCMGCCGCGGTAA"))


V4_515_END = re.compile(R2.iu("CAGCMGCCGCGGTAA"))       # parte 3' do 515F
V4_806_RC = re.compile(R2.iu(revcomp("GGACTACHVGGGTWTCTAAT")[:12]))   # inicio do rc do 806R


def trim_v4_ends(v):
    """Remove restos de primer quando a ASV e mais longa que a V4 (variantes de primer
    nao reconhecidas): 515F no inicio e rc do 806R no fim."""
    if len(v) > 265:
        m = V4_515_END.search(v[:45])
        if m:
            v = v[m.end():]
    if len(v) > 265:
        m = None
        for m in V4_806_RC.finditer(v):
            pass
        if m and m.start() >= len(v) - 45:
            v = v[:m.start()]
    return v


def to_v4(seq, region):
    """Extrai o fragmento 515F-806R (sem primers). V4: a ASV ja e esse fragmento."""
    if region == "V4":
        v = trim_v4_ends(seq)
    else:
        m = V4_515.search(seq)
        if not m:                                 # tolera 1 erro no sitio do 515F
            m = None
            for i in range(120, min(len(seq) - 19, 240)):
                w = seq[i:i + 19]
                if sum(not re.fullmatch(R2.IUPAC[p], b) for p, b in zip("GTGYCAGCMGCCGCGGTAA", w)) <= 1:
                    m = i + 19; break
            if m is None:
                return None
            v = seq[m:]
        else:
            v = seq[m.end():]
        v = trim_v4_ends(v)
    return v if 240 <= len(v) <= 265 else None


def process_study(a, s, runs_meta):
    raw = os.path.join(a.root, "raw", s)
    work = os.path.join(a.work or a.root, "work", s)
    out = os.path.join(a.root, "out", s)
    shutil.rmtree(work, ignore_errors=True); os.makedirs(work); os.makedirs(out, exist_ok=True)
    have = [r for r in runs_meta if os.path.exists(os.path.join(raw, f"{r['run_accession']}_1.fastq.gz"))]
    if len(have) < a.min_samples:
        raise RuntimeError(f"so {len(have)} corridas baixadas")

    # 1. regiao por corrida (confirma a_confirmar; descarta ITS etc.)
    info = {}
    for r in have:
        reg, names, n, per = run_signals(raw, r["run_accession"])
        info[r["run_accession"]] = (reg, names, n)
    # corridas sem regiao reconhecida: tenta pela parte 3' dos primers (primer cortado)
    for r in have:
        run = r["run_accession"]
        if info[run][0] in V4_GROUP or info[run][0] == "ITS":
            continue
        hits, tot = primer_suffix_hits(raw, run)
        fr = lambda p: hits[p] / max(tot, 1)
        f515 = fr("515F"); f34 = max(fr("338F"), fr("341F"), fr("341Fb"))
        r806 = max(fr("806R"), fr("805R"))
        reg = None
        if f515 >= 0.3 and f515 >= f34:
            reg = "V4"
        elif f34 >= 0.3:
            reg = "V3-V4"
        elif r806 >= 0.3:
            reg = "V4|V3-V4?"
        if reg:
            nm = info[run][1] + Counter({f"sufixo:{k}": v for k, v in hits.items() if v})
            info[run] = (reg, nm, info[run][2])
    regs = Counter(v[0] for v in info.values() if v[0] in V4_GROUP)
    if not regs:
        raise RuntimeError("nenhuma corrida com V4")
    firm = Counter({k: v for k, v in regs.items() if k in ("V4", "V3-V4")})
    region = (firm or regs).most_common(1)[0][0]
    region = "V4" if region in ("V4", "V4|V3-V4?") else "V3-V4"
    # corridas da regiao majoritaria (V3-V4|V4 e V4|V3-V4? acompanham a maioria)
    ok_reg = {"V4", "V4|V3-V4?"} if region == "V4" else {"V3-V4", "V3-V4|V4", "V4|V3-V4?"}
    use = [r for r in have if info[r["run_accession"]][0] in ok_reg]
    # unidade amostral: BioSample (lanes somadas) ou a corrida, quando o submetedor
    # reaproveitou BioSamples para varias bibliotecas (>4 corridas por BioSample)
    per = Counter(r["sample_accession"] for r in use)
    unit_run = (max(per.values()) if per else 0) > 4
    unit = (lambda r: r["run_accession"]) if unit_run else \
           (lambda r: r["sample_accession"] or r["run_accession"])
    n_samples = len({unit(r) for r in use})
    with open(os.path.join(out, "runs_region.tsv"), "w") as fh:
        fh.write("run\tregion\tsignals\tusada\n")
        for run, (reg, nm, n) in info.items():
            fh.write(f"{run}\t{reg}\t{';'.join(f'{k}:{v}' for k, v in nm.most_common(5))}\t"
                     f"{int(any(run == u['run_accession'] for u in use))}\n")
    if n_samples < a.min_samples:
        raise RuntimeError(f"so {n_samples} amostras {region} apos confirmacao")

    # 2. primers do estudo e cutadapt
    fwd, rev = pick_primers([info[r["run_accession"]][1] for r in use],
                            sum(info[r["run_accession"]][2] for r in use))
    if not (fwd and rev):           # primers parcialmente cortados no deposito?
        sf, sr = pick_primers_suffix([primer_suffix_hits(raw, r["run_accession"])
                                      for r in use[:10]])
        fwd = fwd or sf; rev = rev or sr
    mixed = {r["run_accession"]: orientation_mixed(raw, r["run_accession"], fwd, rev) for r in use}
    with ThreadPoolExecutor(max(1, a.threads // 2)) as ex:
        list(ex.map(lambda r: cutadapt_run(raw, work, r["run_accession"], fwd, rev, 2,
                                           mixed[r["run_accession"]]), use))

    # 3. DADA2
    d2 = os.path.join(work, "dada2")
    p = subprocess.run(["Rscript", os.path.join(HERE, "dada2_study.R"), work, d2, region,
                        str(a.threads)], capture_output=True, text=True)
    with open(os.path.join(out, "dada2_stderr.log"), "w") as fh:   # para diagnostico
        fh.write(p.stdout + "\n" + p.stderr)
    if p.returncode != 0:
        err = open(os.path.join(d2, "ERRO")).read().strip() if os.path.exists(os.path.join(d2, "ERRO")) \
            else (p.stderr.strip().splitlines() or ["?"])[-1]
        raise RuntimeError(f"DADA2: {err}")
    for f in ("asv_long.tsv", "track.tsv", "params.tsv", "seqtab_nochim.rds"):
        shutil.copy(os.path.join(d2, f), out)

    # modo do DADA2: PE (pares mesclados) ou R1 (so leitura forward, sem sobreposicao)
    with open(os.path.join(out, "params.tsv")) as fh:
        mode = next(csv.DictReader(fh, delimiter="\t")).get("mode", "PE")

    # 4. V4 e colapso por amostra (corridas da mesma amostra sao somadas)
    run2samp = {r["run_accession"]: unit(r) for r in use}
    v4 = Counter(); tot = Counter(); lost = Counter()
    with open(os.path.join(out, "asv_long.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c = int(r["count"]); tot[r["run"]] += c
            v = to_v4(r["seq"], region) if mode == "PE" else \
                (r["seq"] if len(r["seq"]) >= 90 else None)   # fragmento R1, usado inteiro
            if v is None:
                lost[r["run"]] += c; continue
            v4[(run2samp[r["run"]], v)] += c
    with gzip.open(os.path.join(out, "v4_counts.tsv.gz"), "wt") as fh:
        fh.write("sample\tseq_v4\tcount\n")
        for (smp, v), c in v4.items():
            fh.write(f"{smp}\t{v}\t{c}\n")
    st = {"study": s, "region": region, "modo": mode,
          "fragmento": "V4 completa" if mode == "PE" else
                       ("inicio da V4 (R1)" if region == "V4" else "V3 (R1)"),
          "fwd": fwd, "rev": rev, "runs_baixadas": len(have),
          "runs_usadas": len(use), "amostras": n_samples,
          "runs_orientacao_mista": sum(mixed.values()),
          "unidade": "corrida" if unit_run else "biosample",
          "leituras_asv": sum(tot.values()), "frac_perdida_v4": sum(lost.values()) / max(sum(tot.values()), 1),
          "v4_unicas": len({v for _, v in v4}),
          "regioes_corridas": dict(Counter(v[0] for v in info.values()))}
    with open(os.path.join(out, "stats.json"), "w") as fh:
        json.dump(st, fh, indent=1)
    shutil.rmtree(work, ignore_errors=True)
    if a.delete_raw:
        shutil.rmtree(raw, ignore_errors=True)
    return st


def cmd_process(a):
    by = load_runs(a.runs)
    order = study_order(by, a.studies)
    sd = os.path.join(a.root, "state"); os.makedirs(sd, exist_ok=True)
    idle = 0
    while True:
        pending = [s for s in order
                   if os.path.exists(os.path.join(sd, f"{s}.downloaded"))
                   and not os.path.exists(os.path.join(sd, f"{s}.done"))
                   and not os.path.exists(os.path.join(sd, f"{s}.failed"))]
        if not pending:
            remaining = [s for s in order if not os.path.exists(os.path.join(sd, f"{s}.downloaded"))]
            if not remaining or not a.wait:
                log("Nada mais a processar."); return
            if idle == 0:
                log(f"Aguardando downloads ({len(remaining)} estudos faltando)...")
            idle += 1; time.sleep(60); continue
        idle = 0
        s = pending[0]
        t0 = time.time()
        try:
            st = process_study(a, s, by[s])
            open(os.path.join(sd, f"{s}.done"), "w").write(json.dumps(st))
            log(f"{s}: OK {st['region']}{'' if st['modo'] == 'PE' else ' [so R1]'} {st['fwd']}/{st['rev']} | {st['amostras']} amostras | "
                f"{st['v4_unicas']} V4 unicas | perda V4 {st['frac_perdida_v4']:.1%} | "
                f"{(time.time() - t0) / 60:.1f} min")
        except Exception as e:
            open(os.path.join(sd, f"{s}.failed"), "w").write(str(e))
            shutil.rmtree(os.path.join(a.work or a.root, "work", s), ignore_errors=True)
            log(f"{s}: FALHOU ({e})")


def extract_v4_table(out, region, mode, run2samp):
    v4 = Counter(); tot = Counter(); lost = Counter()
    with open(os.path.join(out, "asv_long.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["run"] not in run2samp:
                continue
            c = int(r["count"]); tot[r["run"]] += c
            v = to_v4(r["seq"], region) if mode == "PE" else \
                (r["seq"] if len(r["seq"]) >= 90 else None)
            if v is None:
                lost[r["run"]] += c; continue
            v4[(run2samp[r["run"]], v)] += c
    with gzip.open(os.path.join(out, "v4_counts.tsv.gz"), "wt") as fh:
        fh.write("sample\tseq_v4\tcount\n")
        for (smp, v), c in v4.items():
            fh.write(f"{smp}\t{v}\t{c}\n")
    return v4, tot, lost


def cmd_reextract(a):
    """Refaz so a extracao da V4 (a partir do asv_long.tsv), sem DADA2 nem download."""
    by = load_runs(a.runs)
    for s in a.only_studies.split(","):
        out = os.path.join(a.root, "out", s)
        st = json.load(open(os.path.join(out, "stats.json")))
        used = {l.split("\t")[0] for l in open(os.path.join(out, "runs_region.tsv"))
                if l.rstrip("\n").endswith("\t1")}
        use = [r for r in by[s] if r["run_accession"] in used]
        per = Counter(r["sample_accession"] for r in use)
        unit_run = st.get("unidade") == "corrida"
        run2samp = {r["run_accession"]: (r["run_accession"] if unit_run else
                    (r["sample_accession"] or r["run_accession"])) for r in use}
        v4, tot, lost = extract_v4_table(out, st["region"], st.get("modo", "PE"), run2samp)
        st["frac_perdida_v4"] = sum(lost.values()) / max(sum(tot.values()), 1)
        st["v4_unicas"] = len({v for _, v in v4})
        json.dump(st, open(os.path.join(out, "stats.json"), "w"), indent=1)
        log(f"{s}: perda V4 {st['frac_perdida_v4']:.1%} | {st['v4_unicas']} V4 unicas")


def cmd_status(a):
    by = load_runs(a.runs)
    order = study_order(by, a.studies)
    sd = os.path.join(a.root, "state")
    st = Counter()
    fails = Counter()
    for s in order:
        for k in ("done", "failed", "downloaded"):
            if os.path.exists(os.path.join(sd, f"{s}.{k}")):
                st[k] += 1
                if k == "failed":
                    fails[re.sub(r"\d+", "N", open(os.path.join(sd, f"{s}.failed")).read())[:60]] += 1
                break
        else:
            st["pendente"] += 1
    print(f"estudos: {len(order)} | processados: {st['done']} | falharam: {st['failed']} | "
          f"baixados aguardando: {st['downloaded']} | a baixar: {st['pendente']}")
    print(f"espaco livre: {free_gb(a.root):.0f} GB")
    for k, v in fails.most_common(8):
        print(f"  falha ({v}): {k}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["download", "process", "status", "reextract"])
    ap.add_argument("--only-studies", default="", help="reextract: estudos separados por virgula")
    ap.add_argument("--root", required=True, help="pasta no HD externo")
    ap.add_argument("--runs", required=True, help="runs_download.tsv(.gz)")
    ap.add_argument("--studies", default=None, help="studies_final.tsv (ordem/exclusoes)")
    ap.add_argument("--reads", type=int, default=30000, help="pares por corrida")
    ap.add_argument("--work", default=None,
                    help="pasta temporaria (ex.: no disco interno, para poupar o HD externo)")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("--min-samples", type=int, default=20)
    ap.add_argument("--min-free-gb", type=float, default=50)
    ap.add_argument("--limit", type=int, default=0, help="so os N primeiros estudos (teste)")
    ap.add_argument("--wait", action="store_true", help="process: espera novos downloads")
    ap.add_argument("--delete-raw", action="store_true", help="apaga leituras apos processar")
    ap.add_argument("--max-pending", type=int, default=0,
                    help="download: pausa quando ha N estudos baixados aguardando processamento")
    a = ap.parse_args()
    a.root = os.path.abspath(a.root)
    if a.work:
        a.work = os.path.abspath(a.work)
    os.makedirs(a.root, exist_ok=True)
    {"download": cmd_download, "process": cmd_process, "status": cmd_status,
     "reextract": cmd_reextract}[a.mode](a)


if __name__ == "__main__":
    main()
