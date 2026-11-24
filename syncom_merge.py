#!/usr/bin/env python3
# [S1.11] Merging, quality control and aggregation to genus
# Supplementary 1 of the manuscript. See README.md for the full map.
"""
Uniao dos estudos processados, controle de qualidade e agregacao em genero.

Etapas (cada uma retomavel):
  merge  : junta os v4_counts de todos os estudos .done; gera as sequencias unicas
           (V4 completa e fragmentos so-R1, separadas), a tabela longa de contagens
           e o controle de qualidade por amostra e por estudo
  [taxonomia: Rscript assign_tax.R  -> taxonomy_*.tsv]
  genus  : junta a taxonomia, aplica os filtros e grava uma tabela de generos por
           estudo, mais as tabelas para PRISMA e Figs. 1-2

Uso:
  python3 syncom_merge.py merge --root ~/syncom --studies subset_studies.tsv
  Rscript assign_tax.R ~/syncom/merged <referencia_1genome.fna.gz> 6
  python3 syncom_merge.py genus --root ~/syncom --studies subset_studies.tsv
"""
import argparse, csv, glob, gzip, hashlib, json, os, statistics, sys
from collections import Counter, defaultdict

csv.field_size_limit(10 ** 8)
log = lambda *x: print(*x, file=sys.stderr, flush=True)


def sid(seq, mode):
    return ("v4_" if mode == "PE" else "r1_") + hashlib.sha1(seq.encode()).hexdigest()[:14]


def load_studies(path):
    with open(path) as fh:
        return {r["study_accession"]: r for r in csv.DictReader(fh, delimiter="\t")}


# ------------------------------------------------------------------- merge ---
def cmd_merge(a):
    root, md = a.root, os.path.join(a.root, "merged")
    os.makedirs(md, exist_ok=True)
    done = sorted(os.path.basename(f)[:-5] for f in glob.glob(os.path.join(root, "state", "*.done")))
    meta = load_studies(a.studies)
    done = [s for s in done if s in meta]
    log(f"Estudos concluidos no recorte: {len(done)}")

    seq_info = {}                         # id -> [seq, modo, total, n_amostras, n_estudos]
    seq_studies = defaultdict(set)
    sample_rows, study_rows = [], []
    with gzip.open(os.path.join(md, "counts_long.tsv.gz"), "wt") as fo:
        fo.write("study\tsample\tseq_id\tcount\n")
        for s in done:
            od = os.path.join(root, "out", s)
            st = json.load(open(os.path.join(od, "stats.json")))
            mode = st.get("modo", "PE")
            per_sample = defaultdict(Counter)
            with gzip.open(os.path.join(od, "v4_counts.tsv.gz"), "rt") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    k = sid(r["seq_v4"], mode); c = int(r["count"])
                    per_sample[r["sample"]][k] += c
                    if k not in seq_info:
                        seq_info[k] = [r["seq_v4"], mode, 0, 0]
                    seq_info[k][2] += c; seq_info[k][3] += 1
                    seq_studies[k].add(s)
            for smp, cnt in per_sample.items():
                for k, c in cnt.items():
                    fo.write(f"{s}\t{smp}\t{k}\t{c}\n")
                sample_rows.append({"study": s, "sample": smp, "reads": sum(cnt.values()),
                                    "riqueza": len(cnt)})
            # track: retencao media do estudo
            tr = list(csv.DictReader(open(os.path.join(od, "track.tsv")), delimiter="\t"))
            inp = sum(float(t["input"]) for t in tr if t["input"] not in ("", "NA"))
            fin = sum(float(t["nonchim"]) for t in tr if t["nonchim"] not in ("", "NA"))
            pr = next(csv.DictReader(open(os.path.join(od, "params.tsv")), delimiter="\t"))
            study_rows.append({"study": s, "region": st["region"], "modo": mode,
                               "fragmento": st.get("fragmento", ""), "unidade": st.get("unidade", ""),
                               "amostras": len(per_sample), "retencao": round(fin / max(inp, 1), 3),
                               "maxEE": pr.get("maxEE", ""), "truncF": pr.get("truncF", ""),
                               "truncR": pr.get("truncR", ""), "frac_quimeras": pr.get("frac_chim", "")})
            log(f"  {s}: {len(per_sample)} amostras, modo {mode}")

    # sequencias unicas, separadas por modo; 'classificar' = vale a pena atribuir taxonomia
    n_cls = Counter()
    for mode, fn in (("PE", "seqs_v4.fasta"), ("R1", "seqs_r1.fasta")):
        with open(os.path.join(md, fn), "w") as fh:
            for k, (seq, m, tot, ns) in seq_info.items():
                if m != mode:
                    continue
                if ns >= a.min_samples_seq or tot >= a.min_count_seq:
                    fh.write(f">{k}\n{seq}\n"); n_cls[mode] += 1
    with gzip.open(os.path.join(md, "seq_info.tsv.gz"), "wt") as fh:
        fh.write("seq_id\tmodo\ttotal\tn_amostras\tn_estudos\tclassificada\tseq\n")
        for k, (seq, m, tot, ns) in seq_info.items():
            cl = int(ns >= a.min_samples_seq or tot >= a.min_count_seq)
            fh.write(f"{k}\t{m}\t{tot}\t{ns}\t{len(seq_studies[k])}\t{cl}\t{seq}\n")

    tot_reads = sum(v[2] for v in seq_info.values())
    cov = sum(v[2] for v in seq_info.values() if v[3] >= a.min_samples_seq or v[2] >= a.min_count_seq)
    write_tsv(os.path.join(md, "sample_qc.tsv"), sample_rows)
    write_tsv(os.path.join(md, "study_processing.tsv"), study_rows)
    log(f"\nSequencias unicas: {len(seq_info)} | a classificar: V4 {n_cls['PE']}, R1 {n_cls['R1']} "
        f"({cov / max(tot_reads, 1):.1%} das leituras)")
    log("Proximo passo: Rscript assign_tax.R <root>/merged <referencia.fna.gz> <threads>")


def write_tsv(path, rows, fields=None):
    fields = list(rows[0]) if rows else (fields or [])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fields, delimiter="\t"); w.writeheader(); w.writerows(rows)


# ------------------------------------------------------------------- genus ---
def cmd_genus(a):
    md = os.path.join(a.root, "merged")
    meta = load_studies(a.studies)
    tax = {}
    for fn in ("taxonomy_v4.tsv", "taxonomy_r1.tsv"):
        p = os.path.join(md, fn)
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                tax[r["seq_id"]] = r
    log(f"Sequencias com taxonomia: {len(tax)}")

    # QC por amostra e estudo
    samples = list(csv.DictReader(open(os.path.join(md, "sample_qc.tsv")), delimiter="\t"))
    proc = {r["study"]: r for r in csv.DictReader(open(os.path.join(md, "study_processing.tsv")),
                                                    delimiter="\t")}
    ok_sample = {(r["study"], r["sample"]) for r in samples if int(r["reads"]) >= a.min_depth}
    by_study = defaultdict(list)
    for r in samples:
        by_study[r["study"]].append(r)

    # tabelas de genero por estudo
    sufixo = getattr(a, "sufixo", "") or ""
    gdir = os.path.join(md, "genus" + sufixo); os.makedirs(gdir, exist_ok=True)
    gtab = defaultdict(lambda: defaultdict(Counter))      # study -> sample -> genus -> count
    nonbact = Counter()
    with gzip.open(os.path.join(md, "counts_long.tsv.gz"), "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r["study"], r["sample"]) not in ok_sample:
                continue
            t = tax.get(r["seq_id"])
            c = int(r["count"])
            if t and t.get("Kingdom", "") not in ("Bacteria", "d__Bacteria"):
                if t.get("Kingdom"):
                    nonbact[r["study"]] += c
                    continue
            g = (t or {}).get("Genus", "") or "NA"
            if g == "NA":
                fam = (t or {}).get("Family", "")
                g = f"NA|{fam}" if fam else "NA"
            gtab[r["study"]][r["sample"]][g] += c

    qc_rows, prisma = [], []
    for s, rows in by_study.items():
        n_tot = len(rows)
        kept = [r for r in rows if (s, r["sample"]) in ok_sample]
        rich = statistics.median(int(r["riqueza"]) for r in kept) if kept else 0
        # fracao das leituras com genero atribuido
        tot = sum(sum(c.values()) for c in gtab[s].values())
        assigned = sum(v for c in gtab[s].values() for g, v in c.items() if not g.startswith("NA"))
        frac_gen = assigned / max(tot, 1)
        motivo = ""
        if len(kept) < a.min_samples:
            motivo = f"<{a.min_samples} amostras com >= {a.min_depth} leituras"
        elif rich < a.min_richness:
            motivo = f"riqueza mediana {rich:.0f} < {a.min_richness} (comunidade simplificada?)"
        elif frac_gen < a.min_frac_genus:
            motivo = f"so {frac_gen:.0%} das leituras com genero"
        p = proc.get(s, {})
        m = meta.get(s, {})
        qc_rows.append({"study": s, "cultura": m.get("cultura", ""), "grupo": m.get("grupo", ""),
                        "country": (m.get("country", "") or "").split(":")[0],
                        "region": p.get("region", ""), "modo": p.get("modo", ""),
                        "fragmento": p.get("fragmento", ""), "amostras_proc": n_tot,
                        "amostras_ok": len(kept), "riqueza_mediana": round(rich),
                        "frac_genero": round(frac_gen, 3), "retencao_dada2": p.get("retencao", ""),
                        "status": "incluido" if not motivo else "excluido", "motivo": motivo})
        if not motivo:
            genera = sorted({g for c in gtab[s].values() for g in c})
            with gzip.open(os.path.join(gdir, f"{s}.tsv.gz"), "wt") as fo:
                fo.write("sample\t" + "\t".join(genera) + "\n")
                for smp in sorted(gtab[s]):
                    fo.write(smp + "\t" + "\t".join(str(gtab[s][smp].get(g, 0)) for g in genera) + "\n")
    write_tsv(os.path.join(md, "study_qc" + sufixo + ".tsv"), sorted(qc_rows, key=lambda r: r["study"]))

    # PRISMA: estudos do recorte que nao chegaram ao fim do processamento
    for s in meta:
        f = os.path.join(a.root, "state", f"{s}.failed")
        if os.path.exists(f):
            prisma.append({"study": s, "etapa": "processamento", "motivo": open(f).read().strip()})
        elif s not in proc:
            prisma.append({"study": s, "etapa": "processamento", "motivo": "nao processado"})
    for r in qc_rows:
        if r["status"] == "excluido":
            prisma.append({"study": r["study"], "etapa": "controle de qualidade", "motivo": r["motivo"]})
    write_tsv(os.path.join(md, "prisma_exclusoes" + sufixo + ".tsv"), prisma, ["study", "etapa", "motivo"])

    # ocupacao x abundancia (Fig. S: nucleo)
    inc = [r["study"] for r in qc_rows if r["status"] == "incluido"]
    occ, ab = Counter(), defaultdict(list)
    for s in inc:
        rel = defaultdict(list)
        for smp, c in gtab[s].items():
            t = sum(c.values())
            for g, v in c.items():
                rel[g].append(v / t)
        n = len(gtab[s])
        for g, vals in rel.items():
            if g.startswith("NA"):
                continue
            occ[g] += 1
            ab[g].append(sum(vals) / n)
    write_tsv(os.path.join(md, "genus_occupancy" + sufixo + ".tsv"),
              sorted(({"genus": g, "n_estudos": occ[g], "ocupacao": round(occ[g] / len(inc), 4),
                       "abund_media": f"{statistics.mean(ab[g]):.3e}"} for g in occ),
                     key=lambda r: -r["n_estudos"]))
    ex = len(qc_rows) - len(inc)
    log(f"Estudos incluidos: {len(inc)} | excluidos no controle de qualidade: {ex}")
    log(f"Amostras incluidas: {sum(r['amostras_ok'] for r in qc_rows if r['status'] == 'incluido')}")
    log(f"Generos em >=50% dos estudos: {sum(1 for g in occ if occ[g] >= len(inc) / 2)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["merge", "genus"])
    ap.add_argument("--root", required=True)
    ap.add_argument("--studies", required=True, help="subset_studies.tsv")
    ap.add_argument("--min-samples-seq", type=int, default=2,
                    help="classificar sequencias presentes em >= N amostras...")
    ap.add_argument("--min-count-seq", type=int, default=50, help="...ou com >= N leituras")
    ap.add_argument("--min-depth", type=int, default=2000, help="leituras por amostra")
    ap.add_argument("--min-samples", type=int, default=20)
    ap.add_argument("--min-richness", type=int, default=50, help="riqueza mediana por amostra")
    ap.add_argument("--sufixo", default="",
                    help="genus: sufixo das saidas (ex.: _d5000 -> merged/genus_d5000, study_qc_d5000.tsv)")
    ap.add_argument("--min-frac-genus", type=float, default=0.3)
    a = ap.parse_args()
    a.root = os.path.abspath(os.path.expanduser(a.root))
    {"merge": cmd_merge, "genus": cmd_genus}[a.mode](a)


if __name__ == "__main__":
    main()
