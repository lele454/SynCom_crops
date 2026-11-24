#!/usr/bin/env python3
# [S1.10] Taxonomic assignment, SINTAX fallback
# Supplementary 1 of the manuscript. See README.md for the full map.
"""
Taxonomia com SINTAX (vsearch), pouca memoria. Usa a mesma referencia GTDB no
formato do DADA2 (cabecalho "Dominio;Filo;Classe;Ordem;Familia;Genero;").

Uso:
  python3 assign_tax_sintax.py ~/syncom/merged <referencia.assignTaxonomy.fna.gz> 3 [0.8]
Saidas: taxonomy_v4.tsv e taxonomy_r1.tsv (nivel atribuido quando a confianca >= corte;
         as confiancas de cada nivel ficam nas colunas boot_*)
"""
import csv, gzip, os, re, subprocess, sys

LEVELS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
CODES = ["d", "p", "c", "o", "f", "g"]


def convert_ref(src, dst):
    """DADA2 -> SINTAX: >id;tax=d:X,p:Y,...;"""
    if os.path.exists(dst):
        return
    op = gzip.open if src.endswith(".gz") else open
    n = 0
    with op(src, "rt") as fi, open(dst + ".part", "w") as fo:
        for line in fi:
            if line.startswith(">"):
                parts = [p.strip() for p in line[1:].strip().split(";")]
                parts = [p for p in parts if p]
                # referencias como a SBDI trazem Dominio;Reino;Filo... (ex.: Bacteria;Bacteria;...)
                # remove o nivel "Reino" quando repete o dominio
                if len(parts) >= 2 and parts[1] in ("Bacteria", "Archaea", "Eukaryota") \
                        and parts[1] == parts[0]:
                    parts = [parts[0]] + parts[2:]
                parts = parts[:6]
                tax = ",".join(f"{c}:{re.sub(r'[,:;]', '_', p)}" for c, p in zip(CODES, parts))
                n += 1
                fo.write(f">r{n};tax={tax};\n")
            else:
                fo.write(line)
    os.replace(dst + ".part", dst)
    print(f"referencia convertida: {n} sequencias", flush=True)


def run(md, ref_sintax, name, out, threads, cutoff):
    fa = os.path.join(md, name)
    if not os.path.exists(fa) or os.path.getsize(fa) == 0:
        return
    raw = os.path.join(md, out.replace(".tsv", ".sintax.txt"))
    if not os.path.exists(raw):
        subprocess.run(["vsearch", "--sintax", fa, "--db", ref_sintax, "--tabbedout", raw + ".part",
                        "--sintax_cutoff", str(cutoff), "--strand", "plus",
                        "--threads", str(threads), "--quiet"], check=True)
        os.replace(raw + ".part", raw)
    rows, with_genus = [], 0
    for line in open(raw):
        f = line.rstrip("\n").split("\t")
        q = f[0]; full = f[1] if len(f) > 1 else ""
        row = {"seq_id": q}
        for c, lv in zip(CODES, LEVELS):
            m = re.search(rf"(?:^|,){c}:([^,(]+)\(([\d.]+)\)", full)
            name_, conf = (m.group(1), float(m.group(2))) if m else ("", 0.0)
            row[lv] = name_ if conf >= cutoff else ""
            row["boot_" + lv] = round(conf * 100)
        # niveis abaixo de um nivel nao atribuido ficam vazios
        for i, lv in enumerate(LEVELS):
            if not row[lv]:
                for lv2 in LEVELS[i + 1:]:
                    row[lv2] = ""
                break
        with_genus += bool(row["Genus"])
        rows.append(row)
    with open(os.path.join(md, out), "w", newline="") as fh:
        w = csv.DictWriter(fh, ["seq_id"] + LEVELS + ["boot_" + l for l in LEVELS], delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"{out}: {len(rows)} sequencias | com genero: {100 * with_genus / max(len(rows), 1):.1f}%",
          flush=True)


def main():
    md, ref, threads = os.path.expanduser(sys.argv[1]), sys.argv[2], int(sys.argv[3])
    cutoff = float(sys.argv[4]) if len(sys.argv) > 4 else 0.8
    ref_sintax = os.path.join(md, "ref_sintax.fasta")
    convert_ref(ref, ref_sintax)
    run(md, ref_sintax, "seqs_v4.fasta", "taxonomy_v4.tsv", threads, cutoff)
    run(md, ref_sintax, "seqs_r1.fasta", "taxonomy_r1.tsv", threads, cutoff)


if __name__ == "__main__":
    main()
