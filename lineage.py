#!/usr/bin/env python3
# [S1.15] Genus lineage table
# Supplementary 1 of the manuscript. See README.md for the full map.
"""Tabela genero -> linhagem. Fonte primaria: ref_sintax.fasta (a mesma taxonomia que anotou
os dados); complemento: bac120_metadata do GTDB. Saida: <root>/merged/genus_lineage.tsv"""
import argparse, glob, gzip, os, re
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--root", required=True)
ap.add_argument("--metadata", default=None)
a = ap.parse_args()
root = os.path.abspath(os.path.expanduser(a.root)); md = os.path.join(root, "merged")
vis = {}

fa = os.path.join(md, "ref_sintax.fasta")
n_fa = 0
if os.path.exists(fa):
    op = gzip.open if fa.endswith(".gz") else open
    with op(fa, "rt") as fh:
        for ln in fh:
            if not ln.startswith(">"):
                continue
            t = dict(re.findall(r"([dpcofgs]):([^,;\n]+)", ln.split("tax=")[-1]))
            g = t.get("g", "").strip()
            if g and g not in vis:
                vis[g] = (t.get("p", ""), t.get("c", ""), t.get("o", ""), t.get("f", ""))
                n_fa += 1
print(f"ref_sintax.fasta: {n_fa} generos")

f = a.metadata or next(iter(sorted(glob.glob(os.path.join(root, "genomas", "bac120_metadata*")))), None)
n_md = 0
if f:
    for ch in pd.read_csv(f, sep="\t", usecols=["gtdb_taxonomy"], low_memory=False, chunksize=200000):
        for t in ch.gtdb_taxonomy.dropna().unique():
            p = dict(x.split("__", 1) for x in t.split(";") if "__" in x)
            g = p.get("g", "")
            if g and g not in vis:
                vis[g] = (p.get("p", ""), p.get("c", ""), p.get("o", ""), p.get("f", ""))
                n_md += 1
    print(f"{os.path.basename(f)}: +{n_md} generos")

D = pd.DataFrame([(g, *v) for g, v in vis.items()], columns=["genus", "phylum", "class", "order", "family"])
out = os.path.join(md, "genus_lineage.tsv")
D.sort_values("genus").to_csv(out, sep="\t", index=False)
occ = pd.read_csv(os.path.join(md, "genus_occupancy.tsv"), sep="\t")
falta = sorted(set(occ.genus) - set(D.genus))
print(f"total {len(D)} generos | do dataset ainda sem linhagem: {len(falta)}")
if falta:
    print("exemplos:", ", ".join(falta[:10]))
