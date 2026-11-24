#!/usr/bin/env python3
# [S5.1] Representative genomes of GTDB species clusters
# Supplementary 5 of the manuscript. See README.md for the full map.
"""
Triangulacao funcional, etapa 1: genomas representativos dos generos da shortlist.

1. Baixa a tabela de clusters de especie do GTDB (sp_clusters_r<release>.tsv), ou usa --sp-clusters
2. Para cada genero da shortlist (nome GTDB, ex. Bacillus_A): ate --max-sp representantes de especie,
   priorizando especies com mais genomas no cluster (proxy de especie bem estudada/cultivada)
3. Baixa proteinas (.faa), genoma (.fna) e anotacao (.gff) pelo NCBI datasets

Requer o NCBI datasets CLI no PATH (ver mensagem de erro para instalar).
Saidas: <root>/genomas/<genero>/<acesso>.{faa,fna,gff} e <root>/genomas/selecionados.tsv
Uso: python3 fetch_genomes.py --root ~/syncom [--max-sp 5] [--release 226]
"""
import argparse, gzip, os, shutil, subprocess, sys, urllib.request, zipfile, glob
import pandas as pd

BASES = ["https://data.gtdb.ecogenomic.org/releases/release{r}/{r}.0/",
         "https://data.gtdb.aau.ecogenomic.org/releases/release{r}/{r}.0/"]
SUBS = ["", "auxillary_files/", "auxiliary_files/"]
NAMES = ["sp_clusters_r{r}.tsv", "sp_clusters_r{r}.tsv.gz", "sp_clusters.tsv"]


def download_sp(release, dest_dir):
    for b in BASES:
        for sd in SUBS:
            for nm in NAMES:
                u = (b + sd + nm).format(r=release)
                f = os.path.join(dest_dir, os.path.basename(u))
                try:
                    urllib.request.urlretrieve(u, f)
                    print(f"Baixado: {u}", flush=True)
                    return f
                except Exception:
                    continue
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--release", default="226")
    ap.add_argument("--sp-clusters", default=None, help="arquivo sp_clusters local (pula o download)")
    ap.add_argument("--max-sp", type=int, default=5)
    ap.add_argument("--generos", default=None, help="so estes generos (separados por virgula)")
    ap.add_argument("--so-listar", action="store_true", help="so gera selecionados.tsv, sem baixar")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    gd = os.path.join(root, "genomas"); os.makedirs(gd, exist_ok=True)

    sl = pd.read_csv(os.path.join(root, "merged", "meta", "consorcios", "shortlist.tsv"), sep="\t")
    genera = sorted({g for _, r in sl.iterrows() for g in [r.centro] + r.parceiros.split(" + ")})
    if a.generos:
        genera = [g.strip() for g in a.generos.split(",")]

    f = a.sp_clusters
    if not f:
        prev = sorted(glob.glob(os.path.join(gd, f"sp_clusters*{a.release}*")))
        f = prev[0] if prev else download_sp(a.release, gd)
    if not f:
        sys.exit(f"Nao achei sp_clusters do release {a.release}. Abra https://data.gtdb.ecogenomic.org/releases/"
                 f"release{a.release}/ no navegador, baixe o arquivo sp_clusters e rode com --sp-clusters <arquivo>")
    op = gzip.open if f.endswith(".gz") else open
    with op(f, "rt") as fh:
        sp = pd.read_csv(fh, sep="\t", low_memory=False)
    col = lambda key: next(c for c in sp.columns if key.lower() in c.lower())
    c_rep, c_tax, c_n = col("Representative"), col("GTDB taxonomy"), col("No. clustered")
    sp["genus"] = sp[c_tax].str.extract(r"g__([^;]+)")[0]
    sp["species"] = sp[c_tax].str.extract(r"s__([^;]+)")[0]
    sp["acc"] = sp[c_rep].str.replace(r"^(RS_|GB_)", "", regex=True)

    sel, faltam = [], []
    for g in genera:
        d = sp[sp.genus == g]
        match = "exato"
        if d.empty:                         # nome mudou entre versoes: mesmo nome-base (sem sufixo _X)
            b = g.split("_")[0]
            d = sp[sp.genus.str.split("_").str[0] == b]
            match = "nome_base:" + ",".join(sorted(d.genus.unique())) if len(d) else ""
        d = d.sort_values(c_n, ascending=False).head(a.max_sp)
        if d.empty:
            faltam.append(g); continue
        for _, r in d.iterrows():
            sel.append({"genus": g, "gtdb_genus": r.genus, "match": match, "species": r.species,
                        "acc": r.acc, "n_genomas_cluster": r[c_n]})
    S = pd.DataFrame(sel)
    fsel = os.path.join(gd, "selecionados.tsv")
    if a.generos and os.path.exists(fsel):
        old = pd.read_csv(fsel, sep="\t")
        S = pd.concat([old[~old.genus.isin(genera)], S], ignore_index=True)
    S.to_csv(fsel, sep="\t", index=False)
    aprox = S[S.match.astype(str).str.startswith("nome_base")].drop_duplicates("genus")
    for _, r in aprox.iterrows():
        print(f"Aviso: {r.genus} nao existe no r{a.release}; usados generos de mesmo nome-base ({r.match[10:]})")
    print(f"Generos: {len(genera)} | com genomas: {S.genus.nunique()} | genomas: {len(S)}")
    if faltam:
        print("Sem cluster no GTDB r" + a.release + " (nome mudou entre versoes?): " + ", ".join(faltam))
    print(S.groupby("genus").size().to_string())
    if a.so_listar:
        return

    if not shutil.which("datasets"):
        sys.exit("\nNCBI datasets nao encontrado. Instale com:\n"
                 "  mkdir -p ~/bin && cd ~/bin && curl -o datasets "
                 "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets "
                 "&& chmod +x datasets && export PATH=~/bin:$PATH\n")
    for g, d in S.groupby("genus"):
        od = os.path.join(gd, g); os.makedirs(od, exist_ok=True)
        todo = [x for x in d.acc if not os.path.exists(os.path.join(od, f"{x}.faa"))]
        if not todo:
            continue
        lst = os.path.join(od, "acc.txt"); open(lst, "w").write("\n".join(todo) + "\n")
        z = os.path.join(od, "tmp.zip")
        print(f"{g}: baixando {len(todo)} genomas", flush=True)
        r = subprocess.run(["datasets", "download", "genome", "accession", "--inputfile", lst,
                            "--include", "protein,genome,gff3", "--filename", z])
        if r.returncode != 0 or not os.path.exists(z):
            print(f"  falhou: {g}"); continue
        with zipfile.ZipFile(z) as zf:
            for n in zf.namelist():
                parts = n.split("/")
                if len(parts) < 4:
                    continue
                acc, fn = parts[-2], parts[-1]
                ext = ("faa" if fn == "protein.faa" else "gff" if fn == "genomic.gff"
                       else "fna" if fn.endswith("_genomic.fna") else None)
                if ext:
                    with zf.open(n) as src, open(os.path.join(od, f"{acc}.{ext}"), "wb") as dst:
                        shutil.copyfileobj(src, dst)
        os.remove(z)
    # genomas sem proteinas no NCBI (alguns GenBank/MAG) ficam para predicao com prodigal
    sem = [os.path.basename(p)[:-4] for p in glob.glob(os.path.join(gd, "*", "*.fna"))
           if not os.path.exists(p[:-4] + ".faa")]
    if sem:
        print(f"\nGenomas sem proteinas no NCBI (prever com prodigal na etapa 2): {len(sem)}")
    print("FIM")


if __name__ == "__main__":
    main()
