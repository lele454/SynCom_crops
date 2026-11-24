#!/usr/bin/env python3
# [S3.11] Phylum representation and mediation of the deviation
# Supplementary 3 of the manuscript. See README.md for the full map.
"""
Vies por filo nas redes de coocorrencia de 16S, e teste do mecanismo (copias do 16S).

Achado anterior: Bacillota formam menos associacoes consistentes que generos de abundancia e
ocupacao equivalentes, com ou sem endosporo. Aqui: (a) quantificar isso para todos os filos;
(b) testar se o numero de copias do 16S explica o efeito.

Traços genomicos por genero, do metadata do GTDB (bac120_metadata):
  ssu_count       copias do 16S, SO de genomas completos (em rascunhos e MAGs o rRNA nao
                  monta e a contagem e subestimada) - este e o preditor do mecanismo
  genome_size, gc tamanho do genoma e GC (medianas de genomas de alta qualidade)
Modelos de Poisson quasi, desfecho = parceiros com associacao positiva (FDR < 0,05),
offset = log(pares testados):
  A  base      ~ log_abund + ocupacao                 -> residuo por filo = ranking de filos
  B  mecanismo ~ base + ssu_count + genome_size + gc
  C  filo      ~ B + dummy de Bacillota (e dos filos maiores)
  Mediacao: atenuacao do coeficiente de Bacillota de A+dummy para B+dummy. Se o efeito de
  Bacillota cair para perto de zero ao entrar ssu_count, o mecanismo e o numero de copias.

Saidas: <root>/merged/meta/vies_filos/ (por_genero.tsv, por_filo.tsv, resumo.txt)
Uso: LD_PRELOAD=... python3 phylum_bias.py --root ~/syncom [--metadata arquivo_ja_baixado.tsv.gz]
"""
import argparse, os, re, sys, urllib.request
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inoculant_bias import irls
from spore_bias import tax_map

BASES = ["https://data.gtdb.ecogenomic.org/releases/release{r}/{r}.0/",
         "https://data.gtdb.aau.ecogenomic.org/releases/release{r}/{r}.0/"]
NAMES = ["bac120_metadata_r{r}.tsv.gz", "bac120_metadata.tsv.gz", "bac120_metadata_r{r}.tsv.gz.1"]
COLS = ["gtdb_taxonomy", "ssu_count", "genome_size", "gc_percentage",
        "checkm_completeness", "checkm_contamination", "ncbi_assembly_level"]


def baixar(release, dest):
    for b in BASES:
        for nm in NAMES:
            u = (b + nm).format(r=release)
            f = os.path.join(dest, os.path.basename(u))
            if os.path.exists(f):
                return f
            try:
                print(f"Baixando {u} (algumas centenas de MB)...", flush=True)
                urllib.request.urlretrieve(u, f)
                return f
            except Exception:
                continue
    return None


def traits(path):
    """genero -> copias do 16S (genomas completos), tamanho do genoma e GC (alta qualidade)."""
    use = None
    keep = []
    for ch in pd.read_csv(path, sep="\t", low_memory=False, chunksize=200000):
        if use is None:
            use = [c for c in COLS if c in ch.columns]
            falta = set(COLS) - set(use)
            if falta:
                print("colunas ausentes no metadata: " + ", ".join(sorted(falta)), flush=True)
        keep.append(ch[use])
    d = pd.concat(keep, ignore_index=True)
    d["genus"] = d.gtdb_taxonomy.str.extract(r"g__([^;]+)")[0]
    d = d.dropna(subset=["genus"])
    for c in ("ssu_count", "genome_size", "gc_percentage", "checkm_completeness", "checkm_contamination"):
        if c in d:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    hq = d
    if "checkm_completeness" in d:
        hq = d[(d.checkm_completeness >= 90) & (d.checkm_contamination <= 5)]
    g = hq.groupby("genus").agg(genome_size=("genome_size", "median"), gc=("gc_percentage", "median"),
                                n_genomas_hq=("genus", "size"))
    comp = d
    if "ncbi_assembly_level" in d:
        comp = d[d.ncbi_assembly_level.astype(str).str.contains("Complete", case=False, na=False)]
    s = comp.groupby("genus").agg(ssu_count=("ssu_count", "median"), n_completos=("genus", "size"))
    return g.join(s, how="outer")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--release", default="226")
    ap.add_argument("--metadata", default=None)
    ap.add_argument("--min-filo", type=int, default=8, help="minimo de generos para listar o filo")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, "vies_filos"); os.makedirs(out, exist_ok=True)
    gd = os.path.join(root, "genomas"); os.makedirs(gd, exist_ok=True)
    f = a.metadata or baixar(a.release, gd)
    if not f:
        raise SystemExit(f"metadata do GTDB nao encontrado. Baixe bac120_metadata de "
                         f"https://data.gtdb.ecogenomic.org/releases/release{a.release}/ e use --metadata")

    m = pd.read_csv(os.path.join(mm, "meta_pairs.tsv.gz"), sep="\t",
                    usecols=["genus_a", "genus_b", "mu", "fdr"])
    pos = m[(m.mu > 0) & (m.fdr < 0.05)]
    cnt = lambda d: pd.concat([d.genus_a, d.genus_b]).value_counts()
    G = pd.DataFrame({"n_testados": cnt(m), "grau_pos": cnt(pos)}).fillna(0)
    G = G[G.n_testados >= 10]
    occ = pd.read_csv(os.path.join(md, "genus_occupancy.tsv"), sep="\t").set_index("genus")
    G = G.join(occ[["ocupacao", "abund_media"]]).join(tax_map(md)).join(traits(f))
    G = G.dropna(subset=["ocupacao", "abund_media"])
    G["log_abund"] = np.log10(G.abund_media.clip(lower=1e-8))
    G["log_genoma"] = np.log10(G.genome_size.clip(lower=1e5))
    G["taxa_pos"] = G.grau_pos / G.n_testados
    G["bacillota"] = G.filo.fillna("").str.match(r"^(Bacillota|Firmicutes)").astype(int)
    y = G.grau_pos.values.astype(float); off = G.n_testados.values

    L = [f"Generos testados (>= 10 pares): {len(G)} | com filo atribuido: {G.filo.notna().sum()}",
         f"Com copias do 16S de genomas completos: {G.ssu_count.notna().sum()} | "
         f"com tamanho de genoma: {G.genome_size.notna().sum()}", ""]

    # A. base e ranking de filos
    base_cols = ["log_abund", "ocupacao"]
    tabA, phiA = irls(y, G[base_cols].values, off, True, base_cols)
    L += [f"A. Modelo base (dispersao quasi = {phiA:.1f}):",
          tabA.to_string(float_format=lambda v: f"{v:.3g}"), ""]
    b = tabA.coef.values
    G["esperado"] = np.exp(b[0] + G[base_cols].values @ b[1:]) * off
    P = G.dropna(subset=["filo"]).groupby("filo").apply(lambda d: pd.Series(
        {"n_generos": len(d), "obs": d.grau_pos.sum(), "esp": d.esperado.sum(),
         "ssu_mediano": d.ssu_count.median(), "abund_mediana": d.abund_media.median()}),
        include_groups=False)
    P["razao_obs_esp"] = P.obs / P.esp
    P["p_poisson"] = [2 * min(stats.poisson.cdf(o, e), 1 - stats.poisson.cdf(o - 1, e))
                      for o, e in zip(P.obs, P.esp)]
    P = P[P.n_generos >= a.min_filo].sort_values("razao_obs_esp")
    P.to_csv(os.path.join(out, "por_filo.tsv"), sep="\t", float_format="%.4g")
    L += [f"Ranking de filos (observado/esperado pelo modelo base; filos com >= {a.min_filo} generos):",
          P.to_string(float_format=lambda v: f"{v:.3g}"), ""]

    # B/C. mecanismo e mediacao, so nos generos com os tracos
    H = G.dropna(subset=["ssu_count", "genome_size", "gc"])
    if len(H) < 40:
        L.append(f"Mecanismo: so {len(H)} generos com tracos genomicos; analise nao aplicavel")
    else:
        yh = H.grau_pos.values.astype(float); offh = H.n_testados.values
        cb = base_cols + ["bacillota"]
        t1, p1 = irls(yh, H[cb].values, offh, True, cb)
        cm = base_cols + ["ssu_count", "log_genoma", "gc"]
        t2, p2 = irls(yh, H[cm].values, offh, True, cm)
        cf = cm + ["bacillota"]
        t3, p3 = irls(yh, H[cf].values, offh, True, cf)
        L += [f"(generos com tracos genomicos: {len(H)})",
              f"B1. Base + Bacillota (quasi = {p1:.1f}):", t1.to_string(float_format=lambda v: f"{v:.3g}"), "",
              f"B2. Base + copias do 16S, tamanho do genoma e GC (quasi = {p2:.1f}):",
              t2.to_string(float_format=lambda v: f"{v:.3g}"), "",
              f"B3. B2 + Bacillota (quasi = {p3:.1f}):", t3.to_string(float_format=lambda v: f"{v:.3g}"), ""]
        c1, c3 = t1.loc["bacillota", "coef"], t3.loc["bacillota", "coef"]
        L.append(f"Mediacao: coeficiente de Bacillota {c1:+.3f} -> {c3:+.3f} ao entrar copias do 16S, "
                 f"tamanho e GC (atenuacao de {100 * (1 - c3 / c1) if c1 else float('nan'):.0f}%)")
        L.append("   atenuacao perto de 100% = o traco explica o efeito do filo; perto de 0% = nao explica")
        r = stats.spearmanr(H.ssu_count, H.taxa_pos)
        L.append(f"   correlacao bruta copias do 16S x taxa de associacao: rho = {r.statistic:.3f}, "
                 f"p = {r.pvalue:.3g} (n = {len(H)})")
        L.append(f"   copias do 16S (mediana): Bacillota {H.loc[H.bacillota == 1, 'ssu_count'].median():.1f} | "
                 f"demais {H.loc[H.bacillota == 0, 'ssu_count'].median():.1f}")
        # D. mediacao variavel por variavel
        H2 = H.copy()
        H2["fatia_gtdb"] = pd.Series(H2.index, index=H2.index).str.contains(r"_[A-Z]+$").astype(int)
        cands = {"copias do 16S": ["ssu_count"],
                 "tamanho do genoma": ["log_genoma"],
                 "conteudo GC": ["gc"],
                 "divisao de genero no GTDB": ["fatia_gtdb"],
                 "as tres variaveis genomicas juntas": ["ssu_count", "log_genoma", "gc"]}
        L += ["", "D. Mediacao do efeito de Bacillota, variavel por variavel",
              "   (atenuacao = quanto do coeficiente do filo a variavel absorve)"]
        for nome, cc in cands.items():
            col = base_cols + cc + ["bacillota"]
            try:
                tb, pb = irls(yh, H2[col].values, offh, True, col)
                cx = tb.loc["bacillota", "coef"]
                at = 100 * (1 - cx / c1) if c1 else float("nan")
                L.append(f"   {nome}: {c1:+.3f} -> {cx:+.3f} | atenuacao {at:.0f}% | quasi {pb:.1f}")
            except Exception as ex:
                L.append(f"   {nome}: nao estimavel ({type(ex).__name__}: {ex})")
        n_fatia = int(H2.fatia_gtdb.sum())
        L.append(f"   generos com sufixo de fatia no GTDB: {n_fatia} de {len(H2)} "
                 f"| entre Bacillota: {int(H2.loc[H2.bacillota == 1, 'fatia_gtdb'].sum())}")
        L.append("   endosporos: nao estimavel por mediacao, pois a formacao de endosporo e "
                 "praticamente coextensiva a Bacillota (indicador colinear com o filo)")
    G.to_csv(os.path.join(out, "por_genero.tsv"), sep="\t", float_format="%.5g")
    txt = "\n".join(L)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
