#!/usr/bin/env python3
# [S4.1] Combinations from cliques of confirmed positive associations
# Supplementary 4 of the manuscript. See README.md for the full map.
"""
Consorcios candidatos (SynComs) a partir do nucleo robusto.

Para cada cenario (tier 1; tier 1+2):
  1. rede positiva = pares positivos do nucleo nos tiers do cenario
  2. enumera cliques de 2..max_size generos (todos os pares internos positivos)
  3. descarta grupos com qualquer par interno negativo:
       - consistente_neg (qualquer tier), e
       - associado negativo (FDR < 0,05), salvo --permissivo
  4. mantem grupos com >= 1 genero PGP (lista ampliada, ver PGP_* abaixo)
  5. metricas: elo mais fraco (min r), r medio, min limite inf. do PI, I2 max,
     copresenca (fracao de amostras com todos os membros detectados),
     estudos e culturas onde o grupo esta copresente em >= 50% das amostras
  6. ranking: elo mais fraco, depois copresenca

Saidas: <root>/merged/meta/consorcios/
  consorcios_tier1.tsv, consorcios_tier12.tsv, resumo.txt, lista_pgp.tsv
Uso: LD_PRELOAD=... python3 consortia.py --root ~/syncom
"""
import argparse, gzip, itertools, os, re
from collections import defaultdict
import numpy as np, pandas as pd, networkx as nx

# PGP classicos (mesma lista de robust_core.py / make_figures.py)
PGP_CLASSIC = ["Pseudomonas", "Bacillus", "Paenibacillus", "Streptomyces", "Bradyrhizobium", "Rhizobium",
               "Mesorhizobium", "Sinorhizobium", "Ensifer", "Azospirillum", "Herbaspirillum", "Burkholderia",
               "Paraburkholderia", "Enterobacter", "Pantoea", "Serratia", "Klebsiella", "Azotobacter",
               "Gluconacetobacter", "Arthrobacter", "Microbacterium", "Variovorax", "Brucella",
               "Stenotrophomonas", "Lysobacter", "Flavobacterium", "Achromobacter", "Delftia",
               "Acinetobacter", "Kosakonia", "Priestia", "Peribacillus", "Lysinibacillus", "Rhodococcus",
               "Agrobacterium"]
# Ampliacao: generos com linhagens PGP descritas em rizosfera (conferir referencias antes do artigo)
PGP_EXTRA = ["Chryseobacterium", "Sphingobacterium", "Sphingomonas", "Sphingopyxis", "Novosphingobium",
             "Chryseomonas", "Massilia", "Duganella", "Cupriavidus", "Azoarcus", "Azovibrio",
             "Rhizobacter", "Mitsuaria", "Pedobacter", "Caulobacter", "Phyllobacterium", "Devosia",
             "Ochrobactrum", "Nocardioides", "Micromonospora", "Kitasatospora",
             "Streptosporangium", "Nonomuraea", "Pseudoxanthomonas", "Luteibacter", "Dyella",
             "Rahnella", "Ewingella", "Leclercia", "Cellulomonas", "Curtobacterium", "Kocuria",
             "Brevibacillus", "Bosea", "Methylobacterium", "Azorhizobium", "Allorhizobium",
             "Neorhizobium", "Pararhizobium"]


def pgp_regex(names):
    return re.compile(r"^(" + "|".join(sorted(set(names))) + r")(_[A-Z]+)?$")


def pk(a, b):
    return (a, b) if a < b else (b, a)


def load_presence(md, genera, qc):
    """study -> (cultura, matriz bool amostras x generos do conjunto)"""
    out = {}
    for s, cult in zip(qc.study, qc.cultura):
        f = os.path.join(md, "genus", f"{s}.tsv.gz")
        if not os.path.exists(f):
            continue
        with gzip.open(f, "rt") as fh:
            head = fh.readline().rstrip("\n").split("\t")[1:]
        keep = [g for g in head if g in genera]
        if not keep:
            out[s] = (cult, pd.DataFrame())
            continue
        d = pd.read_csv(f, sep="\t", index_col=0, usecols=["sample"] + keep)
        out[s] = (cult, d > 0)
    return out


_CACHE = {}


def copresence(members, pres):
    """Fracao de amostras com todos os membros; estudos/culturas com copresenca >= 50%. Vetorizado."""
    key = id(pres)
    if key not in _CACHE:
        _CACHE.clear()
        _CACHE[key] = [(cult, {g: j for j, g in enumerate(d.columns)}, d.values.astype(bool), len(d))
                       for cult, d in pres.values()]
    tot = hit = studies = 0
    crops = set()
    for cult, idx, X, n in _CACHE[key]:
        if n == 0:
            continue
        tot += n
        try:
            cols = [idx[m] for m in members]
        except KeyError:
            continue
        x = X[:, cols].all(axis=1)
        hit += int(x.sum())
        if x.mean() >= 0.5:
            studies += 1; crops.add(cult)
    return hit / max(tot, 1), studies, sorted(c for c in crops if isinstance(c, str))


def scenario(name, tiers, core, neg, pgp, pres, max_size, out_dir):
    pos = core[(core.sinal == "positivo") & core.tier.isin(tiers)]
    G = nx.Graph()
    info = {}
    for _, r in pos.iterrows():
        G.add_edge(r.genus_a, r.genus_b)
        info[pk(r.genus_a, r.genus_b)] = r
    rows = []
    for c in nx.enumerate_all_cliques(G):
        if len(c) < 2:
            continue
        if len(c) > max_size:
            break
        c = sorted(c)
        if any(pk(a, b) in neg for a, b in itertools.combinations(c, 2)):
            continue
        npgp = sum(bool(pgp.match(g)) for g in c)
        if npgp == 0:
            continue
        e = [info[pk(a, b)] for a, b in itertools.combinations(c, 2)]
        cp, nst, crops = copresence(c, pres)
        rows.append({"tamanho": len(c), "membros": " + ".join(c),
                     "pgp": " + ".join(g for g in c if pgp.match(g)), "n_pgp": npgp,
                     "r_min": min(x.r_mu for x in e), "r_medio": float(np.mean([x.r_mu for x in e])),
                     "pi_lo_min": min(x.r_pi_lo for x in e), "I2_max": max(x.I2 for x in e),
                     "k_min": min(x.k for x in e), "tier_max": max(x.tier for x in e),
                     "copresenca": cp, "estudos_copresente": nst, "n_culturas": len(crops),
                     "culturas": ";".join(crops)})
    d = pd.DataFrame(rows)
    if len(d):
        d = d.sort_values(["tamanho", "r_min", "copresenca"], ascending=[True, False, False])
    d.to_csv(os.path.join(out_dir, f"consorcios_{name}.tsv"), sep="\t", index=False, float_format="%.4g")

    lines = [f"== Cenario {name} (tiers {','.join(map(str, tiers))}) ==",
             f"Rede positiva: {G.number_of_nodes()} generos, {G.number_of_edges()} pares"]
    if len(d) == 0:
        lines.append("Nenhum consorcio com genero PGP.")
        return lines
    cnt = d.groupby("tamanho").agg(n=("membros", "size"), r_min_mediano=("r_min", "median"),
                                   copres_mediana=("copresenca", "median"))
    lines.append("Consorcios por tamanho (todos os pares positivos, nenhum negativo, >=1 PGP):")
    lines.append(cnt.to_string(float_format=lambda v: f"{v:.3f}"))
    for t, sub in d.groupby("tamanho"):
        lines.append(f"\nTop 5, tamanho {t}:")
        for _, r in sub.head(5).iterrows():
            lines.append(f"  {r.membros}  | r_min={r.r_min:.2f} r_med={r.r_medio:.2f} "
                         f"copres={r.copresenca:.2f} estudos={r.estudos_copresente} culturas={r.n_culturas}")
    return lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--max-size", type=int, default=6)
    ap.add_argument("--permissivo", action="store_true",
                    help="veta so pares consistente_neg (ignora associados negativos)")
    ap.add_argument("--pgp-classica", action="store_true", help="usa so a lista classica de PGP")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, "consorcios"); os.makedirs(out, exist_ok=True)

    names = PGP_CLASSIC + ([] if a.pgp_classica else PGP_EXTRA)
    pgp = pgp_regex(names)
    pd.DataFrame({"genus": sorted(set(names)),
                  "lista": ["classica" if n in PGP_CLASSIC else "ampliada" for n in sorted(set(names))]}
                 ).to_csv(os.path.join(out, "lista_pgp.tsv"), sep="\t", index=False)

    core = pd.read_csv(os.path.join(mm, "nucleo_robusto.tsv"), sep="\t")
    meta = pd.read_csv(os.path.join(mm, "meta_pairs.tsv.gz"), sep="\t",
                       usecols=["genus_a", "genus_b", "mu", "fdr", "classe"])
    isneg = meta.classe.eq("consistente_neg")
    if not a.permissivo:
        isneg |= (meta.mu < 0) & (meta.fdr < 0.05)
    neg = {pk(x, y) for x, y in zip(meta.genus_a[isneg], meta.genus_b[isneg])}

    qc = pd.read_csv(os.path.join(md, "study_qc.tsv"), sep="\t")
    qc = qc[qc.status == "incluido"]
    genera = set(core.genus_a) | set(core.genus_b)
    pres = load_presence(md, genera, qc)

    lines = [f"Pares vetados por sinal negativo: {len(neg)} "
             f"({'so consistentes' if a.permissivo else 'consistentes + associados FDR<0,05'})",
             f"Lista PGP: {len(set(names))} generos ({'classica' if a.pgp_classica else 'classica + ampliada'})",
             f"Estudos para copresenca: {len(pres)}", ""]
    for name, tiers in (("tier1", [1]), ("tier12", [1, 2])):
        lines += scenario(name, tiers, core, neg, pgp, pres, a.max_size, out) + [""]
    txt = "\n".join(lines)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
