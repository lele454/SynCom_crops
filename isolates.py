#!/usr/bin/env python3
# [S5.4] Strain availability and isolation source per genus
# Supplementary 5 of the manuscript. See README.md for the full map.
"""
Disponibilidade de isolados para os generos dos consorcios candidatos.

Do metadata do GTDB (bac120_metadata, ja baixado por phylum_bias.py), por genero:
  genomas_isolados   genomas que NAO sao MAG/SAG (ncbi_genome_category): houve cultivo
  frac_isolados      fracao dos genomas do genero que vem de isolado (0 = so MAG)
  linhagens_tipo     genomas marcados como material-tipo: especie descrita
  isolados_de_solo / isolados_clinicos  fonte de isolamento (so para genomas de isolado)
Classes: acessivel | dificil | sem_cultivado | conferir_fatia (nome ausente nesta versao do
GTDB, mas o nome-base tem isolados: conferir em que fatia caem as linhagens).

Alem do ranking de consorcios, inclui os pares centro-parceiro de cada subpasta passada em
--parceiros-dirs (rizobios -> RZ01.., pgp -> PG01.., fosfato -> PO01..).

Saidas: <root>/genomas/isolados/ (por_genero.tsv, por_consorcio.tsv, resumo.txt)
Uso: python3 isolates.py --root ~/syncom --parceiros-dirs rizobios,pgp
"""
import argparse, glob, os, re, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shortlist import base

COLS = ["gtdb_taxonomy", "ncbi_genome_category", "gtdb_type_designation_ncbi_taxa",
        "gtdb_type_designation", "ncbi_isolation_source", "ncbi_strain_identifiers",
        "checkm_completeness"]
SOLO = re.compile(r"soil|rhizosph|root|plant|sediment|compost|manure|rice|maize|wheat|soybean|"
                  r"sorghum|legume|field|agricultur", re.I)
CLIN = re.compile(r"blood|clinic|patient|human|hospital|urine|sputum|wound|catheter|faec|feces", re.I)
PREFIXO = {"rizobios": "RZ", "pgp": "PG", "fosfato": "PO"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--metadata", default=None)
    ap.add_argument("--min-isolados", type=int, default=3)
    ap.add_argument("--parceiros-dirs", default="rizobios",
                    help="subpastas de merged/meta com parceiros.tsv, separadas por virgula")
    ap.add_argument("--generos", default="", help="generos extras a tabular (sem virar consorcio)")
    ap.add_argument("--sem-rizobios", action="store_true", help="exclui a subpasta rizobios")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    gd = os.path.join(root, "genomas"); out = os.path.join(gd, "isolados"); os.makedirs(out, exist_ok=True)
    f = a.metadata or next(iter(sorted(glob.glob(os.path.join(gd, "bac120_metadata*")))), None)
    if not f:
        raise SystemExit("bac120_metadata nao encontrado em ~/syncom/genomas (rode phylum_bias.py "
                         "ou passe --metadata)")
    mm = os.path.join(root, "merged", "meta", "consorcios")
    fr = os.path.join(mm, "ranking_aptos.tsv")
    R = pd.read_csv(fr if os.path.exists(fr) else os.path.join(mm, "ranking_final.tsv"),
                    sep="\t", index_col=0)
    rankcol = "rank_apto" if "rank_apto" in R else "rank"
    membros = {i: [x.strip() for x in str(m).split("+")] for i, m in R.membros.items()}

    dirs = [x.strip() for x in a.parceiros_dirs.split(",") if x.strip()]
    if a.sem_rizobios:
        dirs = [x for x in dirs if x != "rizobios"]
    novos = {}
    for nd, sub in enumerate(dirs):
        fp = os.path.join(root, "merged", "meta", sub, "parceiros.tsv")
        if not os.path.exists(fp):
            print(f"aviso: {fp} nao existe, ignorado"); continue
        P = pd.read_csv(fp, sep="\t")
        P = P[P.ok_leg] if "ok_leg" in P else P
        pre = PREFIXO.get(sub, sub[:2].upper())
        for i, r in enumerate(P.sort_values("r_leg", ascending=False).itertuples(), 1):
            k = f"{pre}{i:02d}"
            membros[k] = [r.rizobio, r.parceiro]
            novos[k] = 1000 * (nd + 1) + i
    if novos:
        R = pd.concat([R, pd.DataFrame({"membros": [" + ".join(membros[k]) for k in novos],
                                        rankcol: list(novos.values())}, index=list(novos))])
    extra = {x.strip() for x in a.generos.split(",") if x.strip()}
    alvo = {g for v in membros.values() for g in v} | extra

    bases = {base(g) for g in alvo}
    keep = []
    use = None
    for ch in pd.read_csv(f, sep="\t", low_memory=False, chunksize=200000):
        if use is None:
            use = [c for c in COLS if c in ch.columns]
        d = ch[use].copy()
        d["genus"] = d.gtdb_taxonomy.str.extract(r"g__([^;]+)")[0]
        keep.append(d[d.genus.map(base).isin(bases)])   # inclui as outras fatias do mesmo nome-base
    D = pd.concat(keep, ignore_index=True)
    D["species"] = D.gtdb_taxonomy.str.extract(r"s__([^;]+)")[0]
    cat = D.ncbi_genome_category.fillna("").str.lower() if "ncbi_genome_category" in D else pd.Series("", index=D.index)
    D["isolado"] = ~cat.str.contains("metagenome|mag|single cell|sag", na=False)
    tipo = pd.Series("", index=D.index)
    for c in ("gtdb_type_designation", "gtdb_type_designation_ncbi_taxa"):
        if c in D:
            tipo = tipo.str.cat(D[c].fillna("").astype(str), sep=";")
    D["tipo"] = tipo.str.contains("type strain|type species|assembly from type", case=False, na=False)
    src = D.ncbi_isolation_source.fillna("") if "ncbi_isolation_source" in D else pd.Series("", index=D.index)
    D["fonte_solo"] = src.map(lambda v: bool(SOLO.search(str(v)))) & D.isolado
    D["fonte_clinica"] = src.map(lambda v: bool(CLIN.search(str(v)))) & D.isolado

    DA = D[D.genus.isin(alvo)]
    G = DA.groupby("genus").agg(genomas=("genus", "size"), genomas_isolados=("isolado", "sum"),
                                linhagens_tipo=("tipo", "sum"), especies=("species", "nunique"),
                                isolados_de_solo=("fonte_solo", "sum"),
                                isolados_clinicos=("fonte_clinica", "sum"))
    G["frac_isolados"] = G.genomas_isolados / G.genomas
    # linhagem-tipo pressupoe descricao formal e, portanto, cultivo: a ausencia de genoma de
    # isolado nesse caso significa "cultivado sem genoma depositado", nao "nao cultivavel"
    G["classe"] = np.where(G.genomas_isolados == 0,
                  np.where(G.linhagens_tipo > 0, "cultivado_sem_genoma", "sem_cultivado"),
                  np.where((G.genomas_isolados >= a.min_isolados) & (G.linhagens_tipo > 0),
                           "acessivel", "dificil"))
    exemplos = {}
    if "ncbi_isolation_source" in D:
        for g, d in DA[DA.isolado & DA.fonte_solo].groupby("genus"):
            vc = d.ncbi_isolation_source.fillna("").str.lower().value_counts().head(3)
            exemplos[g] = "; ".join(f"{k[:40]} ({v})" for k, v in vc.items() if k)
    G["fontes_solo"] = pd.Series(exemplos)
    G["nome_base"] = [base(g) for g in G.index]
    for g in sorted(alvo - set(G.index)):
        G.loc[g] = {c: (0 if G[c].dtype.kind in "if" else "") for c in G.columns}
        G.loc[g, "classe"] = "sem_cultivado"; G.loc[g, "nome_base"] = base(g)
    porbase = D.assign(nb=D.genus.map(base)).query("isolado").groupby("nb").size()
    G["isolados_no_nome_base"] = [int(porbase.get(base(g), 0)) for g in G.index]
    G.loc[(G.genomas == 0) & (G.isolados_no_nome_base > 0), "classe"] = "conferir_fatia"
    G.to_csv(os.path.join(out, "por_genero.tsv"), sep="\t", float_format="%.3g")

    ordem = {"acessivel": 0, "dificil": 1, "cultivado_sem_genoma": 2,
             "conferir_fatia": 3, "sem_cultivado": 4}
    rows = []
    for i, mem in membros.items():
        cl = [G.classe.get(g, "sem_cultivado") for g in mem]
        pior = max(cl, key=lambda c: ordem[c])
        rows.append({"id": i, "membros": " + ".join(mem), "classe_limitante": pior,
                     "generos_sem_cultivado": "; ".join(g for g, c in zip(mem, cl) if c == "sem_cultivado"),
                     "generos_conferir_fatia": "; ".join(g for g, c in zip(mem, cl) if c == "conferir_fatia"),
                     "generos_cultivado_sem_genoma": "; ".join(g for g, c in zip(mem, cl)
                                                               if c == "cultivado_sem_genoma"),
                     "generos_dificeis": "; ".join(g for g, c in zip(mem, cl) if c == "dificil"),
                     "min_isolados_clinicos": int(max(G.isolados_clinicos.get(g, 0) for g in mem)),
                     "min_isolados_de_solo": int(min(G.isolados_de_solo.get(g, 0) for g in mem)),
                     "rank": R.loc[i, rankcol]})
    C = pd.DataFrame(rows).sort_values(["classe_limitante", "rank"],
                                       key=lambda s: s.map(ordem) if s.name == "classe_limitante" else s)
    C.to_csv(os.path.join(out, "por_consorcio.tsv"), sep="\t", index=False)

    L = [f"Generos dos consorcios: {len(alvo)} | genomas no GTDB: {len(D)} | "
         f"subpastas de parceiros: {', '.join(dirs) if dirs else 'nenhuma'}",
         "Classe por genero: " + ", ".join(f"{k} {v}" for k, v in G.classe.value_counts().items()), "",
         G.sort_values(["classe", "genomas_isolados"], ascending=[True, False])
          [["genomas", "genomas_isolados", "frac_isolados", "linhagens_tipo", "especies",
            "isolados_de_solo", "isolados_clinicos", "classe"]]
          .to_string(float_format=lambda v: f"{v:.2f}"), "",
         "Generos sem genoma no metadata (provavel diferenca de versao do GTDB; "
         "ver isolados_no_nome_base): " +
         (", ".join(sorted(G.index[(G.genomas == 0)])) or "nenhum"), "",
         "Consorcios testaveis in vivo (classe do membro mais limitante):",
         C[["id", "membros", "classe_limitante", "generos_sem_cultivado", "generos_conferir_fatia",
            "generos_cultivado_sem_genoma", "generos_dificeis", "min_isolados_clinicos",
            "min_isolados_de_solo"]].head(40).to_string(index=False)]
    txt = "\n".join(L)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
