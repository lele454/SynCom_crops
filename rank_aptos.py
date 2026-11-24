#!/usr/bin/env python3
# [S4.13] Ranking corrected for sequencing effort
# Supplementary 4 of the manuscript. See README.md for the full map.
"""
Ranking final dos consorcios aptos ao teste in vivo, corrigido pelo vies de sequenciamento.

Problema: generos de interesse clinico (Stenotrophomonas, Enterobacter, Acinetobacter, Serratia...)
tem muito mais genomas depositados e melhor anotados, entao acumulam tracos PGP e BGCs. O escore
funcional premia esforco de sequenciamento, nao potencial biologico.

Correcao: por genero, regredir a funcao anotada (tracos PGP esperados; tipos de BGC) sobre
log10(genomas do genero no GTDB, de sp_clusters). O residuo e a funcao "acima do esperado para
o quanto o genero foi sequenciado". O consorcio recebe a media dos residuos dos membros.

Filtro: so consorcios com apto_in_vivo == "sim" (sem bandeira de biosseguranca e sem membro
provavelmente sem cultivado). O ranking completo continua em ranking_final.tsv.

Componentes (percentis entre os consorcios aptos): associacao, abrangencia, copresenca,
funcao_corrigida, complementar, replicacao. Robustez com --sorteios pesos aleatorios (Dirichlet).
Nota de nomenclatura: Brucella no GTDB inclui Ochrobactrum reclassificado.

Saidas: <root>/merged/meta/consorcios/ (ranking_aptos.tsv, vies_sequenciamento.tsv, ranking_aptos.txt)
Uso: LD_PRELOAD=... python3 rank_aptos.py --root ~/syncom
"""
import argparse, glob, gzip, itertools, os, re
import numpy as np, pandas as pd
from scipy import stats

NOMENCLATURA = {"Brucella": "GTDB: inclui Ochrobactrum reclassificado (rizosferico, linhagens PGP)",
                "Aquipseudomonas": "GTDB: grupo P. alcaligenes/mendocina separado de Pseudomonas",
                "Priestia": "antes Bacillus megaterium/aryabhattai",
                "Pristimantibacillus": "genero recente da familia Paenibacillaceae",
                "Chryseomonas": "GTDB: linhagem separada de Pseudomonas"}


def pct(x):
    return pd.Series(x, dtype=float).rank(pct=True, method="average").fillna(0)


def gtdb_effort(root):
    """genero -> (clusters de especie, genomas totais) no GTDB."""
    f = sorted(glob.glob(os.path.join(root, "genomas", "sp_clusters*")))
    if not f:
        return None
    op = gzip.open if f[0].endswith(".gz") else open
    with op(f[0], "rt") as fh:
        sp = pd.read_csv(fh, sep="\t", low_memory=False)
    col = lambda k: next(c for c in sp.columns if k.lower() in c.lower())
    sp["genus"] = sp[col("GTDB taxonomy")].str.extract(r"g__([^;]+)")[0]
    g = sp.groupby("genus").agg(sp_clusters=(col("GTDB taxonomy"), "size"),
                                genomas_gtdb=(col("No. clustered"), "sum"))
    return g


def residual(y, x):
    """residuo de y ~ x (minimos quadrados), com a inclinacao e o R2."""
    ok = np.isfinite(y) & np.isfinite(x)
    if ok.sum() < 5:
        return pd.Series(y - np.nanmean(y)), np.nan, np.nan
    X = np.column_stack([np.ones(ok.sum()), x[ok]])
    b, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
    pred = b[0] + b[1] * x
    r = stats.pearsonr(x[ok], y[ok])
    return pd.Series(y - pred), b[1], r.statistic ** 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--sorteios", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--todos", action="store_true", help="nao filtrar por aptidao (so corrigir o vies)")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    mm = os.path.join(root, "merged", "meta"); cd = os.path.join(mm, "consorcios")
    R = pd.read_csv(os.path.join(cd, "ranking_final.tsv"), sep="\t", index_col=0)
    S = pd.read_csv(os.path.join(cd, "shortlist.tsv"), sep="\t").set_index("id")

    # ---- funcao anotada por genero
    T = pd.read_csv(os.path.join(root, "genomas", "pgp", "tracos_por_genero.tsv"), sep="\t", index_col=0)
    tr = [c for c in T.columns if c != "n_genomas"]
    G = pd.DataFrame({"tracos_esperados": T[tr].sum(axis=1)})
    fb = os.path.join(root, "genomas", "bgc", "bgc_por_genero.tsv")
    if os.path.exists(fb):
        B = pd.read_csv(fb, sep="\t", index_col=0)
        tipos = [c for c in B.columns if c.startswith("BGC_") and c != "BGC_total" and c != "BGC_Unknown"]
        G["bgc_total"] = B.BGC_total
        G["bgc_tipos"] = (B[tipos] >= 0.5).sum(axis=1)
    eff = gtdb_effort(root)
    if eff is None:
        raise SystemExit("sp_clusters do GTDB nao encontrado em ~/syncom/genomas (rode fetch_genomes.py)")
    G = G.join(eff, how="left")
    G["log_genomas"] = np.log10(G.genomas_gtdb.fillna(1).clip(lower=1))

    lines = ["Vies de sequenciamento (funcao anotada ~ log10 genomas do genero no GTDB):"]
    for c in [x for x in ("tracos_esperados", "bgc_tipos", "bgc_total") if x in G]:
        res, slope, r2 = residual(G[c].values.astype(float), G.log_genomas.values)
        G[c + "_resid"] = res.values
        lines.append(f"  {c}: inclinacao = {slope:+.2f} por decada de genomas | R2 = {r2:.2f}")
    G.to_csv(os.path.join(cd, "vies_sequenciamento.tsv"), sep="\t", float_format="%.3g")

    cols_res = [c for c in G.columns if c.endswith("_resid")]
    G["funcao_corrigida"] = G[cols_res].apply(lambda s: (s - s.mean()) / s.std(ddof=0)).mean(axis=1)

    # ---- consorcios
    D = R.copy()
    if not a.todos:
        D = D[D.apto_in_vivo == "sim"]
    if D.empty:
        raise SystemExit("nenhum consorcio apto; use --todos para ver o ranking corrigido completo")
    membros = {i: [S.loc[i, "centro"]] + S.loc[i, "parceiros"].split(" + ") for i in D.index}
    D["funcao_corrigida"] = [np.mean([G.funcao_corrigida.get(g, np.nan) for g in m]) for m in membros.values()]
    D["genomas_gtdb_mediano"] = [np.nanmedian([G.genomas_gtdb.get(g, np.nan) for g in m]) for m in membros.values()]
    D["nota_nomenclatura"] = ["; ".join(f"{g}: {NOMENCLATURA[g]}" for g in m if g in NOMENCLATURA)
                              for m in membros.values()]
    rep = D.replicacao_temporal.fillna("0/0 de 0").astype(str) if "replicacao_temporal" in D else None
    frac_rep = (rep.str.extract(r"^(\d+)/(\d+)").astype(float).pipe(lambda x: x[0] / x[1]).fillna(0.5)
                if rep is not None else pd.Series(0.5, index=D.index))

    comp = {"associacao": pct(D.rcond_min), "abrangencia": pct(D.n_culturas), "copresenca": pct(D.copresenca),
            "funcao_corrigida": pct(D.funcao_corrigida), "complementar": pct(D.tracos_compl.fillna(0)),
            "replicacao": pct(frac_rep)}
    C = pd.DataFrame(comp, index=D.index)
    for k in C:
        D["c_" + k] = C[k].round(2)
    D["escore_apto"] = C.mean(axis=1)
    rng = np.random.default_rng(a.seed)
    W = rng.dirichlet(np.ones(C.shape[1]), a.sorteios)
    sc = C.values @ W.T
    ranks = (-sc).argsort(0).argsort(0) + 1
    D["frac_top5"] = (ranks <= 5).mean(1)
    D["posicao_mediana"] = np.median(ranks, 1)
    D = D.sort_values("escore_apto", ascending=False)
    D.insert(0, "rank_apto", range(1, len(D) + 1))
    D.to_csv(os.path.join(cd, "ranking_aptos.tsv"), sep="\t", float_format="%.3g")

    show = ["rank_apto", "membros", "escore_apto", "frac_top5", "rcond_min", "n_culturas", "copresenca",
            "tracos_PGP", "funcao_corrigida", "genomas_gtdb_mediano"] + \
           (["replicacao_temporal"] if "replicacao_temporal" in D else [])
    lines += ["", f"Consorcios aptos ao teste in vivo: {len(D)} de {len(R)} "
                  f"({'sem filtro' if a.todos else 'sem bandeira de biosseguranca e sem membro sem cultivado'})",
              "funcao_corrigida: desvios-padrao acima (+) ou abaixo (-) do esperado para o esforco de "
              "sequenciamento do genero", "",
              D[show].to_string(float_format=lambda v: f"{v:.2f}")]
    nn = D[D.nota_nomenclatura != ""]
    if len(nn):
        lines += ["", "Notas de nomenclatura GTDB:"]
        lines += [f"  {i}: {r.nota_nomenclatura}" for i, r in nn.iterrows()]
    txt = "\n".join(lines)
    open(os.path.join(cd, "ranking_aptos.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
