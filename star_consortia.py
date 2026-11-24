#!/usr/bin/env python3
# [S4.2] Star combinations
# Supplementary 4 of the manuscript. See README.md for the full map.
"""
Consorcios em estrela: um genero PGP central + parceiros auxiliares.

Regras:
  - arestas centro-parceiro: pares positivos do nucleo (tiers --tiers) com diagnostico
    "coabundancia" (diag_pares.tsv); --incluir-cooc aceita tambem so_coocorrencia/coabund_fraca
  - centro: genero PGP (lista classica + ampliada de consortia.py)
  - nenhum par interno (centro-parceiro ou parceiro-parceiro) vetado por sinal negativo
    (consistente_neg ou associado negativo FDR < 0,05)
  - tamanho total 2..--max-size
Metricas: elo mais fraco e medio (r condicional e r principal), pares positivos extras entre
parceiros, copresenca (fracao de amostras com todos os membros), estudos e culturas com
copresenca >= 50% das amostras.
Ranking dentro de cada tamanho: elo mais fraco (r condicional), copresenca, n de culturas.

Saidas: <root>/merged/meta/consorcios/estrelas.tsv, estrelas_saturacao.tsv, estrelas_resumo.txt
Uso: LD_PRELOAD=... python3 star_consortia.py --root ~/syncom
"""
import argparse, itertools, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from consortia import PGP_CLASSIC, PGP_EXTRA, pgp_regex, pk, load_presence, copresence


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--tiers", default="1,2")
    ap.add_argument("--incluir-cooc", action="store_true")
    ap.add_argument("--max-size", type=int, default=6)
    ap.add_argument("--diag", default=None, help="tabela de pares (padrao: diagnostico_zeros/diag_pares.tsv)")
    ap.add_argument("--meta", default=None, help="meta-analise para o veto (padrao: meta_pairs.tsv.gz)")
    ap.add_argument("--cultura", default=None, help="copresenca so nos estudos desta cultura")
    ap.add_argument("--prefixo", default="estrelas", help="prefixo dos arquivos de saida")
    ap.add_argument("--max-parceiros", type=int, default=15,
                    help="parceiros considerados por centro (os de maior r condicional)")
    a = ap.parse_args()
    tiers = [int(x) for x in a.tiers.split(",")]
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, "consorcios"); os.makedirs(out, exist_ok=True)
    pgp = pgp_regex(PGP_CLASSIC + PGP_EXTRA)

    dg = pd.read_csv(a.diag or os.path.join(mm, "diagnostico_zeros", "diag_pares.tsv"), sep="\t")
    ok_diag = {"coabundancia"} | ({"so_coocorrencia", "coabund_fraca"} if a.incluir_cooc else set())
    pos_all = dg[dg.sinal == "positivo"]
    edges = pos_all[pos_all.tier.isin(tiers) & pos_all["diag"].isin(ok_diag) & (pos_all.r_mu_cond > 0)]
    E = {pk(r.genus_a, r.genus_b): r for _, r in edges.iterrows()}
    POS = {pk(x, y) for x, y in zip(pos_all.genus_a, pos_all.genus_b)}     # qualquer positivo consistente

    meta = pd.read_csv(a.meta or os.path.join(mm, "meta_pairs.tsv.gz"), sep="\t",
                       usecols=["genus_a", "genus_b", "mu", "fdr", "classe"])
    isneg = meta.classe.eq("consistente_neg") | ((meta.mu < 0) & (meta.fdr < 0.05))
    NEG = {pk(x, y) for x, y in zip(meta.genus_a[isneg], meta.genus_b[isneg])}

    nb = {}
    for (x, y), r in E.items():
        nb.setdefault(x, []).append((y, r)); nb.setdefault(y, []).append((x, r))
    hubs = sorted(g for g in nb if pgp.match(g))

    qc = pd.read_csv(os.path.join(md, "study_qc.tsv"), sep="\t")
    qc = qc[qc.status == "incluido"]
    if a.cultura:
        qc = qc[qc.cultura == a.cultura]
    genera = set(nb)
    pres = load_presence(md, genera, qc)

    rows = []
    for h in hubs:
        parts = sorted(nb[h], key=lambda t: -t[1].r_mu_cond)[:a.max_parceiros]
        for m in range(1, a.max_size):
            for comb in itertools.combinations(parts, m):
                ps = [p for p, _ in comb]
                mem = [h] + ps
                if any(pk(x, y) in NEG for x, y in itertools.combinations(mem, 2)):
                    continue
                rc = [r.r_mu_cond for _, r in comb]; rm = [r.r_mu for _, r in comb]
                extra = sum(pk(x, y) in POS for x, y in itertools.combinations(ps, 2))
                cp, nst, crops = copresence(mem, pres)
                rows.append({"tamanho": len(mem), "centro": h, "parceiros": " + ".join(sorted(ps)),
                             "pgp_parceiros": " + ".join(p for p in sorted(ps) if pgp.match(p)),
                             "rcond_min": min(rc), "rcond_medio": float(np.mean(rc)),
                             "r_min": min(rm), "tier_max": max(int(r.tier) for _, r in comb),
                             "pos_entre_parceiros": extra, "copresenca": cp,
                             "estudos_copresente": nst, "n_culturas": len(crops),
                             "culturas": ";".join(crops)})
    d = pd.DataFrame(rows)
    if d.empty:
        print("Nenhuma estrela encontrada."); return
    # mesma composicao pode surgir com centros diferentes (dois PGP ligados): mantem a de maior rcond_min
    d["membros"] = [" + ".join(sorted([c] + p.split(" + "))) for c, p in zip(d.centro, d.parceiros)]
    d = (d.sort_values(["tamanho", "rcond_min", "copresenca", "n_culturas"], ascending=[True, False, False, False])
           .drop_duplicates("membros"))
    d.to_csv(os.path.join(out, f"{a.prefixo}.tsv"), sep="\t", index=False, float_format="%.4g")

    sat = d.groupby("tamanho").agg(n=("membros", "size"), centros=("centro", "nunique"),
                                   rcond_min_melhor=("rcond_min", "max"),
                                   copres_melhor=("copresenca", "max"),
                                   copres_mediana=("copresenca", "median"),
                                   culturas_max=("n_culturas", "max"),
                                   estudos_max=("estudos_copresente", "max"))
    sat.to_csv(os.path.join(out, f"{a.prefixo}_saturacao.tsv"), sep="\t", float_format="%.4g")

    L = [f"Arestas centro-parceiro: {len(E)} (tiers {a.tiers}; diag: {', '.join(sorted(ok_diag))})",
         f"Centros PGP com parceiros: {len(hubs)} -> {', '.join(hubs)}", "",
         "Saturacao por tamanho:", sat.to_string(float_format=lambda v: f"{v:.3f}"), ""]
    for t, sub in d.groupby("tamanho"):
        L.append(f"Top 5, tamanho {t} (por elo mais fraco):")
        for _, r in sub.head(5).iterrows():
            L.append(f"  [{r.centro}] + {r.parceiros} | rcond_min={r.rcond_min:.2f} "
                     f"copres={r.copresenca:.3f} estudos={r.estudos_copresente} culturas={r.n_culturas}"
                     f"{' +' + str(r.pos_entre_parceiros) + ' pos' if r.pos_entre_parceiros else ''}")
        best = sub.sort_values(["n_culturas", "copresenca"], ascending=False).iloc[0]
        L.append(f"  maior cobertura: [{best.centro}] + {best.parceiros} | culturas={best.n_culturas} "
                 f"estudos={best.estudos_copresente} copres={best.copresenca:.3f}")
        L.append("")
    txt = "\n".join(L)
    open(os.path.join(out, f"{a.prefixo}_resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
