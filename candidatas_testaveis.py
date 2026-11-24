#!/usr/bin/env python3
# [S4.5] Feasibility funnel
# Supplementary 4 of the manuscript. See README.md for the full map.
"""
Combinacoes de 2-3 generos testaveis in planta: todos os membros com linhagem
acessivel, de origem edafica e com menos de --max-clinico de isolados clinicos;
todos os pares internos de nivel 1 com coabundancia; todos replicados no tempo.

Requer a avaliacao de linhagens estendida aos generos dessas combinacoes:
  python3 isolates.py --root <root> --generos "$(cat generos_tier1.txt)"
guardada em <root>/genomas/isolados/por_genero_tier1.tsv (ver --classes).

Saida: <root>/merged/meta/consorcios/candidatas_testaveis.tsv e resumo.txt
"""
import argparse, itertools, os, pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--classes", default=None,
                    help="tabela de cultivabilidade (padrao: genomas/isolados/por_genero_tier1.tsv)")
    ap.add_argument("--max-clinico", type=float, default=0.10)
    ap.add_argument("--tamanhos", default="2,3")
    a = ap.parse_args()
    MM = os.path.join(a.root, "merged", "meta")
    cd = os.path.join(MM, "consorcios")
    cls_f = a.classes or os.path.join(a.root, "genomas", "isolados", "por_genero_tier1.tsv")

    C = pd.read_csv(os.path.join(cd, "consorcios_tier1.tsv"), sep="\t")
    # split literal: pandas trata padrao de varios caracteres como regex
    C["ms"] = C.membros.map(lambda x: str(x).split(" + "))
    lo, hi = (int(x) for x in a.tamanhos.split(","))
    S = C[C.tamanho.between(lo, hi)].copy()
    genes = sorted({x for m in S.ms for x in m})
    open(os.path.join(cd, "generos_tier1.txt"), "w").write(",".join(genes))

    if not os.path.exists(cls_f):
        raise SystemExit(f"falta {cls_f}\nrode: python3 isolates.py --root {a.root} "
                         f'--generos "$(cat {cd}/generos_tier1.txt)" e guarde a saida ali')
    G = pd.read_csv(cls_f, sep="\t", index_col=0)
    G["fc"] = G.isolados_clinicos / G.genomas_isolados.clip(lower=1)
    acess = set(G.index[(G.classe == "acessivel") & (G.fc < a.max_clinico)
                        & (G.isolados_de_solo > 0)])
    D = pd.read_csv(os.path.join(MM, "diagnostico_zeros", "diag_pares.tsv"), sep="\t")
    V = pd.read_csv(os.path.join(MM, "validacao_temporal", "validacao.tsv"), sep="\t")
    L = pd.read_csv(os.path.join(cd, "lista_pgp.tsv"), sep="\t")
    pgp = set(L.genus)
    rep = {frozenset((r.genus_a, r.genus_b)) for _, r in V.iterrows() if bool(r.replicado)}
    coab = {frozenset((r.genus_a, r.genus_b)) for _, r in D.iterrows()
            if r.tier == 1 and r.sinal == "positivo" and r.diag == "coabundancia"}
    tudo = lambda m, s: all(frozenset(p) in s for p in itertools.combinations(sorted(m), 2))

    L1 = [f"generos nas combinacoes de {lo}-{hi} membros: {len(genes)} "
          f"| com cepas acessiveis: {len(acess & set(genes))}",
          f"{'combinacoes enumeradas':34s} {len(S):5d}"]
    for nome, cond in [("todos com linhagem acessivel", lambda m: all(x in acess for x in m)),
                       ("todos os pares coabundancia", lambda m: tudo(m, coab)),
                       ("todos replicados no tempo", lambda m: tudo(m, rep))]:
        S = S[S.ms.map(cond)].copy()
        L1.append(f"{nome:34s} {len(S):5d}")
    S["n_pgp"] = S.ms.map(lambda m: sum(x in pgp for x in m))
    S = S.sort_values("copresenca", ascending=False)
    S.drop(columns=["ms"]).to_csv(os.path.join(cd, "candidatas_testaveis.tsv"),
                                  sep="\t", index=False)
    cols = ["tamanho", "membros", "r_min", "copresenca", "estudos_copresente", "n_culturas", "n_pgp"]
    L1 += ["", "maiores copresencas:", S[cols].head(12).round(3).to_string(index=False)]
    N = S[S.n_pgp == 0]
    L1 += ["", f"sem nenhum genero com promocao relatada: {len(N)}",
           N[cols[:-1]].head(10).round(3).to_string(index=False)]
    # pares entre genero com inoculante registrado e parceiro
    frames = []
    for sub, rot in (("pgp", "promotores"), ("rizobios", "rizobios")):
        f = os.path.join(MM, sub, "parceiros.tsv")
        if os.path.exists(f):
            frames.append(pd.read_csv(f, sep="\t").assign(grupo=rot))
    if frames:
        Pp = pd.concat(frames)
        sim = lambda c: c.astype(str).str.lower().isin(["true", "1", "sim"])
        etapas = [("positivos (FDR < 0,05)", lambda d: d),
                  ("confirmado no hospedeiro", lambda d: d[sim(d.ok_leg)]),
                  ("parceiro acessivel", lambda d: d[[x in acess for x in d.parceiro]]),
                  ("coabundancia", lambda d: d[d.diag == "coabundancia"]),
                  ("replicado no tempo", lambda d: d[sim(d.replicado)]),
                  ("copresenca >= 0,25", lambda d: d[d.copresenca_leg >= 0.25])]
        q = Pp
        for nome, fn in etapas:
            q = fn(q)
            L1.append(f"   {nome:28s} {len(q):4d}  (rizobio: {int((q.grupo=='rizobios').sum())})")
        q = q.sort_values("copresenca_leg", ascending=False)
        q.to_csv(os.path.join(cd, "parceiros_elegiveis.tsv"), sep="\t", index=False)
        L1 += ["", f"pares centro-parceiro elegiveis: {len(q)} "
                   f"(de {len(Pp)} positivos) | com rizobio: {int((q.grupo=='rizobios').sum())}",
               q[["rizobio", "parceiro", "r_leg", "copresenca_leg", "grupo"]]
                .head(10).round(3).to_string(index=False)]

    txt = "\n".join(L1)
    open(os.path.join(cd, "candidatas_testaveis_resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
