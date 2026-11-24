#!/usr/bin/env python3
# [S4.12] Composite ranking
# Supplementary 4 of the manuscript. See README.md for the full map.
"""
Tabela final: ranqueamento dos consorcios candidatos integrando todas as linhas de evidencia.

Componentes (cada um convertido em percentil 0-1 entre os consorcios da shortlist):
  associacao      elo mais fraco (r condicional, coabundancia)
  abrangencia     n de culturas que sustentam o consorcio
  copresenca      fracao de amostras com todos os membros
  funcao_PGP      n de tracos PGP na uniao (pgp_traits.py)
  complementar    tracos PGP + tipos de BGC complementares trazidos pelos parceiros
  metabolitos     tipos de BGC na uniao + sideroforo NIS (bgc_traits.py, se disponivel)
Escore = media dos componentes. Robustez: --sorteios vetores de pesos aleatorios (Dirichlet);
  para cada consorcio, frequencia no top 5 e posicao mediana.
Bandeiras (nao entram no escore, mas definem a coluna apto_in_vivo):
  biosseguranca (generos com patogenos descritos), membro provavelmente sem cultivado.
Contexto (nao entra no escore): estruturacao da comunidade (Z de R2 e de alfa).

Saidas: <root>/merged/meta/consorcios/ranking_final.tsv e ranking_resumo.txt
Uso: LD_PRELOAD=... python3 rank_syncoms.py --root ~/syncom
"""
import argparse, itertools, os
import numpy as np, pandas as pd


def pct(x):
    x = pd.Series(x, dtype=float)
    return x.rank(pct=True, method="average").fillna(0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--sorteios", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    mm = os.path.join(root, "merged", "meta"); cd = os.path.join(mm, "consorcios")
    S = pd.read_csv(os.path.join(cd, "shortlist.tsv"), sep="\t").set_index("id")
    P = pd.read_csv(os.path.join(root, "genomas", "pgp", "tracos_por_consorcio.tsv"), sep="\t").set_index("id")
    fb = os.path.join(root, "genomas", "bgc", "bgc_por_consorcio.tsv")
    B = pd.read_csv(fb, sep="\t").set_index("id") if os.path.exists(fb) else None
    fe = os.path.join(mm, "estruturacao", "consorcios.tsv")
    E = pd.read_csv(fe, sep="\t").set_index("id") if os.path.exists(fe) else None

    n = lambda s: s.fillna("").astype(str).map(lambda v: len([x for x in v.split(";") if x]))
    D = pd.DataFrame(index=S.index)
    D["membros"] = S.centro + " + " + S.parceiros
    D["tamanho"] = S.tamanho
    D["rcond_min"] = S.rcond_min
    D["n_culturas"] = n(S.culturas_sustenta)
    D["culturas"] = S.culturas_sustenta
    D["classe_cultura"] = S.classe_cultura
    D["copresenca"] = S.copresenca
    D["tracos_PGP"] = P.n_tracos_uniao.reindex(D.index)
    D["tracos_compl"] = n(P.tracos_complementares.reindex(D.index))
    comp = {"associacao": pct(D.rcond_min), "abrangencia": pct(D.n_culturas),
            "copresenca": pct(D.copresenca), "funcao_PGP": pct(D.tracos_PGP)}
    if B is not None:
        D["tipos_BGC"] = B.n_tipos_uniao.reindex(D.index)
        D["BGC_compl"] = n(B.tipos_complementares.reindex(D.index))
        D["sideroforo"] = B.sideroforo_NIS_em.reindex(D.index).fillna("")
        comp["complementar"] = pct(D.tracos_compl + D.BGC_compl)
        comp["metabolitos"] = pct(D.tipos_BGC + (D.sideroforo != "").astype(int))
    else:
        comp["complementar"] = pct(D.tracos_compl)
    fv = os.path.join(mm, "validacao_temporal", "validacao.tsv")
    if os.path.exists(fv):                     # replicacao temporal dos pares internos
        V = pd.read_csv(fv, sep="\t", index_col=0)
        fr, nt = [], []
        for i in D.index:
            mem = [S.loc[i, "centro"]] + S.loc[i, "parceiros"].split(" + ")
            ps = ["|".join(p) for p in itertools.combinations(sorted(mem), 2)]
            t = [p for p in ps if p in V.index]
            nt.append(f"{sum(bool(V.loc[p, 'replicado']) for p in t)}/{len(t)} de {len(ps)}")
            fr.append(np.mean([bool(V.loc[p, "replicado"]) for p in t]) if t else np.nan)
        D["replicacao_temporal"] = nt
        comp["replicacao"] = pct(pd.Series(fr, index=D.index).fillna(0.5))   # nao testavel = neutro
    C = pd.DataFrame(comp, index=D.index)
    for k in C:
        D[f"c_{k}"] = C[k].round(2)
    D["escore"] = C.mean(axis=1)

    rng = np.random.default_rng(a.seed)
    W = rng.dirichlet(np.ones(C.shape[1]), a.sorteios)
    sc = C.values @ W.T                                   # consorcios x sorteios
    ranks = (-sc).argsort(0).argsort(0) + 1
    D["frac_top5"] = (ranks <= 5).mean(1)
    D["posicao_mediana"] = np.median(ranks, 1)

    D["biosseguranca"] = S.biosseguranca.fillna("")
    D["sem_cultivado"] = S.sem_cultivado_provavel.fillna("")
    D["apto_in_vivo"] = np.where(D.sem_cultivado != "", "nao (membro sem cultivado)",
                                 np.where(D.biosseguranca != "", "sim, com triagem de linhagem", "sim"))
    if E is not None:
        D["estrut_R2_Z"] = E.R2_Z.reindex(D.index); D["estrut_alfa_Z"] = E.alfa_Z.reindex(D.index)
    D = D.sort_values("escore", ascending=False)
    D.insert(0, "rank", range(1, len(D) + 1))
    D.to_csv(os.path.join(cd, "ranking_final.tsv"), sep="\t", float_format="%.3g")

    cols = ["rank", "membros", "escore", "frac_top5", "posicao_mediana", "rcond_min", "n_culturas",
            "copresenca", "tracos_PGP", "tracos_compl"] + (["replicacao_temporal"] if "replicacao_temporal" in D else []) + (["tipos_BGC", "sideroforo"] if B is not None else []) \
           + ["apto_in_vivo"]
    L = [f"Componentes do escore: {', '.join(C.columns)} (pesos iguais; robustez com {a.sorteios} pesos aleatorios)",
         "", D[cols].to_string(float_format=lambda v: f"{v:.2f}"), "",
         "Bandeiras de biosseguranca:"]
    for i, r in D[D.biosseguranca != ""].iterrows():
        L.append(f"  {i} {r.membros}: {r.biosseguranca}")
    txt = "\n".join(L)
    open(os.path.join(cd, "ranking_resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
