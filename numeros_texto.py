#!/usr/bin/env python3
# [S6.8] Values quoted in the Results text
# Supplementary 6 of the manuscript. See README.md for the full map.
"""Numeros para o texto dos Resultados, a partir das saidas da rodada atual."""
import os, numpy as np, pandas as pd, networkx as nx
R = os.path.expanduser("~/syncom"); MM = f"{R}/merged/meta"
def bloco(t): print("\n" + "=" * 8 + " " + t)

D = pd.read_csv(f"{MM}/diagnostico_zeros/diag_pares.tsv", sep="\t")

bloco("niveis e diagnostico")
print(pd.crosstab(D.tier, D.sinal))
for t in (1, 2, 3):
    s = D[(D.tier == t) & (D.sinal == "negativo")]
    print(f"\nnivel {t} negativos (n={len(s)}):")
    print((s.diag.value_counts(normalize=True) * 100).round(1).to_string())

bloco("maior modulo (nivel 1, positivo, coabundancia)")
e = D[(D.tier == 1) & (D.sinal == "positivo") & (D.diag == "coabundancia")]
G = nx.Graph([(r.genus_a, r.genus_b) for _, r in e.iterrows()])
print(f"grafo: {G.number_of_nodes()} generos, {G.number_of_edges()} arestas")
cl = sorted(nx.find_cliques(G), key=len, reverse=True)[:6]
for c in cl:
    print(f"  {len(c)}: {' + '.join(sorted(c))}")

V = None
for nm in ("validacao.tsv",):
    p = f"{MM}/validacao_temporal/{nm}"
    if os.path.exists(p):
        V = pd.read_csv(p, sep="\t"); print(f"\nvalidacao: {list(V.columns)[:8]}")
if V is not None and cl:
    rep = [c for c in V.columns if "replic" in c.lower()]
    if rep:
        k = rep[0]
        V["par"] = V.genus_a + "|" + V.genus_b
        mp = dict(zip(V.par, V[k]))
        for c in cl[:3]:
            m = sorted(c); tot = ok = 0
            for i in range(len(m)):
                for j in range(i + 1, len(m)):
                    v = mp.get(f"{m[i]}|{m[j]}", mp.get(f"{m[j]}|{m[i]}"))
                    if v is not None:
                        tot += 1; ok += bool(v)
            print(f"  {len(m)} membros: {ok}/{tot} pares replicados no tempo")

bloco("complementaridade por tamanho (so combinacoes totalmente anotadas)")
C = pd.read_csv(f"{MM}/consorcios/complementaridade.tsv", sep="\t")
A = C[C.n_com_tracos == C.tamanho]
g = A.groupby("tamanho").agg(n=("membros", "size"), copres=("copresenca", "mean"),
                             compl=("tracos_compl", "mean"), total=("tracos_total", "mean"),
                             copres_med=("copresenca", "median"))
print(g.round(3).to_string())
print("\nganho marginal de traco complementar / perda marginal de copresenca:")
for t in sorted(g.index)[1:]:
    p = t - 1
    if p in g.index and g.copres[p] > 0:
        dg = (g.compl[t] - g.compl[p]) / max(g.compl[p], 1e-9)
        dp = (g.copres[p] - g.copres[t]) / g.copres[p]
        print(f"  {p}->{t}: ganho {dg:+.2f} | perda {dp:.2f} | razao {dg/dp if dp else float('nan'):.2f}")

bloco("moderacao por cultura (shortlist)")
S = pd.read_csv(f"{MM}/consorcios/shortlist.tsv", sep="\t")
print(S.classe_cultura.value_counts().to_string())
if "p_moderacao_cultura" in S:
    p = pd.to_numeric(S.p_moderacao_cultura, errors="coerce").dropna()
    print(f"p_moderacao_cultura: n={len(p)} | p<0,05: {(p < 0.05).sum()} | mediana {p.median():.3f}")
print(f"\nculturas que sustentam (mediana): {S.culturas_sustenta.median()}")
