#!/usr/bin/env python3
# [S4.8] Core shared by the tied sets
# Supplementary 4 of the manuscript. See README.md for the full map.
"""Nucleo do modulo: empates de clique maximo, frequencia por genero, replicacao."""
import os, pandas as pd, networkx as nx
from collections import Counter
R = os.path.expanduser("~/syncom"); MM = f"{R}/merged/meta"
D = pd.read_csv(f"{MM}/diagnostico_zeros/diag_pares.tsv", sep="\t")
e = D[(D.tier == 1) & (D.sinal == "positivo") & (D.diag == "coabundancia")]
G = nx.Graph([(r.genus_a, r.genus_b) for _, r in e.iterrows()])
cls = list(nx.find_cliques(G)); mx = max(len(c) for c in cls)
top = [c for c in cls if len(c) == mx]
print(f"clique maximo: {mx} generos | empates: {len(top)}")
print(f"nucleo comum a todos os empates: {' + '.join(sorted(set.intersection(*map(set, top))))}")
print("\nfrequencia entre os empates:")
for g, n in Counter(x for c in top for x in c).most_common():
    print(f"  {g}: {n}/{len(top)}")

V = pd.read_csv(f"{MM}/validacao_temporal/validacao.tsv", sep="\t")
rep = [c for c in V.columns if "replic" in c.lower()]
k = rep[0] if rep else None
print(f"\ncoluna de replicacao: {k} | colunas: {list(V.columns)}")
mp = {}
if k:
    for _, r in V.iterrows():
        mp[frozenset((r.genus_a, r.genus_b))] = bool(r[k])
rr = D.set_index(D.apply(lambda r: frozenset((r.genus_a, r.genus_b)), axis=1))
def stats(c):
    m = sorted(c); tot = ok = 0
    for i in range(len(m)):
        for j in range(i + 1, len(m)):
            v = mp.get(frozenset((m[i], m[j])))
            if v is not None:
                tot += 1; ok += v
    return ok, tot
print("\nmaior subconjunto com TODOS os pares replicados no tempo:")
best = None
for c in top:
    for n in range(len(c), 2, -1):
        import itertools
        for s in itertools.combinations(sorted(c), n):
            ok, tot = stats(s)
            if tot == n * (n - 1) // 2 and ok == tot:
                if best is None or n > len(best[0]):
                    best = (s, ok)
                break
        if best and len(best[0]) == n:
            break
if best:
    print(f"  {len(best[0])} membros: {' + '.join(best[0])} ({best[1]} pares, todos replicados)")
else:
    print("  nenhum subconjunto de 3+ membros com todos os pares replicados")
