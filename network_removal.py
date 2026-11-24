#!/usr/bin/env python3
# [S4.9] In silico removal
# Supplementary 4 of the manuscript. See README.md for the full map.
"""
Remocao in silico: a rede de consenso colapsa ou muda quando os consorcios sao retirados?

Rede: pares consistentes da meta-analise principal (meta_pairs.tsv.gz), sem sinal por padrao
(--sinal positivo usa so as associacoes positivas).
Para cada consorcio da shortlist, para a uniao de todos os membros e para os centros PGP:
  remove os nos e mede, na rede restante:
    frac_lcc   fracao dos nos restantes no maior componente
    n_comp     numero de componentes
    eficiencia eficiencia global (media de 1/distancia)
    Q          modularidade (Louvain)
    arestas_perdidas  fracao das arestas da rede ligadas aos nos removidos
  Nulo pareado por grau: --reps conjuntos aleatorios com a mesma distribuicao de grau
  (exclui os membros). Controle positivo: remover os mesmos n nos de maior grau.
  Para cada metrica: z = (obs - media nula) / dp nula e p empirico (dano >= observado).
Centralidade dos membros (rede intacta): grau, intermediacao, coeficiente de participacao
entre modulos, percentis, abundancia media; "keystone" = intermediacao >= percentil 90
e abundancia abaixo da mediana dos generos da rede.

Saidas: <root>/merged/meta/remocao/ (remocao.tsv, centralidade_membros.tsv, resumo.txt)
Uso: LD_PRELOAD=... python3 network_removal.py --root ~/syncom [--reps 100] [--sinal positivo]
"""
import argparse, os, sys
import numpy as np, pandas as pd, networkx as nx
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path, connected_components
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meta_pairs import entrada


def metrics(G, nodes_keep, seed):
    H = G.subgraph(nodes_keep)
    n = H.number_of_nodes()
    if n < 2 or H.number_of_edges() == 0:
        return dict(frac_lcc=0, n_comp=n, eficiencia=0, Q=0)
    idx = {v: i for i, v in enumerate(H.nodes)}
    e = np.array([(idx[a], idx[b]) for a, b in H.edges])
    A = csr_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(n, n))
    nc, lab = connected_components(A, directed=False)
    D = shortest_path(A, directed=False, unweighted=True)
    with np.errstate(divide="ignore"):
        inv = 1 / D
    inv[~np.isfinite(inv)] = 0
    eff = inv.sum() / (n * (n - 1))
    comm = nx.community.louvain_communities(H, seed=seed)
    Q = nx.community.modularity(H, comm)
    return dict(frac_lcc=np.bincount(lab).max() / n, n_comp=nc, eficiencia=eff, Q=Q)


def degree_matched(G, removed, pool_by_deg, degs_sorted, rng):
    out = set()
    for v in removed:
        dv = G.degree(v)
        # grau igual; se nao houver, o mais proximo
        cand = [u for u in pool_by_deg.get(dv, []) if u not in out]
        j = 1
        while not cand:
            cand = [u for d in (dv - j, dv + j) for u in pool_by_deg.get(d, []) if u not in out]
            j += 1
            if j > max(degs_sorted):
                break
        out.add(cand[rng.integers(len(cand))])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--sinal", choices=["todos", "positivo"], default="todos")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, "remocao" + ("" if a.sinal == "todos" else "_positivo")); os.makedirs(out, exist_ok=True)

    m = pd.read_csv(os.path.join(mm, "meta_pairs.tsv.gz"), sep="\t", usecols=["genus_a", "genus_b", "classe", "mu", "fdr"])
    m = m[entrada(m) & ((m.mu > 0) | (a.sinal == "todos"))]
    G = nx.Graph(); G.add_edges_from(zip(m.genus_a, m.genus_b))
    N, E = G.number_of_nodes(), G.number_of_edges()

    sl = pd.read_csv(os.path.join(mm, "consorcios", "shortlist.tsv"), sep="\t")
    sets = [(r.id, [r.centro] + r.parceiros.split(" + ")) for _, r in sl.iterrows()]
    allm = sorted({g for _, s in sets for g in s})
    sets.append(("TODOS_membros", allm))
    sets.append(("TODOS_centros", sorted(set(sl.centro))))
    sets = [(i, [g for g in s if g in G]) for i, s in sets]

    base = metrics(G, list(G.nodes), a.seed)
    pool = [v for v in G.nodes if v not in set(allm)]
    pool_by_deg = {}
    for v in pool:
        pool_by_deg.setdefault(G.degree(v), []).append(v)
    degs_sorted = sorted(pool_by_deg)
    by_deg = sorted(G.nodes, key=lambda v: -G.degree(v))
    keys = ["frac_lcc", "n_comp", "eficiencia", "Q"]
    # dano = queda de frac_lcc e eficiencia; aumento de n_comp; mudanca de Q (aumento = fragmentacao modular)
    dmg = {"frac_lcc": -1, "eficiencia": -1, "n_comp": 1, "Q": 1}

    rows = []
    for sid, s in sets:
        if not s:
            continue
        rest = [v for v in G.nodes if v not in set(s)]
        obs = metrics(G, rest, a.seed)
        lost = sum(1 for u, v in G.edges if u in s or v in s) / E
        null = []
        for _ in range(a.reps):
            rr = degree_matched(G, s, pool_by_deg, degs_sorted, rng)
            null.append(metrics(G, [v for v in G.nodes if v not in rr], a.seed))
        nl = pd.DataFrame(null)
        top = metrics(G, [v for v in G.nodes if v not in set(by_deg[:len(s)])], a.seed)
        row = {"conjunto": sid, "n_nos": len(s), "membros": " + ".join(s), "arestas_perdidas": lost,
               "grau_medio": float(np.mean([G.degree(v) for v in s]))}
        for k in keys:
            mu, sd = nl[k].mean(), nl[k].std(ddof=1)
            row[f"{k}_obs"] = obs[k]; row[f"{k}_delta"] = obs[k] - base[k]
            row[f"{k}_nulo"] = mu; row[f"{k}_z"] = (obs[k] - mu) / sd if sd > 0 else np.nan
            row[f"{k}_p"] = ((dmg[k] * (nl[k] - obs[k]) >= 0).sum() + 1) / (a.reps + 1)
            row[f"{k}_topgrau"] = top[k]
        rows.append(row)
        print(f"{sid}: {len(s)} nos | dLCC={row['frac_lcc_delta']:+.3f} (z={row['frac_lcc_z']:.1f}) "
              f"dEf={row['eficiencia_delta']:+.4f} (z={row['eficiencia_z']:.1f}) "
              f"dQ={row['Q_delta']:+.3f} (z={row['Q_z']:.1f})", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(out, "remocao.tsv"), sep="\t", index=False, float_format="%.4g")

    # centralidade
    bt = nx.betweenness_centrality(G)
    comm = nx.community.louvain_communities(G, seed=a.seed)
    mod = {v: i for i, c in enumerate(comm) for v in c}
    def part(v):
        k = G.degree(v)
        cnt = pd.Series([mod[u] for u in G.neighbors(v)]).value_counts()
        return 1 - ((cnt / k) ** 2).sum()
    occ = pd.read_csv(os.path.join(md, "genus_occupancy.tsv"), sep="\t").set_index("genus")
    ab = occ.abund_media.reindex(list(G.nodes))
    deg = pd.Series(dict(G.degree)); btw = pd.Series(bt)
    C = pd.DataFrame({"genus": allm})
    C = C[C.genus.isin(G.nodes)]
    C["grau"] = C.genus.map(deg); C["grau_pct"] = C.genus.map(deg.rank(pct=True))
    C["intermediacao"] = C.genus.map(btw); C["interm_pct"] = C.genus.map(btw.rank(pct=True))
    C["participacao"] = C.genus.map(part); C["modulo"] = C.genus.map(mod)
    C["abund_media"] = C.genus.map(ab)
    C["keystone"] = (C.interm_pct >= 0.9) & (C.abund_media < ab.median())
    C = C.sort_values("intermediacao", ascending=False)
    C.to_csv(os.path.join(out, "centralidade_membros.tsv"), sep="\t", index=False, float_format="%.4g")

    L = [f"Rede ({a.sinal}): {N} generos, {E} arestas | intacta: frac_lcc={base['frac_lcc']:.3f} "
         f"n_comp={base['n_comp']} eficiencia={base['eficiencia']:.4f} Q={base['Q']:.3f}",
         f"Nulo: {a.reps} conjuntos pareados por grau por consorcio", "",
         "Remocao (delta = depois - antes; z e p contra o nulo pareado por grau; topgrau = remover os n nos de maior grau):"]
    show = R[["conjunto", "n_nos", "arestas_perdidas", "frac_lcc_delta", "frac_lcc_z", "frac_lcc_p",
              "eficiencia_delta", "eficiencia_z", "eficiencia_p", "Q_delta", "Q_z", "Q_p",
              "frac_lcc_topgrau"]]
    L.append(show.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    L += ["", "Centralidade dos membros (rede intacta):",
          C.to_string(index=False, float_format=lambda v: f"{v:.3g}")]
    txt = "\n".join(L)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print("\n" + txt)


if __name__ == "__main__":
    main()
