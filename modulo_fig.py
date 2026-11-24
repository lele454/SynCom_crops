#!/usr/bin/env python3
# [S4.7] Largest mutually associated set, co-presence tie-break
# Supplementary 4 of the manuscript. See README.md for the full map.
"""
Dados da Fig6. Tres conjuntos:
  dez       - o clique maximo da rede de nivel 1 positiva com coabundancia (com a
              frequencia de cada genero entre os empates)
  replicado - maior conjunto cujos pares todos replicam no tempo
  acessivel - idem, restrito a generos com linhagem isolada acessivel
Empates de tamanho maximo sao resolvidos pela maior copresenca do conjunto completo,
e o numero de empates e registrado.
Saidas em <root>/merged/meta/modulo/: frequencia.tsv, copresenca.tsv, pares.tsv, resumo.txt
"""
import argparse, itertools, os, pandas as pd, networkx as nx
from collections import Counter
import consortia as C


def escolhe(cliques, pres):
    """Entre os cliques de tamanho maximo, o de maior copresenca. Retorna (conjunto, empates)."""
    if not cliques:
        return [], 0
    mx = max(map(len, cliques))
    tie = [sorted(c) for c in cliques if len(c) == mx]
    if len(tie) == 1:
        return tie[0], 1
    melhor = max(tie, key=lambda c: (C.copresence(c, pres)[0], "+".join(c)))
    return melhor, len(tie)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    md = os.path.join(a.root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, "modulo"); os.makedirs(out, exist_ok=True)
    L = []

    D = pd.read_csv(os.path.join(mm, "diagnostico_zeros", "diag_pares.tsv"), sep="\t")
    V = pd.read_csv(os.path.join(mm, "validacao_temporal", "validacao.tsv"), sep="\t")
    rep = {frozenset((r.genus_a, r.genus_b)): bool(r.replicado) for _, r in V.iterrows()}
    cls = pd.read_csv(os.path.join(a.root, "genomas", "isolados", "por_genero.tsv"),
                      sep="\t", index_col=0)
    ac = set(cls.index[cls.classe == "acessivel"])

    e = D[(D.tier == 1) & (D.sinal == "positivo") & (D.diag == "coabundancia")]
    G = nx.Graph([(r.genus_a, r.genus_b) for _, r in e.iterrows()])
    R = nx.Graph([(r.genus_a, r.genus_b) for _, r in e.iterrows()
                  if rep.get(frozenset((r.genus_a, r.genus_b)), False)])
    Ra = R.subgraph([n for n in R if n in ac])
    L.append(f"rede nivel 1 positiva com coabundancia: {G.number_of_nodes()} generos, "
             f"{G.number_of_edges()} associacoes")
    L.append(f"restrita as replicadas no tempo: {R.number_of_nodes()} generos, "
             f"{R.number_of_edges()} associacoes")
    L.append(f"restrita tambem a linhagem acessivel: {Ra.number_of_nodes()} generos, "
             f"{Ra.number_of_edges()} associacoes")

    cl = list(nx.find_cliques(G)); mx = max(map(len, cl))
    top = [c for c in cl if len(c) == mx]
    freq = Counter(x for c in top for x in c)
    F = pd.DataFrame({"genus": list(freq), "cliques": [freq[g] for g in freq]})
    F["frac"] = (F.cliques / len(top)).round(3)
    F["acessivel"] = F.genus.isin(ac)
    F = F.sort_values(["cliques", "genus"], ascending=[False, True])
    F.to_csv(os.path.join(out, "frequencia.tsv"), sep="\t", index=False)
    L.append(f"clique maximo: {mx} generos | empates: {len(top)} | "
             f"membros distintos entre os empates: {len(freq)} | "
             f"nucleo comum: {len(set.intersection(*map(set, top)))} generos")

    universo = set(F.genus) | set(R) | set(Ra)
    qc = pd.read_csv(os.path.join(md, "study_qc.tsv"), sep="\t")
    qc = qc[qc.status == "incluido"]
    L.append(f"copresenca calculada sobre {len(qc)} estudos")
    pres = C.load_presence(md, universo, qc)

    dez = list(F.genus.head(mx))
    quin, t_rep = escolhe(list(nx.find_cliques(R)), pres)
    acce, t_acc = escolhe(list(nx.find_cliques(Ra)), pres)
    L.append(f"dez mais frequentes: {' + '.join(dez)}")
    L.append(f"replicado: {len(quin)} generos ({t_rep} empates) -> {' + '.join(quin)}")
    L.append(f"acessivel: {len(acce)} generos ({t_acc} empates) -> {' + '.join(acce)}")

    rows, pr = [], []
    for nome, S in (("dez", dez), ("replicado", quin), ("acessivel", acce)):
        for n in range(2, len(S) + 1):
            vals = [C.copresence(list(s), pres) for s in itertools.combinations(S, n)]
            fr = [v[0] for v in vals]
            rows.append({"conjunto": nome, "tamanho": n, "n_subconjuntos": len(vals),
                         "copres_media": sum(fr) / len(fr), "copres_melhor": max(fr),
                         "estudos_melhor": max(v[1] for v in vals)})
        for x, y in itertools.combinations(S, 2):
            m = D[((D.genus_a == x) & (D.genus_b == y)) | ((D.genus_a == y) & (D.genus_b == x))]
            if len(m):
                r = m.iloc[0]
                pr.append({"conjunto": nome, "genus_a": x, "genus_b": y, "r_mu": r.r_mu,
                           "k": r.k, "tier": r.tier, "diag": r.diag,
                           "replicado": rep.get(frozenset((x, y))),
                           "acessivel_a": x in ac, "acessivel_b": y in ac})
    P = pd.DataFrame(rows); P.to_csv(os.path.join(out, "copresenca.tsv"), sep="\t", index=False)
    Q = pd.DataFrame(pr); Q.to_csv(os.path.join(out, "pares.tsv"), sep="\t", index=False)
    txt = "\n".join(L) + "\n\ncopresenca por tamanho:\n" + P.round(3).to_string(index=False)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print(txt)
    print("\npares do conjunto acessivel:")
    print(Q[Q.conjunto == "acessivel"].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
