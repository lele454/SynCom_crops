#!/usr/bin/env python3
# [S3.4] Two-stage pooling by country and sugarcane exclusion
# Supplementary 3 of the manuscript. See README.md for the full map.
"""
Robustez final: (1) dependencia entre estudos do mesmo pais, (2) exclusao da cana.

1. PAIS (tres niveis). A meta-analise principal trata os estudos como independentes, mas
   estudos do mesmo pais compartilham solo, clima e praticas. Aqui, um ajuste em dois estagios
   que aproxima o modelo de tres niveis: DerSimonian-Laird DENTRO de cada pais (quando ha mais
   de um estudo) e depois DL ENTRE paises. Com isso o pais, e nao o estudo, passa a ser a
   unidade de replicacao. Compara mu, largura do IC e perda de significancia.
   Tambem reporta a concentracao por pais (fracao de estudos do pais dominante) por par.

2. CANA. O recorte diz "culturas anuais/ciclo curto", mas a cana e semiperene (socas por varios
   anos). Refaz a meta-analise sem os estudos de cana e mede a retencao dos pares da rede.

Saidas: <root>/merged/meta/robustez_extra/ (pais_tres_niveis.tsv.gz, sem_cana.tsv.gz, resumo.txt)
Uso: LD_PRELOAD=... python3 robustness_extra.py --root ~/syncom [--pares nucleo|rede]
"""
import argparse, os, sys
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meta_pairs import load, meta, entrada


def dl(z, v):
    """DerSimonian-Laird: mu, var(mu), tau2."""
    z = np.asarray(z, float); v = np.asarray(v, float)
    k = len(z)
    if k == 1:
        return z[0], v[0], 0.0
    w = 1 / v
    mu_f = (w * z).sum() / w.sum()
    Q = (w * (z - mu_f) ** 2).sum()
    C = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    ws = 1 / (v + tau2)
    return (ws * z).sum() / ws.sum(), 1 / ws.sum(), tau2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--pares", choices=["nucleo", "rede"], default="nucleo",
                    help="nucleo: tier 1 do nucleo robusto; rede: todos os pares da rede de entrada")
    ap.add_argument("--min-paises", type=int, default=3)
    ap.add_argument("--min-k", type=int, default=10, help="estudos minimos por par (sem cana)")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, "robustez_extra"); os.makedirs(out, exist_ok=True)

    M = pd.read_csv(os.path.join(mm, "meta_pairs.tsv.gz"), sep="\t")
    M.index = M.genus_a + "|" + M.genus_b
    if a.pares == "nucleo" and os.path.exists(os.path.join(mm, "nucleo_robusto.tsv")):
        C = pd.read_csv(os.path.join(mm, "nucleo_robusto.tsv"), sep="\t")
        alvo = set((C[C.tier == 1].genus_a + "|" + C[C.tier == 1].genus_b))
    else:
        alvo = set(M.index[entrada(M)])
    d, qc = load(root, False)
    pais = dict(zip(qc.study, qc.country.fillna("NA")))
    cult = dict(zip(qc.study, qc.cultura.fillna("")))
    d = d[d.pair.isin(alvo)].copy()
    d["pais"] = d.study.map(pais)
    d["v"] = d.se ** 2
    L = [f"Pares analisados ({a.pares}): {len(alvo)} | estudos: {d.study.nunique()} | "
         f"paises: {d.pais.nunique()}", ""]

    # ---------------------------------------------------------------- 1. pais
    rows = []
    for p, g in d.groupby("pair"):
        por_pais = []
        for pa, gg in g.groupby("pais"):
            mu_p, v_p, _ = dl(gg.z.values, gg.v.values)
            por_pais.append((pa, mu_p, v_p, len(gg)))
        if len(por_pais) < a.min_paises:
            continue
        mu3, v3, tau3 = dl([x[1] for x in por_pais], [x[2] for x in por_pais])
        mu2, v2, tau2 = dl(g.z.values, g.v.values)
        n_dom = max(x[3] for x in por_pais) / len(g)
        rows.append({"pair": p, "k_estudos": len(g), "k_paises": len(por_pais),
                     "frac_pais_dominante": n_dom,
                     "mu_2niveis": mu2, "se_2niveis": np.sqrt(v2), "tau2_estudo": tau2,
                     "mu_3niveis": mu3, "se_3niveis": np.sqrt(v3), "tau2_pais": tau3})
    P = pd.DataFrame(rows)
    if len(P):
        for t in ("2niveis", "3niveis"):
            P[f"p_{t}"] = 2 * stats.norm.sf(np.abs(P[f"mu_{t}"] / P[f"se_{t}"]))
            o = np.argsort(P[f"p_{t}"].values); m = len(P)
            q = np.empty(m)
            q[o] = np.minimum.accumulate((P[f"p_{t}"].values[o] * m /
                                          np.arange(1, m + 1))[::-1])[::-1]
            P[f"fdr_{t}"] = np.minimum(q, 1)
            P[f"r_{t}"] = np.tanh(P[f"mu_{t}"])
        P["mantem"] = (P.fdr_3niveis < 0.05) & (np.sign(P.mu_3niveis) == np.sign(P.mu_2niveis))
        P.to_csv(os.path.join(out, "pais_tres_niveis.tsv.gz"), sep="\t", index=False,
                 float_format="%.5g")
        sig2 = P.fdr_2niveis < 0.05
        L += ["1. Pais como unidade de replicacao (dois estagios: DL dentro do pais, DL entre paises):",
              f"   pares com >= {a.min_paises} paises: {len(P)} | significativos a 2 niveis: "
              f"{int(sig2.sum())} | ainda significativos a 3 niveis: {int((sig2 & P.mantem).sum())} "
              f"({(P.mantem[sig2].mean() if sig2.any() else float('nan')):.0%})",
              f"   correlacao dos efeitos: rho = "
              f"{stats.spearmanr(P.mu_2niveis, P.mu_3niveis).correlation:.3f}",
              f"   |r| mediano: {P.r_2niveis.abs().median():.3f} (2 niveis) vs "
              f"{P.r_3niveis.abs().median():.3f} (3 niveis)",
              f"   IC mediano {P.se_2niveis.median() * 3.92:.3f} (2 niveis) vs "
              f"{P.se_3niveis.median() * 3.92:.3f} (3 niveis) | "
              f"alargamento mediano: {(P.se_3niveis / P.se_2niveis).median():.2f}x",
              f"   concentracao: fracao mediana de estudos do pais dominante = "
              f"{P.frac_pais_dominante.median():.2f}",
              f"   tau2 mediano: estudo {P.tau2_estudo.median():.4f} | entre paises "
              f"{P.tau2_pais.median():.4f}", ""]
    else:
        L += [f"1. Pais: nenhum par com >= {a.min_paises} paises", ""]

    # ---------------------------------------------------------------- 2. sem cana
    d2, _ = load(root, False)
    cana = {s for s, c in cult.items() if c == "cana"}
    d2 = d2[~d2.study.isin(cana)]
    S = meta(d2, a.min_k)
    S.index = S.genus_a + "|" + S.genus_b
    S.to_csv(os.path.join(out, "sem_cana.tsv.gz"), sep="\t", index=False, float_format="%.5g")
    if S.empty:
        L.append(f"2. Sem cana: nenhum par com k >= {a.min_k}")
        txt = "\n".join(L); open(os.path.join(out, "resumo.txt"), "w").write(txt); print(txt); return
    base = M[entrada(M)]
    both = base.index.intersection(S.index)
    ok = (S.loc[both, "fdr"] < 0.05) & (np.sign(S.loc[both, "mu"]) == np.sign(base.loc[both, "mu"]))
    L += [f"2. Sem cana ({len(cana)} estudos excluidos; {d2.study.nunique()} restantes):",
          f"   pares testaveis: {len(S)} | da rede de entrada: {len(both)}/{len(base)} | "
          f"mantidos com mesmo sinal e FDR < 0,05: {int(ok.sum())} ({ok.mean():.0%})",
          f"   correlacao dos efeitos: rho = "
          f"{stats.spearmanr(M.loc[both, 'mu'], S.loc[both, 'mu']).correlation:.3f}"]
    txt = "\n".join(L)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
