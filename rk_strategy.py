#!/usr/bin/env python3
# [S3.14] r/K strategy hypothesis
# Supplementary 3 of the manuscript. See README.md for the full map.
"""
Teste da hipotese r/k: a variabilidade da abundancia explica por que alguns filos formam
menos associacoes consistentes?

Hipotese: copiotroficos de resposta rapida (r) oscilam em pulsos, de modo que a correlacao com
qualquer parceiro muda de estudo para estudo e a meta-analise le isso como inconsistencia.
Oligotroficos de crescimento lento (k) tem abundancia estavel e acoplada ao ambiente.

Metricas por genero (das tabelas de genero, mesmo CLR de effect_sizes.py):
  sd_intra   desvio-padrao do CLR entre amostras, mediana dos estudos onde o genero e testavel
             (prevalencia >= --min-prev): quanto o genero oscila DENTRO de um experimento
  sd_inter   desvio-padrao, entre estudos, do log10 da abundancia media: quanto varia ENTRE
             contextos
  sd_intra_resid   residuo de sd_intra ~ log abundancia (lei de Taylor): oscilacao acima do
                   esperado para o nivel de abundancia
Analises (desfecho = parceiros com associacao positiva; offset = log pares testados):
  1. base + Bacillota                     (efeito do filo, sem controle de variabilidade)
  2. base + sd_intra_resid + sd_inter     (a variabilidade prediz a taxa de associacao?)
  3. 2 + Bacillota                        -> atenuacao = mediacao pela variabilidade
  4. Por filo: variabilidade mediana e razao observado/esperado, com correlacao entre as duas
  5. Coerencia com o eixo r/k: correlacao da variabilidade com as copias do 16S
     (de vies_filos/por_genero.tsv, se existir)

Saidas: <root>/merged/meta/vies_filos/ (rk_por_genero.tsv, rk_por_filo.tsv, rk_resumo.txt)
Uso: LD_PRELOAD=... python3 rk_strategy.py --root ~/syncom [--excluir .../excluir_amostras.tsv]
"""
import argparse, os, sys
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from effect_sizes import load_genus, load_curation, drop_samples
from inoculant_bias import irls
from spore_bias import tax_map

log = lambda *x: print(*x, file=sys.stderr, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--excluir", default=None, help="excluir_amostras.tsv (padrao: o da curadoria)")
    ap.add_argument("--min-prev", type=float, default=0.2)
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, "vies_filos"); os.makedirs(out, exist_ok=True)
    exc = a.excluir or os.path.join(md, "curadoria_final", "excluir_amostras.tsv")
    _, excl = load_curation(None, exc)

    qc = pd.read_csv(os.path.join(md, "study_qc.tsv"), sep="\t", dtype=str)
    inc = sorted(qc.study[qc.status == "incluido"])
    sd_i, ab_s = {}, {}
    for i, s in enumerate(inc, 1):
        samples, genera, M = load_genus(os.path.join(md, "genus", f"{s}.tsv.gz"))
        samples, M = drop_samples(samples, M, excl.get(s))
        if len(samples) < 20:
            continue
        named = np.array([not g.startswith("NA") for g in genera])
        prev = (M > 0).mean(0)
        test = named & (prev >= a.min_prev)
        if test.sum() < 5:
            continue
        rel = M / M.sum(1, keepdims=True)
        Y = np.column_stack([M[:, test], M[:, ~test].sum(1)]) + 0.5
        L = np.log(Y); clr = (L - L.mean(1, keepdims=True))[:, :-1]
        sd = clr.std(0, ddof=1)
        for g, v, am in zip([g for g, t in zip(genera, test) if t], sd, rel[:, test].mean(0)):
            sd_i.setdefault(g, []).append(v)
            ab_s.setdefault(g, []).append(am)
        if i % 25 == 0 or i == len(inc):
            log(f"[{i}/{len(inc)}] estudos")

    V = pd.DataFrame({"sd_intra": {g: float(np.median(v)) for g, v in sd_i.items()},
                      "n_estudos_testavel": {g: len(v) for g, v in sd_i.items()},
                      "sd_inter": {g: float(np.std(np.log10(np.clip(v, 1e-8, None)), ddof=1))
                                   if len(v) > 2 else np.nan for g, v in ab_s.items()}})

    m = pd.read_csv(os.path.join(mm, "meta_pairs.tsv.gz"), sep="\t",
                    usecols=["genus_a", "genus_b", "mu", "fdr"])
    pos = m[(m.mu > 0) & (m.fdr < 0.05)]
    cnt = lambda d: pd.concat([d.genus_a, d.genus_b]).value_counts()
    G = pd.DataFrame({"n_testados": cnt(m), "grau_pos": cnt(pos)}).fillna(0)
    G = G[G.n_testados >= 10]
    occ = pd.read_csv(os.path.join(md, "genus_occupancy.tsv"), sep="\t").set_index("genus")
    G = G.join(occ[["ocupacao", "abund_media"]]).join(tax_map(md)).join(V)
    G = G.dropna(subset=["ocupacao", "abund_media", "sd_intra", "sd_inter"])
    G["log_abund"] = np.log10(G.abund_media.clip(lower=1e-8))
    G["taxa_pos"] = G.grau_pos / G.n_testados
    G["bacillota"] = G.filo.fillna("").str.match(r"^(Bacillota|Firmicutes)").astype(int)
    # lei de Taylor: oscilacao acima do esperado para a abundancia
    X = np.column_stack([np.ones(len(G)), G.log_abund.values])
    b, *_ = np.linalg.lstsq(X, G.sd_intra.values, rcond=None)
    G["sd_intra_resid"] = G.sd_intra - X @ b
    y = G.grau_pos.values.astype(float); off = G.n_testados.values

    L = [f"Generos com variabilidade estimada: {len(G)}",
         f"Lei de Taylor: sd_intra = {b[0]:.2f} {b[1]:+.2f} x log10(abundancia) "
         f"(r = {stats.pearsonr(G.log_abund, G.sd_intra).statistic:.2f})", ""]
    base = ["log_abund", "ocupacao"]
    t1, p1 = irls(y, G[base + ["bacillota"]].values, off, True, base + ["bacillota"])
    c2 = base + ["sd_intra_resid", "sd_inter"]
    t2, p2 = irls(y, G[c2].values, off, True, c2)
    t3, p3 = irls(y, G[c2 + ["bacillota"]].values, off, True, c2 + ["bacillota"])
    L += [f"1. Base + Bacillota (quasi = {p1:.1f}):", t1.to_string(float_format=lambda v: f"{v:.3g}"), "",
          f"2. Base + variabilidade (quasi = {p2:.1f}):", t2.to_string(float_format=lambda v: f"{v:.3g}"), "",
          f"3. 2 + Bacillota (quasi = {p3:.1f}):", t3.to_string(float_format=lambda v: f"{v:.3g}"), ""]
    c1, c3 = t1.loc["bacillota", "coef"], t3.loc["bacillota", "coef"]
    L += [f"Mediacao pela variabilidade: coeficiente de Bacillota {c1:+.3f} -> {c3:+.3f} "
          f"(atenuacao de {100 * (1 - c3 / c1) if c1 else float('nan'):.0f}%)",
          "   perto de 100% = a oscilacao da abundancia explica o efeito do filo", ""]

    bb, *_ = np.linalg.lstsq(np.column_stack([np.ones(len(G)), G[base].values]), np.log(np.maximum(y, .5) / off),
                             rcond=None)
    G["esperado"] = np.exp(bb[0] + G[base].values @ bb[1:]) * off
    P = G.dropna(subset=["filo"]).groupby("filo").apply(lambda d: pd.Series(
        {"n_generos": len(d), "razao_obs_esp": d.grau_pos.sum() / max(d.esperado.sum(), 1e-9),
         "sd_intra_mediano": d.sd_intra.median(), "sd_intra_resid_mediano": d.sd_intra_resid.median(),
         "sd_inter_mediano": d.sd_inter.median()}), include_groups=False)
    P = P[P.n_generos >= 8].sort_values("razao_obs_esp")
    P.to_csv(os.path.join(out, "rk_por_filo.tsv"), sep="\t", float_format="%.4g")
    L += ["4. Por filo (ordenado por razao observado/esperado):",
          P.to_string(float_format=lambda v: f"{v:.3g}")]
    if len(P) >= 5:
        r = stats.spearmanr(P.sd_intra_resid_mediano, P.razao_obs_esp)
        L.append(f"   correlacao entre filos: oscilacao (residuo) x razao obs/esp: rho = {r.statistic:.2f}, "
                 f"p = {r.pvalue:.3g} (n = {len(P)} filos)")

    fv = os.path.join(out, "por_genero.tsv")
    if os.path.exists(fv):
        ssu = pd.read_csv(fv, sep="\t", index_col=0)
        if "ssu_count" in ssu:
            J = G.join(ssu[["ssu_count"]]).dropna(subset=["ssu_count"])
            if len(J) > 30:
                r1 = stats.spearmanr(J.ssu_count, J.sd_intra_resid)
                L += ["", f"5. Coerencia com o eixo r/k: copias do 16S x oscilacao (residuo): "
                          f"rho = {r1.statistic:.2f}, p = {r1.pvalue:.3g} (n = {len(J)})"]
    G.to_csv(os.path.join(out, "rk_por_genero.tsv"), sep="\t", float_format="%.5g")
    txt = "\n".join(L)
    open(os.path.join(out, "rk_resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
