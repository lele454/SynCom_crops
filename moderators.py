#!/usr/bin/env python3
# [S3.7] Moderation by management
# Supplementary 3 of the manuscript. See README.md for the full map.
"""
Moderadores: as associacoes se mantem SOB manejo (adubacao, seca)? - e nao "depois de descontar".

Na analise principal os fatores de desenho entram como covariaveis, para remover o efeito do
tratamento. Aqui eles sao o objeto: para cada estudo que tem contraste, a correlacao e estimada
separadamente DENTRO de cada nivel e os niveis sao comparados no mesmo estudo (pareado).

Harmonizacao dos niveis (valores dos atributos de desenho, por palavra-chave):
  N_mineral   N, NPK, ureia, fertilizante quimico/mineral
  organica    esterco, composto, FYM, biochar, vermicomposto, organico
  controle    CK, control, sem adubacao, nao adubado
  seca        drought, water stress/deficit, limited watering
  hidrico_ok  well-watered, irrigado, optimal watering
Contrastes: N_mineral vs controle | organica vs controle | seca vs hidrico_ok

Por par e contraste: diferenca de Fisher z entre os niveis, dentro do estudo; meta-analise DL
da diferenca (k >= --min-k). Tambem a fracao de pares que seguem significativos em cada nivel.
Recortes: toda a rede de entrada, pares rizobio-parceiro e pares dos consorcios da shortlist.

Saidas: <root>/merged/meta/moderadores/ (niveis.tsv, por_par.tsv.gz, resumo.txt)
Uso: LD_PRELOAD=... python3 moderators.py --root ~/syncom --runs runs_subset.tsv \
        --attrs ~/syncom/merged/curadoria_final/attrs_curado \
        --design ~/syncom/merged/curadoria_final/design_curado.tsv \
        --excluir ~/syncom/merged/curadoria_final/excluir_amostras.tsv
"""
import argparse, csv, os, re, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from effect_sizes import (load_genus, load_attrs, detect_design, design_matrix, load_curation,
                          drop_samples)
from meta_pairs import meta
from shortlist import base

CATS = {   # palavras-chave: so para pre-sugestao no inventario; o mapa manual tem prioridade.
           # a ordem importa: controle e avaliado antes de N_mineral ("0kgN_ha" e controle)
    "controle": r"^\s*(CK|CON|control|controle|untreated|unfertili[sz]ed|check)\b|no[ _]?fertiliz\w*|"
                r"no[ _]?applicat\w*|without[ _]?fertiliz\w*|zero[ _]?N\b|\bN0\b|sem[ _]?adubacao|"
                r"\b0 ?kg ?N",
    "N_mineral": r"\bNPK\b|\bNP\b|\bNK\b|urea|ureia|nitrogen[ _]?fert|chemical[ _]?fert|"
                 r"mineral[ _]?fert|inorganic[ _]?fert|synthetic[ _]?fert|\d+ ?kg ?N|"
                 r"N[ _]?fertiliz|adubacao[ _]?(quimica|mineral)",
    "organica": r"manure|compost|\bFYM\b|biochar|vermicompost|cowdung|organic[ _]?(fert|amend)|"
                r"esterco|straw[ _]?return|residue[ _]?return",
    "seca": r"drought|water[ _]?(stress|deficit|limit)|\bWS\b|\bDS\b|lim[ _]?watering|"
            r"reduced[ _]?irrigat|deficit[ _]?irrigat|estresse[ _]?hidrico|\bseca\b",
    "hidrico_ok": r"well[ _-]?water|optimal[ _]?watering|opt[ _]?watering|full[ _]?irrigat|irrigated|"
                  r"\bWW\b|normal[ _]?water",
}
GRUPOS_OK = {"tratamento", "curado"}   # local, tempo, genotipo e solo nao sao manejo
CONTRASTES = [("N_mineral", "controle"), ("organica", "controle"), ("seca", "hidrico_ok")]
RIZOBIOS = {"Bradyrhizobium", "Rhizobium", "Mesorhizobium", "Sinorhizobium", "Ensifer",
            "Neorhizobium", "Pararhizobium", "Allorhizobium", "Azorhizobium"}
log = lambda *x: print(*x, file=sys.stderr, flush=True)


def classify(v):
    for c, p in CATS.items():
        if re.search(p, str(v), re.I):
            return c
    return None


def resid(Y, X):
    if X is None or X.shape[1] == 0:
        return Y - Y.mean(0)
    Z = np.column_stack([np.ones(Y.shape[0]), X])
    b, *_ = np.linalg.lstsq(Z, Y, rcond=None)
    return Y - Z @ b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--attrs", required=True)
    ap.add_argument("--design", default=None)
    ap.add_argument("--excluir", default=None)
    ap.add_argument("--inventario", action="store_true",
                    help="so lista os niveis de cada fator (para classificacao manual) e sai")
    ap.add_argument("--mapa", default=None,
                    help="TSV study, chave, valor, categoria: classificacao manual dos niveis")
    ap.add_argument("--min-prev", type=float, default=0.2)
    ap.add_argument("--min-n", type=int, default=10, help="amostras minimas por nivel")
    ap.add_argument("--min-k", type=int, default=5, help="estudos minimos por par e contraste")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, "moderadores"); os.makedirs(out, exist_ok=True)
    forced, excl = load_curation(a.design, a.excluir)

    m = pd.read_csv(os.path.join(mm, "meta_pairs.tsv.gz"), sep="\t",
                    usecols=["genus_a", "genus_b", "mu", "fdr", "r_mu"])
    rede = m[m.fdr < 0.05].copy()
    rede["pair"] = rede.genus_a + "|" + rede.genus_b
    want = set(rede.pair)
    rede = rede.set_index("pair")

    run2bs = {}
    with open(a.runs) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            run2bs[r["run_accession"]] = r["sample_accession"]
    qc = pd.read_csv(os.path.join(md, "study_qc.tsv"), sep="\t", dtype=str)
    inc = sorted(qc.study[qc.status == "incluido"])
    cult = dict(zip(qc.study, qc.cultura.fillna("")))

    mapa = {}
    if a.mapa and os.path.exists(a.mapa):
        for r in csv.DictReader(open(a.mapa), delimiter="\t"):
            if r.get("categoria"):
                mapa[(r["study"], r["chave"], r["valor"])] = r["categoria"].strip()

    if a.inventario:                       # inventario dos niveis, para classificar a mao
        inv = []
        for i, s in enumerate(inc, 1):
            samples, _, M = load_genus(os.path.join(md, "genus", f"{s}.tsv.gz"))
            samples, M = drop_samples(samples, M, excl.get(s))
            if len(samples) < 2 * a.min_n:
                continue
            attrs = load_attrs(os.path.join(a.attrs, f"{s}.tsv"))
            amap = {x: attrs.get(run2bs.get(x, x), {}) for x in samples}
            for k, grp in detect_design(amap, samples, forced.get(s)):
                if grp not in GRUPOS_OK:
                    continue
                vc = pd.Series([amap[x].get(k, "") for x in samples]).value_counts()
                for v, n in vc.items():
                    inv.append((s, cult.get(s, ""), k, grp, str(v)[:120], int(n),
                                classify(v) or ""))
            if i % 25 == 0 or i == len(inc):
                log(f"[{i}/{len(inc)}] estudos")
        I = pd.DataFrame(inv, columns=["study", "cultura", "chave", "grupo", "valor",
                                       "n_amostras", "categoria"])
        f = os.path.join(out, "inventario_niveis.tsv")
        I.to_csv(f, sep="\t", index=False)
        print(f"{len(I)} niveis em {I.study.nunique()} estudos -> {f}")
        print(f"ja classificados por palavra-chave: {int((I.categoria != '').sum())}")
        return

    niveis, rows = [], []
    for i, s in enumerate(inc, 1):
        samples, genera, M = load_genus(os.path.join(md, "genus", f"{s}.tsv.gz"))
        samples, M = drop_samples(samples, M, excl.get(s))
        if len(samples) < 2 * a.min_n:
            continue
        attrs = load_attrs(os.path.join(a.attrs, f"{s}.tsv"))
        amap = {x: attrs.get(run2bs.get(x, x), {}) for x in samples}
        keys = detect_design(amap, samples, forced.get(s))   # lista de (chave, grupo)
        if not keys:
            continue
        cat, usada, melhor = {}, None, -1
        for k, _g in keys:                   # a chave que classifica mais amostras em >= 2 categorias
            if _g not in GRUPOS_OK:
                continue
            c = {x: mapa.get((s, k, amap[x].get(k, "")), classify(amap[x].get(k, ""))) for x in samples}
            nc = sum(1 for v in c.values() if v)
            if len({v for v in c.values() if v}) >= 2 and nc > melhor:
                cat, usada, melhor = c, k, nc
        if not usada:
            continue
        for x in samples:
            if cat.get(x):
                niveis.append((s, usada, amap[x].get(usada, ""), cat[x]))
        named = np.array([not g.startswith("NA") for g in genera])
        test = named & ((M > 0).mean(0) >= a.min_prev)
        if test.sum() < 5:
            continue
        G = [g for g, t in zip(genera, test) if t]
        gi = {g: j for j, g in enumerate(G)}
        Y = np.column_stack([M[:, test], M[:, ~test].sum(1)]) + 0.5
        Lg = np.log(Y); clr = (Lg - Lg.mean(1, keepdims=True))[:, :-1]
        outras = [t for t in keys if t[0] != usada]
        X, _ = design_matrix(amap, samples, outras, max(0, len(samples) // 4)) if outras else (None, None)
        prs = [(p, gi[p.split("|")[0]], gi[p.split("|")[1]]) for p in want
               if p.split("|")[0] in gi and p.split("|")[1] in gi]
        if not prs:
            continue
        for c in {v for v in cat.values() if v}:
            idx = [j for j, x in enumerate(samples) if cat.get(x) == c]
            if len(idx) < a.min_n:
                continue
            Xs = X[idx] if X is not None else None
            R = resid(clr[idx], Xs)
            k = (Xs.shape[1] if Xs is not None else 0)
            dfree = len(idx) - 3 - k
            if dfree < 3:
                continue
            C = np.corrcoef(R, rowvar=False)
            for p, ia, ib in prs:
                r = C[ia, ib]
                if np.isfinite(r):
                    rows.append((p, s, c, np.arctanh(np.clip(r, -.999999, .999999)),
                                 1 / np.sqrt(dfree), len(idx)))
        if i % 20 == 0 or i == len(inc):
            log(f"[{i}/{len(inc)}] estudos")

    pd.DataFrame(niveis, columns=["study", "chave", "valor", "categoria"]).drop_duplicates() \
        .to_csv(os.path.join(out, "niveis.tsv"), sep="\t", index=False)
    D = pd.DataFrame(rows, columns=["pair", "study", "cat", "z", "se", "n"])
    if D.empty:
        print("nenhum estudo com contraste harmonizavel"); return

    L = [f"Estudos com contraste harmonizado: {D.study.nunique()}",
         "Amostras por categoria: " + ", ".join(
             f"{c} {int(v)}" for c, v in D.drop_duplicates(['study', 'cat']).groupby('cat').n.sum().items()), ""]
    res_all = []
    for c1, c2 in CONTRASTES:
        A = D[D.cat == c1].set_index(["pair", "study"])
        B = D[D.cat == c2].set_index(["pair", "study"])
        j = A.join(B, lsuffix="_1", rsuffix="_2", how="inner").reset_index()
        if j.empty:
            L.append(f"{c1} vs {c2}: nenhum estudo com os dois niveis"); continue
        j["z"] = j.z_1 - j.z_2
        j["se"] = np.sqrt(j.se_1 ** 2 + j.se_2 ** 2)
        V = meta(j[["pair", "study", "z", "se"]], a.min_k)
        if V.empty:
            L.append(f"{c1} vs {c2}: {j.study.nunique()} estudos, nenhum par com k >= {a.min_k}"); continue
        V["pair"] = V.genus_a + "|" + V.genus_b
        V["contraste"] = f"{c1}_vs_{c2}"
        V = V.join(rede[["r_mu"]].rename(columns={"r_mu": "r_global"}), on="pair")
        # mantem o sinal em cada nivel?
        for tag, sub in ((c1, "_1"), (c2, "_2")):
            g = j.groupby("pair").apply(lambda d: pd.Series(
                {f"r_{tag}": np.tanh(np.average(d["z" + sub], weights=1 / d["se" + sub] ** 2))}),
                include_groups=False)
            V = V.join(g, on="pair")
        res_all.append(V)
        same = (np.sign(V[f"r_{c1}"]) == np.sign(V.r_global)).mean()
        same2 = (np.sign(V[f"r_{c2}"]) == np.sign(V.r_global)).mean()
        mod = (V.fdr < 0.05).mean()
        L += [f"{c1} vs {c2}: {j.study.nunique()} estudos | {len(V)} pares com k >= {a.min_k}",
              f"   mesmo sinal da analise principal: {same:.0%} em {c1}, {same2:.0%} em {c2}",
              f"   pares com moderacao significativa (FDR < 0,05): {mod:.1%}",
              f"   |r| mediano: {V[f'r_{c1}'].abs().median():.3f} em {c1} vs "
              f"{V[f'r_{c2}'].abs().median():.3f} em {c2}",
              f"   diferenca mediana de z (positiva = associacao mais forte em {c1}): "
              f"{V.mu.median():+.3f}"]
            # recortes
        riz = V[[any(base(x) in RIZOBIOS for x in p.split("|")) for p in V.pair]]
        if len(riz):
            L.append(f"   pares rizobio-parceiro: {len(riz)} | mesmo sinal em {c1}: "
                     f"{(np.sign(riz[f'r_{c1}']) == np.sign(riz.r_global)).mean():.0%} | "
                     f"moderacao significativa: {(riz.fdr < 0.05).mean():.1%}")
        L.append("")
    if res_all:
        R = pd.concat(res_all)
        R.to_csv(os.path.join(out, "por_par.tsv.gz"), sep="\t", index=False, float_format="%.4g")
        L += ["Pares com maior moderacao (top 15 por |diferenca|, FDR < 0,05):"]
        top = R[R.fdr < 0.05].reindex(R[R.fdr < 0.05].mu.abs().sort_values(ascending=False).index).head(15)
        for _, r in top.iterrows():
            L.append(f"   {r.pair} [{r.contraste}] dz={r.mu:+.2f} (k={int(r.k)}) r_global={r.r_global:.2f}")
    txt = "\n".join(L)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
