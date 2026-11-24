#!/usr/bin/env python3
# [S3.9] Temporal validation
# Supplementary 3 of the manuscript. See README.md for the full map.
"""
Validacao temporal: os pares descobertos em estudos antigos se replicam em estudos novos?

1. Data de publicacao de cada estudo (ENA, campo first_public), salva em merged/study_dates.tsv
2. Descoberta: estudos publicados ate --corte (padrao 2024-12-31); meta-analise (k >= --min-k-desc)
3. Validacao: estudos publicados depois do corte; meta-analise (k >= --min-k-val) dos pares
   consistentes na descoberta
4. Replicacao = mesmo sinal e p unilateral < 0,05 na direcao prevista; comparada com
   pares "associados" (FDR < 0,05, intervalo de predicao cruza zero) e pares ns da descoberta
5. Pares dos consorcios da shortlist: status na descoberta e na validacao

Saidas: <root>/merged/meta/validacao_temporal/ (descoberta.tsv.gz, validacao.tsv, resumo.txt)
Uso: LD_PRELOAD=... python3 temporal_validation.py --root ~/syncom [--corte 2024-12-31]
"""
import argparse, itertools, json, os, sys, time, urllib.parse, urllib.request
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meta_pairs import load, meta

API = "https://www.ebi.ac.uk/ena/portal/api/search"


def ena_date(acc):
    for result in ("read_study", "study"):
        q = urllib.parse.urlencode({"result": result, "query": f'study_accession="{acc}"',
                                    "fields": "first_public", "format": "json"})
        for _ in range(3):
            try:
                with urllib.request.urlopen(f"{API}?{q}", timeout=60) as r:
                    js = json.loads(r.read() or b"[]")
                if js and js[0].get("first_public"):
                    return js[0]["first_public"]
                break
            except Exception:
                time.sleep(3)
    # alternativa: data mais antiga entre as corridas
    q = urllib.parse.urlencode({"result": "read_run", "query": f'study_accession="{acc}"',
                                "fields": "first_public", "format": "json", "limit": 0})
    try:
        with urllib.request.urlopen(f"{API}?{q}", timeout=120) as r:
            js = json.loads(r.read() or b"[]")
        ds = sorted(x["first_public"] for x in js if x.get("first_public"))
        return ds[0] if ds else ""
    except Exception:
        return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--corte", default="2024-12-31")
    ap.add_argument("--effects", default="effects")
    ap.add_argument("--min-k-desc", type=int, default=10)
    ap.add_argument("--min-k-val", type=int, default=5)
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, "validacao_temporal"); os.makedirs(out, exist_ok=True)

    d, qc = load(root, False, a.effects)
    fdates = os.path.join(md, "study_dates.tsv")
    dates = pd.read_csv(fdates, sep="\t", dtype=str) if os.path.exists(fdates) else \
        pd.DataFrame(columns=["study", "first_public"])
    falta = sorted(set(qc.study) - set(dates.study))
    for i, s in enumerate(falta, 1):
        dates.loc[len(dates)] = [s, ena_date(s)]
        if i % 20 == 0 or i == len(falta):
            print(f"datas ENA: {i}/{len(falta)}", flush=True)
            dates.to_csv(fdates, sep="\t", index=False)
    dt = dict(zip(dates.study, dates.first_public.fillna("")))
    sem = [s for s in qc.study if not dt.get(s)]
    desc = [s for s in qc.study if dt.get(s) and dt[s] <= a.corte]
    val = [s for s in qc.study if dt.get(s) and dt[s] > a.corte]
    L = [f"Corte: {a.corte} | descoberta: {len(desc)} estudos | validacao: {len(val)} estudos | sem data: {len(sem)}"]

    D = meta(d[d.study.isin(desc)], a.min_k_desc)
    D["pair"] = D.genus_a + "|" + D.genus_b
    D.to_csv(os.path.join(out, "descoberta.tsv.gz"), sep="\t", index=False, float_format="%.5g")
    vc = D.classe.value_counts()
    L.append(f"Descoberta: {len(D)} pares testados | consistente_pos {vc.get('consistente_pos', 0)} | "
             f"consistente_neg {vc.get('consistente_neg', 0)}")

    groups = {"consistente": D.classe.str.startswith("consistente"),
              "associado": D.classe.str.startswith("associado"),
              "ns": D.classe.eq("ns")}
    test = D[groups["consistente"] | groups["associado"]].copy()
    ns = D[groups["ns"]].sample(min(20000, int(groups["ns"].sum())), random_state=1)
    test = pd.concat([test, ns])
    fsl = os.path.join(mm, "consorcios", "shortlist.tsv")
    if os.path.exists(fsl):                               # garante os pares dos consorcios
        S0 = pd.read_csv(fsl, sep="\t")
        want = {"|".join(p) for _, r in S0.iterrows()
                for p in itertools.combinations(sorted([r.centro] + r.parceiros.split(" + ")), 2)}
        test = pd.concat([test, D[D.pair.isin(want - set(test.pair))]])
    dv = d[d.study.isin(val) & d.pair.isin(set(test.pair))]
    V = meta(dv, a.min_k_val)
    V.index = V.genus_a + "|" + V.genus_b
    j = test.set_index("pair").join(V[["k", "mu", "se_mu", "r_mu"]].add_suffix("_val"), how="inner")
    sgn = np.sign(j.mu)
    zdir = sgn * j.mu_val / j.se_mu_val
    j["p_unilateral"] = stats.norm.sf(zdir)
    j["mesmo_sinal"] = np.sign(j.mu_val) == sgn
    j["replicado"] = j.mesmo_sinal & (j.p_unilateral < 0.05)
    j["grupo"] = np.where(j.classe.str.startswith("consistente"), "consistente",
                          np.where(j.classe.str.startswith("associado"), "associado", "ns"))
    j.to_csv(os.path.join(out, "validacao.tsv"), sep="\t", float_format="%.4g")

    rho = stats.spearmanr(j.mu, j.mu_val).correlation
    L += [f"Pares testaveis na validacao (k >= {a.min_k_val}): {len(j)} | rho descoberta x validacao = {rho:.3f}", "",
          "Replicacao por classe na descoberta (mesmo sinal; p unilateral < 0,05):"]
    tab = j.groupby(["grupo", np.where(j.mu > 0, "positivo", "negativo")]).agg(
        n=("replicado", "size"), mesmo_sinal=("mesmo_sinal", "mean"), replicado=("replicado", "mean"),
        r_desc=("r_mu", "median"), r_val=("r_mu_val", "median"))
    L.append(tab.to_string(float_format=lambda v: f"{v:.3f}"))
    sig = j[j.grupo != "ns"]
    L.append(f"Rede de entrada (FDR < 0,05 na descoberta): {len(sig)} pares testaveis | mesmo sinal "
             f"{sig.mesmo_sinal.mean():.3f} | replicados {sig.replicado.mean():.3f}")

    sl = os.path.join(mm, "consorcios", "shortlist.tsv")
    if os.path.exists(sl):
        S = pd.read_csv(sl, sep="\t")
        L += ["", "Pares dos consorcios (classe na descoberta -> validacao):"]
        Di = D.set_index("pair")
        for _, r in S.iterrows():
            mem = [r.centro] + r.parceiros.split(" + ")
            parts = []
            for x, y in itertools.combinations(sorted(mem), 2):
                k = f"{x}|{y}"
                if k not in Di.index:
                    continue
                c = Di.loc[k, "classe"]
                if k in j.index:
                    rr = j.loc[k]
                    st = f"r_val={rr.r_mu_val:.2f} k={int(rr.k_val)} {'REPLICADO' if rr.replicado else 'nao replicado'}"
                else:
                    st = "nao testavel na validacao"
                parts.append(f"{x}-{y} [{c}, r={Di.loc[k, 'r_mu']:.2f}] {st}")
            if parts:
                L.append(f"  {r.id}: " + " | ".join(parts))
    txt = "\n".join(L)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
