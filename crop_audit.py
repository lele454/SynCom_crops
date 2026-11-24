#!/usr/bin/env python3
# [S1.3] Audit of the assigned host crop
# Supplementary 1 of the manuscript. See README.md for the full map.
"""Auditoria da cultura atribuida a cada estudo: compara a cultura do recorte com a inferida
do campo de hospedeiro e, na falta dele, do titulo. Saida: Tabelas/crop_audit.tsv"""
import argparse, os, re
import pandas as pd

PAD = [
    ("maize",      r"\bzea\s*mays|\bmaize\b|\bcorn\b|\bmilho\b"),
    ("wheat",      r"\btriticum|\bwheat\b|\btrigo\b"),
    ("rice",       r"\boryza|\brice\b|\bpaddy\b|\barroz\b"),
    ("sorghum",    r"\bsorghum\b|sorgo\b"),
    ("sugarcane",  r"\bsaccharum|sugarcane|cana-de"),
    ("soybean",    r"glycine\s*(max|soja)|\bsoybean\b|\bsoja\b"),
    ("common bean", r"phaseolus|common bean|\bfeijao\b|lima bean"),
    ("pea",        r"\bpisum|\bpea\b|\bpeaola\b|\bervilha\b"),
    ("peanut",     r"arachis|peanut|groundnut|amendoim"),
    ("alfalfa",    r"medicago|alfalfa|lucerne"),
    ("clover",     r"trifolium|clover|berseem"),
    ("lupin",      r"lupinus|lupin"),
    ("faba bean",  r"vicia\s*faba|faba bean|fava"),
    ("mung bean",  r"vigna|mung"),
    ("chickpea",   r"\bcicer|chickpea|\bgrao-de-bico"),
    ("lentil",     r"\blens\s|lentil"),
    ("pigeon pea", r"cajanus|pigeon pea"),
    # fora do escopo declarado
    ("buckwheat",  r"fagopyrum|buckwheat"),
    ("millet",     r"panicum|millet|broomcorn"),
    ("brassica",   r"brassica|pak\s*choi|canola|rapeseed|mustard"),
    ("cucumber",   r"cucumis|cucumber|pepino"),
    ("rye",        r"secale|\brye\b"),
    ("barley",     r"hordeum|barley"),
    ("oat",        r"\bavena|\boat\b"),
    ("timothy",    r"\bphleum|timothy"),
]
GRUPO = {"maize": "grass", "wheat": "grass", "rice": "grass", "sorghum": "grass",
         "sugarcane": "grass", "millet": "grass", "rye": "grass", "barley": "grass",
         "oat": "grass", "timothy": "grass",
         "soybean": "legume", "common bean": "legume", "pea": "legume", "peanut": "legume",
         "alfalfa": "legume", "clover": "legume", "lupin": "legume", "faba bean": "legume",
         "mung bean": "legume", "chickpea": "legume", "lentil": "legume", "pigeon pea": "legume"}
PT = {"milho": "maize", "trigo": "wheat", "arroz": "rice", "sorgo": "sorghum",
      "cana": "sugarcane", "soja": "soybean", "feijao/leguminosas": "other legumes"}
LEG = {"common bean", "pea", "peanut", "alfalfa", "clover", "lupin", "faba bean",
       "mung bean", "chickpea", "lentil", "pigeon pea"}


def acha(txt):
    t = str(txt or "").lower()
    return [n for n, p in PAD if re.search(p, t)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--studies", default="subset_studies.tsv")
    ap.add_argument("--out", default="Tabelas")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root)); md = os.path.join(root, "merged")
    out = os.path.abspath(os.path.expanduser(a.out)); os.makedirs(out, exist_ok=True)
    qc = pd.read_csv(os.path.join(md, "study_qc.tsv"), sep="\t", dtype=str).fillna("")
    sub = pd.read_csv(a.studies, sep="\t", dtype=str).fillna("")
    sub = sub.rename(columns={"study_accession": "study"})
    cols = [c for c in ("study", "host", "study_title") if c in sub.columns]
    T = qc[["study", "cultura", "grupo", "status"]].merge(sub[cols], on="study", how="left").fillna("")

    linhas = []
    for _, r in T.iterrows():
        declarada = PT.get(r.cultura, r.cultura)
        do_host = acha(r.get("host", ""))
        do_tit = acha(r.get("study_title", ""))
        cands = do_host or do_tit
        fonte = "host" if do_host else ("titulo" if do_tit else "")
        if not cands:
            vered = "sem informacao"
        elif declarada == "other legumes":
            vered = "ok" if any(c in LEG for c in cands) else "divergente"
        elif declarada in cands:
            vered = "ok"
        elif len(cands) > 1:
            vered = "ambiguo"
        else:
            vered = "divergente"
        grupos = {GRUPO.get(c, "outro") for c in cands}
        linhas.append({"study": r.study, "status": r.status, "declared_crop": declarada,
                       "declared_group": {"gramínea": "grass", "leguminosa": "legume"}.get(r.grupo, r.grupo),
                       "inferred": "; ".join(cands), "inferred_group": "; ".join(sorted(grupos)),
                       "source": fonte, "verdict": vered,
                       "host": r.get("host", ""), "study_title": str(r.get("study_title", ""))[:120]})
    D = pd.DataFrame(linhas)
    ordem = {"divergente": 0, "ambiguo": 1, "sem informacao": 2, "ok": 3}
    D = D.sort_values(["verdict", "study"], key=lambda s: s.map(ordem) if s.name == "verdict" else s)
    D.to_csv(os.path.join(out, "crop_audit.tsv"), sep="\t", index=False)
    print(D.verdict.value_counts().to_string(), "\n")
    rev = D[D.verdict.isin(["divergente", "ambiguo"]) & (D.status == "incluido")]
    print(f"A revisar entre os estudos incluidos: {len(rev)}")
    if len(rev):
        print(rev[["study", "declared_crop", "inferred", "source", "host"]].to_string(index=False))
    print(f"\nTabela completa em {os.path.join(out, 'crop_audit.tsv')}")


if __name__ == "__main__":
    main()
