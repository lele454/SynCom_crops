#!/usr/bin/env python3
# [S4.4] Shortlist and per-crop analysis
# Supplementary 4 of the manuscript. See README.md for the full map.
"""
Shortlist de consorcios + analise por cultura (universal vs especifico).

1. Seleciona de estrelas.tsv: top --n por elo mais fraco e top --n por cobertura (culturas),
   nos tamanhos 2 e 3 (sem duplicatas)
2. Para cada par centro-parceiro, meta-analise DL por cultura (k >= --min-k-cultura estudos)
   e teste de moderacao entre culturas (Q entre subgrupos, efeitos aleatorios)
3. Copresenca por cultura
4. Classe do consorcio:
     universal : sustentado em >= 4 culturas, nenhuma contraria, sem moderacao (Q_entre p >= 0,05)
     amplo     : sustentado em >= 3 culturas, nenhuma contraria
     restrito  : demais casos
   cultura "sustenta" = todos os pares com k >= minimo e r > 0;
   "contraria" = algum par com r < 0 e IC 95% excluindo zero
5. Anotacoes: cultivabilidade (nome GTDB provisorio = provavelmente sem cultivado) e
   alerta de biosseguranca (generos com patogenos/oportunistas descritos)

Saidas (<root>/merged/meta/consorcios/): shortlist.tsv, shortlist_por_cultura.tsv, shortlist_resumo.txt
Uso: LD_PRELOAD=... python3 shortlist.py --root ~/syncom
"""
import argparse, glob, itertools, os, re, sys
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from consortia import pk, load_presence

PLACEHOLDER = re.compile(r"^(\d|[A-Z]{2,}[A-Z0-9]*\d|UBA\d|[A-Za-z]+-\d|Gp\d|PALSA|Palsa|SG8|RBG|ELB|AG\d|AR\d|"
                         r"CF-|DP-|ER\d|WS-|R-H|SXO|SHN)")
BIOSAFETY = {
    "Agrobacterium": "linhagens fitopatogenicas (plasmideo Ti/Ri)",
    "Burkholderia": "complexo B. cepacia / B. pseudomallei",
    "Pseudomonas": "P. aeruginosa (oportunista); P. syringae (fitopatogeno)",
    "Stenotrophomonas": "S. maltophilia (oportunista)",
    "Chryseobacterium": "C. indologenes (oportunista raro)",
    "Sphingobacterium": "oportunistas raros",
    "Nocardia": "oportunistas (nocardiose)",
    "Klebsiella": "patogenos humanos", "Enterobacter": "oportunistas", "Serratia": "oportunistas",
    "Acinetobacter": "A. baumannii", "Achromobacter": "oportunistas", "Brucella": "patogenos",
    "Ochrobactrum": "oportunistas", "Pantoea": "fitopatogenos/oportunistas",
    "Ralstonia": "fitopatogeno", "Xanthomonas": "fitopatogeno", "Clavibacter": "fitopatogeno",
    "Curtobacterium": "fitopatogeno (C. flaccumfaciens)", "Streptomyces": "S. scabiei (sarna)",
    "Bacillus": "B. cereus/anthracis (grupo cereus)", "Rhodococcus": "R. fascians (fitopatogeno)",
    "Mycobacterium": "micobacterias nao tuberculosas oportunistas (complexo M. tuberculosis no genero)",
    "Mycolicibacterium": "micobacterias nao tuberculosas oportunistas",
    "Flavobacterium": "patogenos de peixe (F. psychrophilum, F. columnare)",
    "Elizabethkingia": "oportunistas", "Cupriavidus": "oportunistas raros",
    "Aquipseudomonas": "grupo P. alcaligenes/mendocina (oportunistas raros)",
}


def base(g):
    return re.sub(r"_[A-Z]+$", "", g)


def dl(z, v):
    """DerSimonian-Laird: mu, se, tau2"""
    w = 1 / v; k = len(z)
    mu_f = (w * z).sum() / w.sum()
    Q = (w * (z - mu_f) ** 2).sum()
    C = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - (k - 1)) / C) if k > 1 and C > 0 else 0.0
    ws = 1 / (v + tau2)
    mu = (ws * z).sum() / ws.sum()
    return mu, np.sqrt(1 / ws.sum()), tau2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--n", type=int, default=10, help="top N por criterio e tamanho")
    ap.add_argument("--tamanhos", default="2,3,4,5", help="tamanhos de consorcio considerados")
    ap.add_argument("--min-k-cultura", type=int, default=3)
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta"); out = os.path.join(mm, "consorcios")

    st = pd.read_csv(os.path.join(out, "estrelas.tsv"), sep="\t")
    sel = []
    for t in [int(x) for x in a.tamanhos.split(",")]:
        s = st[st.tamanho == t]
        sel.append(s.sort_values(["rcond_min", "copresenca"], ascending=False).head(a.n))
        sel.append(s.sort_values(["n_culturas", "estudos_copresente", "copresenca"], ascending=False).head(a.n))
    sl = pd.concat(sel).drop_duplicates("membros").reset_index(drop=True)
    sl["id"] = [f"C{i + 1:02d}" for i in range(len(sl))]

    qc = pd.read_csv(os.path.join(md, "study_qc.tsv"), sep="\t")
    qc = qc[qc.status == "incluido"]
    crop = dict(zip(qc.study, qc.cultura.fillna("NA")))
    crops = sorted(set(crop.values()))

    edges = {}
    for _, r in sl.iterrows():
        edges[r.id] = [pk(r.centro, p) for p in r.parceiros.split(" + ")]
    want = {e for v in edges.values() for e in v}
    parts = []
    for f in sorted(glob.glob(os.path.join(md, "effects", "pairs_*.tsv.gz"))):
        s = os.path.basename(f)[6:-7]
        if s not in crop:
            continue
        d = pd.read_csv(f, sep="\t", usecols=["genus_a", "genus_b", "z", "se"])
        d = d[[pk(x, y) in want for x, y in zip(d.genus_a, d.genus_b)]]
        if len(d):
            d["study"] = s; parts.append(d)
    eff = pd.concat(parts, ignore_index=True)
    eff["pair"] = [pk(x, y) for x, y in zip(eff.genus_a, eff.genus_b)]
    eff["cultura"] = eff.study.map(crop)

    # efeito por par x cultura e moderacao
    pc, mod = {}, {}
    for pr, d in eff.groupby("pair"):
        groups = []
        for c, g in d.groupby("cultura"):
            if len(g) >= a.min_k_cultura:
                mu, se, _ = dl(g.z.values, g.se.values ** 2)
                pc[(pr, c)] = (len(g), mu, se)
                groups.append((mu, se))
        if len(groups) >= 2:
            m = np.array([x for x, _ in groups]); w = 1 / np.array([s for _, s in groups]) ** 2
            Qb = (w * (m - (w * m).sum() / w.sum()) ** 2).sum()
            mod[pr] = stats.chi2.sf(Qb, len(groups) - 1)

    allg = set(g for v in edges.values() for e in v for g in e)
    pres = load_presence(md, allg, qc)

    rows_c, rows = [], []
    for _, r in sl.iterrows():
        mem = [r.centro] + r.parceiros.split(" + ")
        sup, contra = [], []
        for c in crops:
            ks, rs, flag_contra = [], [], False
            for e in edges[r.id]:
                k, mu, se = pc.get((e, c), (0, np.nan, np.nan))
                ks.append(k); rs.append(np.tanh(mu) if k else np.nan)
                if k and mu + 1.96 * se < 0:
                    flag_contra = True
            tot = hit = 0
            for s, (cu, dd) in pres.items():
                if cu != c or len(dd) == 0:
                    continue
                tot += len(dd)
                if all(m in dd.columns for m in mem):
                    hit += int(dd[mem].all(axis=1).sum())
            ok = all(k >= a.min_k_cultura for k in ks) and all(x > 0 for x in rs)
            if ok:
                sup.append(c)
            if flag_contra:
                contra.append(c)
            rows_c.append({"id": r.id, "membros": " + ".join(mem), "cultura": c, "k_min": min(ks),
                           "r_min_cultura": np.nanmin(rs) if any(k for k in ks) else np.nan,
                           "copresenca": hit / tot if tot else np.nan,
                           "status": "sustenta" if ok else ("contraria" if flag_contra else
                                                            ("sem_dados" if min(ks) < a.min_k_cultura else "fraco"))})
        pmod = min((mod.get(e, 1.0) for e in edges[r.id]), default=1.0)
        if len(sup) >= 4 and not contra and pmod >= 0.05:
            classe = "universal"
        elif len(sup) >= 3 and not contra:
            classe = "amplo"
        else:
            classe = "restrito"
        cult = [m for m in mem if PLACEHOLDER.match(m)]
        bio = [f"{m}: {BIOSAFETY[base(m)]}" for m in mem if base(m) in BIOSAFETY]
        rows.append({"id": r.id, "tamanho": r.tamanho, "centro": r.centro, "parceiros": r.parceiros,
                     "rcond_min": r.rcond_min, "copresenca": r.copresenca,
                     "estudos_copresente": r.estudos_copresente, "classe_cultura": classe,
                     "culturas_sustenta": ";".join(sup), "culturas_contraria": ";".join(contra),
                     "p_moderacao_cultura": pmod,
                     "sem_cultivado_provavel": ";".join(cult),
                     "biosseguranca": " | ".join(bio)})
    S = pd.DataFrame(rows); C = pd.DataFrame(rows_c)
    S.to_csv(os.path.join(out, "shortlist.tsv"), sep="\t", index=False, float_format="%.4g")
    C.to_csv(os.path.join(out, "shortlist_por_cultura.tsv"), sep="\t", index=False, float_format="%.4g")

    L = [f"Culturas: {', '.join(crops)}", f"Consorcios na shortlist: {len(S)}", ""]
    for _, r in S.iterrows():
        L.append(f"{r.id} [{r.centro}] + {r.parceiros}  (tam {r.tamanho}, rcond_min={r.rcond_min:.2f})")
        L.append(f"     classe: {r.classe_cultura} | sustenta: {r.culturas_sustenta or '-'} | "
                 f"contraria: {r.culturas_contraria or '-'} | p_moderacao={r.p_moderacao_cultura:.3f}")
        if r.sem_cultivado_provavel:
            L.append(f"     sem cultivado (provavel): {r.sem_cultivado_provavel}")
        if r.biosseguranca:
            L.append(f"     biosseguranca: {r.biosseguranca}")
    L += ["", "Status por cultura (S=sustenta, x=contraria, .=fraco, -=sem dados):"]
    sym = {"sustenta": "S", "contraria": "x", "fraco": ".", "sem_dados": "-"}
    tab = C.pivot(index="id", columns="cultura", values="status").replace(sym)
    L.append(tab.to_string())
    txt = "\n".join(L)
    open(os.path.join(out, "shortlist_resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
