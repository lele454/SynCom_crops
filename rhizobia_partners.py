# [S4.6] Partners of genera with a registered inoculant
# Supplementary 4 of the manuscript. See README.md for the full map.
"""
Parceiros auxiliares de um grupo de generos "centro": quem coocorre com o inoculante.

Centros (--centros): rizobios | pgp | fosfato | lista livre de nomes-base separada por virgula.
As fatias _A/_B/_E do GTDB entram juntas via base(). --foco define o grupo de culturas que
confirma o par e onde a copresenca e medida; r_leg/k_leg/ok_leg = grupo foco (padrao leguminosa,
identico a versao anterior); ok_gram = "tambem vale no grupo de contraste".

Saidas: <root>/merged/meta/<saida>/ (parceiros.tsv, estrelas_rizobio.tsv, resumo.txt)
Uso: LD_PRELOAD=... python3 rhizobia_partners.py --root ~/syncom [--centros pgp]
"""
import argparse, itertools, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meta_pairs import load, meta
from consortia import pk, pgp_regex, load_presence, copresence, PGP_CLASSIC, PGP_EXTRA
from shortlist import BIOSAFETY, PLACEHOLDER, base

RIZOBIOS = ["Bradyrhizobium", "Rhizobium", "Mesorhizobium", "Sinorhizobium", "Ensifer",
            "Neorhizobium", "Pararhizobium", "Allorhizobium", "Azorhizobium"]
PGP_HUBS = ["Pseudomonas", "Bacillus", "Priestia", "Paenibacillus", "Azospirillum",
            "Burkholderia", "Paraburkholderia", "Streptomyces", "Azotobacter",
            "Herbaspirillum", "Gluconacetobacter"]
FOSFATO = ["Enterobacter", "Pantoea", "Serratia", "Klebsiella", "Acinetobacter",
           "Rahnella", "Rhizobium"]
PRESETS = {"rizobios": (RIZOBIOS, "leguminosa", "rizobios"),
           "pgp": (PGP_HUBS, "graminea", "pgp"),
           "fosfato": (FOSFATO, "graminea", "fosfato")}

LEGUMINOSAS = {"soja", "feijao/leguminosas"}
GRAMINEAS = {"milho", "trigo", "sorgo", "cana", "arroz"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--centros", default="rizobios")
    ap.add_argument("--foco", choices=["leguminosa", "graminea"], default=None)
    ap.add_argument("--saida", default=None)
    ap.add_argument("--min-k", type=int, default=5)
    ap.add_argument("--max-size", type=int, default=4)
    ap.add_argument("--max-parceiros", type=int, default=12)
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))

    if a.centros in PRESETS:
        centros, foco_def, saida_def = PRESETS[a.centros]
    else:
        centros = [x.strip() for x in a.centros.split(",") if x.strip()]
        foco_def, saida_def = "graminea", "centros_" + centros[0].lower()
    foco = a.foco or foco_def
    contraste = "graminea" if foco == "leguminosa" else "leguminosa"
    saida = a.saida or saida_def

    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, saida); os.makedirs(out, exist_ok=True)
    pgp = pgp_regex(PGP_CLASSIC + PGP_EXTRA)

    m = pd.read_csv(os.path.join(mm, "meta_pairs.tsv.gz"), sep="\t")
    m.index = m.genus_a + "|" + m.genus_b
    rz = {g for g in set(m.genus_a) | set(m.genus_b) if base(g) in set(centros)}
    isrz = m.genus_a.isin(rz) | m.genus_b.isin(rz)
    pos = m[isrz & (m.mu > 0) & (m.fdr < 0.05)]
    NEG = {pk(x, y) for x, y in zip(m.genus_a[(m.mu < 0) & (m.fdr < 0.05)],
                                    m.genus_b[(m.mu < 0) & (m.fdr < 0.05)])}
    ausentes = sorted(set(centros) - {base(g) for g in rz})
    L = [f"Centros pedidos: {a.centros} | grupo foco: {foco} | contraste: {contraste}",
         f"Generos-centro presentes na rede: {len(rz)} -> {', '.join(sorted(rz))}"]
    if ausentes:
        L.append(f"Centros sem nenhuma fatia na rede: {', '.join(ausentes)}")
    L.append(f"Pares centro-parceiro positivos (FDR < 0,05): {len(pos)}")

    d, qc = load(root, False)
    grp = {}
    for s, c in zip(qc.study, qc.cultura.fillna("")):
        grp[s] = "leguminosa" if c in LEGUMINOSAS else ("graminea" if c in GRAMINEAS else "outra")
    d["grupo"] = d.study.map(grp)
    d = d[d.pair.isin(set(pos.index))]
    M = {}
    for g in ("leguminosa", "graminea"):
        v = meta(d[d.grupo == g], a.min_k)
        v.index = v.genus_a + "|" + v.genus_b
        M[g] = v
        L.append(f"{g}{' (foco)' if g == foco else ''}: {d[d.grupo == g].study.nunique()} estudos | "
                 f"{len(v)} pares testaveis (k >= {a.min_k}) | "
                 f"positivos FDR<0,05: {int(((v.mu > 0) & (v.fdr < 0.05)).sum()) if len(v) else 0}")

    fd = os.path.join(mm, "diagnostico_zeros", "diag_pares.tsv")
    diag = pd.read_csv(fd, sep="\t") if os.path.exists(fd) else None
    if diag is not None:
        diag.index = diag.genus_a + "|" + diag.genus_b
    fv = os.path.join(mm, "validacao_temporal", "validacao.tsv")
    val = pd.read_csv(fv, sep="\t", index_col=0) if os.path.exists(fv) else None

    qc2 = qc[qc.study.map(lambda s: grp.get(s) == foco)]
    pres = load_presence(md, set(pos.genus_a) | set(pos.genus_b), qc2)

    ft = "legum" if foco == "leguminosa" else "gram"
    ct = "gram" if foco == "leguminosa" else "legum"
    rows = []
    for k, r in pos.iterrows():
        hub, par = (r.genus_a, r.genus_b) if r.genus_a in rz else (r.genus_b, r.genus_a)
        if par in rz:
            continue
        cp, nst, _ = copresence([hub, par], pres)
        row = {"rizobio": hub, "parceiro": par, "r_global": r.r_mu, "k_global": r.k,
               "copresenca_leg": cp, "estudos_leg_copresente": nst,
               "pgp": bool(pgp.match(par)),
               "sem_cultivado_provavel": bool(PLACEHOLDER.match(par)),
               "biosseguranca": BIOSAFETY.get(base(par), "")}
        for g, tag in (("leguminosa", "legum"), ("graminea", "gram")):
            v = M[g]
            if k in v.index:
                row[f"k_{tag}"] = int(v.loc[k, "k"]); row[f"r_{tag}"] = v.loc[k, "r_mu"]
                row[f"fdr_{tag}"] = v.loc[k, "fdr"]
                row[f"ok_{tag}"] = bool(v.loc[k, "mu"] > 0 and v.loc[k, "fdr"] < 0.05)
            else:
                row[f"k_{tag}"] = 0; row[f"r_{tag}"] = np.nan; row[f"fdr_{tag}"] = np.nan
                row[f"ok_{tag}"] = False
        row["k_leg"] = row[f"k_{ft}"]; row["r_leg"] = row[f"r_{ft}"]
        row["fdr_leg"] = row[f"fdr_{ft}"]; row["ok_leg"] = row[f"ok_{ft}"]
        row["ok_gram"] = row[f"ok_{ct}"]
        row["diag"] = diag.loc[k, "diag"] if diag is not None and k in diag.index else ""
        row["replicado"] = (bool(val.loc[k, "replicado"]) if val is not None and k in val.index else None)
        rows.append(row)
    if not rows:
        L.append("\nNenhum par centro-parceiro positivo: nada a analisar.")
        txt = "\n".join(L); open(os.path.join(out, "resumo.txt"), "w").write(txt); print(txt); return
    P = pd.DataFrame(rows).sort_values(["rizobio", "r_leg"], ascending=[True, False])
    P.to_csv(os.path.join(out, "parceiros.tsv"), sep="\t", index=False, float_format="%.4g")

    L += ["", f"Parceiros confirmados em {foco}s, por centro (r_leg = r no grupo foco; "
              f"* = tambem em {contraste}s; PGP; repl = replicado no tempo):"]
    conf = P[P.ok_leg]
    for hub, g in conf.groupby("rizobio"):
        L.append(f"  {hub} ({len(g)} parceiros confirmados de {int((P.rizobio == hub).sum())} testados):")
        for _, r in g.head(12).iterrows():
            tags = []
            if r.ok_gram: tags.append(f"*{contraste}s")
            if r.pgp: tags.append("PGP")
            if r.replicado: tags.append("repl")
            if r["diag"]: tags.append(r["diag"])
            if r.biosseguranca: tags.append("biosseg")
            if r.sem_cultivado_provavel: tags.append("sem cultivado")
            L.append(f"    {r.parceiro}: r_leg={r.r_leg:.2f} (k={int(r.k_leg)}) "
                     f"copres={r.copresenca_leg:.2f} em {int(r.estudos_leg_copresente)} estudos"
                     + (" | " + ", ".join(tags) if tags else ""))
    if conf.empty:
        L.append(f"  nenhum parceiro confirmado em {foco}s com os criterios atuais")

    est = []
    for hub, g in conf.groupby("rizobio"):
        cand = list(g.sort_values("r_leg", ascending=False).parceiro[:a.max_parceiros])
        for n in range(1, a.max_size):
            for comb in itertools.combinations(cand, n):
                mem = [hub] + list(comb)
                if any(pk(x, y) in NEG for x, y in itertools.combinations(mem, 2)):
                    continue
                sub = g[g.parceiro.isin(comb)]
                cp, nst, crops = copresence(mem, pres)
                est.append({"tamanho": len(mem), "rizobio": hub, "parceiros": " + ".join(sorted(comb)),
                            "r_leg_min": sub.r_leg.min(), "k_leg_min": int(sub.k_leg.min()),
                            "copresenca_leg": cp, "estudos_leg": nst,
                            "n_pgp": int(sub.pgp.sum()), "todos_em_gramineas": bool(sub.ok_gram.all()),
                            "n_replicados": int(sub.replicado.fillna(False).sum()),
                            "biosseguranca": "; ".join(x for x in sub.biosseguranca if x),
                            "sem_cultivado": "; ".join(sub.parceiro[sub.sem_cultivado_provavel])})
    E = pd.DataFrame(est)
    if len(E):
        E["apto"] = np.where(E.sem_cultivado != "", "nao (sem cultivado)",
                             np.where(E.biosseguranca != "", "sim, com triagem", "sim"))
        E = E.sort_values(["tamanho", "r_leg_min", "copresenca_leg"], ascending=[True, False, False])
        E.to_csv(os.path.join(out, "estrelas_rizobio.tsv"), sep="\t", index=False, float_format="%.4g")
        L += ["", f"Consorcios com o centro no meio (top 5 por tamanho, elo mais fraco em {foco}s):"]
        for t, sub in E.groupby("tamanho"):
            L.append(f"  tamanho {t} ({len(sub)} combinacoes):")
            for _, r in sub.head(5).iterrows():
                L.append(f"    [{r.rizobio}] + {r.parceiros} | r_min={r.r_leg_min:.2f} "
                         f"copres={r.copresenca_leg:.2f} ({int(r.estudos_leg)} estudos) | {r.apto}"
                         + (f" | tambem em {contraste}s" if r.todos_em_gramineas else ""))
            best = sub[sub.apto == "sim"]
            if len(best):
                b = best.iloc[0]
                L.append(f"    melhor sem bandeira: [{b.rizobio}] + {b.parceiros} "
                         f"(r_min={b.r_leg_min:.2f}, copres={b.copresenca_leg:.2f})")
    txt = "\n".join(L)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
