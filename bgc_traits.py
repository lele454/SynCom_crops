#!/usr/bin/env python3
# [S5.3] Secondary-metabolite gene clusters and NIS siderophores
# Supplementary 5 of the manuscript. See README.md for the full map.
"""
Triangulacao funcional, etapa 3: clusters biossinteticos (BGCs) e sideroforos.

- GECCO (Carroll et al.) em cada genoma: BGCs e tipo previsto (NRP, Polyketide, RiPP, Terpene,
  Saccharide, Alkaloid, Unknown)
- Sideroforos independentes de NRPS: proteinas com dominio IucA_IucC (PF04183, sintetase NIS:
  aerobactina, desferrioxamina, achromobactina...), da anotacao de dominios do proprio GECCO
- Por genero: media de BGCs por tipo e fracao de genomas com sintetase NIS
- Por consorcio: tipos de BGC do centro, da uniao e complementares; sideroforo em algum membro

Requer: pip install --user gecco-tool
Saidas: <root>/genomas/bgc/ (bgc_por_genoma.tsv, bgc_por_genero.tsv, bgc_por_consorcio.tsv, resumo.txt)
Uso: python3 bgc_traits.py --root ~/syncom [--cpus 4] [--por-genero 5]
"""
import argparse, glob, os, shutil, subprocess, sys
import pandas as pd

TYPES = ["NRP", "Polyketide", "RiPP", "Terpene", "Saccharide", "Alkaloid", "Unknown"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument("--por-genero", type=int, default=5, help="genomas por genero (menos = mais rapido)")
    a = ap.parse_args()
    if not shutil.which("gecco"):
        sys.exit("gecco nao encontrado: pip install --user gecco-tool  (e ~/.local/bin no PATH)")
    root = os.path.abspath(os.path.expanduser(a.root))
    gd = os.path.join(root, "genomas"); od = os.path.join(gd, "bgc"); os.makedirs(od, exist_ok=True)
    sel = pd.read_csv(os.path.join(gd, "selecionados.tsv"), sep="\t")
    sel = sel.groupby("genus", group_keys=False).head(a.por_genero)

    rows = []
    for i, (_, r) in enumerate(sel.iterrows(), 1):
        fna = glob.glob(os.path.join(gd, r.genus, f"{r.acc}*.fna"))
        if not fna:
            continue
        out = os.path.join(od, r.acc)
        done = os.path.join(out, "FEITO")
        if not os.path.exists(done):
            print(f"[{i}/{len(sel)}] GECCO {r.genus} {r.acc}", flush=True)
            res = subprocess.run(["gecco", "-q", "run", "--genome", fna[0], "-o", out, "-j", str(a.cpus),
                                  "--force-tsv"], capture_output=True, text=True)
            if res.returncode != 0:
                print("  falhou:", res.stderr.strip()[-300:]); continue
            open(done, "w").close()
        cl = glob.glob(os.path.join(out, "*.clusters.tsv"))
        ft = glob.glob(os.path.join(out, "*.features.tsv"))
        row = {"genus": r.genus, "acc": r.acc}
        c = pd.read_csv(cl[0], sep="\t") if cl and os.path.getsize(cl[0]) > 0 else pd.DataFrame()
        row["BGC_total"] = len(c)
        for t in TYPES:
            row[f"BGC_{t}"] = int(c["type"].fillna("Unknown").str.split(";").map(lambda x: t in x).sum()) if len(c) else 0
        nis = 0
        if ft and os.path.getsize(ft[0]) > 0:
            f = pd.read_csv(ft[0], sep="\t")
            hit = f[(f.domain.astype(str).str.contains("PF04183|IucA_IucC")) & (f.i_evalue < 1e-5)]
            nis = hit.protein_id.nunique()
        row["sideroforo_NIS"] = nis
        rows.append(row)
    B = pd.DataFrame(rows)
    B.to_csv(os.path.join(od, "bgc_por_genoma.tsv"), sep="\t", index=False)

    num = [c for c in B.columns if c.startswith("BGC_")]
    G = B.groupby("genus")[num].mean()
    G["frac_sideroforo_NIS"] = B.assign(s=B.sideroforo_NIS > 0).groupby("genus").s.mean()
    G.insert(0, "n_genomas", B.groupby("genus").size())
    G.to_csv(os.path.join(od, "bgc_por_genero.tsv"), sep="\t", float_format="%.2f")

    sl = pd.read_csv(os.path.join(root, "merged", "meta", "consorcios", "shortlist.tsv"), sep="\t")
    types = lambda g: {t for t in TYPES if t != "Unknown" and g in G.index and G.loc[g, f"BGC_{t}"] >= 0.5}
    C = []
    for _, r in sl.iterrows():
        hub, parts = r.centro, r.parceiros.split(" + ")
        th = types(hub); tp = set().union(*[types(p) for p in parts])
        sid = [g for g in [hub] + parts if g in G.index and G.loc[g, "frac_sideroforo_NIS"] >= 0.5]
        C.append({"id": r.id, "centro": hub, "parceiros": r.parceiros,
                  "BGC_medio_total": sum(G.loc[g, "BGC_total"] for g in [hub] + parts if g in G.index),
                  "tipos_centro": ";".join(sorted(th)), "tipos_uniao": ";".join(sorted(th | tp)),
                  "n_tipos_uniao": len(th | tp), "tipos_complementares": ";".join(sorted(tp - th)),
                  "sideroforo_NIS_em": ";".join(sid)})
    C = pd.DataFrame(C).sort_values(["n_tipos_uniao", "BGC_medio_total"], ascending=False)
    C.to_csv(os.path.join(od, "bgc_por_consorcio.tsv"), sep="\t", index=False, float_format="%.1f")

    L = [f"Genomas analisados: {len(B)}", "", "BGCs medios por genoma, por genero (GECCO):",
         G.to_string(float_format=lambda v: f"{v:.2f}"), "", "Consorcios:",
         C.to_string(index=False, float_format=lambda v: f"{v:.1f}")]
    txt = "\n".join(L)
    open(os.path.join(od, "resumo.txt"), "w").write(txt)
    print("\n" + txt)


if __name__ == "__main__":
    main()
