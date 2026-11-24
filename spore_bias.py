#!/usr/bin/env python3
# [S3.15] Endospore hypothesis
# Supplementary 3 of the manuscript. See README.md for the full map.
"""
Hipotese dos endosporos: generos que formam endosporos aparecem menos nas redes de coocorrencia
de 16S porque a leitura conta esporos dormentes junto com celulas ativas, de modo que a
abundancia deixa de acompanhar a atividade e as correlacoes se diluem.

Grupos (por familia, de merged/taxonomy_v4.tsv + taxonomy_r1.tsv):
  endosporo       Bacillota: Bacillaceae, Paenibacillaceae, Planococcaceae, Clostridiaceae...
  bacillota_nao   Bacillota sem endosporo (Lactobacillaceae, Staphylococcaceae...) - contraste
                  que separa "ser Bacillota" de "formar endosporo"
  esporo_aereo    Actinomycetota com esporos aereos (Streptomycetaceae, Micromonosporaceae...)
                  - esporos menos persistentes: a hipotese preve efeito menor
  outros          o resto

Analises (desfecho = parceiros com associacao positiva significativa; offset = pares testados):
  1. Poisson quasi: grau_pos ~ endosporo + esporo_aereo + log_abund + ocupacao
  2. Idem com dummies dos --filos filos mais comuns (o efeito sobrevive ao controle de filo?)
  3. So dentro de Bacillota: endosporo vs nao endosporo
  4. Pareado por estratos (abundancia x ocupacao): razao da taxa observada/esperada
  5. Forca tipica da associacao (|r| mediano dos pares do genero) e ocupacao/abundancia por grupo
     - se os esporos diluem o sinal, |r| cai sem que a prevalencia caia

Saidas: <root>/merged/meta/vies_esporos/ (por_genero.tsv, resumo.txt)
Uso: LD_PRELOAD=... python3 spore_bias.py --root ~/syncom
"""
import argparse, os, re, sys
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inoculant_bias import irls, base

ENDOSPORO = ["Bacillaceae", "Paenibacillaceae", "Planococcaceae", "Sporolactobacillaceae",
             "Alicyclobacillaceae", "Thermoactinomycetaceae", "Clostridiaceae", "Peptostreptococcaceae",
             "Sporomusaceae", "Symbiobacteriaceae", "Desulfitobacteriaceae", "Heliobacteriaceae",
             "Syntrophomonadaceae", "Caldicoprobacteraceae", "Thermoanaerobacteraceae",
             "Sulfobacillaceae", "Bacillaceae_H", "Bacillaceae_G", "Bacillaceae_I",
             "DTU015", "Desulfuribacillaceae", "Pasteuriaceae"]
BACILLOTA_NAO = ["Lactobacillaceae", "Staphylococcaceae", "Streptococcaceae", "Enterococcaceae",
                 "Erysipelotrichaceae", "Carnobacteriaceae", "Leuconostocaceae", "Listeriaceae",
                 "Acholeplasmataceae", "Veillonellaceae", "Selenomonadaceae"]
ESPORO_AEREO = ["Streptomycetaceae", "Micromonosporaceae", "Pseudonocardiaceae", "Streptosporangiaceae",
                "Nocardiopsaceae", "Thermomonosporaceae", "Kineosporiaceae", "Frankiaceae",
                "Actinoplanaceae", "Glycomycetaceae"]
BACILLOTA = re.compile(r"^(Bacillota|Firmicutes)", re.I)


def tax_map(md):
    """genero -> (filo, familia) pelo valor mais frequente nas tabelas de taxonomia."""
    fs = [os.path.join(md, f) for f in ("taxonomy_v4.tsv", "taxonomy_r1.tsv")]
    d = pd.concat([pd.read_csv(f, sep="\t", usecols=["Phylum", "Family", "Genus"], dtype=str)
                   for f in fs if os.path.exists(f)], ignore_index=True)
    d = d.dropna(subset=["Genus"])
    mode = lambda s: s.dropna().mode().iloc[0] if s.notna().any() else np.nan
    return d.groupby("Genus").agg(filo=("Phylum", mode), familia=("Family", mode))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--filos", type=int, default=8, help="filos com dummy no modelo 2")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    md = os.path.join(root, "merged"); mm = os.path.join(md, "meta")
    out = os.path.join(mm, "vies_esporos"); os.makedirs(out, exist_ok=True)

    m = pd.read_csv(os.path.join(mm, "meta_pairs.tsv.gz"), sep="\t",
                    usecols=["genus_a", "genus_b", "mu", "fdr", "r_mu"])
    m["absr"] = m.r_mu.abs()
    pos = m[(m.mu > 0) & (m.fdr < 0.05)]
    cnt = lambda d: pd.concat([d.genus_a, d.genus_b]).value_counts()
    longo = pd.concat([m[["genus_a", "absr"]].rename(columns={"genus_a": "g"}),
                       m[["genus_b", "absr"]].rename(columns={"genus_b": "g"})])
    G = pd.DataFrame({"n_testados": cnt(m), "grau_pos": cnt(pos)}).fillna(0)
    G["absr_mediano"] = longo.groupby("g").absr.median()
    G = G[G.n_testados >= 10]
    occ = pd.read_csv(os.path.join(md, "genus_occupancy.tsv"), sep="\t").set_index("genus")
    G = G.join(occ[["ocupacao", "abund_media"]]).join(tax_map(md))
    G = G.dropna(subset=["ocupacao", "abund_media"])
    G["log_abund"] = np.log10(G.abund_media.clip(lower=1e-8))
    fam = G.familia.fillna("").map(base)
    G["grupo"] = np.where(fam.isin({base(x) for x in ENDOSPORO}), "endosporo",
                 np.where(fam.isin({base(x) for x in ESPORO_AEREO}), "esporo_aereo",
                 np.where(fam.isin({base(x) for x in BACILLOTA_NAO}), "bacillota_nao",
                 np.where(G.filo.fillna("").str.match(BACILLOTA), "bacillota_outras", "outros"))))
    G["endosporo"] = (G.grupo == "endosporo").astype(int)
    G["aereo"] = (G.grupo == "esporo_aereo").astype(int)
    G["taxa_pos"] = G.grau_pos / G.n_testados
    G.to_csv(os.path.join(out, "por_genero.tsv"), sep="\t", float_format="%.5g")

    L = [f"Generos testados (>= 10 pares): {len(G)}",
         "Grupos: " + " | ".join(f"{k} {v}" for k, v in G.grupo.value_counts().items()), "",
         "Taxa de associacao positiva (grau/pares testados) por grupo:",
         G.groupby("grupo").agg(n=("taxa_pos", "size"), taxa_mediana=("taxa_pos", "median"),
                                absr_mediano=("absr_mediano", "median"),
                                ocupacao_mediana=("ocupacao", "median"),
                                abund_mediana=("abund_media", "median"))
          .to_string(float_format=lambda v: f"{v:.4g}"), ""]

    cols = ["endosporo", "aereo", "log_abund", "ocupacao"]
    tab, phi = irls(G.grau_pos.values.astype(float), G[cols].values, G.n_testados.values, True, cols)
    L += [f"1. Poisson com offset log(pares testados) (dispersao quasi = {phi:.1f}):",
          tab.to_string(float_format=lambda v: f"{v:.3g}"), ""]

    top = G.filo.value_counts().head(a.filos).index
    dum = pd.get_dummies(G.filo.where(G.filo.isin(top), "outro_filo"), prefix="f", drop_first=True).astype(float)
    cols2 = cols + list(dum.columns)
    X2 = np.column_stack([G[cols].values, dum.values])
    tab2, phi2 = irls(G.grau_pos.values.astype(float), X2, G.n_testados.values, True, cols2)
    L += [f"2. Idem, com dummies de filo (dispersao quasi = {phi2:.1f}):",
          tab2.loc[[i for i in tab2.index if not i.startswith("f_")]].to_string(float_format=lambda v: f"{v:.3g}"),
          f"   (dummies de filo omitidas; {len(dum.columns)} filos no modelo)", ""]

    B = G[G.grupo.isin(["endosporo", "bacillota_nao", "bacillota_outras"])]
    L.append("3. So dentro de Bacillota (endosporo vs demais Bacillota):")
    if B.endosporo.nunique() == 2 and len(B) >= 15:
        c3 = ["endosporo", "log_abund", "ocupacao"]
        t3, p3 = irls(B.grau_pos.values.astype(float), B[c3].values, B.n_testados.values, True, c3)
        L += [f"   n = {len(B)} (endosporo {int(B.endosporo.sum())}) | dispersao quasi = {p3:.1f}",
              t3.to_string(float_format=lambda v: f"{v:.3g}")]
    else:
        L.append(f"   contraste insuficiente (n = {len(B)}, endosporo = {int(B.endosporo.sum())})")
    L.append("")

    G["estrato"] = (pd.qcut(G.log_abund.rank(method="first"), 5, labels=False).astype(str) + "_" +
                    pd.qcut(G.ocupacao.rank(method="first"), 3, labels=False).astype(str))
    e = G.groupby("estrato").taxa_pos.mean()
    G["esperado"] = G.estrato.map(e)
    agg = G.groupby("grupo").apply(lambda d: pd.Series(
        {"n": len(d), "taxa_obs": d.taxa_pos.mean(), "taxa_esperada": d.esperado.mean()}),
        include_groups=False)
    agg["razao_obs_esp"] = agg.taxa_obs / agg.taxa_esperada
    L += ["4. Pareado por estratos (abundancia x ocupacao):",
          agg.to_string(float_format=lambda v: f"{v:.4g}")]
    G["pct_no_estrato"] = G.groupby("estrato").taxa_pos.rank(pct=True)
    for g in ("endosporo", "esporo_aereo"):
        p = G.loc[G.grupo == g, "pct_no_estrato"]
        if len(p) >= 6:
            L.append(f"   {g}: percentil mediano = {p.median():.2f} | Wilcoxon vs 0,5: "
                     f"p = {stats.wilcoxon(p - 0.5).pvalue:.3g} (n = {len(p)})")
    L += ["", "5. Generos endosporulantes com mais pares testados (os 25 primeiros):",
          G[G.grupo == "endosporo"].nlargest(25, "n_testados")[
              ["familia", "n_testados", "grau_pos", "taxa_pos", "pct_no_estrato", "absr_mediano",
               "ocupacao", "abund_media"]].to_string(float_format=lambda v: f"{v:.4g}")]
    G.to_csv(os.path.join(out, "por_genero.tsv"), sep="\t", float_format="%.5g")
    txt = "\n".join(L)
    open(os.path.join(out, "resumo.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
