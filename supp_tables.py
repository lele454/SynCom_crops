#!/usr/bin/env python3
# [S6.7] Supplementary tables
# Supplementary 6 of the manuscript. See README.md for the full map.
"""Supplementary study tables, in English (S1: studies; S2: processing).

S1 cobre o funil inteiro: os 945 estudos do levantamento, com o motivo de saida de
cada um, e detalhe completo para os que chegaram ao processamento. O status separa
os estudos que produziram estimativas de efeito dos que entram apenas nas medidas
descritivas, para que a contagem da tabela bata com a do texto.

Uso: python3 supp_tables.py --root ~/syncom --studies subset_studies.tsv --out Tabelas"""
import argparse, os, re
import numpy as np
import pandas as pd

COLS = {"study": "study", "country": "country", "cultura": "crop", "grupo": "host_group",
        "region": "region_16S", "tipo_de_leitura": "read_type", "fragmento": "fragment",
        "amostras_proc": "samples_processed", "amostras_ok": "samples_passing_qc",
        "amostras_na_analise": "samples_analysed", "riqueza_mediana": "median_richness",
        "frac_genero": "fraction_reads_with_genus", "retencao_dada2": "read_retention",
        "covariaveis": "n_covariates", "fatores_de_desenho": "design_factors",
        "origem_dos_fatores": "design_factor_source", "generos_testados": "genera_tested",
        "pares": "pairs_tested", "status": "status", "motivo": "exclusion_reason",
        "n_smp": "samples_submitted", "regiao": "region_submitted",
        "tipos_tratamento": "treatment_types", "fatores_tratamento": "treatment_fields",
        "host": "host", "study_title": "study_title", "first_public": "first_public",
        "modo": "read_mode", "unidade": "sample_unit", "amostras": "samples",
        "frac_quimeras": "chimera_fraction", "retencao": "read_retention"}
CROP = {"milho": "maize", "trigo": "wheat", "arroz": "rice", "sorgo": "sorghum",
        "cana": "sugarcane", "soja": "soybean", "feijao/leguminosas": "other legumes"}
GRP = {"gramínea": "grass", "graminea": "grass", "leguminosa": "legume"}
STATUS = {"incluido": "included", "excluido": "excluded"}
EFEITOS = "included (effect sizes)"
DESCRITIVO = "included (descriptive measures only)"
TRIAGEM = "screened out"
SEM_EFEITO = "fewer than 20 samples after sample curation; contributes to descriptive measures only"


def motivo_triagem(region):
    """Por que o estudo do levantamento nao chegou ao processamento."""
    r = str(region) if region is not None else ""
    if "ITS" in r:
        return "ITS amplicon (fungal), not 16S"
    if r == "V5-V7":
        return "16S region without overlap with the V4 window"
    if r.startswith("erro") or r in ("sem_primer", "nao_inferido"):
        return "amplified region could not be determined from the reads"
    return "host crop not among the seven studied, or study not retained at crop curation"
READ = {"V4 completa": "complete V4", "fragmento (leitura direta)": "forward-read fragment"}
FRAG = {"V4 completa": "complete V4", "V3 (R1)": "V3 (forward read)",
        "inicio da V4 (R1)": "V4 start (forward read)"}
SRC = {"automatico": "automatic", "curado": "curated"}
MODE = {"PE": "paired", "R1": "forward only"}
TRT = {"aduba_mineral": "mineral fertilisation", "aduba_organica": "organic amendment",
       "seca_agua": "water/drought", "genotipo_planta": "plant genotype", "patogeno": "pathogen",
       "rotacao_cultivo": "crop rotation", "contaminacao": "contamination",
       "manejo_org_conv": "organic vs. conventional", "salinidade": "salinity",
       "temperatura_co2": "temperature/CO2", "preparo_solo": "tillage", "nenhum": "none"}
WORD = {"tratamento": "treatment", "genotipo": "genotype", "tempo": "time", "local": "site",
        "solo": "soil", "curado": "curated", "especie": "species", "seca": "drought",
        "nivel": "level", "regime": "regime", "bloco": "block", "subbloco": "subblock",
        "parcela": "plot", "subpopulacao": "subpopulation", "nitrogenio": "nitrogen",
        "inoculo": "inoculum", "cultura": "crop", "cultivo": "cropping", "corrida": "run",
        "preservacao": "preservation", "compartimento": "compartment", "ciclo": "cycle",
        "grupo": "group", "sistema": "system", "dose": "dose", "cultivar": "cultivar",
        "salinidade": "salinity", "biosample": "biosample", "nenhum": "none",
        "a_confirmar": "to confirm"}


def traduz_fatores(v):
    if not isinstance(v, str) or not v:
        return v
    v = re.sub(r"\bcur_([a-zA-Z]+)", lambda m: "curated_" + WORD.get(m.group(1), m.group(1)), v)
    v = re.sub(r"\[([a-z]+)\]", lambda m: "[" + WORD.get(m.group(1), m.group(1)) + "]", v)
    return WORD.get(v, v)


def traduz_motivo(v):
    if not isinstance(v, str) or not v:
        return v
    v = re.sub(r"<(\d+) amostras com >= (\d+) leituras", r"<\1 samples with >= \2 reads", v)
    v = re.sub(r"riqueza mediana ([\d.]+) < (\d+).*", r"median richness \1 < \2 (simplified community)", v)
    v = re.sub(r"so (\d+)% das leituras com genero", r"only \1% of reads assigned to a genus", v)
    return v


def ler(p, **kw):
    return pd.read_csv(p, sep="\t", **kw) if os.path.exists(p) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--studies", default="subset_studies.tsv")
    ap.add_argument("--survey", default=os.path.join("survey", "studies.tsv"),
                    help="levantamento completo (945 estudos), para as linhas de triagem")
    ap.add_argument("--out", default="Tabelas")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root)); md = os.path.join(root, "merged")
    out = os.path.abspath(os.path.expanduser(a.out)); os.makedirs(out, exist_ok=True)

    qc = ler(os.path.join(md, "study_qc.tsv"), dtype=str)
    if qc is None:
        raise SystemExit("study_qc.tsv nao encontrado")
    sub = ler(a.studies, dtype=str)
    proc = ler(os.path.join(md, "study_processing.tsv"), dtype=str)
    dates = ler(os.path.join(md, "study_dates.tsv"), dtype=str)
    des = ler(os.path.join(md, "effects", "study_design.tsv"), dtype=str)

    T = qc.copy()
    if sub is not None:
        extra = [c for c in sub.columns if c not in T.columns and c != "study_accession"]
        T = T.merge(sub[["study_accession"] + extra].rename(columns={"study_accession": "study"}),
                    on="study", how="left")
    if dates is not None and "study" in dates.columns:
        T = T.merge(dates, on="study", how="left")
    if des is not None:
        keep = [c for c in ("n", "k", "generos_testados", "pares", "fatores", "fonte") if c in des.columns]
        T = T.merge(des[["study"] + keep].rename(columns={
            "n": "amostras_na_analise", "k": "covariaveis", "fatores": "fatores_de_desenho",
            "fonte": "origem_dos_fatores"}), on="study", how="left")
    if "modo" in T.columns:
        T["tipo_de_leitura"] = T["modo"].map(READ_PT := {"PE": "V4 completa",
                                                         "R1": "fragmento (leitura direta)"}).fillna("")

    for c, m in (("cultura", CROP), ("grupo", GRP), ("status", STATUS),
                 ("tipo_de_leitura", READ), ("fragmento", FRAG), ("origem_dos_fatores", SRC)):
        if c in T.columns:
            T[c] = T[c].map(lambda v: m.get(v, v))
    if "motivo" in T.columns:
        T["motivo"] = T["motivo"].map(traduz_motivo)
    if "fatores_de_desenho" in T.columns:
        T["fatores_de_desenho"] = T["fatores_de_desenho"].map(traduz_fatores)
    if "tipos_tratamento" in T.columns:
        T["tipos_tratamento"] = T["tipos_tratamento"].map(
            lambda v: "; ".join(TRT.get(x, x) for x in str(v).split(";")) if isinstance(v, str) and v else v)
    for c in ("region", "regiao"):
        if c in T.columns:
            T[c] = T[c].map(lambda v: WORD.get(v, v))

    ordem = ["study", "country", "cultura", "grupo", "host", "region", "tipo_de_leitura", "fragmento",
             "amostras_proc", "amostras_ok", "amostras_na_analise", "riqueza_mediana", "frac_genero",
             "retencao_dada2", "covariaveis", "fatores_de_desenho", "origem_dos_fatores",
             "tipos_tratamento", "generos_testados", "pares", "status", "motivo",
             "study_title", "first_public"]
    cols = [c for c in ordem if c in T.columns] + [c for c in T.columns
                                                   if c not in ordem and c != "modo"]
    T = T[cols].rename(columns=COLS)

    # status em tres niveis: o texto reporta 123 estudos com estimativas de efeito e
    # 125 aprovados; sem esta distincao a tabela parece contradizer o texto
    if "status" in T.columns:
        com_efeito = pd.Series(False, index=T.index)
        if des is not None and "study" in des.columns:
            com_efeito = T.study.isin(set(des.study.dropna()))
        elif "samples_analysed" in T.columns:
            com_efeito = T.samples_analysed.notna() & (T.samples_analysed.astype(str) != "")
        inc = T.status == "included"
        T.loc[inc & com_efeito, "status"] = EFEITOS
        T.loc[inc & ~com_efeito, "status"] = DESCRITIVO
        if "exclusion_reason" in T.columns:
            vazio = T.exclusion_reason.isna() | (T.exclusion_reason.astype(str).str.strip() == "")
            T.loc[(T.status == DESCRITIVO) & vazio, "exclusion_reason"] = SEM_EFEITO

    # linhas de triagem: os estudos do levantamento que nao chegaram ao processamento
    S = ler(a.survey, dtype=str)
    if S is not None and "study_accession" in S.columns:
        fora = S[~S.study_accession.isin(set(T.study))].copy()
        R = pd.DataFrame({
            "study": fora.study_accession,
            "country": fora.get("countries", pd.Series(dtype=str)),
            "host": fora.get("hosts", pd.Series(dtype=str)),
            "region_submitted": fora.get("region", pd.Series(dtype=str)),
            "samples_submitted": fora.get("n_runs", pd.Series(dtype=str)),
            "study_title": fora.get("study_title", pd.Series(dtype=str)),
            "first_public": fora.get("first_public", pd.Series(dtype=str)),
            "status": TRIAGEM,
            "exclusion_reason": fora.get("region", pd.Series(dtype=str)).map(motivo_triagem),
        })
        T = pd.concat([T, R.reindex(columns=T.columns)], ignore_index=True)

    # referencias dos trabalhos originais, se fetch_refs.py ja rodou
    refs = ler(os.path.join(out, "study_refs.tsv"), dtype=str)
    if refs is not None and "study" in refs.columns:
        keep = [c for c in ("doi", "pmid", "ref_title", "ref_journal", "ref_year", "ref_source")
                if c in refs.columns]
        T = T.merge(refs[["study"] + keep], on="study", how="left")

    ordem_status = {EFEITOS: 0, DESCRITIVO: 1, "excluded": 2, TRIAGEM: 3}
    T = T.assign(_o=T.status.map(lambda v: ordem_status.get(v, 9))).sort_values(
        ["_o", "study"]).drop(columns="_o")
    T.to_csv(os.path.join(out, "TableS1_studies.tsv"), sep="\t", index=False)

    if proc is not None:
        P = proc.copy()
        if "modo" in P.columns:
            P["modo"] = P["modo"].map(lambda v: MODE.get(v, v))
        if "fragmento" in P.columns:
            P["fragmento"] = P["fragmento"].map(lambda v: FRAG.get(v, v))
        if "unidade" in P.columns:
            P["unidade"] = P["unidade"].map(lambda v: WORD.get(v, v))
        cols2 = [c for c in ("study", "region", "modo", "fragmento", "unidade", "amostras",
                             "maxEE", "truncF", "truncR", "frac_quimeras", "retencao")
                 if c in P.columns]
        P[cols2].sort_values("study").rename(columns=COLS).to_csv(
            os.path.join(out, "TableS2_processing.tsv"), sep="\t", index=False)

    print(f"TableS1: {len(T)} studies")
    if "status" in T.columns:
        for v, n in T.status.value_counts().items():
            print(f"  {v}: {n}")
        if (T.status == TRIAGEM).any():
            print("  motivos de triagem:")
            for v, n in T.loc[T.status == TRIAGEM, "exclusion_reason"].value_counts().items():
                print(f"    {n:4d}  {v}")
    print("columns: " + ", ".join(T.columns))
    if proc is not None:
        print(f"TableS2: {len(proc)} studies")
    print(f"Output in {out}")


if __name__ == "__main__":
    main()
