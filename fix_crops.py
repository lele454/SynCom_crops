#!/usr/bin/env python3
# [S1.4] Removal of out-of-scope crops: reassignment
# Supplementary 1 of the manuscript. See README.md for the full map.
"""Correcao das culturas: exclui estudos e amostras fora do escopo e registra o hospedeiro
como fator de desenho nos estudos que misturam duas culturas-alvo.

  --aplicar  grava as mudancas; sem ele, so mostra o que faria.
"""
import argparse, csv, os, shutil
from collections import defaultdict

ESTUDOS_FORA = {
    "PRJNA1147423": "trigo-sarraceno (Fagopyrum), fora das culturas do estudo",
    "PRJNA1028264": "milheto (Panicum miliaceum), fora das culturas do estudo",
    "PRJNA1184995": "milheto (Panicum miliaceum), fora das culturas do estudo",
    "PRJNA1514812": "metade pak choi; restariam 15 amostras de milho, abaixo do minimo",
}
AMOSTRAS_FORA = {"PRJNA506457": ["rye"]}          # valores de host a excluir (minusculas)
MISTOS = {"PRJNA1056533": "misto_cereal", "PRJNA1366529": "misto_cereal"}
FATOR_MISTO = "host"


def attrs_dir():
    for d in ("stage2/attrs", os.path.expanduser("~/syncom/stage2/attrs")):
        if os.path.isdir(d):
            return d
    raise SystemExit("pasta de atributos nao encontrada")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="~/syncom")
    ap.add_argument("--runs", default="runs_subset.tsv")
    ap.add_argument("--studies", default="subset_studies.tsv")
    ap.add_argument("--design", default="design_curado.tsv")
    ap.add_argument("--aplicar", action="store_true")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root)); md = os.path.join(root, "merged")
    ad = attrs_dir()
    exc_f = os.path.join(md, "curadoria_final", "excluir_amostras.tsv")

    bs2runs = defaultdict(list)
    por_estudo = defaultdict(list)
    with open(a.runs) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            bs2runs[r["sample_accession"]].append(r["run_accession"])
            st = r.get("study_accession") or r.get("study") or ""
            por_estudo[st].append(r["run_accession"])

    novos = []
    for s, motivo in ESTUDOS_FORA.items():
        runs = por_estudo.get(s, [])
        novos += [(s, x, motivo) for x in runs]
        print(f"{s}: estudo inteiro, {len(runs)} corridas | {motivo}")

    for s, valores in AMOSTRAS_FORA.items():
        f = os.path.join(ad, f"{s}.tsv")
        alvo = set()
        for r in csv.DictReader(open(f), delimiter="\t"):
            if r["key"].lower() == "host" and r["value"].strip().lower() in valores:
                alvo.add(r["sample"])
        runs = [x for b in alvo for x in bs2runs.get(b, [])]
        novos += [(s, x, f"hospedeiro fora do escopo: {', '.join(valores)}") for x in runs]
        print(f"{s}: {len(alvo)} amostras ({len(runs)} corridas) com hospedeiro {valores}")

    ja = set()
    if os.path.exists(exc_f):
        for r in csv.DictReader(open(exc_f), delimiter="\t"):
            ja.add((r["study"], r["run"]))
    add = [t for t in novos if (t[0], t[1]) not in ja]
    print(f"\nA acrescentar em excluir_amostras.tsv: {len(add)} corridas "
          f"(ja havia {len(ja)})")

    D = {}
    if os.path.exists(a.design):
        for r in csv.DictReader(open(a.design), delimiter="\t"):
            D[r["study"]] = r["fatores"]
    for s in MISTOS:
        atual = [x.strip() for x in D.get(s, "").split(";") if x.strip()]
        if FATOR_MISTO not in atual:
            atual.append(FATOR_MISTO)
        D[s] = "; ".join(atual)
        print(f"{s}: fatores curados -> {D[s]}")

    if not a.aplicar:
        print("\n(ensaio; use --aplicar para gravar)")
        return

    os.makedirs(os.path.dirname(exc_f), exist_ok=True)
    novo = not os.path.exists(exc_f)
    if not novo:
        shutil.copy(exc_f, exc_f + ".bak")
    with open(exc_f, "a", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        if novo:
            w.writerow(["study", "run", "motivo"])
        for t in add:
            w.writerow(t)

    if os.path.exists(a.design):
        shutil.copy(a.design, a.design + ".bak")
    with open(a.design, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["study", "fatores"])
        for s, v in sorted(D.items()):
            w.writerow([s, v])

    shutil.copy(a.studies, a.studies + ".bak")
    linhas = list(csv.DictReader(open(a.studies), delimiter="\t"))
    campos = list(linhas[0])
    for r in linhas:
        if r["study_accession"] in MISTOS and "cultura" in r:
            r["cultura"] = MISTOS[r["study_accession"]]
    with open(a.studies, "w", newline="") as fh:
        w = csv.DictWriter(fh, campos, delimiter="\t"); w.writeheader(); w.writerows(linhas)
    print("\nGravado. Copias de seguranca com sufixo .bak")


if __name__ == "__main__":
    main()
