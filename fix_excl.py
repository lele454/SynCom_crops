#!/usr/bin/env python3
# [S1.5] Removal of out-of-scope crops: exclusion list
# Supplementary 1 of the manuscript. See README.md for the full map.
"""Retira estudos reprovados de runs_subset.tsv e converte exclusoes SRR->BioSample."""
import csv, shutil, os
RS = "runs_subset.tsv"
EX = os.path.expanduser("~/syncom/merged/curadoria_final/excluir_amostras.tsv")
FORA = {"PRJNA1147423", "PRJNA1028264", "PRJNA1184995", "PRJNA1514812"}

rows = list(csv.DictReader(open(RS), delimiter="\t"))
cols = rows[0].keys()
mapa = {r["run_accession"]: r["sample_accession"] for r in rows}
shutil.copy(RS, RS + ".bak2")

keep = [r for r in rows if r["study_accession"] not in FORA]
with open(RS, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(cols), delimiter="\t")
    w.writeheader(); w.writerows(keep)
print(f"runs_subset: {len(rows)} -> {len(keep)} corridas ({len(rows)-len(keep)} retiradas)")

ex = list(csv.DictReader(open(EX), delimiter="\t"))
shutil.copy(EX, EX + ".bak2")
out, conv, drop, falta = [], 0, 0, []
for r in ex:
    if r["study"] in FORA:
        drop += 1; continue
    run = r["run"]
    if run[:3] in ("SRR", "ERR", "DRR"):
        s = mapa.get(run)
        if not s:
            falta.append(run); continue
        out.append({"study": r["study"], "run": s, "origem": "auto_cultura",
                    "motivo": r.get("origem") or r.get("motivo") or ""})
        conv += 1
    else:
        out.append({"study": r["study"], "run": run,
                    "origem": r.get("origem", ""), "motivo": r.get("motivo", "")})
with open(EX, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["study", "run", "origem", "motivo"], delimiter="\t")
    w.writeheader(); w.writerows(out)
print(f"exclusoes: {len(ex)} -> {len(out)} | convertidas {conv} | removidas (estudo fora) {drop}")
if falta:
    print(f"ATENCAO sem mapeamento ({len(falta)}): {falta[:5]}")
