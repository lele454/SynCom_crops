#!/usr/bin/env bash
# [S3.6] Sequencing depth
# Supplementary 3 of the manuscript. See README.md for the full map.
# =====================================================================================
# Sensibilidade: profundidade minima de 5.000 leituras por amostra (padrao: 2.000).
# Refaz tabelas de genero, QC, efeitos e meta-analise em arquivos com sufixo _d5000,
# SEM tocar nos resultados principais. No fim, compara as duas analises.
#
# Retomavel: rode de novo o mesmo comando para continuar da etapa que faltou.
#   bash run_sens5000.sh
# Rodar de dentro da pasta dos scripts (~/Área de trabalho/SynComs), depois do run_all.sh.
# =====================================================================================
set -u
cd "$(dirname "$0")"
R=~/syncom
M=$R/merged
CUR=$M/curadoria_final
EST=$R/estado_sens5000; LOG=$R/logs_sens5000
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6
mkdir -p "$EST" "$LOG"
[ -e subset_studies.tsv ] || { echo "Falta subset_studies.tsv"; exit 1; }
[ -e $M/meta/meta_pairs.tsv.gz ] || { echo "Rode o run_all.sh antes"; exit 1; }

step() {
    local id=$1 desc=$2; shift 2
    if [ -e "$EST/$id.ok" ]; then echo "[ok]    $id $desc (ja feito)"; return; fi
    echo "[rodando] $id $desc  ($(date +%H:%M))"
    if "$@" > "$LOG/$id.log" 2>&1; then touch "$EST/$id.ok"; echo "[ok]    $id $desc  ($(date +%H:%M))"
    else echo "[ERRO]  $id $desc -- ultimas linhas de $LOG/$id.log:"; tail -15 "$LOG/$id.log"; exit 1; fi
}

step D1 "tabelas de genero com >= 5000 leituras" \
    python3 syncom_merge.py genus --root $R --studies subset_studies.tsv --min-depth 5000 --sufixo _d5000
step D2 "tamanhos de efeito" \
    python3 effect_sizes.py --root $R --runs runs_subset.tsv --attrs $CUR/attrs_curado \
        --design $CUR/design_curado.tsv --excluir $CUR/excluir_amostras.tsv \
        --sufixo _d5000 --outdir effects_d5000
step D3 "meta-analise" \
    python3 meta_pairs.py --root $R --effects effects_d5000 --tag d5000

python3 - <<'PY'
import os, numpy as np, pandas as pd
from scipy import stats
mm = os.path.expanduser("~/syncom/merged/meta")
key = lambda d: d.genus_a + "|" + d.genus_b
A = pd.read_csv(f"{mm}/meta_pairs.tsv.gz", sep="\t"); A.index = key(A)
B = pd.read_csv(f"{mm}/meta_pairs_d5000.tsv.gz", sep="\t"); B.index = key(B)
qc0 = pd.read_csv(os.path.expanduser("~/syncom/merged/study_qc.tsv"), sep="\t")
qc1 = pd.read_csv(os.path.expanduser("~/syncom/merged/study_qc_d5000.tsv"), sep="\t")
both = A.index.intersection(B.index)
rede = A[A.fdr < 0.05]
inter = rede.index.intersection(B.index)
ok = (B.loc[inter, "fdr"] < 0.05) & (np.sign(B.loc[inter, "mu"]) == np.sign(rede.loc[inter, "mu"]))
L = [f"Estudos incluidos: {(qc0.status=='incluido').sum()} (2000 leituras) vs "
     f"{(qc1.status=='incluido').sum()} (5000 leituras)",
     f"Pares testados: {len(A)} vs {len(B)} | em comum: {len(both)}",
     f"rho dos efeitos: {stats.spearmanr(A.loc[both,'mu'], B.loc[both,'mu']).correlation:.3f}",
     f"Rede de entrada (FDR<0,05) na analise principal: {len(rede)} | testaveis a 5000: {len(inter)} | "
     f"mantidos com mesmo sinal e FDR<0,05: {int(ok.sum())} ({ok.mean():.0%})"]
txt = "\n".join(L)
open(f"{mm}/sensibilidade_d5000.txt", "w").write(txt)
print("\n" + txt)
PY
echo; echo "FIM. Resumo: ~/syncom/merged/meta/sensibilidade_d5000.txt"
