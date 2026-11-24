#!/usr/bin/env bash
# [S6.1] Resumable driver for the complete analysis
# Supplementary 6 of the manuscript. See README.md for the full map.
# =====================================================================================
# Rodada final: analise curada como principal + rede de entrada por FDR < 0,05.
#
# 1. Arquiva os resultados antigos em ~/syncom/arquivo/v1_<data>/ (nada e apagado)
# 2. A analise curada (effects_curado) passa a ser a principal (merged/effects)
# 3. Refaz tudo, da meta-analise as figuras, nos locais de sempre
#
# Retomavel: cada etapa concluida fica marcada; se algo falhar ou o notebook desligar,
# rode de novo o mesmo comando e ele continua da etapa que faltou.
#   bash run_all.sh               roda (ou continua)
#   bash run_all.sh --recomecar   refaz todas as etapas (o arquivamento nao se repete)
#
# Rodar de dentro da pasta dos scripts (~/Área de trabalho/SynComs).
# =====================================================================================
set -u
cd "$(dirname "$0")"
R=~/syncom
M=$R/merged
CUR=$M/curadoria_final
RUNS=runs_subset.tsv
EST=$R/run_all_estado
LOG=$R/logs_run_all
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6
export PATH=~/.local/bin:$PATH
export SYNCOM_ENTRADA=fdr
mkdir -p "$EST" "$LOG"
[ "${1:-}" = "--recomecar" ] && find "$EST" -name 'S[0-9]*.ok' ! -name 'S00*' -delete
if [ "${1:-}" = "--desde" ]; then          # refaz a partir de uma etapa: bash run_all.sh --desde S07
    for f in "$EST"/S*.ok; do b=$(basename "$f" .ok); [[ "$b" > "${2}" || "$b" == "${2}" ]] && rm -f "$f"; done
fi

for f in "$RUNS" "$CUR/attrs_curado" "$CUR/design_curado.tsv" "$CUR/excluir_amostras.tsv"; do
    [ -e "$f" ] || { echo "Falta $f (rode apply_curation.py antes)"; exit 1; }
done
CURARGS="--attrs $CUR/attrs_curado --design $CUR/design_curado.tsv --excluir $CUR/excluir_amostras.tsv"

step() {   # step <id> <descricao> <comando...>
    local id=$1 desc=$2; shift 2
    if [ -e "$EST/$id.ok" ]; then echo "[ok]    $id $desc (ja feito)"; return; fi
    echo "[rodando] $id $desc  ($(date +%H:%M))"
    if "$@" > "$LOG/$id.log" 2>&1; then
        touch "$EST/$id.ok"; echo "[ok]    $id $desc  ($(date +%H:%M))"
    else
        echo "[ERRO]  $id $desc -- ultimas linhas de $LOG/$id.log:"; tail -15 "$LOG/$id.log"; exit 1
    fi
}

# ---------------------------------------------------------------- 0. arquivar o antigo
archive() {
    local A=$R/arquivo/v1_$(date +%F)
    mkdir -p "$A/merged" "$A/genomas_pgp" "$A/genomas_bgc"
    [ -d "$M/meta" ] && mv "$M/meta" "$A/merged/"
    if [ -d "$M/effects_curado" ]; then
        [ -d "$M/effects" ] && mv "$M/effects" "$A/merged/effects_v1"
        mv "$M/effects_curado" "$M/effects"
    fi
    for d in "$M"/effects_perm* "$M/curadoria_auto"; do [ -e "$d" ] && mv "$d" "$A/merged/"; done
    [ -d Figuras ] && mv Figuras "$A/Figuras"
    for f in ko_por_genoma.tsv tracos_por_genoma.tsv tracos_por_genero.tsv tracos_por_consorcio.tsv resumo.txt; do
        [ -e "$R/genomas/pgp/$f" ] && mv "$R/genomas/pgp/$f" "$A/genomas_pgp/"
    done
    for f in bgc_por_genoma.tsv bgc_por_genero.tsv bgc_por_consorcio.tsv resumo.txt; do
        [ -e "$R/genomas/bgc/$f" ] && mv "$R/genomas/bgc/$f" "$A/genomas_bgc/"
    done
    [ -e "$R/genomas/selecionados.tsv" ] && cp "$R/genomas/selecionados.tsv" "$A/"
    cat > "$A/LEIAME.txt" <<EOF
Resultados da versao 1 (arquivados em $(date +%F) por run_all.sh).
v1 = analise com deteccao automatica de fatores, sem exclusao de amostras fora da rizosfera,
     rede de entrada = pares "consistentes" (intervalo de predicao exclui zero).
A versao atual (em ~/syncom/merged) usa fatores curados, exclui amostras fora da rizosfera
e define a rede por FDR < 0,05. Os genomas baixados e as saidas brutas do GECCO foram mantidos.
EOF
    echo "Arquivado em $A"
}
step S00 "arquivar resultados antigos" archive
[ -d "$M/effects" ] || { echo "merged/effects nao existe"; exit 1; }

# ---------------------------------------------------------------- 1. meta-analise e robustez
step S01 "meta-analise principal"            python3 meta_pairs.py --root $R
step S02 "sensibilidade"                     python3 sensitivity.py --root $R
step S03 "nucleo robusto"                    python3 robust_core.py --root $R
step S04 "diagnostico de zeros"              python3 zero_diagnostic.py --root $R --runs $RUNS $CURARGS
# ---------------------------------------------------------------- 2. consorcios
step S05 "consorcios (cliques)"              python3 consortia.py --root $R
step S06 "consorcios em estrela"             python3 star_consortia.py --root $R --max-parceiros 10
step S06b "complementaridade funcional"      python3 complementarity.py --root $R
step S06c "dados do modulo (Fig6)"           python3 modulo_fig.py --root $R
step S07 "shortlist e culturas"              python3 shortlist.py --root $R
step S07b "candidatas testaveis in planta"   python3 candidatas_testaveis.py --root $R
step S08 "cobertura por cultura"             python3 crop_coverage.py --root $R
step S09 "arroz"                             python3 rice_meta.py --root $R --runs $RUNS $CURARGS
step S10 "remocao in silico"                 python3 network_removal.py --root $R
step S11 "estruturacao da comunidade"        python3 community_structure.py --root $R --runs $RUNS $CURARGS
# ---------------------------------------------------------------- 3. funcao (genomas)
step S12 "genomas de generos novos"          python3 fetch_genomes.py --root $R --release 232
step S13 "tracos PGP"                        python3 pgp_traits.py --root $R
step S14 "clusters biossinteticos (GECCO)"   python3 bgc_traits.py --root $R --por-genero 3
# ---------------------------------------------------------------- 4. validacao e ranking
step S15 "validacao temporal"                python3 temporal_validation.py --root $R
step S16 "teste de replicacao"               python3 replication_test.py --root $R
step S17 "ranking final"                     python3 rank_syncoms.py --root $R
# ---------------------------------------------------------------- 5. nulo e figuras
step S18 "modelo nulo (permutacao)"          python3 effect_sizes.py --root $R --runs $RUNS $CURARGS --permute 1 --outdir effects_perm1
step S19 "meta-analise do nulo"              python3 meta_pairs.py --root $R --effects effects_perm1 --tag nulo1
step S20 "figuras"                           python3 make_figures.py --root $R --studies subset_studies.tsv --world world_110m.geojson --out Figuras

# ---------------------------------------------------------------- resumo
{
  echo "=========== RESUMO DA RODADA FINAL ($(date +%F)) ==========="
  echo "--- meta-analise";        cat $M/meta/resumo.txt
  echo "--- nucleo robusto";      head -8 $LOG/S03.log
  echo "--- diagnostico zeros";   sed -n 1,30p $M/meta/diagnostico_zeros/resumo.txt
  echo "--- estrelas";            sed -n 1,20p $M/meta/consorcios/estrelas_resumo.txt
  echo "--- validacao temporal";  sed -n 1,20p $M/meta/validacao_temporal/resumo.txt
  echo "--- ranking";             cat $M/meta/consorcios/ranking_resumo.txt
} > $R/RESUMO_FINAL.txt
echo; echo "FIM. Resumo em $R/RESUMO_FINAL.txt | logs em $LOG"
