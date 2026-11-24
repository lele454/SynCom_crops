#!/usr/bin/env Rscript
# [S1.8] DADA2 denoising of one study
# Supplementary 1 of the manuscript. See README.md for the full map.
# DADA2 de um estudo, na regiao nativa (V4 ou V3-V4).
# Uso: Rscript dada2_study.R <dir_fastq_aparados> <dir_saida> <regiao> <threads>
# Entrada: <run>_1.fastq.gz / <run>_2.fastq.gz (primers ja removidos)
# Saida:   asv_long.tsv (run, seq, count), track.tsv, params.tsv, seqtab_nochim.rds

suppressPackageStartupMessages({ library(dada2); library(ShortRead) })
a <- commandArgs(trailingOnly = TRUE)
indir <- a[1]; outdir <- a[2]; region <- a[3]; thr <- as.integer(a[4])
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(100)

fF <- sort(list.files(indir, "_1\\.fastq\\.gz$", full.names = TRUE))
fR <- sub("_1\\.fastq\\.gz$", "_2.fastq.gz", fF)
ok <- file.exists(fR) & file.size(fF) > 100 & file.size(fR) > 100
fF <- fF[ok]; fR <- fR[ok]
runs <- sub("_1\\.fastq\\.gz$", "", basename(fF))
if (length(fF) < 2) stop("menos de 2 corridas com leituras")

# ---- truncagem automatica pelo perfil de qualidade -------------------------
# comprimento do amplicon sem primers (E. coli): V4 ~253; V3-V4 ~400-430
amp <- if (region == "V4") 253 else 430

prof <- function(files, n = 3000) {
  q <- list(); w <- integer(0); L <- 0
  for (f in head(sample(files), 8)) {
    s <- FastqSampler(f, n = n); fq <- yield(s); close(s)
    if (length(fq) == 0) next
    w <- c(w, width(fq))
    m <- as(quality(fq), "matrix")
    q[[length(q) + 1]] <- colMeans(m, na.rm = TRUE)
    L <- max(L, ncol(m))
  }
  if (!length(q)) return(list(mean = numeric(0), len = 0, l10 = 0))
  mm <- sapply(q, function(v) { length(v) <- L; v })
  # l10: comprimento que 90% das leituras alcancam (leituras mais curtas que o
  # truncLen sao descartadas pelo filterAndTrim)
  list(mean = rowMeans(as.matrix(mm), na.rm = TRUE), len = L,
       l10 = as.integer(quantile(w, 0.10)))
}
qcap <- function(p, qmin) {
  v <- stats::filter(p$mean, rep(1 / 10, 10), sides = 1)   # media movel 10 pb
  bad <- which(!is.na(v) & v < qmin)
  qpos <- if (length(bad)) max(bad[1] - 1, 50) else p$len
  min(qpos, p$l10)
}
pF <- prof(fF); pR <- prof(fR)
# alvo: sobreposicao de ~30 pb sobre o amplicon; truncar mais curto = menos erros
target <- amp + 30; minimo <- amp + 12
capF <- qcap(pF, 25); capR <- qcap(pR, 25)
if (capF + capR < minimo) { capF <- qcap(pF, 20); capR <- qcap(pR, 20) }
if (capF + capR < minimo) { capF <- pF$l10; capR <- pR$l10 }
mode <- "PE"
if (capF + capR < minimo) {
  # leituras curtas demais para o merge (ex.: 2x150 na V4): usa so o R1,
  # truncado num comprimento fixo; gera um fragmento do inicio do amplicon
  if (pF$l10 < 90) {
    writeLines(sprintf("ERRO sem_sobreposicao_e_R1_curto: leituras %d+%d < %d", pF$l10, pR$l10, minimo),
               file.path(outdir, "ERRO"))
    quit(status = 3)
  }
  mode <- "R1"
  tF <- max(min(qcap(pF, 20), pF$l10), min(pF$l10, 100)); tR <- 0
} else {
  tF <- min(capF, ceiling(target * 0.55))
  tR <- min(capR, target - tF)
  if (tF + tR < target) tF <- min(capF, target - tR)
}

# ---- filtragem (com tentativas: PE maxEE 2/2 -> PE 2/5 -> so R1) -------------
filtF <- file.path(outdir, "filt", paste0(runs, "_F.fq.gz"))
filtR <- file.path(outdir, "filt", paste0(runs, "_R.fq.gz"))
run_filter <- function(mode, ee) {
  unlink(file.path(outdir, "filt"), recursive = TRUE)
  if (mode == "PE") {
    filterAndTrim(fF, filtF, fR, filtR, truncLen = c(tF, tR), maxN = 0, maxEE = ee,
                  truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = thr)
  } else {
    filterAndTrim(fF, filtF, truncLen = tF, maxN = 0, maxEE = ee[1], truncQ = 2,
                  rm.phix = TRUE, compress = TRUE, multithread = thr)
  }
}
ret <- function(o) sum(o[, 2]) / max(sum(o[, 1]), 1)
ee_used <- "2"
if (mode == "PE") {
  out <- run_filter("PE", c(2, 2)); ee_used <- "2,2"
  if (ret(out) < 0.2) { out <- run_filter("PE", c(2, 5)); ee_used <- "2,5" }
  if (ret(out) < 0.2) {                       # pares ruins: tenta so o R1
    mode <- "R1"; tR <- 0
    tF <- max(min(qcap(pF, 20), pF$l10), min(pF$l10, 100))
    out <- run_filter("R1", 2); ee_used <- "2"
  }
} else {
  out <- run_filter("R1", 2)
}
keep <- out[, 2] >= 500                    # corridas com >=500 leituras filtradas
if (sum(keep) < 2) stop("menos de 2 corridas apos filtragem")
filtF <- filtF[keep]; filtR <- filtR[keep]; runsK <- runs[keep]

# ---- erros, denoising, merge, quimeras ---------------------------------------
errF <- learnErrors(filtF, nbases = 1e8, multithread = thr, randomize = TRUE, verbose = 0)
dF <- dada(filtF, err = errF, multithread = thr, verbose = 0)
if (length(runsK) == 1) dF <- list(dF)
if (mode == "PE") {
  errR <- learnErrors(filtR, nbases = 1e8, multithread = thr, randomize = TRUE, verbose = 0)
  dR <- dada(filtR, err = errR, multithread = thr, verbose = 0)
  if (length(runsK) == 1) dR <- list(dR)
  mg <- mergePairs(dF, filtF, dR, filtR, minOverlap = 12, maxMismatch = 1)
  names(mg) <- runsK
  n_mg <- sum(sapply(mg, function(x) sum(x$abundance)))
  n_dn <- sum(sapply(dF, function(x) sum(getUniques(x))))
  if (n_mg < 0.2 * n_dn) {
    # merge falhou (sobreposicao real menor que a prevista): usa so o R1 ja processado
    mode <- "R1"; tR <- 0
    names(dF) <- runsK
    st <- makeSequenceTable(dF)
  } else {
    st <- makeSequenceTable(mg)
  }
} else {
  names(dF) <- runsK
  st <- makeSequenceTable(dF)
}
stn <- removeBimeraDenovo(st, method = "consensus", multithread = thr, verbose = FALSE)
saveRDS(stn, file.path(outdir, "seqtab_nochim.rds"))

# ---- saidas --------------------------------------------------------------------
idx <- which(stn > 0, arr.ind = TRUE)
long <- data.frame(run = rownames(stn)[idx[, 1]], seq = colnames(stn)[idx[, 2]],
                   count = stn[idx])
write.table(long, file.path(outdir, "asv_long.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
gn <- function(x) sum(getUniques(x))
track <- data.frame(run = runs, input = out[, 1], filtered = out[, 2],
                    denoisedF = NA, merged = NA, nonchim = NA)
track$denoisedF[keep] <- sapply(dF, gn)
track$merged[keep] <- if (mode == "PE") sapply(mg, function(x) sum(x$abundance)) else NA
if (mode == "R1" && exists("dR")) track$denoisedR <- NA
track$nonchim[keep] <- rowSums(stn)[runsK]
write.table(track, file.path(outdir, "track.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(region = region, mode = mode, maxEE = ee_used, retencao_filtro = round(ret(out), 3),
                       truncF = tF, truncR = tR, lenF = pF$len, lenR = pR$len,
                       n_asv = ncol(stn), frac_chim = 1 - sum(stn) / sum(st)),
            file.path(outdir, "params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
unlink(file.path(outdir, "filt"), recursive = TRUE)
cat("OK\n")
