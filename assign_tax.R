#!/usr/bin/env Rscript
# [S1.9] Taxonomic assignment, naive Bayes classifier
# Supplementary 1 of the manuscript. See README.md for the full map.
# Taxonomia (classificador bayesiano do DADA2) para as sequencias unicas.
# Uso: Rscript assign_tax.R <dir merged> <referencia assignTaxonomy .fna.gz> <threads> [tamanho_bloco]
# Processa em blocos (padrao 5000 sequencias), retomavel (blocos prontos sao pulados).
# Saidas: taxonomy_v4.tsv e taxonomy_r1.tsv (com bootstrap de cada nivel)
suppressPackageStartupMessages(library(dada2))
a <- commandArgs(trailingOnly = TRUE)
md <- a[1]; ref <- a[2]; thr <- as.integer(a[3])
bloco <- if (length(a) >= 4) as.integer(a[4]) else 5000L
set.seed(100)
levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

run_set <- function(fasta, out) {
  if (!file.exists(fasta) || file.size(fasta) == 0) return(invisible())
  s <- getSequences(fasta)
  if (!length(s)) return(invisible())
  cdir <- file.path(md, paste0("tax_chunks_", sub("\\.tsv$", "", basename(out))))
  dir.create(cdir, showWarnings = FALSE)
  idx <- split(seq_along(s), ceiling(seq_along(s) / bloco))
  for (i in seq_along(idx)) {
    f <- file.path(cdir, sprintf("chunk_%04d.tsv", i))
    if (file.exists(f)) next
    t0 <- Sys.time()
    x <- s[idx[[i]]]
    gc()
    tx <- assignTaxonomy(unname(x), ref, multithread = thr, outputBootstraps = TRUE,
                         minBoot = 50, tryRC = FALSE, verbose = FALSE)
    tb <- tx$tax; bt <- tx$boot
    nl <- min(ncol(tb), length(levels))
    colnames(tb)[1:nl] <- levels[1:nl]; colnames(bt)[1:nl] <- paste0("boot_", levels[1:nl])
    d <- data.frame(seq_id = names(x), tb, bt, row.names = NULL, check.names = FALSE)
    write.table(d, paste0(f, ".part"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    file.rename(paste0(f, ".part"), f)
    cat(sprintf("%s: bloco %d/%d (%.1f min)\n", basename(out), i, length(idx),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  fs <- sort(list.files(cdir, "^chunk_\\d+\\.tsv$", full.names = TRUE))
  all <- do.call(rbind, lapply(fs, read.delim, check.names = FALSE, colClasses = "character"))
  # limpa prefixos do GTDB (d__, p__, g__...) se existirem
  for (l in intersect(levels, colnames(all))) all[[l]] <- sub("^[a-z]__", "", all[[l]])
  write.table(all, file.path(md, out), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  cat(sprintf("%s: %d sequencias | com genero: %.1f%%\n", out, nrow(all),
              100 * mean(!is.na(all$Genus) & all$Genus != "")))
}
run_set(file.path(md, "seqs_v4.fasta"), "taxonomy_v4.tsv")
run_set(file.path(md, "seqs_r1.fasta"), "taxonomy_r1.tsv")
