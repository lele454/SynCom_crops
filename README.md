# To what extent does the co-occurrence of bacteria in the rhizosphere repeat across different soils, crops and years? A meta-analysis of 123 reprocessed studies

Analysis code for the manuscript above: a meta-analysis of bacterial co-occurrence in
the rhizosphere, built from 123 reprocessed public 16S rRNA amplicon studies
(6,557 samples, seven crops) under a single pipeline.

## What is here

Every script carries, in its header, the identifier of the supplementary section that
describes the corresponding step, so a reviewer reading one section can find the exact
script that produced it. The table below is that correspondence in full.

The layout is deliberately flat: several scripts import helpers from one another by
module name, so moving them into subfolders would break those imports.

## Running the analysis

```bash
bash run_all.sh            # runs every step in order; resumable
bash run_all.sh --recomecar  # clears the markers and recomputes everything
```

Each completed step leaves a marker file, so repeating the same command continues from
where it stopped.

## Requirements

Python 3 and R. Python dependencies are listed in `requirements.txt`; the R steps use
DADA2 and its dependencies. On some Linux installations where Python comes from conda,
the system C++ runtime has to be preloaded for the scientific stack to import:

```bash
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6
```

## Data

No sequence data is redistributed here. Every study is public: the accession of each of
the 123 studies, with its crop assignment, processing parameters and linked publication,
is listed in the supplementary tables of the manuscript. `rhizo_survey.py` (S1.1)
reproduces the survey from the ENA and the remaining scripts run from what it returns.

## Map of scripts to supplementary sections

### Supplementary 1 — Study survey, selection and sequence reprocessing

| ID | Script | Step |
| --- | --- | --- |
| S1.1 | `rhizo_survey.py` | ENA survey of candidate studies |
| S1.2 | `rhizo_stage2.py` | Amplicon region per run and BioSample attributes |
| S1.3 | `crop_audit.py` | Audit of the assigned host crop |
| S1.4 | `fix_crops.py` | Removal of out-of-scope crops: reassignment |
| S1.5 | `fix_excl.py` | Removal of out-of-scope crops: exclusion list |
| S1.6 | `make_subset.py` | Stratified subsampling of runs |
| S1.7 | `syncom_pipe.py` | Per-study processing driver |
| S1.8 | `dada2_study.R` | DADA2 denoising of one study |
| S1.9 | `assign_tax.R` | Taxonomic assignment, naive Bayes classifier |
| S1.10 | `assign_tax_sintax.py` | Taxonomic assignment, SINTAX fallback |
| S1.11 | `syncom_merge.py` | Merging, quality control and aggregation to genus |
| S1.12 | `apply_curation.py` | Sample curation |
| S1.13 | `auto_design.py` | Recovery of design factors from metadata |
| S1.14 | `compare_curado.py` | Curated against uncurated analysis |
| S1.15 | `lineage.py` | Genus lineage table |
| S1.16 | `fetch_refs.py` | Publication linked to each deposit |

### Supplementary 2 — Effect sizes and meta-analysis

| ID | Script | Step |
| --- | --- | --- |
| S2.1 | `effect_sizes.py` | Per-study effect sizes: genus selection, partial correlation, Fisher z |
| S2.2 | `meta_pairs.py` | Random-effects pooling, both intervals, FDR and the null model |

### Supplementary 3 — Stability and nature of the associations

| ID | Script | Step |
| --- | --- | --- |
| S3.1 | `robust_core.py` | Confirmation tiers |
| S3.2 | `sensitivity.py` | Sensitivity subsets |
| S3.3 | `zero_diagnostic.py` | Zero diagnostic |
| S3.4 | `robustness_extra.py` | Two-stage pooling by country and sugarcane exclusion |
| S3.5 | `frag_sens.py` | Read type |
| S3.6 | `run_sens5000.sh` | Sequencing depth |
| S3.7 | `moderators.py` | Moderation by management |
| S3.8 | `rice_meta.py` | Rice-only network |
| S3.9 | `temporal_validation.py` | Temporal validation |
| S3.10 | `replication_test.py` | Entry criterion tested against replication |
| S3.11 | `phylum_bias.py` | Phylum representation and mediation of the deviation |
| S3.12 | `q6_tracos.py` | Trait richness as a predictor of the positive-partner rate |
| S3.13 | `listas_q6.py` | Contrast group matched on abundance and occupancy |
| S3.14 | `rk_strategy.py` | r/K strategy hypothesis |
| S3.15 | `spore_bias.py` | Endospore hypothesis |
| S3.16 | `split_bias.py` | GTDB taxonomic splitting |
| S3.17 | `inoculant_bias.py` | Representation of commercial inoculant genera |
| S3.18 | `tau2_pm_check.py` | Tau-squared re-estimated by Paule-Mandel with Hartung-Knapp |

### Supplementary 4 — Recurrent combinations and feasibility

| ID | Script | Step |
| --- | --- | --- |
| S4.1 | `consortia.py` | Combinations from cliques of confirmed positive associations |
| S4.2 | `star_consortia.py` | Star combinations |
| S4.3 | `complementarity.py` | Functional complementarity |
| S4.4 | `shortlist.py` | Shortlist and per-crop analysis |
| S4.5 | `candidatas_testaveis.py` | Feasibility funnel |
| S4.6 | `rhizobia_partners.py` | Partners of genera with a registered inoculant |
| S4.7 | `modulo_fig.py` | Largest mutually associated set, co-presence tie-break |
| S4.8 | `modulo_nucleo.py` | Core shared by the tied sets |
| S4.9 | `network_removal.py` | In silico removal |
| S4.10 | `community_structure.py` | Independent test of community structuring |
| S4.11 | `crop_coverage.py` | Sampling against rarity as the cause of absence in a crop |
| S4.12 | `rank_syncoms.py` | Composite ranking |
| S4.13 | `rank_aptos.py` | Ranking corrected for sequencing effort |

### Supplementary 5 — Functional traits and strain availability

| ID | Script | Step |
| --- | --- | --- |
| S5.1 | `fetch_genomes.py` | Representative genomes of GTDB species clusters |
| S5.2 | `pgp_traits.py` | Growth-promotion traits from KOfam HMM profiles |
| S5.3 | `bgc_traits.py` | Secondary-metabolite gene clusters and NIS siderophores |
| S5.4 | `isolates.py` | Strain availability and isolation source per genus |

### Supplementary 6 — Software, figures and tables

| ID | Script | Step |
| --- | --- | --- |
| S6.1 | `run_all.sh` | Resumable driver for the complete analysis |
| S6.2 | `make_figures.py` | Figures: shared palette and plotting helpers |
| S6.3 | `figures_new.py` | Figures: panel layout and crop labels |
| S6.4 | `figures_extra.py` | Figures: data loading and saving helpers |
| S6.5 | `fig_candidatas.py` | Figures of the candidate combinations |
| S6.6 | `fig_trios.py` | Figures of the trios |
| S6.7 | `supp_tables.py` | Supplementary tables |
| S6.8 | `numeros_texto.py` | Values quoted in the Results text |
| S6.9 | `figs_final.py` | Final body and supplementary figures |
| S6.10 | `build_tables.py` | Body and supplementary tables |

## Citation

If you use this code, please cite the manuscript: *To what extent does the co-occurrence of bacteria in the rhizosphere repeat across different soils, crops and years? A meta-analysis of 123 reprocessed studies*
