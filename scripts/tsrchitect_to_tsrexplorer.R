setwd("/scratch/rraborn/Czo_results/tsr")

library("TSRchitect")
library("farver")
library("tidyr", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.6")
library("tsrexplorer")
library("ggplot2")
library("viridis")
library("GenomicRanges")
library("dplyr")

load("PdSTRIPE_complete.RData")

CzSTRIPE <- PdSTRIPE #to keep the naming predicable/consistent

#creating the annotation and assembly files
Cz.annot <- "/scratch/rraborn/Czo_results/CzGENOME/Czo_gene_exons.gtf"
Cz.assembly <- "/scratch/rraborn/Czo_results/CzGENOME/Czofingiensis_461_v5.fa"

#writing the tss files to the workspace
tss.1 <- CzSTRIPE@tssCountData[[1]]
tss.2 <- CzSTRIPE@tssCountData[[2]]
tss.3 <- CzSTRIPE@tssCountData[[3]]
tss.4 <- CzSTRIPE@tssCountData[[4]]
tss.5 <- CzSTRIPE@tssCountData[[5]]
tss.6 <- CzSTRIPE@tssCountData[[6]]

#making granges files from tss data frames
tss.1.gr <- makeGRangesFromDataFrame(tss.1,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.2.gr <- makeGRangesFromDataFrame(tss.2,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.3.gr <- makeGRangesFromDataFrame(tss.3,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

#writing the tsr files to the workspace
tsr.1 <- CzSTRIPE@tsrData[[1]]
tsr.2 <- CzSTRIPE@tsrData[[2]]
tsr.3 <- CzSTRIPE@tsrData[[3]]
tsr.4 <- CzSTRIPE@tsrData[[4]]
tsr.5 <- CzSTRIPE@tsrData[[5]]
tsr.6 <- CzSTRIPE@tsrData[[6]]

#making granges files from tss data frames
tsr.1.gr <- makeGRangesFromDataFrame(tsr.1,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="start",
                                     end.field="end",
                                     strand.field = "strand"
)

tsr.2.gr <- makeGRangesFromDataFrame(tsr.2,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="start",
                                     end.field="end",
                                     strand.field = "strand"
)

tsr.3.gr <- makeGRangesFromDataFrame(tsr.3,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="start",
                                     end.field="end",
                                     strand.field = "strand"
)

tsr.4.gr <- makeGRangesFromDataFrame(tsr.4,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="start",
                                     end.field="end",
                                     strand.field = "strand"
)

tsr.5.gr <- makeGRangesFromDataFrame(tsr.5,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="start",
                                     end.field="end",
                                     strand.field = "strand"
)

tsr.6.gr <- makeGRangesFromDataFrame(tsr.6,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="start",
                                     end.field="end",
                                     strand.field = "strand"
)

#making the exp files with all three replicates
Cz.tss <- list(tss.1.gr, tss.2.gr, tss.3.gr)
names(Cz.tss) <- c("Cz_r1", "Cz_r2", "Cz_r3")

TSRs_cz <- list(tsr.1.gr, tsr.2.gr, tsr.3.gr)
names(TSRs_cz) <- c("Cz_r1,", "Cz_r2", "Cz_r3")

exp <- tsr_explorer(Cz.tss, TSRs_cz)
exp <- count_normalization(exp, data_type = "tss", threshold = 3, n_samples = 1)
exp <- count_normalization(exp, data_type = "tsr", threshold = 5, n_samples = 1)

#making the exp files with two replicates
Cz.tss <- list(tss.1.gr, tss.2.gr)
names(Cz.tss) <- c("Cz_r1", "Cz_r2")

TSRs_cz <- list(tsr.1.gr, tsr.2.gr)
names(TSRs_cz) <- c("Cz_r1,", "Cz_r2")

exp2 <- tsr_explorer(Cz.tss, TSRs_cz)
exp2 <- count_normalization(exp2, data_type = "tss", threshold = 3, n_samples = 1)

#making the exp files with a single replicate
Cz.tss.2 <- list(tss.2.gr)
names(Cz.tss.2) <- c("Cz_r2")

TSRs_cz.2 <- list(tsr.2.gr)
names(TSRs_cz.2) <- c("cz_r2")

exp3 <- tsr_explorer(Cz.tss.2, TSRs_cz.2)
exp3 <- count_normalization(exp3, data_type = "tss", threshold = 3, n_samples = 1)

###

p <- plot_correlation(exp, data_type = "tss") +
  ggplot2::theme_bw() +
  ggplot2::theme(text = element_text(size = 6))

p

ggsave("tss_correlation_Czo_single.png", plot = p, device = "png", type = "cairo", height = 6, width = 6)

###
# annotate TSSs
exp <- annotate_features(exp, annotation_file = Cz.annot, data_type = "tss", feature_type = "gene")
exp2 <- annotate_features(exp2, annotation_file = Cz.annot, data_type = "tss", feature_type = "gene")
exp3 <- annotate_features(exp3, annotation_file = Cz.annot, data_type = "tss", feature_type = "gene")

tss_distribution <- genomic_distribution(exp, data_type = "tss", threshold = 5)

p <- plot_genomic_distribution(tss_distribution) +
  ggplot2::theme(text = element_text(size = 6))

p

ggsave("tss_genomic_distribution_singlerep.png", plot = p, device = "png", type = "cairo", height = 2, width = 6)

### making feature plots

features <- detect_features(exp, data_type = "tss", feature_type = "gene", threshold = 5)

p <- plot_detected_features(features, ncol = 3) +
  ggplot2::theme(text = element_text(size = 5))

p

ggsave("tss_feature_plot_Pdec_all_reps.png", plot = p, device = "png", type = "cairo", height = 6, width = 10)
