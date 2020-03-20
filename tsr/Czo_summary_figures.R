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

ggsave("tss_correlation_Czo_single.png", plot = p, device = "png", type = "cairo", height = 3, width = 3)

###
# annotate TSSs
exp <- annotate_features(exp, annotation_file = Cz.annot, data_type = "tss", feature_type = "gene")
exp2 <- annotate_features(exp2, annotation_file = Cz.annot, data_type = "tss", feature_type = "gene")
exp3 <- annotate_features(exp3, annotation_file = Cz.annot, data_type = "tss", feature_type = "gene")

tss_distribution <- genomic_distribution(exp, data_type = "tss", threshold = 5)

p <- plot_genomic_distribution(tss_distribution) +
  ggplot2::theme(text = element_text(size = 6))

ggsave("tss_genomic_distribution_singlerep.png", plot = p, device = "png", type = "cairo", height = .5, width = 4)

### making feature plots

features <- detect_features(exp, data_type = "tss", feature_type = "gene", threshold = 5)

p <- plot_detected_features(features, ncol = 3) +
  ggplot2::theme(text = element_text(size = 5))

ggsave("tss_feature_plot_Pdec_all_reps.png", plot = p, device = "png", type = "cairo", height = 2, width = 4)

## plotting dinucleotide frequencies
# this generates an error; requires the truncated set

frequencies <- dinucleotide_frequencies(exp, genome_assembly = Cz.assembly, threshold = 3)

p <- plot_dinucleotide_frequencies(frequencies, ncol = 3) +
  ggplot2::theme(text = element_text(size = 6))

ggsave("tss_dinucleotide_frequencies_1_rep.png", plot = p, device = "png", type = "cairo", height = 2, width = 5)

## plotting the 5'UTR lengths

max <- max_utr(exp, threshold = 5, max_upstream = 250, max_downstream=100, feature_type = "geneId")
p <- plot_max_utr(max, ncol = 3, upstream = 250, downstream=100)
p
ggsave("max_utr_pdec_3reps.png", plot = p, device = "png", type = "cairo", height = 2, width = 6)

### TSR analyses

exp <- annotate_features(exp, annotation_file = annotation, data_type = "tsr", feature_type = "gene")

exp <- count_normalization(exp, data_type = "tsr")
exp2 <- count_normalization(exp2, data_type = "tsr")

# TSR correlation heatmap and scatter plot
# this is currently throwing up an error
p <- plot_correlation(exp, data_type = "tsr") +
  ggplot2::theme_bw() +
  ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_correlation_pdec_3reps.png", plot = p, device = "png", type = "cairo", height = 9, width = 9)

### tss location heatmap

count_matrix <- tss_heatmap_matrix(exp, threshold = 3, anno_type = "geneId", upstream = 100, downstream = 100)

p <- plot_heatmap(count_matrix, ncol = 3) +
  ggplot2::theme(text = element_text(size = 6))

ggsave("tss_heatmap.png", plot = p, device = "png", type = "cairo", height = 2, width = 4)

### tsr distribution

tsr_distribution <- genomic_distribution(exp, data_type = "tsr")

p <- plot_genomic_distribution(tsr_distribution) +
  ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_genomic_distribution.png", plot = p, device = "png", type = "cairo", height = 1.5, width = 4)

### plotting selected TSR metrics

p <- plot_tsr_metric(exp, tsr_metrics = c("nTAGs", "nTSSs"), log2_transform = TRUE, ncol = 2) +
  ggplot2::theme(text = element_text(size = 6))

ggsave("tsr_metrics.png", plot = p, device = "png", type = "cairo", width = 4, height = 2)

#### Making some figures in yeast

load("ScSTRIPE-wf3.RData")

## setting up the yeast tsrexplorer object

#making the exp files with all three replicates
wt.tss.1 <- ScSTRIPE@tssCountData[[4]]
wt.tss.2 <- ScSTRIPE@tssCountData[[5]]
wt.tss.3 <- ScSTRIPE@tssCountData[[6]]
dia.tss.1 <- ScSTRIPE@tssCountData[[10]]
dia.tss.2 <- ScSTRIPE@tssCountData[[11]]
dia.tss.3 <- ScSTRIPE@tssCountData[[12]]

names(wt.tss.1) <- c("seq", "start", "strand", "score", "isreal")
names(wt.tss.2) <- c("seq", "start", "strand", "score", "isreal")
names(wt.tss.3) <- c("seq", "start", "strand", "score", "isreal")
names(dia.tss.1) <- c("seq", "start", "strand", "score", "isreal")
names(dia.tss.2) <- c("seq", "start", "strand", "score", "isreal")
names(dia.tss.3) <- c("seq", "start", "strand", "score", "isreal")

wt.1.gr <- makeGRangesFromDataFrame(wt.tss.1, ignore.strand = FALSE, seqnames.field = "seq", start.field="start", end.field = "start", strand.field = "strand", keep.extra.columns = TRUE)
wt.2.gr <- makeGRangesFromDataFrame(wt.tss.2, ignore.strand = FALSE, seqnames.field = "seq", start.field="start", end.field = "start", strand.field = "strand", keep.extra.columns = TRUE)
wt.3.gr <- makeGRangesFromDataFrame(wt.tss.3, ignore.strand = FALSE, seqnames.field = "seq", start.field="start", end.field = "start", strand.field = "strand", keep.extra.columns = TRUE)
dia.1.gr <- makeGRangesFromDataFrame(dia.tss.1, ignore.strand = FALSE, seqnames.field = "seq", start.field="start", end.field = "start", strand.field = "strand", keep.extra.columns = TRUE)
dia.2.gr <- makeGRangesFromDataFrame(dia.tss.2, ignore.strand = FALSE, seqnames.field = "seq", start.field="start", end.field = "start", strand.field = "strand", keep.extra.columns = TRUE)
dia.3.gr <- makeGRangesFromDataFrame(dia.tss.3, ignore.strand = FALSE, seqnames.field = "seq", start.field="start", end.field = "start", strand.field = "strand", keep.extra.columns = TRUE)

Sc_tss <- list(wt.1.gr, wt.2.gr, wt.3.gr, dia.1.gr, dia.2.gr, dia.3.gr)
names(Sc_tss) <- c("wt_r1", "wt_r2", "wt_r3", "dia_r1", "dia_r2", "dia_r3")

### tsrs
wt_1_df <- ScSTRIPE@tsrData[[4]]
wt_2_df <- ScSTRIPE@tsrData[[5]]
wt_3_df <- ScSTRIPE@tsrData[[6]]
dia100_1_df <- ScSTRIPE@tsrData[[10]]
dia100_2_df <- ScSTRIPE@tsrData[[11]]
dia100_3_df <- ScSTRIPE@tsrData[[12]]

wt_1_df <- wt_1_df[,-12]
wt_2_df <- wt_2_df[,-12]
wt_3_df <- wt_3_df[,-12]

dia100_1_df <- dia100_1_df[,-12]
dia100_2_df <- dia100_2_df[,-12]
dia100_3_df <- dia100_3_df[,-12]

#wt.1.df.t <- data.frame(wt_1_df$seq, wt_1_df$start, wt_1_df$end, wt_1_df$strand, wt_1_df$score)
#wt.2.df.t <- data.frame(wt_2_df$seq, wt_2_df$start, wt_2_df$end, wt_2_df$strand, wt_2_df$score)
#wt.3.df.t <- data.frame(wt_3_df$seq, wt_3_df$start, wt_3_df$end, wt_3_df$strand, wt_3_df$score)

#names(wt.1.df.t) <- c("seq", "start", "end", "strand", "score")
#names(wt.2.df.t) <- c("seq", "start", "end", "strand", "score")
#names(wt.3.df.t) <- c("seq", "start", "end", "strand", "score")

#dia100_1_df.t <- data.frame(dia100_1_df$seq, dia100_1_df$start, dia100_1_df$end, dia100_1_df$strand, dia100_1_df$score)
#dia100_2_df.t <- data.frame(dia100_2_df$seq, dia100_2_df$start, dia100_2_df$end, dia100_2_df$strand, dia100_2_df$score)
#dia100_3_df.t <- data.frame(dia100_3_df$seq, dia100_3_df$start, dia100_3_df$end, dia100_3_df$strand, dia100_3_df$score)

#names(dia100_1_df.t) <- c("seq", "start", "end", "strand", "score")
#names(dia100_2_df.t) <- c("seq", "start", "end", "strand", "score")
#names(dia100_3_df.t) <- c("seq", "start", "end", "strand", "score")

wt.1.tsr.gr <- makeGRangesFromDataFrame(wt_1_df, ignore.strand = FALSE, seqnames.field = "seq", start.field="start", end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
wt.2.tsr.gr <- makeGRangesFromDataFrame(wt_2_df, ignore.strand = FALSE, seqnames.field = "seq", start.field="start", end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
wt.3.tsr.gr <- makeGRangesFromDataFrame(wt_3_df, ignore.strand = FALSE, seqnames.field = "seq", start.field="start", end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)

dia.1.tsr.gr <- makeGRangesFromDataFrame(dia100_1_df, ignore.strand = FALSE, seqnames.field = "seq", start.field="start", end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
dia.2.tsr.gr <- makeGRangesFromDataFrame(dia100_2_df, ignore.strand = FALSE, seqnames.field = "seq", start.field="start", end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
dia.3.tsr.gr <- makeGRangesFromDataFrame(dia100_1_df, ignore.strand = FALSE, seqnames.field = "seq", start.field="start", end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)

Sc_tsr <- list(wt.1.tsr.gr, wt.2.tsr.gr , wt.3.tsr.gr , dia.1.tsr.gr, dia.2.tsr.gr, dia.3.tsr.gr)
names(Sc_tsr) <- c("wt_1", "wt_2", "wt_3", "dia_1", "dia_2", "dia_3")

Sc <- tsr_explorer(Sc_tss, Sc_tsr)
Sc <- count_normalization(Sc, data_type = "tss", threshold = 3, n_samples = 2)
Sc <- count_normalization(Sc, data_type = "tsr", threshold = 3, n_samples = 2)

### Loading yeast TSSs and TSRs from the system file

TSSs <- system.file("extdata", "yeast_TSSs.RDS", package = "tsrexplorer")
TSSs <- readRDS(TSSs)

TSRs <- system.file("extdata", "yeast_TSRs.RDS", package = "tsrexplorer")
TSRs <- readRDS(TSRs)

annotation <- system.file("extdata", "yeast_annotation.gtf", package="tsrexplorer")
assembly <- system.file("extdata", "yeast_assembly.fasta", package="tsrexplorer")
#create tsr object

Sc_exp <- tsr_explorer(TSSs, TSRs)

Sc_exp <- count_normalization(Sc_exp, data_type = "tss", threshold = 3, n_samples = 2)
Sc_exp <- count_normalization(Sc_exp, data_type = "tsr", threshold = 5, n_samples = 2)

## Differential TSRs

edger_model <- fit_edger_model(
  Sc,
  data_type = "tss", 
  samples = c(
    "dia_r1", 
    "dia_r2",
    "dia_r3",
    "wt_r1",
    "wt_r2",
    "wt_r3"
  ),
  groups = c(1, 1, 1, 2, 2, 2)
)

#Alternative promoter usage
test.AP <- distanceToNearest(wt.1.tsr.gr,dia.1.tsr.gr)
summary(as.data.frame(test.AP)$distance)

AP.index <- which(as.data.frame(test.AP)$distance>50)
AP.table <- as.data.frame(test.AP[AP.index,])

### stacked barplot of AP usage
ggplot(Sc_distance, aes(x=type, y=promoter_status, fill=promoter_status)) +
  geom_bar(stat="identity", width = 0.5) +
  scale_fill_viridis(discrete=TRUE, name="") +
  ylab("Count") + coord_flip() + scale_x_discrete()
ggsave(file="AP_usage_stacked_bar.png", width=4, height = 1)

#making a histogram (excluding values greater than 1k)
p <- ggplot(Sc_distance, aes(Distance_to_TSR))
p + geom_histogram(colours="black",fill="purple",binwidth=25) + scale_x_continuous(limits = c(-100,1000)) +
  geom_vline(xintercept = 50, linetype="dotted", 
             color = "black", size=0.5)
ggsave(file="AP_usage_histogram.png", width=4, height = 3)

## making a gbrowse-style plot

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Gviz", version = "release")