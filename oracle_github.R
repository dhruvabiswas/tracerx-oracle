##  R code for the Nature Medicine publication Biswas D., Birkbak N.J., et al (2019) "A clonal expression biomarker associates with lung cancer mortality"

################
##  load packages
################

library(tidyverse)
library(gridExtra)
library(ggdendro)
library(ggrepel)
library(ggsignif)
library(cowplot)
library(readxl)
library(forcats)
library(reshape2)
library(survival)
library(survminer)
library(rmeta)
library(ReactomePA)
library(ComplexHeatmap)
library(dendsort)
library(circlize)
library(viridis)
library(viridisLite)

################
##  load data
################

##  FIGURE 1

# re-derived shukla signature
shukla_sig <- read.csv("shukla.JNCI.2017_rederived.sig.csv", row.names = 1, stringsAsFactors = F)

# tracerx FPKM matrix 
tracerx_fpkm <- read_rds("tracerx.fpkm.gene.matrix.rds")

##  FIGURE 2

# RNA heterogeneity quadrants for LUAD
ith_comparison_luad <- read_excel("Supp_Table_2.xlsx", "LUAD")

# collated prognostic signatures: 9 LUAD signatures, 280 total genes, 255 unique genes
prog <- read_excel("Supp_Table_3.xlsx")
# fix col-name
colnames(prog)[2] <- colnames(prog)[2] %>% tolower()
# join quadrant col
prog <- dplyr::left_join(x=prog, y=ith_comparison_luad[, c("Gene", "quadrant")], by="Gene")
# make quadrant col a factor
prog$quadrant <- factor(prog$quadrant, levels=c("top_left", "bottom_left", "top_right", "bottom_right"))
# fix signature names
prog$signature <- prog$signature %>% gsub(pattern=" ", replacement="\n")

# TCGA Cox UVA results: n=469 TCGA LUAD patients (stage I-III only)
tcga_res_luad <- read_rds("tcga.res.luad.RDS")

# microarray random signatures 
load("randomSigs.MA.venet.RData", verbose = T)

##  FIGURE 3

# published and clonal prognostic signature pipelines tested in Uppsala data
signature_pval <- read_csv("signature_perfomance_uppsala.csv")

# ORACLE Cox regression hazard ratios and P-values for Uppsala + microarray cohorts
oracle_meta <- read.csv("oracle_metaanalysis.csv", row.names = 1)

# uppsala ORACLE risk-scores
test_rs <- read_csv("ORACLE_riskscore_Uppsala.csv")
test_rs$bin <- ifelse(test_rs$RiskScore > 10.19941, "High", "Low")
test_rs$cohort <- "UPPSALA"
# uppsala clinical data
test_clinical <- read_csv("djureinovic.clinical.csv")   #V6

##  FIGURE 4

# load PRECOG values
precog_z <- read.delim("PRECOG-metaZ.pcl")

# load ith values for NSCLC
ith_comparison_nsclc <- read_excel("Supp_Table_2.xlsx", "NSCLC")

#  supplementary table 6: PRECOG HR and P-vals by RNA heterogeneity quadrant in NSCLC
precog <- read_excel("Supp_Table_6.xlsx", "PRECOG")

# rna-ith and dna-ith metrics
rna_dna_corr <- read.csv("rna_dna_corr.csv")

# scna by rna heterogeneity quadrant
quadrant_scna <- read.csv("quadrant_scna.csv")

##  EXTENDED DATA 1

#  load patient cohort composition data: tracerx, tcga, uppsala, der, okayama, rousseaux, shedden
load("patient_cohort_composition.Rdata", verbose = T)

##  EXTENDED DATA 2

# tracerx normalised counts, filtered for top 500 variably expressed genes
load("tracerx.counts_vsd_filter.Rdata", verbose = T)

# tracerx clinical information
tracerx.clinical <- read_rds("tracerx.clinical.rds")

# tracerx normalised counts, for all genes
tracerx.counts_vsd <- read_rds("tracerx.counts_vsd.rds")

##  EXTENDED DATA 3

# rna-ith metrics: standard deviation, median absolute deviation, coefficient of variation
intra_var <- read_rds("rnaith_scores.rds")

# inter-tumour rna heterogeneity metric: standard deviation in tcga and tracerx
inter_var <- read_rds("inter_var.rds")

##  EXTENDED DATA 4

# clustering concordance scores per gene (calculated in tracerx luad)
hclust_auc <- read_rds("hclust.auc_all.genes.rds")

##  EXTENDED DATA 5

# cox univariate regression results for published luad prognostic signatures in microarray cohorts
load("microarray_cox.Rdata", verbose = T)

##  EXTENDED DATA 6

# oracle risk-score threshold calibration in training cohort
riskscore_calibration <- read.csv("oracle_riskscore_calibration.csv")

# clustering concordance threshold calibration in training cohort
cv_error <- read.csv("oracle_cross_validation.csv")

# oracle risk-score in tracerx luad patients
tx_rs <- read.csv("tracerx_oracle_rs.csv", stringsAsFactors = F)
tx_rs$class <- factor(tx_rs$class, levels = c("Low", "Discordant", "High"))

##  EXTENDED DATA 7

##  load TNM-V8 data
test_clinical_tnm8 <- read.csv("test_clinical_tnm8.csv")

# ORACLE risk-scores in MET500 dataset
met_rs <- read.csv("ORACLE_riskscore_MET500.csv")

# tracerx nature supplementary data (Abbosh et al Nature 2017), includes ki67 staining percentage
tx.nature <- read_excel("nature22364-s2_table1.xlsx")
colnames(tx.nature)[1] <- "PublicationID"

##  EXTENDED DATA 8

# correlations between oracle risk-scores and danaher immune infiltration scores
oracle_immune_corr <- read.csv("oracle_immune_corr.csv", stringsAsFactors = F)
oracle_immune_corr$immune_subset <- factor(oracle_immune_corr$immune_subset, levels = oracle_immune_corr$immune_subset[order(oracle_immune_corr$rho, decreasing = T)])

# ascat purity scores and oracle risk-scores in tracerx luad patients
ascat <- read.csv("tracerx_ascat.csv", stringsAsFactors = F)

# Lambrechts et al (Nature Med 2018) supplementary table: single-cell RNAseq expression data for NSCLC stromal cell types summarised as 52 clusters
sc_expression <- read_excel("41591_2018_96_MOESM3_ESM.xls")
# make df
sc_expression <- sc_expression %>% as.data.frame()
# drop first row
sc_expression <- sc_expression[-1,]
# update colnames
colnames(sc_expression) <- colnames(sc_expression) %>% gsub(pattern=" ", replacement="")
# drop duplicate gene-names
which(duplicated(sc_expression[,1])) %>% length() # QC
sc_expression <- sc_expression[-which(duplicated(sc_expression[,1])),]
# set rownames as gene-name col
rownames(sc_expression) <- sc_expression[,1]
sc_expression <- sc_expression[,-1]
# order rows alphabetically
sc_expression <- sc_expression[order(rownames(sc_expression)), ]
# make values numeric 
sc_expression[] <- sapply(X=sc_expression, FUN=as.numeric)

#  ORACLE signature
oracle <- read_excel("Supp_Table_5.xlsx")
colnames(oracle)[1] <- "Gene"

# correlation between expression and copy-number state of oracle genes in tracerx luad patients
oracle_scna <- read.csv("oracle_scna.csv")

##  EXTENDED DATA 9

# no. of regions per patient
Mseq <- read.csv("Mseq.csv", stringsAsFactors = F)

# ith observed for down-sampled no. of regions
observed_ith <- read.csv("observed_ith.csv", stringsAsFactors = F)
# factor PatientID col by Nregions
observed_ith$PublicationID <- factor(observed_ith$PublicationID , levels = as.character(dplyr::arrange(Mseq, Nregions)$PublicationID))

# mean ith value observed for down-sampled no. of regions
output_mean <- read.csv("mseq_mean.csv", stringsAsFactors = F)
# factor PatientID col by Nregions
output_mean$PublicationID <- factor(output_mean$PublicationID , levels = as.character(dplyr::arrange(Mseq, Nregions)$PublicationID))

# standard deviation of ith values observed for down-sampled no. of regions
output_sd <- read.csv("mseq_sd.csv", stringsAsFactors = F)
# factor PatientID col by Nregions
output_sd$PublicationID <- factor(output_sd$PublicationID , levels = as.character(dplyr::arrange(Mseq, Nregions)$PublicationID))

# rna-ith scores and Danaher immune infiltration scores for tracerx NSCLC patients
rnaith_immune <- read.csv("rnaith_immune.csv", stringsAsFactors = F)

# rna-ith scores and ASCAT purity scores for tracerx NSCLC patients
rnaith_purity <- read.csv("rnaith_purity.csv")

################
##  load functions
################

##  core function for quadrant box-plots
boxplot_DB <- function(data, x, y, fill_var, fill_scale=c("firebrick1", "darkorchid1", "gold1", "turquoise1"), title="", subtitle="", x_axis="x", y_axis="y", test="wilcox", aspect_ratio=1, add_violin=T, exclude_NA=T, log_axis=T, zero_axis=T) {
  
  # tidy input data
  data <- data %>% as.data.frame()
  if(exclude_NA) {data <- data %>% drop_na()}
  
  # extract y-variable summary stats
  if(add_violin){ y_max <- ceiling(max(data[,y, drop=T], na.rm = T)) + 1 }
  if(!add_violin){ y_max <- ceiling(fivenum(data[,y, drop=T], na.rm = T)[4]) + 1 }
  
  # perform statistical test
  if(test=="wilcox") {
    
    ##  pairwise wilcox test
    # calculate p-values for all comparisons
    pv <- pairwise.wilcox.test(x = data[, y, drop=T], g = data[, x, drop=T], p.adjust.method="none", paired=F)$p.value
    # tidy: create "quad_1" and "quad_2" cols to show comparison for each test
    pv <- pv %>% as.data.frame() %>% rownames_to_column("quad_1")
    pv <- melt(pv, id.vars="quad_1")
    colnames(pv)[2:3] <- c("quad_2", "p_val")
    # QC: output to console
    print(pv)
    # filter significant values, and label with asterisks
    pv <- pv[which(pv$p_val < 0.05), ]
    pv$symbol <- ifelse(pv$p_val > 0.01, "*", ifelse(pv$p_val > 0.001, "**", "***"))
    # fix class of "quad_2" col
    pv$quad_2 <- pv$quad_2 %>% as.character()
    
  }
  
  # plot
  gg <- ggplot(data, aes_string(x, y, fill=fill_var)) 
  gg <- gg + scale_fill_manual(values=fill_scale)
  # add violin
  if(add_violin) {gg <- gg + geom_violin(alpha=0.25, scale="width", width=0.8) }
  # add boxplot
  gg <- gg + geom_boxplot(width=0.15, outlier.shape = NA)
  # add significance bars
  if(test=="wilcox") {
    if(nrow(pv)>0) {
      tmp <- as.data.frame(t(pv[, 1:2]))
      tmp[] <- sapply(tmp, as.character)
      gg <- gg + geom_signif(comparisons = as.list(tmp), annotation=pv$symbol)
    }
  }
  # add themes: black border, origin at 0, inward y tickmarks
  gg <- gg + theme_classic() + theme(legend.position = "none") + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = aspect_ratio, axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")), axis.text.x = element_text(size = 8, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),  axis.ticks.length = unit(-1.4, "mm")) 
  
  if(log_axis) {
    sig_values <- c(1, 0.05, 0.01, 10^(-seq(3, y_max)))
    gg <- gg + scale_y_continuous(expand = c(0,0), limits = c(0, y_max), breaks = sapply(X=sig_values, FUN=function(x) {-log10(x)}), labels=sig_values)
  }
  if(zero_axis) {
    gg <- gg + scale_y_continuous(expand = c(0,0), limits = c(0, y_max))
  }
  
  # add plot labels
  gg <- gg + ggtitle(label=title, subtitle=subtitle)
  gg <- gg + xlab(x_axis) + ylab(y_axis) + theme(plot.title=element_text(hjust=0.5, face="bold"))
  # output
  gg %>% return()
  
}

## wrapper function for gene-wise survival association ~ RNA heterogeneity quadrant as box-plots
quadrant_hr_boxplot_DB <- function(prog, gene_hazard_ratios, histology="LUAD", cohort_name="") {
  ##  UPDATE: 17 SEPTEMBER 2018, gene_hazard_ratios = output from rapid Cox UVA, unique gene names
  
  df <- prog
  
  ##  PLOT: p-value ~ signature, coloured by heterogeneity quadrant
  
  ##  prep HR info for joining
  tmp <- dplyr::select(gene_hazard_ratios, Gene, P)
  colnames(tmp)[which(colnames(tmp)=="P")] <- "pval"
  
  ##  prep prog: filter histology & unique gene names
  df <- prog %>% dplyr::filter(Histology==histology) %>% dplyr::select(Gene, quadrant) %>% dplyr::distinct()
  
  ##  inner join: drop NA values (quadrant, pval)
  df <- dplyr::inner_join(x=df, y=tmp, by="Gene")
  
  ##  no. signatures
  n_sigs <- nrow(unique(prog[prog$Histology==histology, "signature"]))
  
  ##prep data for plotting
  # log-transform pvals
  df$log_p <- -log10(df$pval)
  # rename quadrants
  df$q <- c("top_left"="Q1", "bottom_left"="Q2", "top_right"="Q3", "bottom_right"="Q4")[df$quadrant]
  
  ##  box-plot
  # boxplot
  gg <- boxplot_DB(
    data = df, 
    x = "q", 
    y = "log_p", 
    fill_var = "q", 
    title = "", 
    subtitle = paste(n_sigs, histology, "prognostic signatures"), 
    x_axis = "", 
    y_axis = paste0("Gene-Wise Prognostic Value in\n", cohort_name, " cohort (", gene_hazard_ratios$n_patients, " patients)\n-log10(Cox UVA P-value)"), zero_axis = F
  ) 
  gg + geom_hline(yintercept = -log10(0.05), lty="dashed") + geom_hline(yintercept = -log10(0.01), lty="dotted")
  
}

##  scatter plot and correlation function
corr_plot_DB <- function(data, x, y, point_colour="deepskyblue", colour_var=NA, point_border_col = "deepskyblue", title="", x_axis="x", y_axis="y", corr_method="spearman", pch_setting=21, point_size=5, point_alpha=0.5, x_limits = NA, y_limits = NA, aspect_ratio=1, best_fit_line=FALSE, line_colour="black", line_type="dashed") {
  #import
  corr <- data
  # spearman correlation
  if(corr_method=="spearman") {
    #calculate spearman's rho
    gg_cor <- cor.test(as.numeric(corr[,x]), as.numeric(corr[,y]), method = "spearman")$estimate %>% as.numeric() %>% signif(digits=3)
    #calculate spearman's p-value  
    gg_cor_p <- cor.test(corr[,x], corr[,y], method = "spearman")$p.value %>% as.numeric() %>% signif(digits=3)
  }
  # spearman correlation
  if(corr_method=="pearson") {
    #calculate PMCC
    gg_cor <- cor.test(as.numeric(corr[,x]), as.numeric(corr[,y]), method = "pearson")$estimate %>% as.numeric() %>% signif(digits=3)
    #calculate pearson's p-value  
    gg_cor_p <- cor.test(corr[,x], corr[,y], method = "pearson")$p.value %>% as.numeric() %>% signif(digits=3)
  }
  # fix P-value if v small
  if(gg_cor_p == 0) {gg_cor_p <- "< 2.2e-16"}
  # draw plot
  gg <- ggplot(corr, aes_string(x, y)) 
  # add points
  if(is.na(colour_var)) { 
    gg <- gg + geom_point(pch=pch_setting, col=point_border_col, fill=point_colour, size=point_size, alpha=point_alpha) 
  }
  if(!is.na(colour_var)) { 
    gg <- gg + geom_point(aes_string(col=colour_var), size=point_size, alpha=point_alpha) 
  }
  
  # add themes
  gg <- gg + theme_classic() 
  gg <- gg + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) #border
  gg <- gg + theme(aspect.ratio = aspect_ratio)
  # add spearman/pearson correlation outputs
  if(corr_method=="spearman") {
    #gg <- gg + ggtitle(label=title, subtitle=paste("Spearman's Rho = ", gg_cor, "\nSpearman's P-value = ", gg_cor_p, sep=""))
    subtitle <- paste0("Rs = ", gg_cor, "\nP = ", gg_cor_p)
  }
  if(corr_method=="pearson") { 
    #gg <- gg + ggtitle(label=title, subtitle=paste("PMCC = ", gg_cor, "\nPearson's P-value = ", gg_cor_p, sep=""))  
    subtitle <- paste0("PMCC = ", gg_cor, ", P = ", gg_cor_p)
  }
  
  # add title and subtitle
  gg <- gg + ggtitle(label=title, subtitle=subtitle)  
  
  # add best-fit line
  if(best_fit_line) {
    gg <- gg + geom_smooth(method = "lm", se = F, colour=line_colour, linetype=line_type)
  }
  
  gg <- gg + xlab(x_axis) + ylab(y_axis) + theme(plot.title=element_text(hjust=0.5, face="bold"))
  
  if(!is.na(x_limits)) {gg <- gg + xlim(c(x_limits[1], x_limits[2]))}
  if(!is.na(y_limits)) {gg <- gg + ylim(c(y_limits[1], y_limits[2]))}
  
  print(gg)
}

##  function to plot the composition of a cohort of patients as a bar-plot (# patients ~ stage, coloured by therapy status)
surv_object_barplot_DB <- function(surv_object, title) {
  library(plyr)
  library(tidyverse)
  
  ## calculate n-numbers
  #Histology
  n_histology <- data.frame(table(surv_object$Histology), stringsAsFactors = F)
  colnames(n_histology) <- c("Histology", "n_histology")
  n_histology$key <- paste0(n_histology$Histology, " (n=", n_histology$n_histology, ")")
  surv_object$Histology <- n_histology$key[match(surv_object$Histology, n_histology$Histology)]
  #therapy
  n_therapy <- data.frame(table(surv_object$Adjuvant.therapy), stringsAsFactors = F)
  colnames(n_therapy) <- c("Therapy", "n_therapy")
  n_therapy$key <- paste0(n_therapy$Therapy, " (n=", n_therapy$n_therapy, ")")
  surv_object$Adjuvant.therapy <- n_therapy$key[match(surv_object$Adjuvant.therapy, n_therapy$Therapy)]
  
  ## extract plotting parameters
  y_max <- round_any(max(table(surv_object$Stage, surv_object$Histology)), accuracy = 10, f = ceiling)
  y_breaks <- seq(0, y_max, 10) %>% fivenum()
  
  ##  make therapy col a factor
  surv_object$Adjuvant.therapy <- factor(surv_object$Adjuvant.therapy, levels = names(sort(table(surv_object$Adjuvant.therapy), decreasing = T)))
  
  ## bar-plot
  gg <- ggplot(surv_object, aes(x=Stage)) + geom_bar(aes(fill=Adjuvant.therapy)) + facet_wrap(~Histology, nrow=1, scales="free") 
  # labels
  gg <- gg + ggtitle(label = title)  + xlab("Tumour Stage") + ylab("No. of Patients")
  # themes
  gg <- gg + theme_classic() + theme(aspect.ratio=0.75) + labs(fill=NULL) + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + scale_y_continuous(expand = c(0,0), breaks = y_breaks, limits = c(0, y_max)) + scale_fill_manual(values= c("#8DA0CB", "indianred1", "wheat2", "gray75", "#FC8D62", "#E78AC3"))
  gg %>% print()
  
}

## function to calculate and plot gene-wise clustering concordance scores, quantifying the number of patients with all regions in the same cluster against the number of clusters (2 - # of patients)
gene_cluster_DB <- function(
  # specify genes of interest
  signature = oracle$Gene,
  # specify samples of interest
  sample_idx = luad_idx,
  # specify clustering method
  z_score = TRUE,
  hclust_method= "ward.D2",
  dist_method = "manhattan",
  # draw plots
  draw_plots = TRUE
) {
  
  ##   mat1
  # extract signature genes
  mat1 <- tracerx.counts_vsd[which(rownames(tracerx.counts_vsd) %in% signature), ] %>% t() %>% as.data.frame() 
  # extract Histology-specific and stage-specific patients
  mat1 <- mat1[which(rownames(mat1) %in% sample_idx), ]
  # calculate expression Z-scores
  if(z_score) {
    tmp <- sapply(X=mat1, FUN=scale)
    tmp <- as.data.frame(tmp)
    rownames(tmp) <- rownames(mat1)
    mat1 <- tmp
    #sort rows
    mat1 <- mat1[sort(rownames(mat1)), ]  
  }
  
  # initialize data-frame to store AUC per gene
  store_auc <- matrix(nrow=length(signature), ncol=4, dimnames = list(NULL, c("Gene", "AUC", "AUC_rel", "n_patient_concordance"))) %>% as.data.frame()
  
  ##  loop through one gene at a time
  for(i in 1:length(signature)) {
    
    ##  cut_dend data-frame
    # store gene identity
    gene_idx <- colnames(mat1)[i]
    # perform hierarchical clustering: Ward method on Manhattan metric
    mat_cluster_cols <- mat1 %>% dplyr::select(gene_idx) %>% dist(method=dist_method) %>% hclust(method = hclust_method)
    # cut dendrogram into k groups (no. of clusters = 1:total no. samples)
    cut_dend <- cutree(mat_cluster_cols, k=1:nrow(mat1))
    # fix colnames
    colnames(cut_dend) <- paste0("nclust_", colnames(cut_dend))
    # add PatientID col 
    cut_dend <- cut_dend %>% as.data.frame() %>% rownames_to_column("SampleID")
    cut_dend$PatientID <- sapply(X=cut_dend$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split=":"))[1]})
    
    ##  signature_clusters data-frame
    # identify how many clusters a subject's samples are spread across
    signature_clusters <- lapply(X=cut_dend[, -c(1,ncol(cut_dend))], FUN=function(x) aggregate(x=x, by=list(cut_dend$PatientID), FUN=function(x) {length(unique(x))})) %>% as.data.frame()
    # tidy data-frame: drop extra PatientID cols, rename clustering cols
    colnames(signature_clusters)[1] <- "PatientID"
    signature_clusters <- signature_clusters[, -grep(x=colnames(signature_clusters), pattern="Group.1")]
    colnames(signature_clusters) <- gsub(x=colnames(signature_clusters), pattern=".x", replacement="")
    # binarize data-frame: all samples in same cluster (1) or not (0)
    signature_clusters <- data.frame(PatientID=signature_clusters[,1], sapply(X=signature_clusters[,-1], FUN=function(x) {ifelse(x<=1, 1, 0)}))
    # calculate proportion of subjects with all samples in the same cluster
    tmp <- colSums(signature_clusters[,-1])/nrow(signature_clusters)
    signature_clusters <- data.frame(no_clusters=as.numeric(gsub(x=names(tmp), pattern="nclust_", replacement="")), prop_same_cluster=tmp)
    
    ## calculate AUC
    auc <- sum(diff(unique(signature_clusters$no_clusters)) * (head(signature_clusters$prop_same_cluster,-1)+tail(signature_clusters$prop_same_cluster,-1)))/2
    
    ## calculate prop_same_cluster at no_clusters == n_patients
    n_patient_concordance <- subset(signature_clusters, no_clusters==length(unique(cut_dend$PatientID)))$prop_same_cluster
    
    ## draw plots
    if (draw_plots) {
      
      ##  dendrogram
      # plot dendrogram
      gg <- ggdendrogram(mat_cluster_cols, rotate = FALSE, size = 2) + labs(title=gene_idx)
      # add themes
      gg <- gg + theme(aspect.ratio = 1, axis.line.y = element_line(), axis.text.y = element_text(hjust=0.5, size=8), axis.text.x = element_blank()) 
      gg <- gg + ggtitle(label=gene_idx, subtitle = paste0(hclust_method, ", ", dist_method)) + theme(plot.title=element_text(hjust=0.5, face="bold"))
      # print
      p1 <- gg
      
      ##  line-plot
      gg  <- ggplot(signature_clusters, aes(x=no_clusters, y=prop_same_cluster)) 
      # add data line
      gg <- gg + geom_line(col=plot_colour) 
      # add themes
      gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + theme(aspect.ratio = 1, panel.grid.minor = element_blank()) 
      # add plot labels
      gg <- gg + xlab("Number of clusters") + ylab("Proportion of patients with\nall samples in the same cluster") + ggtitle(label=gene_idx, subtitle = paste0("Absolute AUC = ", round(auc, digits=1), "\nRelative AUC = ", round(auc_rel, digits=1))) + theme(plot.title=element_text(hjust=0.5, face="bold"))
      # add vertical lines
      gg <- gg + scale_x_continuous(breaks = c(2, round(length(unique(cut_dend$PatientID))/2), length(unique(cut_dend$PatientID))), limits = c(0, ceiling(length(unique(cut_dend$PatientID))/10)*10)) + theme(panel.grid.major.x = element_line(colour="gray50", linetype="dashed"))
      #gg <- gg + scale_x_continuous(limits = c(0, ceiling(length(unique(cut_dend$PatientID))/10)*10)) + theme(panel.grid.major.x = element_line(colour="gray50"))
      # print
      p2 <- gg
      
      ## print plots
      cowplot::plot_grid(plotlist = list(p1, p2), nrow = 1, rel_widths = c(1,1)) %>% print()
    }
    
    ## store hclust metrics for gene
    store_auc[i,] <- c(gene_idx, auc, n_patient_concordance)
    # output progress to console 
    print(i)
  }
  
  ## output hclust metrics for each gene in signature
  return(store_auc)  
}

################
##  Figure 1
################

##  1c)
{
  ##  input: title
  title <- "Shukla (2017) JNCI signature\nin TRACERx LUAD cohort"
  
  ## input: Shukla signature 
  prognostic_signature <- shukla_sig
  #define riskscore_thresh
  riskscore_thresh <- prognostic_signature$RiskScore %>% unique()
  
  ## input: tracerx FPKM matrix
  gene_matrix <- tracerx_fpkm
  
  ##  select signature genes
  tmp <- gene_matrix[,which(colnames(gene_matrix) %in% prognostic_signature$Gene)] %>% t() %>% as.data.frame()
  colnames(tmp) <- gene_matrix$SampleID
  gene_matrix <- tmp
  
  ##  calculate RiskScore
  RiskScore <- colSums(gene_matrix[which(rownames(gene_matrix) %in% prognostic_signature$Gene), ] * prognostic_signature$beta)
  RiskScore <- data.frame(SampleID=names(RiskScore), RiskScore=RiskScore)
  RiskScore$PatientID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split=":"))[1]})
  RiskScore$RegionID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split=":"))[2]})
  # classify patients as high-risk or low-risk
  RiskScore$RiskScore_bin <- ifelse(RiskScore$RiskScore > riskscore_thresh, "High", "Low")
  
  ## get risk classes
  risk_class <- table(RiskScore$PatientID, RiskScore$RiskScore_bin)
  risk_class <- data.frame(High=as.matrix(risk_class)[,"High"], Low=as.matrix(risk_class)[,"Low"]) %>% rownames_to_column("PatientID")
  # assign as "low", "discordant" or "high"
  risk_class$class <- NA
  risk_class$class <- ifelse(risk_class$High > 0, paste(risk_class$class, "High", sep=""), risk_class$class)
  risk_class$class <- ifelse(risk_class$Low > 0, paste(risk_class$class, "Low", sep=""), risk_class$class)
  risk_class$class <- gsub(x=risk_class$class, pattern="NA", replacement="")
  risk_class$class <- gsub(x=risk_class$class, pattern="HighLow", replacement="Discordant")
  # join to RiskScore df
  RiskScore <- dplyr::left_join(x=RiskScore, y=risk_class[,c("PatientID", "class")], by="PatientID")
  RiskScore$class <- factor(RiskScore$class, levels = c("Low", "Discordant", "High"))
  # no. per class
  tmp <- RiskScore %>% dplyr::select(PatientID, class) %>% distinct()
  
  ## scatter-plot
  gg <- ggplot(RiskScore, aes(x=fct_reorder(PatientID, RiskScore + as.numeric(class), .fun=min), y=RiskScore)) 
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")), axis.text.x = element_text(size = 8, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),  axis.ticks.length = unit(-1.4, "mm")) + scale_y_continuous(expand = c(0,0), limits=c(floor(min(RiskScore$RiskScore)), ceiling(max(RiskScore$RiskScore))))
  gg <- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5))
  gg <- gg + geom_hline(yintercept = riskscore_thresh, col="black", lty="dotted")
  gg <- gg + geom_line(col="black")
  gg <- gg + ggtitle(label="", subtitle = paste0(length(unique(RiskScore$PatientID)), " TRACERx LUAD patients = ", table(tmp$class)["Low"], " low + ", table(tmp$class)["High"], " high + ", table(tmp$class)["Discordant"], " discordant")) + theme(plot.title = element_text(hjust=0.5, face="bold")) + xlab("PatientID") + ylab("Shukla JNCI (2017) RiskScore")
  gg <- gg + theme(legend.position = "bottom") + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) #border
  gg <- gg + theme(aspect.ratio = 0.5)
  gg <- gg + geom_point(pch=16, aes(col=class), size=3, alpha=0.5) + scale_color_manual(values = c("#3B4992FF", "azure4", "#EE0000FF")) + theme(legend.position = "none")
  gg %>% print()
}

##  1d - left) 
{
  ##  calculate percentages for risk classes
  tmp <- RiskScore %>% dplyr::select(PatientID, class) %>% dplyr::distinct()
  tmp <- data.frame(table(tmp$class))
  tmp$class <- tmp$Var1 %>% as.character()
  tmp$class[c(1,3)] <- paste0("Concordant\n", tmp$class[c(1,3)], " Risk")
  tmp <- tmp[c(1,3,2),]
  tmp$class <- factor(tmp$class, levels=tmp$class)
  tmp$Perc <- round(tmp$Freq/sum(tmp$Freq)*100)
  
  ##  bar-plot
  gg <- ggplot(tmp, aes(x=class, y=Perc, fill=class)) + geom_bar(stat="identity")
  gg <- gg + scale_fill_manual(values = c("#3B4992FF", "#EE0000FF", "gray75"), guide = guide_legend(title=NULL))
  gg <- gg + geom_text(size=5, aes(y = (Perc+2.5), label = paste0(Perc, "%"))) 
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.y = element_text(margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")), axis.text.x = element_text(margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),  axis.ticks.length = unit(-1.4, "mm")) + scale_y_continuous(expand = c(0,0), limits=c(0,50)) + theme(legend.position = "none", aspect.ratio = 1)
  gg <- gg + xlab("") + ylab("Survival risk classification (%)")
  gg <- gg + ggtitle(label="Shukla JNCI (2017) signature")
  print(gg)
}

##  1d - right) 
{
  ##  input: title
  title <- "Li (2017) JAMA Onc. signature\nin TRACERx LUAD cohort"
  
  ## input: Li signature 
  li.jamaonc.2017 <- dplyr::filter(prog, signature=="Li\nJAMA\nonc\n2017")$Gene
  ##  DB EDIT REQUIRED: use RiskScore threshold defined in paper
  #RiskScore_median <- TRUE
  
  ## input: tracerx FPKM matrix
  gene_matrix <- tracerx_fpkm
  
  #  select signature genes
  tmp <- gene_matrix[,which(colnames(gene_matrix) %in% li.jamaonc.2017)] %>% t() %>% as.data.frame()
  colnames(tmp) <- gene_matrix$SampleID
  gene_matrix <- tmp
  
  # initialize data-frame to store IGRPI scores
  tmp <- paste(li.jamaonc.2017_IRGPI$IRG_1, li.jamaonc.2017_IRGPI$IRG_2, sep="_")
  IGRPI_score <- matrix(nrow=length(tmp), ncol=ncol(gene_matrix), dimnames = list(tmp, colnames(gene_matrix))) %>% as.data.frame() #edited line to debug (2 SEPT 2018)
  
  # pairwise comparison to generate a score for each IRGP
  for(sample_i in 1:ncol(gene_matrix)) {
    # cycle through tumour samples
    for(genepair_i in 1:nrow(li.jamaonc.2017_IRGPI)) {
      # extract expression value for immune related gene pair 
      irg_1 <- gene_matrix[li.jamaonc.2017_IRGPI$IRG_1[genepair_i], sample_i]
      irg_2 <- gene_matrix[li.jamaonc.2017_IRGPI$IRG_2[genepair_i], sample_i]
      # "An IRGP score of 1 was assigned if IRG 1 was less than IRG 2; otherwise the IRGP score was 0"
      IGRPI_score[genepair_i, sample_i] <- ifelse(irg_1 < irg_2, 1, 0)  
    }
  }
  
  ## calculate RiskScore
  # multiply binary scores by gene-pair coeff
  tmp <- sapply(X=IGRPI_score, FUN=function(x) {x * li.jamaonc.2017_IRGPI$beta}) %>% as.data.frame()
  rownames(tmp) <- rownames(IGRPI_score)
  IGRPI_score.coeff <- tmp
  # calculate RiskScore
  RiskScore <- colSums(IGRPI_score.coeff)
  
  # make RiskScore data-frame
  RiskScore <- data.frame(SampleID=names(RiskScore), RiskScore=RiskScore)
  # using median RiskScore
  riskscore_thresh <- median(RiskScore$RiskScore)  
  # classify patients as high-risk or low-risk
  RiskScore$RiskScore_bin <- ifelse(RiskScore$RiskScore > riskscore_thresh, "High", "Low")
  # add PatientID and RegionID cols
  RiskScore$PatientID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split=":"))[1]})
  RiskScore$RegionID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split=":"))[2]})
  
  # add class col: High, Low, Discordant
  risk_class <- table(RiskScore$PatientID, RiskScore$RiskScore_bin)
  risk_class <- data.frame(High=as.matrix(risk_class)[,"High"], Low=as.matrix(risk_class)[,"Low"]) %>% rownames_to_column("PatientID")
  # assign as "low", "discordant" or "high"
  risk_class$class <- NA
  risk_class$class <- ifelse(risk_class$High > 0, paste(risk_class$class, "High", sep=""), risk_class$class)
  risk_class$class <- ifelse(risk_class$Low > 0, paste(risk_class$class, "Low", sep=""), risk_class$class)
  risk_class$class <- gsub(x=risk_class$class, pattern="NA", replacement="")
  risk_class$class <- gsub(x=risk_class$class, pattern="HighLow", replacement="Discordant")
  # join to RiskScore df
  RiskScore <- dplyr::left_join(x=RiskScore, y=risk_class[,c("PatientID", "class")], by="PatientID")
  RiskScore$class <- factor(RiskScore$class, levels = c("Low", "Discordant", "High"))
  # no. per class
  tmp <- RiskScore %>% dplyr::select(PatientID, class) %>% distinct()
  
  ##  calculate percentages for risk classes
  tmp <- RiskScore %>% dplyr::select(PatientID, class) %>% dplyr::distinct()
  tmp <- data.frame(table(tmp$class))
  tmp$class <- tmp$Var1 %>% as.character()
  tmp$class[c(1,3)] <- paste0("Concordant\n", tmp$class[c(1,3)], " Risk")
  tmp <- tmp[c(1,3,2),]
  tmp$class <- factor(tmp$class, levels=tmp$class)
  tmp$Perc <- round(tmp$Freq/sum(tmp$Freq)*100)
  
  ##  bar-plot
  gg <- ggplot(tmp, aes(x=class, y=Perc, fill=class)) + geom_bar(stat="identity")
  gg <- gg + scale_fill_manual(values = c("#3B4992FF", "#EE0000FF", "gray75"), guide = guide_legend(title=NULL))
  gg <- gg + geom_text(size=5, aes(y = (Perc+2.5), label = paste0(Perc, "%"))) 
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.y = element_text(margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")), axis.text.x = element_text(margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),  axis.ticks.length = unit(-1.4, "mm")) + scale_y_continuous(expand = c(0,0), limits=c(0,50)) + theme(legend.position = "none", aspect.ratio = 1)
  gg <- gg + xlab("") + ylab("Survival risk classification (%)")
  gg <- gg + ggtitle(label="Li JAMA Oncology (2017) signature")
  print(gg)
  
}

################
##  Figure 2
################

##  2a)
{
  ##  input: rna-ith scores
  rnaith_dhruva <- ith_comparison_luad %>% as.data.frame()
  colnames(rnaith_dhruva)[2:4] <- paste("LUAD", colnames(rnaith_dhruva)[2:4], sep="_")
  colnames(rnaith_dhruva) <- colnames(rnaith_dhruva) %>% gsub(pattern="intratumour_heterogeneity_score", replacement="within.heterogeneity") %>% gsub(pattern="intertumour_heterogeneity_score", replacement="between.heterogeneity")
  
  ## Make new quadrant type figure - 2 density plots, one with four to scale circles
  a <- table(factor(rnaith_dhruva[,'LUAD_quadrant'],levels=c("top_left", "bottom_left", "top_right", "bottom_right")))
  b <- a/min(a)
  
  plot(c(100,100,400,400),c(400,100,400,100), pch=16, cex=b*2, ylim=c(0,500),xlim=c(0,500), col=c("firebrick1", "darkorchid2", "gold1", "turquoise2"),axes=F,ylab='',xlab='')
  abline(h=250,v=250,col=1,lty=2,lwd=4)
  text(c(100,100,400,400),c(450,200,450,200), paste(a,'genes'))
  box()
  
  hist(rnaith_dhruva[,'LUAD_between.heterogeneity'], breaks=100, prob=TRUE,col='#6baed699',axes=F,ylab='', xlab='',main='Between tumour, LUAD')
  lines(density(rnaith_dhruva[,'LUAD_between.heterogeneity']), lwd = 3,col = '#08519c')
  abline(v=mean(rnaith_dhruva[,'LUAD_between.heterogeneity']), col=1,lwd=4, lty=2)
  hist(rnaith_dhruva[,'LUAD_within.heterogeneity'], breaks=100, prob=TRUE,col='#fd8d3c99',axes=F,ylab='', xlab='',main='Within tumour, LUAD')
  lines(density(rnaith_dhruva[,'LUAD_within.heterogeneity']), lwd = 3,col = '#a63603')
  abline(v=mean(rnaith_dhruva[,'LUAD_within.heterogeneity']), col=1,lwd=4, lty=2)
  
}

##  2b)
{
  ##  input: prognostic signatures
  prog_histology <- prog
  
  ##  input: heterogeneity scores
  heterogeneity_scores <- ith_comparison_luad
  colnames(heterogeneity_scores) <- colnames(heterogeneity_scores) %>% gsub(pattern="intratumour_heterogeneity_score", replacement="within.tumour.ith") %>% gsub(pattern="intertumour_heterogeneity_score", replacement="between.tumour.ith")
  
  ##  Calculate percentages of expected & observed genes per heterogeneity quadrant
  ##  expected genes per quadrant
  # counts genes per quadrant in full expression data-frame
  q_expected <- data.frame(table(heterogeneity_scores$quadrant))
  colnames(q_expected) <- c("quadrant", "count")
  # convert to percent
  q_expected$percent <- q_expected$count / sum(q_expected$count) * 100
  # prep colnames for joining
  colnames(q_expected)[2:3] <- paste0("expected_", colnames(q_expected)[2:3])
  ##  observed genes per quadrant
  # count genes per quadrant
  q_observed <- data.frame(table(prog_histology$quadrant))
  colnames(q_observed) <- c("quadrant", "count")
  # convert to percent
  q_observed$percent <- q_observed$count / sum(q_observed$count) * 100
  # prep colnames for joining
  colnames(q_observed)[2:3] <- paste0("observed_", colnames(q_observed)[2:3])
  ## join
  q_comparison <- dplyr::full_join(q_expected, q_observed, by="quadrant")
  # order quadrant factor
  q_comparison$quadrant <- factor(q_comparison$quadrant, levels=c("top_right", "bottom_right", "bottom_left", "top_left"))
  # observed % / expected %
  q_comparison$OE <- q_comparison$observed_percent / q_comparison$expected_percent
  
  ##  Fisher's exact test
  # construct contingency matrix
  fisher_matrix <- matrix(c(
    # no. O genes in other 3 quadrants
    (sum(q_comparison$expected_count) - q_comparison[2, "expected_count"]),
    # no. O genes in Q4
    q_comparison[2, "expected_count"],
    # no. E genes in Q4
    (sum(q_comparison$observed_count) - q_comparison[2, "observed_count"]),
    # no. E genes in Q4
    q_comparison[2, "observed_count"]
  ),2,2)
  # perform Fisher's exact test
  fisher_p <- fisher.test(fisher_matrix)$p.value
  
  ##  coxcomb plot
  # order colours
  cols <- c("gold1", "turquoise2", "darkorchid2", "firebrick1")
  ##  make plot
  p <- ggplot(q_comparison, aes(x=quadrant)) 
  # Observed percent / Expected percent
  p <- p + geom_bar(aes(y=OE, fill=quadrant), stat="identity", width=1) 
  # overlay dashed black line: 1-1 observed-expected
  p <- p + geom_bar(aes(y=1), stat="identity", width=1, alpha=0, lty="dashed", colour="black") 
  p <- p + geom_bar(aes(y=ceiling(max(q_comparison$OE))), stat="identity", width=1, alpha=0,  colour="white") 
  # add data labels
  #p <- p + geom_label(aes(y=ceiling(max(q_comparison$OE)), label=paste0(round(q_comparison$OE, digits=1), "x"), colour=quadrant))
  p <- p + geom_text(aes(y=ceiling(max(q_comparison$OE))+0.5, label=paste0(round(q_comparison$OE, digits=1), "x")))
  # add quadrant borders
  p <- p + geom_vline(xintercept = c(0.5, 1.5, 2.5, 3.5), linetype="dashed")
  # make circular plot
  p <- p + coord_polar()
  # colours
  p <- p + scale_color_manual(values=cols) + scale_fill_manual(values=cols)
  # themes
  p <- p + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + theme(aspect.ratio=1) + theme(axis.text = element_blank(), axis.ticks = element_blank()) + theme(legend.position = "none")
  # plot labels
  p <- p + ggtitle(label = "O/E Prognostic Signature Genes ~ Heterogeneity Quadrant", subtitle = paste0(length(unique(prog_histology$signature)), " ", ui_histology, " signatures (n genes = ", sum(q_comparison$observed_count, na.rm = T), ")\nQ4 enrichment: Fisher's exact P = ", signif(fisher_p, digits=3))) + theme(plot.title = element_text(hjust=0.5, face="bold")) + xlab("Prognostic signature genes per heterogeneity quadrant\n(Observed % / Expected %)") + ylab("")
  p %>% print()
}

##  2c)
{
  # input
  prog_histology <- prog
  prog_histology$Histology <- "LUAD"
  
  #  run analysis
  gg <- quadrant_hr_boxplot_DB(prog_histology, tcga_res_luad, "LUAD", "TCGA")
  
  ##  re-plot with greater y-axis upper limit
  y_max <- 9
  sig_values <- c(1, 0.05, 0.01, 10^(-seq(3, y_max)))
  gg <- gg + scale_y_continuous(expand = c(0,0), limits = c(0, y_max), breaks = sapply(X=sig_values, FUN=function(x) {-log10(x)}), labels=sig_values)
  gg %>% print()
}

##  2d)
{
  # count number of significant random signatures
  a <- apply(output_MA[,2:5],1, function(x){sum(x < 0.05)})
  b <- apply(output_MA[,8:11],1, function(x){sum(x < 0.05)})
  aa <- apply(output_MA[,14:17],1, function(x){sum(x < 0.05)})
  ab <- apply(output_MA[,20:23],1, function(x){sum(x < 0.05)})
  
  # bar-plot
  tmp<-cbind(table(factor(a,levels=0:4))/10,table(b)/10,table(aa)/10,table(ab)/10)
  barplot(t(tmp),beside=T,col=oracle_cols,main='1000 random 20-gene signatures', ylab='Percent of rounds', xlab='Number of cohorts significant')
  legend('top',legend=c('Q1 genes', 'Q2 genes','Q3 genes', 'Q4 genes'),col=oracle_cols,pch=15,bty='n')
}

################
##  Figure 3
################

##  3a)
{
  # prep data for plotting
  output <- signature_pval
  sig_values <- c(1, 0.05, 0.01, 0.001)
  metric <- c("logrank_p", "ci_p", "uva_p")
  my_cols <- c(rep(x = c("darkorange1","turquoise3"), (length(unique(output$signature))-1)/2), "turquoise3")
  my_pch <- c(21, 24, 22)
  output$signature <- factor(output$signature, levels=c("p_stepwise_aic", "c_stepwise_aic", "p_stepwise_bic", "c_stepwise_bic", "p_tree", "c_tree", "p_forest", "c_forest", "p_lasso", "c_lasso", "ORACLE"))
  colnames(output)[2] <- "p"
  output$p_adj <- ifelse(output$p < min(sig_values), min(sig_values), output$p)
  
  # lollipop plot
  gg <- ggplot(output, aes(x=signature, y=-log10(p_adj), colour=signature, fill=signature))
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + theme(aspect.ratio = 0.5) + theme(legend.position = "none")
  gg <- gg + xlab("Prognostic Signature") + ylab(paste0("Cox UVA P")) 
  gg <- gg + geom_linerange(aes(ymin = 0, ymax = -log10(p_adj))) + scale_color_manual(values=rep(x = "black", length(output$signature)))
  gg <- gg + geom_hline(yintercept = -log10(0.05), linetype="dashed") + geom_hline(yintercept = -log10(0.01), linetype="dotted")
  gg <- gg + geom_point(position = position_dodge(width=0.5), pch=21, stroke=0, size=4) + scale_fill_manual(values=my_cols)
  gg <- gg + scale_y_continuous(breaks = sapply(X=sig_values, FUN=function(x) {-log10(x)}), labels=sig_values, limits = c(-log10(max(sig_values)), -log10(min(sig_values))), expand = c(0,0))
  gg <- gg + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
  gg %>% print()
}

##  3b)
{
  # input: Cox regression hazard ratios and P-values for Uppsala + microarray cohorts
  tmp <- oracle_meta
  
  # use log values 
  tmp <- tmp[,2:3]
  sum.estimate <- meta.summaries(tmp[-1,1],tmp[-1,2],logscale=TRUE)
  
  # meta-analysis plot
  metaplot(tmp[,1],tmp[,2], labels=rownames(tmp), xlab='Log Hazard Ratio', main='ORACLE', summn = sum.estimate$summary, sumse = sum.estimate$se.summary, sumnn= 1/sum.estimate$se.summary^2, xlim=c(-1, 3), zero=0, colors=meta.colors(box="#3182bd",lines="#a50f15", zero="red", summary="black",text="black"),xaxt='n')
  axis(1, at=log(c(0.5,1,2,4,8,16)), labels=c(0.5,1,2,4,8,16))
  
  # p-value, : https://www.bmj.com/content/343/bmj.d2304
  #If the upper and lower limits of a 95% CI are u and l respectively:
  #1 calculate the standard error: SE = (u − l)/(2×1.96)
  # 2 calculate the test statistic: z = Est/SE
  # 3 calculate the P value: P = exp(−0.717×z − 0.416×z2).
  Est <- 3.1#sum.estimate$summary
  l <- 2.3
  u <- 4.19
  SE <- (u-l)/(2*1.96)
  z = Est/SE
  P = exp(-0.717*z - 0.416*z^2)
  # P= 3.381466e-10
}

##  3c)
{
  ##  input: risk-scores in uppsala cohort
  cohort_riskscore <- test_rs
  clinical_data <- test_clinical
  
  #  make title
  tmp <- subset(clinical_data, PatientID %in% test_rs$PatientID)
  title_full <-  paste("ORACLE\nValidation Data-set", " (", nrow(test_rs), " patients, stage ", min(tmp$Stage_numeric, na.rm = T), "-", max(tmp$Stage_numeric, na.rm = T), ")", sep="")
  
  #stratify cohort
  cohort_riskscore$PatientID <- cohort_riskscore$PatientID %>% as.character()
  #join survival analysis
  cohort_riskscore <- left_join(x=cohort_riskscore, y=clinical_data, by="PatientID")
  # convert bin col to factor
  cohort_riskscore$bin <- factor(cohort_riskscore$bin, levels=c("Low", "High"))
  
  # re-code clinical parameters as factor
  # re-code Stage_numeric as factor
  cohort_riskscore$Stage_numeric <- factor(cohort_riskscore$Stage_numeric)
  # re-code Adjuvant.therapy as factor
  cohort_riskscore$therapy <- ifelse(is.na(cohort_riskscore$Adjuvant.therapy), "No adjuvant treatment", "Adjuvant treatment")
  cohort_riskscore$therapy <- factor(cohort_riskscore$therapy, levels=c("No adjuvant treatment", "Adjuvant treatment"))
  # re-code smoking_history as factor
  cohort_riskscore$smoking_history <- factor(cohort_riskscore$smoking_history, levels=unique(cohort_riskscore$smoking_history))
  # re-code gender as factor
  cohort_riskscore$gender <- factor(cohort_riskscore$gender, levels=unique(cohort_riskscore$gender))
  # re-code performance_status as factor
  cohort_riskscore$performance_status <- factor(cohort_riskscore$performance_status, levels=sort(unique(cohort_riskscore$performance_status)))
  
  # multivariate Cox regression: RiskScore, Stage, Adjuvant.therapy, age, performance status, smoking history, gender, KI67
  multivar_cox.res <- coxph(formula = Surv(cohort_riskscore$OS, cohort_riskscore$status==1) ~ cohort_riskscore$RiskScore + cohort_riskscore$Stage_numeric + cohort_riskscore$therapy + cohort_riskscore$age + cohort_riskscore$performance_status + cohort_riskscore$smoking_history + cohort_riskscore$gender + cohort_riskscore$KI67_percentage)
  
  # extract HR and confidence interval
  tmp <- summary(multivar_cox.res)
  tmp <- data.frame(HR=tmp$coefficients[,2], lower_ci=tmp$conf.int[,3], upper_ci=tmp$conf.int[,4], P=tmp$coefficients[,5]) %>% rownames_to_column("Predictor")
  # fix Predictor col
  tmp$Predictor <- sapply(tmp$Predictor, function(x) unlist(strsplit(x, split="\\$"))[2])
  # add reference rows: Stage_numeric1
  idx <- min(grep("Stage", tmp$Predictor)) - 1
  newRow <- c("Stage_numeric1", 1, NA, NA, NA)
  tmp <- rbind(tmp[1:idx, ], newRow, tmp[-(1:idx), ])
  # add reference rows: therapyNo adjuvant treatment
  idx <- min(grep("therapy", tmp$Predictor)) - 1
  newRow <- c("therapyNo adjuvant treatment", 1, NA, NA, NA)
  tmp <- rbind(tmp[1:idx, ], newRow, tmp[-(1:idx), ])
  # add reference rows: smoking_history
  idx <- min(grep("smoking_history", tmp$Predictor)) - 1
  newRow <- c("smoking_historynever", 1, NA, NA, NA)
  tmp <- rbind(tmp[1:idx, ], newRow, tmp[-(1:idx), ])
  # add reference rows: gender
  idx <- min(grep("gender", tmp$Predictor)) - 1
  newRow <- c("gendermale", 1, NA, NA, NA)
  tmp <- rbind(tmp[1:idx, ], newRow, tmp[-(1:idx), ])
  # add reference rows: performance_status
  idx <- min(grep("performance_status", tmp$Predictor)) - 1
  newRow <- c("performance_status0", 1, NA, NA, NA)
  tmp <- rbind(tmp[1:idx, ], newRow, tmp[-(1:idx), ])
  
  # convert to numeric, and 3 signif figures
  tmp[,-1] <- sapply(tmp[,-1], function(x) signif(as.numeric(x), digits=3))
  
  ## forest plot
  hr_min <- floor(min(tmp$lower_ci, na.rm = T)*100)/100
  hr_max <- ceiling(ceiling(max(tmp$upper_ci, na.rm = T))/10)*10
  ##  forest plot
  gg <- ggplot(tmp, aes(x=Predictor, y=HR)) + coord_flip()
  gg <- gg + geom_hline(yintercept = 1, linetype="dotted")
  gg <- gg + geom_errorbar(aes(x=Predictor, ymin=lower_ci, ymax=upper_ci), width=0.5) 
  gg <- gg + geom_point(pch=15, size=3)
  gg <- gg + theme_minimal() + theme(legend.position = "none") + theme(aspect.ratio = 0.75) + theme(panel.grid = element_blank()) + theme(axis.line.x = element_line(), axis.ticks.x = element_line())
  gg <- gg + scale_y_continuous(trans = "log", expand = c(0,0), breaks = c(0.5,1,2,4,8,16), limits = c(hr_min, hr_max)) + scale_x_discrete(limits=rev(tmp$Predictor)) 
  gg <- gg + xlab("Predictor") + ylab("Hazard Ratio") + ggtitle(title_full)
  gg <- gg + geom_text(aes(x=Predictor, y=hr_max-0.5, label=signif(P, digits=3)), size=3)
  gg %>% print()
}

##  3d)
{
  ##  input: filter stage 1 patients
  tmp <- dplyr::left_join(x=test_rs, y=test_clinical, by="PatientID")
  tmp <- tmp %>% dplyr::filter(Stage_numeric==1)
  cohort_riskscore <- tmp
  
  ## sub-stage test
  # log-rank P
  sdf <- survdiff(Surv(cohort_riskscore$OS, cohort_riskscore$status==1) ~ cohort_riskscore$Stage, data = cohort_riskscore)
  log.rank_p_stage <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1) %>% signif(digits=3)
  # KM plot: substage
  ggsurvplot(survfit(Surv(cohort_riskscore$OS, cohort_riskscore$status==1) ~ cohort_riskscore$Stage, data = cohort_riskscore), risk.table = TRUE, tables.theme = theme_cleantable(), palette = "aaas", title=paste("Substaging\nlog-rank p = ", log.rank_p_stage, sep=""), legend="none", xlab="Time (Years)", ylab="Overall Survival Probability") %>% print()
  
  ##  oracle risk-score test
  # log-rank P 
  sdf <- survdiff(Surv(cohort_riskscore$OS, cohort_riskscore$status==1) ~ cohort_riskscore$bin, data = cohort_riskscore)
  log.rank_p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1) %>% signif(digits=3)
  # KM plot: ORACLE
  ggsurvplot(survfit(Surv(cohort_riskscore$OS, cohort_riskscore$status==1) ~ cohort_riskscore$bin, data = cohort_riskscore), risk.table = TRUE, legend.labs = c("High", "Low"), tables.theme = theme_cleantable(), palette = "aaas", title=paste("ORACLE\nlog-rank p = ", log.rank_p, sep=""), legend="none", xlab="Time (Years)", ylab="Overall Survival Probability") %>% print()
}


################
##  Figure 4
################

## 4a)
{
  # prep pan-cancer z-score values from PRECOG 
  data <- dplyr::inner_join(x=precog_z[, c("Gene", "Unweighted_meta.Z_of_all_cancers")], y=ith_comparison_nsclc[, c("Gene", "quadrant")], by="Gene")
  colnames(data)[2] <- "pancan_z"
  # join NSCLC quadrants
  data$quadrant <- data$quadrant %>% gsub(pattern="top_left", replacement="Q1") %>% gsub(pattern="bottom_left", replacement="Q2") %>% gsub(pattern="top_right", replacement="Q3") %>% gsub(pattern="bottom_right", replacement="Q4")
  data$quadrant <- factor(data$quadrant, levels = c("Q1", "Q2", "Q3", "Q4"))
  
  # box-plot
  boxplot_DB(data = data, x = "quadrant", y = "pancan_z", fill_var = "quadrant", x_axis = "Quadrant", y_axis = "Pan-cancer prognostic value\n(z score)")
}

## 4b)
{
  # input: PRECOG HR and P-vals by RNA heterogeneity quadrant in NSCLC
  precog <- precog %>% as.data.frame()
  rownames(precog) <- precog$`Cancer Type`
  
  # prep data
  quadrant_hr <- precog[, grep("_OR$", colnames(precog))]
  quadrant_pval <- precog[, grep("_p$", colnames(precog))]
  colnames(quadrant_hr) <- colnames(quadrant_pval) <- c("Top_left", "Bottom_left", "Top_right", "Bottom_right")
  
  # Volcano plots
  par(mfrow=c(2,2),las=1, mar=c(3,4.1,2,2.1))
  for(i in c(1,3,2,4)){
    cols <- c('#3182bd','#de2d26')[match(quadrant_hr[,i] > 1, c('FALSE','TRUE'))]
    cols[quadrant_pval[,i] > 0.05] <- '#bdbdbd'
    a <- quadrant_pval[,i]
    a[a < 0.000001] <- 0.000001
    plot(log(quadrant_hr[,i]), -log10(a), xlab='OR', main=colnames(quadrant_hr)[i],col=cols, pch=16,ylab='P-value', axes=F,xlim=c(-1,1), )
    axis(1)
    axis(2, at=-log10(c(1,0.05,0.01,0.001,0.0001,0.00001,0.000001)),labels=c(1,0.05,0.01,0.001,0.0001,0.00001,0.000001))
    box()
    abline(v=0,h=-log10(0.05),lty=3,col='#bdbdbd')
  }
}

## 4c)
{
  # input: rna and dna ith values, filtered for LUAD histology
  tmp <- dplyr::filter(rna_dna_corr, Histology=="Adeno")
  # scatter-plot 
  corr_plot_DB(data = tmp, x = "SCNA_ith", y = "RNA_ith", point_colour = "dodgerblue", point_border_col = "NA", point_alpha = 0.5, point_size = 7, title = "", x_axis = "Copy-Number ITH \n(% subclonal SCNA)", y_axis = "Gene Expression ITH", best_fit_line = T, x_limits = c(0,100), y_limits = c(0,0.5))
}

##  4d)
{
  #test genes
  gene_clonality <- gene_clonality_list[[2]]
  gene_clonality[,16] <- c(4,1,3,2)[match(gene_clonality[,16], c('bottom_right','top_left','top_right','bottom_left'))]
  gene_clonality <- apply(gene_clonality, 2, as.numeric)
  gene_clonality[is.na(gene_clonality)] <- 0
  rownames(gene_clonality) <- rownames(gene_clonality_list[[2]])
  #add annotation about ITH
  gene_clonality <- cbind(gene_clonality, rnaith_dhruva[rownames(gene_clonality),c('LUAD_within.heterogeneity','LUAD_between.heterogeneity','mean_expr')])
  #remove genes with few samples
  gene_clonality <- gene_clonality[rowSums(gene_clonality[,4:9]) > 30,]
  #turn into percentages
  gene_clonality[,4:9] <- gene_clonality[,4:9]/rowSums(gene_clonality[,4:9])
  anyGain <- rowSums(gene_clonality[,c('ClonalGain','SubclonalGain')])
  gene_clonality <- cbind(gene_clonality,anyGain)
  par(las=1,mfrow=c(2,3))
  for(j in c(4:8,20)){
    sumClonal <- setNames(gene_clonality[,j], rownames(gene_clonality))
    #subset to genes in upper or lower quadrant of clonal gains
    b <- rownames(gene_clonality[(sumClonal < quantile(sumClonal)[2] | sumClonal > quantile(sumClonal)[4]),])
    #b <- rownames(gene_clonality[(sumClonal < quantile(sumClonal)[3] | sumClonal > quantile(sumClonal)[3]),])
    #determine fraction in in higher quarter
    a <- b %in% rownames(gene_clonality[(sumClonal > quantile(sumClonal)[4]),])
    #a <- b %in% rownames(gene_clonality[(sumClonal > quantile(sumClonal)[3]),])
    #fishers test
    ftests_clonalgains <- list()
    for(k in 1:4){
      ftests_clonalgains[[k]] <- fisher.test(table(a, gene_clonality[b,'quadrant'] %in% k))
      
    }
    tmp <- unlist(lapply(ftests_clonalgains, function(x){ x$estimate}))
    tmp1 <- unlist(lapply(ftests_clonalgains, function(x){ x$p.value}))
    
    a <-barplot(tmp,main=colnames(gene_clonality)[j], col=c("firebrick1", "darkorchid2", "gold1", "turquoise2"),names.arg=c('Q1','Q2','Q3','Q4'),ylim=c(0,2),ylab='OR',xlab='')
    abline(h=1, lty=2)
    text(a, tmp+0.1, paste('P=',signif(tmp1,3),sep=''))
  }
}

##  4e)
{
  ##  input: scna by rna heterogeneity quadrant
  tmp <- quadrant_scna
  
  ##  prep data
  # log10-transform p-vals
  tmp$log_p_val <- -log10(tmp$p_val)
  # prep y-axis 
  y_max <- tmp$log_p_val %>% max() %>% ceiling()
  sig_values <- c(1, 0.05, 0.01, 10^(-seq(3, y_max)))
  
  ## volcano plot
  gg <- ggplot(tmp, aes(x=odds_ratio, y=log_p_val)) 
  # add lines
  gg <- gg + geom_hline(yintercept = -log10(0.05), linetype="dotted") 
  gg <- gg + geom_vline(xintercept = 1, linetype="dotted")
  # axes
  gg <- gg + scale_y_continuous(expand = c(0,0), limits = c(0, y_max), breaks = sapply(X=sig_values, FUN=function(x) {-log10(x)}), labels=sig_values)
  gg <- gg + scale_x_continuous( limits = c(0.5, 2), breaks = c(0.5, 1, 2), trans="log")
  # themes
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1, axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")), axis.text.x = element_text(size = 8, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),  axis.ticks.length = unit(-1.4, "mm"))
  # axis labels
  gg <- gg + xlab("Odds ratio") + ylab("P value")
  gg <- gg + geom_point(data = tmp, aes(colour=quadrant, fill=quadrant, pch=scna), alpha=1, size=1) + scale_shape_manual(values=c(22, 24, 2, 25, 6)) + scale_colour_manual(values=c("firebrick1", "darkorchid2", "gold1", "turquoise2")) + scale_fill_manual(values=c("firebrick1", "darkorchid2", "gold1", "turquoise2"))
  gg %>% print()
}

##  4f)
{
  #get all reactome pathways
  tmp <- read.table('ReactomePathways.20181008.gmt', sep='\t', as.is=T, fill=T)
  reactome_pathways <- list()
  for(i in 1:nrow(tmp)){
    reactome_pathways[[i]] <- as.character(tmp[i,-c(1:3)])
  }
  names(reactome_pathways) <- tmp[,2]
  
  entrez_id <- names(map2entrez(rownames(rnaith_dhruva)))
  rnaith_dhruva <- cbind(rnaith_dhruva, entrez_id)
  
  tmp <- setNames(c('quadrant', 'LUAD_quadrant','LUSC_quadrant'), c('All','LUAD','LUSC'))
  for(i in 1:3){
    cat(i)
    lower_right_q <- enrichPathway(gene=as.character(rnaith_dhruva[rnaith_dhruva[,tmp[i]] %in% 'bottom_right','entrez_id']),pvalueCutoff=0.05,  readable=T, pAdjustMethod='bonferroni', universe=unique(rnaith_dhruva[,'entrez_id']))
    upper_left_q <- enrichPathway(gene=as.character(rnaith_dhruva[rnaith_dhruva[,tmp[i]] %in% 'top_left','entrez_id']),pvalueCutoff=0.05, readable=T, pAdjustMethod='bonferroni', universe=unique(rnaith_dhruva[,'entrez_id']))
    upper_right_q <- enrichPathway(gene=as.character(rnaith_dhruva[rnaith_dhruva[,tmp[i]] %in% 'top_right','entrez_id']),pvalueCutoff=0.05, readable=T, pAdjustMethod='bonferroni', universe=unique(rnaith_dhruva[,'entrez_id']))
    lower_left_q <- enrichPathway(gene=as.character(rnaith_dhruva[rnaith_dhruva[,tmp[i]] %in% 'bottom_left','entrez_id']),pvalueCutoff=0.05, readable=T, pAdjustMethod='bonferroni', universe=unique(rnaith_dhruva[,'entrez_id']))
    
  }
  
  cnetplot(lower_right_q, showCategory = nrow(lower_right_q), categorySize="geneNum",vertex.label.cex=0.5)
  enrichMap(lower_right_q, vertex.label.cex=0.5)
  
  max_cat <- 5
  g <- cnetplot(lower_right_q, showCategory = max_cat, categorySize="pval",vertex.label.cex=0.5)
  vertex_attr(g)$color <- c(rep('#fed976',max_cat),rep('#9ecae1',length(vertex_attr(g)$size)-max_cat))
  vertex_attr(g)$size[1:max_cat] <- 15
  vertex_attr(g)$size[(max_cat+1):length(vertex_attr(g)$size)] <- 10
  
  plot(g, layout=layout.kamada.kawai(g),vertex.label.color='#000000',vertex.label.cex=0.5, edge.color='#bdbdbd')
}

################
##  Extended Data 1
################

##  ED-1a)
{
  #  TRACERx - RNASEQ
  surv_object_barplot_DB(surv_object = tmp_tracerx, title="TRACERx Lung - RNAseq")
}

##  ED-1b)
{
  #  TCGA - RNASEQ
  surv_object_barplot_DB(surv_object = tmp_tcga, title="The Cancer Genome Atlas - RNAseq")
  #  UPPSALA - RNASEQ
  surv_object_barplot_DB(surv_object = tmp_uppsala, title="Uppsala II NSCLC Cohort - RNAseq\n(Djureinovic JCI Insight 2016)")
}

##  ED-1c)
{
  #  DER - MICROARRAY
  surv_object_barplot_DB(surv_object = tmp_der, title="GSE50081 - microarray\n(Der Journal of Thoracic Oncology 2014)")
  #  OKAYAMA - MICROARRAY
  surv_object_barplot_DB(surv_object = tmp_okayama, title="GSE31210 - microarray\n(Okayama Cancer Research 2012)")
  #  ROUSSEAUX - MICROARRAY
  surv_object_barplot_DB(surv_object = tmp_rousseaux, title="GSE30219 - microarray\n(Rousseaux Science Translational Medicine 2013)")
  #  SHEDDEN - MICROARRAY
  surv_object_barplot_DB(surv_object = tmp_shedden, title="GSE68465 - microarray\n(Shedden Nature Medicine 2008)")
}

################
##  Extended Data 2
################

##  ED-2a)
{
  # input: mat1, tracerx normalised count
  mat1 <- tracerx.counts_vsd_filter %>% t() %>% as.data.frame() 
  # calculate Z-score
  tmp <- sapply(X=mat1, FUN=scale)
  tmp <- as.data.frame(tmp)
  rownames(tmp) <- rownames(mat1)
  mat1 <- tmp
  
  # make mat2
  mat2 <- mat1 %>% rownames_to_column("SampleID")
  #add PatientID and RegionID cols
  mat2$PatientID <- sapply(X=mat2$SampleID, FUN=function(x) {unlist(strsplit(x, split=":"))[1]})
  mat2$RegionID <- sapply(X=mat2$SampleID, FUN=function(x) {unlist(strsplit(x, split=":"))[2]})
  #join histology data
  mat2 <- dplyr::left_join(x = mat2, y=tracerx.clinical, by="PatientID")
  #class histology data
  mat2$Histology_class <- ifelse(mat2$Histology=="Invasive adenocarcinoma", "LUAD", ifelse(mat2$Histology=="Squamous cell carcinoma", "LUSC", "Other"))
  #select cols
  mat2 <- dplyr::select(mat2, SampleID, PatientID, Histology_class)
  #save order arranged by Histology_class
  idx <- dplyr::arrange(mat2, Histology_class)
  #create Colour data-frame
  col <- idx %>% dplyr::select(PatientID, Histology_class) %>% dplyr::distinct()
  tmp <- ifelse(col$Histology_class=="LUAD", "green", ifelse(col$Histology_class=="LUSC", "blue", ifelse(col$Histology_class=="Other", "gold", NA)))
  names(tmp) <- col$PatientID
  col <- tmp
  #spread
  mat2 <- spread(mat2, SampleID, Histology_class) 
  #transpose
  tmp <- mat2[,-1] %>% t() %>% as.data.frame()
  colnames(tmp) <- mat2$PatientID
  mat2 <- tmp
  mat2 <- mat2[, rev(tmp)]
  
  ##  ComplexHeatmap functions
  
  #ht1: expression heatmap
  # specify colours
  colors = structure(c("green4", "darkorchid3", "darkorange"), names = c("LUAD", "LUSC", "Other"))
  #reorder cols and rows
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  mat_cluster_cols <- hclust(dist(t(mat1)))
  mat_cluster_rows <- sort_hclust(hclust(dist(mat1)))
  # make heatmap-1
  ht1 = Heatmap(
    mat1, 
    name = "ht1", 
    col = colorRamp2(seq(-2,4,1.5), viridis(5)),
    column_dend_height = unit(15, "mm"),
    row_dend_width = unit(15, "mm"),
    row_dend_reorder = rev(mat_cluster_rows$order), 
    column_dend_reorder = rev(mat_cluster_cols$order),
    show_column_names = FALSE,
    show_row_names = FALSE,
    width = 1,
    heatmap_legend_param = list(title = NULL, color_bar = "continuous")
  )
  
  #ht2: tumour region heatmap
  # specify colours
  colors = structure(c("green4", "darkorchid3", "darkorange"), names = c("LUAD", "LUSC", "Other"))
  # make heatmap-2
  ht2 = Heatmap(
    matrix=mat2,
    col = colors,
    name = "ht2",
    na_col = "white",
    rect_gp = gpar(col="white"),
    show_row_names = FALSE, 
    show_column_names = FALSE,
    cluster_columns = FALSE,
    heatmap_legend_param = list(title = NULL, color_bar = "discrete"),
    column_names_gp = gpar(fontsize = 4),
    width = 1
  )
  
  # plot ComplexHeatmaps 
  ht1 + ht2
  #ht1: add border
  decorate_heatmap_body("ht1", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.75))})
  #ht2: add border, and lines
  tmp <- dplyr::select(idx, PatientID, Histology_class) %>% dplyr::distinct()
  decorate_heatmap_body("ht2", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.75))})
}

##  ED-2b)
{
  # split patients into 2 groups
  bin <- cutree(mat_cluster_rows, k=2)
  
  # make RiskScore data-frame
  RiskScore <- data.frame(SampleID=names(bin), bin=bin, stringsAsFactors = FALSE)
  # add PatientID col
  RiskScore$PatientID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split=":"))[1]})
  # collapse rows onto PatientID and bin
  RiskScore <- RiskScore %>% dplyr::distinct(PatientID, bin)
  nrow(RiskScore) # QC: should have 48 TRACERx patients
  # join survival data
  RiskScore <- dplyr::left_join(x=RiskScore, y=tracerx.clinical, by="PatientID")
  
  # log-rank test
  sdf <- survdiff(Surv(RiskScore$OS, RiskScore$status==1) ~ RiskScore$bin, data = RiskScore)
  log.rank_p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1) 
  log.rank_p <- log.rank_p %>% signif(digits=3)
  
  #  km plot: discordant_low
  #tmp_title <- paste0("Top 500 Variably Expressed Genes\nTRACERx (n=", nrow(RiskScore), " NSCLC Patients)\nlog-rank P value = ",log.rank_p)
  ggsurvplot(survfit(Surv(RiskScore$OS, RiskScore$status==1) ~ RiskScore$bin, data = RiskScore), risk.table = TRUE, legend.labs = c("High", "Low"), tables.theme = theme_cleantable(), palette = "aaas", title=paste0("Top 500 Variably Expressed Genes\nTRACERx (n=", nrow(RiskScore), " NSCLC Patients)\nlog-rank P value = ",log.rank_p), legend="none", xlab="Time (Months)", ylab="Overall Survival Probability", break.time.by=10)  %>% print()
}

##  ED-2c)
{
  ##  input: 
  # specify Pervenio Lung RS signature
  signature <- dplyr::filter(prog, signature=="Kratz_Lancet_2012")$Gene
  signature_name <- "PERVENIO LUNG RS"
  # font size
  font_size <- 6
  
  ##   mat1
  mat1 <- tracerx.counts_vsd[which(rownames(tracerx.counts_vsd) %in% signature), ] %>% t() %>% as.data.frame() 
  #filter cols: Histology==LUAD patients
  tmp <- dplyr::filter(tracerx.clinical, Histology=="Invasive adenocarcinoma")
  mat1 <- mat1[sapply(X=tmp$PatientID, FUN=function(x) {grep(x, rownames(mat1))}) %>% unlist(), ]
  #Z-score
  tmp <- sapply(X=mat1, FUN=scale)
  tmp <- as.data.frame(tmp)
  rownames(tmp) <- rownames(mat1)
  mat1 <- tmp
  #sort rows
  mat1 <- mat1[sort(rownames(mat1)), ]
  
  ##  mat2
  mat2 <- mat1 %>% rownames_to_column("SampleID")
  #add PatientID and RegionID cols
  mat2$PatientID <- sapply(X=mat2$SampleID, FUN=function(x) {unlist(strsplit(x, split=":"))[1]})
  mat2$RegionID <- sapply(X=mat2$SampleID, FUN=function(x) {unlist(strsplit(x, split=":"))[2]})
  #join histology data
  #mat2 <- dplyr::left_join(x = mat2, y=tx100.manifest)
  mat2 <- dplyr::left_join(x = mat2, y=tracerx.clinical, by="PatientID")
  #class histology data
  mat2$Histology_class <- ifelse(mat2$Histology=="Invasive adenocarcinoma", "LUAD", ifelse(mat2$Histology=="Squamous cell carcinoma", "LUSC", "Other"))
  #select cols
  mat2 <- dplyr::select(mat2, SampleID, PatientID, Histology_class)
  #save order arranged by Histology_class
  idx <- dplyr::arrange(mat2, Histology_class)
  #create Colour data-frame
  col <- idx %>% dplyr::select(PatientID, Histology_class) %>% dplyr::distinct()
  tmp <- ifelse(col$Histology_class=="LUAD", "green", ifelse(col$Histology_class=="LUSC", "blue", ifelse(col$Histology_class=="Other", "gold", NA)))
  names(tmp) <- col$PatientID
  col <- tmp
  #spread
  mat2 <- spread(mat2, SampleID, Histology_class) 
  #transpose
  tmp <- mat2[,-1] %>% t() %>% as.data.frame()
  colnames(tmp) <- mat2$PatientID
  mat2 <- tmp
  mat2 <- mat2[, tmp]
  
  #ht1: expression heatmap
  colors = structure(c("green4", "darkorchid3", "darkorange"), names = c("LUAD", "LUSC", "Other"))
  # make heatmap
  ht1 = Heatmap(
    mat1, 
    name = "ht1", 
    col = colorRamp2(seq(-3,3), viridis(7, option = "magma")),
    column_dend_height = unit(15, "mm"),
    row_dend_width = unit(15, "mm"),
    show_row_names = FALSE,
    width = 1,
    heatmap_legend_param = list(title = NULL, color_bar = "continuous"),
    row_title = signature_name,
    column_names_gp = gpar(fontsize = font_size), 
    cluster_columns = F
  )
  
  #ht2: tumour region heatmap
  colors = structure(c("green4", "darkorchid3", "darkorange"), names = c("LUAD", "LUSC", "Other"))
  #reorder rows
  mat_cluster_rows <- hclust(dist(mat1))
  tmp <- data.frame(SampleID=mat_cluster_rows$labels, order=mat_cluster_rows$order) %>% dplyr::arrange(order)
  tmp$PatientID <- sapply(X=tmp$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split=":"))[1]})
  tmp <- tmp %>% dplyr::group_by(PatientID) %>% dplyr::top_n(n=1, wt=order)
  # make heatmap-2
  ht2 = Heatmap(
    matrix=mat2,
    col = colors,
    name = "ht2",
    na_col = "white",
    rect_gp = gpar(col="white"),
    show_row_names = FALSE, 
    show_column_names = FALSE,
    cluster_columns = FALSE,
    heatmap_legend_param = list(title = NULL, color_bar = "discrete"),
    column_names_gp = gpar(fontsize = font_size),
    width = 1
  )
  
  ##  plot ComplexHeatmaps 
  ht1 + ht2
  #ht1: add border
  decorate_heatmap_body("ht1", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.75))})
  #ht2: add border, and lines
  tmp <- dplyr::select(idx, PatientID, Histology_class) %>% dplyr::distinct()
  decorate_heatmap_body("ht2", { for (i in 1:nrow(tmp)-1) { grid.lines(x = unit(c(i/nrow(tmp), i/nrow(tmp)), "npc"), y = unit(c(0, 1), "npc"), gp = gpar(lty="dotted", col="gray75", lwd=0.75)) } ; grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.75)) })
}

##  ED-2d)
{
  # input: sample-ids for luad patients
  luad_idx <- dplyr::filter(tracerx.clinical, Histology=="Invasive adenocarcinoma")$PatientID
  luad_idx <- sapply(X=luad_idx, FUN=function(x) {grep(paste0("^", x, ":"), colnames(tracerx.counts_vsd), value = T)}) %>% unlist() %>% as.character()
  
  ##  prep tsb_luad for plotting
  # add n_genes
  for(i in 2:length(colnames(tsb_luad))) {
    colnames(tsb_luad)[i] <- paste0(colnames(tsb_luad)[i], " (n=", n_genes$n_genes[which(n_genes$signature == colnames(tsb_luad)[i])], ")")  
  }
  # melt
  tsb_luad <- melt(tsb_luad, "no_clusters")
  colnames(tsb_luad)[2:3] <- c("signature", "prop_same_cluster")
  
  # plot
  # calculate no. patients
  no_patients <- luad_idx %>% sapply(FUN=function(x) {unlist(strsplit(x, split=":"))[1]}) %>% unique() %>% length()
  #make plot
  gg  <- ggplot(tsb_luad, aes(x=no_clusters, y=prop_same_cluster, col=signature)) 
  gg <- gg + geom_vline(xintercept = c(2, 3, round(no_patients/2), no_patients), linetype="dashed")
  # add data line
  gg <- gg + geom_line() 
  # add themes
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + theme(aspect.ratio = 1, panel.grid.minor = element_blank()) + theme(legend.position = "bottom") + guides(colour=guide_legend(title=NULL))
  # add plot labels
  gg <- gg + xlab("No. of clusters") + ylab("Proportion of patients with\nall samples in the same cluster") + ggtitle(label = "Hclust assessment of tumour sampling bias", subtitle = paste0(length(unique(tsb_luad$signature)), " LUAD signatures")) + theme(plot.title=element_text(hjust=0.5, face="bold"))
  p1 <- gg + xlim(0, ceiling(no_patients/10)*10)
  # print
  gg %>% print()
  
  # print median discordance rate
  tmp <- tsb_luad %>% dplyr::filter(no_clusters==28)
  med_discordance <- round((1-median(tmp$prop_same_cluster))*100, digits=1)
  print(paste0("tsb_luad median discordance rate = ", med_discordance, "%"))
}

################
##  Extended Data 3
################

##  ED-3b)
{
  # correlate gene-wise MAD scores with SD scores
  corr_plot_DB(data = intra_var, x = "sd", y = "mad", title = "", point_size = 2, x_axis = "", y_axis = "Median absolute deviation", best_fit_line = T, aspect_ratio = NULL)
  # correlate gene-wise CV scores with SD scores
  corr_plot_DB(data = intra_var, x = "sd", y = "cv", title = "", point_size = 2, x_axis = "Standard deviation", y_axis = "Coefficient of variation", best_fit_line = T, aspect_ratio = NULL)
}

##  ED-3c)
{
  # correlate SD scores, TCGA ~ TRACERx
  corr_plot_DB(data = inter_var, x = "TRACERx.SD", y = "TCGA.SD", corr_method = "pearson", point_size = 2, best_fit_line = T, x_axis = "TRACERx\nintertumour heterogeneity scores", y_axis = "TCGA\nintertumour heterogeneity scores", point_colour = "gold", point_border_col = "gold")
}

################
##  Extended Data 4
################

##  ED-4a)
{
  ## identify Tukey's 5 genes
  # sort by AUC
  tmp <- hclust_auc[, c("Gene", "AUC")]
  tmp <- tmp[order(tmp$AUC), ]
  n <- nrow(tmp)
  tmp <- tmp[c(1, ceiling(n/4), ceiling(n/2), ceiling(3*n/4), n), ]
  
  ## scatter-plot for all genes
  gg <- ggplot(hclust_auc, aes(x=reorder(Gene, AUC), y=AUC)) + geom_point(pch=21, fill="gray75", colour="gray75", alpha=0.5) 
  gg <- gg + geom_point(data = tmp, pch=21, fill="red", colour="red") + geom_text_repel(data = tmp, aes(label=Gene), colour="red")
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(paste0("Gene", "\n(", nrow(hclust_auc), " in total)")) + ylab("Hierarchical clustering concordance\n(integral of curve)") + scale_y_continuous(expand = c(0,0), limits = c(0, 40), breaks = seq(0, 40, 10)) + theme(aspect.ratio = 0.5)
  gg %>% print()
  
  ## line-plot: AUC curves for Tukey's 5 genes
  gene_cluster_DB(signature=tmp$Gene, draw_plots = T) 
}

##  ED-4b)
{
  # join hclust auc scores with heterogeneity scores (LUAD)
  gene_list <- dplyr::inner_join(x=hclust_auc, y=ith_comparison_luad, by="Gene") %>% as.data.frame()
  nrow(gene_list) #QC
  # order quadrants
  gene_list$quadrant <- factor(gene_list$quadrant, levels=c("top_left", "bottom_left", "top_right", "bottom_right"))
  
  # wilcox test
  data = gene_list
  x = "quadrant"
  y = "AUC"
  # calculate p-values for all comparisons
  pv <- pairwise.wilcox.test(x = data[, y, drop=T], g = data[, x, drop=T], p.adjust.method="none", paired=F)$p.value
  # tidy: create "quad_1" and "quad_2" cols to show comparison for each test
  pv <- pv %>% as.data.frame() %>% rownames_to_column("quad_1")
  pv <- melt(pv, id.vars="quad_1")
  colnames(pv)[2:3] <- c("quad_2", "p_val")
  # filter significant values, and label with asterisks
  pv <- pv[which(pv$p_val < 0.05), ]
  pv$symbol <- ifelse(pv$p_val > 0.01, "*", ifelse(pv$p_val > 0.001, "**", "***"))
  # fix class of "quad_2" col
  pv$quad_2 <- pv$quad_2 %>% as.character()
  tmp <- as.data.frame(t(pv[, 1:2]))
  tmp[] <- sapply(tmp, as.character)
  
  #box-plot
  gg <- ggplot(gene_list, aes(x=quadrant, y=AUC, fill=quadrant))
  gg <- gg + geom_boxplot(position = "dodge", lwd=0.5, outlier.size = NA) + scale_fill_manual(values=c("firebrick1", "darkorchid2", "gold1", "turquoise2"))
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + theme(aspect.ratio = 1) + theme(legend.position = "none")
  gg <- gg + xlab("Quadrant") + ylab("Clustering concordance score per gene")
  gg <- gg + ylim(0,20)
  gg %>% print()
  gg <- gg + geom_signif(comparisons = as.list(tmp), annotation=pv$symbol, y_position = seq(15,(15+ncol(tmp)),1), tip_length = 0)
  gg %>% print()
}

################
##  Extended Data 5
################

##  ED-5a)
{
  # first run code for Fig. 2b, need q_comparison and fisher_p objects
  
  # count no. of genes per signature
  n_genes <- data.frame(variable=c("expected", "observed"), n_genes=c(sum(q_comparison$expected_count), sum(q_comparison$observed_count)))
  
  # specify quadrant colours
  cols <- c("firebrick1", "darkorchid2", "gold1", "turquoise2")
  
  # extract signature info
  tmp <- table(prog_histology$quadrant, prog_histology$signature) %>% data.frame()
  colnames(tmp) <- c("quadrant", "signature", "count")
  # calculate percentages
  tmp <- tmp %>% dplyr::group_by(signature) %>% dplyr::mutate(percent=count/sum(count)*100) %>% ungroup
  # make quadrant a factor
  tmp$quadrant <- factor(tmp$quadrant, levels=c("top_left", "bottom_left", "top_right", "bottom_right"))
  
  # tidy signature names for plotting
  tmp$signature <- tmp$signature %>% as.character() %>% gsub(pattern="_", replacement="\n") %>% gsub(pattern="\\.", replacement="\n")
  # add n-genes col
  n_genes <- table(prog_histology$signature) %>% data.frame()
  colnames(n_genes) <- c("signature", "n_genes")
  tmp <- dplyr::left_join(x=tmp, y=n_genes, by="signature")
  # add n_genes to signature name
  tmp$signature <- paste0(tmp$signature, "\n(n=", tmp$n_genes, ")")
  # order columns by Q4
  tmp$signature <- factor(tmp$signature, levels=dplyr::arrange(dplyr::filter(tmp, quadrant=="bottom_right"), percent)$signature)
  
  ##  stacked bar-plot per signature
  gg <- ggplot(tmp, aes(x=signature, y=percent, fill=quadrant)) 
  # underlay horizontal lines
  gg <- gg + geom_hline(yintercept = c(25, 50, 75), colour="azure4", lty="dashed")
  # add bars
  gg <- gg + geom_bar(stat="identity", alpha=0.8) + theme_classic()
  # colours
  gg <- gg + scale_fill_manual(values=cols)
  # plot labels
  gg <- gg + ggtitle(label = "", subtitle = paste0(length(unique(prog_histology$signature)), " signatures (", sum(q_comparison$observed_count, na.rm = T), " total genes)")) + theme(plot.title = element_text(hjust=0.5, face="bold")) + ylab("Prognostic signature genes per\nheterogeneity quadrant (%)") + xlab("Prognostic signature")
  # themes
  gg <- gg  + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + theme(aspect.ratio=1) + theme(legend.position = "none") + theme(aspect.ratio = 0.5)
  gg %>% print()
}

##  ED-5b)
{
  # first run code for Fig. 2b, need q_comparison and fisher_p objects
  
  # count no. of genes per signature
  n_genes <- data.frame(variable=c("expected", "observed"), n_genes=c(sum(q_comparison$expected_count), sum(q_comparison$observed_count)))
  
  # specify quadrant colours
  cols <- c("firebrick1", "darkorchid2", "gold1", "turquoise2")
  
  ##  prep data for ggplot2
  # select percent cols
  q_comparison_melt <- dplyr::select(q_comparison, quadrant, expected_percent, observed_percent)
  # re-order quadrant factor
  q_comparison_melt$quadrant <- factor(q_comparison_melt$quadrant, levels=c("top_left", "bottom_left", "top_right", "bottom_right"))
  # melt
  q_comparison_melt <- melt(q_comparison_melt, id.vars = "quadrant")
  # fix names
  q_comparison_melt$variable <- gsub(x=q_comparison_melt$variable, pattern="_percent", replacement="")
  colnames(q_comparison_melt)[3] <- "percent"
  
  ##  stacked bar-plot  
  gg <- ggplot(q_comparison_melt, aes(x=variable, y=percent, fill=quadrant)) + geom_bar(stat='identity', col="black", width=0.5) 
  # colours
  gg <- gg + scale_fill_manual(values=cols) + scale_colour_manual(values=cols)
  # themes
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + theme(aspect.ratio=1) + theme(legend.position = "none")
  # plot labels
  gg <- gg + ggtitle(label = "", subtitle = paste0(length(unique(prog_histology$signature)), " signatures (", sum(q_comparison$observed_count, na.rm = T), " total genes)\nQ4 enrichment: Fisher's exact P = ", signif(fisher_p, digits=3))) + theme(plot.title = element_text(hjust=0.5, face="bold")) + ylab("Prognostic signature genes\nper heterogeneity quadrant (%)") + xlab("")
  # x axis labels
  gg <- gg + scale_x_discrete(labels=c("expected"=paste0("Expected\n(n=", n_genes$n_genes[which(n_genes$variable=="expected")], ")"), "observed"=paste0("Observed\n(n=", n_genes$n_genes[which(n_genes$variable=="observed")], ")")))
  # add expected labels
  tmp <- dplyr::filter(q_comparison_melt, variable=="expected")
  tmp <- tmp[order(tmp$quadrant, levels(tmp$quadrant), decreasing = T),]
  tmp$cumsum <- tmp$percent %>% cumsum()
  gg <- gg + geom_label(data = tmp, aes(y=(cumsum-percent*0.5), label=paste0(round(percent), "%"), colour=quadrant), nudge_x = -0.4, fill="white")
  gg <- gg + geom_text(data = tmp, aes(y=(cumsum-percent*0.5), label=paste0(round(percent), "%")), nudge_x = -0.4, colour="black")
  # add observed labels
  tmp <- dplyr::filter(q_comparison_melt, variable=="observed")
  tmp <- tmp[order(tmp$quadrant, levels(tmp$quadrant), decreasing = T),]
  tmp$cumsum <- tmp$percent %>% cumsum()
  gg <- gg + geom_label(data = tmp, aes(y=(cumsum-percent*0.5), label=paste0(round(percent), "%"), colour=quadrant), nudge_x = 0.4, fill="white")
  gg <- gg + geom_text(data = tmp, aes(y=(cumsum-percent*0.5), label=paste0(round(percent), "%")), nudge_x = 0.4, colour="black")
  gg <- gg + guides(colour="none")
  gg %>% print()
}

##  ED-5c-f)
{
  # make prog_histology df
  prog_histology <- prog
  prog_histology$Histology <- "LUAD"
  
  # shedden
  quadrant_hr_boxplot_DB(prog_histology, shedden_res_luad, "LUAD", "Shedden Nat Med 2008")
  # okayama
  quadrant_hr_boxplot_DB(prog_histology, okayama_res_luad, "LUAD", "Okayama CR 2012")
  # der
  quadrant_hr_boxplot_DB(prog_histology, der_res_luad, "LUAD", "Der JTO 2014")
  # rousseaux
  quadrant_hr_boxplot_DB(prog_histology, rousseaux_res_luad, "LUAD", "Rousseaux STM 2013")
}

################
##  Extended Data 6
################

##  ED-6c)
{
  # line-plot: cross-validation error in training cohort ~ no. of genes (ordered by clustering concordance)
  gg <- ggplot(cv_error, aes(x=AUC, y=cv_error)) 
  gg <- gg + geom_vline(xintercept = 90, colour="red", linetype="dashed")
  gg <- gg + geom_line(colour="gray75") #+ geom_point()
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + theme(aspect.ratio = 1)
  gg <- gg + xlab("No. of genes (shortlisted by clustering concordance score)") + ylab("Cross-Validation Error in Training Cohort(glmnet 10-fold CV mean error)") 
  gg <- gg + scale_y_continuous(expand = c(0,0), limits = c(12.52, 12.62)) + scale_x_continuous(expand = c(0,0), breaks = seq(10,100,10), limits=c(0,110))
  gg %>% print()
}

##  ED-6d)
{
  # specify significance values for y-axis 
  sig_values <- c(1, 0.01, 1e-5, 1e-10, 1e-15, 1e-20) 
  
  # extract optimal risk-score: average (median) RiskScore threshold amongst significant (log-rank P < 0.01) splits
  tmp <- dplyr::filter(riskscore_calibration, log_rank_p.value < 0.01)
  cutoff <- tmp$riskscore_threshold %>% median()
  
  # line-plot: log-rank p-value in training cohort ~ riskscore cut-off
  gg <- ggplot(riskscore_calibration, aes(x=riskscore_threshold, y=p_adj) )
  gg <- gg + geom_hline(yintercept = -log10(0.01), linetype="dotted")
  gg <- gg + geom_vline(xintercept = cutoff, colour="red", linetype="dashed")
  gg <- gg + geom_line(colour="gray75") #+ geom_point()
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + theme(aspect.ratio = 1)
  gg <- gg + xlab("ORACLE riskscore cut off") + ylab("Log rank p value") 
  gg <- gg + scale_y_continuous(breaks = sapply(X=sig_values, FUN=function(x) {-log10(x)}), labels=sig_values, limits = c(-log10(max(sig_values)), -log10(min(sig_values))), expand = c(0,0))
  gg %>% print()
}

##  ED-6e)
{
  # input: oracle risk-score in tracerx
  RiskScore <- tx_rs
  
  # input: oracle risk-score cut-off 
  riskscore_thresh <- 10.19941
  
  # count no. of patients per class
  #tmp <- RiskScore %>% dplyr::select(PublicationID, class) %>% distinct()
  tmp <- RiskScore %>% dplyr::select(PublicationID, class) %>% distinct()
  
  # scatter-plot
  gg <- ggplot(RiskScore, aes(x=fct_reorder(PublicationID, RiskScore + as.numeric(class), .fun=min), y=RiskScore)) 
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")), axis.text.x = element_text(size = 8, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),  axis.ticks.length = unit(-1.4, "mm")) + scale_y_continuous(expand = c(0,0), limits=c(floor(min(RiskScore$RiskScore)), ceiling(max(RiskScore$RiskScore))))
  gg <- gg + theme(axis.text.x = element_text(angle=90, vjust=0.5))
  gg <- gg + geom_hline(yintercept = riskscore_thresh, col="black", lty="dotted")
  gg <- gg + geom_line(col="black")
  gg <- gg + ggtitle(label="", subtitle = paste0(length(unique(RiskScore$PublicationID)), " TRACERx LUAD patients = ", table(tmp$class)["Low"], " low + ", table(tmp$class)["High"], " high + ", table(tmp$class)["Discordant"], " discordant")) + theme(plot.title = element_text(hjust=0.5, face="bold")) + xlab("PatientID") + ylab("ORACLE RiskScore")
  gg <- gg + theme(legend.position = "bottom") + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) #border
  gg <- gg + theme(aspect.ratio = 0.5)
  gg <- gg + geom_point(pch=16, aes(col=class), size=3, alpha=0.5) + scale_color_manual(values = c("#3B4992FF", "azure4", "#EE0000FF")) + theme(legend.position = "none")
  gg %>% print()
}

################
##  Extended Data 7
################

##  ED-7a)
{
  ##  input: risk-scores in uppsala cohort
  cohort_riskscore <- test_rs
  clinical_data <- test_clinical
  
  #  make title
  tmp <- subset(clinical_data, PatientID %in% test_rs$PatientID)
  title_full <-  paste("ORACLE\nValidation Data-set", " (", nrow(test_rs), " patients, stage ", min(tmp$Stage_numeric, na.rm = T), "-", max(tmp$Stage_numeric, na.rm = T), ")", sep="")
  
  #stratify cohort
  cohort_riskscore$PatientID <- cohort_riskscore$PatientID %>% as.character()
  #join survival analysis
  cohort_riskscore <- left_join(x=cohort_riskscore, y=clinical_data, by="PatientID")
  # convert bin col to factor
  cohort_riskscore$bin <- factor(cohort_riskscore$bin, levels=c("Low", "High"))
  
  # log-rank test
  sdf <- survdiff(Surv(cohort_riskscore$OS, cohort_riskscore$status==1) ~ cohort_riskscore$bin, data = cohort_riskscore)
  log.rank_p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1) 
  log.rank_p <- log.rank_p %>% signif(digits=2)
  
  # KM-plot
  ggsurvplot(survfit(Surv(cohort_riskscore$OS, cohort_riskscore$status==1) ~ cohort_riskscore$bin, data = cohort_riskscore), risk.table = TRUE, legend.labs = c("High", "Low"), tables.theme = theme_cleantable(), palette = "aaas", title=paste(title_full, "\nlog-rank p = ", log.rank_p, sep=""), legend="none", xlab="Time (Years)", ylab="Overall Survival Probability") %>% print()
}

##  ED-7b-c)
{
  ##  input: filter stage 1 patients
  tmp <- dplyr::left_join(x=test_rs, y=test_clinical_tnm8, by="PatientID")
  tmp <- tmp %>% dplyr::filter(Stage_numeric==1)
  cohort_riskscore <- tmp
  
  ##  oracle risk-score test
  # log-rank P
  sdf <- survdiff(Surv(cohort_riskscore$OS, cohort_riskscore$status==1) ~ cohort_riskscore$bin, data = cohort_riskscore)
  log.rank_p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1) %>% signif(digits=3)
  # KM plot: ORACLE
  ggsurvplot(survfit(Surv(cohort_riskscore$OS, cohort_riskscore$status==1) ~ cohort_riskscore$bin, data = cohort_riskscore), risk.table = TRUE, legend.labs = c("High", "Low"), tables.theme = theme_cleantable(), palette = "aaas", title=paste("ORACLE\nlog-rank p = ", log.rank_p, sep=""), legend="none", xlab="Time (Years)", ylab="Overall Survival Probability") %>% print()
  
  ## sub-stage test
  # log-rank P
  sdf <- survdiff(Surv(cohort_riskscore$OS, cohort_riskscore$status==1) ~ cohort_riskscore$Stage, data = cohort_riskscore)
  log.rank_p_stage <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1) %>% signif(digits=3)
  # KM plot: substage
  ggsurvplot(survfit(Surv(cohort_riskscore$OS, cohort_riskscore$status==1) ~ cohort_riskscore$Stage, data = cohort_riskscore), risk.table = TRUE, tables.theme = theme_cleantable(), palette = "aaas", title=paste("Substaging\nlog-rank p = ", log.rank_p_stage, sep=""), legend="none", xlab="Time (Years)", ylab="Overall Survival Probability") %>% print()
}

##  ED-7d)
{
  ##  join cohorts
  # prep test_rs
  test_rs <- dplyr::left_join(x=test_rs, y=dplyr::select(test_clinical, PatientID, Stage, Stage_numeric), by="PatientID")
  test_rs$cohort <- "UPPSALA"
  # row-bind uppsala and met500 data
  tmp <- rbind(test_rs, met_rs)
  #  make cohort col a factor
  tmp$cohort <- factor(tmp$cohort, levels = c("UPPSALA", "MET500"))
  
  ##  get scale values
  y_min <- min(tmp$RiskScore) %>% floor()
  y_max <- max(tmp$RiskScore) %>% ceiling() + 1
  ## calculate n-numbers
  n <- table(tmp$Stage_numeric)
  n <- data.frame(x_old=names(n), x_new=paste0(names(n), "\n(n=", n, ")"), stringsAsFactors = F)
  
  ##  calculate P-values
  pv <- pairwise.wilcox.test(x = tmp$RiskScore, g = tmp$Stage_numeric, p.adjust.method="none", paired=F)$p.value
  # tidy: create "quad_1" and "quad_2" cols to show comparison for each test
  pv <- pv %>% as.data.frame() %>% rownames_to_column("quad_1")
  pv <- melt(pv, id.vars="quad_1")
  colnames(pv)[2:3] <- c("quad_2", "p_val")
  # QC: output to console
  print(pv)
  # filter significant values, and label with asterisks
  pv <- pv[which(pv$p_val < 0.05), ]
  pv$symbol <- ifelse(pv$p_val > 0.01, "*", ifelse(pv$p_val > 0.001, "**", "***"))
  # fix class of "quad_2" col
  pv$quad_2 <- pv$quad_2 %>% as.character()
  # QC: output to console
  print(pv)
  #
  p_val <- as.data.frame(t(pv[, 1:2]))
  p_val[] <- sapply(p_val, as.character)
  
  #  beeswarm plot
  gg <- ggplot(tmp, aes(x=factor(Stage_numeric), y=RiskScore)) 
  gg <- gg + geom_hline(yintercept = 10.19941, linetype="dotted")
  gg <- gg + geom_boxplot(fill=NA, outlier.shape = NA) 
  gg <- gg + geom_beeswarm(aes(colour=cohort), size = 1) + scale_color_manual(values = c("orange", "dodgerblue", "palegreen3"))
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + theme(aspect.ratio = 1)
  gg <- gg + xlab("Disease Stage") + ylab("ORACLE risk score")
  gg <- gg + scale_y_continuous(limits=c(y_min, y_max), expand = c(0,0))
  gg <- gg + scale_x_discrete(labels = n$x_new)
  gg <- gg + geom_signif(comparisons = as.list(p_val), annotation=pv$symbol, tip_length = 0)
  gg %>% print()
}

##  ED-7e)
{
  #  join oracle risk-scores and ki67 staining percentage
  tmp <- join_all(dfs=list(tx_rs, tx.nature[, c("PublicationID", "Ki67_percentage")]), by="PublicationID")
  # scatter plot
  corr_plot_DB(data = tmp, x="Ki67_percentage", y="RiskScore", x_axis = "Ki67 percentage", y_axis = "ORACLE RiskScore", title = "", best_fit_line = TRUE)
}

################
##  Extended Data 8
################

##  ED-8a)
{
  ##  heatmap
  gg <- ggplot(oracle_immune_corr, aes(x=immune_subset, y=0, fill=rho)) 
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.25)) + labs(fill="Spearman Rho") + coord_equal()
  gg <- gg + xlab("Immune Subset") + theme(axis.ticks=element_blank()) + theme(axis.text.x = element_text(angle=90, hjust=1))
  gg <- gg + ylab("") + theme(axis.text.y = element_blank())
  gg <- gg + geom_tile(color="white", size=1)
  gg <- gg + geom_text(aes(label=sig), color="white", size=5, y=-0.1)
  gg %>% print()
}

##  ED-8b)
{
  # scatter plot
  corr_plot_DB(data = ascat, x="ASCAT.purity", y="RiskScore", x_axis = "WES tumour purity", y_axis = "ORACLE risk score", title = "", best_fit_line = TRUE)
}

##  ED-8c)
{
  # specify marker genes for 7 stromal cell-types: "Alveolar", "Bcell", "Epithelial", "Fibroblast", "Myeloid", "Tcell", "Vascular" 
  marker_genes <- c("AGER", "MS4A1", "EPCAM", "COL6A2", "CD68", "CD3D", "FLT1")
  # include ORACLE genes
  genes <- c(marker_genes, oracle$Gene)
  
  # extract these genes from Lambrechts et al cluster data
  db <- sc_expression %>% rownames_to_column("Gene")
  db <- dplyr::filter(db, Gene %in% genes)
  # make Gene col a factor
  db$Gene <- factor(db$Gene, levels = genes)
  # melt
  db <- melt(db, id.vars = "Gene")
  # add celltype col
  db$variable <- db$variable %>% gsub(pattern="EC", replacement="Vascular")
  db$celltype <- gsub("[0-9]", "", db$variable) 
  
  # bar-plot, faceted by gene and coloured by stromal expression cluster
  ggplot(db, aes(x=variable, y=value, fill=celltype)) + facet_wrap(~Gene, ncol = length(marker_genes)) + geom_bar(stat="identity", position="dodge") + theme_bw() + scale_fill_brewer(palette="Dark2") + ggtitle("Lambrechts Nature Medicine 2018") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + xlab("Cluster") + ylab("Gene Expression (relative units)") + theme(aspect.ratio=0.5) + theme(legend.position="bottom")
}

##  ED-8d)
{
  # bar-plot: pearson's r ~ oracle genes, coefficients of correlation between expression and copy-number state of oracle genes
  ggplot(oracle_scna, aes(x=reorder(Gene, pmcc), y=pmcc, fill=sig)) + geom_bar(stat="identity") + xlab("ORACLE genes") + ylab("Correlation of gene expression\nwith copy-number status\n(Pearson r)") + theme_classic() + theme(legend.position = "none") + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 0.5)
}

################
##  Extended Data 9
################

##  ED-9a)
{
  ## plot 
  gg <- ggplot(observed_ith, aes(x=biopsy_number, y=ith_score)) + facet_wrap(~PublicationID, nrow=5) 
  # theme
  gg <- gg + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + theme(aspect.ratio=1) + theme(strip.background = element_rect(colour="white"))
  # labels
  gg <- gg + xlab("No. of Biopsies") + ylab("Gene Expression ITH") 
  # scales
  gg <- gg + scale_x_continuous(expand = c(0,0), limits = c(0, 7), breaks = seq(0,7,2)) + scale_y_continuous(expand = c(0,0), limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2))
  # add points
  gg <- gg + geom_point(pch=21, fill="black", stroke=0, alpha=0.5, size=3)
  # add mean
  gg <- gg + geom_line(data=output_mean, aes(x=biopsy_number, y=mean), col="firebrick1")
  # add lower sd
  gg <- gg + geom_line(data=output_sd, aes(x=biopsy_number, y=sd_lo), col="dodgerblue")
  # add upper sd
  gg <- gg + geom_line(data=output_sd, aes(x=biopsy_number, y=sd_up), col="dodgerblue")
  gg %>% print()
}

##  ED-9b)
{
  # prep for loop
  plots <- vector("list", (ncol(rnaith_immune)-2))
  # run for loop
  for (i in 1:(ncol(rnaith_immune)-2)) {
    db <- corr_plot_DB(data = rnaith_immune, x = colnames(rnaith_immune)[i+2], y="ITH", corr_method = "spearman", title = colnames(rnaith_immune)[i+2], x_axis = "Immune infiltrate", y_axis = "RNA ITH", best_fit_line = T, point_alpha = 0.75)
    plots[[i]] <- db$plot
  }
  # correlation plots: RNA-ITH ~ immune infiltrate
  cowplot::plot_grid(plotlist=plots)
}

##  ED-9c)
{
  # correlation plot: RNA-ITH ~ tumour purity
  corr_plot_DB(data = rnaith_purity, x="purity", y="ith", x_axis = "WES tumour purity", y_axis = "RNA ITH", title = "", best_fit_line = TRUE)
}

################
##  Extended Data 10
################

{
  #get all reactome pathways
  tmp <- read.table('ReactomePathways.20181008.gmt', sep='\t', as.is=T, fill=T)
  reactome_pathways <- list()
  for(i in 1:nrow(tmp)){
    reactome_pathways[[i]] <- as.character(tmp[i,-c(1:3)])
  }
  names(reactome_pathways) <- tmp[,2]
  
  entrez_id <- names(map2entrez(rownames(rnaith_dhruva)))
  rnaith_dhruva <- cbind(rnaith_dhruva, entrez_id)
  
  tmp <- setNames(c('quadrant', 'LUAD_quadrant','LUSC_quadrant'), c('All','LUAD','LUSC'))
  for(i in 1:3){
    cat(i)
    lower_right_q <- enrichPathway(gene=as.character(rnaith_dhruva[rnaith_dhruva[,tmp[i]] %in% 'bottom_right','entrez_id']),pvalueCutoff=0.05,  readable=T, pAdjustMethod='bonferroni', universe=unique(rnaith_dhruva[,'entrez_id']))
    upper_left_q <- enrichPathway(gene=as.character(rnaith_dhruva[rnaith_dhruva[,tmp[i]] %in% 'top_left','entrez_id']),pvalueCutoff=0.05, readable=T, pAdjustMethod='bonferroni', universe=unique(rnaith_dhruva[,'entrez_id']))
    upper_right_q <- enrichPathway(gene=as.character(rnaith_dhruva[rnaith_dhruva[,tmp[i]] %in% 'top_right','entrez_id']),pvalueCutoff=0.05, readable=T, pAdjustMethod='bonferroni', universe=unique(rnaith_dhruva[,'entrez_id']))
    lower_left_q <- enrichPathway(gene=as.character(rnaith_dhruva[rnaith_dhruva[,tmp[i]] %in% 'bottom_left','entrez_id']),pvalueCutoff=0.05, readable=T, pAdjustMethod='bonferroni', universe=unique(rnaith_dhruva[,'entrez_id']))
    
    par(las=1)
    barplot(lower_right_q, showCategory=nrow(lower_right_q)) # Q4
    barplot(upper_left_q, showCategory=nrow(upper_left_q)) # Q2
    barplot(upper_right_q, showCategory=nrow(upper_right_q)) # Q3
    barplot(lower_left_q, showCategory=nrow(lower_left_q)) # Q1
    
    dotplot(lower_right_q, showCategory=nrow(lower_right_q))
    dotplot(upper_left_q, showCategory=nrow(upper_left_q))
    dotplot(upper_right_q, showCategory=nrow(upper_right_q))
    dotplot(lower_left_q, showCategory=nrow(lower_left_q))
    
  }
  
}

