# Auto-generated R script from notebook BSBM_project_notebook2.ipynb
# Conversion date: 2025-06-05

# ---- markdown ----
# # Baltic Sea BioMed summer school: PROJECT WORK, DAY 2 
# This notebook contains analysis of RNA-seq data generated from human liver 
# cancer cell line HEPG2 (https://en.wikipedia.org/wiki/Hep_G2), that is commonly 
# used model for liver toxicity. Cells have been treated with 
# DMSO (control), PCB118 or PCB153 for 24 hours. PCB118 and PCB153 are 
# polychlorinated biphenyls (https://en.wikipedia.org/wiki/Polychlorinated_biphenyl). 
# RNA-seq libraries are produced from mRNAs, mapped to human genome hg38 and 
# differential gene expression (DGE) is analyzed with 
# DESeq2 (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
# and using DMSO as a reference.  
#
#   DAY 1: We'll get familiar with the data and produce various common visualizations.  
#   DAY 2: We'll do gene-set enrichment analysis (GSEA) and predict upstream transcription regulators using Lisa online tool. 

# ---- markdown ----
# ## 1. Load libraries and import RNA-seq data frame

# ---- code ----
## load libraries
# these are used in GSEA analysis
library(clusterProfiler)
library(enrichplot)

# these are used for plotting
library(ggplot2)
library(ggrepel)
library(gridExtra)

# ---- markdown ----
# Next, we'll use `readRDS()` to read-in the data from yesterday (remember that 
# we stored it as R object?). We also use `set.seed()` to set seed value for 
# random number generators. This is done to retain consistent results. You can choose your own seed number :)

# ---- code ----
## read RNA-seq data frame from day 1
# define home directory
home_dir <- getwd()

# read saved data frame as df and check that it's ok
df <- readRDS(paste0(home_dir, "/df.rds"))
head(df)

# set random generator seed for consistent results (you can choose your own)
set.seed(38642)

# TPM cutoff might be needed 
cutoff_tpm <- 1   # cutoff value for expressed genes

# ---- markdown ----
# ## 2. Gene-set enrichment analysis (GSEA)
# More comprehensive instructions and theoretical background can be found at 
# Biomedical knowledge mining book (https://yulab-smu.top/biomedical-knowledge-mining-book/index.html)

# ---- markdown ----
# **OPTIONAL STEP:** download Hallmark gene set from mSigDB 

# ---- code ----
## OPTIONAL: load pathways from mSigDB ##
#########################################
## you can browse mSigDB signature pathways from their webpage:
# https://www.gsea-msigdb.org/gsea/msigdb/index.jsp

## or using msigdbr tool:
library(msigdbr) # MSigDB - molecular signature data base
as.data.frame(msigdbr_collections()) # check pathway categories and subcategories. 

# load gene set names and gene symbols belonging to hallmark category
pathway <- msigdbr(species = "Homo sapiens", 
                    collection = "H", 
                    subcollection  = NULL)


# older version of msigdbr uses this format
#pathway <- msigdbr(species = "Homo sapiens", category = "H", subcategory = "")


# An other example to download reactome pathways
#pathway2 <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") # new way
# pathway2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") # older msigdbr versions

# this downloads table as a "tibble" format with various gene annotations
head(pathway)

# for GSEA, we only need columns for gene set name (gs_name) and associated gene symbols (gene_symbols)
pathway <- pathway[,c("gs_name", "gene_symbol")]
head(pathway)

# write pathway to file for later use
write.table(x = pathway, 
           file = paste0(home_dir,"/Hallmark_pathway_mSigDB.tsv"),
           quote = FALSE,
           row.names = FALSE,
           col.names = TRUE,
           sep = "\t"
           )

# ---- code ----
## GSEA of PCB118 treatment responses
# save DMSO vs PCB118 logFC values as a matrix and add gene symbols as rownames
ma2 <- as.matrix(df$DMSO.50nM_vs_PCB118.50microM_Log2_Fold_Change)
rownames(ma2) <- df$symbol
colnames(ma2) <- "PCB118"

# order according to logFC 
ma2 <- ma2[order(ma2, decreasing = TRUE),]

# check if there are NAs
table(is.na(ma2))

# remove NAs
ma2 <- na.omit(ma2)

# check first and last entries
head(ma2)
tail(ma2)

# upload hallmarks pathway (if not done already)
pathway <- read.delim(paste0(home_dir,"/Hallmark_pathway_mSigDB.tsv"))

## run GSEA analysis (set pvalueCutoff to 1 if you want to have all results)
gsea_results <- GSEA(geneList = ma2, 
                     TERM2GENE = pathway, 
                     pvalueCutoff = 0.05)

# you can look results as a data frame 
gsea_table <- as.data.frame(gsea_results)
gsea_table[, 1:(ncol(gsea_table)-2)] # last columns are a bit too big for screen..

# ---- code ----
# classical GSEA plots
gseaplot(gsea_results, by = "all", title = gsea_results$Description[1], geneSetID = 1)
#gseaplot(gsea_results, by = "all", title = gsea_results$Description[2], geneSetID = 2)

# ---- code ----
## dotplot (from enrichplot package)
dotplot(gsea_results, showCategory=15, split=".sign") + facet_grid(.~.sign)

# ---- code ----
## network plot for GSEA results
# select how many GSEA categories to plot
top_n_category <- 5

# plot network graph showing shared genes and fold change from matrix as color
#pdf(paste0(home_dir, "cnet_plot_of_top_", top_n_category ,"_GSEA_categories.pdf"))
cnetplot(gsea_results, 
         categorySize = "geneNum", 
         foldChange = ma2, 
         showCategory = top_n_category, 
         node_label = "category")
#dev.off()

# ---- code ----
## PCB153 GSEA analysis
# copy and modify the above code for this

# ---- markdown ----
# ## 3. Upstream regulator prediction using Lisa
# Next we'll use online tool Lisa (documentation: http://lisa.cistrome.org/doc), 
# to predict upstream transcription regulators. With this tool we can predict 
# which transcription factors and coregulators are responsible for inducing and 
# repressed genes in response to PCB118 or PCB153 treatments. 

# ---- markdown ----
# ### 3.1 Prepare gene lists 
# Lisa want's to have 50-500 genes to make reliable predictions. Our DEG categories 
# are all bigger than 500 genes. We'll order genes based on adjusted p-value and 
# take top 500 induced (up) or repressed (dn) genes for our gene lists.

# ---- code ----
## write table of induced gene symbols for LISA analysis: http://lisa.cistrome.org/

# sort according to adj. pvalue for PCB118 analysis
df <- df[order(df$DMSO.50nM_vs_PCB118.50microM_adj_pvalue),]

# create directory for data
dir.create(paste0(home_dir, "/data"))

# write table with top 500 PCB118-induced gene symbols for LISA analysis
write.table(x = head(df[df$PCB118_up == 1, ]$symbol, n = 500),
          file = paste0(home_dir, "/data/PCB118_up_top500.txt"),
          quote = FALSE, row.names = FALSE, col.names = FALSE)
          
# write table with top 500 PCB118-induced gene symbols for LISA analysis
write.table(x = head(df[df$PCB118_dn == 1, ]$symbol, n = 500),
          file = paste0(home_dir, "/data/PCB118_dn_top500.txt"),
          quote = FALSE, row.names = FALSE, col.names = FALSE)

# sort according to adj. pvalue for PCB153 analysis
df <- df[order(df$DMSO.50nM_vs_PCB153.50microM_adj_pvalue),]

# write table with top 500 PCB153-induced gene symbols for LISA analysis
write.table(head(df[df$PCB153_up == 1, ]$symbol, n = 500),
          file = paste0(home_dir, "/data/PCB153_up_top500.txt"),
          quote = FALSE, row.names = FALSE, col.names = FALSE)
          
# write table with top 500 PCB153-induced gene symbols for LISA analysis
write.table(x = head(df[df$PCB153_dn == 1, ]$symbol, n = 500),
          file = paste0(home_dir, "/data/PCB153_dn_top500.txt"),
          quote = FALSE, row.names = FALSE, col.names = FALSE)

# continue with these gene lists on Lisa online tool

# ---- markdown ----
# ### 3.2 run Lisa online tool
# Open Lisa (http://lisa.cistrome.org/) online tool. Select Species: Human 
# and insert gene symbols (copy-paste from text editor) of induced (up) genes to Gene set 1, 
# and repressed (dn) gene lists to Gene set 2. Add Job name: PCB118 or PCB153 and 
# optionally email address. Click Run to initiate analysis. If you added email 
# address, you will be notified when the run has been completed.  
# Once the runs are complete, press Dowload all samples, dowload and unzip the 
# results .zip file to home_dir/lisa_results/.  
# Analysis produces ~60 files, but we'll only use two that contain 
# de-duplicated transcription regulatory predictions based on ChIP-seq data for 
# induced (gene set 1 = gs1) and repressed (gs2) genes. 

# ---- markdown ----
# ### 3.3 Prepare figures from Lisa results
# Let's build a similar scatter-plot of lisa results as is shown on Lisa webpage 
# and label top-n most significant transcription factors (TFs)!  
#
# This time well make a vector of file names that contain 
# "chipseq_cauchy_combine_dedup" and use `lapply()` function ("list apply") 
# to read in result files as a list. Scatter plot images of -log(p-values) are 
# generated in `lapply()` loop and lables are added to top_n transcription 
# regulators with lowest p-value. Images are collected to a list called p3, 
# which can be accessed later on.

# ---- code ----
# Simple scatter plot
## define how many TFs you would like to label in plot
top_n <- 20

## upload Lisa results
# set home directory to point project folder with Lisa results
lisa_res_dir <- paste0(home_dir, "/lisa_results/")

# select files with results from ChIP-seq model
in_files <- list.files(lisa_res_dir, pattern = "chipseq_cauchy_combine_dedup", full.names = TRUE)
# in_files

# if you followed instructions while uploading data to Lisa, your...
# predictions for up genes have gs1 in file name, eg. PCB153_2024_07_12_0739250.215__hg38_gs1.txt_chipseq_cauchy_combine_dedup.csv
# predictions for dn genes have gs2 in file name, eg. PCB153_2024_07_12_0739250.215__hg38_gs2.txt_chipseq_cauchy_combine_dedup.csv

# read in data from the files as a list
lisa_results <- lapply(in_files, function(i){
  read.csv(i)
})

# name list items (check that the naming corresponds the order of your files!)
names(lisa_results) <- c("PCB118_up", "PCB118_dn", "PCB153_up", "PCB153_dn")

# i=1
# collect results from different treatments to a list
p3 <- lapply(1:(length(lisa_results)/2), function(i){
  ## predictions for repressed genes
  # i*2 is the dn group and i*2-1 is the up group
  dn_data <- lisa_results[[i*2]]
  # rename p-value column to contain "dn_" prefix
  names(dn_data)[2] <- paste0("dn_", names(dn_data)[2])
  
  ## predictions for induced genes
  up_data <- lisa_results[[i*2-1]]
  # rename p-value column to contain "up_" prefix
  names(up_data)[2] <- paste0("up_", names(up_data)[2])
  
  # check that both up and dn data contain same TFs before merging data frames 
  print(paste0("dn and up data contain the same TFs? ", identical(sort(dn_data$TF), sort(up_data$TF))))
  # not same TF data is kept for both though
  # identical(sort(dn_data$X), sort(up_data$X))
  
  # merge data frames for plotting
  data <- merge.data.frame(x = up_data, y = dn_data, by = "TF")
  
  # add negative log10 p-values  for plotting
  data$neg_log_pval_up <- -log10(data$up_pval)
  data$neg_log_pval_dn <- -log10(data$dn_pval)
  
  # add rank data
  data$rank_up <- rank(data$up_pval)
  data$rank_dn <- rank(data$dn_pval)
  data$is_top_n <- data$rank_up < top_n | data$rank_dn < top_n
  data$label <- ifelse(data$is_top_n, data$TF, "")
  
  # treatment name
  treatment_name <- gsub("_dn", "", names(lisa_results[i*2]))
  
  # plot 
  ggplot(data, aes(x=neg_log_pval_up, y = neg_log_pval_dn, label = label))+
    theme_classic()+
    geom_point(shape = 19, col=rgb(1,0,0,0.5))+
    geom_text_repel(max.overlaps = 50, size = 4, color = "black")+
    labs(title = paste0("Lisa predictions for ", treatment_name, " treatment (top ",top_n," TFs labeled)"),
         x = paste0(names(lisa_results[i*2 -1]), "  [-log10(p-value)]"),
         y = paste0(names(lisa_results[i*2]), " [-log10(pvalue)]"))
  
})

# images were colletecd to a list named p3.
# show first image
p3[1]

# show second image
p3[2]

## show/print both images in a grid
# uncomment "pdf()" and "dev.off()" lines to print pdf of the results
# pdf(paste0(home_dir, "/Lisa_results_simple_top", top_n, ".pdf"), height = 16, width = 8)
# do.call(grid.arrange, p3)
# dev.off()

# ---- markdown ----
# Next we'll generate a bit more detailed scatter plot that combines Lisa results 
# data with the RNA-seq information.  
# Here we'll combine two data frames using `merge.data.frame()` function. 
# Then we'll remove Lisa predictions that are not expressed in the RNA-seq 
# data (ie. do not meet TPM cutoff value), adjust the circle sizes to represent 
# median expression in RNA-seq, and add logFC from the RNA-seq data as a color. 

# ---- code ----
# a more detailed scatter plot with expression data from RNA-seq
## define how many TFs you would like to label in plot
top_n <- 20

## upload Lisa results
# set home directory to point project folder with Lisa results
lisa_res_dir <- paste0(home_dir, "/lisa_results/")

# select files with results from ChIP-seq model
in_files <- list.files(lisa_res_dir, pattern = "chipseq_cauchy_combine_dedup", full.names = TRUE)
# in_files

# if you followed instructions while uploading data to Lisa, your...
# predictions for up genes have gs1 in file name, eg. PCB153_2024_07_12_0739250.215__hg38_gs1.txt_chipseq_cauchy_combine_dedup.csv
# predictions for dn genes have gs2 in file name, eg. PCB153_2024_07_12_0739250.215__hg38_gs2.txt_chipseq_cauchy_combine_dedup.csv

# read in data from the files as a list
lisa_results <- lapply(in_files, function(i){
  read.csv(i)
})

# name list items (check that the naming corresponds the order of your files!)
names(lisa_results) <- c("PCB118_up", "PCB118_dn", "PCB153_up", "PCB153_dn")

# collect results from different treatments to a list
p4 <- lapply(1:(length(lisa_results)/2), function(i){
  ## predictions for repressed genes
  # i*2 is the dn group and i*2-1 is the up group
  dn_data <- lisa_results[[i*2]]
  # rename p-value column to contain "dn_" prefix
  names(dn_data)[2] <- paste0("dn_", names(dn_data)[2])
  
  ## predictions for induced genes
  up_data <- lisa_results[[i*2-1]]
  # rename p-value column to contain "up_" prefix
  names(up_data)[2] <- paste0("up_", names(up_data)[2])
  
  # check that both up and dn data contain same TFs before merging data frames 
  print(paste0("dn and up data contain the same TFs? ", identical(sort(dn_data$TF), sort(up_data$TF))))
  # not same TF data is kept for both though
  # identical(sort(dn_data$X), sort(up_data$X))
  
  # merge data frames for plotting
  data <- merge.data.frame(x = up_data, y = dn_data, by = "TF")
  
  # merge expression data from df
  data <- merge.data.frame(x = data, y = df, by.x = "TF", by.y = "symbol")
  
  # remove non-expressed
  data <- data[data$is_exp,]
  
  # add negative log10 p-values  for plotting
  data$neg_log_pval_up <- -log10(data$up_pval)
  data$neg_log_pval_dn <- -log10(data$dn_pval)
  
  # add rank data
  data$rank_up <- rank(data$up_pval)
  data$rank_dn <- rank(data$dn_pval)
  data$is_top_n <- data$rank_up < top_n | data$rank_dn < top_n
  data$label <- ifelse(data$is_top_n, data$TF, "")
  
  # treatment name
  treatment_name <- gsub("_dn", "", names(lisa_results[i*2]))
  
  # factorize expression values for size plotting
  data$tpm_fact <- cut(data$median_exp, 
                       right = FALSE, # right border open-ended: [1,2)
                       breaks = c(0, 2.5, 5, 10, 20, 30, 1000), 
                       labels = c("<2.5", "<5", "<10", "<20", "<30","30-"),
                       ordered = TRUE)
  # grep fold change column name for later use
  col_name_logfc <- grep("Log2_Fold_Change", grep(treatment_name, names(df), value = T), value = T)
  
  # plot 
  ggplot(data, aes(x=neg_log_pval_up, 
                   y = neg_log_pval_dn, 
                   label = label, 
                   size = tpm_fact,
                   # col = grep("Log2_Fold_Change", grep(treatment_name, names(df), value = T), value = T)))+
                   # col = col_name_logfc))+
                   col = .data[[col_name_logfc]]))+
    theme_classic()+
    geom_point()+
    scale_color_gradientn(colors = rev(c("#A50026", "#F46D43", "#FDAE61", "#FFFFBF", "#E0F3F8",  "#74ADD1", "#313695")),
                          values = scales::rescale(c(min(data[,col_name_logfc], na.rm = TRUE), -1, -0.5, 0, 0.5 , 1, max(data[,col_name_logfc], na.rm = TRUE))),
                          na.value = "gray")+
    # 
    # scale_size_discrete(range = c(2,6), na.value = 1)+
    # scale_size_continuous(trans = "log10", range = c(1, 10)) + 
    geom_text_repel(max.overlaps = 50, size = 4, color = "black")+
    labs(title = paste0("Lisa predictions for ", treatment_name, " treatment (top ",top_n," expressed TFs labeled). TPM > ", cutoff_tpm, "."),
         x = paste0(names(lisa_results[i*2 -1]), "  [-log10(p-value)]"),
         y = paste0(names(lisa_results[i*2]), " [-log10(pvalue)]"),
         col = paste0(treatment_name, " vs DMSO\n(Log2FC)"),
         size = "Expression\n(median TPM)")
  
})

# show 1st image
p4[1]

# show 2nd image
p4[2]

## show ggplot images in grid
# uncomment "pdf()" and "dev.off()" lines to print pdf of the results
# pdf(paste0(home_dir, "/Lisa_results_pretty_top", top_n, ".pdf"), height = 16, width = 8)
# do.call(grid.arrange, p4)
# dev.off()

# ---- code ----
names(df)

