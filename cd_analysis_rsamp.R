setwd("~/Documents/r_saves/R_sessions/final/data")
library(dplyr)

# Read the text file into a data frame
file_path <- "~/Documents/r_saves/R_sessions/final/data/GSE67250_RnaSeqCounts.txt" # Replace with your file path
data <- read.table(file_path, header = TRUE)

# Create the seqdata data frame
seqdata <- data

# Extract the sample names from the column names (excluding the GeneID column)
sample_names <- colnames(data)[-1]

# Create the sampleInfo data frame
sampleInfo <- tibble(
  Samples = sample_names,
  Condition = ifelse(grepl("^NL", sample_names), "control", "crohn's disease")
)

# Display the seqdata and sampleInfo data frames
sampleInfo <- as.data.frame(sampleInfo)


#bar plots of library sizes:
library(ggplot2)
library(dplyr)

countdata <- seqdata[,-1]
rownames(countdata) <- seqdata[,1]
COl_counts <- colSums(countdata)

data_for_counts_bar <- cbind(sampleInfo, COl_counts)
counts_bargraph <- barplot(height = data_for_counts_bar$COl_counts, names.arg = data_for_counts_bar$Samples, ylab = "Library Size", xlab = "Sample")
# NO SAMPLES HAVE LIBRARY SIZES OF UNDER 2 MILLION, CAN KEEP ALL

countdata <- rbind(countdata, COl_counts)

# FILTERING LOW COUNT GENES:
library("edgeR")
myCPM <- cpm(countdata)
# make a data frame with CPM values added:
my.df <- data.frame()
for (sample.name in colnames(countdata)) {
  mat <- cbind(CPM = myCPM[,sample.name], Counts = countdata[,sample.name], Sample = rep(sample.name, nrow(countdata)))
  my.df <- rbind(my.df, mat)
}
my.df$CPM <- as.numeric(my.df$CPM)
freq_vs_cpm <- ggplot(data=my.df, aes(x=CPM)) +
  geom_histogram(binwidth=0.5) +
  labs(title = "Distribution of CPM",x="CPM",y="Frequency") +
  theme_minimal() +
  xlim(0,600) + ylim(0,100)
freq_vs_cpm

subset.df <- my.df[my.df$Sample == "NL2_RNASEQ",]
subset.df$CPM <- as.numeric(subset.df$CPM)
subset.df$Counts <- as.numeric(subset.df$Counts)
scat1 <- ggplot(data=subset.df,aes(x=CPM,y=Counts)) + geom_point() + xlim(0,2.5) + ylim(0,100)
scat1
# CPM thresh = 0.125
AboveThresh <- (myCPM > 0.125)
num.samples.aboveThresh <- rowSums(AboveThresh)
my.df2 <- as.data.frame(num.samples.aboveThresh)
aboveThresh_histo <- ggplot(data=my.df2,aes(x=num.samples.aboveThresh)) + geom_histogram()
aboveThresh_histo

# filter the genes that are not above the threshold in all 6 samples (10654 genes):
length(which(num.samples.aboveThresh < 6))
keep.genes <- names(which(num.samples.aboveThresh >= 6))
counts_filtered <- countdata[keep.genes,]

# Make filtered data a DGE object:
dgeObj <- DGEList(counts_filtered)


# QUALITY CONTROL:
log.cpm <- cpm(dgeObj, log=TRUE)
log.cpm.df <- data.frame()
for (mysample in colnames(counts_filtered)) {
  Mat <- cbind(Sample = rep(mysample,nrow(log.cpm)), Log.CPM = log.cpm[,mysample])
  log.cpm.df <- rbind(log.cpm.df, Mat)
}
log.cpm.df$Sample <- as.factor(log.cpm.df$Sample)
log.cpm.df$Log.CPM <- as.numeric(log.cpm.df$Log.CPM)
cpm_boxplot <- ggplot(data=log.cpm.df,aes(x=Sample,y=Log.CPM)) + geom_boxplot() +
  labs(x="Samples",y="Log CPM",title="Distribution of Log CPM across samples")
cpm_boxplot
# Conclusion: No real outliers or strange samples

sampleInfo$Condition <- as.factor(sampleInfo$Condition)
col.condition <- c("purple","cyan")[sampleInfo$Condition]
data.frame(sampleInfo$Condition,col.condition)
plotMDS(dgeObj,col=col.condition)
legend("topleft",fill=c("purple","cyan"), legend = levels(sampleInfo$Condition))
title("Condition Type")


# NORMALIZATION:
dgeObj <- calcNormFactors(dgeObj)
par(mfrow=c(1,2))
plotMD(log.cpm,column=2)
abline(h=0,col="blue")
plotMD(log.cpm,column=4)
abline(h=0,col="blue")

logcounts <- cpm(dgeObj,log=TRUE)
plotMD(logcounts,column=2)
abline(h=0,col="red")
plotMD(logcounts,column=4)
abline(h=0,col="red")
# conclusion: minor shift/change after normalization

# DIFF EXPRESSION:
mat_design <- as.formula(~ Condition)
model_matrix <- model.matrix(mat_design, data=sampleInfo)

dgeObj <- estimateCommonDisp(dgeObj)
dgeObj <- estimateGLMTrendedDisp(dgeObj)
dgeObj <- estimateTagwiseDisp(dgeObj)
# NOTE: ADD OTHERS HERE?
plotBCV(dgeObj)

fit <- glmFit(y=dgeObj,design = model_matrix)
lrt.CvsC <- glmLRT(fit, coef = 2)
topTags(lrt.CvsC)
# lowest false discovery rate is 40% (all of them are around here), none are significant at FDR <= 0.05

diff_ex_results <- as.data.frame(topTags(lrt.CvsC,n=Inf))
de <- decideTestsDGE(lrt.CvsC,adjust.method = "none",p.value=0.05, lfc=0)
# 560 downregulated, 379 upreg
# Log FC vs log CPM, significant in red
de2 <- cbind(de, diffExp=as.logical(de))
detags <- rownames(dgeObj)[as.logical(de)]
plotSmear(lrt.CvsC,de.tags = detags)
# decreasing p value to 0.01:
de <- decideTestsDGE(lrt.CvsC,adjust.method = "none",p.value=0.01, lfc=0)
detags <- rownames(dgeObj)[as.logical(de)]
plotSmear(lrt.CvsC,de.tags = detags)

# annotating with info from org.Hs.eg.db
library(org.Hs.eg.db)
all.genes <- rownames(diff_ex_results)

status <- all.genes %in% keys(org.Hs.eg.db, keytype="SYMBOL")
annotations <- select(org.Hs.eg.db,keys=rownames(diff_ex_results),keytype="SYMBOL",columns = c("GENENAME","ENTREZID","SYMBOL"))
annotations <- annotations[-12443,]
anno_results <- cbind(diff_ex_results,annotations)
head(anno_results)

#volcano plot of logP vs log FC
negLogPValue <- -log10(anno_results$PValue)
absLogFC <- abs(anno_results$logFC)
thresh <- -log10(0.01)
Significant <- (negLogPValue >= thresh & absLogFC > 3)
anno_results <- cbind(diff_ex_results,negLogPValue,Significant,annotations)
anno_results$logFC <- as.numeric(anno_results$logFC)
anno_results$Significant <- as.factor(anno_results$Significant)
volc <- ggplot(data=anno_results,aes(logFC,negLogPValue)) +
  geom_point(color=1+Significant) +
  theme_light() +
  labs(x="Log Fold Change",y="-log10(P value)")
volc


# HEATMAP CREATION:
library(pheatmap)
de_genes_expr <- logcounts[rownames(logcounts) %in% detags,]
# Scale the expression values by row (gene)
scaled_expr <- t(scale(t(de_genes_expr)))
# Create the color map: blue for downregulated, red for upregulated
color_map <- colorRampPalette(c("blue", "white", "red"))(100)
# Create the heatmap
# png("heatmap.png", width = 800, height = 800)

# Modify the sampleInfo data frame
sampleInfo_mod <- sampleInfo[, "Condition", drop = FALSE] # Keep only the "Condition" column
rownames(sampleInfo_mod) <- sampleInfo$Samples # Set row names to the sample names

pheatmap(scaled_expr,
         color = color_map,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = sampleInfo_mod)
# dev.off()
# browseURL("heatmap.png")

# downregulated only:

downregulated_genes <- rownames(diff_ex_results[diff_ex_results$logFC < 0, ])
downregulated_expr <- scaled_expr[rownames(scaled_expr) %in% downregulated_genes, ]
color_map <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap(downregulated_expr,
         color = color_map,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = sampleInfo_mod)

# Sort the results by log fold change (logFC) in descending order
sorted_results <- diff_ex_results[order(diff_ex_results$logFC, decreasing = TRUE),]

# Extract the top 10 most upregulated genes
top10_upregulated <- sorted_results[1:10,]

# Sort the results by log fold change (logFC) in ascending order
sorted_results <- diff_ex_results[order(diff_ex_results$logFC),]

# Extract the top 10 most downregulated genes
top10_downregulated <- sorted_results[1:10,]

# Combine the top 10 most upregulated and downregulated genes
top10_genes <- rbind(top10_upregulated, top10_downregulated)

# Print the results
print(top10_genes)

# NEW HEATMAPS BUT WITH ONLY TOP 20 THIS TIME:
# Extract the top 20 most upregulated genes
top20_upregulated <- sorted_results[1:20,]
# Extract the top 20 most downregulated genes
top20_downregulated <- sorted_results[(nrow(sorted_results)-19):nrow(sorted_results),]
# Get the expression values for the top 20 most upregulated genes
upregulated_genes_expr <- logcounts[rownames(logcounts) %in% rownames(top20_upregulated),]
# Get the expression values for the top 20 most downregulated genes
downregulated_genes_expr <- logcounts[rownames(logcounts) %in% rownames(top20_downregulated),]
# Scale the expression values by row (gene)
scaled_upregulated_expr <- t(scale(t(upregulated_genes_expr)))
scaled_downregulated_expr <- t(scale(t(downregulated_genes_expr)))

# Create the heatmaps
color_map <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap(scaled_upregulated_expr,
         color = color_map,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = sampleInfo_mod)

pheatmap(scaled_downregulated_expr,
         color = color_map,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = sampleInfo_mod)


# ATTEMPT WIKIPATHWAYS ANALYSIS:
# Create a new data frame with rownames as gene symbols
annotations_new <- annotations
rownames(annotations_new) <- annotations$SYMBOL

# Get the Entrez IDs for the differentially expressed genes
de_symbols <- detags
de_annotations <- annotations_new[de_symbols, ]
de_entrez <- de_annotations$ENTREZID

library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(DOSE)

data(de_entrez, package="DOSE")
enrich_result <- enrichWP(de_entrez, organism = "Homo sapiens") 
dotplot(enrich_result)

theirDownreg <- c("PTGIS","ATP10A","GDF10","DOK5","CGREF1","LOC388630","NAP1L2","LINC00116","FBLN1","FAM86DP","SPOCD1","FMO3","RASCRP1","ARHGEF4",
                  "SORBS1","KCNH2","AOC3","BMPER","IL33","IQSEC3","HSPB6","STEAP1","LOC100506409","CLIC3","COL5A3","CNTNAP3","EFNA5","CDH3","WNT2B","ARHGAP42")
commonDownreg <- intersect(theirDownreg, row.names(top20_upregulated))

theirUpreg <- c("SNX10","OBSCN","CADM1","LINC00256B","PREX2","NOX4","TSPAN11","CCNJL","ITGA8","TLR4","PLXDC2","KALRN","IGF1","PAG1","FAM117B","RORB","ATP6VOD2",
                "TANC1","ARL4C","VASH2","RAPGEF5","CLCN4","RAB38","ADAMTS5","SRD5A1","COL24A1","CARD16","C7orf69","ALCAM","DCDC2","GRAMD1C","AR","C10orf10","SCN3A",
                "CLIC6","MAFB","TMEM26","ANO4","PLCL1","IFI27")
dommonUpreg <- intersect(theirUpreg,row.names(top20_downregulated))
