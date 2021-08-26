library("DESeq2")

# Variables to replace
    # EXP: the name of the dataset you wish to analyze (likely age_brainregion)
    # DATE: date of the analysis
    # X, Y, Z: number of replicates for HT, KO, WT

# Load EXP count sheet for EXP
filepath <- "parentDir/totalExpression_EXP_DATE.csv"
rsem.in <- as.matrix(read.csv(filepath, row.names=1))
head(rsem.in)

# Set up metadata
    # model = genotype 
    # sample = individual sample name
colData = data.frame(model=c(rep(c("EXP_HT"), X), rep(c("EXP_KO"), Y), rep(c("EXP_WT"), Z)),sample=colnames(rsem.in))
colData$model = relevel(colData$model, "EXP_WT")

# Optional: find genes below expression threshhold. Create 'NA' data frame to re-add lowly expressed 
# genes for easier comparison between analyses.
drop <- rsem.in
drop$abundance <-
    drop$abundance[apply(drop$length,
                             1,
                             function(row) !all(row > 5)),]                              
drop$counts <-
  drop$counts[apply(drop$length,
                            1,
                             function(row) !all(row > 5)),]
drop$length <-
  drop$length[apply(drop$length,
                             1,
                             function(row) !all(row > 5)),]
                    
dropped <- drop$counts
dropRes <- data.frame(row.names=make.unique(row.names(drop$counts)), matrix(ncol=6, nrow=dim(drop$counts)[1]))
cols = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
colnames(dropRes) <- cols
EXPName1 <- data.frame(do.call('rbind', strsplit(as.character(row.names(dropRes)), '_', fixed=TRUE)))
EXPName1 <- EXPName1["X2"]               
dropRes$symbol <- EXPName1$X2

# Filter lowly expressed genes found above out of the main dataset for the DE analysis
rsem.in$abundance <-
  rsem.in$abundance[apply(rsem.in$length,
                             1,
                             function(row) all(row > 5 )),]
rsem.in$counts <-
  rsem.in$counts[apply(rsem.in$length,
                             1,
                             function(row) all(row > 5 )),]
rsem.in$length <-
  rsem.in$length[apply(rsem.in$length,
                             1,
                             function(row) all(row > 5 )),]

# Convert data structure to proper format and run DESeq
dds <- DESeqDataSetFromMatrix(countData = rsem.in,
                              colData = colData,
                              design = ~model)
dds <- DESeq(dds)
check <- estimateSizeFactors(dds)
norm <- counts(check, normalized=TRUE)
                       
# Combine normalized count dataframe with the 'NA' dataframe and save
all <- rbind(norm, dropped)                                          
EXPName <- data.frame(do.call('rbind', strsplit(as.character(row.names(all)), '_', fixed=TRUE)))
EXPName <- EXPName["X2"]               
row.names(all) <- EXPName$X2
all <- all[order(row.names(all)),]
head(all)
write.csv(all, "EXP_normalizedCounts.csv", row.names=TRUE)
                       
# Plot pca using different metadata conditions to see relationships between samples
pdf(file="EXP_pca.pdf")
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="sample")
dev.off()
pdf(file="EXP_modelpca.pdf")
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="model")
dev.off()

# Extract DE results from different genotype comparisions
    # Currently set up such that the LFC directionality is in relation to the changes in the first 
    # genotype listed: i.e. + lfc == upregulation in that geno, - lfc == downregulation
    # KOHT relationship is optional, sometimes interesting, primarily for QC purposes
resHT <- results(dds, contrast=c("model", "EXP_HT","EXP_WT"))
resKO <- results(dds, contrast=c("model", "EXP_KO","EXP_WT"))
resKOHT <- results(dds, contrast=c("model", "EXP_KO","EXP_HT"))
head(results(dds, tidy=TRUE), n=6)
summary(dds)

# Combine unscaled HTvsWT results and 'NA' dropped gene dataframe then write to .csv
EXPNameHT <- data.frame(do.call('rbind', strsplit(as.character(row.names(resHT)), '_', fixed=TRUE)))
EXPNameHT <- EXPNameHT["X2"]               
resHT$symbol <- EXPNameHT$X2
resMergeHT <- rbind(resHT, dropRes)
dim(resMergeHT)
resMergeHT <- resMergeHT[order(resMergeHT$pvalue),]
write.csv(resMergeHT, "EXP_HTWT-DESeq2_allgenes.csv")
# Scale HTvsWT results using lfcShrink(): https://rdrr.io/bioc/DESeq2/man/lfcShrink.html
   # Combine post-shrinkage results with 'NA' dropped dataframe then write to .csv
resHT <- lfcShrink(dds=dds, contrast=c("model", "EXP_HT", "EXP_WT"), res=resHT, type="normal")
EXPNameHT <- data.frame(do.call('rbind', strsplit(as.character(row.names(resHT)), '_', fixed=TRUE)))
EXPNameHT <- EXPNameHT["X2"]               
resHT$symbol <- EXPNameHT$X2
resMergeHT <- rbind(resHT, dropRes)
dim(resMergeHT)
resMergeHT <- resMergeHT[order(resMergeHT$pvalue),]
write.csv(resMergeHT, "EXP_HTWT-DESeq2_LFCshrink.csv")

# Repeat above for KOvsWT
EXPNameKO <- data.frame(do.call('rbind', strsplit(as.character(row.names(resKO)), '_', fixed=TRUE)))
EXPNameKO <- EXPNameKO["X2"]               
resKO$symbol <- EXPNameKO$X2
resMergeKO <- rbind(resKO, dropRes)
dim(resMergeKO)
resMergeKO <- resMergeKO[order(resMergeKO$pvalue),]
write.csv(resMergeKO, "EXP_KOWT-DESeq2_allgenes.csv")
resKO <- lfcShrink(dds=dds, contrast=c("model", "EXP_KO", "EXP_WT"), res=resKO, type="normal")
EXPNameKO <- data.frame(do.call('rbind', strsplit(as.character(row.names(resKO)), '_', fixed=TRUE)))
EXPNameKO <- EXPNameKO["X2"]               
resKO$symbol <- EXPNameKO$X2
resMergeKO <- rbind(resKO, dropRes)
dim(resMergeKO)
resMergeKO <- resMergeKO[order(resMergeKO$pvalue),]
write.csv(resMergeKO, "EXP_KOWT-DESeq2_LFCshrink.csv")

# Repeat above for KOvsHT
EXPNameKOHT <- data.frame(do.call('rbind', strsplit(as.character(row.names(resKOHT)), '_', fixed=TRUE)))
EXPNameKOHT <- EXPNameKOHT["X2"]               
resKOHT$symbol <- EXPNameKOHT$X2
resMergeKOHT <- rbind(resKOHT, dropRes)
dim(resMergeKOHT)
resMergeKOHT <- resMergeKOHT[order(resMergeKOHT$pvalue),]
write.csv(resMergeKOHT, "EXP_KOHTWT-DESeq2_allgenes.csv")
resKOHT <- lfcShrink(dds=dds, contrast=c("model", "EXP_KO", "EXP_HT"), res=resKOHT, type="normal")
EXPNameKOHT <- data.frame(do.call('rbind', strsplit(as.character(row.names(resKOHT)), '_', fixed=TRUE)))
EXPNameKOHT <- EXPNameKOHT["X2"]               
resKOHT$symbol <- EXPNameKOHT$X2
resMergeKOHT <- rbind(resKOHT, dropRes)
dim(resMergeKOHT)
resMergeKOHT <- resMergeKOHT[order(resMergeKOHT$pvalue),]
write.csv(resMergeKOHT, "EXP_KOHT-DESeq2_LFCshrink.csv")
