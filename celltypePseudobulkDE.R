library("DESeq2")
filepath <- "parentDir/cellType_totalExpression_genoWT_DATE.csv"

rsem.in <- as.matrix(read.csv(filepath, row.names=1))
head(rsem.in)
colData = data.frame(row.names=colnames(rsem.in),model=c(rep(c("GENE_geno"), 4), rep(c("GENE_WT"), 4)),sample=colnames(rsem.in))
colData$model = relevel(colData$model, "GENE_WT")

dds <- DESeqDataSetFromMatrix(countData = rsem.in,
                              colData = colData,
                              design = ~model)
dds
dds <- DESeq(dds)

res <- results(dds)
head(res,6)
summary(dds)

res <- lfcShrink(dds=dds, coef=2, res=res, type="normal")
res <- res[order(res$pvalue),]
head(res, 6)

pdf(file="GENE_PROJECT_cellType_genoWT_volcano_DATE.pdf")              
high <- row.names(res)[1:10]
high
                         
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="GENE_PROJECT_cellType_genoWT"))
with(subset(res, padj<0.05), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, row.names(res) %in% high), text(log2FoldChange, -log10(pvalue), labels=high, cex=0.7, font=1, col="red"))
dev.off()                   

res["GENE",]
write.csv(res, "GENE_PROJECT_cellType_genoWT-DESeq2_DATE.csv")