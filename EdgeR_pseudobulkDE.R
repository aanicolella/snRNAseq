library(edgeR)
library(fdrtool)

###### DE BY GENOTYPE ######
#read in pseudobulk sheets
filepath <- "/path/to/pseudobulkTable.csv"
rsem.in <- as.matrix(read.csv(filepath, row.names=1))
head(rsem.in)

# Set up metadata
# model = genotype - adjust numbers and order to match your data
# sample = individual sample name
colData = data.frame(model=c(rep(c("KO"), 5), rep(c("HT"), 5), rep(c("WT"), 5)),sample=colnames(rsem.in))
colData$model = as.factor(colData$model)
colData$model = relevel(colData$model, "WT")
# Add batch number if data comes from multiple experiments, example below
#colData$batch = as.factor(c(1,1,1,1,2,1,1,1,1,2,1,1,2,2,2))

dge = DGEList(counts = rsem.in, group = colData$model)

#filter out lowly expressed genes
keep = filterByExpr(dge)
dge = dge[keep, , keep.lib.sizes = FALSE]

#normalize library sizes
dge = calcNormFactors(dge)

#### Generalized linear models (GLM) for more complex designs
#create model matrix, can add covariates if desired
design = model.matrix(~model, data = colData)
# Design matrix when multiple batches are used
#design = model.matrix(~model + batch, data = colData)
#estimate dispersions
dge = estimateDisp(dge, design)

### Likelihood ratio test (LRT) - more power, less control of false positives
fit = glmFit(dge, design)
#HT vs WT
lrtHT = glmLRT(fit, coef = 2)
#KO vs WT
lrtKO = glmLRT(fit, coef = 3)

##explore the results
topTags(lrtHT)
topTags(lrtKO)
#see how many FDR < 0.05
summary(decideTests(lrtHT))
summary(decideTests(lrtKO))

# Create dropped gene df--dropped genes are added back as NA rows to allow for easier merging of results from multiple experiments for comparison
rem <- as.data.frame(keep[keep==FALSE])
dropRes <- data.frame(row.names=row.names(rem), matrix(ncol=5, nrow=length(row.names(rem))))
cols = c("logFC","logCPM","LR","PValue","FDR")
colnames(dropRes) = cols
#output results
resHT = lrtHT$table
resHT$FDR = p.adjust(resHT$PValue, "fdr")
# Optional recalculation of nominal and adjusted p-values using fdrtool
#norm_HT=-sign(resHT$logFC)*qnorm(resHT$PValue/2)
#pval_HT=fdrtool(norm_HT,plot=F)$pval
#resHT["pval_fdrtools"]=pval_HT
#resHT["padj_fdrtools"]=p.adjust(pval_HT,"fdr")
resHT_merge = rbind(resHT,dropRes)
dir.create("HT",showWarnings=FALSE)
write.csv(resHT_merge, "HT/PROJECT_HTWT_EdgeR_DATE.csv")

resKO = lrtKO$table
resKO$FDR = p.adjust(resKO$PValue, "fdr")
# Optional recalculation of nominal and adjusted p-values using fdrtool
#norm_KO=-sign(resKO$logFC)*qnorm(resKO$PValue/2)
#pval_KO=fdrtool(norm_KO,plot=F)$pval
#resKO["pval_fdrtools"]=pval_KO
#resKO["padj_fdrtools"]=p.adjust(pval_KO,"fdr")
resKO_merge = rbind(resKO,dropRes)
dir.create("KO",showWarnings=FALSE)
write.csv(resKO_merge, "KO/PROJECT_KOWT_EdgeR_DATE.csv")

##### DE within cell types ######
#input desired cell types 
celltypes <- c()
for (i in celltypes){
  #load data from pseudobulk sheets
  fp <- paste0("/path/to/celltype/pseudobulkSheets/", i, "_totalExpression_DATE.csv")
  rsem.in <- as.matrix(read.csv(fp, row.names = 1))
  
  #make sure this matches your data
  colData = data.frame(model=c(rep(c("KO"), 5), rep(c("HT"), 5), rep(c("WT"), 5)),sample=colnames(rsem.in))
  colData$model = as.factor(colData$model)
  colData$model = relevel(colData$model, "WT")
  # Add batch number if data comes from multiple experiments, example below
  #colData$batch = as.factor(c(1,1,1,1,2,1,1,1,1,2,1,1,2,2,2))
  
  dge = DGEList(counts = rsem.in, group = colData$model)
  
  keep = filterByExpr(dge)
  dge = dge[keep, , keep.lib.sizes = FALSE]
  
  dge = calcNormFactors(dge)
  
  ### GLM (LRT)
  design = model.matrix(~model, data = colData)
  # Design matrix when multiple batches are used
  #design = model.matrix(~model + batch, data = colData)
  dge = estimateDisp(dge, design)
  
  ##LRT
  fit = glmFit(dge, design)
  lrtHT = glmLRT(fit, coef = 2)
  lrtKO = glmLRT(fit, coef = 3)
  
  rem <- as.data.frame(keep[keep==FALSE])
  dropRes <- data.frame(row.names=row.names(rem), matrix(ncol=5, nrow=length(row.names(rem))))
  cols = c("logFC","logCPM","LR","PValue","FDR")
  colnames(dropRes) = cols
  
  resHT = lrtHT$table
  resHT$FDR = p.adjust(resHT$PValue, "fdr")
  # Optional recalculation of nominal and adjusted p-values using fdrtool
  #norm_HT=-sign(resHT$logFC)*qnorm(resHT$PValue/2)
  #pval_HT=fdrtool(norm_HT,plot=F)$pval
  #resHT["pval_fdrtools"]=pval_HT
  #resHT["padj_fdrtools"]=p.adjust(pval_HT,"fdr")
  resHT_merge = rbind(resHT,dropRes)
  dir.create("HT/byCellType",showWarnings=FALSE)
  dir.create(paste0("HT/byCellType/", i),showWarnings=FALSE)
  write.csv(resHT_merge, paste0("HT/byCellType/",i,"/",i,"_HTWT_EdgeR_DATE.csv"))
  
  resKO = lrtKO$table
  resKO$FDR = p.adjust(resKO$PValue, "fdr")
  # Optional recalculation of nominal and adjusted p-values using fdrtool
  #norm_KO=-sign(resKO$logFC)*qnorm(resKO$PValue/2)
  #pval_KO=fdrtool(norm_KO,plot=F)$pval
  #resKO["pval_fdrtools"]=pval_KO
  #resKO["padj_fdrtools"]=p.adjust(pval_KO,"fdr")
  resKO_merge = rbind(resKO,dropRes)
  dir.create("KO/byCellType",showWarnings=FALSE)
  dir.create(paste0("KO/byCellType/", i),showWarnings=FALSE)
  write.csv(resKO_merge, paste0("KO/byCellType/",i,"/",i,"_KOWT_EdgeR_DATE.csv"))
  
}