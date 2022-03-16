library(edgeR)
##### DE within cell types ######
#input desired cell subtypes 
celltypes <- c()
for (i in celltypes){
  #load data from pseudobulk sheets
  fp <- paste0("/path/to/subtype/pseudobulkSheets/", i, "_totalExpression_DATE.csv")
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
  dir.create("HT/byCellType/targetCellType/targetSubtypes",showWarnings=FALSE)
  write.csv(resHT_merge, paste0("HT/byCellType/targetCellType/targetSubtypes/",i,"_HTWT_EdgeR_DATE.csv"))
  
  resKO = lrtKO$table
  resKO$FDR = p.adjust(resKO$PValue, "fdr")
  # Optional recalculation of nominal and adjusted p-values using fdrtool
  #norm_KO=-sign(resKO$logFC)*qnorm(resKO$PValue/2)
  #pval_KO=fdrtool(norm_KO,plot=F)$pval
  #resKO["pval_fdrtools"]=pval_KO
  #resKO["padj_fdrtools"]=p.adjust(pval_KO,"fdr")
  resKO_merge = rbind(resKO,dropRes)
  dir.create("KO/byCellType/targetCellType/targetSubtypes/",showWarnings=FALSE)
  write.csv(resKO_merge, paste0("KO/byCellType/targetCellType/targetSubtypes/",i,"_KOWT_EdgeR_DATE.csv"))
  
}