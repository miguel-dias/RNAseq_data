# Make a basic volcano plot
with(tt$table, plot(logFC, -log10(PValue), pch=20, main="Volcano_plot_edgeR", xlim=c(-10,10)))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano_plot_DESeq2", xlim=c(-10,10)))
# Add colored points: blue if -log10(padj)>10 and log2FC>2n if both)
#with(subset(res, padj<.000000000000001 ), points(log2FoldChange, -log10(padj), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(tt$table, PValue<0.0000000001 & abs(logFC)>=1), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(res, pvalue<0.0000000001 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#library(calibrate) #to identify genes in the volcano plot
#with(resfinal, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano_plot_DESeq2", xlim=c(-10,10)))
#with(subset(resfinal, pvalue<0.00000000000000000000000000000000001 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
#with(subset(resfinal, pvalue<0.00000000000000000000000000000000001 & abs(log2FoldChange)>=1), textxy(log2FoldChange, -log10(pvalue), labs=gene, cex=.8))


#                     VOLCANO PLOT OF EDGER TABLE(TT$TABLE)
## Draw a volcano plot
plot(tt$table$logFC, xlab="logFC",
     -log10(tt$table$PValue), ylab="-log10(PValue)",
     main = "Volcano plot",
     col=densCols(tt$table$logFC, -log10(tt$table$PValue)),
     xlim=c(-10,10),
     pch=20)

grid(col="grey", lty="solid")
## Line marking the null effect
abline(v=0, col="black") 
## Lines marking arbitrary effect sizes of 1 and -1 respectively, 
## which correspond to fold-changes of 2 and 1/2 in the raw microarray 
## intensities. 
abline(v=c(-1,1), col="red") 
## Line marking an arbitrary threshold of 1% on the p-value
abline(h=-log10(0.01), col="red")
## Draw a legend with the number of genes declared positive
legend("topleft", 
       legend=c(paste("Pvalue<=", 0.01, ":", sum(tt$table$PValue <= 0.01 & abs(tt$table$logFC)>=1))))

#                     VOLCANO PLOT OF DESEQ2 TABLE (res)
plot(res$log2FoldChange, xlab="logFC",
     -log10(res$pvalue), ylab="-log10(PValue)",
     main = "Volcano plot",
     col=densCols(res$log2FoldChange, -log10(res$pvalue)),
     xlim=c(-11,11), ylim=c(0,91),
     pch=20)
grid(col="grey", lty="solid")
## Line marking the null effect
abline(v=0, col="black") 
## Lines marking arbitrary effect sizes of 1 and -1 respectively, 
## which correspond to fold-changes of 2 and 1/2 in the raw microarray 
## intensities. 
abline(v=c(-1,1), col="red") 
## Line marking an arbitrary threshold of 1% on the p-value
abline(h=-log10(0.01), col="red")
## Remove NA values by replacing them with 1
res$padj=ifelse(is.na(res$padj), 1, res$padj)
res$pvalue=ifelse(is.na(res$pvalue), 1, res$pvalue)
## Draw a legend with the number of genes declared positive
legend("topleft", 
       legend=c(paste("pvalue<=", 0.01, ":", sum(res$pvalue<=0.01 & abs(res$log2FoldChange)>=1))))
