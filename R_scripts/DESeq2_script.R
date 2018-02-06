
library("DESeq2")

#Create a DESeqDataset(dds) from the count table and corresponding metadata
dds=DESeqDataSetFromHTSeqCount(sampleTable = Tcell, design= ~ Subset) 


#Prefiltering option. Remove rows in which there are very few reads, reducing the memory size of 
#the dds data object, and increase the speed of the transformation and testing functions within 
#DESeq2. More strict filtering to increase power is automatically applied via independent filtering 
#on the mean of normalized counts within the results function.

#keep = rowSums(counts(dds)) >= 10
#dds = dds[keep,]

#--------Differential expression analysis-------------##

#The standard differential expression analysis steps are wrapped into a single function, DESeq.
dds=DESeq(dds)

#Results tables are generated using the function results, which extracts a results table with 
#log2 fold changes, p values and adjusted p values. The comparison can be ordered using the 
#contrast arguments of results. In this case, logFC=log2(Treg/Tconv) 

res = results(dds, contrast=c("Subset","Treg","Tconv"))
summary(res) #para ver summary com outro pvalue: summary(res,alpha=x)

out of 32530 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 1347, 4.1% 
LFC < 0 (down)   : 885, 2.7% 
outliers [1]     : 106, 0.33% 
low counts [2]   : 9671, 30% 
(mean count < 4)

#To see how many many adjusted p-values were less than 0.1?
#sum(res$padj < 0.1, na.rm=TRUE)
#[1] 2232


#Results function automatically performs independent filtering based on the mean of normalized 
#counts for each gene, optimizing the number of genes which will have an adjusted p value below 
#a given FDR cutoff, alpha. By default the argument alpha is set to 0.1. To adjust pvalue cutoff, 
#alpha must be set to that value:
#res05 = results(dds, contrast=c("Subset", "Treg", "Tconv"), alpha=0.05)

#----------------SAMPLE CLUSTERING AND VISUALIZATION--------##


#             Heatmap of the sample-to-sample distances


#Another use of the transformed data is sample clustering. Here, dist function is applied to  
#transpose the transformed count matrix to get sample-to-sample distances, which requires vsd(variance-stabilizing transformation)
#A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between
#samples:

vsd = varianceStabilizingTransformation(dds)
sampleDists = dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) =colnames(dds)
colnames(sampleDistMatrix) = colnames(dds)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
library(pheatmap)
pheatmap (sampleDistMatrix, 
          clustering_distance_rows=sampleDists,
          clustering_distance_cols=sampleDists, 
          col=colors)  


# Principal component analysis (PCA) plot of the samples

#Related to the distance matrix is the PCA plot, which shows the samples in the 2D plane spanned by their
#first two principal components (subset and replicates). This type of plot is useful for visualizing 
#the overall effect of experimental covariates and batch effects.
par(mar=c(5,10,1,2))
plotPCA(vsd, intgroup = c("Subset"))

#----------------EXPLORING AND EXPORTING RESULTS-------------##

#MA-plot shows the log2 fold changes attributable to a given variable over the mean of normalized counts 
#for all samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. 
#Points which fall out of the window are plotted as open triangles pointing either up or down.

plotMA(res, ylim=c(-10,10), main="Tconv vs. Treg")
#plot(res$baseMean, res$log2FoldChange,
     #main="MA plot",
     #xlab= "Mean of normalized counts",
     #ylab = "Log2FC (Treg/Tconv)",
     #col=densCols(res$log2FoldChange, res$baseMean),
     #grid(lty="solid",col="grey"),
     #pch=20)
#abline(col="black", h=0)
#abline(col="red", h=c(1, -1))

#MA-plot for the shrunken log2 fold changes: Remove the noise associated with log2 fold changes from 
#low count genes without requiring arbitrary filtering thresholds.

resLFC = lfcShrink(dds, coef=2)
plotMA(resLFC, ylim=c(-6,6), main="Tconv vs. Treg")


#Plotting the dispersion estimates is a useful diagnostic. The dispersion plot below is 
#typical, with the final estimates shrunk from the gene-wise estimates towards the fitted 
#estimates. Some gene-wise estimates are flagged as outliers and not shrunk towards the 
#fitted value. The amount of shrinkage can be more or less than seen here, depending on the 
#sample size, the number of coefficients, the row mean and the variability of the gene-wise estimates.

plotDispEsts(dds)


#Examine the counts of reads for a single gene across the groups: plotCounts normalizes counts by sequencing 
#depth and adds a pseudocount of 1/2 to allow for log scale plotting. The counts are grouped by the 
#variables in intgroup, where more than one variable can be specified. which.min specify the gene which 
#had the smallest p adj value from the results table (res) created above. You can select the gene to plot 
#by rowname or by numeric index.

plotCounts(dds, gene=which.min(res$padj), intgroup="Subset", main="FoxP3")
#plotCounts(dds, gene=rownames(res, "ENSG00000112149.9"), intgroup="Subset", main="FoxP3")



#To convert  the Ensembl ID to gene name:
ensembl.all = read.table('gencode.v26.annotation.gtf', sep='\t', header =F);
ensembl.gene.names = data.frame(sapply(strsplit(sapply(strsplit(as.character(ensembl.all[,9]),'gene_name '), "[", 2), ";"), "[", 1));
ensembl.id.names = data.frame(sapply(strsplit(sapply(strsplit(as.character(ensembl.all[,9]),'gene_id '), "[", 2), ";"), "[", 1));
ensembl.code = cbind(ensembl.gene.names, ensembl.id.names);
ensembl.code = ensembl.code[!duplicated(ensembl.code[,2]),];
rownames(ensembl.code) = ensembl.code[,2];
colnames(ensembl.code) = c('gene', 'id');
res.names = merge(as.data.frame(res), ensembl.code, by="row.names", all.x=F);

#Guardar tabela com resultados (ordenado por pvalue)	
resOrderedPValue = res.names[order(res.names$pvalue),]
write.table(resOrderedPValue, file="toptags_Deseq2.txt",sep="\t",quote=F)


#END
##############################################################

              Outros grÃ¡ficos:		
                
                pdf(file = "test.pdf")
              plotCounts(dds, gene=which.min(res$padj), intgroup="Subset") 
              d <- plotCounts(dds, gene=which.min(res$padj), intgroup="Subset",
                              returnData=TRUE)
              library("ggplot2")
              ggplot(d, aes(x=condition, y=count)) +
                geom_point(position=position_jitter(w=0.1,h=0)) +
                scale_y_log10(breaks=c(25,100,400))
              dev.off()
              
              rld <- rlog(cds, blind=FALSE)
              vsd <- varianceStabilizingTransformation(cds, blind=FALSE)
              vsd.fast <- vst(cds, blind=FALSE)
              
              library("pheatmap")
              select <- order(rowMeans(counts(cds,normalized=TRUE)),decreasing=TRUE)[1:20]
              nt <- normTransform(cds) # defaults to log2(x+1)
              log2.norm.counts <- assay(nt)[select,]
              df <- as.data.frame(colData(cds)[,c("Subset","LibraryLayout")])
              sampleDists <- dist(t(assay(rld)))
              library("RColorBrewer")
              sampleDistMatrix <- as.matrix(sampleDists)
              rownames(sampleDistMatrix) <- paste(rld$Subset)
              colnames(sampleDistMatrix) <- NULL
              colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
              
              pdf(file="heatmap_deseq2.pdf", width=10, height=10)
              pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
              dev.off()
              
              pdf(file="PCA_deseq2.pdf", width=10,height=10)
              data <- plotPCA(vsd, intgroup=c("Subset"), returnData=TRUE)
              percentVar <- round(100 * attr(data, "percentVar"))
              ggplot2(data, aes(PC1, PC2, color=condition1, shape=type)) +
                geom_point(size=3) +
                xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                ylab(paste0("PC2: ",percentVar[2],"% variance"))
              dev.off()
              
              
              pdf(file = "vulcano_plots_DESeq2.pdf", width =10, height =10)
              with(res, plot(res$log2FoldChange, -log10(res$padj), pch=20, main="Tconv vs. Treg", xlim=c(-10,10)))
              #with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
              #with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
              dev.off()
              