library(gplots)
palette<-colorRampPalette(c("darkblue","white","darkred"))(256)

cor_tCOL <- 1 - cor(data.matrix(Tabela))
#cor_tCOL <- 1 - cor(data.matrix(counts)
distancet <- as.dist(cor_tCOL)
hclust_complete <- hclust(distancet, method = "complete")
dendcomplete <- as.dendrogram(hclust_complete)

cor_tROW <- 1 - cor(t(data.matrix(counts)))
#cor_tROW <- 1 - cor(t(data.matrix(tabelaTPMtx[1:1000,:9])))
distancetROW <- as.dist(cor_tROW)
hclust_completeROW <- hclust(distancetROW, method = "complete")
dendcompleteROW <- as.dendrogram(hclust_completeROW)

pdf(file = "TPMdoubledendo_correlation_correlation.pdf", width =400, height =400)

heatmap.2(as.matrix(counts),dendrogram="both",scale="row",
col=palette,margins=c(10,10),cexCol=1,trace='none',main="teste",
#labCol=c("sample1", ...,"sample6")
,Rowv=FALSE
,labRow = NULL
,Colv=dendcomplete)

dev.off()

