#Calculate the Mean centered fold expression level (depth adjusted reads per million) across samples
data1=nc
data1=cbind(data1, meancentered=rowMeans(data1))        
data2=data1/data1[ ,7]       
data2<-data2[,-7]
data2=cbind(data2, Row.names=row.names(data2))

#Heatmap of all the DE genes in the FinalTable 
data2_df=as.data.frame(data2) #transform data2 matrix to a dataframe
common=merge(Finaltable, data2_df, by="Row.names")#merge the two tables
row.names(common)=common$Row.names
common=common[,-c(1:16)]
heattable=common[order(common[ ,1], decreasing = T),]
write.table(heattable, "table4heatmap.txt", sep="\t", quote=F)
heattable_2=read.table(file="table4heatmap.txt", header=T)
library(gplots)
heatmap.2(as.matrix(heattable_2), dendrogram ="column"
          , col=redgreen(75), scale="row",key=TRUE
          , Rowv=F, Colv=T, symkey=FALSE, labRow = F, density.info="none", trace="none", cexRow=0.5)





#Heatmap od the 300 genes with the smallest pvalue (from edger values) 
#head(sort(Finaltable$pvalue ,decreasing=FALSE), n = 300)
#list1<-head(Finaltable[order(Finaltable$pvalue, decreasing= F),], n = 300)
#data2_df=as.data.frame(data2) #transform data2 matrix to a dataframe
#common=merge(list1,data2_df, by="Row.names")
#row.names(common)=common$Row.names
#common=common[,-c(1:18)]
#final=common[order(common[ ,1], decreasing = T),]
#write.table(final, "final.txt", sep="\t", quote=F)
#test=read.table(file="final.txt", header=T)
#library(gplots)
#heatmap.2(as.matrix(test), dendrogram ="column"
 #         , col=redgreen(75), scale="row",key=TRUE
  #        , Rowv=F, Colv=T, symkey=FALSE, labRow = F, density.info="none", trace="none", cexRow=0.5)
