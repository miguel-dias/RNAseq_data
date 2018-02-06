##------------------------------DATA PROCESSING
#DESEQ2 TABLE
library(DESeq2)
resUP<-subset(res.names,res.names$log2FoldChange>0) #get the genes that are up-regulated
resDOWN<-subset(res.names,res.names$log2FoldChange<0) #get the genes that are down-regulated

resUP$FCabs=2^(resUP$log2FoldChange) #add a column with FC absolute values in ttUP table
resDOWN$FCabs=-1/(2^(resDOWN$log2FoldChange)) #add a column with FC absolute values in ttDown table
resFCabs<-rbind(resUP,resDOWN) #bind the two tables 

resFCabs$padj=ifelse(is.na(resFCabs$padj), 1, resFCabs$padj)#replace all padj NA values with 1 (outliers)
resFCabs$pvalue=ifelse(is.na(resFCabs$pvalue), 1, resFCabs$pvalue)#replace all pvalue NA values with 1 (outliers)
hist(resFCabs$padj)

respvalue01=subset(resFCabs,resFCabs$pvalue<0.01) #Create a table with genes whose pvalue is < 0.01
resP01FC2=subset(respvalue01,abs(respvalue01$FCabs)>=2)
write.csv(resP01FC2$Row.names,file="toptags_P01FC2_DESeq2.csv")
nrow(resP01FC2[resP01FC2$FCabs>=2,]) # How many genes are upregulated
nrow(resP01FC2[resP01FC2$FCabs<=2,]) #  How many genes are downregulated
write.csv(resP01FC2, file="toptags_P01FC2_DESeq2.csv", quote=F) #export csv file
#EDGER TABLE
library(edgeR)
ttUP<-subset(tt.names,tt.names$logFC>0)#get the genes that are up-regulated
ttDOWN<-subset(tt.names, tt.names$logFC<0)#get the genes that are down-regulated
ttUP$FCabs=2^(ttUP$logFC)#add a column with FC absolute values in ttUP table
ttDOWN$FCabs=-1/(2^(ttDOWN$logFC))#add a column with FC absolute values in ttDown table
ttFCabs<-rbind(ttUP,ttDOWN)#bind the two tables 
#subset(ttFCabs,ttFCabs[,4]=NA)
#ttFCabs[is.na(ttFCabs)] = 1 #replace NA with 1 (outliers)
hist(ttFCabs$FDR) #inspect FDR values
ttpvalue01=subset(ttFCabs,ttFCabs$PValue<0.01) #Create a table with genes whose pvalue is < 0.01
ttP01FC2=subset(ttpvalue01,abs(ttpvalue01$FCabs)>=2)#Create a table with genes whose abs(FCabs)>=2
nrow(ttP01FC2[ttP01FC2$FCabs>=2,]) #  How many genes are upregulated: 
nrow(ttP01FC2[ttP01FC2$FCabs<=2,]) #  How many genes are downregulated:
write.csv(ttP01FC2,file="toptags_P01FC2_edger.csv", quote=F) #export a csv file 

#merge das duas tabelas por gene id
Finaltable=merge(as.data.frame(resP01FC2),ttP01FC2,by="Row.names")
Finaltable$FCcheck=(Finaltable$FCabs.x)*(Finaltable$FCabs.y) #check FCabs values
Finaltable=Finaltable[ -c(9,16)]# Remove duplicated columns
write.csv(Finaltable,file="HClist_DEgenes.csv", quote = F)
nrow(Finaltable[Finaltable$FCabs.x>=2,])
nrow(Finaltable[Finaltable$FCabs.x<=2,])

#Geneontology bar plot of all DE genes
#Run gene ontology of Biological processes and download the results
BP_all=read.table(file="Gene_ontology/WebGestalt/Biological_Process_P01/enrichment_results_wg_result1516018401.txt", header=T, sep = "\t")
BP_all1=BP_all[,-c(1,3,4,6:11)]
par(mar=c(3,21,1,1))
barplot(BP_all1[c(15:1),"O"], horiz = TRUE,
        main = "Biological Process of all DE genes", 
        names.arg = BP_all1[c(15:1),"description"], 
        space=1, las=2, col = "Red", border = F,
        xlim=c(0,200))

#Geneontology bar plot of Down-regulated genes
Down=subset(Finaltable, Finaltable$FCabs.x<0)
write.csv(Down$Row.names, file = "Down4GO.csv", sep=",", quote=F)
BP_Down=read.table(file="Gene_ontology/WebGestalt/Biological_Processes_P01_Down/enrichment_results_wg_result1516019732.txt", header=T, sep = "\t")
par(mar=c(3,19,1,1))
barplot(BP_Down[c(15:1),"O"], horiz = TRUE,
        main = "Biological Process of Down-regulated genes",
        names.arg = BP_Down[c(15:1),"description"], space=1, las=2,
        col = "Red", border = F,
        xlim=c(0,140))

#Geneontology bar plot of UP-regulated genes
UP=subset(Finaltable, Finaltable$FCabs.x>0)
write.csv(UP$Row.names, file = "UP4GO.csv", sep=",", quote=F)
BP_UP=read.table(file="Gene_ontology/WebGestalt/Biological_Processes_P01_Up/enrichment_results_wg_result1516019899.txt", header=T, sep = "\t")
par(mar=c(3,17,1,1))
barplot(BP_UP[c(15:1),"O"], horiz = TRUE,
        main = "Biological Process of Up-regulated genes", 
        names.arg = BP_UP[c(15:1),"description"], space=1, las=2, 
        col = "Red", border = F,
        xlim=c(0,140))

#Pathways analyses bar plot of all DE genes
Path_all=read.table(file="Gene_ontology/WebGestalt/Pathway_allDEgenes_P01/enrichment_results_wg_result1516027648.txt", header=T, sep = "\t")
Path_all=Path_all[,-c(1,3,4,6:11)]
par(mar=c(3,18,1,1))
barplot(BP_UP[c(15:1),"O"], horiz = TRUE, 
        main = "Enriched pathways of all DE genes", 
        names.arg = Path_all[c(15:1),"description"], space=1, 
        las=2, col = "Red", border = F,
        xlim=c(0,140))
#Pathways analyses bar plot of UP DE genes
Path_Up=read.table(file="Gene_ontology/WebGestalt/Pathways_UP/enrichment_results_wg_result1516624090.txt", header=T, sep = "\t")
Path_Up=Path_Up[,-c(1,3,4,6:11)]
par(mar=c(3,18,1,1))
barplot(Path_Up[c(10:1),"O"], horiz = TRUE, 
        main = "Enriched pathways of up-regulated genes", 
        names.arg = Path_Up[c(10:1),"description"], space=1, 
        las=2, col = "Red", border = F,
        xlim=c(0,40))
#Pathways analyses bar plot of DOWN DE genes
Path_Down=read.table(file="Gene_ontology/WebGestalt/Pathways_DOWN/enrichment_results_wg_result1516624482.txt", header=T, sep = "\t")
Path_Down=Path_Down[,-c(1,3,4,6:11)]
par(mar=c(3,21,1,1))
barplot(Path_Down[c(10:1),"O"], horiz = TRUE, 
        main = "Enriched pathways of Down-regulated genes", 
        names.arg = Path_Down[c(10:1),"description"], space=1, 
        las=2, col = "Red", border = F,
        xlim=c(0,20))
#Scatter plot Tconv vs Treg using counts table
counts_1=cbind(counts, Tconv_mean=rowMeans(counts[ ,1:3]))
counts_2=cbind(counts_1, Treg_mean=rowMeans(counts[ ,4:6]))
counts_3=counts_2[ ,-1:-2:-3:-4:-5:-6]
counts_4=log2(1+counts_3)

axis.range <- c(0, ceiling(max(counts_4[,c("Tconv_mean","Treg_mean")])))
par(mar=c(5,8,2,1))
plot(counts_4,
     main="",
     xlab="Tconv 
     log2(read counts + 1)",
     ylab="Treg 
     log2(read counts + 1)",
     xlim=axis.range, ylim=axis.range,
     pch=20,
     col = densCols(x=counts_4[ ,2], y=counts_4[ ,1]))
grid(lty="solid", col="darkgray")

abline(h=0, col="black")
abline(v=0, col="black")
abline(a=0, b=1)
abline(a=1, b=1, col="red")
abline(a=-1, b=1, col="red")

#Check for outliers???
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
 


#Bar plot of Downregulated genes with the number of genes after each bar
BP_Down$log10Pvalue=(-log10(BP_Down$PValue))
a <- BP_Down[c(15:1),"log10Pvalue"]
b <- BP_Down[c(15:1),"O"]
par(mar=c(4,4,1,1), oma=c(0,0,0,0), las=1)
posbar <- barplot(a, horiz = TRUE, col="darkolivegreen3", 
                  main="Biological Process of Down-regulated genes",
                  space=1, xlim=c(0,10),
                  names.arg = BP_Down[c(15:1),"description"])
mtext(side=1, text="-Log10(Pvalue)", line=2.5)
text(y=posbar, x=a, pos=4,labels=b)
#Bar plot of Upregulated genes with the number of genes after each bar
BP_UP$log10Pvalue=(-log10(BP_UP$PValue))
c <- BP_UP[c(15:1),"log10Pvalue"]
d <- BP_UP[c(15:1),"O"]
par(mar=c(4,4,1,1), oma=c(0,0,0,0), las=1)
posbar <- barplot(c, horiz = TRUE, col="darkolivegreen3", 
                  main="Biological Process of Up-regulated genes",
                  ylab = "", xlab="", xlim=c(0,20),
                  names.arg = BP_Down[c(15:1),"description"])
mtext(side=1, text="-Log10(Pvalue)", line=2.5)
text(y=posbar, x=c, pos=4,labels=d)

#####################################################
#TF_uniprot=read.table("../../../Downloads/Transcription_factors_DNAbinding_uniprot.txt", sep = "\t", header = T)
#TF_UPtable=merge(TF_uniprot, UP, by="gene.x")
#TF_Downtable=merge(TF_uniprot, Down, by="gene.x")
#TF=merge(TF_Down_GTRD, TF_Downtable, all.x=T, all.y=T) #merge the two tables from uniprot and GTRD databases
#PM_uniprot=read.table("../../../Downloads/Plasma_membrane_uniprot.txt", sep = "\t", header = T)
#PM_UPtable=merge(PM_uniprot, UP, by="gene.x")
#PM_Downtable=merge(PM_uniprot, Down, by="gene.x")
#write.table(PM_UPtable, file = "Plasma_membrane_UP_table.csv", sep = ";", quote = F)
#write.table(PM_Downtable, file = "Plasma_membrane_Down_table.csv", sep = ";", quote = F)

###
TF_GTRD=read.table("TF_from_GTRD.csv", sep = ";", header = T)
colnames(TF_GTRD)=gsub("ï..", "", colnames(TF_GTRD))
TF_UP_GTRD=merge(TF_GTRD, UP, by="Row.names")
TF_Down_GTRD=merge(TF_GTRD, Down, by="Row.names")
#Inspect which of these are regulated by Foxp3
Foxp3_regulated_genes=read.table("FoxP3-regulated_genes.txt", header=T)
Foxp3UP=merge(TF_UP_GTRD, Foxp3_regulated_genes, by="Row.names")
Foxp3DOWN=merge(TF_Down_GTRD, Foxp3_regulated_genes, by="Row.names")

#TFDOWN=merge(TF_Down_GTRD, TF_Downtable, all.x=T, all.y=T) #merge the two tables from uniprot and GTRD databases
#TFUP=merge(TF_UP_GTRD, TF_UPtable, all.x=T, all.y=T) #merge the two tables from uniprot and GTRD databases
#####################################################
#Scatter plot Tconv vs Treg using nc table
#nc_1=cbind(nc, Tconv_mean=rowMeans(nc[ ,1:3]))
#nc_2=cbind(nc_1, Treg_mean=rowMeans(nc[ ,4:6]))
#nc_3=nc_2[ ,-1:-2:-3:-4:-5:-6]
#nc_4=log2(1+nc_3)

#axis.range <- c(0, ceiling(max(nc_4[,c("Tconv_mean","Treg_mean")])))
#par(mar=c(5,8,2,1))
#plot(nc_4,
     #main="",
     #xlab="Tconv 
     #log2(read counts + 1)",
     #ylab="Treg 
    #log2(read counts + 1)",
    #xlim=axis.range, ylim=axis.range,
     #pch=20,
     #col = densCols(x=nc_4[ ,2], y=nc_4[ ,1]))
#grid(lty="solid", col="darkgray")

#abline(h=0, col="black")
#abline(v=0, col="black")
#abline(a=0, b=1)
#abline(a=1, b=1, col="red")
#abline(a=-1, b=1, col="red")
