
#Load the edgeR package and use the utility function, readDGE, to read in the COUNT files created from 
#htseq-count:
#Criar tabela no excel -> gravar como  csv -> abrir edgeR -> importar:

library(edgeR)
Tcell=read.csv ("Count_files/table.csv",header=T,sep=";") 

head(Tcell)

Name         Countf Subset
Tconv274 Tconv274.count  Tconv
Tconv276 Tconv276.count  Tconv
Tconv277 Tconv277.count  Tconv
Treg274  Treg274.count   Treg
Treg276  Treg276.count   Treg
Treg277  Treg277.count   Treg

#Read the count files created from HTseq-counts:
counts = readDGE(Tcell$Countf)$counts  

#visualize:
head(counts)

#Samples
#Tags                 Tconv274 Tconv276 Tconv277 Treg274 Treg276 Treg277
#ENSG00000000005.5         0        0        0       0       0       0
#ENSG00000000419.12     2938     2426     4221    2024    2528    3453
#ENSG00000000457.13     2273     2405     2823    1793    2524    3161
#ENSG00000000460.16      487      553      671     524     555     598
#ENSG00000000938.12       12       28       50      10       0       4
#ENSG00000000971.15       37       90      105       0       0       6


#Filter weakly expressed and noninformative(e.g., non-aligned) features:
  
noint = rownames(counts) %in% c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
cpms = cpm(counts)
keep=rowSums(cpms>1)>=3 & !noint #remove featres without at least 1 read per million in n (3) of the samples, where n is the size of smallest group of replicates.
counts=counts[keep,]

  
#Visualize and inspect the count table as follows:
head(counts) 
tail(counts)
dim(counts) 

#Create a DGEList object (edgeRs container for RNA-seq count data):
  
d = DGEList(counts = counts, group = Tcell$Subset)

#Estimate normalization factors using:
  
d = calcNormFactors(d)

#Inspect the relationships between all pairs ofsamples using a multidimensional scaling (MDS) plot:
#dim1=relação entre subsets
#dim2=relaçao intra subsets
par(mar=c(5,5,1,2)) 
plotMDS(d, labels = Tcell$shortname,col = c( "darkgreen","blue")[factor(Tcell$Subset)],
        xlim=c(-5/2,5/2))

#Estimate tagwise dispersion:
d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)

#Create a visual representation of the mean-variance relationship using the plotMeanVar and plotBCV functions:
  
plotMeanVar(d, show.tagwise.vars = TRUE, NBline = TRUE) #explore the mean-variance relationship, where each
#dot represents the estimated mean and variance of each gene, with binned variances as well as the trended 
#common dispersion overlaid.
plotBCV(d) #illustrated the relationship of biological coefficient of variaion vs mean log CPM.

#Test for differential expression as follows(Treg/Tconv):
de = exactTest(d, pair = c("Tconv","Treg"))
#summary(decideTests(de, p.value=0.01))

#Tconv+Treg
#-1        364
#0       11962
#1         619

#Use the topTags function to present a tabular summary of the differential expression
#statistics, ordered by pvalue:
tt = topTags(de, n = nrow(d))
head(tt$table)

#Inspect the depth-adjusted reads per milion for the 5 top differentially expressed genes:
nc = cpm(d, normalized.lib.sizes = TRUE)
rn = rownames(tt$table)
head(nc[rn,order(Tcell$Subset)])

#Samples
#Tags                  Tconv274   Tconv276   Tconv277   Treg274   Treg276  Treg277
#ENSG00000049768.14  6.808078 11.7342478 11.0705458 491.80501 435.00555 436.1491
#ENSG00000134460.15  1.630973  0.8220139  0.9434554 109.17294  97.35805 164.7392
#ENSG00000163599.14 17.767724 22.6362065 22.7332604 368.23805 344.24201 379.5911
#ENSG00000115414.18  1.309721  0.8939401  0.8129775  48.31303  40.86810  59.3879
#ENSG00000112149.9   6.004947  4.7676803  7.1461728 106.53221 112.43629 118.7859

#Create a graphical summary such as an M (log-fold change) versus A (log-average expression) 
#plot here showing the genes selected as differentially expressed (with a 10% false discovery rate:
                                              
deg = rn[tt$table$PValue<0.01] #para restringirmos os resultados em funçao do p-value:  substituir FDR por PValue]
plotSmear(d, de.tags = deg, xlim=c(-3,15), ylim=c(-10,9))  #plots the logFC against the log counts per million.
plot(tt$table$logCPM, tt$table$logFC,
     main="MA plot",
     xlab= "Average logCPM",
     ylab = "LogFC (Treg/Tconv)",
col=densCols(tt$table$logFC,tt$table$logCPM),
grid(lty="solid",col="grey"),
pch=20,
xlim=c(-3,15), ylim=c(-10,10))
abline(col="black", h=0)
abline(col="red", h=c(1, -1))

#To convert the Ensembl ID to gene name:
ensembl.all = read.table('gencode.v26.annotation.gtf', sep='\t', header =F);
ensembl.gene.names = data.frame(sapply(strsplit(sapply(strsplit(as.character(ensembl.all[,9]),'gene_name '), "[", 2), ";"), "[", 1));
ensembl.id.names = data.frame(sapply(strsplit(sapply(strsplit(as.character(ensembl.all[,9]),'gene_id '), "[", 2), ";"), "[", 1));
ensembl.code = cbind(ensembl.gene.names, ensembl.id.names);
ensembl.code = ensembl.code[!duplicated(ensembl.code[,2]),];
rownames(ensembl.code) = ensembl.code[,2];
colnames(ensembl.code) = c('gene', 'id');
tt.names = merge(tt$table, ensembl.code, by="row.names", all.x=F);

#Order table by logFC:
#ttnamesOrderedlogFC = tt.names[order(tt.names$logFC),]
#Save the result table as a CSV file (alternative formats are possible) :
write.csv(tt.names,file="toptags_edgeR.csv", sep="/t", quote = F)

#   END
###################################################### 
#Make a basic volcano plot
#with(tt$table, plot(logFC, -log10(PValue), pch=20, main="Volcano_plot_edgeR", xlim=c(-10,10)))
# Add red points if -log10(pvalue)>2 and abs log2FC>=2)
#with(subset(tt$table, PValue<0.01 & abs(logFC)>=1), points(logFC, -log10(PValue), pch=20, col="red"))


#[para guardar em pdf:
#    >pdf(file = "nome do pdf.pdf", width =10, height =10)â€‹
#  >plotSmear(d, de.tags = deg, pair = c("Tconv","Treg"), main="Tconv vs. Treg")
#  >dev.off()
#  ]
#####################################################################