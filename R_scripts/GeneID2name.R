##http://guertinlab.org/wp-content/uploads/2015/02/ChIP_RNA_seq_analysis.pdf

These genes all have Ensembl gene indentifiers, but it is useful to have gene names. To retrieve a ’key’
to convert as the Ensembl ID to gene name do the following.

ensembl.all = read.table('gencode.v26.annotation.gtf', sep='\t', header =F);
ensembl.gene.names = data.frame(sapply(strsplit(sapply(strsplit(as.character(ensembl.all[,9]),'gene_name '), "[", 2), ";"), "[", 1));
ensembl.id.names = data.frame(sapply(strsplit(sapply(strsplit(as.character(ensembl.all[,9]),'gene_id '), "[", 2), ";"), "[", 1));
ensembl.code = cbind(ensembl.gene.names, ensembl.id.names);
ensembl.code = ensembl.code[!duplicated(ensembl.code[,2]),];
rownames(ensembl.code) = ensembl.code[,2];
colnames(ensembl.code) = c('gene', 'id');
res.names = merge(resOrdered, ensembl.code, by="row.names", all.x=F);
#res.names = res.names[,c(1:8)];
resnamesOrdered = res.names[order(res$pvalue),]
write.csv(resnamesOrdered, file="Deseq2_results.csv")


Next we may want to select specific rows. We can seach our favorite gene, or find genes that pass our
adjusted p-value threshold. Note, that all NA values will be returned unless we specify ’!is.na()’.

esr1.gene = res.names[res.names$gene == 'ESR1',];
print(esr1.gene);
low.padj = res.names[res.names$padj < 0.1 & !is.na(res.names$padj),];
print(low.padj);
