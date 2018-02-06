library("tximport")
dir="Kallisto_results"
list.files(dir)
samples=read.table(file.path(dir, "abudance_files.csv"), sep=";", header=TRUE)
files=file.path(dir, samples$File)
all(file.exists(files))
#[1] TRUE
library(rhdf5)
txi.kallisto=tximport(files, type=c("kallisto"),
                      txIn=TRUE, txOut=FALSE, countsFromAbundance = c("no"),
                      tx2gene=NULL)
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = FALSE)
