#Import transcript-level estimates
BiocManager::install("tximport")
library(tximport)
library(tximportData)
library(rhdf5)
library(biomaRt)
library(readr)


#modify header names in *.h5 files
filh5 <- list.files("first_project/kallisto_results", pattern=".h5", recursive=TRUE, full.names=TRUE)
for (currentFile in filh5) {
  oldids <- h5read(currentFile, "/aux/ids")
  #newids <- gsub("\\|.*", "", oldids)
  h5write(oldids, currentFile, "/aux/ids")
}

#View(h5ls(filh5[[1]]))
#upload metadata 
metadata <- read_csv("first_project/metadata.txt")

#rename with sample with run
names(filh5)= metadata$BioSample

# change kallisto output to input for deseq2
txi.kallisto <- tximport(filh5, type = "kallisto", txOut = TRUE)
head(txi.kallisto$counts)
counts= txi.kallisto[["counts"]]
abundance= txi.kallisto[["abundance"]]
write.csv(abundance, file="kalisto_abundance.csv")
write.csv(counts, file="kalisto_counts.csv")

#_____________________________DeseQ2___________________

library(DESeq2)

rownames(metadata)= metadata$BioSample

cond1="Healthy control" 
cond2="Diabetic Patient"

dds <- DESeqDataSetFromTximport(txi.kallisto, metadata, ~disease_state)
dds <- DESeq(dds)
res=results(dds, contrast = c("disease_state",cond2 ,cond1))
res=as.data.frame(res)
res=res[complete.cases(res), ]# remove nulls
write.csv(res, file = "res.csv", row.names = T,quote = F)
# get the significant
#DEGs.fc1.1=res[res$padj< 0.05 & abs(res$log2FoldChange)>1.1,]


# change transcripts to gene Ids
#name mapping
library(stringr)
ids <- as.list(rownames(res))
ensemble.id=sapply(ids, function(x) strsplit(as.character(x),"\\.")[[1]][1])
rownames(res)= ensemble.id

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
genemap <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","hgnc_symbol"),
                 filters = "ensembl_transcript_id",
                 values = ensemble.id,
                 mart = mart )

symbols <- tapply(genemap$hgnc_symbol, genemap$ensembl_transcript_id, paste, collapse="; ")
res$symbol <- symbols[ rownames(res) ]
rownames(res) = NULL
res=res[ ! duplicated(res$symbol),]
res=res[ ! is.na(res$symbol),]
rownames(res) = res$symbol
res = res[,-7]

write.csv(res, file="res.csv")

#####################################
##gprofiler###

#install.packages("gprofiler2")
#library(gprofiler2)
#gp_up = gost(list(row.names(res)), organism = "hsapiens")










