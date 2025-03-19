# Single Cell RNA seq project
# Cuvier's Team
# Schaak - Heurteau - Depierre* 
# 2017

# get a GRANGES of only CODING genes dm6 release i.e exlude non coding RNA

`%ni%` <- Negate(`%in%`) 
library(GenomicRanges)

# Get GRanges genes from dm6 version
source("https://bioconductor.org/biocLite.R")
biocLite("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
TxDb_dm6 = TxDb.Dmelanogaster.UCSC.dm6.ensGene

# Chromosome names =  Ensembl Style
seqlevelsStyle(TxDb_dm6) <- "Ensembl"
seqlevels(TxDb_dm6)[7]='Y'
genes_dm6 = genes(TxDb_dm6)
cds.gr.list = cdsBy(TxDb_dm6, by="gene")
coding_genes_dm6 = genes_dm6[names(genes_dm6) %in% names(cds.gr.list)]
seqlevels(coding_genes_dm6)[1:8]

# Focus on main chromosomes : 
#~ 1                        2L
#~ 2                        2R
#~ 3                        3L
#~ 4                        3R
#~ 5                         4
#~ 6                         X
#~ 7                         Y
#~ 8 dmel_mitochondrion_genome
coding_genes_dm6 = coding_genes_dm6[seqnames(coding_genes_dm6) %in% seqlevels(coding_genes_dm6)[1:8]]
seqlevels(coding_genes_dm6) = seqlevels(coding_genes_dm6)[1:8]


length(coding_genes_dm6)
#~ [1] 13899

# delete genes (1 gene here) which have TSS position >5000pb from chromosome border (because profiles are done on this window, only 1 gene is deleted)
coding_genes_dm6=coding_genes_dm6[names(coding_genes_dm6) %ni% names(coding_genes_dm6[start(coding_genes_dm6)<5000])]

length(coding_genes_dm6)
#~ [1] 13898

saveRDS(coding_genes_dm6, "/home/belhocine/Bureau/work_ddepierre/David_M2_K9K27DE/DHMG/r6.13/coding_genes_dm6_Granges.rds")
saveRDS(coding_genes_dm6, "/home/belhocine/Bureau/work_ddepierre/David_M2_K9K27DE/DEAnalysis/r6.13/coding_genes_dm6_Granges.rds")

saveRDS(cod_genes, "/home/belhocine/Bureau/work_ddepierre/David_M2_K9K27DE/DEAnalysis/r6.13/coding_genes_dm6_Granges.rds")
saveRDS(cod_genes, "/home/belhocine/Bureau/work_ddepierre/David_M2_K9K27DE/DHMG/r6.13/coding_genes_dm6_Granges.rds")










