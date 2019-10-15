##/export/bio/R/R-3.2.1/bin/R

#source("http://bioconductor.org/biocLite.R")
#biocLite("DEXSeq")
library("DEXSeq", lib.loc="/netapp/home/singh/R/x86_64-unknown-linux-gnu-library/3.2")
#biocLite("BiocParallel")
library("BiocParallel", lib.loc="/netapp/home/singh/R/x86_64-unknown-linux-gnu-library/3.2")
#biocLite("statmod")
#library(statmod)
#library(parallel)
#library("Rmpi", lib.loc="/netapp/home/singh/R/x86_64-unknown-linux-gnu-library/3.2")
sessionInfo()
packageVersion("DEXSeq")

BPPARAM = MulticoreParam(workers=20)
#register(MulticoreParam(5))
#place with python scripts /export/bio/R/R-3.0.1/lib64/R/library/DEXSeq/python_scripts/

#1 covert gtf to DEXseq specific gff

#/export/bio/python/python2.7 /export/bio/R/R-3.0.1/lib64/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py Oreochromis_niloticus.BROADON2.gtf Oreochromis_niloticus.BROADON2.DEXseq.gff

#/export/bio/python/python2.7 /export/bio/R/R-3.0.1/lib64/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py Metriaclima_zebra.BROADMZ2.gtf Metriaclima_zebra.BROADMZ2.DEXseq.gff

#2 get counts for bam

#counts_dexseq.sh ## this didnt work on the server because python failed to import HTseq when qsubbed

#qsub -q all.q@romulus -V -pe mpi 8 -l h_vmem=8G -b y -cwd -N count -e counts_dexseq.err -o counts_dexseq.log /bin/sh counts_dexseq.sh # works only on romulus

#3 dexseq in R

#library(DEXSeq)
setwd("/netapp/home/singh/pooja_phd/data/pharyngeal_jaw_transcriptomes/raw_demultiplexed_bam/paper1/map2Onil/all_bam/alternative_splicing/")


# reads files
countFiles = list.files(".", pattern="DEXseqcount.txt$", full.names=TRUE)
countFiles

#read gff for DEXseq
flattenedFile = list.files("/netapp/home/singh/pooja_phd/data/cichlid_reference_genomes/new_Oreochromis/gtf_brawand/", pattern="seq.gff$", full.names=TRUE)
flattenedFile

# import sample table
sampleTable <- read.table("paper1_speciesinfo_Ld.txt", header=T, row.names=1)

#start analysis
dxd = DEXSeqDataSetFromHTSeq(countFiles, sampleData=sampleTable, design= ~ sample + exon + condition:exon, flattenedfile=flattenedFile)


# normlisation

dxd = estimateSizeFactors(dxd)

dxd = estimateDispersions(dxd, BPPARAM = BPPARAM)
#dxd = estimateDispersions(dxd)


#estimate dispersion


png("DEXseq_dispersion_estimates_Ld.png")
plotDispEsts(dxd)
dev.off()

#test of diff exp of exons

dxd = testForDEU(dxd, BPPARAM = BPPARAM)
#dxd = testForDEU(dxd)

#remove totalt gene expression effects

dxd = estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM = BPPARAM)
#dxd = estimateExonFoldChanges(dxd, fitExpToVar="condition")

#results

dxr1 = DEXSeqResults(dxd)
write.table(dxr1, file = "DEXseq_fullresults_Ld.txt", sep = "\t", quote=F)

#plot MA plot

png("DEXseq_MAplot_Ld.png")
plotMA(dxr1, cex=0.8 )
dev.off()

#DE isoforms at FDR 1%

genes <- table(dxr1$padj < 0.01)
write.table(genes, file = "DEXseq_sig_results_0.01_Ld.txt", sep = "\t", quote=F)


#DE genes affected at FDR 1%

genesN <- table ( tapply( dxr1$padj < 0.01, dxr1$groupID, any) )
write.table(genesN, file = "DEXseq_sig_results_0.01_genesAffected_Ld.txt", sep = "\t", quote=F)


#DE isoforms at FDR 5%

genes <- table(dxr1$padj < 0.05)
write.table(genes, file = "DEXseq_sig_results_0.05_Ld.txt", sep = "\t", quote=F)


#DE genes affected at FDR 5%

genesN <- table ( tapply( dxr1$padj < 0.05, dxr1$groupID, any) )
write.table(genesN, file = "DEXseq_sig_results_0.05_genesAffected_Ld.txt", sep = "\t", quote=F)
#netapp/home/tmp/RtmpAeSB6R/downloaded_packages


save(dxr1, file="./dexClust.RData")
