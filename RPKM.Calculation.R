library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)
library(GenomicAlignments)
dev.plot2pdf("file.pdf")
source('~/Sequencing Short Read Analysis/Count.R')
library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)
library(GenomicAlignments)
gtf <- import("refseq_hg19.gtf.gz", asRangedData=FALSE)
class(gtf)
idx <- mcols(gtf)$type == "exon"
exons <- gtf[idx]
genes <- split(exons, mcols(exons)$gene_id)
geneID
length(unique(gtf$gene_id))
idx <- mcols(gtf)$type == "exon"
exons <- gtf[idx]
genes <- split(exons, mcols(exons)$gene_id)
genes
params <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE))
bam <- readGAlignments("LY1_ATAC_chr18.bam", param=params)
counts <- summarizeOverlaps( features=genes,
reads=bam,
mode="Union",
ignore.strand=TRUE,
SingleEnd=FALSE,
param=params)
counts
count_table <- assays(counts, withDimnames=TRUE)$counts
count_table
count_table["BCL2",]
1000000*count_table/sum(count_table)
count_table["BCL2",]
reduce(gene)
reduce(genes)
genes[[1]]
width(reduce(genes[[1]]))
sum(width(reduce(genes[[1]])))
sum(width(reduce(genes)))
gene[[1]]
reduce(genes[[1]])
#reduce the exons that you want into overlapping regions
width(reduce(genes[[1]])
width(reduce(genes[[1]])
width(reduce(genes[[1]]
sum(widtch(reduce(genes)))
#normalize count by size of genes
count_table/sum(width(reduce(genes)))
count_table/sum(width(reduce(genes)))*1000
#normalize based on total number of genes
9count_table/sum(width(reduce(genes)))*1000)/sum(count_table)
(count_table/sum(width(reduce(genes)))*1000)/sum(count_table)
1000000*(count_table/sum(width(reduce(genes)))*1000)/sum(count_table)
write.table(count_table,
file="counts_table.txt", sep = "\t",
row.names=TRUE, col.names=FALSE, quote=F)
#counts per nucleotide: count_table/sum(width(reduce(genes)))*1000
