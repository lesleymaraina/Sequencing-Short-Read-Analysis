library(GenomicRanges) 
library(Rsamtools) 
library(rtracklayer) 
library(GenomicAlignments)

#setwd("TEACHING")

# import GTF 
gtf <- import("refseq_hg19.gtf.gz", asRangedData=FALSE)

# what did we just get ?
class(gtf)

# EXERCISE: how many genes in GTF ?

# which of these rows are exons?
# find where exons are in gtf
# get rid of anything thats not an exon in a gene
idx <- mcols(gtf)$type == "exon"

# get exons only
exons <- gtf[idx]

# genes
# for each entiity collect index and exon
genes <- split(exons, mcols(exons)$gene_id)


# map reads to exon
# set up BAM scanning parameters
# only look at certain reads, and add flag
# dont add in memory the reads that are not mapped
params <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE))
  
# scan BAM file
bam <- readGAlignments("LY1_ATAC_chr18.bam", param=params)

# calculate counts
# map reads to diff genes
counts <- summarizeOverlaps( features=genes,
                             reads=bam,
                             mode="Union", #define what we mean by overlap
                             ignore.strand=TRUE,
                             SingleEnd=FALSE,
                             param=params) 

# what did we get ?
counts

# make count table
count_table <- assays(counts, withDimnames=TRUE)$counts

# EXERCISE 1: how many counts for BCL2 ? 247

# EXERCISE 2: convert to count per million reads (cpm)

# EXERCISE 3: calculate RPKM ?

# write count table to disk
write.table(count_table,
            file="counts_table.txt", sep = "\t",
            row.names=TRUE, col.names=FALSE, quote=F)
