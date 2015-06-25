library(Rsamtools)

# set wd 
#setwd("TEACHING")

bamfile <- 'LY1_ATAC_chr18.bam'

# optional
# sortBam(bamfile, "bam.sorted") # sort entries in bam file
# indexBam("bam.sorted.bam") # index reads in bam file

bam <- scanBam(bamfile)

# what did we get ?
class(bam)

#flag: how to store 0/1 information
#look at the first names in a list
#pos: chromosola position
# fields available

names(bam[[1]])

# let's look at first read in the alignment file
bam[[1]]$seq[1]

# let's count nt frequency
# for each sequence, calculate the number of nt withing each sequence
alphabetFrequency(bam[[1]]$seq[1])

# lets's plot GC content
#give fraction of GC in a set
gcFunction = function(x) { 
  alf <- alphabetFrequency(x, as.prob=TRUE)
  rowSums(alf[,c("G", "C")])
}
readGC <- gcFunction(bam[[1]][["seq"]])
hist(readGC)

#assess mapping quality of reads
# same scale as quality score
#quality score of 0: read didn't map well to the genome
# MAPQ
bam[[1]]$mapq

# histogram of mapq
hist(bam[[1]]$mapq)

# EXERCISE: calculate % mapped reads with qual below 20?
mean(bam[[1]]$mapq < 20, na.rm=T)*100

# let's identify covered regions (peaks)

# read in 
aln <- readGAlignments('LY1_ATAC_chr18.bam')

# what is this ?
class(aln)

# get coverage
# converts alignments into coverage
#coverage: for each nt, how many reads overlap each nt
cover <- coverage(aln)

# look at chr18 coverage
#first 10,000 reads are telomeres and have read of 0
cover$chr18

# identify islands of coverage >= 1 (view in IRange)
# regions of the genome above the threshhold
islands <- slice(cover, lower=1)

# identify max height in each island
# view maximum coverage
islandPeakHeight <- viewMaxs(islands)

# identify width of each island
islandWidth <- width(islands)

# find peaks 
# regions of coverage with atleast 100 bp and at least peak coverage 2
peaks <- islands[islandPeakHeight >= 2 & islandWidth >=100]

# look at them
peaks$chr18

# how many we have
length(peaks$chr18)
