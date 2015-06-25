# QC suing Shortreads package

# install
#source("http://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
#biocLite("GenomicAlignments")

library(Rsamtools)
library(ShortRead)
library(lattice)
#setwd("TEACHING")

# load reads
fq <- readFastq("reads.fastq.gz")

# look at data
fq 

# how many reads
length(fq)

# read #50000 
fq[50000]

#show actual sequence
# get sequence
sread(fq[50000])

quality(fq[50000])

# get quality scores for 1 read
# convert quality score into number, and stores it in a matrix
quals<-as(head(quality(fq[50000])), "matrix")
# convert quality score into numbers
head(quality(fq[50000])
# look
quals

# translate into proba
#translate quality score into probabilities
#probability of not making the right call
10 ^ (-quals / 10)

# make boxplot
# same as: PhredQuality
#gives the box plot function for each column in the matrix
#boxplot: median
#outliers: avg +/- std deviation
#outline = F; remove the outliers
boxplot(as.matrix(PhredQuality(quality(fq))), outline=F, ylim=c(0,41), main="Per Cycle Read Quality", xlab="Cycle", ylab="Phred Quality", col=(c("gold","darkgreen")))

# plot nucleotide content
abc<-alphabetByCycle(sread(fq), alphabet=c("A","C","T","G"))
abc.df <- data.frame(Nucleotide=rownames(abc), Cycle=as.vector(col(abc)), Count=as.vector(abc))
print(xyplot(Count ~ Cycle, abc.df, group=Nucleotide, type="l", auto.key=list(lines=TRUE, points=FALSE))) 

# super fast way to get full report
qas <- qa("reads.fastq.gz")

# visualize it
browseURL(report(qas))

# trim reads from left using a sliding window (2 out 5 reads below "4" score)
#trim the reads from the right side
#how nt have a base score below 20, and look starting at least 2 nt; 4 is the asking code which is 20
fqt <- trimTailw(fq, 2, "4", halfwidth=2)

# check
fqt

# check new distribution of read length
hist(width(fqt))

# write trimmed reads to disk
#trim and write to new file
writeFastq(fq, "treads.fastq.gz", "w")


