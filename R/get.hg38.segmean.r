# get.hg38.segmean.r
# This is a modifed version GetGeneLevelCNA.r from
# https://github.com/ZhenyuZ/eqtl/blob/master/cnv/GetGeneLevelCNA.r
#
# Input GDC segmentation data, window size and calculate segment mean
# in genomic regions every window.size bases 

# split chromosome by window.size to create region
# input "chromosome" contains two columns in order: chr and len
# window.size is any interger > 1
make.region = function(chromosomes, window.size) {
  chromosome$nseg = ceiling(chromosome$len / window.size)
  region = data.frame(chr = unlist(apply(chromosome, 1, function(x) rep(x[1], x[3]))), 
            start = unlist(apply(chromosome, 1, function(x) (seq(as.numeric(x[3]))-1)*window.size + 1)), 
            end = unlist(apply(chromosome, 1, function(x) c((seq(as.numeric(x[3]))-1)*window.size, as.numeric(x[2]))[-1]))
    )
  region$name = with(region, paste0(chr, ":", start, "-", end))
  return(region)
}


# read segmentation file and convert to GRanges object (only autosomes)
seg2gr = function(file){
  seg = read.table(file, h=T, stringsAsFactors=F)
  # load into granges object
  gr = with(seg[seg$Chromosome %in% 1:22, ], GRanges( seqnames = paste0("chr", Chromosome),
                                                    ranges = IRanges(start = Start, end = End),
                                                    Segment_Mean = Segment_Mean))
  return(gr)
}

# calculate mean of segment_mean in gr by region 
# region is a data frame containing the following 4 columns in order: chr, start, end, name
get.segmean = function(gr, region) {
  window.mean = NULL
  for (chr in chromosome$chr) {
    my.gr = gr[seqnames(gr) == chr]
    my.region = region[region$chr == chr, ]
    window.mean = c(window.mean, apply(my.region, 1, 
                  function(x){y = restrict(my.gr, start  = as.numeric(x[2]), end = as.numeric(x[3])); 
                  weighted.mean(y$Segment_Mean, width(y))})  )
  }
  names(window.mean) = region$name
  return(window.mean)
}


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
segfile = args[1]
window.size = as.numeric(args[2])
outfile = args[3]


# load library and init
options(stringsAsFactors=F)
options(scipen=999)

library("GenomicRanges")

# reconstruct autosomal chromosome information for GRCh38
chromosome = data.frame(chr=paste0("chr", 1:22), len=c(248956422, 242193529, 198295559, 190214555, 181538259, 
              170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 
              101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468))

# subset by window.size
region = make.region(chromosome, window.size)

# read segment file
gr = seg2gr(segfile)

# calculate segment mean 
segmean = get.segmean(gr, region)

# write.table
write.table(segmean, outfile, col.names=F, row.names=T, sep="\t", quote=F)


# transform and collect output together in bash
# ls | grep TCGA > f
# while read -r line; do awk -v line=$line 'BEGIN {printf line}; {printf "\t"$2}; END{printf "\n"}' $line >> aggregate;  done <"f"


library(data.table)

aggregate.file = "aggregate"
aggregate = fread(aggregate.file, h=F)

data = t(as.matrix(aggregate[, 2: ncol(aggregate), with=F]))
colnames(data) = aggregate$V1
rownames(data) = fread(sample.file, h=F)$V1

write.table(data, "TCGA.cnv.window_1M.txt", col.names=T, row.names=T, sep="\t", quote=F)


















 