(objs <- load("full.strand.list.RData"))
strand <- "both"
chipseq.dt <- full.strand.list[[strand]]
dim(chipseq.dt)
# [1] 15126283        4
chipseq.df <- data.frame(chipseq.dt)
windows.df <- data.frame(windows)

# R version 3.1.2
require(data.table)    # latest commit from github
require(GenomicRanges) # "1.18.4"
GR_fun <- function(query, windows){
  windows.gr <- with(windows, GRanges(chrom, IRanges(chromStart, chromEnd)))
  query.gr <- with(query, {
    GRanges(chrom, IRanges(chromStart, chromEnd), coverage=query$count)
  })
  findOverlaps(query.gr, windows.gr)
}

## Note that Arun's timings use a newer version of
## GenomicRanges=1.18.4 which is slower than the old version that I
## was using in
## http://www.bioconductor.org/packages/2.14/BiocViews.html#___Software
system.time(foverlaps(chipseq.dt, windows, nomatch=0L, which=TRUE))
#    user  system elapsed  
#   1.020   0.049   1.079  

system.time(GR_fun(chipseq.df, windows.df))
#    user  system elapsed  
#   2.351   0.398   2.752

