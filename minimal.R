## 1. <link to data>
RData.url <- "http://cbio.ensmp.fr/~thocking/data/full.strand.list.RData"

## 2. <code to load data from link>
if(!file.exists("full.strand.list.RData")){
  download.file(RData.url, "full.strand.list.RData")
}
(objs <- load("full.strand.list.RData"))
strand <- "both"
chipseq.dt <- full.strand.list[[strand]]
chipseq.df <- data.frame(chipseq.dt)
windows.df <- data.frame(windows)

## 3. GR_fun <- { … } # self contained function that creates ranges from the data and runs findOverlaps.
require(GenomicRanges)#"1.16.4"
GR_fun <- function(query, windows){
  windows.gr <- with(windows, GRanges(chrom, IRanges(chromStart, chromEnd)))
  query.gr <- with(query, {
    GRanges(chrom, IRanges(chromStart, chromEnd), coverage=query$count)
  })
  findOverlaps(query.gr, windows.gr)
}

## 4. DT_fun <- { … } # self contained function that runs foverlaps (or use setDT to convert from data.frame to data.table).
require(data.table) ###"Rdatatable/data.table@84ba1151299ba49e833e68a2436630216b653306"
DT_fun <- function(query, windows){
  query.dt <- data.table(query)
  ##setkey(query.dt, chrom, chromStart, chromEnd)
  windows.dt <- data.table(windows)
  setkey(windows.dt, chrom, chromStart, chromEnd)
  foverlaps(query.dt, windows.dt, nomatch=0L, which=TRUE)
}

## 5. so that I can call them as: GR_fun(query1, window); DT_fun(query1, window) etc..
library(microbenchmark)
microbenchmark(data.table={
  foverlaps(chipseq.dt, windows, nomatch=0L, which=TRUE)
}, data.table.setkey=DT_fun(chipseq.df, windows.df),
  GenomicRanges=GR_fun(chipseq.df, windows.df),
  times=2)
