works_with_R("3.1.2",
             GenomicRanges="1.16.4",
             dplyr="0.4.0",
             microbenchmark="1.3.0",
             "Rdatatable/data.table@84ba1151299ba49e833e68a2436630216b653306")

RData.url <- "http://cbio.ensmp.fr/~thocking/data/full.strand.list.RData"
if(!file.exists("full.strand.list.RData")){
  download.file(RData.url, "full.strand.list.RData")
}

load("full.strand.list.RData")

time.list <- list()
for(strand in names(full.strand.list)){
  chipseq <- full.strand.list[[strand]]
  windows.gr <- with(windows, GRanges(chrom, IRanges(chromStart, chromEnd)))
  chipseq.gr <- with(chipseq, {
    GRanges(chrom, IRanges(chromStart, chromEnd), coverage=chipseq$count)
  })
  ## times <- microbenchmark(`GenomicRanges::findOverlaps`={
  ##   hits.gr <- findOverlaps(chipseq.gr, windows.gr)
  ##   overlap.chipseq <- as.data.frame(chipseq.gr[queryHits(hits.gr), ])
  ##   overlap.windows <- as.data.frame(windows.gr[subjectHits(hits.gr), ])
  ## }, `data.table::foverlaps`={
  ##   overlap.dt <- foverlaps(chipseq, windows, nomatch=0L)
  ## }, times=5)
  ## print(times)

  ## stopifnot(nrow(overlap.dt) == nrow(overlap.chipseq))
  ## stopifnot(overlap.dt$chrom == overlap.chipseq$seqnames)
  ## stopifnot(overlap.dt$chromStart.i == overlap.chipseq$start)
  ## stopifnot(overlap.dt$chromEnd.i == overlap.chipseq$end)
  
  ## stopifnot(nrow(overlap.dt) == nrow(overlap.windows))
  ## stopifnot(overlap.dt$chrom == overlap.windows$seqnames)
  ## stopifnot(overlap.dt$chromStart == overlap.windows$start)
  ## stopifnot(overlap.dt$chromEnd == overlap.windows$end)

  write.table(chipseq, "chipseq.bedGraph",
              quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)
  windows.bed <- data.frame(windows)[, c("chrom", "chromStart", "chromEnd")]
  write.table(windows.bed, "windows.bed",
              quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)
  cmd <- "intersectBed -wb -a windows.bed -b chipseq.bedGraph > overlap.bedGraph"

  chipseq.read <- fread("chipseq.bedGraph")
  setnames(chipseq.read, c("chrom", "chromStart", "chromEnd", "coverage"))
  setkey(chipseq.read, chrom, chromStart, chromEnd)
  stopifnot(nrow(chipseq.read) == nrow(chipseq))
  stopifnot(chipseq$chrom == chipseq.read$chrom)
  stopifnot(chipseq$chromStart == chipseq.read$chromStart)
  stopifnot(chipseq$chromEnd == chipseq.read$chromEnd)

  windows.read <- fread("windows.bed")
  setnames(windows.read, c("chrom", "chromStart", "chromEnd"))
  setkey(windows.read, chrom, chromStart, chromEnd)
  stopifnot(nrow(windows.read) == nrow(windows))
  stopifnot(windows$chrom == windows.read$chrom)
  stopifnot(windows$chromStart == windows.read$chromStart)
  stopifnot(windows$chromEnd == windows.read$chromEnd)

  chipseq.read.gr <- with(chipseq.read, {
    GRanges(chrom, IRanges(chromStart, chromEnd), coverage=coverage)
  })
  windows.read.gr <- with(windows.read, {
    GRanges(chrom, IRanges(chromStart, chromEnd))
  })
  ## .windows.read means that we use the data read from the
  ## windows.bed file, which for some reason speeds up foverlaps!
  times <- microbenchmark(`GenomicRanges::findOverlaps`={
    hits.gr <- findOverlaps(chipseq.gr, windows.gr)
    overlap.chipseq <- as.data.frame(chipseq.gr[queryHits(hits.gr), ])
    overlap.windows <- as.data.frame(windows.gr[subjectHits(hits.gr), ])
  }, `data.table::foverlaps`={
    overlap.dt <- foverlaps(chipseq, windows, nomatch=0L)
  }, `GenomicRanges::findOverlaps.windows.read`={
    hits.gr <- findOverlaps(chipseq.gr, windows.read.gr)
    overlap.chipseq.read <- as.data.frame(chipseq.read.gr[queryHits(hits.gr), ])
    overlap.windows.read <- as.data.frame(windows.gr[subjectHits(hits.gr), ])
  }, `data.table::foverlaps.windows.read`={
    overlap.dt <- foverlaps(chipseq, windows.read, nomatch=0L)
  }, intersectBed={
    system(cmd)
  }, times=5)
  print(times)
  
  stopifnot(nrow(overlap.dt) == nrow(overlap.chipseq))
  stopifnot(overlap.dt$chrom == overlap.chipseq$seqnames)
  stopifnot(overlap.dt$chromStart.i == overlap.chipseq$start)
  stopifnot(overlap.dt$chromEnd.i == overlap.chipseq$end)
  
  stopifnot(nrow(overlap.dt) == nrow(overlap.windows))
  stopifnot(overlap.dt$chrom == overlap.windows$seqnames)
  stopifnot(overlap.dt$chromStart == overlap.windows$start)
  stopifnot(overlap.dt$chromEnd == overlap.windows$end)
}
