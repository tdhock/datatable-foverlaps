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

all.diff.list <- list()
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
  cmd <- "intersectBed -wa -wb -a windows.bed -b chipseq.bedGraph > overlap.bedGraph"

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
  }, fread.foverlaps.write={
    ## Is intersectBed slower simply because it needs to read/write
    ## the files from/to disk? Try read/write in R to compare.
    CSR <- fread("chipseq.bedGraph")
    setnames(CSR, c("chrom", "chromStart", "chromEnd", "coverage"))
    setkey(CSR, chrom, chromStart, chromEnd)
    WR <- fread("windows.bed")
    setnames(WR, c("chrom", "chromStart", "chromEnd"))
    setkey(WR, chrom, chromStart, chromEnd)
    ODT <- foverlaps(CSR, WR, nomatch=0L)
    write.table(ODT, file="overlap-R.bedGraph",
                quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)
  },times=5)
  print(times)
  time.list[[strand]] <-
    data.table(times,
               chipseq.rows=nrow(chipseq),
               windows.rows=nrow(windows),
               overlap.rows=nrow(overlap.chipseq),
               strand)

  ## The following two code blocks ensure that the results of
  ## data.table::foverlaps and GenomicRanges::findOverlaps are
  ## identical.
  stopifnot(nrow(overlap.dt) == nrow(overlap.chipseq))
  stopifnot(overlap.dt$chrom == overlap.chipseq$seqnames)
  stopifnot(overlap.dt$chromStart.i == overlap.chipseq$start)
  stopifnot(overlap.dt$chromEnd.i == overlap.chipseq$end)
  
  stopifnot(nrow(overlap.dt) == nrow(overlap.windows))
  stopifnot(overlap.dt$chrom == overlap.windows$seqnames)
  stopifnot(overlap.dt$chromStart == overlap.windows$start)
  stopifnot(overlap.dt$chromEnd == overlap.windows$end)

  ## Also check consistency of R and intersectBed.
  intersectBed <- fread("overlap.bedGraph")
  bed.names <-
    c("chrom", "windowStart", "windowEnd",
      "chrom2", "chipseqStart", "chipseqEnd", "coverage")
  setnames(intersectBed, bed.names)
  data.list <-
    list(R=ODT %>%
         mutate(start=i.chromStart,
                end=i.chromEnd,
                method="R",
                region=sprintf("%s:%d-%d", chrom, chromStart, chromEnd)) %>%
         select(start, end, method, region, coverage),
         intersectBed=intersectBed %>%
         mutate(start=chipseqStart,
                end=chipseqEnd,
                method="intersectBed",
                region=sprintf("%s:%d-%d", chrom, windowStart, windowEnd)) %>%
         select(start, end, method, region, coverage))
  region.list <- list()
  for(method in names(data.list)){
    dt <- data.list[[method]]
    method.list <- split(dt, dt$region)
    for(region.name in names(method.list)){
      region.list[[region.name]][[method]] <- method.list[[region.name]]
    }
  }
  diff.list <- list()
  for(region.name in names(region.list)){
    dt.list <- region.list[[region.name]]
    dt.rows <- sapply(dt.list, nrow)
    size.1 <- nrow(dt.list[[1]])
    size.2 <- nrow(dt.list[[2]])
    criteria <-
      c(size=size.1 - size.2,
        starts=sum(dt.list[[1]]$start != dt.list[[2]]$start),
        ends=sum(dt.list[[1]]$end != dt.list[[2]]$end),
        coverage=sum(dt.list[[1]]$coverage != dt.list[[2]]$coverage))
    category <- if(sum(criteria) == 0){
      "identical"
    }else{
      if(criteria[["size"]] == 0){
        shifted <-
          c(sum(dt.list[[1]]$start[-1] != dt.list[[2]]$start[-size.2]),
            sum(dt.list[[1]]$start[-size.1] != dt.list[[2]]$start[-1]))
        if(any(shifted == 0)){
          "samesize.shift1"
        }else{
          diff.list[[region.name]] <- region.list[[region.name]]
          coverage <- do.call(rbind, dt.list)
          ggplot()+
            geom_step(aes(start/1e3, coverage),
                      data=coverage, color="grey")+
            theme_bw()+
            theme(panel.margin=grid::unit(0, "cm"))+
            facet_grid(method ~ .)
          "samesize.shift>1"
        }
      }
    }
    if(is.null(category)){
      category <- "diffsize"
    }
    all.diff.list[[paste(strand, region.name)]] <- 
      data.table(strand, region.name, category)
  }
  
}

differences <- do.call(rbind, all.diff.list)

full.strand.times <- do.call(rbind, time.list)

save(full.strand.times, differences, file="full.strand.times.RData")
