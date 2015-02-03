works_with_R("3.1.2",
             GenomicRanges="1.16.4",
             dplyr="0.4.0",
             microbenchmark="1.3.0",
             "Rdatatable/data.table@84ba1151299ba49e833e68a2436630216b653306")

downloads <- "http://tare.medisin.ntnu.no/chipseqbenchmark/downloads/"

colClasses <-
  c("chrom"="factor",
    "chromStart"="integer",
    "chromEnd"="integer",
    "count"="integer")

## Use the *nix wc program to quickly determine the number of lines
## of a file.
wc <- function(f){
  stopifnot(is.character(f))
  stopifnot(length(f)==1)
  if(file.exists(f)){
    cmd <- sprintf("wc -l '%s'",f)
    as.integer(sub(" .*","",system(cmd,intern=TRUE)))
  }else{
    0L
  }
}

read.bench <- function(filename){
  times <- microbenchmark(fread={
    dt <- fread(filename)
    setnames(dt, names(colClasses))
    setkey(dt, chrom, chromStart, chromEnd)
  }, read.table={
    df <-
      read.table(filename,
                 sep="\t",
                 colClasses = as.character(colClasses),
                 nrows = wc(filename))
  }, times=2)
  stopifnot(nrow(dt) == nrow(df))
  list(times=times, fread=dt, read.table=df)
}

overlap.times.list <- read.times.list <-
  site.list <- region.list <- coverage.list <- list()
for(tf.name in c("max", "nrsf", "srf")){

  ## First read manually annotated peak regions.
  ann.file <-
    paste0("Manually_classified_peaks/classified_peaks_", tf.name, ".txt")
  ann.url <- paste0(downloads, ann.file)
  nrsf.peaks <- read.table(ann.url)
  bg.names <- if(tf.name == "max"){
    c("background", "balanced", "replicate")
  }else{
    c("background1",
      "background2",
      "replicate",
      "subset")    
  }
  names(nrsf.peaks) <-
    c("chrom",
      "chromStart",
      "chromEnd",
      "bases",
      "region.index",
      "initial",
      bg.names,
      "overall")

  ## Also read TF binding site annotations.
  sites.file <-
    paste0("Manually_classified_sites/classified_sites_", tf.name, ".txt")
  sites.url <- paste0(downloads, sites.file)
  nrsf.sites <- read.table(sites.url)
  names(nrsf.sites) <-
    c("chrom",
      "chromStart",
      "chromEnd",
      "pwm.score",
      "tags.per.bp.within.site.sample",
      if(tf.name != "max")"tags.per.bp.within.site.replicate" else NULL,
      "region.index",
      "classification")

  ann.code <-
    c("-1"="noPeaks",
      "0"="ambiguous",
      "1"="onePeak")

  pos.neg <- data.table(tf.name, nrsf.peaks) %>%
    mutate(annotation=ann.code[as.character(overall)])

  ## Windows to show.
  window.cols <- c("region.index", "chrom", "chromStart", "chromEnd")
  sites.and.peaks <-
    rbind(data.table(nrsf.sites[, window.cols], what="site"),
          data.table(nrsf.peaks[, window.cols], what="region"))
  expand.bases <- 3000L
  windows <- sites.and.peaks %>%
    group_by(region.index, chrom) %>%
    summarise(chromStart=min(chromStart)-expand.bases,
              chromEnd=max(chromEnd)+expand.bases,
              sites=sum(what=="site"),
              regions=sum(what=="region"))
  str(windows)
  windows[, chrom := as.character(chrom) ]
  str(windows)
  setkey(windows, chrom, chromStart, chromEnd)

  ## Database of visual annotations.
  pos.neg.regions <- pos.neg %>%
    select(tf.name, region.index, chrom, chromStart, chromEnd, annotation)

  tf.dir <- file.path("/home/thocking/genomecov", tf.name)
  tf.files <- Sys.glob(file.path(tf.dir, "*.bedGraph"))
  bg.files <- sub(tf.name, "backgr", tf.files)
  bedgraph.files <- c(tf.files, bg.files[file.exists(bg.files)])
  sample.list <- list()
  for(bg.file.i in seq_along(bedgraph.files)){
    bg.file <- bedgraph.files[[bg.file.i]]
    sample.id <- sub(".bedGraph", "", basename(bg.file))
    experiment <- basename(dirname(bg.file))
    cat(sprintf("%4d / %4d %s %s\n",
                bg.file.i, length(bedgraph.files),
                sample.id, experiment))

    bg.list <- read.bench(bg.file)

    full.strand.list <- list(both=bg.list$fread)
    full.strand.times <-
      list(both=data.table(strand="both",
             rows=nrow(bg.list$fread), bg.list$times))
    for(strand in c("+", "-")){
      plus.file <- paste0(bg.file, strand, "strand")
      strand.list <- read.bench(plus.file)
      full.strand.list[[strand]] <- strand.list$fread
      full.strand.times[[strand]] <-
        data.table(strand, rows=nrow(strand.list$fread), strand.list$times)
    }
    sample.times <- do.call(rbind, full.strand.times) %>%
      mutate(seconds=time/1e9)
    read.times.list[[bg.file]] <-
      data.table(sample.id, tf.name, experiment, sample.times)

    ##save(full.strand.list, windows, file="full.strand.list.RData")

    window.strand.list <- list()
    for(strand in names(full.strand.list)){
      one.strand <- full.strand.list[[strand]]
      big.gr <- with(one.strand, GRanges(chrom, IRanges(chromStart, chromEnd)))
      small.gr <- with(windows, GRanges(chrom, IRanges(chromStart, chromEnd)))
      times.indices <- microbenchmark(`GenomicRanges::findOverlaps`={
        hits.gr <- findOverlaps(big.gr, small.gr)
      }, `data.table::foverlaps`={
        hits.dt <- foverlaps(one.strand, windows, nomatch=0L, which=TRUE)
      }, times=2)
      stopifnot(length(hits.gr) == nrow(hits.dt))
      stopifnot(queryHits(hits.gr) == hits.dt$xid)
      stopifnot(subjectHits(hits.gr) == hits.dt$yid)
      
      times.tables <- microbenchmark(`GenomicRanges::findOverlaps`={
        big.gr <- with(one.strand, {
          GRanges(chrom, IRanges(chromStart, chromEnd))
        })
        small.gr <- with(windows, {
          GRanges(chrom, IRanges(chromStart, chromEnd))
        })
        hits.gr <- findOverlaps(big.gr, small.gr)
        df.join <- 
          cbind(as.data.frame(big.gr[queryHits(hits.gr), ]),
                as.data.frame(small.gr[subjectHits(hits.gr), ]))
      }, `data.table::foverlaps`={
        one.join <- foverlaps(one.strand, windows, nomatch=0L)
      }, times=2)
      stopifnot(nrow(one.join) == nrow(df.join))
      stopifnot(one.join$chrom == df.join$seqnames)
      stopifnot(one.join$chromStart == df.join$chromStart)
      stopifnot(one.join$chromEnd == df.join$chromEnd)

      write.table(one.strand, "chipseq.bedGraph",
                  quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)
      windows.bed <- data.frame(windows)[, c("chrom", "chromStart", "chromEnd")]
      write.table(windows.bed, "windows.bed",
                  quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)
      cmd <- "intersectBed -wa -wb -a windows.bed -b chipseq.bedGraph > overlap.bedGraph"
      times.IO <- microbenchmark(intersectBed={
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
      }, times=2)
      system("wc -l overlap-R.bedGraph overlap.bedGraph")

      meta <- 
        data.table(sample.id, tf.name, experiment, strand,
                   query.rows=nrow(one.strand),
                   subject.rows=nrow(windows),
                   overlap.rows=nrow(one.join))
      overlap.times.list[[paste(bg.file, strand)]] <-
        rbind(data.table(meta, what="indices", times.tables),
              data.table(meta, what="tables", times.indices),
              data.table(meta, what="IO", times.IO))
    }
  }
}

TF.benchmark <-
  list(overlap=do.call(rbind, overlap.times.list),
       read=do.call(rbind, read.times.list))

save(TF.benchmark, file="TF.benchmark.RData")
