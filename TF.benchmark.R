works_with_R("3.1.2",
             GenomicRanges="1.18.4",
             dplyr="0.4.0",
             microbenchmark="1.3.0",
             ##"Rdatatable/data.table@84ba1151299ba49e833e68a2436630216b653306",
             data.table="1.9.4")

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

## This is the largest data set that we consider.
tf.name <- "nrsf"
bg.file <- "/home/thocking/genomecov/nrsf/huds_k562_rep2.bedGraph"
strand <- "+"

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
  ##str(windows)
  windows[, chrom := as.character(chrom) ]
  ##str(windows)
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


    full.strand.list <- list()
    
    ## bg.list <- read.bench(bg.file)
    ## full.strand.list$both <- bg.list$fread
    full.strand.times <- list()
    ## full.strand.times$both <-
    ##   data.table(strand="both",
    ##              rows=nrow(bg.list$fread),
    ##              bg.list$times)
    
    strand.files <- list(both=bg.file)
    ##for(strand in c("+", "-")){
      plus.file <- strand.files[[strand]] <-
        paste0(bg.file, strand, "strand")
      strand.list <- read.bench(plus.file)
      full.strand.list[[strand]] <- strand.list$fread
      full.strand.times[[strand]] <-
        data.table(strand, rows=nrow(strand.list$fread), strand.list$times)
    ##}
    sample.times <- do.call(rbind, full.strand.times) %>%
      mutate(seconds=time/1e9)
    read.times.list[[bg.file]] <-
      data.table(sample.id, tf.name, experiment, sample.times)

    ##save(full.strand.list, windows, file="full.strand.list.RData")

    window.strand.list <- list()
    ##for(strand in names(full.strand.list)){
      one.strand <- full.strand.list[[strand]]
      big.df <- data.frame(one.strand)
      small.df <- data.frame(windows)
      times.tables <- microbenchmark(`GenomicRanges::findOverlaps`={
        big.gr <- with(big.df, {
          GRanges(chrom, IRanges(chromStart, chromEnd), coverage=count)
        })
        small.gr <- with(small.df, {
          GRanges(chrom, IRanges(chromStart, chromEnd))
        })
        hits.gr <- findOverlaps(big.gr, small.gr)
        hits.big <- as.data.frame(big.gr[queryHits(hits.gr), ])
        hits.small <- as.data.frame(small.gr[subjectHits(hits.gr), ])
      }, `data.table::foverlaps`={
        big.dt <- data.table(big.df)
        setkey(big.dt, chrom, chromStart, chromEnd)
        small.dt <- data.table(small.df)
        setkey(small.dt, chrom, chromStart, chromEnd)
        one.join <- foverlaps(big.dt, small.dt, nomatch=0L)
      }, times=2)
      stopifnot(all.equal(as.character(one.join$chrom),
                          as.character(hits.big$seqnames)))
      stopifnot(all.equal(as.character(one.join$chrom),
                          as.character(hits.small$seqnames)))
      stopifnot(all.equal(one.join$chromStart, hits.small$start))
      stopifnot(all.equal(one.join$chromEnd, hits.small$end))
      stopifnot(all.equal(one.join$count, hits.big$coverage))
      stopifnot(all.equal(one.join$i.chromStart, hits.big$start))
      stopifnot(all.equal(one.join$i.chromEnd, hits.big$end))

      times.indices <- microbenchmark(`GenomicRanges::findOverlaps`={
        hits.gr <- findOverlaps(big.gr, small.gr)
      }, `data.table::foverlaps`={
        hits.dt <- foverlaps(big.dt, small.dt, nomatch=0L, which=TRUE)
      }, times=2)
      stopifnot(all.equal(queryHits(hits.gr), hits.dt$xid))
      stopifnot(all.equal(subjectHits(hits.gr), hits.dt$yid))
      
      windows.bed <- small.df[, c("chrom", "chromStart", "chromEnd")]
      write.table(windows.bed, "windows.bed",
                  quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)

      strand.file <- strand.files[[strand]]
      bedtools <- function(intersectBed){
        cmd <- 
          paste(intersectBed,
                "-wa -wb -a windows.bed -b",
                strand.file,
                "> overlap.bedGraph")
        system(cmd)
      }
      R.cmd <-
        paste("R --no-save --args",
              strand.file,
              "windows.bed overlap-startup.bedGraph < intersect.R")
      
      times.IO <- microbenchmark(`intersectBed-2.14.3`={
        bedtools("/usr/bin/intersectBed")
      }, `intersectBed-2.17.0`={
        bedtools("bedtools-2.17.0/bin/intersectBed -sorted")
      }, `intersectBed-2.22.1`={
        bedtools("bedtools2/bin/intersectBed -sorted")
      }, startup.fread.foverlaps.write={
        system(R.cmd)
      }, fread.foverlaps.write={
        ## Is intersectBed slower simply because it needs to read/write
        ## the files from/to disk? Try read/write in R to compare.
        one.strand <- fread(strand.file)
        setnames(one.strand, c("chrom", "chromStart", "chromEnd", "coverage"))
        setkey(one.strand, chrom, chromStart, chromEnd)
        windows <- fread("windows.bed")
        setnames(windows, c("chrom", "chromStart", "chromEnd"))
        setkey(windows, chrom, chromStart, chromEnd)
        one.join <- foverlaps(one.strand, windows, nomatch=0L)
        write.table(one.join, file="overlap-R.bedGraph",
                    quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)
      }, times=2)
      system("wc -l overlap-R.bedGraph overlap-startup.bedGraph overlap.bedGraph")

      meta <- 
        data.table(sample.id, tf.name, experiment, strand,
                   query.rows=nrow(one.strand),
                   subject.rows=nrow(windows),
                   overlap.rows=nrow(one.join))
      overlap.times.list[[paste(bg.file, strand)]] <-
        rbind(data.table(meta, what="tables", times.tables),
              data.table(meta, what="indices", times.indices),
              data.table(meta, what="IO", times.IO))
    ##}#strand
  }#bg.file.i
}#tf.name


TF.benchmark <-
  list(overlap=do.call(rbind, overlap.times.list),
       read=do.call(rbind, read.times.list))

save(TF.benchmark, file="TF.benchmark.RData")
