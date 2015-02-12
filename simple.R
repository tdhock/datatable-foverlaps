works_with_R("3.1.2",
             GenomicRanges="1.18.4",
             data.table="1.9.4")

bedGraph <- data.table(chrom="chr1",
                       chromStart=as.integer(c(0, 200)),
                       chromEnd=as.integer(c(200, 300)),
                       coverage=as.integer(c(1, 2)))
bg.gr <- with(bedGraph, {
  GRanges(chrom, IRanges(chromStart, chromEnd), coverage=coverage)
})
bg.gr.1 <- with(bedGraph, {
  GRanges(chrom, IRanges(chromStart+1L, chromEnd), coverage=coverage)
})
bg.1 <- with(bedGraph, {
  data.table(chrom,
             chromStart=chromStart+1L,
             chromEnd,
             coverage)
})
test <- function(chromStart, chromEnd, expected){
  list(window=data.table(chrom="chr1",
         chromStart=as.integer(chromStart),
         chromEnd=as.integer(chromEnd)),
       expected=as.integer(expected))
}
tests <-
  list(test(200, 1000, 2),
       test(199, 1000, c(1, 2)),
       test(0, 200, 1),
       test(0, 201, c(1, 2)))
methods <-
  list("intersectBed-2.22.1"=function(win){
    write.table(win, "win.bed",
                quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)
    write.table(bedGraph, "bg.bedGraph",
                quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)
    cmd <- 
    paste("bedtools2/bin/intersectBed",
          "-sorted",
          "-wa -wb",
          "-a win.bed",
          "-b bg.bedGraph",
          "> overlap.bedGraph")
    system(cmd)
    overlap.df <- read.table("overlap.bedGraph")
    names(overlap.df) <-
      c("chrom.win", "chromStart.win", "chromEnd.win",
        "chrom.bg", "chromStart.bg", "chromEnd.bg", "coverage.bg")
    overlap.df$coverage.bg
  }, foverlaps=function(win){
    setkey(win, chrom, chromStart, chromEnd)
    overlap.dt <- foverlaps(bedGraph, win, nomatch=0L)
    overlap.dt$coverage
  }, findOverlaps=function(win){
    win.gr <- with(win, GRanges(chrom, IRanges(chromStart, chromEnd)))
    hits.gr <- findOverlaps(bg.gr, win.gr)
    overlap.df <- as.data.frame(bg.gr[queryHits(hits.gr), ])
    overlap.df$coverage
  }, "foverlaps+1"=function(win){
    win$chromStart <- win$chromStart+1L
    setkey(win, chrom, chromStart, chromEnd)
    overlap.dt <- foverlaps(bg.1, win, nomatch=0L)
    overlap.dt$coverage
  }, "findOverlaps+1"=function(win){
    win.gr <- with(win, GRanges(chrom, IRanges(chromStart+1L, chromEnd)))
    hits.gr <- findOverlaps(bg.gr.1, win.gr)
    overlap.df <- as.data.frame(bg.gr[queryHits(hits.gr), ])
    overlap.df$coverage
  })
result.list <- list()
for(test.i in seq_along(tests)){
  test.list <- tests[[test.i]]
  for(method in names(methods)){
    fun <- methods[[method]]
    computed <- fun(test.list$window)
    result <- all.equal(computed, test.list$expected)
    status <- ifelse(isTRUE(result), "correct", result)
    result.list[[paste(test.i, method)]] <-
      data.table(test.i, method, status,
                 expected=paste(test.list$expected, collapse=","),
                 computed=paste(computed, collapse=","),
                 chromStart=test.list$window$chromStart,
                 chromEnd=test.list$window$chromEnd)
  }
}
simple <- do.call(rbind, result.list)
save(simple, file="simple.RData")
