works_with_R("3.1.2",
             GenomicRanges="1.16.4",
             dplyr="0.4.0",
             microbenchmark="1.3.0",
             "Rdatatable/data.table@84ba1151299ba49e833e68a2436630216b653306")

tgz.url <- "http://cbio.ensmp.fr/~thocking/data/overlap-benchmark.tgz"
if(!file.exists(tgz <- "overlap-benchmark.tgz")){
  download.file(tgz.url, tgz)
  system(paste("tar xf", tgz))
}
windows <- fread("overlap-benchmark/windows.bed")
setnames(windows, c("chrom", "chromStart", "chromEnd"))
chip.seq <- fread("overlap-benchmark/chip-seq.bedGraph")
setnames(chip.seq, c("chrom", "chromStart", "chromEnd", "coverage"))
windows.gr <- with(windows, GRanges(chrom, IRanges(chromStart, chromEnd)))
chip.seq.gr <- with(chip.seq, {
  GRanges(chrom, IRanges(chromStart, chromEnd))
})
setkey(windows, chrom, chromStart, chromEnd)
setkey(chip.seq, chrom, chromStart, chromEnd)
windows.extra <- windows %>%
  mutate(sites=rnorm(n()),
         regions=rnorm(n()))
times <- microbenchmark(`GenomicRanges::findOverlaps`={
  hits.gr <- findOverlaps(chip.seq.gr, windows.gr)
}, `data.table::foverlaps`={
  hits.dt <- foverlaps(chip.seq, windows, nomatch=0L, which=TRUE)
}, `data.table::foverlaps.extra`={
  extra.dt <- foverlaps(chip.seq, windows.extra, nomatch=0L, which=TRUE)
}, times=5)
print(times)

load("TF.benchmark.RData")
TF.benchmark$overlap %>%
  filter(sample.id == "yale_k562_rep1",
         experiment=="max",
         strand=="-")
