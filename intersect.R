works_with_R("3.1.2",
             "Rdatatable/data.table@84ba1151299ba49e833e68a2436630216b653306")

argv <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argv) == 3)
stopifnot(is.character(argv))
big.file <- argv[1]
small.file <- argv[2]
out.file <- argv[3]

big <- fread(big.file)
big.names <- names(big)
big.names[1:3] <- c("chrom", "chromStart", "chromEnd")
setnames(big, big.names)
setkey(big, chrom, chromStart, chromEnd)

small <- fread(small.file)
small.names <- names(small)
small.names[1:3] <- c("chrom", "chromStart", "chromEnd")
setnames(small, small.names)
setkey(small, chrom, chromStart, chromEnd)

overlap <- foverlaps(big, small, nomatch=0L)
write.table(overlap, out.file,
            quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)

