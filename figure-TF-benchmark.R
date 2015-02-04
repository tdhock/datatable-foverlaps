works_with_R("3.1.2",
             dplyr="0.4.0",
             "hadley/tidyr@3de52c46f12d0d4cba603aa4e63e39f2370a9cfd",
             ##"Rdatatable/data.table@84ba1151299ba49e833e68a2436630216b653306",
             directlabels="2014.6.13",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7")

load("TF.benchmark.RData")

algo.colors <-
  c("intersectBed-2.22.1\ngithub"="#f8091f",
    "intersectBed-2.14.3\nubuntu"="#f87b87",
    "intersectBed-2.17.0\nguillimin"="#f8b3ba",
    "data.table::foverlaps"="#f37036",
    "fread.foverlaps.write"="#f37036",
    "startup.fread.foverlaps.write"="black",
    "data.table::fread"="#f37036",
    "read.table"="black",
    "GenomicRanges::findOverlaps"="#1892aa")

over.end <- list(dl.trans(x=x+0.1), "last.qp")

## from https://github.com/Rdatatable/data.table/wiki/talks/EARL2014_OverlapRangeJoin_Arun.pdf
minutes.wide <-
  data.frame(query.rows=c(80e6, 65e6),
             subject.rows=c(33e3, 36e3),
             `data.table::foverlaps`=c(2, 4),
             `GenomicRanges::findOverlaps`=c(16, 28* 24),
             check.names=FALSE)
minutes.tall <- 
gather(minutes.wide, method, minutes, -c(query.rows, subject.rows))

refs <-
  data.table(unit=c("1 second"),
             seconds=c(1),
             vjust=c(-0.5))

bed.labs <-
  c("intersectBed-2.17.0"="intersectBed-2.17.0\nguillimin",
    "intersectBed-2.14.3"="intersectBed-2.14.3\nubuntu",
    "intersectBed-2.22.1"="intersectBed-2.22.1\ngithub")
ov <- TF.benchmark$overlap %>%
  mutate(method=ifelse(expr %in% names(bed.labs),
           bed.labs[paste(expr)],
           paste(expr)),
         seconds=time/1e9) %>%
  filter(seconds < 50)

overlab.df <-
  data.table(rows=3.5e7,
             method=c("GenomicRanges::findOverlaps", "data.table::foverlaps"),
             seconds=c(6.5, 15))

with.labels <- 
  ggplot()+
  ## geom_text(aes(rows, seconds, label=method, color=method),
  ##           data=overlab.df, size=3)+
  ##geom_hline(aes(yintercept=seconds), data=refs, color="grey50")+
  ## geom_text(aes(3.75e7, seconds, label=unit, vjust=vjust),
  ##           data=refs, color="grey50", size=3)+
  geom_point(aes(query.rows, seconds, color=method),
             data=ov, pch=1)+
  xlab("rows in bedGraph file")+
  guides(color="none")+
  theme_bw()+
  scale_y_continuous("seconds")+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ ., scales="free")
direct.label(with.labels, "over.end")

with.labels+
  geom_point(aes(query.rows, minutes*60, color=method),
             data=minutes.tall)+
  scale_y_log10()+
  scale_x_log10()

ov %>%
  filter(query.rows==max(query.rows))

ov.list <- split(ov, ov$what)
for(comparison in names(ov.list)){
  o <- ov.list[[comparison]]
  rows.lab <- mean(sort(unique(o$query.rows), decreasing=TRUE)[1:2])
  expr.list <- split(o, o$expr, drop=TRUE)
  lab.list <- list()
  for(expr in names(expr.list)){
    seconds <- with(expr.list[[expr]], {
      approx(query.rows, seconds, rows.lab)$y
    })
    lab.list[[expr]] <- data.table(method=expr, rows=rows.lab, seconds)
  }
  overlab.df <- do.call(rbind, lab.list)

  no.labels <- 
  ggplot()+
  ##geom_hline(aes(yintercept=seconds), data=refs, color="grey50")+
  ## geom_text(aes(3.75e7, seconds, label=unit, vjust=vjust),
  ##           data=refs, color="grey50", size=3)+
  geom_point(aes(query.rows, seconds, color=method),
             data=o, pch=1)+
  scale_x_continuous("rows in bedGraph file", limits=c(NA, 4e7),
                     breaks=c(1, 2)*1e7)+
  guides(color="none")+
  theme_grey()+
    scale_color_manual(values=algo.colors)+
  scale_y_continuous("seconds")

  manual.labels <-
    no.labels +
    geom_text(aes(rows, seconds, label=method, color=method),
              data=overlab.df, size=3)

  dl <- direct.label(no.labels, "over.end")

  pdf.name <- sprintf("figure-TF-benchmark-%s.pdf", comparison)
  pdf(pdf.name, 5, 3)
  print(dl)
  dev.off()

  dlog <-
    dl+
      scale_x_log10("rows in bedGraph file", limits=c(NA, 4e7))+
      scale_y_log10()
  log.pdf.name <- sprintf("figure-TF-benchmark-%s-log.pdf", comparison)
  pdf(log.pdf.name, 5, 3)
  print(dlog)
  dev.off()
}

o <- ov.list$IO %>%
  filter(!grepl("guillimin|ubuntu", method))

no.labels <- 
  ggplot()+
  geom_point(aes(query.rows, seconds, color=method),
             data=o, pch=1)+
    scale_color_manual(values=algo.colors)+
  scale_x_continuous("rows in bedGraph file", breaks=c(1, 2)*1e7,
                     limits=c(NA, 4e7))+
  guides(color="none")+
  theme_grey()+
  scale_y_continuous("seconds")

dl <-
  direct.label(no.labels, "over.end")

pdf("figure-TF-benchmark-IO-some.pdf", 5, 3)
print(dl)
dev.off()

refs <-
  data.table(unit=c("1 second", "1 minute"),
             seconds=c(1, 60),
             vjust=c(-0.5, 1.5))
read.times <- TF.benchmark$read %>%
  mutate(method=ifelse(expr=="read.table", "read.table", "data.table::fread"))
lab.df <-
  data.table(rows=3.6e7,
             method=c("data.table::fread", "read.table"),
             seconds=c(10, 50))
with.labels <-
  ggplot()+
  ## geom_text(aes(rows, seconds, label=method, color=method),
  ##           data=lab.df, size=3)+
  ## geom_hline(aes(yintercept=seconds), data=refs, color="grey50")+
  ## geom_text(aes(3.75e7, seconds, label=unit, vjust=vjust),
  ##           data=refs, color="grey50", size=3)+
    scale_color_manual(values=algo.colors)+
  geom_point(aes(rows, time/1e9, color=method),
             data=read.times, pch=1)+
  scale_x_continuous("rows in bedGraph file", limits=c(NA, 3.5e7),
                     breaks=c(1, 2)*1e7)+
  ylab("seconds")+
  guides(color="none")
dl <- direct.label(with.labels, "over.end")

pdf("figure-TF-benchmark.pdf", 5, 3)
print(dl)
dev.off()
