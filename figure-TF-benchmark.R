works_with_R("3.1.2",
             dplyr="0.4.0",
             "hadley/tidyr@3de52c46f12d0d4cba603aa4e63e39f2370a9cfd",
             "Rdatatable/data.table@84ba1151299ba49e833e68a2436630216b653306",
             directlabels="2014.6.13",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7")

load("TF.benchmark.RData")

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

ov <- TF.benchmark$overlap %>%
  mutate(result=rep(c("indices", "tables"), each=4),
         method=sub(".tables|indices$", "", expr)) %>%
  filter(result == "tables")
stopifnot(nrow(ov) == nrow(TF.benchmark$overlap)/2)

overlab.df <-
  data.table(rows=3.5e7,
             method=c("GenomicRanges::findOverlaps", "data.table::foverlaps"),
             seconds=c(6.5, 15))
with.labels <- 
  ggplot()+
  geom_text(aes(rows, seconds, label=method, color=method),
            data=overlab.df, size=3)+
  geom_hline(aes(yintercept=seconds), data=refs, color="grey50")+
  geom_text(aes(3.75e7, seconds, label=unit, vjust=vjust),
            data=refs, color="grey50", size=3)+
  geom_point(aes(query.rows, time/1e9, color=method),
             data=ov, pch=1)+
  xlab("rows in bedGraph file")+
  guides(color="none")

with.labels+
  geom_point(aes(query.rows, minutes*60, color=method),
             data=minutes.tall)+
  scale_y_log10()+
  scale_x_log10()

pdf("figure-TF-benchmark-overlap.pdf", 5, 3)
print(with.labels)
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
  geom_text(aes(rows, seconds, label=method, color=method),
            data=lab.df, size=3)+
  geom_hline(aes(yintercept=seconds), data=refs, color="grey50")+
  geom_text(aes(3.75e7, seconds, label=unit, vjust=vjust),
            data=refs, color="grey50", size=3)+
  geom_point(aes(rows, time/1e9, color=method),
             data=read.times, pch=1)+
  xlab("rows in bedGraph file")+
  guides(color="none")

pdf("figure-TF-benchmark.pdf", 5, 3)
print(with.labels)
dev.off()
