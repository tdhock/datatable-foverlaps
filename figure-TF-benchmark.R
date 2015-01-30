works_with_R("3.1.2",
             directlabels="2014.6.13",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7")


load("TF.benchmark.RData")

## from https://github.com/Rdatatable/data.table/wiki/talks/EARL2014_OverlapRangeJoin_Arun.pdf
minutes.wide <-
  data.frame(query.rows=c(80e6, 65e6),
             subject.rows=c(33e3, 36e3),
             foverlaps=c(2, 4),
             findOverlaps=c(16, 28* 24))

refs <-
  data.table(unit=c("1 second", "1 minute"),
             seconds=c(1, 60),
             vjust=c(-0.5, 1.5))
lab.df <-
  data.table(rows=3.5e7,
             expr=c("fread", "read.table"),
             seconds=c(10, 50))

with.labels <-
  ggplot()+
  geom_text(aes(rows, seconds, label=expr, color=expr), data=lab.df)+
  geom_hline(aes(yintercept=seconds), data=refs, color="grey50")+
  geom_text(aes(3.75e7, seconds, label=unit, vjust=vjust),
            data=refs, color="grey50")+
  geom_point(aes(rows, time/1e9, color=expr),
             data=TF.benchmark$read, pch=1)+
  guides(color="none")

pdf("figure-TF-benchmark.pdf")
print(with.labels)
dev.off()
