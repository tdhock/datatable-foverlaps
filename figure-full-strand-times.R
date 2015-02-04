works_with_R("3.1.2",
             ##"Rdatatable/data.table@84ba1151299ba49e833e68a2436630216b653306",
             dplyr="0.4.0",
             reshape2="1.2.2",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7")

load("full.strand.times.RData")

full.strand.times %>%
  mutate(seconds=time/1e9) %>%
  select(expr, strand, chipseq.rows, seconds) %>%
  arrange(chipseq.rows, expr)

ggplot()+
  geom_point(aes(time/1e9, expr), data=full.strand.times)+
  facet_grid(. ~ strand)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))

ggplot()+
  scale_x_log10()+
  scale_y_log10()+
  geom_point(aes(chipseq.rows, time/1e9, color=expr),
             data=full.strand.times, pch=1)

levs <-
  c("intersectBed",
    "fread.foverlaps.write",
    "data.table::foverlaps",
    "GenomicRanges::findOverlaps",
    "data.table::foverlaps\nwindows.read",
    "GenomicRanges::findOverlaps\nwindows.read")
times <- full.strand.times %>%
  arrange(strand, expr) %>%
  mutate(expr.fac=reorder(expr, time),
         expr.chr=sub(".windows.read", "\nwindows.read", expr),
         expr.fac2=factor(expr.chr, levs),
         seconds=time/1e9,
         size=sprintf("%.3f", chipseq.rows/1e6))
mean.times <- times %>%
  group_by(strand, expr) %>%
  summarise(mean=mean(seconds))
dcast(mean.times, expr ~ strand)

dots <-
  ggplot()+
  scale_x_log10("seconds")+
  geom_point(aes(seconds, expr.fac2), data=times, pch=1)+
  facet_grid(. ~ size, labeller=label_both)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))

pdf("figure-full-strand-times.pdf", 5, 3)
print(dots)
dev.off()
