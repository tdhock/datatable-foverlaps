works_with_R("3.1.2",
             dplyr="0.4.0",
             reshape2="1.2.2",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7")

load("full.strand.times.RData")

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

times <- full.strand.times %>%
  arrange(chipseq.rows, expr) %>%
  mutate(expr.fac=reorder(expr, time),
         seconds=time/1e9)
mean.times <- times %>%
  group_by(strand, expr) %>%
  summarise(mean=mean(seconds))
dcast(mean.times, expr ~ strand)

dots <- ggplot()+
  scale_x_log10("seconds")+
  geom_point(aes(seconds, expr.fac), data=times, pch=1)+
  facet_grid(. ~ strand)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))

pdf("figure-full-strand-times.pdf", 5, 3)
print(dots)
dev.off()
