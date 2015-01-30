works_with_R("3.1.2",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7")

load("piecewise.constant.RData")

ggplot()+
  geom_point(aes(ref.rows, time/1e9, color=expr),
             data=piecewise.constant, pch=1)+
  facet_wrap("query.rows")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))

ggplot()+
  geom_point(aes(ref.rows, time/1e9, color=expr),
             data=piecewise.constant, pch=1)
