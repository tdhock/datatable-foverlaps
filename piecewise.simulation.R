works_with_R("3.1.2",
             microbenchmark="1.3.0",
             ##"Rdatatable/data.table@84ba1151299ba49e833e68a2436630216b653306")

set.seed(1)
for(n.breaks in c(1, 10, 100)){
  breaks <- sort(runif(n.breaks))
  fun.dt <-
    data.table(min=c(0, breaks),
               max=c(breaks, 1),
               value=c(0, seq_along(breaks)))
  ggplot()+
    geom_segment(aes(min, value, xend=max, yend=value), data=fun.dt)
  point.vals <- c(1, 10, 100)
  x.all <- runif(max(point.vals))
  for(n.points in point.vals){
    x <- x.all[1:n.points]
    a <- approx(c(0, fun.dt$max), c(0,fun.dt$value,
  }
}
