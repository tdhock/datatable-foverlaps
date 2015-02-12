works_with_R("3.1.2",
             reshape2="1.2.2",
             dplyr="0.4.0",
             xtable="1.7.3")

load("simple.RData")
simple %>%
  group_by(method) %>%
  summarise(correct=sum(status=="correct"),
            tests=n())
short <- simple %>%
  mutate(status=ifelse(status=="correct", "ok", "incorrect"),
         test.name=paste0("test", test.i))
wide <- dcast(short, method ~ test.name, value.var="computed")
test.info <- short %>%
  filter(method=="foverlaps")
test.dt <- test.info[, .(chromStart, chromEnd, expected)]
test.mat <- t(as.matrix(test.dt))
test.methods <- cbind(rownames(test.mat), test.mat)

xt <- xtable(rbind(test.methods, as.matrix(wide)))
print(xt, include.rownames=FALSE, include.colnames=FALSE,
      hline.after=c(0, 3, nrow(xt)),
      file="table-simple.tex")
