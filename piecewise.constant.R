works_with_R("3.1.2",
             microbenchmark="1.3.0",
             "Rdatatable/data.table@84ba1151299ba49e833e68a2436630216b653306")

load("oracle.regularized.RData")
load("oracle.intervals.RData")
load("oracle.optimal.RData")
load("dp.peaks.sets.RData")

result.list <- list()
for(set.name in names(dp.peaks.sets)){
  train.sets <- dp.peaks.sets[[set.name]]
  interval.list <- oracle.intervals[[set.name]]
  optimal.list <- oracle.optimal[[set.name]]
  fit.list <- oracle.regularized$models[[set.name]]
  for(set.i in seq_along(train.sets)){
    train.chunks <- train.sets[[set.i]]
    test.chunks <- names(optimal.list)[!names(optimal.list) %in% train.chunks]
    fit <- fit.list[[set.i]]$L1.reg
    testSet <- names(fit.list)[[set.i]]
    for(test.chunk in test.chunks){
      test.features <- interval.list[[test.chunk]]$all.features
      pred.log.lambda <- fit$predict(test.features)
      optimal.samples <- optimal.list[[test.chunk]]
      ## How fast can we turn these log.lambda values into the optimal
      ## number of segments for each sample?
      peaks <- data.frame(pred.log.lambda, foverlaps=NA, for.list=NA)
      times <- microbenchmark(foverlaps.setkey2={
        optimal.dt.list <- list()
        for(sample.id in names(optimal.samples)){
          optimal.dt.list[[sample.id]] <-
            data.table(sample.id, optimal.samples[[sample.id]])
        }
        optimal.dt <- do.call(rbind, optimal.dt.list)
        setkey(optimal.dt, sample.id, min.log.lambda, max.log.lambda)
        pred.dt <-
          data.table(sample.id=rownames(pred.log.lambda),
                     log.lambda=as.numeric(pred.log.lambda),
                     log.lambda2=as.numeric(pred.log.lambda))
        setkey(pred.dt, sample.id, log.lambda, log.lambda2)
        overlap.dt <- foverlaps(pred.dt, optimal.dt)
      }, for.list={
        for(sample.id in names(optimal.samples)){
          log.lambda <- pred.log.lambda[sample.id, ]
          peaks[sample.id, "for.list"] <- 
            subset(optimal.samples[[sample.id]],
                   min.log.lambda < log.lambda &
                   log.lambda < max.log.lambda)$peaks
        }
      }, times=10)
      ## check results.
      stopifnot(nrow(overlap.dt) == nrow(pred.dt))
      peaks[pred.dt$sample.id, "foverlaps"] <- overlap.dt$peaks
      with(peaks, stopifnot(foverlaps == for.list))
      ## What if we start with the dt?
      times2 <- microbenchmark(foverlaps.setkey={
        pred.dt <-
          data.table(sample.id=rownames(pred.log.lambda),
                     log.lambda=as.numeric(pred.log.lambda),
                     log.lambda2=as.numeric(pred.log.lambda))
        setkey(pred.dt, sample.id, log.lambda, log.lambda2)
        overlap.dt <- foverlaps(pred.dt, optimal.dt)
      }, foverlaps={
        overlap.dt <- foverlaps(pred.dt, optimal.dt)
      }, times=10)
      result.list[[paste(set.name, set.i, test.chunk)]] <-
        data.table(set.name, set.i, test.chunk,
                   ref.rows=nrow(optimal.dt),
                   query.rows=nrow(pred.dt),
                   overlap.rows=nrow(overlap.dt),
                   rbind(times, times2))
    }
  }
}
## foverlaps is not faster than a simple vector scan in these cases
## where the reference data.frames are not very big.
piecewise.constant <- do.call(rbind, result.list)

save(piecewise.constant, file="piecewise.constant.RData")
