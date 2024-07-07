#' Simulation of adjustment to -.5 SE
#' Over time 75% are are below the new mean with a huge number at -0.5 SD
#' But the effect on the overall mean is positive, with the minimum
#' salary being asymptotically 0.5 SDs higher than where the mean would have
#' been with no anomaly adjustements. 

#' Effect of adjustment to -.5 SE
#' with shrink
#'
target <- -.5
shrink <- .5
#'
zd <- list()
set.seed(123)
zd[[1]] <- list()
zd[[1]]$y <- rnorm(1000)
zd[[1]]$mean <- mean(zd[[1]]$y)
zd[[1]]$sd <- sd(zd[[1]]$y)
for(step in 2:20) {
  y <- zd[[step-1]]$y
  m <- zd[[step-1]]$mean
  s <- zd[[step-1]]$sd
  y <- y + (y < m+ target*s) * shrink * (-y + (m+target*s))
  zd[[step]] <- list()
  zd[[step]]$y <- y
  zd[[step]]$mean <- mean(y)
  zd[[step]]$sd <- sd(y)
}
zdd <- data.frame(step = 1:10, mean = sapply(zd, function(l) l$mean),
                  sd = sapply(zd, function(l) l$sd))
zdd

zall <- data.frame(y = unlist(sapply(zd, function(l) l$y)))
dim(zall)
head(zall)
zall <- tolong(zall, sep = '.', timevar = 'step')
names(zall)
zall$median <- with(zall, capply(y, step, median))
zall$q1 <- with(zall, capply(y, step, quantile, .25))
zall$q3 <- with(zall, capply(y, step, quantile, .75))
zall$mean <- with(zall, capply(y, step, mean))
zall$sd <- with(zall, capply(y, step, sd))
zall$medsd <- with(zall, (median-mean)/sd)
zall <- sortdf(zall, ~ step)
xyplot(y ~ step, zall) +
  layer_(panel.bwplot(..., horizontal = FALSE))+
  xyplot(mean ~ step, zall, type = 'l')+
  xyplot(q1 ~ step, zall, type = 'l', col = 'blue')+
  xyplot(q3 ~ step, zall, type = 'l', col = 'blue')+
  xyplot(I(mean+sd) ~ step, zall, type = 'l', col = 'red')+
  xyplot(I(mean-sd) ~ step, zall, type = 'l', col = 'red')
zall %>% lapply(class)
tab(zall, I(y<mean) ~ step)  
tab(zall, I(y<mean) ~ step, pct = 0)  
xyplot(medsd ~ step, zall)
