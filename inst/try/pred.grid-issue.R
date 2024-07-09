#' Experiment with pred.grid to work with predict.lm
#'
#' Problem: if the original data has factors with non-default contrasts
#' predict.lm (and perhaps others) don't work with newdata even if
#' the new data frame has correctly specified contrasts.
#'
#' We explore the following questions:
#'
#' - if a factor has treatment contrasts but the base level is not the
#'   first lexicographically,
#'   - will things work with expand.grid generated with levels of the factors
#'
{
  set.seed(23453)   # rank deficient model
  set.seed(23451)   # full rank model
  dd <- data.frame(
  sex = factor(c('m','f','a')[sample(3,20,T)], levels = c('m' ,'f','a')),
  rank = factor(c('asp','acp','p')[sample(3,20,T)], levels = c('asp' ,'acp','p')),
  age = 51:(51+19),
  y = 100 + rnorm(20)

)
fit <- glm(y ~ age + sex*rank, dd, family = 'gaussian')   #

}
summary(fit)
dd$fit <- predict(fit)


fun1 <- function (...)   # spida2::pred.grid on 2024-07-08
{
  nams <- as.character(as.list(substitute(list(...)))[-1L])
  x <- list(...)
  if (is.null(names(x)))
    names(x) <- nams
  else if (any(names(x) == ""))
    names(x)[names(x) == ""] <- nams[names(x) == ""]
  ret <- lapply(x, unique)
  ret <- lapply(ret, sort)
  ret <- do.call(expand.grid, c(ret, stringsAsFactors = FALSE))
  for (nn in names(ret)) {
    if (is.factor(ret[[nn]]))
      contrasts(ret[[nn]]) <- contrasts(x[[nn]])
  }
  ret
}

pred1 <- with(dd, fun1(sex, rank, age))
pred1 <- with(dd, fun1(sex, rank, age))
pred1$fit <- predict(fit, newdata = pred1)   # contrasts dropped



############
fun2 <- function (...)  { # spida2::pred.grid on 2024-07-08
  nams <- as.character(as.list(substitute(list(...)))[-1L])
  x <- list(...)
  if (is.null(names(x)))
    names(x) <- nams
  else if (any(names(x) == ""))
    names(x)[names(x) == ""] <- nams[names(x) == ""]
  uniq <- function(x) {
    # seems to preserve contrasts, etc. for factors
    # but keeps only unique values of numerical and other inputs
    x <- sort(x)
    ind <- c(1, diff(as.numeric(factor(x))))
    x[ind == 1]
  }
  x <- lapply(x, uniq)
  ret <- do.call(expand.grid, c(x, stringsAsFactors = FALSE)
  ret
}





colnames(helm) <- c('first','second')
contrasts(dd$sex) <- helm
pred <- with(dd, fun2(sex, age = 51, rank) )
fit <- lm(y ~ age + sex * rank, dd)
pred$fit <- predict(fit, newdata = pred)
# different approach
#
pred1 <- subset(dd, select = c(sex,rank))
pred1 <- pred1[!duplicated(pred1),]
pred1$age <- 51
predict(fit, newdata = dd)
predict.lm


dd$fitd <- predict(fit, newdata = dd)

dd$fit <- predict(fit)

head(dd)
contrasts(dd$sex)
contrasts(pred$sex)

  for (nn in names(ret)) {
    if (is.factor(ret[[nn]]))
      contrasts(ret[[nn]]) <- contrasts(x[[nn]])
  }
  ret
  x
}
dd$sex <- relevel(dd$sex, 'm')

with(dd, fun2(sex, age, rank)) $sex |> contrasts()    # correct
with(dd, fun2(sex, age, rank)) $sex |> unique() |> contrasts()    # XXX in              correct

contrasts(dd$sex) <- cbind(c(-1,1,0)/2,c(1,1,-2)/2 )

with(dd, fun2(sex, age = age, rank))
with(dd, fun2(sex, age = age, rank))$sex |> contrasts()    # correct


pred1$fit <- predict(fit, newdata = pred1)   # contrasts dropped

}

