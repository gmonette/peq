#' Graphical comparison of AIC and BIC
#'
#' Compares AIC and BIC for a named list of models. The models should have been fitted
#' to the same data.
#'
#' Updated with version in UTFA
#'
#' @param ll a list of models using the same data.
#' @param rot,srt angle to rotate labels
#' @param prefix for labels
#'
#' @examples
#' fit1 <- lm(mpg ~ wt, mtcars)
#' fit2 <- lm(mpg ~ wt + cyl, mtcars)
#' fit3 <- lm(mpg ~ wt * hp + cyl , mtcars)
#' fitl <- list('wt' = fit1, 'wt+cyl' = fit2, fit3=fit3)
#' aic(fitl)
#' ics(fitl)
#' icp(fitl)
#' @export
icp <-
  function(ll, ..., rot = 0, srt = rot, prefix = ""){
    if(is.null(names(ll))) {
      # names(ll) <- paste0(prefix,seq_along(ll))
      names(ll) <- sapply(ll, function(x) as.character(formula(x)[-c(1,2)]))
    }
    a <- ics(ll)
    a$model <- sapply(ll, formula)
    a <- sortdf(a, ~ edf)
    tobj <- xyplot(AIC +BIC~ edf, a, type = 'l', outer = T,...,
                   labs = rownames(a), fonts = 2,
                   subscripts = TRUE,
                   layout = c(1,2),
                   srt = srt,
                   scales = list(y = list(relation = 'free')))+
      layer(panel.text(..., labels = labs, fonts = 2, srt = srt))
    print(tobj)
    print(a)
    print(sortdf(a, ~ AIC))
    invisible(ll)
  }
#' AIC and BIC for a list of functions
#'
#' @param ll a list of models using the same data.
#' @describeIn icp compute AIC for a list of models.
#' @export
aic <- function(ll, k =2){
  if(is.null(names(ll))) names(ll) <- as.character(seq_along(ll))
  # ret <- as.data.frame(do.call(rbind,lapply(ll, extractAIC, k=2)))
  ret <- as.data.frame(do.call(rbind,Map(extractAIC,ll, k=k)))
  names(ret) <- c('edf','AIC')
  ret

  #cmd <- paste('with(ll, AIC(',paste(names(ll), collapse = ','),'))',collapse = '')
  #eval(str2lang(cmd))
}
#' @describeIn icp compute BIC for a list of models.
#' @export
bic <- function(ll, k = sapply(ll,function(f)log(nobs(f)))){
  ret <- aic(ll,k)
  names(ret)[2] <- "BIC"
  ret
}
#' Combine and scale AIC and BIC
#'
#' @describeIn icp compute BIC for a list of models
#' @export
ics <- function(ll){
  ret <- aic(ll)
  ret$BIC <- bic(ll)$BIC
  min_bic <- min(ret$BIC)
  min_aic <- min(ret$AIC)
  ret <- within(ret,
                {
                  BIC <- BIC - min_bic
                  AIC <- AIC - min_aic
                })
  ret$fmla <- sapply(ll, formula)
  ret
}
