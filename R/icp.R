#' AIC for a list of functions
#'
#' @param ll a list of models using the same data.
#' @export
aic <- function(ll){
  cmd <- paste('with(ll, AIC(',paste(names(ll), collapse = ','),'))',collapse = '')
  eval(str2lang(cmd))
}
#' @export
bic <- function(ll){
  cmd <- paste('with(ll, BIC(',paste(names(ll), collapse = ','),'))',collapse = '')
  eval(str2lang(cmd))
}
#' @export
icp <- function(ll, ...){
  a <- aic(ll)
  b <- bic(ll)
  a$BIC <- b$BIC

  a <- sortdf(a, ~ df)
  pl <- xyplot(AIC +BIC~ df, a, type = 'l',
               outer = T,
               ...,
               labs = rownames(a), fonts = 2,
               subscripts = TRUE,
               layout = c(1,2),
               scales = list(y = list(relation = 'free')))+
    layer(panel.text(..., labels = labs, fonts = 2))
  print(pl)
  print(a)
  print(sortdf(a, ~ AIC))
  invisible(ll)
}
