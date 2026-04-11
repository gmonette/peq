#' Gender gap decompositions without SEs
#'
#' @param fitl a named list of models
#' @param g character string with the name of the variable used for grouping, e.g. "Gender"
#' @param comp character string giving the name of the level used as a comparator, e.g. "M"
#' @param comon data set from which all models can be fit
#' @param cond optional vector with names of variables used for conditioning
#' @param refit logical (default TRUE) should models be refit using 'data'
#'
#' @returns an object of class 'decomps' consisting of a list with two data frames: the first, called 'gaps', consists of the
#' rbind of 'data' for each model
#' together with additional variables:
#' \itemize{
#'   \item res_: ordinary residual
#'   \item pred_: ordinary predicted value
#'   \item resp_: response
#'   \item predg_: value predicted if comparator
#'   \item resg_: residual from value predicted if comparator
#'   \item model: identifying the model using 'names(fitl)'
#' }
#'
#' The second component, called 'summ', is the mean summary of
#' of 'gaps' keeping variables indicating each model, the
#' grouping ('g') and conditioning ('cond') variables, if any.
#' @export
decomps <- function(fitl, g = "Gender", comp = "M", data = na.omit(getD(full)), cond = NULL, refit = TRUE) {
  disp <- function(...) {
    NULL
  }
  # 2024-08-18: simplified version with just residuals
  along <- function(x) {
    ret <- seq_along(x)
    names(ret) <- names(x)
    ret
  }
  addrownames <- function(x) {
    rownames(x) <- as.character(seq_len(nrow(x)))
    x
  }

  incmat <- function(d1, d2) {
    # incidence matrix for position in d2 for combinations of values in d1
    # d1 and d2 are data frames with variables in d1 a subset of variables in d2
    imats <- list()
    for(v in names(d1)) {
      imats[[v]] <- 1* outer(d1[[v]], d2[[v]], '==')
    }
    Reduce('*',imats)
  }
  pred.grid_ <- function(x) {
    x <- lapply(x, unique)
    # to ensure that factors are in their internal order, not the order of their appearance in the data
    x <- lapply(x, sort)
    expand.grid(x, stringsAsFactors = FALSE)
  }
  #
  library(car)
  library(spida2)

  ret <- list()
  ret[['names']] <- list(gname = g, gcomplevel = comp)
  ret[['names']][['cond']] <- cond

  #
  # Refit with data set used for full model to check for common data set
  #

  full <- fitl[[length(fitl)]]
  ret[['data']] <- data
  d <- data
  if(refit) {
    fit_up <- lapply(fitl, update, data = d)
  }
  else {
    fit_up <- fitl
  }

  fitl <- fit_up

  #
  # Create id variable to link between models in long form for dout
  #

  idname <- 'id.'
  while(idname %in% names(d)) { # generate idname not in dataset
    idname <- paste0(idname,'.')
  }
  ret[['names']][['idname']] <- idname
  d[[idname]]  <- 1:nrow(d)

  #
  # Create comparator data set
  #

  dg <- d
  dg[[g]][] <- comp   # set to level of comparator group

  #
  # Generate data sets with mean residuals names coef
  # using different strategies
  #

  res <- lapply(fitl, resid)
  pred <- lapply(fitl, predict)
  resp <- lapply(fitl, function(f) resid(f)  + predict(f))
  predg <- lapply(fitl, predict, newdata = dg)

  gaps  <- d[rep(1:nrow(d), length(fitl)), ]
  gaps$model <- factor(rep(names(fitl), each = nrow(d)), levels = names(fitl))

  gaps$res_ <- unlist(res)
  gaps$pred_ <- unlist(pred)
  gaps$resp_ <- unlist(resp)
  gaps$predg_ <- unlist(predg)
  gaps$resg_ <- with(gaps, resp_ - predg_)
  gaps$diffg_ <- with(gaps, pred_ - predg_)


  ret[['gaps']] <- gaps

  fmup <- paste('~ model +', g)
  if(!is.null(cond)) fmup <- paste(c(fmup, cond), collapse = '+')

  summ <- up(gaps, as.formula(fmup), agg = ~ res_ + pred_ + resp_ + predg_ + resg_ + diffg_)
  summ <- sortdf(summ, ~ model)

  ret[['summ']] <- summ
  class(ret) <- 'decomps'
  ret
}

#' @describeIn decomp experimental conditional version with 'cond' argument
#' @export
decomp3 <- function(fitl, g, comp, data = na.omit(getD(full)), cond = NULL, refit = TRUE) {
  disp <- function(...) {
    NULL
  }
  # 2024-08-18:
  #    add colums to dout for different resids
  #    create output components like 'gaps_each'
  #    for various residuals
  # exportPattern("^[[:alpha:]]+")
  # fitl: list of fitted models that should be nested with fullest model last
  # g; character: names of group variable
  # comp: character: comparator level name
  #
  # Might consider including a full model beyond but last
  # but not clear what happens if some contrasts are estimable in
  # smaller models but not in fullest model
  #
  # What does this do?
  #
  # 1. Takes the regression with groups
  # 2. Sets group to comparator for all
  # 3. Computes the difference between the average predicted value
  #    (which is equal to the mean value if comparator with
  #    the predicted value from model
  #
  # Last fit is full fit
  # g is name of group variable and comp is comparator level
  #
  # - produce data frame with individual discrepancies
  #   should be long data frame replicating original data for each model
  # - and 2 wald object with all gaps for gap plot
  # - and for between gap p-values (proportion explained by ....)
  # - one with discrepancies and one
  # - changes in discrepancies
  # - return with a class and write a plot methods to produce
  #   a variety of plots
  #   and a print method
  #
  # substitute for lsfit that produces numerical errors with singular models
  #
  # lssvd <- function(x, y, zero = 10^(-14)) {
  #   xp <- svd(x,nu=ncol(x), nv = ncol(x))
  #   uy <- t(xp$u)%*%cbind(y)
  #   dinv <- 1/xp$d
  #   dinv[abs(xp$d) < zero] <- 0
  #   coef <- t(t(xp$v)*dinv) %*% uy
  #   resid <- y - x %*% coef
  #   list(coef = coef, resid = resid, sse = sum(resid^2))
  # }
  #
  along <- function(x) {
    ret <- seq_along(x)
    names(ret) <- names(x)
    ret
  }
  addrownames <- function(x) {
    rownames(x) <- as.character(seq_len(nrow(x)))
    x
  }

  incmat <- function(d1, d2) {
    # incidence matrix for position in d2 for combinations of values in d1
    # d1 and d2 are data frames with variables in d1 a subset of variables in d2
    imats <- list()
    for(v in names(d1)) {
      imats[[v]] <- 1* outer(d1[[v]], d2[[v]], '==')
    }
    Reduce('*',imats)
  }
  pred.grid_ <- function(x) {
    x <- lapply(x, unique)
    # to ensure that factors are in their internal order, not the order of their appearance in the data
    x <- lapply(x, sort)
    expand.grid(x, stringsAsFactors = FALSE)
  }
  #
  library(car)
  library(spida2)

  ret <- list()
  ret[['names']] <- list(gname = g, gcomplevel = comp)
  ret[['names']][['cond']] <- cond

  #
  # Refit with data set used for full model to check for common data set
  #

  full <- fitl[[length(fitl)]]
  ret[['data']] <- data
  d <- data
  if(refit) {
    fit_up <- lapply(fitl, update, data = d)
  }
  else {
    fit_up <- fitl
  }
  compcoefs <- lapply(seq_along(fitl),
                      function(ii) {
                        cbind(
                          compareCoefs(fitl[[ii]], fit_up[[ii]], print = FALSE),
                          diff = coef(fitl[[ii]])-coef(fit_up[[ii]]))
                      })
  names(compcoefs) <- names(fitl)
  ret[['compare']] <- compcoefs

  fitl <- fit_up

  #
  # Create id variable to link between models in long form for dout
  #

  idname <- 'id.'
  while(idname %in% names(d)) { # generate idname not in dataset
    idname <- paste0(idname,'.')
  }
  ret[['names']][['idname']] <- idname
  d[[idname]]  <- 1:nrow(d)

  #
  # Create comparator data set
  #

  dg <- d
  dg[[g]][] <- comp   # set to level of comparator group

  #
  # Geerate data sets with mean residuals names 'coef'
  # using different strategies
  #

  p_res <- lapply(fitl, resid)
  p_pred <- lapply(fitl, predict)
  p_y <- lapply(fitl, function(fit) predict(fit) + resid(fit))
  p_pred_comp <- lapply(fitl, predict, newdata = dg)

  d_gaps  <- d[rep(1:nrow(d), length(fitl)), ]
  d_gaps$model <- factor(rep(names(fitl), each = nrow(d)), levels = names(fitl))

  d_gaps$p_res <- unlist(p_res)
  d_gaps$p_pred <- unlist(p_pred)
  d_gaps$p_y <- unlist(p_y)
  d_gaps$p_pred_comp <- unlist(p_pred_comp)
  d_gaps$p_res_comp <- with(d_gaps, p_y - p_pred_comp)

  d_gaps$coef <- d_gaps$p_res_comp

  fmup <- paste('~ model +', g)
  if(!is.null(cond)) fmup <- paste(c(fmup, cond), collapse = '+')

  gaps_comp <- up(d_gaps, as.formula(fmup), agg = ~ coef)
  gaps_comp <- sortdf(gaps_comp, ~ model)



  ret[['gaps_comp']] <- gaps_comp %>% addrownames()

  #
  # actual resids from comparator
  #
  d1 <- d[,c(g,cond), drop=FALSE]
  d1 <- pred.grid_(d1)
  groupL <-incmat(d1, d)
  groupL <- groupL/apply(groupL,1,sum)
  dc_groupL <- d1
  dc_groupL$groupL <- groupL
  ret[['groupL']] <- dc_groupL

  ## Individual residuals from each fit

  gresids <- lapply(fitl, function(f) {
    y <- predict(f) + resid(f)
    y - predict(f, newdata = dg)
  })
  ggaps <- lapply(fitl, function(f) {
    y <- predict(f)
    y - predict(f, newdata = dg)
  })
  inds <- rep(1:nrow(d), length(gresids))
  dout <- d[inds,]
  dout$gresids <- unlist(gresids)
  dout$ggaps <- unlist(ggaps)
  dout$model <- factor(rep(names(fitl), each = nrow(d)), levels = names(fitl))
  ret[['dout']] <- dout

  #
  # gaps calculated from full model
  #

  mfs <- lapply(fitl, model.matrix, data = d)
  mgs <- lapply(fitl, model.matrix, data = dg)
  mffull <- mfs[[length(mfs)]]
  Bs_qr <- lapply(seq_along(mfs), function(ii) lsfit(mfs[[ii]], mffull, intercept =FALSE )$coef)
  Bs <- lapply(seq_along(mfs), function(ii) lssvd(mfs[[ii]], mffull)$coef)
  resids <- lapply(
    along(mfs),
    function(ii) {
      sum(abs(lsfit(mfs[[ii]], mffull, intercept =FALSE )$residuals))
    }
  )
  resids2 <- lapply(
    along(mfs),
    function(ii) {
      sum(abs(lsfit(mffull, mfs[[ii]], intercept =FALSE )$residuals))
    }
  )
  names(resids) <- names(fitl)
  Ls <- lapply(seq_along(mfs), function(ii) {
    groupL %*% ((mfs[[ii]] - mgs[[ii]]) %*% Bs[[ii]])
  })
  Lpred <- lapply(seq_along(mfs), function(ii) {
    groupL %*% mfs[[ii]]%*% Bs[[ii]]
  })
  Ls <- do.call(rbind, Ls)
  Lpred <- do.call(rbind, Lpred)
  # data <- expand.grid(zork = levels(d[[g]]), model = names(fitl))
  data <- cbind(
    d1[rep(1:nrow(d1), length(fitl)),],
    model = rep(factor(names(fitl), levels = names(fitl)), each = nrow(d1))
  )
  # names(data)[1] <- g
  attr(Ls, 'data') <- data
  gaps <- waldf(full,Ls)

  attr(Lpred, 'data') <- data
  pred <- waldf(full, Lpred)

  ret[['diags']][['Bs_qr']] <- Bs_qr
  ret[['diags']][['Bs']] <- Bs
  ret[['diags']][['B_diffs']] <- lapply(along(Bs),
                                        function(ii) {
                                          Bs[[ii]] - Bs_qr[[ii]]
                                        })

  ret[['gaps']] <- gaps %>% addrownames
  ret[['resids']] <- resids
  ret[['resids2']] <- resids2
  # disp(resids2)
  # disp(unlist(resids2) > 10^(-11))
  if(any(unlist(resids2) > 10^(-11))) warning('Some models may not be nested in last model. See resids2')

  ret[['pred']] <- pred %>% addrownames

  #
  # gaps calculated from individual models
  #

  # from above: mfs <- lapply(fitl, model.matrix, data = d)
  # from above: mgs <- lapply(fitl, model.matrix, data = dg)
  # not needed: mffull <- mfs[[length(mfs)]]
  # not needed: Bs <- lapply(seq_along(mfs), function(ii) lsfit(mfs[[ii]], mffull, intercept =FALSE )$coef)
  Ls_each <- lapply(
    along(mfs),
    function(ii) {
      Lmat <- groupL %*% (mfs[[ii]] - mgs[[ii]])    # not needed:  %*% Bs[[ii]]
      attr(Lmat, 'data') <- subset(data, model == names(fitl)[ii])
      waldf(fitl[[ii]], Lmat)
    })
  ret[['gaps_each']] <-Ls_each %>% lapply(subset, select = -L) %>% do.call(rbind,.) %>% addrownames

  ## Comparing fitted values from 'predict' and from X %*% beta ------ is okay

  fitted_values <- lapply(along(fitl), function(ii) {
    list() %>%
      within({
        yhat <- predict(fitl[[ii]])
        yhatmat <- mfs[[ii]] %*% coef(fitl[[ii]])
        yhatg <- predict(fitl[[ii]], newdata = dg)
        yhatmatg <- mgs[[ii]] %*% coef(fitl[[ii]])
        Lyhat <- groupL %*% (yhat - yhatg)
        Lymat <- groupL %*% (yhatmat - yhatmatg)
        wald <- waldf(fitl[[ii]], groupL %*% (mfs[[ii]] - mgs[[ii]]))

      })
  })
  ret[['fitted_values']] <- fitted_values

  #
  # Differences between models
  #

  diffmat <- function(n) {
    cbind(-diag(n-1),0) + cbind(0, diag(n-1))
  }
  # Ldiffs <- kronecker(diffmat(length(fitl)), diag(length(levels(d[[g]])))) %*% Ls
  Ldiffs <- kronecker(diffmat(length(fitl)), diag(nrow(d1))) %*% Ls
  ret[['Ldiffs']] <- Ldiffs
  gapdiffs <- paste(names(fitl)[-1], '-', names(fitl)[-length(names(fitl))])
  #data <- expand.grid(zork = levels(d[[g]]), gapdiffs = paste(names(fitl)[-1], '-', names(fitl)[-length(names(fitl))]))
  data <- cbind(
    d1[rep(1:nrow(d1),length(gapdiffs)),],
    gapdiffs = rep(gapdiffs, each = nrow(d1))
  )

  # names(data)[1] <- g
  attr(Ldiffs, 'data') <- data
  gapdiffs <- waldf(full, Ldiffs)
  ret[['gapdiffs']] <- gapdiffs

  #
  # Combine residuals and 'penalties'
  #

  geach <- ret[['gaps_each']]
  geach$Type <- 'Each'
  gfull <- ret[['gaps']]
  gfull$Type <- 'Full'
  gpen  <- ret[['gapdiffs']]
  gpen$Type <- 'Penalty'
  gpen$model <- gpen$gapdiffs
  allgaps <- Rbind(geach,gfull,gpen)

  ret[['allgaps']] <- allgaps

  class(ret) <- 'decomp'
  ret
}
