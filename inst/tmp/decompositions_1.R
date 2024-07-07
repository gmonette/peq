
decomp <- function(fitl, g, comp, data = na.omit(getD(full)), refit = TRUE) {
  disp <- function(...) {
    NULL
  }
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
  # USE spida2::lssvd
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

  #
  library(car)
  library(spida2)
  ret <- list()
  ret[['names']] <- list(gname = g, gcomplevel = comp)

  #
  # Refit with data set used for full model to check for common data set
  #

  full <- fitl[[length(fitl)]]
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
  # actual resids from comparator
  #

  groupL <- t( 1* outer(d[[g]], levels(d[[g]]), '=='))  # indicator matrix for groups
  groupL <- groupL / apply(groupL, 1, sum)
  ret[['groupL']] <- groupL

  ## Individual residuals from each fit

  gresids <- lapply(fitl, function(f) {
    y <- predict(f) + resid(f)
    y - predict(f, newdata = dg)
  })
  inds <- rep(1:nrow(d), length(gresids))
  dout <- d[inds,]
  dout$gresids <- unlist(gresids)
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
  data <- expand.grid(zork = levels(d[[g]]), model = names(fitl))
  names(data)[1] <- g
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


  ret[['gaps']] <- gaps
  ret[['resids']] <- resids
  ret[['resids2']] <- resids2
  # disp(resids2)
  # disp(unlist(resids2) > 10^(-11))
  if(any(unlist(resids2) > 10^(-11))) warning('Some models may not be nested in last model. See resids2')

  ret[['pred']] <- pred

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
  ret[['gaps_each']] <-Ls_each %>% lapply(subset, select = -L) %>% do.call(rbind,.)



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
  Ldiffs <- kronecker(diffmat(length(fitl)), diag(length(levels(d[[g]])))) %*% Ls
  ret[['Ldiffs']] <- Ldiffs
  data <- expand.grid(zork = levels(d[[g]]), gapdiffs = paste(names(fitl)[-1], '-', names(fitl)[-length(names(fitl))]))
  names(data)[1] <- g
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
}                                                           ## END decomp ####
