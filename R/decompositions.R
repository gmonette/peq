# PROBLEMS: fixed?
# - output in z$gaps obtained through wald test is not wholly consistent
#   with regression gaps while taking means of gresid appears to
#   be. This might have something to do with estimation from a full
#   model and numerical errors. Needs investigating. In the meantime
#   gapplot uses gresids in z$dout
#   Note that resplot already does
if(FALSE){

  X <- cbind(1, 1:100, (1:100)^2, 1:100)
  # X <- cbind(1, 1:100, (1:100)^2)
  y <- cbind(1:100)
  Y <- cbind(1:100) %*% rbind(rep(1,1000))
  Y[] <- Y + rnorm(100*1000)
  dim(Y)

  Bf <- lsfit(X,Y, intercept =FALSE)$coef

  Bs <- peq::lssvd(X,Y,zero = 10^(-7))$coef
  B1ret <- lssvd(X,Y)
  B1ret[c('d','p','rank','zero')]
  B1 <- B1ret$coef
  dim(Bf)
  dim(Bs)
  var(t(Bf)) %>% svd(nu=0,nv=0)
  var(t(Bs)) %>% svd(nu=0,nv=0)
  var(t(B1)) %>% svd(nu=0,nv=0)
  Bs %>% t %>% head
  plot(t(Bs)[,c(2,4)])
  plot(t(Bf)[,c(2,4)])


}

if(FALSE){                                             ## Sample data ----------
  {
    library(spida2)
    library(latticeExtra)
    library(nlme)
    library(car)
    set.seed(123)
    lnt <- function(x) 100 * log(x)
    et <- function(x) exp(x/100)
    expand.grid(age= seq(30,70), g = c('A','B','C'), area = c('a','b','c')) %>%
      within({
        p <- 1/(1 + exp(-scale(age*(as.numeric(g)-2)*(as.numeric(area)-2))))
        n <- rpois(length(p), abs(p+3))
        y <- age -.0002*age^2 + as.numeric(g) + as.numeric(area) + rnorm(n)
        ages <- (age -50)/10
      }) -> d1
    dd <- with(d1, d1[rep(1:nrow(d1), n),])
    dd$id <- 1:nrow(dd)
    dim(dd)
    head(dd)
    sq2 <- function(x) gsp(x, c(-1,1), c(1,2,1), c(1,1,1))[,-1]
    list() %>%
      within(
        {
          fit1  <- lm(lnt(y) ~ g, dd)
          fit2 <- lm (lnt(y) ~ g * ages + sq2(ages), dd)
          fit3 <- lm (lnt(y) ~ g + ages * area + sq2(ages), dd)
          full <- lm (lnt(y) ~ g * ages * area * sq2(ages),dd)

        }
      ) %>% rev -> fitl
  }                                      ## RUN sample data ----
}
#' Graphical comparison of AIC and BIC
#'
#' Compares AIC and BIC for a list of models. The models should have been fitted
#' to the same data.
#'
#' @param ll a list of models using the same data.
#' @param rot angle to rotate labels
#'
#' @export
icp <- function(ll, ..., rot = 0, srt = rot){
  if(is.null(names(ll))) {
    names(ll) <- paste0("model_",seq_along(ll))
  }
  a <- aic(ll)
  b <- bic(ll)
  a$BIC <- b$BIC
  a <- sortdf(a, ~ df)
  tobj <- xyplot(AIC +BIC~ df, a, type = 'l', outer = T,...,
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
#' Decompose pay gaps suing a sequence of models
#'
#' Pay equity gaps between equity-seeking groups and a comparator
#'
#' @param fitl a list of fitted models. All models should be nested in the last.
#'        Usually each model is nested in the next model.
#' @param g name (character) of group variable.
#' @param comp level of variable 'g' to use as comparator.
#' @param data used for fit model on common data. Needed for nlme models
#' @param cond vector of names of variables for conditioning (decomp2 only)
#' @examples
#' \dontrun{
#' library(peq)
#' library(spida2)
#' }
#' @export
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

#' @describeIn decomp conditional version with 'cond' argument
#' @export
decomp2 <- function(fitl, g, comp, data = na.omit(getD(full)), cond = NULL, refit = TRUE) {
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
}                                                           ## END decomp2 ####


#' svd version of lsfit
#'
#' Works with rank-deficient models (currently not centering)
#'
#' @param x predictor matrix, should include intercept term if needed
#' @param y response vector or matrix
#' @param zero value below which a latent value of the data matrix is considered to be 0, default 10^(-7) in parallel with \code{\link{lsfit}}.
#'
#' @export
lssvd <- function(x, y, zero = 1e-07,...) {
   spida2::lssvd(x, y, zero = zero, ...)
}
#' @describeIn lssvd oldversion
#' @export
lssvd_old <- function(x, y, zero = 10^(-7), has_intercept = all(cbind(x)[,1] ==1)) {
  lssvd_nc <- function(x,y, zero) {
    ret <- list()
    # x and y should be matrices
    xp <- svd(x, nu = ncol(x), nv = ncol(x))
    uy <- t(xp$u)%*%y
    dinv <- 1/xp$d
    dinv[abs(xp$d) < zero] <- 0
    ret$coefficients <- t(t(xp$v)*dinv) %*% uy
    ret$d <- xp$d
    ret$zero <- zero
    ret$p <- ncol(x)
    ret$rank <- sum(xp$d > zero)
    ret
  }
  disp <- function(...) NULL
  # if(!has_intercept) stop('has_intercept FALSE not yet implemented.')
  # if(!all(x[,1] == 1)) stop('first column must be intercept term.')
  x <- cbind(x)
  y <- cbind(y)
  xorig <- x
  xn <- colnames(x)
  yn <- colnames(y)
  # if(has_intercept) {
  if(FALSE) {   ## NEEDS FIXING    -----
    x <- x[,-1, drop = FALSE]
    xc <- scale(x)
    yc <- scale(y)
    yc[is.na(yc)] <- 0
    xm <- apply(x,2, mean)
    ym <- apply(y,2, mean)
    xs <- apply(x,2, sd)
    ys <- apply(y,2, sd)
    B <- lssvd_nc(xc, yc, zero = zero)
    disp(B)
    B <- ys * t(t(B)/c(xs))
    disp(B)
    disp(xm)
    disp(ym)
    Beta <- rbind( ym - rbind(xm) %*% B, B)
  }
  else ret <- lssvd_nc(x, y, zero = zero)
  Beta <- ret$coefficients
  colnames(Beta) <- yn
  rownames(Beta) <- xn
  ret$coefficients <- Beta
  ret$resids <-  y - xorig%*%Beta
  ret$sse <- apply(ret$resids, 2, function(x) sum(x^2))
  ret
}
if(FALSE) {
  zy <- cbind(1, 1:100, (1:100)^2-((1:100))^3)
  zx <- cbind(1,-50:49,(-50:49)^2,(-50:49)^3)
  lssvd(zx, zy)
  (zy - zx %*% lssvd(zx, zy, zero = 10^(-18))$coef) %>% head
  (zy - zx %*% lssvd(zx, zy, zero = 10^(-16))$coef) %>% head
  (zy - zx %*% lssvd(zx, zy, zero = 10^(-16), has_intercept=F)$coef) %>% head
  lsfit(zx, zy, intercept = FALSE)$resid %>% head
}
#' Plotting a decomp object
#'
#' @param decomp a decomp object creaate by \code{\link{decomp}}.
#'
#' @export
gapplot <- function(obj, data = obj$gaps_each, log = FALSE, rot = 45,
                    at = seq(-200000,100000,10000),
                    auto.key = list(space = 'right'), ...) {
# gapplot <- function(obj, data = obj$dout, log = FALSE, rot = 45,
#                     at = seq(-200000,100000,10000),...) {
  library(latticeExtra)
  library(spida2)
  disp <- function(...) {
    NULL
  }
  # Used gresids but then reverted to gaps_each after checking consistency
  dollar_gap <- function(x, reference) {
    et <- function(x) exp(x/100)
    lnt <- function(x)  100 * log(x)
    et(x + reference) - et(reference)
  }
#  data$coef <- with(data, capply(gresids, data[,c(obj$names$gname,'model')], mean, na.rm = T))
#  form <- as.formula(paste('~', obj$names$gname, "+ model"))
#  disp(dim(data))
#  data <- up(data, form)
#  data <- sortdf(data, ~ model)
#  print_data <- data
#  print_data$coef <- fmt(print_data$coef)
#  disp(print_data)
#  if(log) stop('log might not be working correctly dt dout fix')
  if(log) {
    disp(obj$names$gname)
    disp(names(obj$pred))
    disp( obj$pred[[obj$names$gname]])
    disp(obj$names$complevel)
    sel <- obj$pred[[obj$names$gname]] == obj$names$gcomplevel
    disp(sel)
    ref <- obj$pred[sel,'coef'][1]
    disp(ref)
    disp(data$coef)
    data$coef <- dollar_gap(data$coef, ref)
    disp(data$coef)
  }
  disp(data$coef)
  data$gr_ <- with(data, reorder(factor(data[[obj$names$gname]]), - coef))
  fmla <- "coef ~ model"
  if(!is.null(obj$names$cond)) {
    fmla <- paste(fmla, "|" , paste(obj$names$cond, collapse = '*'))
  }
  fmla <- as.formula(fmla)

  xyplot(fmla,
         data, ...,
#         groups = data[[obj$names$gname]], type = 'b',
         groups = gr_, type = 'b',
         scales = list(y =
                         list(at = at,
                              labels = fmt(at,0)),
                       x = list(rot = rot)),
         ylab = 'Group-weighted adjusted salary gap\nfrom comparator group',
         xlab = "Cumulatively adjusted factors",
         auto.key = auto.key)+
    layer_(panel.grid(v=-1,h=-1))

}
#' boxplot statistics for resplot
#'
#' @export
    qstats_resplot <- function (x, coef = 0.05, do.conf = TRUE, do.out = TRUE, ends = .05) {
      # adapted from grDevices::boxplot.stats to return quantile whiskers
      disp <- function(...) NULL
      disp('in qstats')
      nna <- !is.na(x)
      n <- sum(nna)
      stats <- stats::fivenum(x, na.rm = TRUE)
      stats[c(1,5)] <- quantile(x, c(ends, 1-ends), na.rm = TRUE)
      big <- na2f(x > max(stats))
      small <- na2f(x < min(stats))
      disp(big)
      out <- x[big | small]
      # iqr <- diff(stats[c(2, 4)])
      #out <- (x > stats[5]) | (x < stats[1])
      list(stats = stats, n = n, conf = FALSE, out = if (length(out)>0) out else numeric())
    }

#' Plot individual residuals
#'
#'
#' @export
resplot <- function(obj, data = obj$dout, log = FALSE, which = 1,
                    at = seq(-200000,100000,10000),..., rot = 45,
                    par.strip.text = list(cex = 1)) {
    disp <- function(...) NULL
    library(latticeExtra)
    library(spida2)
    env <- environment()

  # qstats <- boxplot.stats
  dollar_gap <- function(x, reference) {   # also in gapplot
    et <- function(x) exp(x/100)
    lnt <- function(x)  100 * log(x)
    et(x + lnt(reference)) - reference
  }
  if(log) {
    disp(obj$names$gname)
    disp(names(obj$pred))
    disp( obj$pred[[obj$names$gname]])
    disp(obj$names$complevel)
    sel <- obj$pred[[obj$names$gname]] == obj$names$gcomplevel
    disp(sel)
    ref <- obj$pred[sel,'coef'][1]
    disp(ref)
    dim(head(data))
    #disp(data$gresids)
    #disp(dollar_gap(data$gresids, ref))
    data$gresids <- dollar_gap(data$gresids, ref)
  }
if(which ==1)  xyplot(gresids ~ model | data[[obj$names$gname]],
         data, ..., type = 'p', pch = '.', cex = 1, #alpha = .4,
         par.strip.text = par.strip.text,
         # groups = data[[obj$names$idname]],
         scales = list(y =
                         list(at = at,
                              labels = fmt(at,0),
                              alternating = 2),
                       x = list(alternating = F, rot = rot)),
         ylab = 'Individual salary differences\nfrom comparator prediction',
         xlab = "Cumulatively adjusted factors",
         auto.key = list(space = 'right'))+
    # layer_(panel.grid(v=-1,h=-1))+
    layer(panel.bwplot(..., horizontal = FALSE, fill = 'grey90', lwd = 2,
                       varwidth = T,
                       pch = '|',coef = 6, lty = 1,
                       stats = qstats_resplot)) +
    layer(panel.abline(h=0))
else  xyplot(gresids ~ data[[obj$names$gname]] | model,
         data, ..., type = 'p', pch = '.', cex = 1, #alpha = .4,
         # groups = data[[obj$names$idname]],
         par.strip.text = par.strip.text,
         scales = list(y =
                         list(at = at,
                              labels = fmt(at,0),
                              alternating = 2),
                       x = list(alternating = F, rot = rot)),
         ylab = 'Individual salary differences\nfrom comparator prediction',
         xlab = "Groups",
         auto.key = list(space = 'right'))+
    layer(panel.bwplot(...,
                       horizontal = FALSE,
                       fill = 'grey90',
                       lwd = 2,
                       varwidth = T,
                       pch = '|',
                       coef = .1,
                       lty = 1,
                       stats = qstats.resplot)) +
    layer_(panel.grid(v=-1,h=-1))+
    layer(panel.abline(h=0))
} # ;  resplot(z, at = seq(-100,110,10),log = T,which = 2, data = subset(z$dout, !model %in%c('full','fit1'))) # end of resplot                                        ## RUN --------------

## testing ####

# test decomp2



