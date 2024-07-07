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

##
## decomp:  ####
##
#'
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
#' library(latticeExtra)
#' library(kableExtra)
#' simdata <- read.table(header = TRUE, stringsAsFactors = TRUE,text ="
#' sal  gender  area   age  rank
#' 100  F       A      30   a   # group A: gap: -10 and -20 for A
#' 110  F       A      40   b
#' 120  F       A      50   b
#' 120  M       A      40   a
#' 130  M       A      50   b
#' 140  M       A      60   b
#' 100  A       A      40   a
#' 110  F       B      30   a # group B: gap: -20 and -30 for A
#' 130  M       B      30   a
#' 140  M       B      40   a
#' 140  M       B      40   b
#' 149  M       B      50   b
#' 150  M       B      50   c
#' 151  M       B      50   c
#' 150  M       B      50   b
#' 150  M       B      50   b
#' 150  M       B      50   c
#' 150  M       B      50   c
#' 150  M       B      50   c
#' 150  M       B      50   b
#' 120  A       B      50   b
#' ")
#' simdata$gender <- relevel(simdata$gender, 'M')
#' tps(pch=16)
#' xyplot(sal ~ age | area, simdata, groups = gender, auto.key = TRUE)
#'
#'
#' contr.helmert(3)
#' contr.rot <- function(n) {
#'    ret <- contr.helmert(n)
#'    disp(ret)
#'    ret <- t( t(ret)/sqrt(apply(ret^2, 2, sum)))
#'    ret
#' }
#' contr.rot(3)
#' contr.rot(3) %>% crossprod
#' ?C
#'
#'
#' contrasts(simdata$gender)
#'
#' (fit0 <- lm(sal ~ gender, simdata)) %>% summary
#' (fit1 <- lm(sal ~ gender + age, simdata)) %>% summary
#' (fit2 <- lm(sal ~ gender + age + area, simdata)) %>% summary
#' (fit3 <- lm(sal ~ area/gender + age -1 , simdata)) %>% summary
#' (fit4 <- lm(sal ~ area/gender + age + rank -1 , simdata)) %>% summary
#'
#' fitl <- list(gender = fit0, age = fit1, area = fit2, g_by_a = fit3, rank = fit4)
#'
#' pred <- with(simdata, spida2::pred.grid(gender,age,area))
#' pred$fit3 <- predict(fit3, newdata = pred)
#'
#' xyplot(fit3 ~ age | area, pred, groups = gender, type = 'l',auto.key = TRUE) +
#' xyplot(jitter(sal,10) ~ age | area, simdata, groups = gender, auto.key = TRUE)
#'
#' # Using fit1: no effect of area
#' L1 <- rbind(Agap = c(0,1,0,0),
#'             Fgap = c(0,0,1,0))
#' waldf(fit1, L1)
#'
#' # Using fit2: additive area
#' L2 <- rbind(Agap = c(0,1,0,0,0),
#'             Fgap = c(0,0,1,0,0))
#' waldf(fit2, L2)
#'
#' # Using fit3: no effect of area
#' L3 <- rbind("Agap in A" = c(0,0,0,1,0,0,0),
#'             "Fgap in A" = c(0,0,0,0,0,1,0),
#'             "Agap in B" = c(0,0,0,0,1,0,0),
#'             "Fgap in B" = c(0,0,0,0,0,0,1))
#' waldf(fit3, L3)
#'
#' SEs <-   waldf(fit3, L3)$se
#' wts <- 1/SEs^2   # precision weights
#'
#' # Type III average over Areas
#'
#' W3 <- rbind( A3 = std(c(1, 0 , 1, 0)),
#'              F3 = std(c(0, 1 , 0, 1)))
#' W3
#'
#' # Average by size of Areas  (dubious)
#'
#' nArea <- tab__(simdata, ~ area)
#' nArea <- nArea[c(1,1,2,2)]
#'
#' Wsize <- rbind(
#'     AnA = std(nArea * c(1,0,1,0)),
#'     FnA = std(nArea * c(0,1,0,1))
#' )
#' Wsize
#'
#' # Average by size of equity group
#'
#' ## Compare relative weights of different methods
#' ## Idea: Using additive model may give reasonable weights if no group is very small
#' ## but bad idea for very small groups were the gap might not reflect the
#' ## average **experience** of members of the small group.
#' ## Weights will tend to look like (1/nF + 1/nM)^1 which will look like nF if small ..... so ....
#' ## ... but maybe not so simple with other adjustment variables like age, etc....
#'
#' ns <- tab__(simdata, ~ gender + area)
#' ns
#' Weqn <- rbind(
#'     AnA = std(c(ns['A','A'],0,ns['A','B'],0)),
#'     FnA = std(c(0,ns['F','A'],0,ns['F','B']))
#' )
#' Weqn
#'
#' waldf(fit3, W3 %*% L3)
#' waldf(fit3, Wsize %*% L3)
#' waldf(fit3, Weqn %*% L3)
#'
#' decomp(fitl, 'gender', 'M') %>% gapplot
#' decomp2(fitl, 'gender', 'M', cond = 'area') %>% gapplot
#'
#' decomp(fitl, 'gender', 'M') %>% decomp_table(log=FALSE)
#' decomp(fitl, 'gender', 'M') %>% decomp_table(log=TRUE)
#' decomp2(fitl, 'gender', 'M', cond = 'area') %>% gapplot
#' decomp2(fitl, 'gender', 'M', cond = 'rank') %>% gapplot
#'
#' # TODO: Use SEs for fit3 to show how precision weighted means are same as
#' # gaps obtained from non-interaction model
#'
#' decomp2(fitl, 'gender','M', cond = c('area','rank')) %>% gapplot
#'
#' #
#' #
#' # Simple example with age and gender show
#' # possibility of a reversal depending
#' # on the choice of comparator group
#' #
#' dage <- read.table(header = TRUE, stringsAsFactors = TRUE,text ="
#' Gender    Age   Salary
#' M         30    30000
#' M         50    50000
#' M         60    60000
#' M         65    65000
#' M         65    65000
#' M         70    69000
#' M         70    70000
#' M         70    71000
#' F         30    30000
#' F         40    35000
#' F         40    36000
#' F         40    34000
#' F         50    40000
#' ")
#'
#'
#' tps(pch=16:17)
#' xyplot(Salary ~ Age, dage, groups = Gender, alpha = .5)
#' fit1 <- lm(Salary ~ Gender, dage)
#' fit2 <- lm(Salary ~ Gender + Age, dage)
#' fit3 <- lm(Salary ~ Gender * Age, dage)
#' fitl <- list(Gender=fit1, 'Gender + Age'=fit2, 'Gender * Age'= fit3)
#' icp(fitl)
#' decomp(fitl, 'Gender', 'M')
#' decomp(fitl, 'Gender', 'M') %>% decomp_table
#' decomp(fitl, 'Gender', 'F') %>% gapplot
#' decomp(fitl, 'Gender', 'M') %>% gapplot
#'
#'
#' #
#' # The following is example of a reversal paradox in which
#' # the direction of the gap is reversed depending on the
#' # choice of reference groups. This can only occur with
#' # interaction with Gender in the model and
#' #
#'
#'
#' dage2 <- within(
#'   dage,
#'   {
#'     Salary2 <- ifelse(Gender == 'F', Salary + 10000, Salary)
#'   }
#' )
#' fit21 <- lm(Salary2 ~ Gender, dage2)
#' fit22 <- lm(Salary2 ~ Gender + Age, dage2)
#' fit23 <- lm(Salary2 ~ Gender * Age, dage2)
#' fit2l <- list(Gender=fit21, 'Gender + Age'=fit22, 'Gender * Age'= fit23)
#' decomp(fit2l, 'Gender', 'F') %>% decomp_table
#' decomp(fit2l, 'Gender', 'M') %>% decomp_table
#' decomp(fit2l, 'Gender', 'F') %>% gapplot
#' decomp(fit2l, 'Gender', 'M') %>% gapplot
#'
#'
#'
#'
#' }
#' @export
decomp <- function(fitl, g, comp, data = na.omit(getD(full)), refit = TRUE,
                   verbose = FALSE) {
  if(verbose == FALSE) disp <- function(...) {
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

  addrownames <- function(x) {
    rownames(x) <- as.character(seq_len(nrow(x)))
    x
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
  ret[['gaps_each']] <- Ls_each %>% lapply(subset, select = -L)  %>%
    do.call(rbind,.) %>% addrownames


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
}

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
#' @import utils
#' @import stats
#' @import lattice
#' @import latticeExtra
#' @export
gapplot <- function(obj, data = obj$gaps_each, log = FALSE, rot = 45,
                    ylab = 'Group-weighted adjusted salary gaps',
                    # ylab = 'Group-weighted adjusted salary gaps\nfrom comparator group',
                    xlab = 'Cumulatively adjusted factors',
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
  orig_gap <- function(gap, pred) {
    # gap and pred in 100 x ln(salary)
    # pred is predicted for group, not comparator
    exp(pred/100) * (1 - exp(-gap/100))
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
  if(FALSE) {  # won't work if subset selected
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

  if(log) {  # should work with subset if log
    disp(dim(data))
    disp(head(data))
    disp(dim(obj$pred))
    disp(head(obj$pred))
    pred <- obj$pred[rownames(data),]
    data$coef <- orig_gap(data$coef,pred$coef)
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
                         list(at = at, alternating = 1,
                              labels = fmt(at,0)),
                       x = list(rot = rot, alternating = 1)),
         ylab = ylab,
         xlab = xlab,
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
} #  resplot(z, at = seq(-100,110,10),log = T,which = 2, data = subset(z$dout, !model %in%c('full','fit1'))) # end of resplot                                        ## RUN --------------

## testing ####

# test decomp2
#
#' Table of gaps, percent gaps from model fitting log salaries
#'
#' @param z output of decomp2 or decomp1 (untested)
#'
#' @export
decomp_table <- function(x, log = TRUE, p = TRUE, n_min = 1) {
  if(log) {
    decomp_table_log(x, p = p, n_min = n_min)
  }
  else {
    decomp_table_raw(x, p = p, n_min = n_min)
  }
}
#' @export
decomp_table_log <- function(z, p = TRUE, n_min = 1, reduction = "Change"){
  #
  # This is the version for log models with p-values for gaps
  #
  # z is the result of using decomp or decomp2 on a list of models
  #
  orig_gap <- function(gap, pred) {
    # gap and pred in 100 x ln(salary)
    # pred is predicted for group, not comparator
    exp(pred/100) * (1 - exp(-gap/100))
  }
  sel <- function(arr, ind) {
    do.call(`[`,c(list(arr), ind))
  }
  fmtpc <- function(x, nsmall) {
    ret <- fmt(x, nsmall)
    ret <- paste0(ret,"\\percent" )
    ret
  }
  clean <- function(x) {
    #where <- x
    x[] <- gsub("-Inf","  NA", x)
    where <- grepl('Inf|NaN', x)
    x[where] <- ''
    where <- grepl('Inf|NaN', x)
    x[where] <- ''
    x
  }
  #
  # need to merge predictions and gaps
  # allgaps is function of g, model, cond vars and type
  # pred is function of g, model, cond vars
  #
  # Easier to work on gaps each and take differences
  #
  # If we need p-values then we'll need to use gapdiffs
  #
  #
  gaps <- z$gaps_each
  pred <- z$pred
  fmla <- as.formula(paste('coef ~ ', paste(c(z$names$gname,"model", z$names$cond), collapse = '+')))
  fmla_n <-  as.formula(paste(' ~ ', paste(c(z$names$gname, z$names$cond), collapse = '+')))
  fmla_p <- as.formula(paste('`p-value` ~ ', paste(c(z$names$gname,"model", z$names$cond), collapse = '+')))  #new

  gap_tab <- tab__(fmla, gaps)
  pred_tab <- tab__(fmla, pred)
  n_tab <- tab__(fmla_n, z$data)
  p_tab <- tab__(fmla_p, gaps)   # new

  value_tab <- gap_tab
  value_tab[] <- orig_gap(gap_tab, pred_tab)    # percent to dollar value
  p_tab[] <- pfmt(p_tab)
  p_tab[] <- paste0("p:",p_tab)

  # Work out reductions by using difference

  dims <- dim(value_tab)
  inds1 <- indslast <- inds <- lapply(dims, seq_len)
  indslast[[2]] <- indslast[[2]][-1]
  inds1[[2]] <- inds1[[2]][-length(inds1[[2]])]

  gap_diffs <- do.call('[', c(list(gap_tab), indslast)) -
    do.call('[', c(list(gap_tab), inds1))

  val_diffs <- do.call('[', c(list(value_tab), indslast)) -
    do.call('[', c(list(value_tab), inds1))

  # set up indices for large table combining gaps and reductions

  dims_ret <- dims
  dims_ret[2] <- 2* dims_ret[2] # one column for N, dims_ret[2] for gaps
  # and dim_ret[2] - 1 for reductions
  inds <- lapply(dims_ret, seq_len)

  ret1 <- array('', dims_ret)    # row 1 of each cell

  # add N in row 1 first column
  #
  ii <- inds
  ii[[2]] <- 1

  ret1 <- do.call("[<-", c(list(ret1),ii, list(fmt(n_tab,0))))

  # add gap pcts in row 1 gap columns

  ii[[2]] <- seq(2,dim(ret1)[2],2)
  ret1 <- do.call("[<-", c(list(ret1),ii, list(fmtpc(gap_tab,1))))

  # add reductions in row 1 diff columns

  ii[[2]] <- seq(3,dim(ret1)[2],2)
  ret1 <- do.call("[<-", c(list(ret1),ii, list(fmtpc(gap_diffs,1))))


  # Row 2 add values

  ret2 <- array('', dims_ret)    # row 2 of each cell

  # Add gap dollar value in row gap columns

  ii[[2]] <- seq(2,dim(ret2)[2],2)
  ret2 <- do.call("[<-", c(list(ret2),ii, list(fmt(value_tab,0))))

  # Add reductions in dollar value in row 2 reduction columns

  ii[[2]] <- seq(3,dim(ret1)[2],2)
  ret2 <- do.call("[<-", c(list(ret2),ii, list(fmt(val_diffs,0))))

  if(p) {

    # make row 3 with p-values

    ret3 <- array(' ', dims_ret)     # row 3 of each cell (no-empty for prevent vertical centering)
    ii[[2]] <- seq(2,dim(ret3)[2],2)
    ret3 <- do.call("[<-", c(list(ret3),ii, list(p_tab)))

  }
  ret <- ret1
  # disp(ret1[,1])
  # disp(n_tab)
  if(p) ret[] <- kbind(ret1, '\n', ret2, '\n', ret3)
  else ret[] <- kbind(ret1, '\n', ret2)


  dimnames(ret)[[1]] <- dimnames(gap_tab)[[1]]
  if(length(ii) > 2){
    dimnames(ret)[seq(3,length(ii))] <- dimnames(gap_tab)[seq(3,length(ii))]
  }

  droplast <- function(x) x[-length(x)]
  dimnames(ret)[[2]] <- droplast(c("N",rbind(dimnames(gap_tab)[[2]],reduction)))
  names(dimnames(ret)) <- names(dimnames(gap_tab))
  # class(ret) <- 'decomp_table_log'
  ret <- as.table(ret)

  # blank if n < n_min

  if(TRUE) {

  retkeep <- ret
  np <- seq_along(dim(ret))
  np <- c(np[-2], 2)
  npinv <- np
  npinv[np] <- seq_along(np)
  retkeep <- aperm(retkeep, np)
  retkeep[] <- c(n_tab)
  retkeep <- aperm(retkeep, npinv)
  retkeep <- (retkeep >= n_min) | (slice.index(retkeep, 2) == 1)
  ret[!retkeep] <- ' '
  }


  clean(ret)
}

# decomp_table(z, log = T)

# %>% k_(' ') %>% ksf


# decomp_table_log(z, p = F) %>% k_(' ') %>% ksf


# new version with p-value option
#' @export
decomp_table_raw <- function(z, p = TRUE, n_min = 1, reduction = 'Change'){
  #
  # This is the version for models were raw gaps are not transformed to dollar gaps
  #
  # z is the result of using decomp or decomp2 on a list of models
  #
  # Not needed:
  orig_gap <- function(gap, pred) {
    # gap and pred in 100 x ln(salary)
    # pred is predicted for group, not comparator
    exp(pred/100) * (1 - exp(-gap/100))
  }
  sel <- function(arr, ind) {
    do.call(`[`,c(list(arr), ind))
  }
  # Not needed:
  #
  # fmtpc <- function(x, nsmall) {
  #   ret <- fmt(x, nsmall)
  #   ret <- paste0(ret,"\\percent" )
  #   ret
  # }
  clean <- function(x) {
    #where <- x
    x[] <- gsub("-Inf","  NA", x)
    where <- grepl('Inf|NaN', x)
    x[where] <- ''
    where <- grepl('Inf|NaN', x)
    x[where] <- ''
    x
  }
  #
  # need to merge predictions and gaps
  # allgaps is function of g, model, cond vars and type
  # pred is function of g, model, cond vars
  #
  # Easier to work on gaps each and take differences
  #
  # If we need p-values then we'll need to use gapdiffs
  #
  # Still need to work out 'value' because it's used later for dimensions, etc.
  #
  gaps <- z$gaps_each
  pred <- z$pred
  fmla <- as.formula(paste('coef ~ ', paste(c(z$names$gname,"model", z$names$cond), collapse = '+')))
  fmla_n <-  as.formula(paste(' ~ ', paste(c(z$names$gname, z$names$cond), collapse = '+')))
  fmla_p <- as.formula(paste('`p-value` ~ ', paste(c(z$names$gname,"model", z$names$cond), collapse = '+')))  #new

  gap_tab <- tab__(fmla, gaps)
  pred_tab <- tab__(fmla, pred)
  n_tab <- tab__(fmla_n, z$data)
  p_tab <- tab__(fmla_p, gaps)   # new

  value_tab <- gap_tab
  value_tab[] <- orig_gap(gap_tab, pred_tab)    # percent to dollar value
  p_tab[] <- pfmt(p_tab)
  p_tab[] <- paste0("p:",p_tab)

  # Work out reductions by using difference

  dims <- dim(value_tab)
  inds1 <- indslast <- inds <- lapply(dims, seq_len)
  indslast[[2]] <- indslast[[2]][-1]
  inds1[[2]] <- inds1[[2]][-length(inds1[[2]])]

  gap_diffs <- do.call('[', c(list(gap_tab), indslast)) -
    do.call('[', c(list(gap_tab), inds1))

  val_diffs <- do.call('[', c(list(value_tab), indslast)) -
    do.call('[', c(list(value_tab), inds1))

  # set up indices for large table combining gaps and reductions

  dims_ret <- dims
  dims_ret[2] <- 2* dims_ret[2] # one column for N, dims_ret[2] for gaps
  # and dim_ret[2] - 1 for reductions
  inds <- lapply(dims_ret, seq_len)

  ret1 <- array('', dims_ret)    # row 1 of each cell

  # add N in row 1 first column
  #
  ii <- inds
  ii[[2]] <- 1

  ret1 <- do.call("[<-", c(list(ret1),ii, list(fmt(n_tab,0))))

  # add gap raw in row 1 gap columns

  ii[[2]] <- seq(2,dim(ret1)[2],2)
  ret1 <- do.call("[<-", c(list(ret1),ii, list(fmt(gap_tab,0))))

  # add reductions in row 1 diff columns

  ii[[2]] <- seq(3,dim(ret1)[2],2)
  ret1 <- do.call("[<-", c(list(ret1),ii, list(fmt(gap_diffs,0))))


  # # Row 2 add values
  #
  # ret2 <- array('', dims_ret)    # row 2 of each cell
  #
  # # Add gap dollar value in row gap columns
  #
  # ii[[2]] <- seq(2,dim(ret2)[2],2)
  # ret2 <- do.call("[<-", c(list(ret2),ii, list(fmt(value_tab,0))))
  #
  # # Add reductions in dollar value in row 2 reduction columns
  #
  # ii[[2]] <- seq(3,dim(ret1)[2],2)
  # ret2 <- do.call("[<-", c(list(ret2),ii, list(fmt(val_diffs,0))))

  if(p) {

    # make row 3 with p-values

    ret3 <- array(' ', dims_ret)     # row 3 of each cell (no-empty for prevent vertical centering)
    ii[[2]] <- seq(2,dim(ret3)[2],2)
    ret3 <- do.call("[<-", c(list(ret3),ii, list(p_tab)))

  }
  ret <- ret1
  if(p) ret[] <- kbind(ret1, '\n', ret3)
  else ret[] <- kbind(ret1)


  dimnames(ret)[[1]] <- dimnames(gap_tab)[[1]]
  if(length(ii) > 2){
    dimnames(ret)[seq(3,length(ii))] <- dimnames(gap_tab)[seq(3,length(ii))]
  }

  droplast <- function(x) x[-length(x)]
  dimnames(ret)[[2]] <- droplast(c("N",rbind(dimnames(gap_tab)[[2]],reduction)))
  names(dimnames(ret)) <- names(dimnames(gap_tab))
  # class(ret) <- 'decomp_table_log'
  ret <- as.table(ret)

  # blank if n < n_min

  if(TRUE){
    retkeep <- ret
    np <- seq_along(dim(ret))
    np <- c(np[-2], 2)
    npinv <- np
    npinv[np] <- seq_along(np)
    retkeep <- aperm(retkeep, np)
    retkeep[] <- c(n_tab)
    retkeep <- aperm(retkeep, npinv)
    retkeep <- (retkeep >= n_min) | (slice.index(retkeep, 2) == 1)
    ret[!retkeep] <- ' '
  }

  clean(ret)
}



decomp_table_log_original <- function(z){
  # z is the result of using decomp or decomp2 on a list of models
  orig_gap <- function(gap, pred) {
    # gap and pred in 100 x ln(salary)
    # pred is predicted for group, not comparator
    exp(pred/100) * (1 - exp(-gap/100))
  }
  sel <- function(arr, ind) {
    do.call(`[`,c(list(arr), ind))
  }
  fmtpc <- function(x, nsmall) {
    ret <- fmt(x, nsmall)
    ret <- paste0(ret,"\\percent" )
    ret
  }
  #
  # need to merge predictions and gaps
  # allgaps is function of g, model, cond vars and type
  # pred is function of g, model, cond vars
  #
  # Easier to work on gaps each and take differences
  #
  # If we need p-values then we'll need to use gapdiffs
  #
  #
  gaps <- z$gaps_each
  pred <- z$pred
  fmla <- as.formula(paste('coef ~ ', paste(c(z$names$gname,"model", z$names$cond), collapse = '+')))
  fmla_n <-  as.formula(paste(' ~ ', paste(c(z$names$gname, z$names$cond), collapse = '+')))
  fmla_p <- as.formula(paste('`p-value` ~ ', paste(c(z$names$gname,"model", z$names$cond), collapse = '+')))

  gap_tab <- tab__(fmla, gaps)
  pred_tab <- tab__(fmla, pred)
  n_tab <- tab__(fmla_n, z$data)
  p_tab <- tab__(fmla_p, gaps)

  value_tab <- gap_tab
  value_tab[] <- orig_gap(gap_tab, pred_tab)

  dims <- dim(value_tab)
  inds1 <- indslast <- inds <- lapply(dims, seq_len)
  indslast[[2]] <- indslast[[2]][-1]
  inds1[[2]] <- inds1[[2]][-length(inds1[[2]])]

  gap_diffs <- do.call('[', c(list(gap_tab), indslast)) -
    do.call('[', c(list(gap_tab), inds1))

  val_diffs <- do.call('[', c(list(value_tab), indslast)) -
    do.call('[', c(list(value_tab), inds1))

  dims_ret <- dims
  dims_ret[2] <- 2* dims_ret[2] # one column for N, dims_ret[2] for gaps
                                # and dim_ret[2] - 1 for reductions
  inds <- lapply(dims_ret, seq_len)

  ret1 <- array('', dims_ret)    # row 1 of each cell

  # add N in row 1 first column
  ii <- inds
  ii[[2]] <- 1

  ret1 <- do.call("[<-", c(list(ret1),ii, list(fmt(n_tab,0))))

  # add gap pcts in row 1 gap columns

  ii[[2]] <- seq(2,dim(ret1)[2],2)
  ret1 <- do.call("[<-", c(list(ret1),ii, list(fmtpc(gap_tab,1))))
  ii[[2]] <- seq(3,dim(ret1)[2],2)
  ret1 <- do.call("[<-", c(list(ret1),ii, list(fmtpc(gap_diffs,1))))


  # Row 2 add values

  ret2 <- array('', dims_ret)    # row 2 of each cell

  ii[[2]] <- seq(2,dim(ret2)[2],2)
  ret2 <- do.call("[<-", c(list(ret2),ii, list(fmt(value_tab,0))))
  ii[[2]] <- seq(3,dim(ret1)[2],2)
  ret2 <- do.call("[<-", c(list(ret2),ii, list(fmt(val_diffs,0))))
  dim(ret1)
  dim(ret2)
  ret <- ret1
  ret[] <- kbind(ret1,'\n',ret2)

  dimnames(ret)[[1]] <- dimnames(gap_tab)[[1]]
  if(length(ii) > 2){
    dimnames(ret)[seq(3,length(ii))] <- dimnames(gap_tab)[seq(3,length(ii))]
  }

  droplast <- function(x) x[-length(x)]
  dimnames(ret)[[2]] <- droplast(c("N",rbind(dimnames(gap_tab)[[2]],"Reduction")))
  names(dimnames(ret)) <- names(dimnames(gap_tab))
  # class(ret) <- 'decomp_table_log'
  as.table(ret)
}

# new version with p-value option

decomp_table_log_original2 <- function(z, p = TRUE){
  # z is the result of using decomp or decomp2 on a list of models
  orig_gap <- function(gap, pred) {
    # gap and pred in 100 x ln(salary)
    # pred is predicted for group, not comparator
    exp(pred/100) * (1 - exp(-gap/100))
  }
  sel <- function(arr, ind) {
    do.call(`[`,c(list(arr), ind))
  }
  fmtpc <- function(x, nsmall) {
    ret <- fmt(x, nsmall)
    ret <- paste0(ret,"\\percent" )
    ret
  }
  #
  # need to merge predictions and gaps
  # allgaps is function of g, model, cond vars and type
  # pred is function of g, model, cond vars
  #
  # Easier to work on gaps each and take differences
  #
  # If we need p-values then we'll need to use gapdiffs
  #
  #
  gaps <- z$gaps_each
  pred <- z$pred
  fmla <- as.formula(paste('coef ~ ', paste(c(z$names$gname,"model", z$names$cond), collapse = '+')))
  fmla_n <-  as.formula(paste(' ~ ', paste(c(z$names$gname, z$names$cond), collapse = '+')))
  fmla_p <- as.formula(paste('`p-value` ~ ', paste(c(z$names$gname,"model", z$names$cond), collapse = '+')))  #new

  gap_tab <- tab__(fmla, gaps)
  pred_tab <- tab__(fmla, pred)
  n_tab <- tab__(fmla_n, z$data)
  p_tab <- tab__(fmla_p, gaps)   # new

  value_tab <- gap_tab
  value_tab[] <- orig_gap(gap_tab, pred_tab)
  p_tab[] <- pfmt(p_tab)


  dims <- dim(value_tab)
  inds1 <- indslast <- inds <- lapply(dims, seq_len)
  indslast[[2]] <- indslast[[2]][-1]
  inds1[[2]] <- inds1[[2]][-length(inds1[[2]])]

  gap_diffs <- do.call('[', c(list(gap_tab), indslast)) -
    do.call('[', c(list(gap_tab), inds1))

  val_diffs <- do.call('[', c(list(value_tab), indslast)) -
    do.call('[', c(list(value_tab), inds1))

  dims_ret <- dims
  dims_ret[2] <- 2* dims_ret[2] # one column for N, dims_ret[2] for gaps
  # and dim_ret[2] - 1 for reductions
  inds <- lapply(dims_ret, seq_len)

  ret1 <- array('', dims_ret)    # row 1 of each cell

  # add N in row 1 first column
  ii <- inds
  ii[[2]] <- 1

  ret1 <- do.call("[<-", c(list(ret1),ii, list(fmt(n_tab,0))))

  # add gap pcts in row 1 gap columns

  ii[[2]] <- seq(2,dim(ret1)[2],2)
  ret1 <- do.call("[<-", c(list(ret1),ii, list(fmtpc(gap_tab,1))))
  ii[[2]] <- seq(3,dim(ret1)[2],2)
  ret1 <- do.call("[<-", c(list(ret1),ii, list(fmtpc(gap_diffs,1))))


  # Row 2 add values

  ret2 <- array('', dims_ret)    # row 2 of each cell

  ii[[2]] <- seq(2,dim(ret2)[2],2)
  ret2 <- do.call("[<-", c(list(ret2),ii, list(fmt(value_tab,0))))
  ii[[2]] <- seq(3,dim(ret1)[2],2)
  ret2 <- do.call("[<-", c(list(ret2),ii, list(fmt(val_diffs,0))))
  dim(ret1)
  dim(ret2)
  ret <- ret1
  ret[] <- kbind(ret1,'\n',ret2)

  ###########################3  ADD ROW 3 HERE


  dimnames(ret)[[1]] <- dimnames(gap_tab)[[1]]
  if(length(ii) > 2){
    dimnames(ret)[seq(3,length(ii))] <- dimnames(gap_tab)[seq(3,length(ii))]
  }

  droplast <- function(x) x[-length(x)]
  dimnames(ret)[[2]] <- droplast(c("N",rbind(dimnames(gap_tab)[[2]],"Reduction")))
  names(dimnames(ret)) <- names(dimnames(gap_tab))
  # class(ret) <- 'decomp_table_log'
  ret <- as.table(ret)
  # disp(dim(ret))
  # ret[n_tab == 0, -1] <-''
  clean(ret)
}
#' Multiple rows for kable
#'
#' @param x element in first row
#' @param ...  additional elements
#' @param prefix add before first element
#' @param align 'r', 'l' or 'c', default 'r'
#' @export
kbind <-
function(x,..., prefix = '', align = 'r') {
  # pasting tables together to possibly have multiple lines
  # in cells by using: kbind(x,'\n',y)
  # Note: preformat numbers to control format
  #
  # x can have a structure such as a table, matrix, etc.
  # If the first element does not have the desired structure
  # for the result, provide it with the 'prefix' argument.
  x[] <- kableExtra::linebreak(paste0(prefix,x,...), align = align)
  x
}

