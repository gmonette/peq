lssvd <- function(x,y, zero = 10^(-14)) {
  xp <- svd(x,nu=ncol(x), nv = ncol(x))
  uy <- t(xp$u)%*%cbind(y)
  dinv <- 1/xp$d
  dinv[abs(xp$d) < zero] <- 0
  coef <- t(t(xp$v)*dinv) %*% uy
  resid <- y - x %*% coef
  list(coef = coef, resid = resid, sse = sum(resid^2), p = ncol(x), rank = sum(abs(xp$d)>=zero), zero = zero)
}
zx <- cbind(1, 1:10, (1:10)^2, 10:1)
zy <- cbind(10:1)
lsfit(zx, zy, intercept = FALSE)
[c('coef','resid')]
lssvd(zx, zy)



# log fits
{
  library(spida2)
  library(latticeExtra)
  library(nlme)
  library(car)
}

sq3 <- function(x) gsp(x, c(-1,0,1.5), c(1,2,2,1), c(1,1,1))[,-1]
d9$Rank_ <- sub(" .Conditional.",'', d9$rank) %>% sub(', Teaching Stream','',.) 
d9$Rank_ <- reorder(factor(d9$Rank_), d9$salary)
tab(d9, ~ Rank_)
subset(d9, !is.na(Rank_) & Rank_ != 'Sr Lecturer') %>% droplevels -> d9_

within(
  list(),
  {
    raw <-         lm((salary) ~ Group6, d9_)
    age <-         lm((salary) ~ Group6 * ages + sq3(ages), d9_)
    term <- lm((salary) ~ (Term + Group6 + ages)^2 - Term:Group6 , d9_)
    stream <- lm((salary) ~ (Stream + Term + Group6 + ages)^2 
                 - Term:Group6 
                 + Stream*sq3(ages), d9_)
    area <-        lm((salary) ~ (Stream + Term + Group6 + Area +ages)^2
                      - Term:Group6 -Area:Group6 + Stream*sq3(ages), 
                      d9_)
    rank <-        lm((salary) ~ (Stream + Term + Group6 + Area + ages + Rank_)^2 
                      - Term:Group6 - Area:Group6 - Group6:Rank_ + Stream*sq3(ages)
                      - Stream:ages - Stream:Rank_ - ages:Rank_ , d9_)
  }) %>% rev -> fitlist
fitlist %>% icp  
fitlist %>% lapply(Anova)  
fitlist %>% lapply(summary)  


within(
  list(),
  {
    raw <-         lm((salary) ~ Group6, d9_)
    age <-         lm((salary) ~ Group6 * ages , d9_)
    term <- lm((salary) ~ (Term + Group6 + ages)^2 - Term:Group6 , d9_)
    stream <- lm((salary) ~ (Stream + Term + Group6 + ages)^2 - Term:Group6 , d9_)
    area <-        lm((salary) ~ (Stream + Term + Group6 + Area +ages)^2 - Term:Group6 -Area:Group6 , d9_)
    rank <-        lm((salary) ~ (Stream + Term + Group6 + Area + ages + Rank_)^2 - Term:Group6 - Area:Group6 - Group6:Rank_ 
                      - Stream:ages - Stream:Rank_ - ages:Rank_ , d9_)
    full <-  lm((salary) ~ (Stream + Term + Group6 + Area + ages + Rank_)^2 - Term:Group6 - Area:Group6 - Group6:Rank_ 
                - Stream:ages - Stream:Rank_ - ages:Rank_ , d9_)
  }) %>% rev -> fitlist_ns # no spline


for(nn in names(fitlist)){
  d9_[[paste0('gap___',nn)]] <- d9_$salary - exp(waldf(fitlist[[nn]], pred = within(d9_, Group6[] <- 'Male Comparator'))$coef/100)
}

names(d9_) %>% grepv('___',.)
d9l <- tolong(d9_, sep = '___', timevar = 'adjust')
with(d9l, tapply(gap, list(Group6, adjust), mean)) -> mat

mat <- mat - outer(rep(1,nrow(mat)),mat[1,] )
mat <- mat[, c('raw','age','term','stream','area','rank')]
dmat <- rbind(diag(6), cbind(-diag(5),0) + cbind(0,diag(5)))
rownames(dmat) <- c(names(fitlist), paste(names(fitlist)[-1], '-', names(fitlist)[-length(names(fitlist))]))

(mat %*% t(dmat) ->fmat)  %>% fmt(0) %>% print(quote = F) -> pmat

mat %>% as.table %>% as.data.frame %>% 
  rename(c("Var1",'Var2'),c('Group','Adjust')) %>% 
  within(
    {
      Adjustment <- factor(Adjust, levels = c('raw','age','term','stream','area','rank'))
      Adjustment <- tr(Adjustment,'raw','unadjusted')
    }
  ) -> dgaps
{
  tps(pch = 21:25, cex = 1, col = COLS, fill = paste0(COLS,'33'), alpha = 1)
  subset(dgaps, !grepl('Surveyed', Group)) %>% 
    droplevels %>% 
    within({
      Group <- reorder(Group, - Freq)
      Adjustment_ind <- as.numeric(Adjustment)
    }) -> z
  spint <- function(x) gsp(x, sort(rep((1:6),each =2) + c(-.1,.1)), c(0,0,3,0,3,0,3,0,3,0,3,0,0), 1)
  
  fintra <- lm(Freq ~ Group *spint(Adjustment_ind), z) 
  filldf <- expand.grid(Adjustment_ind = seq(.9,6.1,.01), Group = levels(z$Group)) 
  head(filldf)
  filldf$y <- predict(fintra, newdata = filldf)
  head(z)
  
    xyplot( Freq ~ Adjustment, z, group = Group, type = 'p', 
            scales = list(y = list(at = seq(-60000,10000, 10000), labels = fmt(seq(-60000,10000, 10000),0))),
            subset = !grepl('Surveyed', Group),
            ylab = 'Average adjusted salary difference',
            xlab = 'Cumulatively adjusted factors',
            auto.key = list(space='right')) +
   layer_(panel.abline(h = seq(-60000,10000,1000), col = 'gray95')) +
      xyplot(y ~ Adjustment_ind, filldf, group =  Group, type = 'l', lwd = 2, lty = 1, alpha = .4)
    
    xyplot( Freq ~ Adjustment, z, group = Group, type = 'b', 
            scales = list(y = list(at = seq(-60000,10000, 10000), labels = fmt(seq(-60000,10000, 10000),0))),
            subset = !grepl('Surveyed', Group),
            ylab = 'Average adjusted salary difference',
            xlab = 'Cumulatively adjusted factors',
            auto.key = list(space='right')) +
      layer_(panel.abline(h = seq(-60000,10000,1000), col = 'gray95')) 
    
}


#'
#' Add commentary on interpretation of Graph
#' Maybe show in different orders
#' maybe sepearate Term and Stream
#' - Note deviation fro Oaxaca-Blinder: interpretation of 'endowments' not appropriate here
#' - Much more a question of opportunities possibly unevenly distributed
#' - Also age and experience not equatable since age conflates age with cohorts. It reflects historical opportunities.
#' 
#'
#'
z <- lm(lnt(salary) ~ (Stream + Term + Group6 + Area + ages + Rank_)^2 - Term:Group6 - Area:Group6 - Group6:Rank_ + Stream*sq3(ages)
   - Stream:ages - Stream:Rank_ - ages:Rank_ , d9_)
summary(z)

#'
#' Builing a hypothesis matrix
#'
#' Given a sequence of nested models, $M_1, M_2,...,M_F$, and a comparator
#' level, $C$ for a variable, $G$, we want to estimate the mean gap between
#' levels of $G$ with the comparator level within each model.
#' 
#' The individual gaps for model $M_i$ can be obtained by subtituting level
#' $C$ for variable $G$
#' 
library(spida2)
library(latticeExtra)
set.seed(123)
expand.grid(age= seq(30,70), g = c('A','B','C'), area = c('a','b','c')) %>% 
  within({
    p <- 1/(1 + exp(-scale(age*(as.numeric(g)-2)*(as.numeric(area)-2))))
    n <- rpois(length(p), abs(p+3))
    y <- age + as.numeric(g) + as.numeric(area) + rnorm(n)
  }) -> d1
dd <- with(d1, d1[rep(1:nrow(d1), n),])
dd$id <- 1:nrow(dd)
dim(dd)
head(dd)

list() %>% 
  within(
    {
  fit1  <- lm(y ~ g, dd)
  fit2 <- lm (y ~ g + age, dd)
  fit3 <- lm (y ~ g + age + area, dd)
  
    }
  ) %>% rev -> fitl


decomp <- function(fitl, g, comp, droplast = NULL) {
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
  lssvd <- function(x, y, zero = 10^(-14)) {
    xp <- svd(x,nu=ncol(x), nv = ncol(x))
    uy <- t(xp$u)%*%cbind(y)
    dinv <- 1/xp$d
    dinv[abs(xp$d) < zero] <- 0
    coef <- t(t(xp$v)*dinv) %*% uy
    resid <- y - x %*% coef
    list(coef = coef, resid = resid, sse = sum(resid^2))
  }
  along <- function(x) {
    ret <- seq_along(x)
    names(ret) <- names(x)
    ret
  }
  
  #    
  library(car)
  ret <- list()
  ret[['names']] <- list(gname = g, gcomplevel = comp)
  
  #
  # Refit with data set used for full model to check for common data set
  #
  
  full <- fitl[[length(fitl)]]
  d <- getD(full)
  fit_up <- lapply(fitl, update, data = d)
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
    seq_along(mfs), 
    function(ii) {
      sum(abs(lsfit(mfs[[ii]], mffull, intercept =FALSE )$residuals))
    }
  )
  resids2 <- lapply(
    seq_along(mfs), 
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
  
  ret[['pred']] <- pred
  
  #
  # gaps calculated from individual models 
  #
  
  # from above: mfs <- lapply(fitl, model.matrix, data = d)
  # from above: mgs <- lapply(fitl, model.matrix, data = dg)
  # not needed: mffull <- mfs[[length(mfs)]]
  # not needed: Bs <- lapply(seq_along(mfs), function(ii) lsfit(mfs[[ii]], mffull, intercept =FALSE )$coef)
  Ls_each <- lapply(
    seq_along(mfs), 
    function(ii) {
      Lmat <- groupL %*% (mfs[[ii]] - mgs[[ii]])    # not needed:  %*% Bs[[ii]]
      attr(Lmat, 'data') <- subset(data, model == names(fitl)[ii]) 
      waldf(fitl[[ii]], Lmat)
    })
  ret[['gaps_each']] <-Ls_each %>% lapply(subset, select = -L) %>% do.call(rbind,.)
  
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
  
  class(ret) <- 'decomp'
  ret
} 

decomp(fitl, 'g', 'C') -> z
z$names
z$gaps %>% class
decomp(fitlist, 'Group6','Male Comparator') -> zz

names(zz)
zz$pred %>% tab(coef ~ model + Group6)  # very good








decomp(fitlist_ns, 'Group6','Male Comparator') -> zz_ns
zz$diags
zz$resids
zz$resids2

merge(zz$gaps,zz$gaps_each, by = c('Group6','model'), all = T) %>% 
  sortdf(~model/Group6)-> gapsboth
gapsboth[,c('Group6','model','coef.x','coef.y','se.x','se.y')]
zz$gapdiffs
zz$dout %>% lapply(class)

xyplot(gresids ~ fact)

{
tps(grid.pars = list(fontfamily = 'Serif'))
gapplot <- function(z) {
  disp(names(z$gaps))
  data <- z$gaps
  xyplot(coef ~ model, 
         z$gaps, 
         #subset = subset,
         groups = data[['Group6']], type = 'b',
        scales = list(y = list(at = seq(-60000,10000, 10000), labels = fmt(seq(-60000,10000, 10000),0))),
        ylab = 'Average adjusted salary difference',
        xlab = "Cumulatively adjusted factors",
         auto.key = list(space = 'right'))
  
}
gapplot(z)
zz$gaps %>% head
gapplot(zz)
}
{
  tps(grid.pars = list(fontfamily = 'Serif'))
  resplot <- function(zz, subset = TRUE) {
    xyplot(gresids ~ model | g, 
           z$dout, 
           subset = subset,
           groups = g, type = 'b',
           scales = list(y = list(at = seq(-60000,10000, 10000), labels = fmt(seq(-60000,10000, 10000),0))),
           ylab = 'Average adjusted salary difference',
           xlab = "Cumulatively adjusted factors",
           auto.key = list(space = 'right'))
    
  }
  gapplot(z)
}

{
  tps(grid.pars = list(fontfamily = 'Serif'))
  resplot <- function(z, subset = TRUE) {
    
    xyplot(gresids ~ factor(model) | g, groups =id.,
           z$dout, 
           subset = subset,
           type = 'b',
           scales = list(y = list(at = seq(-60000,10000, 10000), labels = fmt(seq(-60000,10000, 10000),0))),
           ylab = 'adjusted salary differences',
           xlab = "Cumulatively adjusted factors",
           auto.key = list(space = 'right'))
    
  }
  resplot(z)
}


zz <- decomp(fitlist, 'Group6','Male Comparator')


xyplot( Freq ~ Adjustment, zz, groups = id., type = 'b', 
        scales = list(y = list(at = seq(-60000,10000, 10000), labels = fmt(seq(-60000,10000, 10000),0))),
        subset = !grepl('Surveyed', Group),
        ylab = 'Average adjusted salary difference',
        xlab = 'Cumulatively adjusted factors',
        auto.key = list(space='right')) +
  layer_(panel.abline(h = seq(-60000,10000,1000), col = 'gray95')) 

names(z)
xyplot(y ~ age | area, z, groups = g, auto.key = T)

## lsfit with models not of full rank

xz <- cbind(1, 1:10, 1:10 -5, rnorm(10))
yz <- 1:10 +5
lsfit(xz, yz, intercept = FALSE)
n <- as.numeric
expand.grid(a = letters[1:2],b = letters[1:2])[rep(1:4,20),] %>% 
  within(
    {
      y <- n(a) * n(b) + rnorm(a)
      y[a=='a' & b == 'a'] <- NA
    }
  ) -> zd
fit <- lm(y ~ a*b, zd, na.action = na.exclude)
summary(fit)
zd$pred <- predict(fit)
zd$predn <- predict(fit, newdata = zd)
zd %>% with(sum(abs(pred - predn), na.rm = T))
zd$predw <- waldf(fit, pred = zd)$coef
zd %>% with(sum(abs(pred - predw), na.rm = T))
