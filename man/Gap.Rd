% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decompositions.R
\name{Gap}
\alias{Gap}
\alias{plot.Gap}
\alias{summary.Gap}
\title{Create a Gap object for a single model to plot or make tables}
\usage{
Gap(formula, data, pred, groupvar = "Gender", complevel = "M", ...)

\method{plot}{Gap}(x, fmla, ylim = NULL, ylab = NULL, f = 0.1, type = 1, ...)

\method{summary}{Gap}(object, formula, type = 1, ...)
}
\arguments{
\item{formula}{model formula for fitting with \code{\link{lm}}
or fitted model. The
**grouping variable** should be first so that \code{\link{subrow}}
will differentiate correctly. **Maybe only pred matters**}

\item{data}{data set}

\item{pred}{prediction data set. It should include all variables in the
model formula. The grouping variable should be first so that
\code{\link{subrow}} will differentiate correctly.}

\item{groupvar}{name of grouping variable. Default: Gender}

\item{complevel}{level of groupvar that is used as comparator. Default: 'M'}

\item{...}{other arguments}

\item{x}{Gap object}

\item{fmla}{formula for plotting}

\item{ylim}{limits for y axis. Default NULL to get calculation in code}

\item{ylab}{label for vertical axis. Default NULL}

\item{f}{factor to expand default 'ylim'. Default: .1}

\item{object}{a Gap object created by \code{\link{Gap}}}
}
\value{
A list containing the following components
  - fit: the fitted model
  - data augmented with components
    - Fit: Fitted (predicted) value for individual
    - Line: Fitted value for individual if they belonged to comparator group
    - Gap: Line - Fit
    - Residual: Line - Salary
    - residual: Fit - Salary
  - pred augmented with components
    - Fit: Fitted (predicted) value for individual
    - Line: Fitted value for individual if they belonged to comparator group
    - Gap: Line - Fit
}
\description{
**CHECK THIS:** Fits a model with lm and plots **weighted adjusted gaps** or
simply **adjusted gaps** depending on the inclusion/exclusion of
variables in the plotting formula relative to the fitting formula.
}
\section{Methods (by generic)}{
\itemize{
\item \code{plot(Gap)}: plot method

\item \code{summary(Gap)}: method for printing kable output

}}
\examples{
\dontrun{

# data and initial part copied from 'decomp'

library(peq)
library(spida2)
library(latticeExtra)
library(kableExtra)
simdata <- read.table(header = TRUE, stringsAsFactors = TRUE,text ="
sal  gender  area   age  rank
100  F       A      30   a   # group A: gap: -10 and -20 for A
110  F       A      40   b
120  F       A      50   b
120  M       A      40   a
130  M       A      50   b
140  M       A      60   b
100  A       A      40   a
110  F       B      30   a # group B: gap: -20 and -30 for A
130  M       B      30   a
140  M       B      40   a
140  M       B      40   b
149  M       B      50   b
150  M       B      50   c
151  M       B      50   c
150  M       B      50   b
150  M       B      50   b
150  M       B      50   c
150  M       B      50   c
150  M       B      50   c
150  M       B      50   b
120  A       B      50   b
")
simdata$gender <- relevel(simdata$gender, 'M')
tps(pch=16)
xyplot(sal ~ age | area, simdata, groups = gender, auto.key = TRUE)


contr.helmert(3)
contr.rot <- function(n) {
   ret <- contr.helmert(n)
   disp(ret)
   ret <- t( t(ret)/sqrt(apply(ret^2, 2, sum)))
   ret
}
contr.rot(3)
contr.rot(3) \%>\% crossprod
contrasts(simdata$gender)

(fit0 <- lm(sal ~ gender, simdata)) \%>\% summary
(fit1 <- lm(sal ~ gender + age, simdata)) \%>\% summary
(fit2 <- lm(sal ~ gender + age + area, simdata)) \%>\% summary
(fit3 <- lm(sal ~ area/gender + age -1 , simdata)) \%>\% summary
(fit4 <- lm(sal ~ area/gender + age + rank -1 , simdata)) \%>\% summary

fitl <- list(gender = fit0, age = fit1, area = fit2, g_by_a = fit3, rank = fit4)

pred <- with(simdata, spida2::pred.grid(gender,age,area))
pred$fit3 <- predict(fit3, newdata = pred)

xyplot(fit3 ~ age | area, pred, groups = gender, type = 'l',auto.key = TRUE) +
xyplot(jitter(sal,10) ~ age | area, simdata, groups = gender, auto.key = TRUE)

# Using fit1: no effect of area
L1 <- rbind(Agap = c(0,1,0,0),
            Fgap = c(0,0,1,0))
waldf(fit1, L1)

# Using fit2: additive area
L2 <- rbind(Agap = c(0,1,0,0,0),
            Fgap = c(0,0,1,0,0))
waldf(fit2, L2)

# Using fit3: no effect of area
L3 <- rbind("Agap in A" = c(0,0,0,1,0,0,0),
            "Fgap in A" = c(0,0,0,0,0,1,0),
            "Agap in B" = c(0,0,0,0,1,0,0),
            "Fgap in B" = c(0,0,0,0,0,0,1))
waldf(fit3, L3)

SEs <-   waldf(fit3, L3)$se
wts <- 1/SEs^2   # precision weights

# Type III average over Areas

W3 <- rbind( A3 = std(c(1, 0 , 1, 0)),
             F3 = std(c(0, 1 , 0, 1)))
W3

# Average by size of Areas  (dubious)

nArea <- tab__(simdata, ~ area)
nArea <- nArea[c(1,1,2,2)]

Wsize <- rbind(
    AnA = std(nArea * c(1,0,1,0)),
    FnA = std(nArea * c(0,1,0,1))
)
Wsize

# Average by size of equity group

## Compare relative weights of different methods
## Idea: Using additive model may give reasonable weights if no group is very small
## but bad idea for very small groups were the gap might not reflect the
## average **experience** of members of the small group.
## Weights will tend to look like (1/nF + 1/nM)^1 which will look like nF if small ..... so ....
## ... but maybe not so simple with other adjustment variables like age, etc....

ns <- tab__(simdata, ~ gender + area)
ns
Weqn <- rbind(
    AnA = std(c(ns['A','A'],0,ns['A','B'],0)),
    FnA = std(c(0,ns['F','A'],0,ns['F','B']))
)
Weqn

waldf(fit3, W3 \%*\% L3)
waldf(fit3, Wsize \%*\% L3)
waldf(fit3, Weqn \%*\% L3)

decomp(fitl, 'gender', 'M') \%>\% gapplot
decomp2(fitl, 'gender', 'M', cond = 'area') \%>\% gapplot

decomp(fitl, 'gender', 'M') \%>\% decomp_table(log=FALSE) |> kb()
decomp(fitl, 'gender', 'M') \%>\% decomp_table(log=TRUE)|> kb()
decomp2(fitl, 'gender', 'M', cond = 'area') \%>\% gapplot
decomp2(fitl, 'gender', 'M', cond = 'rank') \%>\% gapplot

## Using Gap ####

library(spida2)
library(peq)
pred <- with(simdata, pred.grid(gender, area, age = seq(30,60,1) ))

gap <- Gap(sal ~ gender + age + area, simdata, pred , 'gender')
names(gap)
gap$pred |> xyplot(Fit~ age | area, data = _, group=gender, auto.key = T)
gap$pred \%>\% {
xyplot(Gap~ age | area, data = ., group=gender, auto.key = T,
   fit = .$Fit,
   upper = .$upper,
   lower = .$lower,
   subscripts = TRUE) +
glayer(panel.fit(...))
}
}
}
