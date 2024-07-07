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
contr.rot(3) %>% crossprod
?C


contrasts(simdata$gender)

(fit0 <- lm(sal ~ gender, simdata)) %>% summary
(fit1 <- lm(sal ~ gender + age, simdata)) %>% summary
(fit2 <- lm(sal ~ gender + age + area, simdata)) %>% summary
(fit3 <- lm(sal ~ area/gender + age -1 , simdata)) %>% summary
(fit4 <- lm(sal ~ area/gender + age + rank -1 , simdata)) %>% summary

fitl <- list(gender = fit0, age = fit1, area = fit2, g_by_a = fit3, rank = fit4)

pred <- with(simdata, spida2::pred.grid(gender,age,area))
simdata |> Filter(is.factor, x = _) |> lapply(contrasts)
pred |> Filter(is.factor, x = _) |> lapply(contrasts)

debug(model.frame.default)

pred$fit3 <- predict(fit3, newdata = pred)
pred |> Filter(is.factor, x = _) |> lapply(contrasts)


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

waldf(fit3, W3 %*% L3)
waldf(fit3, Wsize %*% L3)
waldf(fit3, Weqn %*% L3)

decomp(fitl, 'gender', 'M') %>% gapplot
decomp2(fitl, 'gender', 'M', cond = 'area') %>% gapplot

decomp(fitl, 'gender', 'M') %>% decomp_table(log=FALSE)
decomp(fitl, 'gender', 'M') %>% decomp_table(log=TRUE)
decomp2(fitl, 'gender', 'M', cond = 'area') %>% gapplot
decomp2(fitl, 'gender', 'M', cond = 'rank') %>% gapplot

fitl |> lapply(getX) |> lapply( function(m) {
  c(Matrix::rankMatrix(m), ncol(m))
})

# TODO: Use SEs for fit3 to show how precision weighted means are same as
# gaps obtained from non-interaction model

decomp2(fitl, 'gender','M', cond = c('area','rank')) %>% gapplot

#
#
# Simple example with age and gender show
# possibility of a reversal depending
# on the choice of comparator group
#
dage <- read.table(header = TRUE, stringsAsFactors = TRUE,text ="
Gender    Age   Salary
M         30    30000
M         50    50000
M         60    60000
M         65    65000
M         65    65000
M         70    69000
M         70    70000
M         70    71000
F         30    30000
F         40    35000
F         40    36000
F         40    34000
F         50    40000
")


tps(pch=16:17)
xyplot(Salary ~ Age, dage, groups = Gender, alpha = .5)
fit1 <- lm(Salary ~ Gender, dage)
fit2 <- lm(Salary ~ Gender + Age, dage)
fit3 <- lm(Salary ~ Gender * Age, dage)
fitl <- list(Gender=fit1, 'Gender + Age'=fit2, 'Gender * Age'= fit3)
icp(fitl)
decomp(fitl, 'Gender', 'M')
decomp(fitl, 'Gender', 'M') %>% decomp_table
decomp(fitl, 'Gender', 'F') %>% gapplot
decomp(fitl, 'Gender', 'M') %>% gapplot


#
# The following is example of a reversal paradox in which
# the direction of the gap is reversed depending on the
# choice of reference groups. This can only occur with
# interaction with Gender in the model and
#


dage2 <- within(
  dage,
  {
    Salary2 <- ifelse(Gender == 'F', Salary + 10000, Salary)
  }
)
fit21 <- lm(Salary2 ~ Gender, dage2)
fit22 <- lm(Salary2 ~ Gender + Age, dage2)
fit23 <- lm(Salary2 ~ Gender * Age, dage2)
fit2l <- list(Gender=fit21, 'Gender + Age'=fit22, 'Gender * Age'= fit23)
decomp(fit2l, 'Gender', 'F') %>% decomp_table
decomp(fit2l, 'Gender', 'M') %>% decomp_table
decomp(fit2l, 'Gender', 'F') %>% gapplot
decomp(fit2l, 'Gender', 'M') %>% gapplot

