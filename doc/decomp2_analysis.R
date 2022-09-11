#' ---
#' author: 'decomp2 analysis'
#' date: "2022-09-10"
#' output: bookdown::pdf_document2
#' ---

#' ```{r setup, include=FALSE}
#' knitr::opts_chunk$set(echo = TRUE)
#' ```
#'
#' # Components of decomp2
#'
#' How this works:
#'
#' 1. refit models with groups set to comparator
#'   - issues:
#'     - Safe (STABLE) prediction: If Group is involved in (interacts with) unstable prediction
#'       then there can be error but if other variables are unstable it's
#'       okay because the full model is the same with respect to those variables
#'   - The result in 'dout' and resids from comparator model are in **dout$gresids**
#'   - We will add the differences between predicted values as gaps
#'
#'   Note: difference between decomp and decomp2
#'   decomp2$groupL is a data frame with group and cond variables as well as
#'   a matrix decomp2$groupL$groupL
#'
#' 2. The group incidence matrix is combined with submodel to full model
#'    Wald L matrices to form hypothesis matrices for the submodel gaps
#'    and for the disparity reductions between models.
#'
#'
library(spida2)
library(peq)
library(latticeExtra)
mtcars %>%
  within(
    {
      Cyl <- factor(cyl)
      Carb <- factor(carb)
    }
  ) -> z
list() %>%
  within(
    {

    fit1 <- lm(mpg ~ Cyl, z)
    fit2 <- lm(mpg ~ Cyl + gear, z)
    fit3 <- lm(mpg ~ Cyl * gear, z)
    fit4 <- lm(mpg ~ Cyl * gear * Carb, z)
    }
  ) %>% rev -> fitlist

library(peq)
fitlist %>% icp

# undebug(decomp2)
zzc <- decomp2(fitlist, "Cyl", "4", z, c("gear","Carb"))
zzm <- decomp(fitlist, "Cyl", '4', z)

zzm$groupL %>% apply(1, sum)
zzc$groupL$groupL %>% na20 %>% apply(1,sum)

zzm$names
zzc$names

names(zzm)

gapplot(zzm)
gapplot(zzc)

zzm$dout %>% all.equal(zzc$dout)

zzc$allgaps %>% dim
zzc$allgaps$coef

zz <- decomp(fitlist, "Cyl", "4", z)
zz %>% gapplot


#' This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

#' When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

#' ```{r cars}
#' summary(cars)
#' ```

#' ## Including Plot

#' You can also embed plots, for example:

#' ```{r pressure, echo=FALSE}
#' plot(pressure)
#' ```

#' Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
