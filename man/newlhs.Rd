% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decompositions.R
\name{newlhs}
\alias{newlhs}
\alias{newdep}
\title{Change the dependent variable of a formula}
\usage{
newlhs(fmla, depvar = "")

newdep(fmla, depvar = "")
}
\arguments{
\item{fmla}{a formula}

\item{depvar}{a character string with the intended name of the LHS of the
formula, default: ''}
}
\value{
a formula with the same RHS as 'fmla' and the value of depvar
         for its LHS
}
\description{
Change the dependent variable of a formula
}
\section{Functions}{
\itemize{
\item \code{newdep()}: deprecated name

}}
