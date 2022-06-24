## 
## 
## From R/gnew/R/utils.R on 2022_06_23
## 
## 
## small utilities for bad typists ------
##
#' getAnywhere
#' 
#' @param ... parameters for getAnywhere
#' 
#' @export
ga <- function(...) {
  getAnywhere(...)
}
##
## Kable utilities ------
##
#' Collect options for kbl
#'
#' @export
kb <- function(
  x, 
  caption = NULL,
  format.args = list(big.mark = ','), 
  ...) {
  kbl(x, 
      booktabs = T, 
      linesep = '', 
      format.args = format.args,
      caption= caption, ...
  ) %>% 
    kable_styling(latex_options = "HOLD_position")
}
#' Turn rownames into labelled column for kable
#'   
#' @export
add_rownames <- function(x, label = names(dimnames(x))[1], ...) {
  x <- cbind(x) # to make it a matrix if it isn't
  vals <- cbind(rownames(x))
  x <- cbind(vals, x)
  dimnames(x)[[2]][1] <- label
  rownames(x) <- NULL
  x
}
#' Generic round
#' 
#' @export
rnd <- function(x,...) {
  if(is.numeric(x))  round(x, ...)
  else if(is.list(x)) {
    x[] <- lapply(x, rnd, ...)
    x
  } else x
}
#' Formatting that thinks it's smart
#' 
#' Chooses number of significant digits semi-intelligently
#' 
#' @param x numeric object to be formatted
#' @param nsmall, nsig, nright, 
#' @param scientific FALSE by default
#' @param big.mark is ',' by default
#' @param ... other arguments to \code{\link{format}}
#' 
#' @export
fmt <- function(x, nsmall = ns_(x, nsig, nright), nsig = 4, nright = 1,
                scientific = FALSE, big.mark = ',',...) {
  if(!is.numeric(x)) return(x)
  if(missing(nsmall) && missing(nright) && all(floor(x) == x)) nsmall <- 0
  disp <- function(x) NULL
  ns_ <- function(x, nsig, nright) {
    ok <- is.finite(x) & (na20(x) != 0)
    disp(ok)
    if(sum(ok) ==0) return(0)
    max(c(-log10(abs(x[ok]))+nsig, nright), na.rm = T)
  }
  x[] <- round(x, nsmall)
  disp(nsmall)
  x[] <-format(x, nsmall = nsmall, big.mark = big.mark,scientific = scientific, ...)
  x
}
##
## Factor utilities ------
##

nuniq <- function(x) {
  length(unique(x))
}

na2na <- function(x, replace , ...) UseMethod('na2na')
na2na.default <- function(x, replace= "No Answer", ... ) {
  x[is.na(x)] <- replace
  x
}
na2na.factor <- function(x, replace= "No Answer", ...) {
  levs <- levels(x)
  x <- as.character(x)
  x <- na2na(x, replace = replace)
  factor(x, levels = unique(c(levs, replace)))
}
# zf <- factor(c(NA,LETTERS[1:5],NA))
# zf
# zf %>% na2na
# zf %>% na2na('not surveyed')
#' @export
relevel_last <- function(x, ref, ...) {
  x <- as.factor(x)
  rotate <- function(x) c(x[-1],x[1])
  if(!(ref %in% levels(x))) return( x )
  ret <- relevel(x, ref, ...) 
  factor(ret, levels = rotate(levels(ret)))
}

#' @export
relevel_first <- function(x, ref, ...) {
  x <- as.factor(x)
  if(!(ref %in% levels(x))) return( x )
  relevel(x, ref, ...) 
}


#' @export
relevel_size <- function(x, ord = -capply(x, x, length), ...) {
  reorder(as.factor(x), ord)
}

#' @export
combine_small_levels <- function(x, size = 6, name = 'Smaller groups') {
  tofac <- is.factor(x)
  if(tofac) lev <- levels(x)
  x <- as.character(x)
  nx <- capply(x, x, length)
  x[nx <= size] <- name
  if(tofac) {
    new_levs <- intersect(lev, unique(x))
    new_levs <- c(new_levs, name)
    x <- factor(x, levels = new_levs)
  }
  x
}

# in spida2: years <- function(x,...) as.numeric(format(x,"%Y"))

##
## Lattice utilities ------
##
#' Auto.key arguments
#' 
#' @export
ak <- function(reverse.rows = T, lines = T, points = T, ...) {
  list(space = 'right', reverse.rows = reverse.rows, lines = lines,
       points = points)
}
#' Set latice parameters for base (not superpose) elements
#' 
#' @export
td_ <- function(...) td(..., superpose = FALSE)
#'
#' Presumably to create a heading
#' 
#' @export
main <- function(x, line2 = NULL, font = 1, cex = 1, ...) {
  if(!is.null(line2)) x <- paste0(x,'\n',line2)
  list(x, font = font, cex = cex, ...)
}
#' tr_gs: special case
#' 
#' @export
tr_gs <- function(x) {
  x <- tr(x, c("M:F1", "M:FA", "F:F1", "F:FA"), c('Male F1', 'Male FA', 'Female F1', 'Female FA') )
  factor(x, levels =  c('Male F1', 'Female F1','Male FA',  'Female FA'))
}
#' add number of occurences in parentheses
#' 
#' @export
addn <- function(x, big.mark = ',') {
  # add frequency in parentheses
  library(knitr)
  library(kableExtra)
  xx <- as.character(x)
  n <- format(capply(xx, xx, length), big.mark= ',', justify = 'none')
  n <- paste0('(',gsub(' ','',n),')')
  ret <- apply(data.frame(xx, n),1, paste, collapse = ' ')
  ret <- reorder(factor(ret), as.numeric(as.factor(x)))
  levels(ret) <- text_spec(levels(ret), format='latex', monospace = F)
  ret
}
# addn(c('a','a','a','b'))
##
## Miscellaneous ------
##
#' Quick grepv
#' 
#' 
#' @export
g <- function(string, data = dd) {
  grepv(string, names(data))
}
#' 
#' 
#' @export
reorder_last <- function(x, last) {
  wh <- which(x %in% last)
  if(length(wh) == 0) x
  else c(x[-wh], x[wh])
}
#' relevel for character vectors
#' 
#' @export
reorder_first <- function(x, first) {
  wh <- which(x %in% first)
  if(length(wh) == 0) x
  else c(x[wh], x[-wh])
}



