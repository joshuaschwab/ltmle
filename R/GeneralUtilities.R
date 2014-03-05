#General utilities

is.equal <- function(...) {
  isTRUE(all.equal(...))
}

#same as model.matrix, but leave NA rows
model.matrix.NA <- function(object, data) {
  MM <- model.matrix(object, data)
  MM <- MM[match(rownames(data), rownames(MM)), , drop=FALSE]
  rownames(MM) <- rownames(data)
  return(MM)
}

#if warning is in ignoreWarningList, ignore it; otherwise post it as usual
SuppressGivenWarnings <- function(expr, warningsToIgnore) {
  h <- function (w) {
    if (w$message %in% warningsToIgnore) invokeRestart( "muffleWarning" )
  }
  withCallingHandlers(expr, warning = h )
}

get.stack <- function(x, ifnotfound=NULL, mode="any") {
  #Look for object x (string) searching through the calling stack, starting at the current parent frame and going back through the parents
  if (!is.character(x)) stop("x must be a character string")
  for (f in (sys.parent(1)):0) {
    if (exists(x, envir=sys.frame(f), mode=mode, inherits=FALSE)) {
      return(get(x, envir=sys.frame(f), mode=mode, inherits=FALSE))
    }
  }
  return(ifnotfound)
}


rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))

# Given row and column numbers of a matrix with num.rows rows, compute the single index
sub2ind <- function(row, col, num.rows) {
  return((col - 1) * num.rows + row)
}

repmat <- function(X,m,n){
  #R equivalent of repmat (matlab)
  mx <- dim(X)[1]
  nx <- dim(X)[2]
  if ((m == 0) || (n == 0)) return(matrix(numeric(0), nrow=mx*m, ncol=nx*n)) #avoids warnings when m or n is 0
  return(matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T))
}

drop3 <- function(x) {
  #if x is an array with 3 dimensions and third dimension has one level, return a matrix with it dropped; otherwise error
  stopifnot(length(dim(x))==3)
  stopifnot(dim(x)[3]==1)
  dn <- dimnames(x)
  dim(x) <- dim(x)[1:2]
  dimnames(x) <- dn[1:2]
  return(x)
}
