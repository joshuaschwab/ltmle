#' @export
MaintainTreatment <- function(data, current.node, nodes) {
  #if the previous Anode is 1, all subsequent Anodes are 1 
  Anodes <- nodes$A
  if (!(current.node %in% Anodes)) return(NULL)
  if (!(any(Anodes < current.node))) return(NULL)
  
  prev.a.node <- max(Anodes[Anodes < current.node])
  is.deterministic <- data[, prev.a.node] == 1
  return(list(is.deterministic=is.deterministic, prob1=1))  
}

#' @export
MaintainControl <- function(data, current.node, nodes) {
  #if the previous Anode is 0, all subsequent Anodes are 0
  Anodes <- nodes$A
  if (!(current.node %in% Anodes)) return(NULL)
  if (!(any(Anodes < current.node))) return(NULL)
  
  prev.a.node <- max(Anodes[Anodes < current.node])
  is.deterministic <- data[, prev.a.node] == 0
  return(list(is.deterministic=is.deterministic, prob1=0))  
}


#' @export
deterministic.g.function_template <- function(data, current.node, nodes) {
  # data: the 'data' data.frame passed to ltmle/ltmleMSM
  # current.node: the column index of data corresponding to the A or C node (see is.deterministic below)
  # nodes: list of column indicies, components: A, C, L, Y, AC (Anodes and Cnodes combined and sorted), 
  #   LY (Lnodes and Ynodes combined, sorted, "blocks" removed - see ?ltmle)
  # Note that nodes may be passed to ltmle as either the names of nodes or numerical column indicies, but they
  #   are all converted to numerical indicies before deterministic.g.function is called
  
  # deterministic.g.function will be called at all Anodes and Cnodes
  # return(NULL) is equivalent to return(list(is.deterministic=rep(FALSE, nrow(data)), prob1=numeric(0)))
  
  #define is.deterministic here: vector of logicals, length=nrow(data)
  #define prob1 here: the probability that data[is.deterministic, current.node] == 1, 
  #  vector of length 1 or length(which(is.deterministic))
  is.deterministic <- stop("replace me!")
  prob1 <- stop("replace me!")
  return(list(is.deterministic=is.deterministic, prob1=prob1))  
}

#' @export
deterministic.Q.function_template <- function(data, current.node, nodes, called.from.estimate.g) {
  # data: the 'data' data.frame passed to ltmle/ltmleMSM
  # current.node: the column index of data corresponding to the A or C node (see is.deterministic below)
  # nodes: list of column indicies, components: A, C, L, Y, AC (Anodes and Cnodes combined and sorted), 
  #   LY (Lnodes and Ynodes combined, sorted, "blocks" removed - see ?ltmle)
  # called.from.estimate.g: TRUE or FALSE - your function will be called with called.from.estimate.g=TRUE during 
  #   estimation of g and called.from.estimate.g=FALSE during estimation of Q. During estimation of g, only
  #   the is.deterministic element of the return list will be used.
  # Note that nodes may be passed to ltmle as either the names of nodes or numerical column indicies, but they
  #   are all converted to numerical indicies before deterministic.Q.function is called
  
  # It is not necessary to specify that deterministic Y events (Y==1) indicate 
  #   a deterministic Q value of 1; this is automatic.
  # deterministic.Q.function will be called at all Lnodes and Ynodes (after removing "blocks") 
  #   and Anodes and Cnodes (see called.from.estimate.g above)
  # return(NULL) is equivalent to return(list(is.deterministic=rep(FALSE, nrow(data)), Q.value=numeric(0)))
  
  #define is.deterministic here: vector of logicals, length=nrow(data)
  #define Q.value here: the iterated expectation of the final Y, vector of length 1 or length(which(is.deterministic))
  is.deterministic <- stop("replace me!")
  Q.value <- stop("replace me!")
  return(list(is.deterministic=is.deterministic, Q.value=Q.value))  
}
