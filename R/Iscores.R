#' Iscores: compute the imputation KL-based scores
#'
#' @param imputations a list of list of imputations matrices containing no missing values of the same size as X.NA
#' @param methods a vector of characters indicating which methods are considered for imputations. It should have the same length as the list imputations.
#' @param X.NA a matrix containing missing values, the data to impute.
#' @param m the number of multiple imputation to consider, defaulting to the number of provided multiple imputations.
#' @param num.proj an integer specifying the number of projections to consider for the score.
#' @param num.trees.per.proj an integer, the number of trees per projection.
#' @param min.node.size the minimum number of nodes in a tree.
#' @param sample.splitting a boolean indicating if sample splitting should be used.
#' @param frac.test a numeric between 0 and 1, the fraction of test points.
#' @param n.cores an integer, the number of cores to use.
#' @param ... additional parameters.
#'
#' @return The scores for each imputation method.
#'
#' @export
Iscores <-function(imputations,
                   methods,
                   X.NA,
                   m = length(imputations[[1]]),
                   num.proj = 100,
                   num.trees.per.proj = 5,
                   min.node.size = 5,
                   sample.splitting = TRUE,
                   frac.test = 0.5,
                   n.cores = 1,
                   ...) {


  if (!is.null(names(imputations))) {
    if (!identical(names(imputations),methods)) {
      stop("imcompatible names between imputations and methods.")
    }
    if (length(imputations)!=length(methods)) {
      stop("different lengths of imputations and methods.")
    }
  } else {
    names(imputations) <- methods
  }

  ## original data
  # candidate missing value points
  ind.candidates <- which(!complete.cases(X.NA))
  ind.candidates <- sort(ind.candidates)
  nrofmissing <- sum(is.na(X.NA))

  # save the metrics
  scores.all.dr.kl<- list()


  # Initialize
  for (method in methods){
    scores.all.dr.kl[[method]] <- rep(NA, m)
  }


  # get the NA patterns
  NA.pat <- X.NA
  NA.pat[!is.na(NA.pat)] <- 1
  NA.pat.unique <- unique(NA.pat)
  NA.pat.groups <- apply(NA.pat.unique,
                         1,
                         function(p) which(apply(NA.pat, 1, function(pp) identical(pp,p))))

  # remove the fully observed patterns from the data
  V<-apply(NA.pat.unique,
           1, sum, na.rm=T)
  NA.pat.unique <- NA.pat.unique[which(V!=ncol(X.NA)),,drop=F]
  NA.pat.groups <- NA.pat.groups[which(V!=ncol(X.NA))]


  for (method in methods){

    scores.all <- list()

    # switch cases
    ## Continue here with checking!!
    for (part in 1:2) {

      # apply on the unique patterns
      scores.all.part <- mclapply(1:nrow(NA.pat.unique), function(j){

        # if sample splitting is chosen, then we separe in two parts
        if (sample.splitting) {
          # in case only one point is there we do not separate
          if (length(NA.pat.groups[[j]])==1) {
            ids.pattern.test <- NA.pat.groups[[j]]
            ids.pattern.train <- 1:nrow(X.NA)
          } else {
            if (part == 1) {
              ids.pattern.test <- NA.pat.groups[[j]][1:(floor(length(NA.pat.groups[[j]])*frac.test))]
            } else {
              ids.pattern.test <- NA.pat.groups[[j]][-c(1:(floor(length(NA.pat.groups[[j]])*frac.test)))]
            }
            ids.pattern.train <- setdiff(1:nrow(X.NA), ids.pattern.test)
          }
        } else {
          ids.pattern.test <- NA.pat.groups[[j]]
          ids.pattern.train <- 1:nrow(X.NA)
        }

        scoredr <- list()

        scores <- lapply(1:m, function(set){


          X.h <-imputations[[method]][[set]] #imputations[[method]][[sample((1:m)[-set],1)]]
          X.h <- as.matrix(X.h)[ids.pattern.train,,drop=F]
          #set.seed(1)
          # Delete Comment: This is the only place where we would get
          # different results from the old code
          object.dr <- densityRatioScore(X = X.NA[ids.pattern.train,,drop=F],
                                         Xhat = X.h,
                                         x =  NA.pat.unique[j,],
                                         X.true = X[ids.pattern.train,,drop=F],
                                         compute.glm = TRUE,
                                         num.proj=num.proj,
                                         num.trees.per.proj = num.trees.per.proj,
                                         min.node.size = min.node.size)


          ## Define the test set!
          Z <-  as.matrix(imputations[[method]][[set]])[ids.pattern.test,, drop=F]
          Z <- unname(Z)
          Z <- as.matrix(Z)



          scoredr.kl <- sum(compute_drScore(object = object.dr, Z = Z, type = "kl.score"))



          return(list(scoredr.kl = scoredr.kl))

        }
        )

        scoredr.kl <- unlist(lapply(scores, function(l) l$scoredr.kl))


        return(list(scoredr.kl =scoredr.kl))
      }, mc.cores = n.cores)

      scores.all[[part]] <- scores.all.part
    }


    dat.scoredr.kl <- do.call(rbind, lapply(unlist(scores.all, recursive = FALSE), function(l) l$scoredr.kl))
    dat.scoredr.kl <- colSums(dat.scoredr.kl) / sum(lengths(NA.pat.groups))

    scores.all.dr.kl[[method]] <- dat.scoredr.kl

  }

  iscores <- do.call(rbind, scores.all.dr.kl)
  colnames(iscores) <- "Iscores"
  return(iscores)

}






#' Combined projection forests
combine2Forests <- function(mod1, mod2) {
  # make several tests, to see if models fit together?
  # (formula, num.independent.variables, mtry, min.node.size, splitrule, treetype, call?, importance.mode, num.samples, replace, ...)

  res = mod1

  res$num.trees = res$num.trees + mod2$num.trees
  res$inbag.counts = c(res$inbag.counts, mod2$inbag.counts)
  res$forest$child.nodeIDs = c(res$forest$child.nodeIDs, mod2$forest$child.nodeIDs)

  for (i in 1:mod1$num.trees) {
    res$forest$split.varIDs[[i]] <- res$var[res$forest$split.varIDs[[i]]+1]
    res$forest$split.varIDs[[i]][res$forest$child.nodeIDs[[i]][[1]]==0] <- 0
  }
  for (i in 1:mod2$num.trees) {
    mod2$forest$split.varIDs[[i]] <- mod2$var[mod2$forest$split.varIDs[[i]]+1]
    mod2$forest$split.varIDs[[i]][mod2$forest$child.nodeIDs[[i]][[1]]==0] <- 0
  }

  res$forest$split.varIDs = c(res$forest$split.varIDs, mod2$forest$split.varIDs)
  res$forest$split.values = c(res$forest$split.values, mod2$forest$split.values)
  if (!is.null(res$forest$terminal.class.counts)) {
    res$forest$terminal.class.counts = c(res$forest$terminal.class.counts, mod2$forest$terminal.class.counts)
    #res$forest$terminal.class.counts = c(res$forest$terminal.class.counts, mod2$forest$terminal.class.counts)
  }
  res$forest$num.trees = res$forest$num.trees + mod2$forest$num.trees
  res$call$num.trees = res$num.trees
  res$num.independent.variables <- length(res$full.vars)
  res$forest$independent.variable.names <- res$full.vars
  res$forest$is.ordered <- rep(TRUE, length(res$full.vars))
  res$var <- 0:(length(res$full.vars)-1)
  res
}

combineForests <- function(list.rf) {
  Reduce(list.rf, f = combine2Forests)
}

#' Prediction for crf
#' @param object a crf object
#' @param Z a matrix of candidate points
#' @param cutoff the cutoff for the ensemble
#' @param return.regions should the individual indicators be returned?
#' @export
predict.crf <- function(object, Z, cutoff = NULL, return.regions = FALSE, alpha=0.05, combined = FALSE) {



  if (!is.matrix(Z)) {
    stop("Z should be a matrix.")
  }

  if(is.null(cutoff)){
    cutoff <- object$optimal.oob.cutoff
  }

  if (is.null(object$list.vars)){

    nodes.Z <- predict(object$forest, data = data.frame(X=Z),
                       type = "terminalNodes", num.threads =1)$predictions
  } else if (combined) {

    nodes.Z <- predict(object$list.rf, data = data.frame(X=Z),
                       type = "terminalNodes")$predictions
  } else {

    nodes.Z <- lapply(1:length(object$list.rf), function(i){

      z <- data.frame(X=na.omit(Z[,object$list.vars[[i]], drop=F]))

      if(nrow(z)==0){
        return(NA)
      }else{
        return(predict(object$list.rf[[i]], z,
                       type="terminalNodes", num.threads =1)$predictions)
      }
    })


    # obs.per.tree <- lapply(1:length(object$list.rf), function(i) which(complete.cases(Z[,object$list.vars[[i]], drop=F])))
    #
    #
    # nodes.Z <- lapply(1:length(object$list.rf), function(i){
    #   nodes.forall.Z <- rep(NA, nrow(Z))
    #   nodes.forall.Z[obs.per.tree[[i]]] <- nodes.Z[[i]]
    #   return(nodes.forall.Z )
    # }
    # )
    #
    nodes.Z <- Reduce(nodes.Z, f = function(x,y) cbind(x,y))# nodes.Z  <- t(do.call(rbind,nodes.Z))#
  }

  indicators <- sapply(1:ncol(nodes.Z), function(i) as.numeric(nodes.Z[,i] %in% object$confidence.nodes[[i]]))


  #tau <- apply(indicators, 1, mean)
  if (is.null(object$weights)){
    object$weights<-rep(1/ncol(indicators), ncol(indicators))

  }

  if (is.null(dim(indicators))){
    indicators <- matrix(indicators, nrow=1)
    tau <- mean(indicators)
    tauweighted <- sum(indicators*object$weights  )
    ids <- tau > cutoff
  }else {


    tau <- apply(indicators, 1, mean)
    ids <-  tau > cutoff


    tauweighted <- apply(indicators*object$weights, 1, sum)
    idsweighted <-  tau > cutoff

  }


  # if (is.null(cutoff)){
  #   # Imputation
  #   ids <- which.max(apply(indicators, 1, mean))
  # } else {
  # ids <- apply(indicators, 1, mean) > cutoff
  # }


  return(list(#Zhat = Z[ids, ],
    tau=tau, tauweighted=tauweighted))#, class = as.numeric(ids), indicators = if (return.regions) indicators else NULL))
}






Iscores <-function(imputations,
                   methods,
                   X.NA,
                   m = length(imputations[[1]]),
                   num.proj = 100,
                   num.trees.per.proj = 5,
                   min.node.size = 5,
                   sample.splitting = TRUE,
                   frac.test = 0.5,
                   n.cores = 1,
                   ...) {


  if (!is.null(names(imputations))) {
    if (!identical(names(imputations),methods)) {
      stop("imcompatible names between imputations and methods.")
    }
    if (length(imputations)!=length(methods)) {
      stop("different lengths of imputations and methods.")
    }
  } else {
    names(imputations) <- methods
  }

  ## original data
  # candidate missing value points
  ind.candidates <- which(!complete.cases(X.NA))
  ind.candidates <- sort(ind.candidates)
  nrofmissing <- sum(is.na(X.NA))

  # save the metrics
  scores.all.dr.kl<- list()


  # Initialize
  for (method in methods){
    scores.all.dr.kl[[method]] <- rep(NA, m)
  }


  # get the NA patterns
  NA.pat <- X.NA
  NA.pat[!is.na(NA.pat)] <- 1
  NA.pat.unique <- unique(NA.pat)
  NA.pat.groups <- apply(NA.pat.unique,
                         1,
                         function(p) which(apply(NA.pat, 1, function(pp) identical(pp,p))))

  # remove the fully observed patterns from the data
  V<-apply(NA.pat.unique,
           1, sum, na.rm=T)
  NA.pat.unique <- NA.pat.unique[which(V!=ncol(X.NA)),,drop=F]
  NA.pat.groups <- NA.pat.groups[which(V!=ncol(X.NA))]


  for (method in methods){

    scores.all <- list()

    # switch cases
    ## Continue here with checking!!
    for (part in 1:2) {

      # apply on the unique patterns
      scores.all.part <- mclapply(1:nrow(NA.pat.unique), function(j){

        # if sample splitting is chosen, then we separe in two parts
        if (sample.splitting) {
          # in case only one point is there we do not separate
          if (length(NA.pat.groups[[j]])==1) {
            ids.pattern.test <- NA.pat.groups[[j]]
            ids.pattern.train <- 1:nrow(X.NA)
          } else {
            if (part == 1) {
              ids.pattern.test <- NA.pat.groups[[j]][1:(floor(length(NA.pat.groups[[j]])*frac.test))]
            } else {
              ids.pattern.test <- NA.pat.groups[[j]][-c(1:(floor(length(NA.pat.groups[[j]])*frac.test)))]
            }
            ids.pattern.train <- setdiff(1:nrow(X.NA), ids.pattern.test)
          }
        } else {
          ids.pattern.test <- NA.pat.groups[[j]]
          ids.pattern.train <- 1:nrow(X.NA)
        }

        scoredr <- list()

        scores <- lapply(1:m, function(set){


          X.h <-imputations[[method]][[set]] #imputations[[method]][[sample((1:m)[-set],1)]]
          X.h <- as.matrix(X.h)[ids.pattern.train,,drop=F]
          #set.seed(1)
          # Delete Comment: This is the only place where we would get
          # different results from the old code
          object.dr <- densityRatioScore(X = X.NA[ids.pattern.train,,drop=F],
                                         Xhat = X.h,
                                         x =  NA.pat.unique[j,],
                                         X.true = X[ids.pattern.train,,drop=F],
                                         compute.glm = TRUE,
                                         num.proj=num.proj,
                                         num.trees.per.proj = num.trees.per.proj,
                                         min.node.size = min.node.size)


          ## Define the test set!
          Z <-  as.matrix(imputations[[method]][[set]])[ids.pattern.test,, drop=F]
          Z <- unname(Z)
          Z <- as.matrix(Z)



          scoredr.kl <- sum(compute_drScore(object = object.dr, Z = Z, type = "kl.score"))



          return(list(scoredr.kl = scoredr.kl))

        }
        )

        scoredr.kl <- unlist(lapply(scores, function(l) l$scoredr.kl))


        return(list(scoredr.kl =scoredr.kl))
      }, mc.cores = n.cores)

      scores.all[[part]] <- scores.all.part
    }


    dat.scoredr.kl <- do.call(rbind, lapply(unlist(scores.all, recursive = FALSE), function(l) l$scoredr.kl))
    dat.scoredr.kl <- colSums(dat.scoredr.kl) / sum(lengths(NA.pat.groups))

    scores.all.dr.kl[[method]] <- dat.scoredr.kl

  }

  iscores <- do.call(rbind, scores.all.dr.kl)
  colnames(iscores) <- "Iscores"
  return(iscores)

}






#' Combined projection forests
combine2Forests <- function(mod1, mod2) {
  # make several tests, to see if models fit together?
  # (formula, num.independent.variables, mtry, min.node.size, splitrule, treetype, call?, importance.mode, num.samples, replace, ...)

  res = mod1

  res$num.trees = res$num.trees + mod2$num.trees
  res$inbag.counts = c(res$inbag.counts, mod2$inbag.counts)
  res$forest$child.nodeIDs = c(res$forest$child.nodeIDs, mod2$forest$child.nodeIDs)

  for (i in 1:mod1$num.trees) {
    res$forest$split.varIDs[[i]] <- res$var[res$forest$split.varIDs[[i]]+1]
    res$forest$split.varIDs[[i]][res$forest$child.nodeIDs[[i]][[1]]==0] <- 0
  }
  for (i in 1:mod2$num.trees) {
    mod2$forest$split.varIDs[[i]] <- mod2$var[mod2$forest$split.varIDs[[i]]+1]
    mod2$forest$split.varIDs[[i]][mod2$forest$child.nodeIDs[[i]][[1]]==0] <- 0
  }

  res$forest$split.varIDs = c(res$forest$split.varIDs, mod2$forest$split.varIDs)
  res$forest$split.values = c(res$forest$split.values, mod2$forest$split.values)
  if (!is.null(res$forest$terminal.class.counts)) {
    res$forest$terminal.class.counts = c(res$forest$terminal.class.counts, mod2$forest$terminal.class.counts)
    #res$forest$terminal.class.counts = c(res$forest$terminal.class.counts, mod2$forest$terminal.class.counts)
  }
  res$forest$num.trees = res$forest$num.trees + mod2$forest$num.trees
  res$call$num.trees = res$num.trees
  res$num.independent.variables <- length(res$full.vars)
  res$forest$independent.variable.names <- res$full.vars
  res$forest$is.ordered <- rep(TRUE, length(res$full.vars))
  res$var <- 0:(length(res$full.vars)-1)
  res
}

combineForests <- function(list.rf) {
  Reduce(list.rf, f = combine2Forests)
}

#' Prediction for crf
#' @param object a crf object
#' @param Z a matrix of candidate points
#' @param cutoff the cutoff for the ensemble
#' @param return.regions should the individual indicators be returned?
#' @export
predict.crf <- function(object, Z, cutoff = NULL, return.regions = FALSE, alpha=0.05, combined = FALSE) {



  if (!is.matrix(Z)) {
    stop("Z should be a matrix.")
  }

  if(is.null(cutoff)){
    cutoff <- object$optimal.oob.cutoff
  }

  if (is.null(object$list.vars)){

    nodes.Z <- predict(object$forest, data = data.frame(X=Z),
                       type = "terminalNodes", num.threads =1)$predictions
  } else if (combined) {

    nodes.Z <- predict(object$list.rf, data = data.frame(X=Z),
                       type = "terminalNodes")$predictions
  } else {

    nodes.Z <- lapply(1:length(object$list.rf), function(i){

      z <- data.frame(X=na.omit(Z[,object$list.vars[[i]], drop=F]))

      if(nrow(z)==0){
        return(NA)
      }else{
        return(predict(object$list.rf[[i]], z,
                       type="terminalNodes", num.threads =1)$predictions)
      }
    })


    # obs.per.tree <- lapply(1:length(object$list.rf), function(i) which(complete.cases(Z[,object$list.vars[[i]], drop=F])))
    #
    #
    # nodes.Z <- lapply(1:length(object$list.rf), function(i){
    #   nodes.forall.Z <- rep(NA, nrow(Z))
    #   nodes.forall.Z[obs.per.tree[[i]]] <- nodes.Z[[i]]
    #   return(nodes.forall.Z )
    # }
    # )
    #
    nodes.Z <- Reduce(nodes.Z, f = function(x,y) cbind(x,y))# nodes.Z  <- t(do.call(rbind,nodes.Z))#
  }

  indicators <- sapply(1:ncol(nodes.Z), function(i) as.numeric(nodes.Z[,i] %in% object$confidence.nodes[[i]]))


  #tau <- apply(indicators, 1, mean)
  if (is.null(object$weights)){
    object$weights<-rep(1/ncol(indicators), ncol(indicators))

  }

  if (is.null(dim(indicators))){
    indicators <- matrix(indicators, nrow=1)
    tau <- mean(indicators)
    tauweighted <- sum(indicators*object$weights  )
    ids <- tau > cutoff
  }else {


    tau <- apply(indicators, 1, mean)
    ids <-  tau > cutoff


    tauweighted <- apply(indicators*object$weights, 1, sum)
    idsweighted <-  tau > cutoff

  }


  # if (is.null(cutoff)){
  #   # Imputation
  #   ids <- which.max(apply(indicators, 1, mean))
  # } else {
  # ids <- apply(indicators, 1, mean) > cutoff
  # }


  return(list(#Zhat = Z[ids, ],
    tau=tau, tauweighted=tauweighted))#, class = as.numeric(ids), indicators = if (return.regions) indicators else NULL))
}




