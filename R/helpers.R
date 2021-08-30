#' Combining projection forests
#' @param mod1 first forest
#' @param mod2 second forest
#' @return a new forest combining the first and the second forest
#' @import stats
combine2Forests <- function(mod1, mod2) {

  res <- mod1

  res$num.trees <- res$num.trees + mod2$num.trees
  res$inbag.counts <- c(res$inbag.counts, mod2$inbag.counts)
  res$forest$child.nodeIDs <- c(res$forest$child.nodeIDs, mod2$forest$child.nodeIDs)

  for (i in 1:mod1$num.trees) {
    res$forest$split.varIDs[[i]] <- res$var[res$forest$split.varIDs[[i]]+1]
    res$forest$split.varIDs[[i]][res$forest$child.nodeIDs[[i]][[1]]==0] <- 0
  }
  for (i in 1:mod2$num.trees) {
    mod2$forest$split.varIDs[[i]] <- mod2$var[mod2$forest$split.varIDs[[i]]+1]
    mod2$forest$split.varIDs[[i]][mod2$forest$child.nodeIDs[[i]][[1]]==0] <- 0
  }

  res$forest$split.varIDs <- c(res$forest$split.varIDs, mod2$forest$split.varIDs)
  res$forest$split.values <- c(res$forest$split.values, mod2$forest$split.values)
  if (!is.null(res$forest$terminal.class.counts)) {
    res$forest$terminal.class.counts <- c(res$forest$terminal.class.counts, mod2$forest$terminal.class.counts)
  }
  res$forest$num.trees <- res$forest$num.trees + mod2$forest$num.trees
  res$call$num.trees <- res$num.trees
  res$num.independent.variables <- length(res$full.vars)
  res$forest$independent.variable.names <- res$full.vars
  res$forest$is.ordered <- rep(TRUE, length(res$full.vars))
  res$var <- 0:(length(res$full.vars)-1)
  return(res)
}

#' Combining a list of forest
#'
#' @param list.rf a list of forests
#' @return a forest combination of the forests stored in list.rf
combineForests <- function(list.rf) {
  Reduce(list.rf, f = combine2Forests)
}

#' Prediction for crf
#' @param object a crf object
#' @param Z a matrix of candidate points
#' @param cutoff the cutoff for the ensemble
#' @param return.regions should the individual indicators be returned?
#' @param alpha the level
#' @param combined the option to have the forests combined or not
#' @return a list of the predictions from a crf object.
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

      z <- data.frame(X=stats::na.omit(Z[,object$list.vars[[i]], drop=F]))

      if(nrow(z)==0){
        return(NA)
      }else{
        return(predict(object$list.rf[[i]], z,
                       type="terminalNodes", num.threads =1)$predictions)
      }
    })

    nodes.Z <- Reduce(nodes.Z, f = function(x,y) cbind(x,y))
  }

  indicators <- sapply(1:ncol(nodes.Z), function(i) as.numeric(nodes.Z[,i] %in% object$confidence.nodes[[i]]))

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




  return(list(tau=tau, tauweighted=tauweighted))
}


#' Computation of the density ratio score
#' @param X a matrix of the observed data containing missing values
#' @param Xhat a  matrix of imputations having same size as X
#' @param x pattern of missing values
#' @param X.true the actual true data if available
#' @param num.proj an integer specifying the number of projections to consider for the score.
#' @param num.trees.per.proj an integer, the number of trees per projection.
#' @param projection.function a function providing the possible projections
#' @param min.node.size the minimum number of nodes in a tree.
#' @param multiclass a boolean
#' @param compute.glm a boolean
#' @param ... additional parameters#'
#' @return a crf object
densityRatioScore <- function(X, # full data with missing values
                              Xhat, # imputation
                              x = NULL, # pattern to input
                              X.true = NULL, # actual true data
                              num.proj = 10,
                              num.trees.per.proj = 1,
                              projection.function = NULL,
                              min.node.size = 1,
                              multiclass = FALSE,
                              compute.glm = FALSE,
                              ...){


  if (multiclass & num.trees.per.proj > 1) {
    stop("more trees per proj not implemented yet for multiclass.")
  }


  M <- X
  M[!is.na(M)] <- 0
  M[is.na(M)] <- 1


  # detect the type of missing value pattern in x
  if (!is.null(x)){
    ids.x.na <- which(is.na(x))
  }else{
    ids.x.na <- 0
  }

  list.rf <- list()
  list.rf.oracle <- list()
  list.glm <- list()
  list.vars <- list()
  list.ids.NA <- list()
  list.trueindices <- list()


  i <- 0

  while (length(list.rf) != num.proj){

    i <- i + 1

    #sampling from the missing coordinates (at least one should be present)

    if (is.null(projection.function)) {
      if (sum(ids.x.na) > 0){


        if (length(ids.x.na) == 1) {
          num.var.na <- 1
          vars.na <- ids.x.na
        } else {
          num.var.na <- sample(1:length(ids.x.na),1)
          vars.na <- sample(ids.x.na, num.var.na, replace=F)
        }
        vars.na <- sort(vars.na, decreasing = F)

        vars <- vars.na
        if(length(ids.x.na)==ncol(X)){
          dim.proj <- 0
        }else{

          if (ncol(X) == 2) {
            dim.proj <- sample(c(0,1),1)
          } else {
            dim.proj <- sample(0:(ncol(X)-length(vars)), size = 1, replace = FALSE)
          }

          avail <- 1:ncol(X)
          avail <- avail[-vars]
          if(length(avail)==1){
            if(dim.proj==1){
              vars <- c(vars, avail)
            }}else{
              vars <- c(vars,sample(avail, size = dim.proj, replace = FALSE))
            }

        }
        vars <- sort(vars, decreasing=F)

      } else {

        # 1) subsampling of variables
        if(ncol(X)==2){
          dim.proj <- 2
        }else{
          dim.proj <- sample(2:ncol(X), size = 1, replace = FALSE)
        }

        vars <- sample(1:ncol(X), size = dim.proj, replace = FALSE)
        vars <- sort(vars, decreasing=F)

      }
    } else {
      vars <- projection.function(x)
      vars <- sort(vars, decreasing=F)
    }




    # from here on we should have a projection made of vars.


    dim.proj <- length(vars)


    # 2) complete case selection + adversarial sample construction
    X.proj.complete <- stats::na.omit(X[,vars,drop=F])


    X.proj.complete <- as.matrix(X.proj.complete)
    colnames(X.proj.complete) <- NULL

    list.ids.NA[[length(list.ids.NA)+1]] <- which(apply(X[,vars,drop=F], 1, function(x) any(is.na(x))))


    # as adversarial we bind the most of the data with missing values + additional with no missing values
    ids.with.missing <- which(apply(M[,vars,drop=F], 1, function(x) sum(x)!=0))


    if (multiclass) {
      Y.proj <- lapply(Xhat, function(l) {
        if (length(ids.with.missing) != 0) {

          Y.proj <- l[ids.with.missing,vars,drop=F]

          if (nrow(Y.proj) >= nrow(X.proj.complete)) {

            Y.proj <- Y.proj[1:nrow(X.proj.complete),,drop=F]

          } else {

            Y.proj <- rbind(Y.proj, l[-ids.with.missing, vars,drop=F][1:(nrow(X.proj.complete)-nrow(Y.proj)),,drop=F])

          }
        } else {
          Y.proj <- X.proj.complete
        }
        colnames(Y.proj) <- NULL
        Y.proj <- as.matrix(Y.proj)
        return(Y.proj)
      })
    } else {
      if (length(ids.with.missing) != 0) {

        Y.proj <- Xhat[ids.with.missing,vars,drop=F]

        if (nrow(Y.proj) >= nrow(X.proj.complete)) {

          Y.proj <- Y.proj[1:nrow(X.proj.complete),,drop=F]

        } else {

          Y.proj <- rbind(Y.proj, Xhat[-ids.with.missing, vars,drop=F][1:(nrow(X.proj.complete)-nrow(Y.proj)),,drop=F])

        }
      } else {
        Y.proj <- X.proj.complete
      }
      colnames(Y.proj) <- NULL
      Y.proj <- as.matrix(Y.proj)
    }




    # storing proj
    list.vars[[length(list.vars)+1]] <- vars
    list.trueindices[[length(list.trueindices)+1]] <- which(stats::complete.cases(X[,vars]))
    # fitting single tree

    colnames(X.proj.complete) <- NULL

    if (multiclass) {
      d <- data.frame(class = rep(1:(length(Xhat)+1), each = nrow(X.proj.complete)),
                      x=rbind(X.proj.complete, Reduce(x = Y.proj, f = rbind)))
    } else {
      d <- data.frame(class = c(rep(1, each=nrow(X.proj.complete)), rep(0, each=nrow(Y.proj))),
                      x=rbind(X.proj.complete, Y.proj))
    }

    st <- tryCatch({obj <- ranger::ranger(probability = TRUE,
                                          formula = class~., data = d,
                                          num.trees = num.trees.per.proj, mtry = dim.proj,
                                          keep.inbag = TRUE, min.node.size = min.node.size)

    obj},
    error = function(e) NA)


    if (compute.glm) {
      gl <- tryCatch({obj <- stats::glm(formula = class~., data = d,family = "binomial")
      obj},
      error = function(e) NA)
    }

    if (!is.null(X.true)) {

      if (multiclass) {
        stop("multiclass not implemented for oracle")
      }

      if (length(ids.with.missing) != 0) {
        Y.proj.oracle <- Xhat[ids.with.missing,vars,drop=F]
        X.true.oracle <- X.true[ids.with.missing,vars,drop=F]
      } else {
        Y.proj.oracle <- Xhat[,vars,drop=F]
        X.true.oracle <- X.true[,vars,drop=F]
      }

      Y.proj.oracle <- as.matrix(Y.proj.oracle)
      colnames(Y.proj.oracle) <- NULL
      X.true.oracle <- as.matrix(X.true.oracle)
      colnames(X.true.oracle) <- NULL

      d.oracle <- data.frame(class = c(rep(1, each=nrow(X.true.oracle)),
                                       rep(0, each=nrow(Y.proj.oracle))),
                             x=rbind(X.true.oracle, Y.proj.oracle))

      st.oracle <- tryCatch({obj <- ranger::ranger(probability = TRUE,
                                                   formula = class~., data = d.oracle,
                                                   num.trees = num.trees.per.proj, mtry = dim.proj,
                                                   keep.inbag = TRUE, min.node.size = min.node.size)
      obj},
      error = function(e) NA)
    }

    if (!any(is.na(st))) {


      # store var corresponding
      st$var <- vars-1
      st$full.vars <- paste("X.", 1:ncol(X),sep="")

      # storing tree
      list.rf[[length(list.rf)+1]] <- st

    }

    if (compute.glm && !any(is.na(gl))) {

      # store var corresponding
      gl$var <- vars-1
      gl$full.vars <- paste("X.", 1:ncol(X),sep="")

      # storing tree
      list.glm[[length(list.glm)+1]] <- gl
    }

    if (!is.null(X.true) && !any(is.na(st.oracle))) {


      # store var corresponding
      st.oracle$var <- vars-1
      st.oracle$full.vars <- paste("X.", 1:ncol(X),sep="")

      # storing tree
      list.rf.oracle[[length(list.rf.oracle)+1]] <- st.oracle

    }


  }



  object <- list(X = X,
                 list.rf = list.rf,
                 list.rf.oracle = list.rf.oracle,
                 list.glm = list.glm,
                 list.vars = list.vars,
                 list.ids.NA = list.ids.NA,
                 list.trueindices = list.trueindices,
                 num.trees.per.proj = num.trees.per.proj,
                 num.proj = num.proj)

  # combine the forests together
  object$list.rf <- combineForests(list.rf = object$list.rf)
  if (!is.null(X.true)) {
    object$list.rf.oracle <- combineForests(list.rf = object$list.rf.oracle)
  }



  class(object) <- "crf"



  return(object)
}



#' compute the density ratio score
#' @param object a crf object
#' @param Z a matrix of candidate points
#' @param multiclass a boolean
#' @param type a character
#' @return a numeric value, the score.s
compute_drScore <- function(object, Z = Z, multiclass = FALSE, type = "simple"){

  if (multiclass & object$num.trees.per.proj > 1) {
    stop("more trees per proj not implemented yet for multiclass.")
  }

  if (multiclass) {
    preds.all.f.h <- lapply(Z, function(ZZ) predict(object$list.rf, data.frame(X=ZZ), predict.all = TRUE)$predictions)
  } else {

    if (type == "kl.score") {
      preds.all.f.h <- predict(object$list.rf, data.frame(X=Z), predict.all = TRUE)$predictions
    } else if (type == "kl.oracle.proj") {
      preds.all.f.h <- predict(object$list.rf.oracle, data.frame(X=Z), predict.all = TRUE)$predictions
    }
  }


  if (length(dim(preds.all.f.h))==2) {
    preds.all.f.h <- apply(preds.all.f.h, 1, mean)
    p.f.h <- matrix(preds.all.f.h,nrow=1)
    dr.f.h <- matrix(apply(p.f.h, 2, function(p) truncProb(p) / (1-truncProb(p))),ncol=1)
    dr.f.h[!is.finite(dr.f.h)] <- 0
    kl.f.h <- log(dr.f.h)
    kl.f.h <- colMeans(kl.f.h)
    scoredr <-  kl.f.h
    return(scoredr)
  }


  if (!object$num.trees.per.proj > 1) {
    preds.all.f.h  <- t(apply(preds.all.f.h,1,cumsum))[,seq(object$num.trees.per.proj,
                                                            object$list.rf$num.trees,
                                                            object$num.trees.per.proj),drop=F]
    if (ncol(preds.all.f.h)<=2) {
      preds.all.f.h <- cbind(preds.all.f.h[,1],
                             matrix(apply(preds.all.f.h,1,diff),ncol=1)) / object$num.trees.per.proj
    } else {
      preds.all.f.h <- cbind(preds.all.f.h[,1],
                             t(apply(preds.all.f.h,1,diff))) / object$num.trees.per.proj
    }
  }
    if (dim(preds.all.f.h)[1] == 1) {
      p.f.h <- matrix(preds.all.f.h[1,1,],nrow=1)
      dr.f.h <- matrix(apply(p.f.h, 2, function(p) truncProb(p) / (1-truncProb(p))),ncol=1)
    } else {
      p.f.h <- preds.all.f.h[,1,]
      dr.f.h  <- t(apply(p.f.h, 2, function(p) truncProb(p) / (1-truncProb(p))))
    }
    dr.f.h[!is.finite(dr.f.h)] <- 0
    kl.f.h <- log(dr.f.h)
    kl.f.h <- colMeans(kl.f.h)
    scoredr <-  kl.f.h


  return(scoredr)
}

#' Truncation of probability
#' @param p a numeric value between 0 and 1 to be truncated
#' @return a numeric value, the truncated probability.
truncProb <- function(p) {
  return(pmin(pmax(p, 10^{-9}), 1-10^{-9}))
}


