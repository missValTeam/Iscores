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
#' @param list.rf a list of forests
#' @return a forest combination of the forests stored in list.rf
combineForests <- function(list.rf) {
  Reduce(list.rf, f = combine2Forests)
}

#' Computation of the density ratio score
#' @param X a matrix of the observed data containing missing values.
#' @param Xhat a  matrix of imputations having same size as X.
#' @param x pattern of missing values.
#' @param num.proj an integer specifying the number of projections.
#' @param num.trees.per.proj an integer, the number of trees per projection.
#' @param projection.function a function providing the user-specific projections.
#' @param min.node.size the minimum number of observations in a leaf of a tree.
#' @param normal.proj, a boolean, if TRUE, sample from the NA of the pattern and
#' additionally from the non NA. If FALSE, sample only from the NA of the pattern.
#' @return a fitted random forest based on random projections
densityRatioScore <- function(X, # full data with missing values
                              Xhat, # imputation
                              x = NULL, # pattern to input
                              num.proj = 10,
                              num.trees.per.proj = 1,
                              projection.function = NULL,
                              min.node.size = 1,
                              normal.proj=T){

  M <- X
  M[!is.na(M)] <- 0
  M[is.na(M)] <- 1


  # detect the type of missing value pattern in x
  if (!is.null(x)){
    ids.x.na <- which(is.na(x))
  }else{
    ids.x.na <- 0
  }


  list.rf <- lapply(1:num.proj, FUN =function(i){

    vars <- sample.vars.proj(ids.x.na = ids.x.na,
                             X=X,
                             projection.function = projection.function,
                             normal.proj=normal.proj)



    dim.proj <- length(vars)
    X.proj.complete <- na.omit(X[,vars,drop=F])
    X.proj.complete <- matrix(X.proj.complete,ncol=length(vars), nrow=nrow(X.proj.complete), byrow = F)
    colnames(X.proj.complete) <- NULL

    ids.with.missing <- which(apply(M[,vars,drop=F], 1, function(x) sum(x)!=0))

    if(nrow(X.proj.complete) <=2){
      return(NA)
    }

    if (length(ids.with.missing) == 0) {
      return(NA)
    }

      if(normal.proj==TRUE){

        patternxA <- matrix(x[vars],ncol=length(vars))
        patternxA[is.na(patternxA)] <- 0
        patternxA <- 1-patternxA

        M.A <- M[ids.with.missing,vars,drop=F]

        ## Better: only choose the pattern that is relevant, i.e. the pattern that occurs in x_A (in the test set)
        kern <- rbfdot(sigma = 0.25)
        B <- kernelMatrix(kern, x=patternxA,
                          y=M.A)
        drawA <- which(B==1)
        Y.proj <- Xhat[ids.with.missing,vars,drop=F][drawA, ,drop=F]
      }else{
        Y.proj <- Xhat[ids.with.missing,vars,drop=F]
        drawA <- c()
      }

    #### class balancing
    cl.bl.output <- class.balancing(X.proj.complete = X.proj.complete,
                                    Y.proj = Y.proj,
                                    drawA = drawA,
                                    Xhat = Xhat,
                                    ids.with.missing =ids.with.missing,
                                    vars=vars)

    X.proj.complete <- cl.bl.output$X.proj.complete
    Y.proj <- cl.bl.output$Y.proj

    if(nrow(Y.proj)==0){
      return(NA)
    }

    colnames(Y.proj) <- NULL
    Y.proj <- as.matrix(Y.proj)
    colnames(X.proj.complete) <- NULL

    #print(vars)
    #print(dim(X.proj.complete)==dim(Y.proj))


    d <- data.frame(class = c(rep(1, each=nrow(X.proj.complete)), rep(0, each=nrow(Y.proj))),
                    X=rbind(X.proj.complete, Y.proj))


    st <- tryCatch({obj <- ranger::ranger(probability = TRUE,
                                          formula = class~., data = d,
                                          num.trees = num.trees.per.proj, mtry = dim.proj,
                                          keep.inbag = TRUE, min.node.size = min.node.size)

    obj},
    error = function(e) NA)

    if (any(is.na(st))){
      warning("Forest for a projection was NA, will redo")

    }



    if (!any(is.na(st))) {

      st$var <- vars-1
      st$full.vars <- paste("X.", 1:ncol(X),sep="")
      return(st)
    }

  })


  inds <- lapply(list.rf,function(l) {
    if(length(l)==1){
      return(FALSE)
    }else{
      return(TRUE)
    }})

  list.rf <- list.rf[unlist(inds)]
  print(paste0("nr of projections ", length(list.rf)))

  object <- list()
  object$list.rf <- combineForests(list.rf = list.rf)

  return(object)
}


#' compute the density ratio score
#' @param object a crf object.
#' @param Z a matrix of candidate points.
#' @param num.proj an integer specifying the number of projections.
#' @param num.trees.per.proj an integer, the number of trees per projection.
#' @return a numeric value, the DR I-Score.
compute_drScore <- function(object, Z = Z, num.trees.per.proj, num.proj){


  preds.all.f.h <- predict(object$list.rf, data.frame(X=Z), predict.all = TRUE)$predictions


  if (length(dim(preds.all.f.h))==2) {
    preds.all.f.h <- apply(preds.all.f.h, 1, mean)
    p.f.h <- matrix(preds.all.f.h,nrow=1)
    dr.f.h <- matrix(apply(p.f.h, 2, function(p) truncProb(p) / (1-truncProb(p))),ncol=1)
    dr.f.h[!is.finite(dr.f.h)] <- 0
    kl.f.h <- log(dr.f.h)
    kl.f.h <- colMeans(kl.f.h)
    scoredr <-  kl.f.h

  }else{

    if (!num.trees.per.proj > 1) {
      preds.all.f.h  <- t(apply(preds.all.f.h,1,cumsum))[,seq(num.trees.per.proj,
                                                              num.trees.per.proj*num.proj,
                                                              num.trees.per.proj),drop=F]
      if (ncol(preds.all.f.h)<=2) {
        preds.all.f.h <- cbind(preds.all.f.h[,1],
                               matrix(apply(preds.all.f.h,1,diff),ncol=1)) / num.trees.per.proj
      } else {
        preds.all.f.h <- cbind(preds.all.f.h[,1],
                               t(apply(preds.all.f.h,1,diff))) / num.trees.per.proj
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
  }

  return(scoredr)
}

#' Truncation of probability
#' @param p a numeric value between 0 and 1 to be truncated
#' @return a numeric value, the truncated probability.
truncProb <- function(p) {
  return(pmin(pmax(p, 10^{-9}), 1-10^{-9}))
}

#' Balancing of Classes
#' @param X.proj.complete matrix with complete projected observations.
#' @param Y.proj matrix with projected imputed observations.
#' @param drawA vector of indices corresponding to current missingness pattern.
#' @param Xhat matrix of full imputed observations.
#' @param ids.with.missing vector of indices of observations with missing values.
#' @param vars vectors of variables in projection.
#' @return a list of new X.proj.complete and Y.proj.
class.balancing <- function(X.proj.complete,Y.proj,
                            drawA, Xhat,
                            ids.with.missing,vars){

  if (nrow(Y.proj) >= nrow(X.proj.complete)) {

      X.proj.complete <- X.proj.complete[sample(1:nrow(X.proj.complete), size=nrow(Y.proj), replace=T),]

  } else {

        Xhat0 <- Xhat[ids.with.missing,vars,drop=F][-drawA,,drop=FALSE]

        if ((nrow(Y.proj) < nrow(X.proj.complete)*0.75) & nrow(Xhat0)>0 ){

          ind <- sample(1:nrow(Xhat0), size=(nrow(X.proj.complete)-nrow(Y.proj)), replace=TRUE)
          Y.proj <- rbind(Y.proj, Xhat0[ind,])

        } else {
          Y.proj <- Y.proj[sample(1:nrow(Y.proj), size=nrow(X.proj.complete), replace=T),,drop=FALSE]
        }
  }

  return(list(X.proj.complete=X.proj.complete, Y.proj=Y.proj))
}

#' Sampling of Projections
#' @param ids.x.na a vector of indices corresponding to NA in the given missingness pattern.
#' @param X a matrix of the observed data containing missing values.
#' @param projection.function a function providing the user-specific projections.
#' @param normal.proj, a boolean, if TRUE, sample from the NA of the pattern and
#' additionally from the non NA. If FALSE, sample only from the NA of the pattern.
#' @return a vector of variables corresponding to the projection.
sample.vars.proj <- function(ids.x.na,
                             X,
                             projection.function = NULL,
                             normal.proj=T){

  if (is.null(projection.function)) {

    if(normal.proj==TRUE){
      if (length(ids.x.na) == 1) {
        num.var.na <- 1
        vars.na <- ids.x.na
      } else {
        num.var.na <- sample(1:length(ids.x.na),1)
        vars.na <- sample(ids.x.na, num.var.na, replace=F)
      }
      vars.na <- sort(vars.na, decreasing = F)

      vars <- vars.na #ids.x.na

      if (ncol(X) == 2) {
        dim.proj <- 1  #sample(c(0,1),1)
      } else {

        dim.proj <-sample(1:(ncol(X)-length(vars)), size = 1, replace = FALSE)
      }

      avail <- 1:ncol(X)
      avail <- avail[-vars]
      if(length(avail)==1){
        if(dim.proj==1){
          vars <- c(vars, avail)
        }}else{
          vars <- c(vars,sample(avail, size = dim.proj, replace = FALSE))
        }


      vars <- sort(vars, decreasing=F)
      #print(vars)
      # sample only from the NAs
    }else{
      if (length(ids.x.na) == 1) {
        num.var.na <- 1
        vars.na <- ids.x.na
      } else {
        num.var.na <- sample(1:length(ids.x.na),1)
        vars.na <- sample(ids.x.na, num.var.na, replace=F)
      }
      vars.na <- sort(vars.na, decreasing = F)

      vars <- vars.na #ids.x.na
    }

  } else {
    vars <- projection.function(X)
    vars <- sort(vars, decreasing=F)
  }
  return(vars)
}


#' doevaluation: compute the imputation KL-based scoring rules
#' @param imputations a list of list of imputations matrices containing no missing values of the same size as X.NA
#' @param methods a vector of characters indicating which methods are considered for imputations. It should have the same length as the list imputations.
#' @param X.NA a matrix containing missing values, the data to impute.
#' @param m the number of multiple imputation to consider, defaulting to the number of provided multiple imputations.
#' @param num.proj an integer specifying the number of projections to consider for the score.
#' @param num.trees.per.proj an integer, the number of trees per projection.
#' @param min.node.size the minimum number of nodes in a tree.
#' @param n.cores an integer, the number of cores to use.
#' @param projection.function a function providing the user-specific projections.
#' @import ranger
#' @import parallel
#' @return a vector made of the scores for each imputation method.

doevaluation <-function(imputations,
                        methods,
                        X.NA,
                        m,
                        num.proj,
                        num.trees.per.proj,
                        min.node.size,
                        n.cores = 1,
                        projection.function = NULL) {

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

  X.NA <- as.matrix(X.NA)

  colnames(X.NA)<-NULL
  rownames(X.NA)<-NULL

  # candidate missing value points
  ind.candidates <- which(!stats::complete.cases(X.NA))
  ind.candidates <- sort(ind.candidates)
  nrofmissing <- sum(is.na(X.NA))

  # Initialize
  scores.all.dr.kl<- list()
  for (method in methods){
    scores.all.dr.kl[[method]] <- rep(NA, m)
  }

  NA.pat <- X.NA
  NA.pat[!is.na(NA.pat)] <- 1
  NA.pat.unique <- unique(NA.pat)
  NA.pat.groups <- apply(NA.pat.unique,
                         1,
                         function(p) list(unname(which(apply(NA.pat, 1, function(pp) identical(pp,p))))))

  NA.pat.groups<-unlist(NA.pat.groups, recursive = F)

  # remove the fully observed patterns from the data
  V<-apply(NA.pat.unique,
           1, sum, na.rm=T)
  NA.pat.unique <- NA.pat.unique[which(V!=ncol(X.NA)),,drop=F]
  NA.pat.groups <- NA.pat.groups[which(V!=ncol(X.NA))]

  lengths.groups <- unlist(lapply(NA.pat.groups, FUN = function(l) length(l)))
  inds.groups.length1 <- which(lengths.groups ==1)
  #print(length(inds.groups.length1))

  if(length(inds.groups.length1)>1){

    # average differently for this new group
    average.diff <- TRUE
    obs.groups.length1 <- unlist(NA.pat.groups[inds.groups.length1])
    pat.groups.length1 <- ifelse(colSums(NA.pat.unique[inds.groups.length1,])>=1,1,NA)

    NA.pat.groups <- NA.pat.groups[-inds.groups.length1]
    NA.pat.unique <- NA.pat.unique[-inds.groups.length1,]

    NA.pat.groups[[length(NA.pat.groups)+1]] <- obs.groups.length1
    NA.pat.unique <- rbind(NA.pat.unique,pat.groups.length1)
    rownames(NA.pat.unique)<-NULL
  }else{
    average.diff <- FALSE
  }


  for (method in methods){

    print(paste0("Evaluating method ", method))

    dat.scoredr.kl <- mclapply(1:nrow(NA.pat.unique), function(j){

      if(average.diff==TRUE){
        #print(paste0("Pattern ", j, " out of ", nrow(NA.pat.unique) ))
        if(j==nrow(NA.pat.unique)){
          parts <- 1
        }else{
          if (length(NA.pat.groups[[j]])==1){
            parts <- 1
          }else{
            parts <- c(1,2)
          }
        }
      }else{
        if (length(NA.pat.groups[[j]])==1){
          parts <- 1
        }else{
          parts <- c(1,2)
        }
      }


      scores.all <- list()
      for (part in parts) {

        if(j!=nrow(NA.pat.unique)){
          # in case only one point is there we do not separate
          if (length(NA.pat.groups[[j]])==1) {
            ids.pattern.test <- NA.pat.groups[[j]]
            ids.pattern.train <- 1:nrow(X.NA)
            ### We actually include everything in this case
          } else {
            if (part == 1) {
              ids.pattern.test <- NA.pat.groups[[j]][1:(floor(length(NA.pat.groups[[j]])*0.5))]

            } else{
              ids.pattern.test <- NA.pat.groups[[j]][-c(1:(floor(length(NA.pat.groups[[j]])*0.5)))]
            }

            ids.pattern.train <- setdiff(1:nrow(X.NA), ids.pattern.test)

          }
        }else{

          if (average.diff==T){
            ids.pattern.train <- 1:nrow(X.NA)
            ids.pattern.test <- NA.pat.groups[[j]]
          }else{

            if (length(NA.pat.groups[[j]])==1) {
              ids.pattern.test <- NA.pat.groups[[j]]
              ids.pattern.train <- 1:nrow(X.NA)
              ### We actually include everything in this case
            } else {
              if (part == 1) {
                ids.pattern.test <- NA.pat.groups[[j]][1:(floor(length(NA.pat.groups[[j]])*0.5))]

              } else{
                ids.pattern.test <- NA.pat.groups[[j]][-c(1:(floor(length(NA.pat.groups[[j]])*0.5)))]
              }

              ids.pattern.train <- setdiff(1:nrow(X.NA), ids.pattern.test)

            }

          }
        }

        scores <- lapply(1:m, function(set){

          X.h <- imputations[[method]][[set]]
          X.h <- as.matrix(X.h)[ids.pattern.train,,drop=F]

          if(any(is.na(X.h))==TRUE){
            stop(paste0("Method ", method, " returns NAs in imputation set ", set,
                        ". Please check the imputation."))
          }

          if(average.diff==TRUE){
            if(j!=nrow(NA.pat.unique)){
              normal.proj <- TRUE
            }else{
              normal.proj <- FALSE
            }
          }else{
            normal.proj <- TRUE
          }

          #tmp1<-Sys.time()
          object.dr <-  tryCatch({obj <- densityRatioScore(X = X.NA[ids.pattern.train,,drop=F],
                                                           Xhat = X.h,
                                                           x =  NA.pat.unique[j,],
                                                           num.proj=num.proj,
                                                           num.trees.per.proj = num.trees.per.proj,
                                                           min.node.size = min.node.size,
                                                           projection.function = projection.function,
                                                           normal.proj=normal.proj)
          obj},
          error = function(e) NA)


          ## Define the test set!
          Z <-  as.matrix(imputations[[method]][[set]])[ids.pattern.test,, drop=F]
          Z <- unname(Z)
          Z <- as.matrix(Z)


          if(any(is.na(object.dr)) || any(is.null(object.dr))){
            scoredr.kl <- NA
          }else{

            if(j==nrow(NA.pat.unique)){
              if(average.diff ==FALSE){
                scoredr.kl <- mean(compute_drScore(object = object.dr, Z = Z,
                                                   num.trees.per.proj = num.trees.per.proj ,
                                                   num.proj = num.proj))
              }else{
                scoredr.kl <- unlist(compute_drScore(object = object.dr, Z = Z,
                                                     num.trees.per.proj = num.trees.per.proj ,
                                                     num.proj = num.proj))
              }

            }else{
              scoredr.kl <- mean(compute_drScore(object = object.dr, Z = Z,
                                                 num.trees.per.proj = num.trees.per.proj ,
                                                 num.proj = num.proj))
            }

          }

          return(list(scoredr.kl = scoredr.kl))

        }
        )

        scores.all.part<- unlist(lapply(scores, function(l) l$scoredr.kl))

        scores.all[[part]] <- scores.all.part

      }

      if(length(parts)==1){
        dat.scoredr.kl <- unlist(scores.all)
      }else{
        dat.scoredr.kl <- mean(unlist(scores.all), na.rm=T)

      }

      return(dat.scoredr.kl)
    }, mc.cores = n.cores)


    scores.all.dr.kl[[method]] <- mean(unlist(dat.scoredr.kl), na.rm=T)

  }

  print(paste0("done scoring"))

  return(evaluation=list(scores.all.dr.kl=scores.all.dr.kl))

}




