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
#' @examples
#' n <- 100
#' X <- cbind(rnorm(n),rnorm(n))
#' X.NA <- X
#' X.NA[,1] <- ifelse(runif(n)<=0.2, NA, X[,1])
#'
#' imputations <- list()
#'
#' imputations[[1]] <- lapply(1:5, function(i) {
#'  X.loc <- X.NA
#'  X.loc[is.na(X.NA[,1]),1] <- mean(X.NA[,1],na.rm=TRUE)
#'  return(X.loc)
#' })
#'
#' imputations[[2]] <- lapply(1:5, function(i) {
#'  X.loc <- X.NA
#'  X.loc[is.na(X.NA[,1]),1] <- sample(X.NA[!is.na(X.NA[,1]),1], size = sum(is.na(X.NA[,1])), replace = TRUE)
#'  return(X.loc)
#' })
#'
#' methods <- c("mean","sample")
#'
#' Iscores(imputations = imputations, methods = methods, X.NA = X.NA, num.proj = 10, num.trees.per.proj = 2)

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
      scores.all.part <- parallel::mclapply(1:nrow(NA.pat.unique), function(j){

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
  return(apply(iscores,1,mean))

}




