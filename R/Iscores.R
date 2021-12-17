#' Iscores: compute the imputation KL-based scoring rules
#' @param imputations a list of list of imputations matrices containing no missing values of the same size as X.NA
#' @param methods a vector of characters indicating which methods are considered for imputations. It should have the same length as the list imputations.
#' @param X.NA a matrix containing missing values, the data to impute.
#' @param m the number of multiple imputation to consider, defaulting to the number of provided multiple imputations.
#' @param num.proj an integer specifying the number of projections to consider for the score.
#' @param num.trees.per.proj an integer, the number of trees per projection.
#' @param min.node.size the minimum number of nodes in a tree.
#' @param n.cores an integer, the number of cores to use.
#' @param projection.function a function providing the user-specific projections.
#' @param rescale, a boolean, TRUE if the scores should be rescaled such that the max score is 0.
#' @import ranger
#' @import parallel
#' @import kernlab
#' @examples
#' n <- 100
#' X <- cbind(rnorm(n),rnorm(n))
#' X.NA <- X
#' X.NA[,1] <- ifelse(stats::runif(n)<=0.2, NA, X[,1])
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
#'  X.loc[is.na(X.NA[,1]),1] <- sample(X.NA[!is.na(X.NA[,1]),1],
#'  size = sum(is.na(X.NA[,1])), replace = TRUE)
#'  return(X.loc)
#' })
#'
#' methods <- c("mean","sample")
#'
#' Iscores(imputations,
#' methods,
#' X.NA,
#' num.proj=5
#' )
#'
#' @return a vector made of the scores for each imputation method.
#' @export
Iscores <- function(imputations,
                   methods,
                   X.NA,
                   m=length(imputations[[1]]),
                   num.proj=100,
                   num.trees.per.proj = 5,
                   min.node.size=10,
                   n.cores = 1,
                   projection.function = NULL,
                   rescale =TRUE
                   ){

  Iscores.values <- doevaluation(imputations=imputations,
                                  methods=methods,
                                  X.NA=X.NA,
                                  m=m,
                                  num.proj=num.proj,
                                  num.trees.per.proj = num.trees.per.proj,
                                  min.node.size=min.node.size,
                                  n.cores = n.cores,
                                  projection.function = projection.function)



  Iscores.values<- do.call(rbind, Iscores.values)

  if(rescale==TRUE){
    names.methods <- colnames(Iscores.values)
    Iscores.values <-lapply(Iscores.values, FUN=function(l){
      l-max(unlist(Iscores.values))
    })

    Iscores.values <- do.call(cbind,Iscores.values)
    colnames(Iscores.values) <- names.methods
  }

  return(Iscores.values)

}
