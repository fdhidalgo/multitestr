free_stepdown <- function(tstats, rep_tstats){
  tstat_order <- rank(abs(tstats), ties.method = "random")
  bs_qvals <- matrix(ncol = ncol(rep_tstats), nrow = nrow(rep_tstats))
  for (i in 1:ncol(bs_qvals)) {
    for (j in 1:length(tstats)) {
      if (j == 1) {
        bs_qvals[which(tstat_order == j), i] <-
          abs(rep_tstats[which(tstat_order == 1), i][[1]])
      } else {
        bs_qvals[which(tstat_order == j), i] <-
          max(
            abs(rep_tstats[which(tstat_order == j), i][[1]]),
            abs(bs_qvals[which(tstat_order == j-1), i][[1]])
          )
      }
    }
  }

  bs_pvals_adjusted <- rep(NA, length(tstats))
  for (i in 1:length(tstats)) {
    bs_pvals_adjusted[i] <-
      (sum(bs_qvals[i, ] >= abs(tstats[i])) + 1) / (ncol(rep_tstats) + 1)
  }

  #Enforce Monotonicity
  for (j in 0:(length(bs_pvals_adjusted)-1)) {
    if (j == 0) {
      bs_pvals_adjusted[which(tstat_order == length(tstats))] <-  bs_pvals_adjusted[which(tstat_order == length(tstats))]
    } else {
      bs_pvals_adjusted[which(tstat_order == length(tstats)-j)] <-
        max(
          bs_pvals_adjusted[which(tstat_order == length(tstats)-j)],
          bs_pvals_adjusted[which(tstat_order == length(tstats) - (j-1))]
        )
    }
  }
  bs_pvals_adjusted
}

#' Free step-down adjusted p-values
#'
#' @description
#' Generates adjusted p-values using the algorithm described in Westfall and
#'   Young (1993).
#'
#' @param full_formulas A list of objects of class \code{formula}, specifying
#'   the full (unrestricted) models.
#' @param null_formulas A list of objects of class \code{formula}, specifying
#'   the null (restricted) models.
#' @param data A data frame containing the variables in the models.
#' @param coef_list A character vector specifying the variable for each model
#'   for which the p-value is to be adjusted.
#' @param weights An optional list of weights to be used in the fitting process.
#'   Each element of the list should be a numeric vector with length equals to
#'   the number of rows in \code{data}. If \code{NULL}, each observation is
#'   given a weight of 1. Weights can be used to specify subgroups, by giving units outside the subgroup of interest a weight of 0.
#' @param cluster An optional character string indicating the column in the
#'   data frame that records the cluster to which an observation belongs (for
#'   computing standard errors that account for clustering).
#' @param nboots The number of bootstrap replicates, a single positive integer.
#' @param boot_type A character string indicating the type of resampling
#'   required. Possible values are \code{wild} or \code{pairs}.
#' @param parallel logical, indicating if parallel operation is to be used. If
#'   \code{TRUE}, a parallel backend should be registered prior to running
#'   \code{boot_stepdown}. See the documentation for \code{foreach} in the
#'   \pkg{foreach} package for more details.
#' @param pb logical. Should a progress bar be displayed if \code{parallel} is
#'   \code{FALSE}? Ignored if \code{parallel} is \code{TRUE}.
#'
#' @return A data frame reporting unadjusted and adjusted p-values for each
#'   hypothesis provided in \code{full_formulas} and \code{null_formulas}.
#' @export
#' @importFrom stats fitted formula lm model.frame na.action na.omit resid
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @references Westfall, Peter H., and S. Stanley Young. 1993. \emph{Resampling-
#' Based Multiple Testing: Examples and Methods for P-Value Adjustment.}
#' John Wiley & Sons, New York.
#'
#' @examples
#' # Replicate Casey et al 2012, Table II Column 3
#' F <- lapply(sprintf("zscore_%d ~ 0 + t + tothhs + road + ward", 1:12),
#'             as.formula)
#' pvals <- boot_stepdown(full_formulas = F,
#'                        null_formulas = lapply(F, update, . ~ . - t),
#'                        data = gobifo,
#'                        coef_list = "t",
#'                        nboots = 100, # small nboots for demonstration only
#'                        parallel = FALSE,
#'                        boot_type = "pairs")
#' dplyr::mutate_if(pvals, is.numeric, round, 3)
boot_stepdown <- function(full_formulas, null_formulas, data, coef_list,
                          weights = NULL, cluster = NULL, nboots = 10000,
                          boot_type = c("wild", "pairs"), parallel = TRUE,
                          pb = TRUE) {

  boot_type <- match.arg(boot_type)
  if (is.null(cluster)) {
    data$cluster <- as.character(1:nrow(data))
  } else {
    if (!(cluster %in% names(data)))
      stop(sprintf("`%s` is not a column in `data`", cluster))
    data$cluster <- as.character(data[, cluster])
  }

  # a 1 column dataframe of cluster indicators from original dataset
  boot_cluster <- data.frame(cluster = as.character(data$cluster),
                             stringsAsFactors = FALSE)

  if (is.null(weights))
    weights <- list(rep(1, nrow(data)))

  # Only keep variables uses to test hypotheses
  # Delete incomplete cases
  # Append weights
  datasets <- mapply(
    null_formulas, full_formulas, list(data), weights,
    FUN = function(x, y, z, w) {
      dat <- na.omit(z[, unique(c(all.vars(x), all.vars(y), "cluster"))])
      if (is.null(na.action(dat)))
        dat$weights <- w
      else
        dat$weights <- w[-na.action(dat)]
      return(dat)
    },
    SIMPLIFY = FALSE)

  # cluster indicators for each dataset (can vary across datasets)
  clusters <- lapply(datasets, function(x) x$cluster)
  weights <- lapply(datasets, function(x) x$weights)

  # Fit full and null models
  full_models <- mapply(x = full_formulas, y = datasets, w = weights,
                        FUN = function(x, y, w) {
                          newEnv <- new.env(parent = environment(x))
                          newEnv$w <- w
                          environment(x) <- newEnv
                          lm(x, data = y, weights = w)
                        },
                        SIMPLIFY = FALSE)
  null_models <- mapply(null_formulas, datasets, weights,
                        FUN = function(x, y, w) {
                          newEnv <- new.env(parent = environment(x))
                          newEnv$w <- w
                          environment(x) <- newEnv
                          lm(x, data = y, weights = w)
                        },
                        SIMPLIFY = FALSE)

  # append cluster indicators
  for(i in 1:length(full_models)){
    full_models[[i]]$cluster <- as.character(datasets[[i]]$cluster)
  }

  if(boot_type == "wild"){
    gen_boot <- gen_wildbs
  } else if(boot_type == "pairs"){
    gen_boot <- gen_pairsbs
  } else {
    stop("boot_type must be one of ``wild'' or ``pairs''")
  }

  tstats <- gen_boot(null_models = null_models,
                     full_models = full_models,
                     cluster_list = boot_cluster,
                     clusters = clusters,
                     coef_list = coef_list,
                     permute = FALSE)

  if (parallel == FALSE) {
    rep_tstats <- matrix(nrow = length(tstats), ncol = nboots)
    if (pb) progress_bar <- txtProgressBar(1, nboots)
    for(i in 1:nboots) {
      rep_tstats[,i] <- gen_boot(null_models = null_models,
                                 full_models = full_models,
                                 cluster_list = boot_cluster,
                                 clusters = clusters,
                                 coef_list = coef_list)
      if (pb) setTxtProgressBar(progress_bar, i)
    }
    if (pb) cat("\n")
  } else {
    `%dopar%` <- foreach::`%dopar%`
    rep_tstats <- foreach::foreach(i = 1:nboots,
                                   .combine = 'cbind',
                                   .packages = c('sandwich', 'lmtest')) %dopar% {
                                     gen_boot(null_models = null_models,
                                              full_models = full_models,
                                              clusters = clusters,
                                              cluster_list = boot_cluster,
                                              coef_list = coef_list)
                                   }
  }

  #Generate bootstrap unadjusted p-values
  bs_pvals <- rep(NA, length(tstats))
  for(i in 1:length(tstats)){
    bs_pvals[i] <- (sum(abs(rep_tstats[i, ]) >= abs(tstats[i])) + 1) / (nboots + 1)
  }
  bs_pvalues_adjusted <- free_stepdown(tstats = tstats,
                                       rep_tstats = rep_tstats)

  data.frame(Hypothesis = paste0("Hypothesis ", 1:length(full_formulas)),
             Variable = coef_list,
             bs_pvalues_unadjusted = bs_pvals,
             bs_pvalues_adjusted = bs_pvalues_adjusted)
}
