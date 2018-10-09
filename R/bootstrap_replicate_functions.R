gen_wildbs <- function(null_models, full_models, cluster_list, clusters,
                       coef_list, permute = TRUE) {
  #Wild Parametric Bootstrap. The null and full models are
  #equal-sized lists of formulas that specify the models to be
  #tested. The argument 'coef_list' specifies the coefficients for
  #which a p-value needs to be generated in the full model (only
  #one per model)

  resids_null <- lapply(null_models, resid)
  fitted_null <- lapply(null_models, fitted)

  ## ------------------------------------------------------------------------
  ## Sampling
  if (permute) {
    # Rademacher distribution
    wild_perm <- dplyr::summarise(dplyr::group_by(cluster_list, cluster),
                                  wild_perm = sample(c(1, -1), 1))
  } else {
    wild_perm <- dplyr::summarise(dplyr::group_by(cluster_list, cluster),
                                  wild_perm = 1)
  }
  ## ------------------------------------------------------------------------

  bs_data <- lapply(full_models, model.frame) # one for each model

  # loop over hypotheses to test
  for (i in 1:length(bs_data)) {
    bs_data[[i]]$cluster <- clusters[[i]]
    bs_data[[i]] <- dplyr::left_join(bs_data[[i]], wild_perm,
                                     by = "cluster")
    bs_data[[i]][, 1] <- fitted_null[[i]] +
      resids_null[[i]] * bs_data[[i]]$wild_perm
  }

  # Estimate models on bootstrapped data
  # The new.env construct is needed because `lm` cannot find weights otherwise
  models_bs <- mapply(full_models, bs_data,
                      FUN = function(x, y) {
                        f <- formula(x)
                        newEnv <- new.env(parent = environment(f))
                        newEnv$w <- y$`(weights)`
                        lm(formula = f, data = y, weights = w)
                      },
                      SIMPLIFY = FALSE)

  for(i in 1:length(models_bs)){
    models_bs[[i]]$cluster <- bs_data[[i]]$cluster
  }
  #Estimate cluster robust VCOV and pvals
  bs_vcovs <- lapply(models_bs, function(x) {
    lmtest::coeftest(x, sandwich::vcovCL(x , x$cluster))
  })
  mapply(bs_vcovs, coef_list, FUN = function(x, y) x[y ,'t value'])
}

gen_pairsbs <- function(null_models, full_models, cluster_list, clusters,
                        coef_list, permute = TRUE){
  #Pairs Non-Parametric Bootstrap. The null and full models are
  #equal-sized lists of formulas that specify the models to be
  #tested. The argument 'coef_list' specifies the coefficients for
  #which a p-value needs to be generated in the full model (only
  #one per model)

  ## cluster_list is a dataframe with 1 column

  if (permute) {
    bs_clusters <- dplyr::sample_frac(unique(cluster_list), replace = TRUE)
    bs_clusters$cluster <- as.character(bs_clusters$cluster)
    bs_clusters$bs_cluster <- as.character(1:nrow(bs_clusters))
    #    indx <- sample(1:nrow(model.frame(full_models[[1]])), nrow(model.frame(full_models[[1]])), replace = TRUE)
  } else {
    bs_clusters <- unique(cluster_list)
    bs_clusters$cluster <- as.character(bs_clusters$cluster)
    bs_clusters$bs_cluster <- as.character(1:nrow(bs_clusters))
    #    indx <- 1:nrow(model.frame(full_models[[1]]))
  }


  bs_data <- lapply(full_models, model.frame)
  # loop over hypotheses to test
  for(i in 1:length(bs_data)){
    bs_data[[i]]$cluster <- as.character(clusters[[i]])
    bs_data[[i]] <- dplyr::left_join(bs_clusters, bs_data[[i]], by = "cluster")
    #    bs_data[[i]] <- bs_data[[i]][indx, ]
  }

  # Estimate models on bootstrapped data
  models_bs <- mapply(full_models, bs_data,
                      FUN = function(x, y) {
                        f <- formula(x)
                        newEnv <- new.env(parent = environment(f))
                        newEnv$w <- y$`(weights)`
                        lm(formula = f,
                           data = y,
                           weights = w)
                      },
                      SIMPLIFY = FALSE)
  for(i in 1:length(models_bs)){
    #    models_bs[[i]]$cluster <- 1:nrow(bs_data[[i]])
    models_bs[[i]]$cluster <- bs_data[[i]]$bs_cluster
  }
  #Estimate cluster robust VCOV and pvals
  bs_vcovs <- lapply(models_bs, function(x) {
    lmtest::coeftest(x, sandwich::vcovCL(x , x$cluster))
  })
  full_vcovs <- lapply(full_models, function(x) {
    lmtest::coeftest(x, sandwich::vcovCL(x , x$cluster))
  })
  if(permute == FALSE){
    return(mapply(bs_vcovs, coef_list, FUN = function(x, y) {
      x[y ,'t value']
    }))
  }
  if(permute == TRUE){
    return(
      (mapply(bs_vcovs, coef_list, FUN = function(x, y) x[y ,'Estimate']) -
         mapply(full_vcovs, coef_list, FUN = function(x, y) x[y ,'Estimate'])) /
        mapply(bs_vcovs, coef_list, FUN = function(x, y) x[y ,'Std. Error'])
    )
  }

}
