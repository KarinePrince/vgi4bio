#' @title Short Title of the Function (One Line)
#'
#' @author Karine Princé
#' @references [Giraud et al. 2015](https://doi.org/10.1111/biom.12431)
#' @description A longer description of what the function does (several lines)
#' @import stats
#' @param x a data.frame containing original species count data. The data.frame should have 4 columns with headers, namely:
#'  "species": the variable that stores the name/code of the species,
#'  "site": the variable that stores the observation sites,
#'  "scheme": the variable that stores the monitoring program to which the observation data are associated
#'  "effort" a measure of effort that depends on the monitoring program (for example, it can be set to the number of survey points-years implemented in the standardardized program for the site stored in the variable Site)
#' #' @return What the function returns
#'
#'
#' @examples
#'
#' #' @export
#'



## estimation of abundance based on opportunistic survey data combined to standardized data

fun_abund_estim <- function(x, family="quasipoisson", datapred=NULL, type = "response"){

  ## Define


  ## run model
  cat("\nRun Abundance Model\n=======================\n")
  mod <- glm(N_obs ~ species*scheme + species*site + site*scheme, offset=log(effort),
             data = d, family= family)
  return(list(summary(mod),confint(mod)))

  ## look at predictive power of the model
  if(datapred == NULL){
    N_comb <- NULL
  } else {
    cat("\nCompute Abundance Predictions\n=======================\n")
    N_comb <- predict(mod, newdata = datapred, type = type)
  }
}


## function to calculate the median, the first and third quartiles of the Pearsons' correlation coefficients between the estimates obtained by a given approach and "reference" estimates
fun_pearsons <- function(N){
  ns <- split(N, tabpr$Species)
  n <- split(N_ref, tabpr$Species)
  quantile(na.omit(sapply(1:length(ns), function(i) {
    ifelse(sum(ns[[i]])>0, cor(ns[[i]], n[[i]]), NA)
  })),
  c(0.25,0.5,0.75))
}

## function to calculate the median as well as the first and third quartiles of the Bray-Curtis dissimilarity between the estimates obtained by a given approach and  "reference" estimates
fun_braycurt <- function(N){
  ns <- lapply(split(N, tabpr$Species), function(x) x/sum(x))
  n <- lapply(split(N_ref, tabpr$Species), function(x) x/sum(x))
  quantile(na.omit(sapply(1:length(ns), function(i){
    ifelse(sum(ns[[i]])>0, sum(abs(ns[[i]]-n[[i]]))/sum(ns[[i]]+n[[i]]), NA)
  })),
  c(0.25,0.5,0.75))
}









