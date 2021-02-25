#' @title Short Title of the Function (One Line)
#'
#' @author Karine Princé
#' @references
#' @description estimation of abundance based on opportunistic survey data combined to standardized data
#' @importFrom stats glm predict
#' @param df a data.frame containing original species count data. The data.frame should contain 5 columns with the following headers:
#'  "species": the variable that stores the name/code of the species,
#'  "site": the variable that stores the observation sites,
#'  "dataset": the variable that stores the dataset (e.g. monitoring program) to which the observation data are associated
#'  "effort": a measure of effort that depends on the monitoring program (for example, it can be set to the number of survey points-years implemented in the standardardized program for the site stored in the variable Site)
#'  "N_obs": the observed abundance/count of each *species* detected at each *site*;
#'  @param family string setting the distribution of the error term in the GLM, default='quasipoisson', but can be 'poisson'.
#'  @param newdata an optional data frame in which to look for variables with which to predict and assess the predictive power of the model. Default is NULL. Is newdata is not provided, predicted values are not computed.
#'  @type the type of prediction required. The default is on the scale of the linear predictors; the alternative "response" is on the scale of the response variable. Thus for a quasipoisson poisson model the default predictions are of log-odds (probabilities on log scale) and type = "response" gives the predicted probabilities. The "terms" option returns a matrix giving the fitted values of each term in the model formula on the linear predictor scale.


#' @examples see the vignette("tutorial") for some examples
#'
#' @export
fun_abund_estim <- function(df, family="quasipoisson", newdata=NULL, type = c("link", "response", "terms")){

  ## Check input data prior modeling
  list_var <- c("species","site","N_obs","dataset","effort")
  if(!all(list_var %in% df))
    stop("variables are missing or data frame headers do not match expected format\n")
  if (any(is.na(df$N_obs)))
    stop("Observation variable contain 'NA' values\n")

  ## run model
  cat("\nRun Abundance Model\n=======================\n")
  mod <- stats::glm(N_obs ~ species*dataset + species*site + site*dataset, offset=log(effort),
             data = df, family= family)

  ## look at predictive power of the model
  if(is.null(datapred)){
    N_pred <- NULL
  } else {
    cat("\nCompute Abundance Predictions\n=======================\n")
    N_pred <- stats::predict(mod, newdata = datapred, type = type)
  }

  return(list(mod=mod,N_pred=N_pred))
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









