#' @title Run many regression models quickly, for many response variables.
#'
#' @description Performs any type of cross-sectional regression analysis quickly for many response variables.
#'
#' @param response Character vector of variable names to use as response of interest in separate regression models.
#' @param response_name Label to display for each response of interest
#' @param Xvar Main predictor to be used in each regression.
#' @param dataset Dataset in which these variables exist.
#' @param adj_vars Character vector of adjustment covariate names. 
#' @param family Family for regression model, e.g. 'gaussian', 'binomial', 'poisson'.
#' @param exponentiate Whether or not to exponentiate the estimate from the model (e.g. for logistic or Poisson regression).
#' @param p_adjust_meth Method by which to provide adjustment of P-values. 
#' @param robust If TRUE, returns robust standard error estimates. 
#' @param sort If TRUE, sorts by p-value. 
#' @param decimals Specifies the number of decimals to which to round estimates.
#' @export
#' @importFrom lmtest "coeftest"
#' @importFrom sandwich "vcovHC"

fast_glm2 <- function(response = NULL,
                      response_name = NULL,
                      Xvar = NULL,
                      dataset = NULL,
                      adj_vars=NULL,
                      family = 'gaussian',
                      exponentiate = FALSE,
                      p_adjust_meth = 'fdr',
                      robust = FALSE,
                      sort = TRUE,
                      decimals = 3){
  
  # regression
  cat('Running regression',length(response),'times..\n')
  reg_list <- list()
  for(i in 1:(length(response))){
    if(length(adj_vars) > 1){
      form <- reformulate(c(Xvar, adj_vars),response[i])
      N <- nrow(na.omit(dataset[, c(response[i], Xvar, adj_vars)]))
    }else{
      form <- reformulate(c(Xvar),response[i])
      N <- nrow(na.omit(dataset[, c(response[i], Xvar)]))
    }
    mod <- glm(form, data=dataset, x=FALSE, y=FALSE, family = family)
    if(robust){
      robmod <- lmtest::coeftest(mod, df = Inf,
                                 vcov=sandwich::vcovHC(mod, type='HC0'))
      coef_se <- robmod[2,]
      coef_lci <- coef_se[1]-1.96*coef_se[2]
      coef_uci <- coef_se[1]+1.96*coef_se[2]
      reg_list[[i]] <- c(N, coef_se, coef_lci, coef_uci)
    }
    else reg_list[[i]] <- c(N, summary(mod)$coef[2,],confint.default(mod)[2,])
  }
  # format
  cat('Reformatting output..\n')
  reg_list <- do.call('rbind', reg_list)
  reg_list <- as.data.frame(reg_list)
  if(is.null(response_name)){reg_list$variable <- response[1:(length(response))]}
  else{reg_list$variable <- response_name[1:(length(response_name))]}
  reg_list <- reg_list[,c(8,1,2,6,7,3,4,5)]
  colnames(reg_list)[2:8] <- c('N','estimate','LCI','UCI','SE','t','p')
  if(family == 'binomial'){exponentiate <- TRUE}
  if(exponentiate) {reg_list[,3:5] <- exp(reg_list[,3:5])}
  reg_list$p_adjust <- p.adjust(reg_list$p, method=p_adjust_meth)
  if(sort) reg_list <- reg_list[order(reg_list$p),]
  cat('Variables significantly related after multiple testing correction:',
      length(reg_list$p_adjust[reg_list$p_adjust<0.05]),'\n')
  reg_list[,3:7] <- round(reg_list[,3:7], decimals)
  return(reg_list)
}