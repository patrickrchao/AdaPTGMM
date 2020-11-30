#' Weighted Mean computation
#'
#' @param values Vector of values
#' @param weights Vector of nonnegative weights
#'
#' @return weighted mean
#' @noRd
.weighted_mean <- function(values,weights){
  return(sum(values*weights)/sum(weights))
}

#' Perform checks for valid parameters
#'
#' @noRd
.input_checks <- function(x,pvals,z,testing,model_type,rendpoint,lendpoint,nclasses,niter_fit,niter_ms,nfit,masking_params,masking_shape,alphas,cr){
  alpha_m <- masking_params$alpha_m
  zeta <- masking_params$zeta
  lambda <- masking_params$lambda
  if(is.null(x) | (is.null(pvals) & is.null(z))){
    stop("Invalid inputs for x, pvals, and test_statistics. None were inputted.")
  }
  if(!is.data.frame(x)){
    stop("Invalid input, x must be a dataframe.")
  }
  if(testing == "interval"){
    if(is.null(z)){
      stop("Invalid input for z. Must include test statistics to perform interval null testing.")
    }else if(!is.null(pvals)){
      warning("Inputted p-values will be ignored in interval null testing, they will be computed automatically.")
    }
    if(is.null(rendpoint)){
      stop("Invalid input for rendpoint. Must include right endpoint for interval null testing.")
    }
  }
  #if(length(x) != max(length(pvals), length(z))){
  #  stop("Invalid inputs, x and pvals/test statistics have different length.")
  #}
  if(!is.null(pvals)){
     if((min(pvals<0) | max(pvals) > 1)){
        stop("Invalid p-values, p-values must be in range [0,1].")
     }
  }

  if(alpha_m<0 | alpha_m>1 | alpha_m>lambda | lambda<0 | lambda>1 | zeta<0  | lambda+alpha_m*zeta>1){
    stop("Invalid input for alpha_m and lambda, must all be between 0 and 1 and 0<alpha_m<=lambda<lambda+alpha_m*zeta<=1.")
  }
  if(masking_shape!="tent" & masking_shape != "comb"){
    stop("Invalid masking shape inputted, must be `tent` or `comb`.")
  }
  if(min(alphas) < 0 | max(alphas) > 1){
    stop("Invalid alphas inputted, alphas must be in range [0,1].")
  }
  if(niter_fit <= 0 | niter_ms <= 0 | nfit <= 0){
    stop("Invalid number of iterations for niter_fit, niter_ms, or nfit, must be an integer greater than 0.")
  }
  if(min(nclasses) < 2){
    stop("Invalid number of classes, minimum 2 classes.")
  }
  if(!is.null(rendpoint) & !is.null(lendpoint)){
    if(lendpoint > rendpoint){
      stop("Invalid input for endpoints. Right endpoint must be greater than left endpoint.")
    }
  }

  #if(!(cr  %in% c("AIC","BIC","HIC","AICC","spread"))){
  #  stop("Invalid criterion inputted. Valid criterion are AIC, BIC, HIC, AICC, spread")
  #}
}

#' Helper function to find inverse of a function
#' Used to find test statistics that correspond to mirrored p-values
#' From Mike Axiak, https://stackoverflow.com/questions/10081479/solving-for-the-inverse-of-a-function-in-r
#' @param f function
#'
#' @noRd
.inverse = function (f, lower = 0, upper = 50) {
  return(function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper, tol=.Machine$double.eps^2)[1])
}

clean_beta_formulas <- function(beta_formulas,include_intercept_model){
  if(include_intercept_model & !("+1" %in% beta_formulas)){
    beta_formulas <- c("+1",beta_formulas)
  }
  return(unlist(lapply(beta_formulas,complete_pkg)))
}





# From Lihua Lei's AdaPTMT Package
# https://github.com/lihualei71/adaptMT/blob/master/R/helpers.R
complete_pkg <- function(formula){
  formula <- as.character(formula)
  formula <- tail(formula, 1)
  formula <- tail(strsplit(formula, "~")[[1]], 1)
  formula <- paste0(" ", formula)
  if (grepl("ns\\(", formula)){
    if (!requireNamespace("splines", quietly = TRUE)){
      stop("package \'splines\' not found. Please install.")
    }
    formula <- gsub("([^:])ns\\(", "\\1splines::ns\\(", formula)
  }
  if (grepl("[^a-z]s\\(", formula)){
    if (!requireNamespace("mgcv", quietly = TRUE)){
      stop("package \'mgcv\' not found. Please install.")
    }
  }
  if (grepl("([^:a-z])mgcv::s\\(", formula)){
    formula <- gsub( "([^:a-z])mgcv::s\\(","\\1s\\(", formula)
  }
  formula <- paste0("class  ~ ",formula)
  formula <- as.formula(formula)
  return(formula)
}

set_default_target <- function(target_alpha_level,alphas,default_value=0.05){
  if(is.null(target_alpha_level)){
    if(length(alphas) == 1){
      target_alpha_level = alphas
    }else{
      target_alpha_level = default_value
    }
  }
  return(target_alpha_level)
}

select_initialization <- function(masking_params, initialization){
  if(masking_params$zeta == 1 &
     ((0.5 - masking_params$alpha_m ) == (masking_params$lambda - 0.5)) &
     initialization == "kmeans"){

    warning("Symmetric masking function found, possible unstable performance with k-means. Setting initialization scheme to `random`.")
    initialization = "random"
  }
  return(initialization)
}

#' #' Helper function to ensure all formulas are valid
#' #'
#' #' @param x dataframe of covariates
#' #' @param beta_formulas list of all beta formulas
#' .check_formulas <- function(x,beta_formulas){
#'   new_formulas <- lapply(beta_formulas,function(formula){
#'   tryCatch({
#'         .evaluate_formula(x,formula)
#'         return(formula)
#'       }, error = function(e) {
#'         warning(paste0("Invalid beta_formula found: ",formula,". Ignoring formula."))
#'         return("invalid")
#'       })
#'   })
#'
#'   new_formulas <- unlist(new_formulas[new_formulas!="invalid"])
#'   if(length(new_formulas)==0){
#'     stop("All formulas are invalid. Stopping.")
#'   }else if(length(new_formulas) == 1 & new_formulas[1] == "+1"){
#'     warning("Only using intercept model.")
#'   }
#'   return(new_formulas)
#' }
#'
#' .evaluate_formula <- function(x,formula){
#'   new_x <- eval(parse(text=formula),x)
#' }
