#' Perform checks for valid parameters
#'
#' @noRd
.input_checks <- function(x,pvals,z,se,testing,model_type,rendpoint,lendpoint,nclasses,niter_fit,niter_ms,nfits,masking_params,masking_shape,alphas,cr,symmetric_modeling){
  alpha_m <- masking_params$alpha_m
  zeta <- masking_params$zeta
  lambda <- masking_params$lambda
  if(is.null(x) | (is.null(pvals) & is.null(z))){
    stop("Invalid inputs for x, pvals, and test_statistics. None were inputted.")
  }
  if(!is.data.frame(x)){
    stop("Invalid input, x must be a dataframe.")
  }
  if(! (testing %in% c("one_sided","interval","two_sided"))){
    stop("Invalid testing input, valid options are one_sided, interval, two_sided")
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
  if(!is.null(se)){
    if((min(se<=0))){
      stop("Invalid standard errors, standard errors must be at least 0.")
    }
  }
  if(!is.null(se)){
    if(length(se) != nrow(x) & length(se) != 1){
      stop("Invalid standard errors, the standard errors are the wrong length.")
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
  if(niter_fit <= 0 | niter_ms <= 0 | nfits <= 0){
    stop("Invalid number of iterations for niter_fit, niter_ms, or nfits, must be an integer greater than 0.")
  }
  if(min(nclasses) < 2){
    stop("Invalid number of classes, minimum 2 classes.")
  }
  if(!is.null(rendpoint) & !is.null(lendpoint)){
    if(lendpoint > rendpoint){
      stop("Invalid input for endpoints. Right endpoint must be greater than left endpoint.")
    }
  }

  if(symmetric_modeling & testing == "one_sided"){
    stop("One sided testing should not use a symmetric model.")
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
  return(sapply(beta_formulas,complete_pkg))
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
