
#' FLXMRhyreg_het
#'
#' @description Function used in flexmix M-Step to estimate hybrid model accounting for heteroscedastisity
#'
#' @param formula Model formula
#' @param formula_sigma formula for estimation of sigma to account for heteroscedasticity
#' @param data dataset to be used
#' @param family default = "hyreg"
#' @param type a vector containing an indicator wheter that datapoint contains to TTO or DCE Data
#' @param type_cont indicator for continous data
#' @param type_dich indicator for dichotoums data
#' @param variables_both character vactor; variables to be fitted on TTO and DCE data, if not specified all variables from formula are used
#' @param variables_cont character vactor; variables to be fitted only on TTO data
#' @param variables_dich character vactor; variables to be fitted only on DCE data
#' @param stv named vector or list of named vactors containing start values,
#'            has to be a vector if the same start values should be used for all latent classes,
#'            has to be a list of named vectores if different start values are assumed for the latent classes
#'            has to include start values for sigma and theta as well
#' @param stv_sigma vector containing start values for sigma estimates
#' @param offset offset as in flexmix
#' @param optimizer optimizer to be used in bbmle::mle2, default = "optim"
#' @param opt_method optimization method to be used in optimizer, default = "BFGS"
#' @param lower opt_method must be set to "L-BFGS-B", lower bound for censored data
#' @param upper opt_method must be set to "L-BFGS-B", upper bound for censored data
#' @param ...  additional arguments for flexmix or  bbmle::mle2
#'

#'
#'
#' @return a model with hybrid likelihood
#'
#' @author Kim Rand & Svenja Elkenkamp
#' @examples
#'
#' formula <- y ~  -1 + x1 + x2 + x3
#' formula_sigma <- y ~ x1 + x2 + x3
# 'k <- 1
#' stv <- setNames(c(0.2,0,1,1),c(colnames(simulated_data_norm)[3:5],c("theta")))
#' stv_sigma <- setNames(c(0.2,0.2,0.1,1),c(colnames(simulated_data_norm)[3:5],c("(Intercept)")))
#'
# 'rm(counter)
#'
#' x <- model.matrix(formula,simulated_data_norm)
#' y <- simulated_data_norm$y
#' w <- 1

#'model <- FLXMRhyreg_het(formula = formula,
#'                      formula_sigma = formula_sigma,
#'                     family=c("hyreg"),
#'                     type =  simulated_data_norm$type,
#'                     stv = stv,
#'                     stv_sigma = stv_sigma,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = list(iter.max = 1000, verbose = 4),
#'                     offset = NULL,
#'                     optimizer = "optim",
#'                     variables_both =  names(stv)[!is.element(names(stv),c("theta"))],
#'                     variables_cont = NULL,
#'                     variables_dich = NULL,
#'                     lower = -Inf,
#'                     upper = Inf,
#')

#' @importFrom flexmix flexmix
#' @export




### FLXMRhyreg ###

#Error in validObject(.Object) :
#  ungültiges Objekt der Klasse “FLXcomponent”:
# invalid object for slot "df" in class "FLXcomponent": got class "NULL",
# should be or extend class "numeric"



FLXMRhyreg_het <- function(formula= . ~ .,
                           formula_sigma = formula_sigma,
                           data,
                           family=c("hyreg"),
                           type = NULL,
                           type_cont = NULL,
                           type_dich = NULL,
                           variables_both = NULL,
                           variables_cont = NULL,
                           variables_dich = NULL,
                           stv = NULL, # has to include starting value for theta but not sigma
                           stv_sigma = NULL,
                           offset = NULL,
                           opt_method = "BFGS",
                           optimizer = "optim",
                           lower = -Inf,
                           upper = Inf,
                           #                       non_linear = FALSE,
                           #                       formula_orig = formula_orig,
                           ...
)
{
  # rm(counter)
  family <- match.arg(family)

  # refit function has to depend on x,y,w.

  hyregrefit <- function(x, y, w) {
    warning(paste0("Not defined", "Please try xreg2:::refit"))
    return(NA)
  }

  z <- new("FLXMRglm", weighted=TRUE, formula=formula,
           name=paste("FLXMRhyreg"), offset = offset,
           family="hyreg", refit=hyregrefit)

  z@preproc.y <- function(x){
    if (ncol(x) > 1)
      stop(paste("for the", family, "family y must be univariate"))
    x
  }

  if(family=="hyreg"){
    z@defineComponent <- function(para) {  # we get para from the first estimation with start values

      ### NEW ###
      predict <- function(x, type, ...) {
        dotarg = list(...)
        if("offset" %in% names(dotarg)) offset <- dotarg$offset

        if(type == type_cont){
          # change for non-linear functions
          p <- x %*% para$coef[is.element(names(para$coef),c(variables_cont,variables_both))]  # Xb in xreg
        }
        if(type == type_dich){
          # change for non-linear functions
          p <- (x %*% para$coef[is.element(names(para$coef),c(variables_dich,variables_both))]) * para$theta
        }

        if (!is.null(offset)) p <-  p + offset
        p
      }

      #
      #
      logLik <- function(x, y, sigma = para$sigma, theta = para$theta, return_vector = TRUE, ...){
        #para$sigma must be a vector

        #prepare sigma
        xsigma <-  model.matrix(formula_sigma,data[type == type_cont,])
        colnames(xsigma) <- paste0(colnames(xsigma),"_h")
        sigma <- exp(xsigma %*% sigma[colnames(xsigma)])


        # prepare data
        # change for non-linear functions
        # use formula_orig
        x1 <- x[type == type_cont,c(variables_cont,variables_both)]
        x2 <-  x[type == type_dich,c(variables_dich,variables_both)]
        y1 <- y[type == type_cont]
        y2 <-  y[type == type_dich]


        # linear predictor
        # change for non-linear functions
        Xb1 <- x1 %*% para$coef[colnames(x1)] # only cont and both variables
        Xb2 <- (x2 %*% para$coef[colnames(x2)]) * exp(theta)  # only dich and both variables


        # pvals and likelihood
        logistic_tmp <- .5+.5*tanh(Xb2/2)
        pvals2 <- log(y2 *logistic_tmp + (1-y2)* (1-logistic_tmp)) # dich_logistic

        pvals1 <- dnorm(y1, mean=Xb1, sd=sigma, log=TRUE) # cont_normal


        # for box constraints
        #in xreg:
        if(upper != Inf ){
          censV <- y1 == upper
          pvals1[which(censV)] <- log(pnorm((Xb1[censV]-upper)/sigma[censV],0,1) - pnorm((Xb1[censV]-Inf)/sigma,0,1))
        }

        if(lower != -Inf ){
          censV <- y1 == lower
          pvals1[which(censV)] <- log(pnorm((Xb1[censV]-(-Inf))/sigma[censV],0,1) - pnorm((Xb1[censV]-lower)/sigma,0,1))
        }


        pvals <- c(pvals1,pvals2)

        pvals[pvals == -Inf] <- log(.Machine$double.xmin)
        pvals[pvals == Inf] <- log(.Machine$double.xmax)

        if(return_vector == TRUE){
          return(pvals)
        }else{
          return(-sum(pvals))         # get neg log L for optimizer
        }

      }

      new("FLXcomponent",
          parameters=list(coef=para$coef,
                          sigma=para$sigma,
                          theta = para$theta,
                          #  stderror = para$stderror,
                          #  pvalue = para$pvalue,
                          fit_mle = para$fit_mle),
          # counter = counter),
          #minLik = para$minLik),
          logLik=logLik, predict=predict,
          df=para$df)
    }


    z@fit <- function(x, y, w, component, ...){


      # function to use in mle, same as logLik but depending on stv and giving out the neg logL directly
      logLik2 <- function(stv){

        # prepare sigma
        stv_sigma <- stv[is.element(names(stv),names(stv_sigma))]
        stv <-  stv[!is.element(names(stv),names(stv_sigma))]

        xsigma <-  model.matrix(formula_sigma,data[type == type_cont,])
        colnames(xsigma) <- paste0(colnames(xsigma),"_h")
        sigma <- exp(xsigma %*% stv_sigma[colnames(xsigma)])

        # prepare data
        x1 <- x[type == type_cont,c(variables_cont,variables_both)]
        x2 <-  x[type == type_dich,c(variables_dich,variables_both)]
        y1 <- y[type == type_cont]
        y2 <-  y[type == type_dich]


        #prepare stv
        theta <- exp(stv[is.element(names(stv),c("theta"))][[1]])
        stv_cont <- stv[!is.element(names(stv),c("sigma","theta", variables_dich))]
        stv_dich <- stv[!is.element(names(stv),c("sigma","theta", variables_cont))]


        # linear predictors
        # use formula_orig for non-linear functions
        Xb1 <- x1 %*% stv_cont[colnames(x1)] # hier könnte man ggf nur TTO spezifische Variablen einfließen lassen, Interaktionen etc beachten
        Xb2 <- x2 %*% stv_dich[colnames(x2)]
        Xb2 <- Xb2*theta


        # pvals and Likelihood calculation
        logistic_tmp <- 0.5 + 0.5*tanh(Xb2/2)
        pvals2 <- log(y2 *logistic_tmp + (1-y2)* (1-logistic_tmp)) # dich_logistic

        pvals1 <- dnorm(y1, mean=Xb1, sd=sigma, log=TRUE) # cont_normal

        # for box constraints:
        #in xreg:
        if(upper != Inf ){
          censV <- y1 == upper
          pvals1[which(censV)] <- log(pnorm((Xb1[censV]-upper)/sigma[censV],0,1) - pnorm((Xb1[censV]-Inf)/sigma,0,1))
        }

        if(lower != -Inf ){
          censV <- y1 == lower
          pvals1[which(censV)] <- log(pnorm((Xb1[censV]-(-Inf))/sigma[censV],0,1) - pnorm((Xb1[censV]-lower)/sigma,0,1))
        }


        pvals <- c(pvals1,pvals2)


        pvals[pvals == -Inf] <- log(.Machine$double.xmin) # or without log?
        pvals[pvals == Inf] <- log(.Machine$double.xmax)


        # use weights for EM algo
        pvals_w <- pvals * w


        return(-sum(pvals_w))   # look up Leisch: FlexMix: A General Framework for Finite Mixture
        #                  Models and Latent Class Regression in R, p. 3
        # w must be posterior class probabilities for each observation, pvals is already the log?
      }


      # changes names of stv_sigma to be able to get estimates for it
      names(stv_sigma) <- paste0(names(stv_sigma),"_h")
      bbmle::parnames(logLik2) <- c(colnames(x),"theta",names(stv_sigma)) # set names of inputs for logLik2



      ### MODEL ESTIMATION ###
      if(!exists("counter")){
        counter <<- 1


        # use different stv for different components
        # implement stv as lits? and ask if it is a list, then
        # use stv[[counter]] as stv
        # for the next iterations of EM its not requried since we use component$coef

        stv_in <- stv # without sigma !
        stv_sigma_in <- stv_sigma

        if(is.list(stv)){
          stv_in <- stv[[counter]]
        }
        if(is.list(stv_sigma)){
          stv_sigma_in <- stv_sigma[[counter]]
        }


        fit_mle <- bbmle::mle2(minuslogl = logLik2,
                               start = c(stv_in,stv_sigma_in),
                               optimizer = optimizer,
                               method = opt_method,
                               lower = lower,
                               upper = upper)

      }else{
        if(counter < k){
          counter <<- counter + 1


          stv_in <- stv # without sigma !
          stv_sigma_in <- stv_sigma

          if(is.list(stv)){
            stv_in <- stv[[counter]]
          }
          if(is.list(stv_sigma)){
            stv_sigma_in <- stv_sigma[[counter]]
          }


          fit_mle <- bbmle::mle2(minuslogl = logLik2,
                                 start = c(stv,stv_sigma),
                                 optimizer = optimizer,
                                 method = opt_method,
                                 lower = lower,
                                 upper = upper)


        }else{
          stv_new <- setNames(c(component$coef,component$theta,component$sigma),c(colnames(x),"theta",names(stv_sigma)))
          fit_mle <- bbmle::mle2(minuslogl = logLik2,
                                 start = stv_new,
                                 optimizer = optimizer,
                                 method = opt_method,
                                 lower = lower,
                                 upper = upper)
        }
      }

      z@defineComponent(para = list(coef = fit_mle@coef[!is.element(names(fit_mle@coef),c(names(stv_sigma),"theta"))],
                                    df = ncol(x) + 1, # wrong? how to change
                                    sigma = fit_mle@coef[is.element(names(fit_mle@coef),names(stv_sigma))],
                                    theta = fit_mle@coef[is.element(names(fit_mle@coef),c("theta"))],
                                    fit_mle = fit_mle,
                                    # counter = counter,
                                    # stderror = summary(fit_mle)@coef[,2],  # trying to get slot "coef" from an object (class "summaryDefault") that is not an S4 object
                                    #  pvalue = summary(fit_mle)@coef[,4],
                                    minLik = fit_mle@min)
      )
    }
  }

  z
}
