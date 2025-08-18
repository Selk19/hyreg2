
#' FLXMRhyreg_het
#'
#' @description Function used in flexmix M-Step to estimate hybrid model accounting for heteroscedastisity
#'
#' @param formula Model formula
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
#' @examples Put Example here

#' @importFrom flexmix flexmix
#' @export




### FLXMRhyreg ###
# not completly working (refit is missing)


# at the moment we get:
# Error in postunscaled + FLXdeterminePostunscaled(model[[n]], components[[n]]) :
# non-conformable arrays
# change flexmix:::FLXdeterminePostunscaled ? and upload an adpted version?


#### WITH Formula for sigma to account for heteroscedastisity ###
# formula for sigma can only be the same as in formula estimated on only the TTO values
# idea: use additional formula_sigma argument to be more flexible
# stv_sigma must be named vector
# names(stv_sigma) have to be same variable names as in formula/data but with _h at the ending

FLXMRhyreg_het <- function(formula= . ~ .,
                        #   formula_sigma = NULL,
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
                       non_linear = FALSE,
                       formula_orig = formula_orig,
                       ...
)
{
 # rm(counter)
  family <- match.arg(family)

  # refit function has to depend on x,y,w.

  hyregrefit <- function(x, y, w) {
    # use mle2?
    fit <- bbmle::mle2() ## or use optim?

    # in FLXMRglm they use:
    #c(glm.fit(x, y, weights=w, offset=offset,
    #         family=get(family, mode="function")()),
    # list(call = sys.call(), offset = offset,
    #      control = eval(formals(glm.fit)$control),
    #      method = "weighted.glm.fit"))
    #  fit$df.null <- sum(w) + fit$df.null - fit$df.residual - fit$rank
    #  fit$df.residual <- sum(w) - fit$rank
    #  fit$x <- x
    #  fit
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


        # choose subset of x and y depending on type
        # at the moment sigma formula can use only same dependet variables as formula itself

        # change for non-linear functions
        # use formula_orig
        x1 <- x[type == type_cont,c(variables_cont,variables_both)]
        x2 <-  x[type == type_dich,c(variables_dich,variables_both)]
        y1 <- y[type == type_cont]
        y2 <-  y[type == type_dich]


        # if sigma uses Intercept include 1s
        # intercept must be given on first position
        # CHECK:
        # hat stv_sigma Namen am Ende ein _h? Wenn Nein, dann Error Meldung entsprechend
        if(!is.element("(Intercept)",colnames(x)) & is.element("(Intercept)_h",names(stv_sigma))){
          sigma <- exp(cbind(rep(1,dim(x1)[1]),x1) %*% sigma) # sigma must be a vector now
        }else{
          sigma <- exp( x1[,unlist(strsplit(names(stv_sigma),"_h"))] %*% sigma)
        }


        # change for non-linear functions
        Xb1 <- x1 %*% para$coef[colnames(x1)] # only cont and both variables
        Xb2 <- (x2 %*% para$coef[colnames(x2)]) * exp(theta)  # only dich and both variables


        logistic_tmp <- .5+.5*tanh(Xb2/2)
        pvals2 <- log(y2 *logistic_tmp + (1-y2)* (1-logistic_tmp)) # dich_logistic


        # Likelihood calculation
        pvals1 <- dnorm(y1, mean=Xb1, sd=sigma, log=TRUE) # cont_normal



        ### from xreg ###
        ### for box constraints


        #in xreg:
        if(upper != Inf ){
          censV <- y1 == upper
          pvals1[which(censV)] <- log(pnorm((Xb1[censV]-upper)/sigma[censV],0,1) - pnorm((Xb1[censV]-Inf)/sigma,0,1))
        }

        if(lower != -Inf ){
          censV <- y1 == lower
          pvals1[which(censV)] <- log(pnorm((Xb1[censV]-(-Inf))/sigma[censV],0,1) - pnorm((Xb1[censV]-lower)/sigma,0,1))
        }

        # if(upper != Inf){
        #   censV <- y1 == upper # adapt for lower
        #   pvals1[which(censV)] <- pnorm(q = Xb1[censV], mean = upper, sd = sigma[censV], lower.tail = T, log.p = T) -  # mean = 2 in xreg
        #     pnorm(q =  Xb1[censV], mean = lower, sd = sigma[censV], lower.tail = T, log.p = T) # mean = Inf in xreg
        #   # or do we have to use log.p = F and than log(pnorm() - pnorm())?
        # }


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
        # variables_cont, variables_both, variables_dich
        # as charachter, names of variables to be fitted for only specific type of data


        x1 <- x[type == type_cont,c(variables_cont,variables_both)]
        x2 <-  x[type == type_dich,c(variables_dich,variables_both)]
        y1 <- y[type == type_cont]
        y2 <-  y[type == type_dich]



        # include vector of 1s for intercept ?
        # if Intercept in stv_sigma but not in stv, than include ones
        # names(stv_sigma) have to be same variable names as in formula but with _h at the ending


        # CHECK:
        # hat stv_sigma Namen am Ende ein _h? Wenn Nein, dann Error Meldung entsprechend
        if(!is.element("(Intercept)",colnames(x)) & is.element("(Intercept)_h",names(stv_sigma))){
          sigma <- exp( as.matrix(cbind(rep(1,dim(x1)[1]),x1)) %*% stv[is.element(names(stv),names(stv_sigma))]) # exp like in xreg?
        }else{
          sigma <- exp( x1[,unlist(strsplit(names(stv_sigma),"_h"))] %*% stv[is.element(names(stv),names(stv_sigma))])
        }

        theta <- exp(stv[is.element(names(stv),c("theta"))][[1]])
        stv_cont <- stv[!is.element(names(stv),c("sigma","theta", variables_dich,names(stv_sigma)))]
        stv_dich <- stv[!is.element(names(stv),c("sigma","theta", variables_cont,names(stv_sigma)))]

        # use formula_orig for non-linear functions
        Xb1 <- x1 %*% stv_cont[colnames(x1)] # hier könnte man ggf nur TTO spezifische Variablen einfließen lassen, Interaktionen etc beachten
        Xb2 <- x2 %*% stv_dich[colnames(x2)]


        Xb2 <- Xb2*theta

        logistic_tmp <- 0.5 + 0.5*tanh(Xb2/2)
        pvals2 <- log(y2 *logistic_tmp + (1-y2)* (1-logistic_tmp)) # dich_logistic

        # Likelihood calculation
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

        # if(upper != Inf){
        #   censV <- y1 == upper # adapt for lower
        #   pvals1[which(censV)] <- pnorm(q = Xb1[censV], mean = upper, sd = sigma[censV], lower.tail = T, log.p = T) -  # mean = 2 in xreg
        #     pnorm(q =  Xb1[censV], mean = lower, sd = sigma[censV], lower.tail = T, log.p = T) # mean = Inf in xreg
        #   # or do we have to use log.p = F and than log(pnorm() - pnorm())?
        # }



        pvals <- c(pvals1,pvals2)


        pvals[pvals == -Inf] <- log(.Machine$double.xmin) # or without log?
        pvals[pvals == Inf] <- log(.Machine$double.xmax)

        # what to do with NaN


        pvals_w <- pvals * w


        return(-sum(pvals_w))   # look up Leisch: FlexMix: A General Framework for Finite Mixture
        #                  Models and Latent Class Regression in R, p. 3
        # w must be posterior class probabilities for each observation, pvals is already the log?
      }


       bbmle::parnames(logLik2) <- c(colnames(x),"theta",names(stv_sigma)) # set names of inputs for logLik2

       if(!exists("counter")){
         counter <<- 1

         # use different stv for different components
         # implement stv as lits? and ask if it is a list, then
         # use stv[[counter]] as stv
         # for the next iterations of EM its not requried since we use component$coef

         if(class(stv) == "list"){
           stv <- stv[[counter]]
           stv_sigma <- stv_sigma[[counter]]
         }
         fit_mle <- bbmle::mle2(minuslogl = logLik2,
                                start = c(stv,stv_sigma),
                                optimizer = optimizer,
                                method = opt_method,
                                lower = lower,
                                upper = upper)

       }else{
         if(counter < k){
           counter <<- counter + 1

           if(class(stv) == "list"){
             stv <- stv[[counter]]
             stv_sigma <- stv_sigma[[counter]]
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
                                    df = ncol(x)+1, # not changed yet
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
