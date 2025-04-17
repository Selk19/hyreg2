
# Write R package
# https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

#' FLXMRhyreg
#'
#' @description Function used in flexmix M-Step to estimate hybrid model
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
#' @param offset offset as in flexmix
#' @param optimizer optimizer to be used in bbmle::mle2, default = "optim"
#' @param opt_method optimization method to be used in optimizer, default = "BFGS"
#' @param lower opt_method must be set to "L-BFGS-B", lower bound for censored data
#' @param upper opt_method must be set to "L-BFGS-B", upper bound for censored data
#' @param ... additional arguments for flexmix or bbmle::mle2
#'
#' @return a model with hybrid likelihood
#'
#' @author Kim Rand & Svenja Elkenkamp
#' @examples Put Example here

#' @importFrom flexmix flexmix
#' @importFrom  bbmle mle2
#' @importFrom  bbmle summary
#' @export




### FLXMRhyreg ###
# not completly working (refit is missing)


#### WITH SIGNA AND THETA INCLUDED IN ESTIMATION ###

FLXMRhyreg <- function(formula= .~. ,
                       family=c("hyreg"),
                       type = NULL,
                       type_cont = NULL,
                       type_dich = NULL,
                       variables_both = NULL,
                       variables_cont = NULL,
                       variables_dich = NULL,
                       stv = NULL, # has to include starting values for sigma and theta as well
                       offset = NULL,
                       opt_method = "BFGS",
                       optimizer = "optim",
                       lower = -Inf,
                       upper = Inf,
                       ...
)
{
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
          p <- x %*% para$coef[is.element(names(para$coef),c(variables_cont,variables_both))]  # Xb in xreg
        }
        if(type == type_dich){
          p <- (x %*% para$coef[is.element(names(para$coef),c(variables_dich,variables_both))]) * para$theta
        }

        if (!is.null(offset)) p <-  p + offset
        p
      }



      logLik <- function(x, y, sigma = para$sigma, theta = para$theta, return_vector = TRUE, ...){


        # choose subset of x and y depending on type
        x1 <- x[type == type_cont,c(variables_cont,variables_both)]
        x2 <-  x[type == type_dich,c(variables_dich,variables_both)]
        y1 <- y[type == type_cont]
        y2 <-  y[type == type_dich]


        sigma <- exp(sigma)

        Xb1 <- x1 %*% para$coef[colnames(x1)] # only cont and both variables
        Xb2 <- (x2 %*% para$coef[colnames(x2)]) * exp(theta)  # only dich and both variables


        logistic_tmp <- .5+.5*tanh(Xb2/2)
        pvals2 <- log(y2 *logistic_tmp + (1-y2)* (1-logistic_tmp)) # dich_logistic


        # Likelihood calculation
        pvals1 <- dnorm(y1, mean=Xb1, sd=sigma, log=TRUE) # cont_normal



        ### from xreg ###
        ### for box constraints ?

        if(upper != Inf){
          censV <- y1 == upper # adadpt for lower
          pvals1[which(censV)] <- pnorm(q = Xb1[censV], mean = upper, sd = sigma, lower.tail = T, log.p = T) -  # mean = 2 in xreg
            pnorm(q =  Xb1[censV], mean = lower, sd = sigma, lower.tail = T, log.p = T) # mean = Inf in xreg
          # or do we have to use log.p = F and than log(pnorm() - pnorm())?
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
        #  counter = counter), # Error: unzulässiger Name für Slot der Klasse “FLXcomponent”; fit_mle
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



        sigma <- exp(stv[is.element(names(stv),c("sigma"))][[1]])
        theta <- exp(stv[is.element(names(stv),c("theta"))][[1]])
        stv_cont <- stv[!is.element(names(stv),c("sigma","theta", variables_dich))]
        stv_dich <- stv[!is.element(names(stv),c("sigma","theta", variables_cont))]


        Xb1 <- x1 %*% stv_cont[colnames(x1)] # hier könnte man ggf nur TTO spezifische Variablen einfließen lassen, Interaktionen etc beachten
        Xb2 <- x2 %*% stv_dich[colnames(x2)]


        Xb2 <- Xb2*theta

        logistic_tmp <- 0.5 + 0.5*tanh(Xb2/2)
        pvals2 <- log(y2 *logistic_tmp + (1-y2)* (1-logistic_tmp)) # dich_logistic

        # Likelihood calculation
        pvals1 <- dnorm(y1, mean=Xb1, sd=sigma, log=TRUE) # cont_normal

        # for box constraints:

        if(upper != Inf){
          censV <- y1 == upper # adapt for lower, anpassen für lower bzw upper und lower zeitgleich
          pvals1[which(censV)] <- pnorm(q = Xb1[censV], mean = upper, sd = sigma, lower.tail = T, log.p = T) -  # mean = 2 in xreg
            pnorm(q =  Xb1[censV], mean = lower, sd = sigma, lower.tail = T, log.p = T) # mean = Inf in xreg
          # or do we have to use log.p = F and than log(pnorm() - pnorm())?

        }



        pvals <- c(pvals1,pvals2)


        pvals[pvals == -Inf] <- log(.Machine$double.xmin) # or without log?
        pvals[pvals == Inf] <- log(.Machine$double.xmax)

        # what to do with NaN


        pvals_w <- pvals * w


        return(-sum(pvals_w))   # look up Leisch: FlexMix: A General Framework for Finite Mixture
        #                  Models and Latent Class Regression in R, p. 3
        # w must be posterior class probabilities for each observation, pvals is already the log?
      }


      bbmle::parnames(logLik2) <- c(colnames(x),"sigma","theta") # set names of inputs for logLik2

      if(!exists("counter")){
        # if(iter == 1){
        counter <<- 1

        # use different stv for different components
        # implement stv as lits? and ask if it is a list, then
        # use stv[[counter]] as stv
        # for the next iterations of EM its not requried since we use component$coef

        if(class(stv) == "list"){
          stv <- stv[[counter]]
        }

        fit_mle <- bbmle::mle2(minuslogl = logLik2,
                               start = stv,
                               optimizer = optimizer,
                               method = opt_method,
                               lower = lower,
                               upper = upper)

      }else{
        if(counter < k){
          counter <<- counter + 1

          if(class(stv) == "list"){
            stv <- stv[[counter]]
          }

          fit_mle <- bbmle::mle2(minuslogl = logLik2,
                                 start = stv,
                                 optimizer = optimizer,
                                 method = opt_method,
                                 lower = lower,
                                 upper = upper)



        }else{
          stv_new <- setNames(c(component$coef,component$sigma,component$theta),c(colnames(x),"sigma","theta"))
          fit_mle <- bbmle::mle2(minuslogl = logLik2,
                                 start = stv_new,
                                 optimizer = optimizer,
                                 method = opt_method,
                                 lower = lower,
                                 upper = upper)
        }
      }


      z@defineComponent(para = list(coef = fit_mle@coef[!is.element(names(fit_mle@coef),c("sigma","theta"))],
                                    df = ncol(x)+1, # not changed yet
                                    sigma = fit_mle@coef[is.element(names(fit_mle@coef),c("sigma"))],
                                    theta = fit_mle@coef[is.element(names(fit_mle@coef),c("theta"))],
                                    fit_mle = fit_mle,
                                  #  counter = counter,
                                   # stderror = summary(fit_mle)@coef[,2],  # trying to get slot "coef" from an object (class "summaryDefault") that is not an S4 object
                                  #  pvalue = summary(fit_mle)@coef[,4],
                                    minLik = fit_mle@min)
      )
    }
  }

  z
}
