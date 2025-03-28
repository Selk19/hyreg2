
# Write R package
# https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

#' FLXMRhyreg
#'
#' @description Function used in flexmix M-Step to estimate hybrid model
#'
#' @param formula: Model formula
#'
#' @return a model with hybrid likelihood
#'
#' @author Kim Rand & Svenja Elkenkamp
#' @examples Put Example here


#' @export


### FLXMRhyreg ###
# not completly working (refit is missing)


#### WITH SIGNA AND THETA INCLUDED IN ESTIMATION ###

FLXMRhyreg <- function(formula= .~. ,
                       family=c("hyreg"),
                       type = NULL,
                       type_cont = NULL,
                       type_dich = NULL,
                       stv = NULL, # has to include starting values for sigma and theta as well
                       offset = NULL,
                       opt_method = "BFGS",
                       optimizer = "optim",
                       lower = -Inf,
                       upper = Inf
)
{
  family <- match.arg(family)

  # fit function has to depend on x,y,w.
  # actual Problem: we need logLik2 but it can´t be found from R
  # The fit() function returns an object of class "FLXcomponent"

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
    #    if (ncol(x) > 1)
    #      stop(paste("for the", family, "family y must be univariate"))
    x
  }

  if(family=="hyreg"){
    z@defineComponent <- function(para) {  # we get para from the first estimation with start values

      ### NEW ###
      predict <- function(x, type, ...) {
        dotarg = list(...)
        if("offset" %in% names(dotarg)) offset <- dotarg$offset

        if(type == type_cont){
          p <- x %*% para$coef  # Xb in xreg
        }
        if(type == type_dich){
          p <- (x %*% para$coef) * para$theta
        }

        if (!is.null(offset)) p <-  p + offset
        p
      }

      #
      #
      logLik <- function(x, y, sigma = para$sigma, theta = para$theta, return_vector = TRUE, ...){


        # choose subset of x and y depending on type
        # use predict?
        x1 <- x[type == type_cont,]
        x2 <-  x[type == type_dich,]
        y1 <- y[type == type_cont]
        y2 <-  y[type == type_dich]

        sigma <- exp(sigma)

        Xb1 <- x1 %*% para$coef
        Xb2 <- (x2 %*% para$coef) * exp(theta)


        logistic_tmp <- .5+.5*tanh(Xb2/2)
        pvals2 <- log(y2 *logistic_tmp + (1-y2)* (1-logistic_tmp)) # dich_logistic


        # Likelihood calculation
        pvals1 <- dnorm(y1, mean=Xb1, sd=sigma, log=TRUE) # cont_normal



        ### from xreg ###
        ### for box constraints ?

        if(upper != Inf){
          censV <- y1 == upper
          pvals1[which(censV)] <- pnorm(q = Xb1[censV], mean = upper, sd = sigma, lower.tail = T, log.p = T) -  # mean = 2 in xreg
            pnorm(q =  Xb1[censV], mean = lower, sd = sigma, lower.tail = T, log.p = T) # mean = Inf in xreg
          # or do we have to use log.p = F and than log(pnorm() - pnorm())?
        }


        pvals <- c(pvals1,pvals2)

        pvals[pvals == -Inf] <- log(.Machine$double.xmin) # or without log?
        pvals[pvals == Inf] <- log(.Machine$double.xmax) # or without log?

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
                          stderror = para$stderror,
                          pvalue = para$pvalue),
          #minLik = para$minLik),
          logLik=logLik, predict=predict,
          df=para$df)
    }


    z@fit <- function(x, y, w, component, ...){


      # function to use in mle, same as logLik but depending on stv and giving out the neg logL directly
      logLik2 <- function(stv){

        x1 <- x[type == type_cont,]
        x2 <-  x[type == type_dich,]
        y1 <- y[type == type_cont]
        y2 <-  y[type == type_dich]



        sigma <- exp(stv[is.element(names(stv),c("sigma"))][[1]]) # exp like in xreg?
        theta <- exp(stv[is.element(names(stv),c("theta"))][[1]]) # exp like in xreg?
        stv <- stv[!is.element(names(stv),c("sigma","theta"))]

        Xb1 <- x1 %*% stv
        Xb2 <- x2 %*% stv


        Xb2 <- Xb2*theta

        logistic_tmp <- 0.5 + 0.5*tanh(Xb2/2)
        pvals2 <- log(y2 *logistic_tmp + (1-y2)* (1-logistic_tmp)) # dich_logistic

        # Likelihood calculation
        pvals1 <- dnorm(y1, mean=Xb1, sd=sigma, log=TRUE) # cont_normal

        # for box constraints ?

        ### from xreg ###
        ### for box constraints ?

        if(upper != Inf){
          censV <- y1 == upper # anpassen für lower bzw upper und lower zeitgleich
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

      if(iter == 1){ # maybe we can ask if the object components already has elements and use them, maybe use k in FLEXmstep for different start values in different classes
        fit <- bbmle::mle2(minuslogl = logLik2,
                           start = stv,
                           optimizer = optimizer,
                           method = opt_method,
                           lower = lower,
                           upper = upper)
      }else{
        stv_new <- setNames(c(component$coef,component$sigma,component$theta),c(colnames(x),"sigma","theta"))
        fit <- bbmle::mle2(minuslogl = logLik2,
                           start = stv_new,
                           optimizer = optimizer,
                           method = opt_method,
                           lower = lower,
                           upper = upper)
      }

      # z@definecomponent(para = list(coef = fit@coef[!is.element(names(fit@coef),c("sigma","theta"))],
      #                               df = ncol(x)+1, # not changed yet
      #                               sigma =  fit@coef[is.element(names(fit@coef),c("sigma"))], # does this has to be calculted outside mle2?
      #                               theta = fit@coef[is.element(names(fit@coef),c("theta"))])) # does this has to be calculted outside mle2?

      z@defineComponent(para = list(coef = fit@coef[!is.element(names(fit@coef),c("sigma","theta"))],
                                    df = ncol(x)+1, # not changed yet
                                    sigma = fit@coef[is.element(names(fit@coef),c("sigma"))],
                                    theta = fit@coef[is.element(names(fit@coef),c("theta"))],
                                    stderror = summary(fit)@coef[,2],
                                    pvalue = summary(fit)@coef[,4],
                                    minLik = fit@min)
      )
    }
  }

  z
}

### REFIT OUTSIDE THE M STEP DRIVER ###
# construct output with standard errors and p values (like refit should normaly do)
# can we include this in the m step driver? change method for refit?

gendf <- function(list){
  df <- list()
  for(j in 1:length(list)){
    df[[j]] <- data.frame("Estimates" = c(list[[j]]@parameters[["coef"]], #  maybe include an other apply function here?
                                          list[[j]]@parameters[["sigma"]],
                                          list[[j]]@parameters[["theta"]]),
                          "Std Error" = list[[j]]@parameters[["stderror"]],
                          "pvalue" = list[[j]]@parameters[["pvalue"]])
    # give AIC, loglik here as additional arguments?
  }

  return(df)
}


refit <- function(object){
  comp <- object@components
  out <- lapply(comp, function(j) {gendf(j)})
  return(out)
}


devtools::document()

