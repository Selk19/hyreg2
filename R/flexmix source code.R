
require(flexmix)

### Code from flexmix ###

RemoveGrouping <- function(formula) {
  lf <- length(formula)
  formula1 <- formula
  if(length(formula[[lf]])>1) {
    if (deparse(formula[[lf]][[1]]) == "|"){
      formula1[[lf]] <- formula[[lf]][[2]]
    }
    else if (deparse(formula[[lf]][[1]]) == "("){
      form <- formula[[lf]][[2]]
      if (length(form) == 3 && form[[1]] == "|")
        formula1[[lf]] <- form[[2]]
    }
  }
  formula1
}

.FLXgetGroupingVar <- function(x)
{
  lf <- length(x)
  while (lf > 1) {
    x <- x[[lf]]
    lf <- length(x)
  }
  x
}

.FLXgetGrouping <- function(formula, data)
{
  group <- factor(integer(0))
  formula1 <- RemoveGrouping(formula)
  if (!identical(formula1, formula))
    group <- factor(eval(.FLXgetGroupingVar(formula), data))
  return(list(group=group, formula=formula1))
}

setMethod("FLXgetModelmatrix", signature(model="FLXM"),
          function(model, data, formula, lhs=TRUE, ...)
          {
            formula <- RemoveGrouping(formula)
            if (length(grep("\\|", deparse(model@formula)))) stop("no grouping variable allowed in the model")
            if(is.null(model@formula))
              model@formula = formula

            ## model@fullformula = update.formula(formula, model@formula)
            ## <FIXME>: ist das der richtige weg, wenn ein punkt in beiden
            ## formeln ist?
            model@fullformula = update(terms(formula, data=data), model@formula)
            ## </FIXME>

            if (lhs) {
              mf <- if (is.null(model@terms)) model.frame(model@fullformula, data=data, na.action = NULL)
              else model.frame(model@terms, data=data, na.action = NULL, xlev = model@xlevels)
              model@terms <- attr(mf, "terms")
              response <- as.matrix(model.response(mf))
              model@y <- model@preproc.y(response)
            }
            else {
              mt1 <- if (is.null(model@terms)) terms(model@fullformula, data=data) else model@terms
              mf <- model.frame(delete.response(mt1), data=data, na.action = NULL, xlev = model@xlevels)
              model@terms<- attr(mf, "terms")
              ## <FIXME>: warum war das da???
              ## attr(mt, "intercept") <- attr(mt1, "intercept")
              ## </FIXME>
            }
            X <- model.matrix(model@terms, data=mf)
            model@contrasts <- attr(X, "contrasts")
            model@x <- model@preproc.x(X)
            model@xlevels <- .getXlevels(model@terms, mf)
            model
          })


initPosteriors <- function(k, cluster, N, groups) {
  if(is(cluster, "matrix")){
    postunscaled <- cluster
    if (!is.null(k)) if (k != ncol(postunscaled)) stop("specified k does not match the number of columns of cluster")
  }
  else{
    if(is.null(cluster)){
      if(is.null(k))
        stop("either k or cluster must be specified")
      else
        cluster <- ungroupPriors(as.matrix(sample(seq_len(k), size = sum(groups$groupfirst), replace=TRUE)),
                                 groups$group, groups$groupfirst)
    }
    else{
      cluster <- as(cluster, "integer")
      if (!is.null(k)) if (k != max(cluster)) stop("specified k does not match the values in cluster")
      k <- max(cluster)
    }
    postunscaled <- matrix(0.1, nrow=N, ncol=k)
    for(K in seq_len(k)){
      postunscaled[cluster==K, K] <- 0.9
    }
  }
  postunscaled
}


ungroupPriors <- function(x, group, groupfirst) {
  if (!length(group)) group <- seq_along(groupfirst)
  if (nrow(x) >= length(group[groupfirst])) {
    x <- x[order(as.integer(group[groupfirst])),,drop=FALSE]
    x <- x[as.integer(group),,drop=FALSE]
  }
  x
}


setMethod("flexmix",
          signature(formula = "formula", model="list"),
          function(formula, data=list(), k=NULL, cluster=NULL,
                   model=NULL, concomitant=NULL, control=NULL, weights=NULL)
          {
            mycall = match.call()
            control = as(control, "FLXcontrol")
            if (!is(concomitant, "FLXP")) concomitant <- FLXPconstant()

            groups <- .FLXgetGrouping(formula, data)
            model <- lapply(model, FLXcheckComponent, k, cluster)
            k <- unique(unlist(sapply(model, FLXgetK, k)))
            if (length(k) > 1) stop("number of clusters not specified correctly")

            model <- lapply(model, FLXgetModelmatrix, data, formula)

            groups$groupfirst <-
              if (length(groups$group)) {flexmix:::groupFirst(groups$group)
              }else {rep(TRUE, FLXgetObs(model[[1]]))}

            if (is(weights, "formula")) {
              weights <- model.frame(weights, data = data, na.action = NULL)[,1]
            }
            ## check if the weights are integer
            ## if non-integer weights are wanted modifications e.g.
            ## for classify != weighted and
            ## plot,flexmix,missing-method are needed
            if (!is.null(weights) & !identical(weights, as.integer(weights)))
              stop("only integer weights allowed")
            ## if weights and grouping is specified the weights within each
            ## group need to be the same
            if (!is.null(weights) & length(groups$group)>0) {
              unequal <- tapply(weights, groups$group, function(x) length(unique(x)) > 1)
              if (any(unequal)) stop("identical weights within groups needed")
            }

            postunscaled <- flexmix:::initPosteriors(k, cluster, FLXgetObs(model[[1]]), groups)

            if (ncol(postunscaled) == 1L)
              concomitant <- FLXPconstant()

            concomitant <- FLXgetModelmatrix(concomitant, data = data,
                                             groups = groups)

            ############################# MODEL FIT ####################################################

            z <- FLXfit(model=model, concomitant=concomitant, control=control,
                        postunscaled=postunscaled, groups=groups, weights = weights)

            z@formula = formula
            z@call = mycall
            z@k0 = as.integer(k)
            z
          })


# FLXfit Source Code

# https://r-forge.r-project.org/scm/viewvc.php/locClass/skel/R/flexmixSVM.R?logsort=cvs&view=markup&root=locclass&pathrev=237


setMethod("FLXfit", signature(model="list"),
          function(model, concomitant, control, postunscaled=NULL, groups, weights)
          {
            ### initialize
            k <- ncol(postunscaled)
            N <- nrow(postunscaled)
            control <- flexmix:::allweighted(model, control, weights)
            if(control@verbose>0)
              cat("Classification:", control@classify, "\n")
            if (control@classify %in% c("SEM", "random")) iter.rm <- 0
            group <- groups$group
            groupfirst <- groups$groupfirst
            if(length(group)>0) postunscaled <- flexmix:::groupPosteriors(postunscaled, group)

            logpostunscaled <- log(postunscaled)
            postscaled <- exp(logpostunscaled - flexmix:::log_row_sums(logpostunscaled))

            llh <- -Inf
            if (control@classify %in% c("SEM", "random")) llh.max <- -Inf
            converged <- FALSE
            components <- rep(list(rep(list(new("FLXcomponent")), k)), length(model))


            # included for manual walk

            ### EM
            for(iter in seq_len(control@iter.max)) {
              iter <<- iter # included for FLXMRhyreg
              ### M-Step
              postscaled = flexmix:::.FLXgetOK(postscaled, control, weights) # prob for each row to be in the different classes

              # calculate probabilities for each class
              prior <- if (is.null(weights)){  # here sth not working with FLXMRhyreg
                ungroupPriors(concomitant@fit(concomitant@x, postscaled[groupfirst,,drop=FALSE]),
                              group, groupfirst)
              }else {ungroupPriors(concomitant@fit(concomitant@x, (postscaled/weights)[groupfirst & weights > 0,,drop=FALSE], weights[groupfirst & weights > 0]),
                                   group, groupfirst)
              }
              # Check min.prior
              nok <- if (nrow(prior) == 1) which(prior < control@minprior) else {
                if (is.null(weights)) which(colMeans(prior[groupfirst,]) < control@minprior)  # how does this work with [groupfirst,,drop = FALSE]???
                else which(colSums(prior[groupfirst,] * weights[groupfirst])/sum(weights[groupfirst]) < control@minprior)
              }
              if(length(nok)) { # what does this do?
                if(control@verbose>0)
                  cat("*** Removing", length(nok), "component(s) ***\n")
                prior <- prior[,-nok,drop=FALSE]  # I don´t understand drop = FALSE here
                prior <- prior/rowSums(prior)
                postscaled <- postscaled[,-nok,drop=FALSE]
                postscaled[rowSums(postscaled) == 0,] <- if (nrow(prior) > 1) prior[rowSums(postscaled) == 0,]
                else prior[rep(1, sum(rowSums(postscaled) == 0)),]
                postscaled <- postscaled/rowSums(postscaled)
                if (!is.null(weights)) postscaled <- postscaled * weights
                k <- ncol(prior)
                if (k == 0) stop("all components removed")
                if (control@classify=="random") {
                  llh.max <- -Inf
                  iter.rm <- iter
                }
                model <- lapply(model, FLXremoveComponent, nok)
                components <- lapply(components, "[", -nok)
              }

              ## ANWENDUNG DES M STEPS
              # components enthält geschätze Parameter für jede Klasse (und noch mehr, s.components@)
              components <- lapply(seq_along(model), function(i) FLXmstep(model[[i]], postscaled, components[[i]])) # gives parameter for model components, can be generateed by mle as well
              postunscaled <- matrix(0, nrow = N, ncol = k)
              for (n in seq_along(model))
                postunscaled <- postunscaled + FLXdeterminePostunscaled(model[[n]], components[[n]]) # using logLik? # https://github.com/cran/flexmix/blob/master/R/flexmix.R
              if(length(group)>0)
                postunscaled <- flexmix:::groupPosteriors(postunscaled, group)

              #FLXdeterminePostunscaled: matrix(sapply(components[[n]], function(x) x@logLik(model[[n]]@x, model[[n]]@y)), nrow = nrow(model[[n]]@y))
              # components[[n]][[1]]@logLik

              ### E-Step
              ## Code changed thanks to Nicolas Picard
              ## to avoid problems with small likelihoods
              postunscaled <- if (nrow(prior) > 1) {postunscaled + log(prior)
              } else {sweep(postunscaled, 2, log(prior), "+")} # why???
              logpostunscaled <- postunscaled
              postunscaled <- exp(postunscaled)
              postscaled <- exp(logpostunscaled - flexmix:::log_row_sums(logpostunscaled))
              ##<FIXME>: wenn eine beobachtung in allen Komonenten extrem
              ## kleine postunscaled-werte hat, ist exp(-postunscaled)
              ## numerisch Null, und damit postscaled NaN
              ## log(rowSums(postunscaled)) ist -Inf
              ##</FIXME>
              if (any(is.nan(postscaled))) {
                index <- which(as.logical(rowSums(is.nan(postscaled))))
                postscaled[index,] <- if(nrow(prior)==1) rep(prior, each = length(index)) else prior[index,]
                postunscaled[index,] <- .Machine$double.xmin
              }
              ### check convergence
              llh.old <- llh
              llh <- if (is.null(weights)) sum(flexmix:::log_row_sums(logpostunscaled[groupfirst,,drop=FALSE]))
              else sum(log_row_sums(logpostunscaled[groupfirst,,drop=FALSE])*weights[groupfirst])
              if(is.na(llh) | is.infinite(llh))
                stop(paste(formatC(iter, width=4),
                           "Log-likelihood:", llh))
              if (abs(llh-llh.old)/(abs(llh)+0.1) < control@tolerance){
                if(control@verbose>0){
                  flexmix:::printIter(iter, llh)
                  cat("converged\n")
                }
                converged <- TRUE
                break
              }
              if (control@classify=="random") {
                if (llh.max < llh) {
                  components.max <- components
                  prior.max <- prior
                  postscaled.max <- postscaled
                  postunscaled.max <- postunscaled
                  llh.max <- llh
                }
              }
              if(control@verbose && (iter%%control@verbose==0))
                flexmix:::printIter(iter, llh)
            } # end of for loop EM Algo?


            ### Construct return object
            if (control@classify=="random") {
              components <- components.max
              prior <- prior.max
              postscaled <- postscaled.max
              postunscaled <- postunscaled.max
              llh <- llh.max
              iter <- control@iter.max - iter.rm
            }

            components <- lapply(seq_len(k), function(i) lapply(components, function(x) x[[i]]))
            names(components) <- paste("Comp", seq_len(k), sep=".")
            cluster <- max.col(postscaled)
            size <-  if (is.null(weights)) tabulate(cluster, nbins=k) else tabulate(rep(cluster, weights), nbins=k)
            names(size) <- seq_len(k)
            concomitant <- flexmix:::FLXfillConcomitant(concomitant, postscaled[groupfirst,,drop=FALSE], weights[groupfirst])
            df <- concomitant@df(concomitant@x, k) + sum(sapply(components, sapply, slot, "df"))
            control@nrep <- 1
            prior <- if (is.null(weights)) colMeans(postscaled[groupfirst,,drop=FALSE])
            else colSums(postscaled[groupfirst,,drop=FALSE] * weights[groupfirst])/sum(weights[groupfirst])

            retval <- new("flexmix", model=model, prior=prior,
                          posterior=list(scaled=postscaled,
                                         unscaled=postunscaled),
                          weights = weights,
                          iter=iter, cluster=cluster, size = size,
                          logLik=llh, components=components,
                          concomitant=concomitant,
                          control=control, df=df, group=group, k=as(k, "integer"),
                          converged=converged)
            retval
          })


# While this is called in line 248
# different parameter for the different classes are estimated
setMethod("FLXmstep", signature(model = "FLXM"), function(model, weights, components, ...) {
  if ("component" %in% names(formals(model@fit)))
    sapply(seq_len(ncol(weights)), function(k) model@fit(model@x, model@y, weights[,k], component = components[[k]]@parameters))
  else
    sapply(seq_len(ncol(weights)), function(k) model@fit(model@x, model@y, weights[,k]))
})

