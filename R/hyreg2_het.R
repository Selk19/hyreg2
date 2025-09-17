
#' hyreg2_het: function for model estimation for EQ5D valueset data accounting for heteroscedasticity in continous data
#'
#'
#' @description Estimation of hybrid model for EQ-5D data
#' @param formula linear model formula, using “| x” will include a grouping variable x, see details
#' @param formula_sigma formula for sigma estimation, if not provided formula is taken without grouping
#'  variable, see details
#' @param data a dataframe containing the data, see details
#' @param type either the name of the column in data containing an indicator whether that datapoint (row)
#'  contains continuous or dichotomous data as character, or the whole vector containing the indicator.
#' @param type_cont value of type referring to continuous data
#' @param type_dich value of type referring to dichotomous data
#' @param k numeric, number of latent classes to be estimated via flexmix::flexmix
#' @param control control vector for flexmix::flexmix
#' @param stv named vector or list of named vectors containing start values for all coefficients
#'  formula, including sigma and theta, see details
#' @param stv_sigma named vector with start values for sigma estimation.
#'  Have to be the same variables as given in formula_sigma. Should NOT end with _h, see details
#' @param offset offset as in flexmix::flexmix
#' @param optimizer optimizer to be used in bbmle::mle2, default = "optim"
#' @param opt_method optimization method to be used in optimizer, default = "BFGS"
#' @param lower default = -INF, lower bound for censored data. If this is used, opt_method must be
#'  set to "L-BFGS-B"
#' @param upper default = INF, upper bound for censored data. If this is used, opt_method must be
#'  set to "L-BFGS-B",
#' @param latent data type to use in component identification, must be one of "both", "cont" or "dich",
#'  default = “both”,  see details
#' @param id_col character, name of the grouping variable, only needed if latent != "both”, see details
#' @param classes_only logical, default FALSE, indicates whether the function should perform only
#'  classification, rather than both classification and model estimation, only possible for latent != "both",
#'   see datails
#' @param variables_both character vector; variables to be fitted on both continuous and dichotomous data.
#'  If not specified, all variables from formula are used. If provided and not all variables from formula
#'  are included, variables_cont and variables_dich must be provided as well, while one of them can be NULL,
#'   see details.
#' @param variables_cont character vector; variables to be fitted only on continous data. If provided,
#' variables_both and variables_dich must be provided as well.
#' @param variables_dich character vactor; variables to be fitted only on dichotomous data, if provided,
#'  variables_both and variables_cont must be provided as well.
#' @param ... additional arguments for flexmix::flexmix or bbmle::mle2
#'
#' @return model object of type flemix, coefficients named ..._h are coefficients for heteroscedasticity
#'
#' @details
#'
#' estimation of sigma:
#' text here
#'
#' formula:  a typical R formula of the form y ~ x1 + x2 + … should be provided.
#' Additionally, it is possible to include a grouping variable for repeated measures
#' by using “| xg” where xg is the column containing the group-memberships.
#'  Formula will look like this: y ~ x1 + x2 +… | xg.  In flexmix, this is called the concomitant
#'  variable specification: the model is fit conditional on grouping, so that all observations with
#'  the same group are treated as belonging together when computing likelihood contributions.
#'  One possible grouping variable can e.g. be an id number to identify answers by the same participants.
#'  We highly recommend using a grouping variable, since otherwise the algorithm tends to classify all
#'  continuous data into one estimated class and all dichotomous data into the other.
#'
#'
#' dataframe: a dataframe having the following columns:
#' all independent and the dependet variable used in formula,
#' one column for the grouping variable if grouping should be used,e.g. id numbers of participants with
#' repeated measurements,
#' one column indicating if the observations belongs to continuous or dichotomous data
#' with the entries type_cont and type_dich, e.g. a column called "type" with the entries "TTO" for continous
#' datapoints and "DCE" for dichotomous datapoints, type_cont will be "TTO" and type_dich will be "DCE" then.
#' One row should match one observation (one datapoint).
#'
#'
#' start values (stv): the given start values must be a named vector if the same start
#'  values should be used for all latent classes. Otherwise a list of named vectors should be used
#'  if different start values are assumed for each latent class. There has to be one entry in the list
#'  for each latent class.  Each start value vector must include start values for sigma and theta.
#'  Currently, it is necessary to use the names "sigma" and "theta" for these values. If users are unsure
#'  for which variables start values must be provided, this can be checked by calling
#'  colnames(model.matrix(formula,data)). Here, the formula should not include the grouping variable part.
#'
#'
#' latent: in some situations, e.g. for controlling preference heterogeneity, it can be useful to
#' estimate the latent classes only on one type of data and afterwards estimate the model parameters on all
#' data. In such cases, input variable latent can be used to specify on which type of data the classification
#' should be done. If “cont” or “dich” is used, the input parameter id_col must be specified and gives the name,
#' i.e. a character string, of the grouping variable for classification. Some groups may be removed from the
#' data, since they have only continuous or only dichotomous observations. Then in a first step,
#' a model is estimated only on the continuous/ dichotomous data and the achieved classification is stored.
#' In a next step model parameters are estimated on both types of data using this classification, i.e. in the
#' second step, no classification is done but only estimating models with k = 1 on all types of data for the
#' classes estimated in step one. Setting input variable classes_only on TRUE, the second step is left out and
#' directly the estimated classes from step one are given as output.
#'
#'
#' Variables_both, variables_cont, variables_dich: It is possible to specify partial coefficients,
#'  which are used only on continuous or dichotomous data.
#'  Example:  Suppose different models should be specified for continuous and dichotomous  data:
#'  • Model continuous data: y ~  x1 + x3
#'  • Model dichotomous data: y ~  x1 + x2
#'  The formula input to hyreg2 must then include all parameters that occur in either model:
#'   y ~ x1 + x2 + x3
#'  The assignment of parameters to data types is then achieved via the input arguments variables_both, variables_cont, and variables_dich:
#'   variables_both = “x1”, variables_cont = “x3” and  variables_dich = “x2”.
#'  Every variable included in the provided formula (except the grouping variable ) must appear in exactly one of these vectors. One of the variables_ vectors can also be NULL, if no variables should be used only on this type of the data.

#'
#' @author Svenja Elkenkamp, Kim Rand & John Grosser
#' @examples
#'
#'formula <- y ~  -1 + x1 + x2 + x3
#'formula_sigma <- y ~  x1 + x2 + x3
#'
#'k <- 1

#'stv <- setNames(c(0.2,0,1,1),c(colnames(simulated_data_norm)[3:5],c("theta")))
#'stv_sigma <- setNames(c(0.2,0,1,1),c(colnames(simulated_data_norm)[3:5],c("(Intercept)")))
#'control = list(iter.max = 1000, verbose = 4)

#'rm(counter)

#'mod <- hyreg2_het(formula = formula,
#'                    formula_sigma = formula_sigma,
#'                     data =  simulated_data_norm,
#'                     type =  simulated_data_norm$type, # or "type"
#'                     stv = stv,
#'                     stv_sigma = stv_sigma,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = control,
#'                     latent = "both",
#'                     id_col = "id"
#')
#'summary_hyreg2(mod)

#' @importFrom flexmix flexmix
#' @importFrom bbmle mle2
#' @export


hyreg2_het <-function(formula,
                  formula_sigma = NULL,
                  data,
                  type,
                  type_cont,
                  type_dich,
                  k = 1,
                  control = NULL,
                  stv = NULL, # has to include starting values for sigma and theta as well
                  stv_sigma = NULL,
                  offset = NULL,
                  opt_method = "BFGS",
                  optimizer = "optim",
                  lower = -Inf,
                  upper = Inf,
                  latent = "both", # one of "both", "cont", "dich"
                  id_col = NULL, # as character
                  classes_only = FALSE,
                  variables_both = NULL,
                  variables_dich = NULL,
                  variables_cont = NULL,
               #   non_linear = FALSE,
                  # additional arguments for flexmix or optimizer ?

                  # MISSING:
                  # non linear regression not implemented yet
                  #   Xb in FLXMRhyreg must be computed differently for that
                  # FLXMRhyreg_het is not implemnted here, maybe better write new function hyreg2_het?

                  ...){

  dotarg <- list(...)

  if(is.character(type) & length(type) == 1 ){
    type <- data[,type]
  }

  # prepare formula handling
  formula_string <- paste(deparse(formula), collapse = "")
  formula_parts <- strsplit(formula_string, "\\|")[[1]]
  formula_short <- as.formula(formula_parts[1])

  # set formula_sigma if not provided
  if(is.null(formula_sigma)){
    formula_sigma <- formula_short
  }



  ### STV Check ###
  # check stv for names in x (model.matrix(formula,data))

  # CHECK FOR STV
  # no stv
  if(!is.list(stv)){
    if(is.null(stv)){
      warning(paste0("No stv provided. Set all stv from data to 0.1 and theta = 1"))
      stv <- setNames(c(rep(0.1,dim(model.matrix(formula_short,data))[2]),1), c(colnames(model.matrix(formula_short,data)),c("theta")))
    }else{

      # one or more stv missing
      if(any(!is.element(colnames(model.matrix(formula_short,data)),names(stv)))){
        miss <- colnames(model.matrix(formula_short,data))[!is.element(colnames(model.matrix(formula_short,data)),names(stv))]
        stop(paste0("start values missing for ", paste(miss, collapse = ", ") ," Please provide stv values for all relevant variables."
        ))
      }

      # theta missing
      if(!is.element("theta",names(stv))){
        stv <- c(stv,setNames(1,"theta"))
        warning(paste0("No stv for theta provided, set to 1"))
      }

      # sigma in stv must be deleted (since we want to estimate it with stv_sigma)
      if(is.element("sigma",names(stv))){
        stv <- stv[names(stv) != "sigma"]
        warning(paste0("sigma deleted from stv"))
      }

      # stv for variables not in formula given, d.h. zu viele angegeben
      if(any(!is.element(names(stv), c(colnames(model.matrix(formula_short,data)),"theta")))){
        much <- names(stv)[!is.element(names(stv), c(colnames(model.matrix(formula_short,data)),"theta"))]
        stop(paste0("Too many stv provided. ",paste(much, collapse = ", "), " is/are no variable in formula and does not need a stv"))
      }
      #  check order in FLXMRhyreg
    }

  }else{ # stv is a list

    # warning(paste0("stv is a list. Please ensure that variables_both, variables_dich and variables_cont are provided"))

    for(i in 1:length(stv)){
      # one or more stv missing
      if(any(!is.element(colnames(model.matrix(formula_short,data)),names(stv[[i]])))){
        miss <- colnames(model.matrix(formula_short,data))[!is.element(colnames(model.matrix(formula_short,data)),names(stv[[i]]))]
        stop(paste0("start values missing for ", paste(miss, collapse = ", ") ," Please provide stv values for all relevant variables."
        ))
      }

      # theta missing
      if(!is.element("theta",names(stv[[i]]))){
        stv[[i]] <- c(stv[[i]],setNames(1,"theta"))
        warning(paste0("No stv for theta provided, set to 1"))
      }

      # sigma in stv must be deleted (since we want to estimate it with stv_sigma)
      if(is.element("sigma",names(stv))){
        stv <- stv[names(stv) != "sigma"]
        warning(paste0("sigma deleted from stv"))
      }

      # stv for variables not in formula given, d.h. zu viele angegeben
      if(any(!is.element(names(stv[[i]]), c(colnames(model.matrix(formula_short,data)),"theta")))){
        much <- names(stv[[i]])[!is.element(names(stv[[i]]), c(colnames(model.matrix(formula_short,data)),"theta"))]
        stop(paste0("Too many stv provided. ",paste(much, collapse = ", "), " is/are no variable in formula and does not need a stv"))
      }

    }
  }


  # CHECK FOR STV_SIGMA
  # no stv_sigma
  if(!is.list(stv_sigma)){
    if(is.null(stv_sigma)){
      warning(paste0("No stv_sigma provided. Set all stv from formula_sigma to 0.1"))
      stv_sigma <- setNames(c(rep(0.1,dim(model.matrix(formula_sigma,data))[2])), c(colnames(model.matrix(formula_sigma,data))))
    }else{

      # one or more stv_sigma missing
      if(any(!is.element(colnames(model.matrix(formula_sigma,data)),names(stv_sigma)))){
        miss <- colnames(model.matrix(formula_sigma,data))[!is.element(colnames(model.matrix(formula_sigma,data)),names(stv_sigma))]
        stop(paste0("sigma start values missing for ", paste(miss, collapse = ", ") ," Please provide stv_sigma values for all relevant variables."
        ))
      }

      # stv_sigma for variables not in formula given, d.h. zu viele angegeben
      if(any(!is.element(names(stv_sigma), c(colnames(model.matrix(formula_sigma,data)))))){
        much <- names(stv_sigma)[!is.element(names(stv_sigma), c(colnames(model.matrix(formula_sigma,data))))]
        stop(paste0("Too many stv_sigma provided. ",paste(much, collapse = ", "), " is/are no variable in formula_sigm and does not need a stv_sigma"))
      }
      #  check order in FLXMRhyreg_het
    }

  }else{ # stv_sigma is a list


    for(i in 1:length(stv_sigma)){
      # one or more stv missing
      if(any(!is.element(colnames(model.matrix(formula_sigma,data)),names(stv_sigma[[i]])))){
        miss <- colnames(model.matrix(formula_sigma,data))[!is.element(colnames(model.matrix(formula_sigma,data)),names(stv_sigma[[i]]))]
        stop(paste0("sigma start values missing for ", paste(miss, collapse = ", ") ," Please provide stv_sigma values for all relevant variables."
        ))
      }

      # stv for variables not in formula given, d.h. zu viele angegeben
      if(any(!is.element(names(stv_sigma[[i]]), c(colnames(model.matrix(formula_sigma,data)))))){
        much <- names(stv_sigma[[i]])[!is.element(names(stv_sigma[[i]]), c(colnames(model.matrix(formula_sigma,data))))]
        stop(paste0("Too many stv_sigma provided. ",paste(much, collapse = ", "), " is/are no variable in formula_sigma and does not need a stv_sigma"))
      }
    }
  }




  ### TYPE Check ###
  if(is.null(type) | is.null(type_dich) | is.null(type_cont)){
    stop(paste0("inputs for type, type_dich and typ_cont needed"))
  }else{
    if(!is.element(type_dich,unique(type))){
      warning(paste0("Provided type_dich is not part of type"))
    }
    if(!is.element(type_cont,unique(type))){
      warning(paste0("Provided type_cont is not part of type"))
    }
  }



  ### VARIABALES Check ###
  if(!is.list(stv)){
    if(is.null(variables_both) & is.null(variables_dich) & is.null(variables_cont)){
      variables_both <- names(stv)[!is.element(names(stv),c("sigma","theta"))]
    }else{
      if(any(!is.element(names(stv)[!is.element(names(stv),c("sigma","theta"))],
                         c(variables_both,variables_dich,variables_cont)))){
        # check if all variables are included
        stop(paste0("all variables from stv have to be part of exactly one of the vectors variables_both, variables_dich or variables_cont"))
        # alternative: do not provide any of the vectors
        # than all relevant variables are set to variables_both automatically
      }
      if(any(table(c(variables_both,variables_cont,variables_dich))>1)){
        stop(paste0("variables can only be part in one of the vectors variables_both, variables_dich and variables_cont"))
      }
    }

  }else{

    # if stv is a list, both classes have to depend on the same variable set,
    # stv of different classes have same names (hence using stv[[1]] is okay here)
    # different set of variables for each class not supported yet!

    if(is.null(variables_both) & is.null(variables_dich) & is.null(variables_cont)){
      variables_both <- names(stv[[1]])[!is.element(names(stv[[1]]),c("sigma","theta"))]
    }else{
      if(any(!is.element(names(stv[[1]])[!is.element(names(stv[[1]]),c("sigma","theta"))],
                         c(variables_both,variables_dich,variables_cont)))){
        # check if all variables are included
        stop(paste0("all variables from stv have to be part of exactly one of the vectors variables_both, variables_dich or variables_cont"))
        # alternative: do not provide any of the vectors
        # than all relevant variables are set to variables_both automatically
      }
      if(any(table(c(variables_both,variables_cont,variables_dich))>1)){
        stop(paste0("variables can only be part in one of the vectors variables_both, variables_dich and variables_cont"))
      }
    }
  }




  ### NON LINEAR FUNCTIONS ###
  # NOT IMPLEMENTED YET

  formula_orig <- formula
  # if(non_linear == TRUE){
  #   # formula <- function to keep only names of data columns
  # }
  # # for linear functoins formula and formula_orig are the same




  ### ESTIMATION ###
  if(latent == "both"){

    model <- list(FLXMRhyreg_het( data = data,
                              formula_sigma = formula_sigma,
                             type= type,
                             stv = stv,
                             stv_sigma = stv_sigma,
                             type_cont = type_cont,
                             type_dich = type_dich,
                             variables_both = variables_both,
                             variables_cont = variables_cont,
                             variables_dich = variables_dich,
                             opt_method = opt_method,
                             optimizer = optimizer,
                             lower = lower,
                             upper = upper,
                             non_linear = non_linear))
    #formula_orig = formula



    fit <- flexmix::flexmix(formula = formula, data = data, k = k, model = model, control = control)
    rm(counter, envir = .GlobalEnv) # counter will be created during the M-step driver

    return(fit)

  }else{
    # latent = "cont" or "dich"

    # Prepare first step
    if(is.null(id_col)){
      stop("id_col needed")
    }

    idframe <- data.frame(id = data[,id_col],type)
    idcount <- as.data.frame(table(unique(idframe)))
    #as.character(idcount[idcount$Freq == 0,"id"])
    data <- data[!is.element(as.character(data[,id_col]), as.character(idcount[idcount$Freq == 0,"id"])),]
    type <- idframe[!is.element(as.character(idframe[,"id"]), as.character(idcount[idcount$Freq == 0,"id"])),"type"]

    if(any(idcount$Freq == 0)){
      miss <- idcount[idcount$Freq == 0,"id"]
      warning(paste0( "IDs ",paste(miss, collapse = ", "), " were removed, since they were only part in one type of data"))
    }


    # FIRST STEP: GET LATENT CLASSES
    if(latent == "cont"){
      data_cont <- data[type == type_cont,]
      model <- list(FLXMRhyreg( data = data_cont,
                                type= type[type == type_cont],
                               stv = stv,
                               type_cont = type_cont,
                               type_dich = type_dich,
                               variables_both = variables_both,
                               variables_cont = variables_cont,
                               variables_dich = variables_dich,
                               opt_method = opt_method,
                               optimizer = optimizer,
                               lower = lower,
                               upper = upper,
                               non_linear = non_linear,
                               formula_orig = formula_orig))


      mod <- flexmix::flexmix(formula = formula, data = data_cont, k = k, model = model, control = control)
      rm(counter, envir = .GlobalEnv)

      data_cont$mod_comp <- mod@cluster

      data$roworder <- 1:nrow(data)
      data <- merge(data, unique(data_cont[,c(id_col,"mod_comp")]), by = id_col)
      data <- data[order(data$roworder), ]

      # später auch ausgeben können, welche ID zu welcher Klasse zugeordnet wurde
      id_classes <- data_cont[,c(id_col,"mod_comp")]

    }

    if(latent == "dich"){
      data_dich <- data[type == type_dich,]
      model <- list(FLXMRhyreg( data = data_dich,
                                type= type[type == type_dich],
                               stv = stv,
                               type_cont = type_cont,
                               type_dich = type_dich,
                               variables_both = variables_both,
                               variables_cont = variables_cont,
                               variables_dich = variables_dich,
                               opt_method = opt_method,
                               optimizer = optimizer,
                               lower = lower,
                               upper = upper,
                               non_linear = non_linear,
                               formula_orig = formula_orig))


      mod <- flexmix::flexmix(formula = formula, data = data_dich, k = k, model = model, control = control)
      rm(counter, envir = .GlobalEnv) # counter will be created during the M-step driver

      data_dich$mod_comp <- mod@cluster

      data$roworder <- 1:nrow(data)
      data <- merge(data, unique(data_dich[,c(id_col,"mod_comp")]), by = id_col)
      data <- data[order(data$roworder), ]

      # später auch ausgeben können, welche ID zu welcher Klasse zugeordnet wurde
      id_classes <- data_dich[,c(id_col,"mod_comp")]
    }

    # return only estimated classes without model new model coefficients
    if(classes_only == TRUE){
      return(id_classes)
    }

    # SECOND STEP: GET MODEL ESTIMATES
    data_list <- list()
    for(i in unique(mod@cluster)){
      data_list[[i]] <- data[data$mod_comp == i,]
    }


    mod_list <- lapply(data_list, function(xy){
      if(is.null(xy)){
        mod <- NULL
        warning( paste("one or more components are empty. Set mod to NULL"))
      }else{
        model <- list(FLXMRhyreg(type= type[data$mod_comp == unique(xy$mod_comp)],
                                 #type = type,
                                 stv = stv, # stv can be a list
                                 type_cont = type_cont,
                                 type_dich = type_dich,
                                 variables_both = variables_both,
                                 variables_cont = variables_cont,
                                 variables_dich = variables_dich,
                                 opt_method = opt_method,
                                 optimizer = optimizer,
                                 lower = lower,
                                 upper = upper,
                                 non_linear = non_linear,
                                 formula_orig = formula_orig))

        mod <- flexmix::flexmix(formula = formula, data = xy, k = 1, model = model, control = control)
        rm(counter, envir = .GlobalEnv) # counter will be created during the M-step driver
      }
      return(mod)
    })

    mod_list$id_classes <- unique(id_classes)
    return(mod_list)
  }
}

