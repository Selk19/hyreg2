

#' hyreg2: function for model estimation for EQ5D valueset data
#'
#' @description Estimation of hybrid model for EQ-5D data, not implemnted for het yet
#'
#' @param formula Model formula
#' @param data a dataframe containing the data
#' @param type a vector containing an indicator wheter that datapoint contains to TTO or DCE Data
#' @param type_cont indicator for continous data
#' @param type_dich indicator for dichotoums data
#' @param k numeric, number of latent classes to be estimated via flexmix::flexmix
#' @param control control vector for flexmix::flexmix
#' @param stv named vector or list of named vactors containing start values,
#'            has to be a vector if the same start values should be used for all latent classes,
#'            has to be a list of named vectores if different start values are assumed for the latent classes
#'            has to include start values for sigma and theta as well
#'            Using colnames(model.matrix(formula,data)) (formula without |) you can check, which variables need a stv value.
#'
#' @param offset offset as in flexmix
#' @param optimizer optimizer to be used in bbmle::mle2, default = "optim"
#' @param opt_method optimization method to be used in optimizer, default = "BFGS"
#' @param lower opt_method must be set to "L-BFGS-B", lower bound for censored data
#' @param upper opt_method must be set to "L-BFGS-B", upper bound for censored data
#' @param latent one of "both","cont" or "dich", see details
#' @param id_col character, column name containing participant ids, needed if latent != "both
#' @param classes_only logical, default FALSE, indicate if only classification should be done without estimating model parameters
#'                    only possible for latent = "cont" or "dich"
#' @param variables_both character vactor; variables to be fitted on TTO and DCE data, if not specified all variables from formula are used
#' @param variables_cont character vactor; variables to be fitted only on TTO data
#' @param variables_dich character vactor; variables to be fitted only on DCE data
#' @param ... additional arguments for flexmix::flexmix or bbmle::mle2
#'
#' @return model of type flemix
#'
#' @author Svenja Elkenkamp & Kim Rand
#' @examples
#'
#'formula <- y ~  -1 + x1 + x2 + x3 | id
#'
#'k <- 2

#'stv <- setNames(c(0.2,0,1,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
#'control = list(iter.max = 1000, verbose = 4)

#'rm(counter)

#'mod <- hyreg2(formula = formula,
#'                     data =  simulated_data_norm,
#'                     type =  simulated_data_norm$type,
#'                     stv = stv,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = control,
#'                     latent = "cont",
#'                     id_col = "id"
#')
#'summary_hyreg2(mod)

#' @importFrom flexmix flexmix
#' @importFrom bbmle mle2
#' @export


hyreg2 <-function(formula,
                  data,
                  type,
                  type_cont,
                  type_dich,
                  k = 1,
                  control = NULL,
                  stv = NULL, # has to include starting values for sigma and theta as well
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
              #    non_linear = FALSE,
                  # additional arguments for flexmix or optimizer ?

                  # MISSING:
                  # non linear regression not implemented yet
                  #   Xb in FLXMRhyreg must be computed differently for that
                  # FLXMRhyreg_het is not implemnted here, maybe better write new function hyreg2_het?

                  ...){

  dotarg <- list(...)


  # prepare formula handling
  formula_string <- paste(deparse(formula), collapse = "")
  formula_parts <- strsplit(formula_string, "\\|")[[1]]
  formula_short <- as.formula(formula_parts[1])


  ### STV Check ###
  # check stv for names in x (model.matrix(formula,data))

  # no stv
  if(!is.list(stv)){
    if(is.null(stv)){
      warning(paste0("No stv provided. Set all stv from data to 0.1 and sigma = 1 and theta = 1"))
      stv <- setNames(c(rep(0.1,dim(model.matrix(formula_short,data))[2]),1,1), c(colnames(model.matrix(formula_short,data)),c("sigma","theta")))
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

      # sigma missing
      if(!is.element("sigma",names(stv))){
        stv <- c(stv,setNames(1,"sigma"))
        warning(paste0("No stv for sigma provided, set to 1"))
      }

      # stv for variables not in formula given, d.h. zu viele angegeben
      if(any(!is.element(names(stv), c(colnames(model.matrix(formula_short,data)),"theta","sigma")))){
        much <- names(stv)[!is.element(names(stv), c(colnames(model.matrix(formula_short,data)),"theta","sigma"))]
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

        # sigma missing
        if(!is.element("sigma",names(stv[[i]]))){
          stv[[i]] <- c(stv[[i]],setNames(1,"sigma"))
          warning(paste0("No stv for sigma provided, set to 1"))
        }

        # stv for variables not in formula given, d.h. zu viele angegeben
        if(any(!is.element(names(stv[[i]]), c(colnames(model.matrix(formula_short,data)),"theta","sigma")))){
          much <- names(stv[[i]])[!is.element(names(stv[[i]]), c(colnames(model.matrix(formula_short,data)),"theta","sigma"))]
          stop(paste0("Too many stv provided. ",paste(much, collapse = ", "), " is/are no variable in formula and does not need a stv"))
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
#  if(non_linear == TRUE){
    # formula <- function to keep only names of data columns
#  }
  # for linear functoins formula and formula_orig are the same



 ### ESTIMATION ###
 if(latent == "both"){

   model <- list(FLXMRhyreg(type= type,
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
     model <- list(FLXMRhyreg(type= type[type == type_cont],
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
   # SECOND STEP: GET MODEL ESTIMATES
   if(latent == "dich"){
     data_dich <- data[type == type_dich,]
     model <- list(FLXMRhyreg(type= type[type == type_dich],
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


### new summary function ###

#' summary_hyreg2
#'
#' @description get model parameters of model generated by hyreg2 oder hyreg2_het
#' @param object modelobject generated with hyreg2 or hyreg2_het
#' @return summary object of bbmle::mle2 model
#'
#' @author Kim Rand & Svenja Elkenkamp
#' @examples
#'formula <- y ~  -1 + x1 + x2 + x3 | id
#'k <- 2
#'stv <- setNames(c(0.2,0,1,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
#'control = list(iter.max = 1000, verbose = 4)

#'rm(counter)

#'mod <- hyreg2(formula = formula,
#'                     data =  simulated_data_norm,
#'                     type =  simulated_data_norm$type,
#'                     stv = stv,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = control,
#'                     latent = "both",
#'                     id_col = "id"
#')
#'summary_hyreg2(mod)
#'
#' @importFrom bbmle mle2
#' @export


summary_hyreg2 <- function(object){
  if(class(object) == "list"){
    object <- object[1:(length(object) - 1)]
    out <- lapply(object, function(k){
      comp <- k@components
      out_list <- lapply(comp, function(j){bbmle::summary(j[[1]]@parameters[["fit_mle"]])})
    })
  }else{
    comp <- object@components
    out <- lapply(comp, function(j){bbmle::summary(j[[1]]@parameters[["fit_mle"]])})
  }
  return(out)
}


##############################
#' getstv
#'
#' @description function to export coefficent values and names from a model fitted with k = 1 and latent = "both".
#' These values can be used as stv for a new model with k > 1
#' @param mod modeloutput from hyreg2 oder hyreg2_het using k = 1 and latent = "both"
#' @return named vector of parameter estimates from mod. Can be used as stv for additional model estimations
#'
#' @author Kim Rand & Svenja Elkenkamp
#' @examples
#'formula <- y ~  -1 + x1 + x2 + x3 | id
#'
#'k <- 1
#'stv <- setNames(c(0.2,0,1,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
#'control = list(iter.max = 1000, verbose = 4)

#'rm(counter)

#'mod <- hyreg2(formula = formula,
#'                     data =  simulated_data_norm,
#'                     type =  simulated_data_norm$type,
#'                     stv = stv,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = control,
#'                     latent = "both",
#'                     id_col = "id"
#')
#'new_stv <- getstv(mod)
#'
#'### use new_stv in new model ###
#'

# random numbers from normal dist
#'formula <- y ~  -1 + x1 + x2 + x3 | id

#'k <- 2

#'hyflex_mod <- hyreg2(formula = formula,
#'                     data =  simulated_data_norm,
#'                     type =  simulated_data_norm$type,
#'                     stv = new_stv,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = control,
#'                     latent = "both",
#'                     id_col = "id"
#')

#' @export

getstv <- function(mod){
  stv <- c(mod@components[["Comp.1"]][[1]]@parameters$coef,
           mod@components[["Comp.1"]][[1]]@parameters$sigma,
           mod@components[["Comp.1"]][[1]]@parameters$theta)
  return(stv)
}


##################################



# for use of fit_mle
# internal function, not for direct user use
gendf <- function(list){
  df <- list()
  for(j in 1:length(list)){
    df[[j]] <- data.frame("Estimates" = c(list[[j]]@parameters[["coef"]], #  maybe include an other apply function here?
                                          list[[j]]@parameters[["sigma"]],
                                          list[[j]]@parameters[["theta"]]),
                          "Std Error" = bbmle::summary(list[[j]]@parameters[["fit_mle"]])@coef[,2],
                          "pvalue" =bbmle::summary(list[[j]]@parameters[["fit_mle"]])@coef[,4])
    # give AIC, loglik here as additional arguments?
  }

  return(df)
}


### refit function ###

# get parametervalues of models
# object is outcome of hyreg2
refit <- function(object){
  if(class(object) == "list"){
    out <- lapply(object, function(k){
      comp <- k@components
      out_list <- lapply(comp, function(j){gendf(j)})
    })
  }else{
    comp <- object@components
    out <- lapply(comp, function(j){gendf(j)})
  }
  return(out)
}



