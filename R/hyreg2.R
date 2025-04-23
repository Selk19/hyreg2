

#' hyreg2
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
#' @param offset offset as in flexmix
#' @param optimizer optimizer to be used in bbmle::mle2, default = "optim"
#' @param opt_method optimization method to be used in optimizer, default = "BFGS"
#' @param lower opt_method must be set to "L-BFGS-B", lower bound for censored data
#' @param upper opt_method must be set to "L-BFGS-B", upper bound for censored data
#' @param latent one of "both","cont" or "dich", see details
#' @param id_col character, column name containing participant ids
#' @param variables_both character vactor; variables to be fitted on TTO and DCE data, if not specified all variables from formula are used
#' @param variables_cont character vactor; variables to be fitted only on TTO data
#' @param variables_dich character vactor; variables to be fitted only on DCE data
#' @param ... additional arguments for flexmix::flexmix or bbmle::mle2
#'
#' @return model of type flemix
#'
#' @author Kim Rand & Svenja Elkenkamp
#' @examples Put Example here

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
                  variables_both = NULL,
                  variables_dich = NULL,
                  variables_cont = NULL,
                  # additional arguments for flexmix or optimizer ?
                  #cluster = NULL,
                  #concomitant=NULL,
                  #weights=NULL,

                  # non linear regression not implemented yet
                  #   Xb in FLXMRhyreg must be computed differently for that

                  # FLXMRhyreg_het is not implemnted here, maybe better write new function hyreg2_het?


                  ...){



if(is.null(variables_both)){
  variables_both <- names(stv)[!is.element(names(stv),c("sigma","theta"))]
}

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
                            upper = upper))


   fit <- flexmix::flexmix(formula = formula, data = data, k = k, model = model, control = control)
   rm(counter, envir = .GlobalEnv) # counter wird während flexmix erstellt

   return(fit)

 }else{

   idframe <- data.frame(id = data[,id_col],type)
   idcount <- as.data.frame(table(unique(idframe)))
   #as.character(idcount[idcount$Freq == 0,"id"])
   data <- data[!is.element(as.character(data[,id_col]), as.character(idcount[idcount$Freq == 0,"id"])),]
   type <- idframe[!is.element(as.character(idframe[,"id"]), as.character(idcount[idcount$Freq == 0,"id"])),"type"]
   # type anpassen ohne die entsprechenden Zeilen zu den IDs, die raus sind


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
                              upper = upper))


     mod <- flexmix::flexmix(formula = formula, data = data_cont, k = k, model = model, control = control)
     rm(counter, envir = .GlobalEnv) # counter wird während flexmix erstellt

     data_cont$mod_comp <- mod@cluster

     data$roworder <- 1:nrow(data)
     data <- merge(data, unique(data_cont[,c(id_col,"mod_comp")]), by = id_col)
     data <- data[order(data$roworder), ]

   }

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
                              upper = upper))


     mod <- flexmix::flexmix(formula = formula, data = data_dich, k = k, model = model, control = control)
     rm(counter, envir = .GlobalEnv) # counter wird während flexmix erstellt

     data_dich$mod_comp <- mod@cluster

     data$roworder <- 1:nrow(data)
     data <- merge(data, unique(data_dich[,c(id_col,"mod_comp")]), by = id_col)
     data <- data[order(data$roworder), ]
   }


   data_list <- list()
   for(i in unique(mod@cluster)){
     data_list[[i]] <- data[data$mod_comp == i,]
   }

   mod_list <- lapply(data_list, function(xy){
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
                              upper = upper))
     mod <- flexmix::flexmix(formula = formula, data = xy, k = 1, model = model, control = control)
     rm(counter, envir = .GlobalEnv)
     return(mod)
   })

   return(mod_list)
 }
}

##############################

#' summary_hyreg2
#'
#' @description get model parameters of hyreg2 model
#' @param object modelobject generated with hyreg2
#' @return summary object of bbmle::mle2 model
#'
#' @author Kim Rand & Svenja Elkenkamp
#' @examples Put Example here
#' @importFrom bbmle mle2
#' @export


summary_hyreg2 <- function(object){
  if(class(object) == "list"){
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



##################################


## TEST ####


library(flexmix)
library(EQ5Ddata)
TTOonly <- hyregdata[hyregdata$method == "TTO" & hyregdata$fb_flagged == 0 & hyregdata$state_id > 0,]
DCEonly <- hyregdata[hyregdata$method == "DCE_A" & hyregdata$state_id < 197,]

TTOonly <- TTOonly[1:250,]
DCEonly <- DCEonly[1:150,]
data <- rbind(TTOonly,DCEonly)

formula <- value ~ -1 + mo2 + sc2 + ua2 + pd2 + ad2 + mo3 + sc3 + ua3 + pd3 + ad3 +
  mo4 + sc4 + ua4 + pd4 + ad4 + mo5 + sc5 + ua5 + pd5 + ad5 | id

k <- 2

#cluster <- NULL
#concomitant=NULL
#control=NULL
control = list(iter.max = 500, verbose = 5)
#weights=NULL
stv <- setNames(c(rep(0.1,20),1,1),c(colnames(data)[17:36],c("sigma","theta")))


mod1 <- hyreg2(formula = formula,
       data = data,
       type = data$method,
       stv = stv,
       k = k,
       type_cont = "TTO",
       type_dich = "DCE_A",
       control = control,
       latent = "cont", # "dich" not working yet
       id_col = "id"
#      variables_cont = c("mo5","sc5"),
#      variables_both = c("mo2","sc2","ua2","pd2","ad2","mo3","sc3","ua3","pd3","ad3",
#      "mo4","sc4","ua4","pd4","ad4","ua5","pd5")
      )

# if you get an Error like this:
# Error in names(object) <- nm : attempt to set an attribute on NULL
# use rm(counter) and try again

summary(mod1)
summary_hyreg2(mod1)
#xreg2:::refit(mod1)
