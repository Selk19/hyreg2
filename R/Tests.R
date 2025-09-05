

## TESTs ####

library(flexmix)
library(EQ5Ddata)

#################
### USE DATA1 ###
#################


# random numbers from normal dist
formula <- y ~  -1 + x1 + x2 + x3 | ID

k <- 2

stv <- setNames(c(0.2,0,1,1,1),c(colnames(data1)[3:5],c("sigma","theta")))
control = list(iter.max = 1000, verbose = 4)

rm(counter)

hyflex_mod <- hyreg2(formula = formula,
                     data =  data1,
                     type =  data1$type,
                     stv = stv,
                     k = k,
                     type_cont = "TTO",
                     type_dich = "DCE_A",
                     opt_method = "L-BFGS-B",
                     control = control,
                     latent = "cont",
                     id_col = "ID"
)


summary(hyflex_mod)
summary_hyreg2(hyflex_mod)
parameters(hyflex_mod, component=1)
parameters(hyflex_mod, component=2)

(sum(hyflex_mod@cluster == data1$c))/dim(data1)[1]

# if latent was "cont" or "dich"
proof <- merge(unique(data1[,c("ID","c")]),hyflex_mod[["id_classes"]], by = "ID")
sum((proof$c == proof$mod_comp)/dim(proof)[1])


# with latent == "dich" we get warning
# In bbmle::mle2(minuslogl = logLik2, start = stv_new, optimizer = optimizer,  :
#                  convergence failure: code=52 (ERROR: ABNORMAL_TERMINATION_IN_LNSRCH)



### RESULT ###

# WITH latent = "both" and | ID in formula #
# 74 % of datapoints are classified to the correct class
# prop 1: 67% (402 of 600), prop 2: 33 % (198 of 600)
# parameter estimation is well for comp 2 of the data but bad for comp 1


# WITH latent = "cont" and | ID in formula #
# 95 % of datapoints are classified to the correct class
# # prop 1: 50% (300 of 600), prop 2: 50 % (300 of 600)
# parameter estimation is well ans close to the true values, except for theta in one class ( estimated as 3, true value is 5)

# --> code seem to work fine with latent = "cont"


############################
### USING EQ5D DATA SETS ###
############################


# data
TTOonly <- hyregdata[hyregdata$method == "TTO" & hyregdata$fb_flagged == 0 & hyregdata$state_id > 0,]
DCEonly <- hyregdata[hyregdata$method == "DCE_A" & hyregdata$state_id < 197,]

TTOonly <- TTOonly[1:250,]
DCEonly <- DCEonly[1:150,]
data <- rbind(TTOonly,DCEonly)


# model
formula <- value ~ -1 + mo2 + sc2 + ua2 + pd2 + ad2 + mo3 + sc3 + ua3 + pd3 + ad3 +
  mo4 + sc4 + ua4 + pd4 + ad4 + mo5 + sc5 + ua5 + pd5 + ad5 | id


k <- 1

control = list(iter.max = 5000, verbose = 5)
stv2 <- setNames(c(rep(0.1,20),1,1),c(colnames(data)[17:36],c("sigma","theta")))
stv1 <- setNames(c(rep(0.1,20),1,1),c(colnames(data)[17:36],c("sigma","theta")))
#stv <- list(stv1,stv2) # not working


mod1 <- hyreg2(formula = formula,
               data = data[data$method == "TTO",],
               type = data[data$method == "TTO",]$method,
               stv = stv1,
               # upper = 2,
               # lower = 0,
               k = k,
               type_cont = "TTO",
               type_dich = "DCE_A",
               opt_method = "L-BFGS-B",
               control = control,
               latent = "both",
               id_col = "id"
               #      variables_cont = c("mo5","sc5"),
               #      variables_both = c("mo2","sc2","ua2","pd2","ad2","mo3","sc3","ua3","pd3","ad3",
               #      "mo4","sc4","ua4","pd4","ad4","ua5","pd5")
)

# if you get an Error like this:
# Error in names(object) <- nm : attempt to set an attribute on NULL
# use rm(counter) and try again


### SUMMARY ###

summary(mod1)
summary_hyreg2(mod1)



# check if 1 = class 1 in data or not
# proportion of correct classification
(sum(mod1@cluster == simulated_data$class))/dim(simulated_data)[1]



# TTO lm #
# lm(value ~ -1 + mo2 + sc2 + ua2 + pd2 + ad2 + mo3 + sc3 + ua3 +
#      pd3 + ad3 + mo4 + sc4 + ua4 + pd4 + ad4 + mo5 + sc5 + ua5 +
#      pd5 + ad5,
#    data =data[data$method == "TTO",])



# compare to xreg#
library(xreg)

modformula <- value ~
  mo2 * MO2 + sc2 * SC2 + ua2 * UA2 + pd2 * PD2 +  ad2 * AD2 +
  mo3 * MO3 + sc3 * SC3 + ua3 * UA3 + pd3 * PD3 + ad3 * AD3 +
  mo4 * MO4 + sc4 * SC4 + ua4 * UA4 + pd4 * PD4 + ad4 * AD4 +
  mo5 * MO5 + sc5 * SC5 + ua5 * UA5 + pd5 * PD5 + ad5 * AD5




hybC <- hyreg(modformula, data, datatype = "d_method", ll = 0, ul = 2)
hybC



### RESULT ###
# for k= 1 we get the same estimates as with xreg
# censoring (using lower and upper) also leads to the same estimates as with xreg function hyreg

# k = 2 works, but the results can not be evaluated, since we do not know the true parameter values


####################################
#### USING SIMULATED DATA (EQ5D) ###
####################################

### with simulated_data ###

formula <- y ~ -1 + mo2 + sc2 + ua2 + pd2 + ad2 + mo3 + sc3 + ua3 + pd3 + ad3 +
  mo4 + sc4 + ua4 + pd4 + ad4 + mo5 + sc5 + ua5 + pd5 + ad5 | id


k <- 1


control = list(iter.max = 5000, verbose = 5)
stv2 <- setNames(c(rep(0.1,20),1,1),c(colnames(data)[17:36],c("sigma","theta")))
stv1 <- setNames(c(rep(0.1,20),1,1),c(colnames(data)[17:36],c("sigma","theta")))
#stv <- list(stv1,stv2) # not working


mod1 <- hyreg2(formula = formula,
               data = simulated_data,
               type = simulated_data$type,
               stv = stv1,
               # upper = 2,
               # lower = 0,
               k = k,
               type_cont = "TTO",
               type_dich = "DCE_A",
               opt_method = "L-BFGS-B",
               control = control,
               latent = "both",
               id_col = "id"
               #      variables_cont = c("mo5","sc5"),
               #      variables_both = c("mo2","sc2","ua2","pd2","ad2","mo3","sc3","ua3","pd3","ad3",
               #      "mo4","sc4","ua4","pd4","ad4","ua5","pd5")
)

# if you get an Error like this:
# Error in names(object) <- nm : attempt to set an attribute on NULL
# use rm(counter) and try again


### SUMMARY ###

summary(mod1)
summary_hyreg2(mod1)

# check if 1 = class 1 in data or not
# proportion of correct classification
(sum(mod1@cluster == simulated_data$class))/dim(simulated_data)[1]



# TTO lm #
# lm(value ~ -1 + mo2 + sc2 + ua2 + pd2 + ad2 + mo3 + sc3 + ua3 +
#      pd3 + ad3 + mo4 + sc4 + ua4 + pd4 + ad4 + mo5 + sc5 + ua5 +
#      pd5 + ad5,
#    data =simulated_data[simulated_data$type == "TTO",])



# compare to xreg#

library(xreg)

### with simulated data ###
modformula <- y ~
  mo2 * MO2 + sc2 * SC2 + ua2 * UA2 + pd2 * PD2 +  ad2 * AD2 +
  mo3 * MO3 + sc3 * SC3 + ua3 * UA3 + pd3 * PD3 + ad3 * AD3 +
  mo4 * MO4 + sc4 * SC4 + ua4 * UA4 + pd4 * PD4 + ad4 * AD4 +
  mo5 * MO5 + sc5 * SC5 + ua5 * UA5 + pd5 * PD5 + ad5 * AD5

simulated_data$type_num <- simulated_data$type == "TTO"

hyb <- hyreg(modformula, simulated_data, datatype = "type_num")
hyb



#########################
### SIMULATED_DATA_MO ###
#########################

#### Using simulated_data_mo ####

formula <- y ~ -1 + mo2 + mo3 + mo4 +  mo5


k <- 2

control = list(iter.max = 5000, verbose = 5)
stv_mo <- setNames(c(rep(0.1,4),1,1),c(colnames(simulated_data_mo)[3:6],c("sigma","theta")))


modMO <- hyreg2(formula = formula,
               data = simulated_data_mo,
               type = simulated_data_mo$type,
               stv = stv_mo,
               # upper = 2,
               # lower = 0,
               k = k,
               type_cont = "TTO",
               type_dich = "DCE_A",
               opt_method = "L-BFGS-B",
               control = control,
               latent = "cont",
               id_col = "id"
)

summary(modMO)
summary_hyreg2(modMO)


(sum(modMO@cluster == simulated_data_mo$class))/dim(simulated_data_mo)[1]



# if latent was "cont" or "dich"
proof <- merge(unique(simulated_data_mo[,c("id","class")]),modMO[["id_classes"]], by = "id")
sum((proof$class == proof$mod_comp)/dim(proof)[1])




