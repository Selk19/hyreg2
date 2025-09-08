
## TESTs ####

library(flexmix)
library(EQ5Ddata)

#################
### USE DATA1 ###
#################


# random numbers from normal dist
formula <- y ~  -1 + x1 + x2 + x3 | id

k <- 2

stv <- setNames(c(0.2,0.2,0.2,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
control = list(iter.max = 1000, verbose = 4)

rm(counter)

hyflex_mod <- hyreg2(formula = formula,
                     data =  simulated_data_norm,
                     type =  simulated_data_norm$type,
                     stv = stv,
                     k = k,
                     type_cont = "TTO",
                     type_dich = "DCE_A",
                     opt_method = "L-BFGS-B",
                     control = control,
                     latent = "both",
                     id_col = "id"
)


summary(hyflex_mod)
summary_hyreg2(hyflex_mod)
#parameters(hyflex_mod, component=1) # works only for latent = "both"
#parameters(hyflex_mod, component=2)  # works only for latent = "both"

# if latent was "both"
(sum(hyflex_mod@cluster == simulated_data_norm$class))/dim(simulated_data_norm)[1]

# if latent was "cont" or "dich"
proof <- merge(unique(simulated_data_norm[,c("id","class")]),hyflex_mod[["id_classes"]], by = "id")
sum((proof$class == proof$mod_comp)/dim(proof)[1])


# with latent == "dich" we get warning
# In bbmle::mle2(minuslogl = logLik2, start = stv_new, optimizer = optimizer,  :
#                  convergence failure: code=52 (ERROR: ABNORMAL_TERMINATION_IN_LNSRCH)



### RESULT ###

# WITH latent = "both" and | ID in formula #
# 74 % of datapoints are classified to the correct class
# prop 1: 67% (402 of 600), prop 2: 33 % (198 of 600)
# parameter estimation is well for comp 2 of the data but bad for comp 1


# WITH latent = "cont" and | id in formula #
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


k <- 2

control = list(iter.max = 5000, verbose = 5)
stv2 <- setNames(c(rep(0.2,20),1,1),c(colnames(data)[17:36],c("sigma","theta")))
stv1 <- setNames(c(rep(0.1,20),1,1),c(colnames(data)[17:36],c("sigma","theta")))
stv <- list(stv1,stv2) # not working yet (we get back only start values?)

# if formula has an intercept, use this
stvint <- setNames(c(rep(0.1,20),1,1,1),c(colnames(data)[17:36],c("sigma","theta","(Intercept)")))


mod1 <- hyreg2(formula = formula,
               data = data,
               type = data$method,
               stv = stv,
               #   upper = 2,
               #   lower = 0,
               k = k,
               type_cont = "TTO",
               type_dich = "DCE_A",
               opt_method = "L-BFGS-B",
               control = control,
               latent = "both",
               id_col = "id",
               # variables_cont = c("mo5","sc3"),
               # variables_both = c("mo2","sc2","ua2","pd2","ad2","mo3","sc3","ua3","pd3","ad3",
               #  "mo4","sc4","ua4","pd4","ad4","ua5","pd5", "ad5")
)

# if you get an Error like this:
# Error in names(object) <- nm : attempt to set an attribute on NULL
# use rm(counter) and try again

# using stv as list:
# estimates are just the start values except sigma
# why does this happen? optimizer seem not to work correct here


### SUMMARY ###

summary(mod1)
summary_hyreg2(mod1)



# check only TTO with lm #
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




hyb <- hyreg(modformula, data, datatype = "d_method")# ll = 0, ul = 2)
hyb



### RESULT ###
# for k= 1 we get the same estimates as with xreg
# censoring (using lower and upper) also leads to the same estimates as with xreg function hyreg

# k = 2 works (with and without censoring),
# but the results can not be evaluated, since we do not know the true parameter values

# partial coefficients using variables_cont, variables_dich and variables_both work
# Error messages occur in situations, where the entries are wrong (that is good)

# sometimes I get
# Warning message:
# In bbmle::mle2(minuslogl = logLik2, start = stv_new, optimizer = optimizer,  :
# convergence failure: code=1 (NEW_X)
# but having a close look we seem not to have any convergence problem ?!
# I do not know whats the problem or what causes this warning...

# not using |id in formula is problematic, because than it says TTO data are the one class and DCE data are the other
# --> this is not what we want, therefore we need to use | id




####################################
#### USING SIMULATED DATA (EQ5D) ###
####################################

### with simulated_data ###

formula <- y ~ -1 + mo2 + sc2 + ua2 + pd2 + ad2 + mo3 + sc3 + ua3 + pd3 + ad3 +
  mo4 + sc4 + ua4 + pd4 + ad4 + mo5 + sc5 + ua5 + pd5 + ad5 | id


k <- 2


control = list(iter.max = 5000, verbose = 5)
stv2 <- setNames(c(rep(0.1,20),1,1),c(colnames(simulated_data)[3:22],c("sigma","theta")))
stv1 <- setNames(c(rep(0.1,20),1,1),c(colnames(simulated_data)[3:22],c("sigma","theta")))
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
               latent = "cont",
               id_col = "id"
)

# if you get an Error like this:
# Error in names(object) <- nm : attempt to set an attribute on NULL
# use rm(counter) and try again


### SUMMARY ###

summary(mod1)
summary_hyreg2(mod1)

# check if 1 = class 1 in data or not
# proportion of correct classification (latent = "both")
(sum(mod1@cluster == simulated_data$class))/dim(simulated_data)[1]


# if latent was "cont" or "dich"
proof <- merge(unique(simulated_data[,c("id","class")]),mod1[["id_classes"]], by = "id")
sum((proof$class == proof$mod_comp)/dim(proof)[1])




### RESULT ###
# for k = 2 and latent = "cont"
#  96 % of datapoints are classified to the correct class
# the estimates of class 2 (from simulated_data) are very well and close to the true values from simulation
# for class 1 it is less accurate, problematic are especially low values close to zero like 0.005 (lower than 0.1 ?)
# but I think in general it is okay



#########################
### SIMULATED_DATA_MO ###
#########################

#### Using simulated_data_mo ####

formula <- y ~ -1 + mo2 + mo3 + mo4 +  mo5 | id


k <- 2

control = list(iter.max = 5000, verbose = 5)
stv_mo <- setNames(c(rep(0.3,4),1,1),c(colnames(simulated_data_mo)[3:6],c("sigma","theta")))


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


# latent was "both"
(sum(modMO@cluster != simulated_data_mo$class))/dim(simulated_data_mo)[1]


# if latent was "cont" or "dich"
proof <- merge(unique(simulated_data_mo[,c("id","class")]),modMO[["id_classes"]], by = "id")
sum((proof$class == proof$mod_comp)/dim(proof)[1])




### RESULT ###
# for k = 2 and latent = "both"
# 62.5 % of datapoints are classified to the correct class
# some estimates seem to be mixed between classes:
# class 1 from data: estimate of mo2 and mo4 are close to true values,
# but estimates of mo3 and mo5 are close to true values of class 2
# all in all estimates and classification is not statisfying here


# for k = 2 and latent = "cont"
# we get many different results between 50 % and 97 % correct classification

#  with loglik 619.5833 in first step
# we get 97 % of datapoints classified to the correct class
# estimates for both classes are very close to the true values then

# --> latent = "cont" can lead to very good results, but when?
# depending on start values? different start values for different components not working yet




### use getstv to generate stv values ###

# use latent = "cont" and k = 1 in code above to generate modMO

modMO2 <- hyreg2(formula = formula,
                data = simulated_data_mo,
                type = simulated_data_mo$type,
                stv = getstv(modMO),
                # upper = 2,
                # lower = 0,
                k = 2,
                type_cont = "TTO",
                type_dich = "DCE_A",
                opt_method = "L-BFGS-B",
                control = control,
                latent = "cont",
                id_col = "id"
)


summary_hyreg2(modMO2)
