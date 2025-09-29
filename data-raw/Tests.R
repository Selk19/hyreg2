
## TESTs ####

#devtools::installgithub("intelligentaccident/EQ5D_data")
library(EQ5Ddata)
library(flexmix)
library(hyreg2)

###############################
### USE SIMULATED_DATA_NORM ###
###############################


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
                     latent = "cont",
                     id_col = "id"
)


summary(hyflex_mod)
summary_hyreg2(hyflex_mod)

# ratio of correct classification:
# if latent was "both"
(sum(hyflex_mod@cluster == simulated_data_norm$class))/dim(simulated_data_norm)[1]

# if latent was "cont" or "dich"
proof <- merge(unique(simulated_data_norm[,c("id","class")]),hyflex_mod[["id_classes"]], by = "id")
sum((proof$class == proof$mod_comp)/dim(proof)[1])




############################
### USING EQ5D DATA SETS ###
############################


# data
TTOonly <- hyregdata[hyregdata$method == "TTO" & hyregdata$fb_flagged == 0 & hyregdata$state_id > 0,]
DCEonly <- hyregdata[hyregdata$method == "DCE_A" & hyregdata$state_id < 197,]

# use only subdataset for faster estimations
#TTOonly <- TTOonly[1:250,]
#DCEonly <- DCEonly[1:150,]

data <- rbind(TTOonly,DCEonly)


# model
formula <- value ~ -1 + mo2 + sc2 + ua2 + pd2 + ad2 + mo3 + sc3 + ua3 + pd3 + ad3 +
  mo4 + sc4 + ua4 + pd4 + ad4 + mo5 + sc5 + ua5 + pd5 + ad5 | id


k <- 1

control = list(iter.max = 5000, verbose = 5)
stv2 <- setNames(c(rep(0.2,20),1,1),c(colnames(data)[17:36],c("sigma","theta")))
stv1 <- setNames(c(rep(0.1,20),1,1),c(colnames(data)[17:36],c("sigma","theta")))
stv <- list(stv1,stv2)

# if formula has an intercept, use this as stv
# stvint <- setNames(c(rep(0.1,20),1,1,1),c(colnames(data)[17:36],c("sigma","theta","(Intercept)")))


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
                variables_cont = c("mo5","sc5"),
                variables_both = c("mo2","sc2","ua2","pd2","ad2","mo3","sc3","ua3","pd3","ad3",
                 "mo4","sc4","ua4","pd4","ad4","ua5","pd5", "ad5")
)

# if you get an Error like this:
# Error in names(object) <- nm : attempt to set an attribute on NULL
# use rm(counter) and try again



### SUMMARY ###

summary(mod1)
summary_hyreg2(mod1)




# compare to xreg#
#devtools::installgithub("intelligentaccident/xreg")
library(xreg)

modformula <- value ~
  mo2 * MO2 + sc2 * SC2 + ua2 * UA2 + pd2 * PD2 +  ad2 * AD2 +
  mo3 * MO3 + sc3 * SC3 + ua3 * UA3 + pd3 * PD3 + ad3 * AD3 +
  mo4 * MO4 + sc4 * SC4 + ua4 * UA4 + pd4 * PD4 + ad4 * AD4 +
  mo5 * MO5 + sc5 * SC5 + ua5 * UA5 + pd5 * PD5 + ad5 * AD5




hyb <- hyreg(modformula, data, datatype = "d_method")# ll = 0, ul = 2)
hyb






####################################
#### USING SIMULATED_DATA (EQ5D) ###
####################################


formula <- y ~ -1 + mo2 + sc2 + ua2 + pd2 + ad2 + mo3 + sc3 + ua3 + pd3 + ad3 +
  mo4 + sc4 + ua4 + pd4 + ad4 + mo5 + sc5 + ua5 + pd5 + ad5 | id


k <- 2


control = list(iter.max = 5000, verbose = 5)
stv2 <- setNames(c(rep(0.1,20),1,1),c(colnames(simulated_data)[3:22],c("sigma","theta")))
stv1 <- setNames(c(rep(0.1,20),1,1),c(colnames(simulated_data)[3:22],c("sigma","theta")))
#stv <- list(stv1,stv2)


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

# proportion of correct classification:
# check if 1 = class 1 in data or not, maybe change == into !=
# latent = "both"
(sum(mod1@cluster == simulated_data$class))/dim(simulated_data)[1]


# if latent was "cont" or "dich"
proof <- merge(unique(simulated_data[,c("id","class")]),mod1[["id_classes"]], by = "id")
sum((proof$class == proof$mod_comp)/dim(proof)[1])




#########################
### SIMULATED_DATA_MO ###
#########################

#### Using simulated_data_mo ####

formula <- y ~ -1 + mo2 + mo3 + mo4 +  mo5 |id


k <- 2

control = list(iter.max = 5000, verbose = 5)
stv_mo <- setNames(c(rep(0.3,4),1,1),c(colnames(simulated_data_mo)[3:6],c("sigma","theta")))

stv_mo1 <- setNames(c(rep(0.1,4),1,1),c(colnames(simulated_data_mo)[3:6],c("sigma","theta")))
stv_mo2 <- setNames(c(rep(0.6,4),1,1),c(colnames(simulated_data_mo)[3:6],c("sigma","theta")))
stvl <- list(stv_mo1,stv_mo2)


modMO <- hyreg2(formula = formula,
                data = simulated_data_mo,
                type = simulated_data_mo$type,
                stv = stvl,
                # upper = 2,
                # lower = 0,
                k = k,
                type_cont = "TTO",
                type_dich = "DCE_A",
                opt_method = "L-BFGS-B",
                control = control,
                latent = "cont",
               # classes_only = TRUE,
                id_col = "id"

)

summary(modMO)
summary_hyreg2(modMO)


# proportion of correct classification:
# latent was "both"
(sum(modMO@cluster != simulated_data_mo$class))/dim(simulated_data_mo)[1]


# if latent was "cont" or "dich"
proof <- merge(unique(simulated_data_mo[,c("id","class")]),modMO[["id_classes"]], by = "id")
sum((proof$class == proof$mod_comp)/dim(proof)[1])



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









####################################################
########## TEST FOR HETEROSCEDASTICITY #############
####################################################


TTOonly <- hyregdata[hyregdata$method == "TTO" & hyregdata$fb_flagged == 0 & hyregdata$state_id > 0,]
DCEonly <- hyregdata[hyregdata$method == "DCE_A" & hyregdata$state_id < 197,]

# use only subdataset for faster estimations
#TTOonly <- TTOonly[1:250,]
#DCEonly <- DCEonly[1:150,]

data <- rbind(TTOonly,DCEonly)

# model
formula <- value ~ -1 + mo2 + sc2 + ua2 + pd2 + ad2 + mo3 + sc3 + ua3 + pd3 + ad3 +
  mo4 + sc4 + ua4 + pd4 + ad4 + mo5 + sc5 + ua5 + pd5 + ad5 | id

k <- 1

control = list(iter.max = 10000, verbose = 5)

# for sigma estimation
formula_sigma <- value ~  mo2 + sc2 + ua2 + pd2 + ad2 + mo3 + sc3 + ua3 + pd3 + ad3 +
  mo4 + sc4 + ua4 + pd4 + ad4 + mo5 + sc5 + ua5 + pd5 + ad5
stvs <- setNames(c(rep(0.1,20),1),c(colnames(data)[17:36],c("theta")))
stv_sigma <- setNames(c(rep(0.1,20),1),c(colnames(data)[17:36],c("(Intercept)")))



mod1 <- hyreg2_het(formula = formula,
                   formula_sigma = formula_sigma, # if not provided, same as formula is taken
                   data = data,
                   type = data$method,
                   stv = stvs,
                   stv_sigma = stv_sigma, # if not provided all variables from formula_sigma set to 0.1
                   #   upper = 2,
                   #   lower = 0,
                   k = k,
                   type_cont = "TTO",
                   type_dich = "DCE_A",
                   opt_method = "L-BFGS-B",
                   control = control,
                   latent = "both",
                   id_col = "id",
                  # variables_cont = c("mo5","sc5"),
                 #  variables_both = c("mo2","sc2","ua2","pd2","ad2","mo3","sc3","ua3","pd3","ad3",
                  #  "mo4","sc4","ua4","pd4","ad4","ua5","pd5", "ad5")
)

# if you get an Error like this:
# Error in names(object) <- nm : attempt to set an attribute on NULL
# use rm(counter) and try again



### SUMMARY ###
summary_hyreg2(mod1)



### compare to xreg ###
library(xreg)

modformula <- value ~
  mo2 * MO2 + sc2 * SC2 + ua2 * UA2 + pd2 * PD2 +  ad2 * AD2 +
  mo3 * MO3 + sc3 * SC3 + ua3 * UA3 + pd3 * PD3 + ad3 * AD3 +
  mo4 * MO4 + sc4 * SC4 + ua4 * UA4 + pd4 * PD4 + ad4 * AD4 +
  mo5 * MO5 + sc5 * SC5 + ua5 * UA5 + pd5 * PD5 + ad5 * AD5


modformula_het <- value ~ INTERCEPT +
  mo2 * HMO2 + sc2 * HSC2 + ua2 * HUA2 + pd2 * HPD2 + ad2 * HAD2 +
  mo3 * HMO3 + sc3 * HSC3 + ua3 * HUA3 + pd3 * HPD3 + ad3 * HAD3 +
  mo4 * HMO4 + sc4 * HSC4 + ua4 * HUA4 + pd4 * HPD4 + ad4 * HAD4 +
  mo5 * HMO5 + sc5 * HSC5 + ua5 * HUA5 + pd5 * HPD5 + ad5 * HAD5

hyb <- hyreg(modformula, data, datatype = "d_method", hetcont = modformula_het)# ll = 0, ul = 2)
hyb



