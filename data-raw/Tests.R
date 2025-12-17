
## TESTs ####

#devtools::installgithub("intelligentaccident/EQ5D_data")
library(EQ5Ddata)
library(flexmix)
library(hyreg2)

###############################
### USE SIMULATED_DATA_NORM ###
###############################

#classic
formula <- y ~  -1 + x1 + x2 + x3 | id

# non-classic
formula = y ~   1/exp( INTERCEPT + x1 * beta1 + x2 * beta2) | id

k <- 2

stv1 <- setNames(c(0.2,0.2,0.2,1,1),c(c("INTERCEPT","beta1","beta2"),c("sigma","theta")))
stv2<- setNames(c(2,1,0.5,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
control = list(iter.max = 1000, verbose = 4)



formula = formula
data =  simulated_data_norm
type =  simulated_data_norm$type
stv = stv1
k = k

type_cont = "TTO"
type_dich = "DCE_A"
opt_method = "L-BFGS-B"
latent = "cont"
id_col = "id"
variables_both = NULL
variables_cont = NULL
variables_dich = NULL
formula_type_classic = FALSE
upper = 2
lower = -Inf



hyflex_mod <- hyreg2(formula = formula,
                     data =  simulated_data_norm,
                     type =  simulated_data_norm$type,
                     stv = stv,
                     k = k,
                     type_cont = type_cont,
                     type_dich = type_dich,
                     opt_method = opt_method,
                     upper = upper,
                     lower = lower,
                     control = control,
                     latent = latent ,
                     id_col = id_col,
                     variables_both =   variables_both,
                     variables_cont =   variables_cont,
                     variables_dich =  variables_dich,
                     formula_type_classic = formula_type_classic
)




summary(hyflex_mod)
summary_hyreg2(hyflex_mod)

plot_hyreg2(data = simulated_data_norm,
            x = "id",
            y= "y",
            id_col = "id",
            class_df_model = give_class(data = simulated_data_norm,
                                        model = hyflex_mod,
                                        id_col = "id"))



### hyreg2_het ###
formula <- y ~  -1 + x1 + x2 + x3 | id
formula_sigma <- y ~   x1 + x2 + x3


stv <- setNames(c(0.2,0.2,0.2,1),c(c("x1","x2","x3"),c("theta")))
stv_sigma <- setNames(c(0.2,0.2,1,3),c("x1","x2","x3","(Intercept)"))

data =  simulated_data_norm
type =  simulated_data_norm$type
k = 2

type_cont = "TTO"
type_dich = "DCE_A"
opt_method = "L-BFGS-B"
latent = "cont"
id_col = "id"
variables_both = NULL
variables_cont = NULL
variables_dich = NULL
upper = Inf
lower = -Inf

hyflex_mod_het <- hyreg2_het(formula = formula,
                         formula_sigma = formula_sigma,
                         stv_sigma = stv_sigma,
                     data =  simulated_data_norm,
                     type =  simulated_data_norm$type,
                     stv = stv,
                     k = k,
                     type_cont = type_cont,
                     type_dich = type_dich,
                     opt_method = opt_method,
                     upper = upper,
                     lower = lower,
                     control = control,
                     latent = latent ,
                     id_col = id_col,
                     variables_both =   variables_both,
                     variables_cont =   variables_cont,
                     variables_dich =  variables_dich
)

summary_hyreg2(hyflex_mod_het)

plot_hyreg2(data = simulated_data_norm,
            x = "id",
            y= "y",
            id_col = "id",
            class_df_model = give_class(data = simulated_data_norm,
                                        model = hyflex_mod_het,
                                        id_col = "id"))


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






####################################################
########## TEST FOR HETEROSCEDASTICITY #############
####################################################


TTOonly <- hyregdata[hyregdata$method == "TTO" & hyregdata$fb_flagged == 0 & hyregdata$state_id > 0,]
DCEonly <- hyregdata[hyregdata$method == "DCE_A" & hyregdata$state_id < 197,]

# use only subdataset for faster estimations
TTOonly <- TTOonly[1:250,]
DCEonly <- DCEonly[1:150,]

data <- rbind(TTOonly,DCEonly)

# model
formula <- value ~ -1 + mo2 + sc2 + ua2 + pd2 + ad2 + mo3 + sc3 + ua3 + pd3 + ad3 +
  mo4 + sc4 + ua4 + pd4 + ad4 + mo5 + sc5 + ua5 + pd5 + ad5 | id

k <- 2

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
                   latent = "cont",
                   id_col = "id"
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



