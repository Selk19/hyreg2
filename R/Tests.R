

## TESTs ####
#

library(flexmix)
library(EQ5Ddata)
TTOonly <- hyregdata[hyregdata$method == "TTO" & hyregdata$fb_flagged == 0 & hyregdata$state_id > 0,]
DCEonly <- hyregdata[hyregdata$method == "DCE_A" & hyregdata$state_id < 197,]

TTOonly <- TTOonly[1:250,]
DCEonly <- DCEonly[1:150,]
data <- rbind(TTOonly,DCEonly)

formula <- value ~ -1 + mo2 + sc2 + ua2 + pd2 + ad2 + mo3 + sc3 + ua3 + pd3 + ad3 +
  mo4 + sc4 + ua4 + pd4 + ad4 + mo5 + sc5 + ua5 + pd5 + ad5 |id

# with simulated_data
formula <- y ~ -1 + mo2 + sc2 + ua2 + pd2 + ad2 + mo3 + sc3 + ua3 + pd3 + ad3 +
  mo4 + sc4 + ua4 + pd4 + ad4 + mo5 + sc5 + ua5 + pd5 + ad5 | id


k <- 2


#cluster <- NULL
#concomitant=NULL
#control=NULL
control = list(iter.max = 5000, verbose = 5)
#weights=NULL
stv2 <- setNames(c(rep(0.1,20),1,1),c(colnames(data)[17:36],c("sigma","theta")))
stv1 <- setNames(c(rep(0.1,20),1,1),c(colnames(data)[17:36],c("sigma","theta")))
#stv <- list(stv1,stv2) # not working


mod1 <- hyreg2(formula = formula,
              # data = data,
              # type = data$method,
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
               latent = "both", # "dich" not working yet
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
#xreg2:::refit(mod1)

# proportion of correct classification
(sum(mod1@cluster == simulated_data$class))/dim(simulated_data)[1]

x <- 4
mod1@cluster[c(x,x+120)]



# compare to xreg#
modformula <- value ~
  mo2 * MO2 + sc2 * SC2 + ua2 * UA2 + pd2 * PD2 +  ad2 * AD2 +
  mo3 * MO3 + sc3 * SC3 + ua3 * UA3 + pd3 * PD3 + ad3 * AD3 +
  mo4 * MO4 + sc4 * SC4 + ua4 * UA4 + pd4 * PD4 + ad4 * AD4 +
  mo5 * MO5 + sc5 * SC5 + ua5 * UA5 + pd5 * PD5 + ad5 * AD5

hyb <- hyreg(modformula, data, datatype = "d_method")
hyb

hybU <- hyreg(modformula, data, datatype = "d_method",init = hyb,  ul = 2)
hybU

hybC <- hyreg(modformula, data, datatype = "d_method",init = hyb, ll = 0, ul = 2)
hybC




### with simulated data ###
modformula <- y ~
  mo2 * MO2 + sc2 * SC2 + ua2 * UA2 + pd2 * PD2 +  ad2 * AD2 +
  mo3 * MO3 + sc3 * SC3 + ua3 * UA3 + pd3 * PD3 + ad3 * AD3 +
  mo4 * MO4 + sc4 * SC4 + ua4 * UA4 + pd4 * PD4 + ad4 * AD4 +
  mo5 * MO5 + sc5 * SC5 + ua5 * UA5 + pd5 * PD5 + ad5 * AD5

simulated_data$type_num <- simulated_data$type == "TTO"

hybC <- hyreg(modformula, simulated_data, datatype = "type_num")
hybC


