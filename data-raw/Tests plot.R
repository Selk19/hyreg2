### TEST PLOT ###

### TEST simulated_data_norm ###

# true classes
plot_hyreg2(data = simulated_data_norm,
           x = "id",
           y = "y",
           id_col = "id",
           id_df_model = simulated_data_norm[,c("id","class")])

# model classification
plot_hyreg2(data = simulated_data_norm,
           x = "x2",
           y = "y",
           id_col = "id",
           id_df_model = give_id(data = simulated_data_norm,
                                 model = hyflex_mod,
                                 id_col = "id"))

# only TTO
# true classes
plot_hyreg2(data = simulated_data_norm,
           x = "x2",
           y = "y",
           id_col = "id",
           id_df_model = simulated_data_norm[,c("id","class")])

# model classification
plot_hyreg2(data = simulated_data_norm,
           x ="id",
           y = "y",
           id_col = "id",
           id_df_model = give_id(simulated_data_norm,hyflex_mod,"id"),
           type_to_plot = list("type","DCE_A"))





#### simulated_data ####


simulated_data$mo <- ifelse(simulated_data$mo2 == 1, 2,
                            ifelse(simulated_data$mo3 == 1, 3,
                                   ifelse(simulated_data$mo4 == 1, 4,
                                          ifelse(simulated_data$mo5 == 1, 5, 1))))

simulated_data$sc <- ifelse(simulated_data$sc2 == 1, 2,
                            ifelse(simulated_data$sc3 == 1, 3,
                                   ifelse(simulated_data$sc4 == 1, 4,
                                          ifelse(simulated_data$sc5 == 1, 5, 1))))

simulated_data$ua <- ifelse(simulated_data$ua2 == 1, 2,
                            ifelse(simulated_data$sc3 == 1, 3,
                                   ifelse(simulated_data$sc4 == 1, 4,
                                          ifelse(simulated_data$sc5 == 1, 5, 1))))

simulated_data$pd <- ifelse(simulated_data$pd2 == 1, 2,
                            ifelse(simulated_data$pd3 == 1, 3,
                                   ifelse(simulated_data$pd4 == 1, 4,
                                          ifelse(simulated_data$pd5 == 1, 5, 1))))

simulated_data$ad <- ifelse(simulated_data$ad2 == 1, 2,
                            ifelse(simulated_data$ad3 == 1, 3,
                                   ifelse(simulated_data$ad4 == 1, 4,
                                          ifelse(simulated_data$ad5 == 1, 5, 1))))



# true classes
plot_hyreg2(data = simulated_data,
           x = "mo",
           y = "y",
           id_col = "id",
           id_df_model = simulated_data[,c("id","class")],
           colors = c("#F8766D","#00BFC4"))

# model classification
plot_hyreg2(data = simulated_data,
           x = "mo",
           y = "y",
           id_col = "id",
           id_df_model = give_id(data = simulated_data,
                                 model = mod1,
                                 id_col = "id"),
           colors = c("#00BFC4","#F8766D"))





#### simulated_data_mo ####


# code dummies back
simulated_data_mo$mo <- ifelse(simulated_data_mo$mo2 == 1, 2,
                               ifelse(simulated_data_mo$mo3 == 1, 3,
                                      ifelse(simulated_data_mo$mo4 == 1, 4,
                                             ifelse(simulated_data_mo$mo5 == 1, 5, 1))))



give_id(simulated_data_mo,modMO)

# if model was estimated without | id
simulated_data_mo$idn <- rownames(simulated_data_mo)

# all ids
plot_hyreg2(simulated_data_mo,
           "mo",
           "y",
           "id",
           give_id(simulated_data_mo,modMO,"id"))

# only TTO ids
plot_hyreg2(simulated_data_mo[simulated_data_mo$type == "TTO",],
           "id",
           "y",
           "id",
           give_id(simulated_data_mo,modMO,"id")[simulated_data_mo$type == "TTO",])


