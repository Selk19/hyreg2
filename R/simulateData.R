
# Set parameters
set.seed(42)
n_samples_tto <- 120
n_samples_dce <- 120
sigma <- 0.1
theta <- 0.8
stv <- c(0.02, 0.15, 0.4, 0.6, # MO
         0.01, 0.2, 0.36, 0.5, # SC
         0.015, 0.25, 0.5, 0.8, # UA
         0.1, 0.3, 0.4, 0.6, # PD
         0.09, 0.19, 0.6, 0.7) # AD



### DUMMIES ###
dummy_mo <- as.data.frame(model.matrix(~ factor(sample(1:5,n_samples_tto,replace = TRUE),levels = 1:5) - 1))
colnames(dummy_mo) <- paste0("MO", 1:5)
dummy_mo <- dummy_mo[,-1]


dummy_sc <- as.data.frame(model.matrix(~ factor(sample(1:5,n_samples_tto,replace = TRUE),levels = 1:5) - 1))
colnames(dummy_sc) <- paste0("SC", 1:5)
dummy_sc <- dummy_sc[,-1]


dummy_ua <- as.data.frame(model.matrix(~ factor(sample(1:5,n_samples_tto,replace = TRUE),levels = 1:5) - 1))
colnames(dummy_ua) <- paste0("UA", 1:5)
dummy_ua <- dummy_ua[,-1]


dummy_pd <- as.data.frame(model.matrix(~ factor(sample(1:5,n_samples_tto,replace = TRUE),levels = 1:5) - 1))
colnames(dummy_pd) <- paste0("PD", 1:5)
dummy_pd <- dummy_pd[,-1]


dummy_ad <- as.data.frame(model.matrix(~ factor(sample(1:5,n_samples_tto,replace = TRUE),levels = 1:5) - 1))
colnames(dummy_ad) <- paste0("AD", 1:5)
dummy_ad <- dummy_ad[,-1]



# alle zusammen
dummy_df <- cbind(dummy_mo, dummy_sc)
dummy_df <- cbind(dummy_df, dummy_ua)
dummy_df <- cbind(dummy_df, dummy_pd)
dummy_df <- cbind(dummy_df, dummy_ad)




### ONE COMPONENT ###

# linear predictor
Xb_tto <- as.matrix(dummy_df) %*% stv
Xb_dce <- (as.matrix(dummy_df) %*% stv) * theta

# y
y_tto <- rnorm(n_samples_tto, mean = Xb_tto, sd = sigma)
logistic_tmp <- 0.5 + 0.5 * tanh(Xb_dce / 2)
y_dce <- rbinom(n_samples_dce, size = 1, prob = logistic_tmp)

### dataset ###
data_tto <- data.frame(
  type = rep("TTO", n_samples_tto),
  y = y_tto,
  mo2 = dummy_df[, 1],
  mo3 = dummy_df[, 2],
  mo4 = dummy_df[, 3],
  mo5 = dummy_df[, 4],
  sc2 = dummy_df[, 5],
  sc3 = dummy_df[, 6],
  sc4 = dummy_df[, 7],
  sc5 = dummy_df[, 8],
  ua2 = dummy_df[, 9],
  ua3 = dummy_df[, 10],
  ua4 = dummy_df[, 11],
  ua5 = dummy_df[, 12],
  pd2 = dummy_df[, 13],
  pd3 = dummy_df[, 14],
  pd4 = dummy_df[, 15],
  pd5 = dummy_df[, 16],
  ad2 = dummy_df[, 17],
  ad3 = dummy_df[, 18],
  ad4 = dummy_df[, 19],
  ad5 = dummy_df[, 20])

# data_tto <- data.frame(
#   type = rep("tto", n_samples_tto),
#   y = y_tto,
#   x1 = x_tto[, 1],
#   x2 = x_tto[, 2],
#   x3 = x_tto[, 3]
# )


data_dce <- data.frame(
  type = rep("DCE_A", n_samples_dce),
  y = y_dce,
  mo2 = dummy_df[, 1],
  mo3 = dummy_df[, 2],
  mo4 = dummy_df[, 3],
  mo5 = dummy_df[, 4],
  sc2 = dummy_df[, 5],
  sc3 = dummy_df[, 6],
  sc4 = dummy_df[, 7],
  sc5 = dummy_df[, 8],
  ua2 = dummy_df[, 9],
  ua3 = dummy_df[, 10],
  ua4 = dummy_df[, 11],
  ua5 = dummy_df[, 12],
  pd2 = dummy_df[, 13],
  pd3 = dummy_df[, 14],
  pd4 = dummy_df[, 15],
  pd5 = dummy_df[, 16],
  ad2 = dummy_df[, 17],
  ad3 = dummy_df[, 18],
  ad4 = dummy_df[, 19],
  ad5 = dummy_df[, 20])

simulated_data <- rbind(data_tto, data_dce)



# check with base R functions for each type
lm(y~ -1 + mo2 + mo3 + mo4 + mo5 +
     sc2 + sc3 + sc4 + sc5 +
     ua2 + ua3 + ua4 + ua5 +
     pd2 + pd3 + pd4 + pd5 +
     ad2 + ad3 + ad4 + ad5,
  # data = data_tto,
    data = simulated_data[simulated_data$type == "TTO",])



glm(y~ -1 + mo2 + mo3 + mo4 + mo5 +
      sc2 + sc3 + sc4 + sc5 +
      ua2 + ua3 + ua4 + ua5,
    data = simulated_data[simulated_data$type == "DCE_A",], family = "binomial")

