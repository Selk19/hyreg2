

# these datasets can be used in Tests file to test the functions


######################################
### Random Numbers from NormalDist ###
######################################

### COMP 1 ###
# Set parameters
set.seed(42)
n_samples_tto <- 200
n_samples_dce <- 100
sigma <- 1.0
theta <- 5
stv <- c(0.5, -0.3, 0.8)  # coef for x

# simulate matrix
x_tto <- matrix(rnorm(n_samples_tto * length(stv)), ncol = length(stv))
x_dce <- matrix(rnorm(n_samples_dce * length(stv)), ncol = length(stv))

# Xb (linear predictors) and y
Xb_tto <- x_tto %*% stv
Xb_dce <- (x_dce %*% stv) * theta

y_tto <- rnorm(n_samples_tto, mean = Xb_tto, sd = sigma)  # continuous outcomes
logistic_tmp <- 0.5 + 0.5 * tanh(Xb_dce / 2)
y_dce <- rbinom(n_samples_dce, size = 1, prob = logistic_tmp)  # binary outcomes


### dataset ###
data_tto <- data.frame(
  type = rep("TTO", n_samples_tto),
  y = y_tto,
  x1 = x_tto[, 1],
  x2 = x_tto[, 2],
  x3 = x_tto[, 3]
)

data_dce <- data.frame(
  type = rep("DCE_A", n_samples_dce),
  y = y_dce,
  x1 = x_dce[, 1],
  x2 = x_dce[, 2],
  x3 = x_dce[, 3]
)


simulated_data1 <- rbind(data_tto, data_dce)



### COMP 2 ###
# Set parameters
set.seed(37)
n_samples_tto <- 200
n_samples_dce <- 100
sigma <- 0.5
theta <- 2
stv <- c(1.4, 2.3, -0.2)  # coef for x


# simulate matrix
x_tto <- matrix(rnorm(n_samples_tto * length(stv)), ncol = length(stv))
x_dce <- matrix(rnorm(n_samples_dce * length(stv)), ncol = length(stv))

# Xb (linear predictors) and y
Xb_tto <- x_tto %*% stv
Xb_dce <- (x_dce %*% stv) * theta

y_tto <- rnorm(n_samples_tto, mean = Xb_tto, sd = sigma)  # continuous outcomes
logistic_tmp <- 0.5 + 0.5 * tanh(Xb_dce / 2)
y_dce <- rbinom(n_samples_dce, size = 1, prob = logistic_tmp)  # binary outcomes


### dataset ###
data_tto <- data.frame(
  type = rep("TTO", n_samples_tto),
  y = y_tto,
  x1 = x_tto[, 1],
  x2 = x_tto[, 2],
  x3 = x_tto[, 3]
)

data_dce <- data.frame(
  type = rep("DCE_A", n_samples_dce),
  y = y_dce,
  x1 = x_dce[, 1],
  x2 = x_dce[, 2],
  x3 = x_dce[, 3]
)


simulated_data2 <- rbind(data_tto, data_dce)




### whole dataset with two comp. ###
simulated_data1$c <- 1
simulated_data2$c <- 2

simulated_data1$ID <- c(1:100, 1:100, 1:100)
simulated_data2$ID <- c(101:200, 101:200, 101:200)


simulated_data <- rbind(simulated_data1,simulated_data2)
simulated_data_too <- simulated_data[simulated_data$type == "TTO",]
simulated_data_dce <- simulated_data[simulated_data$type == "DCE_A",]


data1 <- simulated_data
rm(simulated_data)
rm(simulated_data1)
rm(simulated_data2)




###########################
#### ONLY FOR MOBILITY ####
###########################



### COMP 1 ###

# Set parameters
set.seed(42)
n_samples_tto <- 120
n_samples_dce <- 120
sigma <- 0.001
theta <- 0.2
stv <- c(0.005, 0.01, 0.08, 0.1) # MO




### DUMMIES ###
dummy_mo <- as.data.frame(model.matrix(~ factor(sample(1:5,n_samples_tto,replace = TRUE),levels = 1:5) - 1))
colnames(dummy_mo) <- paste0("MO", 1:5)
dummy_mo <- dummy_mo[,-1]



# linear predictor
Xb_tto <- as.matrix(dummy_mo) %*% stv
Xb_dce <- (as.matrix(dummy_mo) %*% stv) * theta

# y
y_tto <- rnorm(n_samples_tto, mean = Xb_tto, sd = sigma)
logistic_tmp <- 0.5 + 0.5 * tanh(Xb_dce / 2)
y_dce <- rbinom(n_samples_dce, size = 1, prob = logistic_tmp)

### dataset ###
data_tto <- data.frame(
  type = rep("TTO", n_samples_tto),
  y = y_tto,
  mo2 = dummy_mo[, 1],
  mo3 = dummy_mo[, 2],
  mo4 = dummy_mo[, 3],
  mo5 = dummy_mo[, 4])



data_dce <- data.frame(
  type = rep("DCE_A", n_samples_dce),
  y = y_dce,
  mo2 = dummy_mo[, 1],
  mo3 = dummy_mo[, 2],
  mo4 = dummy_mo[, 3],
  mo5 = dummy_mo[, 4])

simulated_data1 <- rbind(data_tto, data_dce)
simulated_data1$class <- 1
simulated_data1$id <- c(1:n_samples_tto,1:n_samples_dce)



# check with base R functions for each type
lm(y~ -1 + mo2 + mo3 + mo4 + mo5,
   # data = data_tto,
   data = simulated_data1[simulated_data1$type == "TTO",])



glm(y~ -1 + mo2 + mo3 + mo4 + mo5,
    data = simulated_data1[simulated_data1$type == "DCE_A",], family = "binomial")





### COMP 2 ###

# Set parameters
set.seed(84)
n_samples_tto <- 120
n_samples_dce <- 120
sigma <- 0.1
theta <- 2
stv <- c(0.2, 0.4, 0.6, 0.8) # MO




### DUMMIES ###
dummy_mo <- as.data.frame(model.matrix(~ factor(sample(1:5,n_samples_tto,replace = TRUE),levels = 1:5) - 1))
colnames(dummy_mo) <- paste0("MO", 1:5)
dummy_mo <- dummy_mo[,-1]



# linear predictor
Xb_tto <- as.matrix(dummy_mo) %*% stv
Xb_dce <- (as.matrix(dummy_mo) %*% stv) * theta

# y
y_tto <- rnorm(n_samples_tto, mean = Xb_tto, sd = sigma)
logistic_tmp <- 0.5 + 0.5 * tanh(Xb_dce / 2)
y_dce <- rbinom(n_samples_dce, size = 1, prob = logistic_tmp)

### dataset ###
data_tto <- data.frame(
  type = rep("TTO", n_samples_tto),
  y = y_tto,
  mo2 = dummy_mo[, 1],
  mo3 = dummy_mo[, 2],
  mo4 = dummy_mo[, 3],
  mo5 = dummy_mo[, 4])



data_dce <- data.frame(
  type = rep("DCE_A", n_samples_dce),
  y = y_dce,
  mo2 = dummy_mo[, 1],
  mo3 = dummy_mo[, 2],
  mo4 = dummy_mo[, 3],
  mo5 = dummy_mo[, 4])

simulated_data2 <- rbind(data_tto, data_dce)
simulated_data2$class <- 2
simulated_data2$id <- c((1+n_samples_tto):(2*n_samples_tto) , (1+n_samples_dce):(2*n_samples_dce))



### all together ###
simulated_data_mo <- rbind(simulated_data1, simulated_data2)

rm(simulated_data1,simulated_data2)




###########################
### ALL DIMENSIONS EQ5D ###
###########################

### COMP 1 ###

# Set parameters
set.seed(42)
n_samples_tto <- 120
n_samples_dce <- 120
sigma <- 0.1
theta <- 3
stv <- c(0.001, 0.05, 0.08, 0.1, # MO
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

simulated_data1 <- rbind(data_tto, data_dce)
simulated_data1$class <- 1
simulated_data1$id <- c(1:n_samples_tto,1:n_samples_dce)



# check with base R functions for each type
lm(y~ -1 + mo2 + mo3 + mo4 + mo5 +
     sc2 + sc3 + sc4 + sc5 +
     ua2 + ua3 + ua4 + ua5 +
     pd2 + pd3 + pd4 + pd5 +
     ad2 + ad3 + ad4 + ad5,
  # data = data_tto,
    data = simulated_data1[simulated_data1$type == "TTO",])



glm(y~ -1 + mo2 + mo3 + mo4 + mo5 +
      sc2 + sc3 + sc4 + sc5 +
      ua2 + ua3 + ua4 + ua5,
    data = simulated_data1[simulated_data1$type == "DCE_A",], family = "binomial")





### COMP 2 ###

# Set parameters
set.seed(84)
n_samples_tto <- 120
n_samples_dce <- 120
sigma <- 0.02
theta <- 2
stv <- c(0.2, 0.4, 0.6, 0.8, # MO
         0.1, 0.3, 0.4, 0.5, # SC
         0.2, 0.25, 0.6, 0.7, # UA
         0.05, 0.2, 0.27, 0.8, # PD
         0.15, 0.35, 0.4, 0.65) # AD



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

simulated_data2 <- rbind(data_tto, data_dce)
simulated_data2$class <- 2
simulated_data2$id <- c((1+n_samples_tto):(2*n_samples_tto) , (1+n_samples_dce):(2*n_samples_dce))



### all together ###
simulated_data <- rbind(simulated_data1, simulated_data2)





