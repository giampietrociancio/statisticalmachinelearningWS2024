# sub task 6 -----------------------------------------
set.seed(25012001)
library(tidyverse)
library(FNN)

simul_epe <- function(sample_size = 500, dimension = 10, sd, n_sim = 1000, type) {
  
  epe <- matrix(0, nrow = n_sim, ncol = dimension)
  colnames(epe) <- paste("p =", 1:dimension)
  
  for (sim in 1:n_sim) {
    for (p in 1:dimension) {
      
      # generate simulation data
      x_t <- as.data.frame(matrix(runif(sample_size * p, min = -1, max = 1), nrow = sample_size, ncol = p))
      y_t <- x_t[, 1] + rnorm(sample_size, mean = 0, sd = sd)  # y = f[X] with f[X] = X[,1]
      
      x_0 <- as.data.frame(matrix(rep(0, p), nrow = 1))
      y_0 <- 0 + rnorm(1, 0, sd)
      
      if (type == "lm") {
        lin_model <- lm(y_t ~ ., data = x_t)
        y_hat <- predict(lin_model, x_0)
      }
      
      if (type == "knn") {
        y_hat <- knn.reg(train = x_t, y = y_t, test = x_0, k = 1)$pred
      }
      
      epe_tmp <- (y_0 - y_hat)^2
      epe[sim, p] <- epe_tmp
    }
  }
  
  return(colMeans(epe))
}

average_epe_lin_sd_0 <- simul_epe(sample_size = 500, dimension = 10, sd = 0, n_sim = 1000, type = "lm")
average_epe_lin_sd_1 <- simul_epe(sample_size = 500, dimension = 10, sd = 1, n_sim = 1000, type = "lm")

average_epe_knn_sd_0 <- simul_epe(sample_size = 500, dimension = 10, sd = 0, n_sim = 1000, type = "knn")
average_epe_knn_sd_1 <- simul_epe(sample_size = 500, dimension = 10, sd = 1, n_sim = 1000, type = "knn")

dim <- 1:10

data_sd_0 <- data.frame(
  dim = rep(dim, 2),
  epe = c(average_epe_lin_sd_0, average_epe_knn_sd_0),
  method = factor(rep(c("lm", "knn"), each = length(dim)))
)

data_sd_1 <- data.frame(
  dim = rep(dim, 2),
  epe = c(average_epe_lin_sd_1, average_epe_knn_sd_1),
  method = factor(rep(c("lm", "knn"), each = length(dim)))
)

plot_sd_0 <- ggplot(data_sd_0, aes(x = dim, y = epe, color = method)) +
  geom_line(linewidth = 1) +
  labs(title = "EPE with sd = 0", x = "p", y = "EPE")

plot_sd_1 <- ggplot(data_sd_1, aes(x = dim, y = epe, color = method)) +
  geom_line(linewidth = 1) +
  labs(title = "EPE with sd = 1", x = "p", y = "EPE")

### SD = 0. The linear model remains effectively zero due to no noise, and with 500 observations, it accurately captures the true linear relationship between X1 and Y. The abundance of data prevents over-reliance on any "useless" predictors in higher dimensions.

# For KNN, the 1-NN method starts with an error close to zero, but as dimensionality increases, the error rises due to the curse of dimensionality.The volume of space grows exponentially, making it harder for KNN to find neighbors, leading to a higher likelihood of selecting observations near the boundary.

### SD = 1. With an SD of 1 across all 10 dimensions, the linear model maintains a relatively constant expected prediction error (EPE) around 1, reflecting the variance of Y. Since the true relationship is linear, the model shows no bias and effectively manages the "useless" predictors, driving their beta coefficients close to zero.

# Given that the true function is linear, linear regression benefits from low bias and stable variance, while 1-NN exhibits low bias but high variance, leading to overfitting.

# sub task 7 -----------------------------------------
set.seed(25012001)
## setup
data('diabetes', package = 'lars')

## functions
train_test_split <- function(df, train_idx){
    if (!is.data.frame(df)){
        df <- data.frame(df)
    }
    train <- df[train_idx,]
    test <- df[-train_idx,]
    return(list(train=train, test=test))
}

mse <- function(fit) mean(fit$residuals^2)

extract_significat_cols <- function(fit, alpha=.05){
    t_test_treshold <- qnorm(1-alpha/2)
    col_names <- dimnames(coef(summary(fit_1)))[[1]][-1]
    is_significant <- (abs(coef(summary(fit_1))[,'t value']) > t_test_treshold)[-1]
    return(col_names[is_significant])
}

## analysis
train_size <- 400
train_idx <- sample(1:nrow(diabetes), train_size) # nolint: seq_linter.
y_list <- train_test_split(diabetes$y, train_idx)
x_list <- train_test_split(diabetes$x, train_idx)

x_train <- x_list$train
y_train <- y_list$train
df_train <- data.frame(cbind(y_train, x_train))
x_test <- x_list$test
y_test <- y_list$test
df_test <- data.frame(cbind(y_test, x_test))

### 1: The order in which data have been collected may not be random: e.g. different groups may be taken into considerations for one study and observations may be collected per group. Not splitting the data randomly would result in training set that usese observations only form the first groups and viceversa for the test test. I.e. the two sample_ would be biased. Randomizing would ensure heterogenous groups in both sample_s.

### 2: Standardizing covariates may represent an issure in terms of interpretation of estimated coefficient: the "betas" would represent the variation in the depedent variable for one unit change in one of the covariates z-score (everything else constant). I.e. variation of "y" for one SD change in one independent variable. Without standardization, the intepretation is more straightforward (unit change in y for unit change in one x) and easy to interpret.

pairs(df_train, pch = 19)
as.dist(round(cor(df_train), digits = 3))
### 3: first column shows the relationship between the dependent variable and the covariates. BMI, MAP, LTG, GLU show remarkable degree of positive linear correlation, while HDL shows negative linear correlation. The presence of important linear relationship shows that the best functional form of the model would be a linear regression. Amonng the covariates, we can notice a very high positive linear correlation between TC and LDL. Collineraity should be take into consideration during the feature selection process to investigate causality and avoid numerical instability issues. 

fit_1 <- lm(y_train ~ ., data=df_train)
mse_1 <- mse(fit_1)
y_hat_1 <- predict(fit_1, newdata=df_test)
mse_1_test <- mean((y_hat_1 - y_test)^2)
### 4: Surprisingly, the model performs (in MSE) OOS than IS.

significant_cols <- extract_significat_cols(fit_1)
df_train_s <- data.frame(cbind(y_train, x_train[,significant_cols]))
df_test_s <- data.frame(cbind(y_test, x_test[,significant_cols]))  

fit_2 <- lm(y_train ~ ., data=df_train_s)
mse_2 <- mse(fit_2)
y_hat_2 <- predict(fit_2, newdata=df_test)
mse_2_test <- mean((y_hat_2 - y_test)^2)

F_1  <- summary(fit_1)$fstatistic
F_2  <- summary(fit_2)$fstatistic
summary(fit_1)
summary(fit_2)
### 5: The F-test is used in regression analysis to test the hypothesis that all model parameters are zero. We can see that F test rejects the null strongly in both models (more "strong" rejetction in the restricted model as expected).

# sub task 8 -----------------------------------------
library(leaps)
library(regclass)
library(tidyverse)

find_best_AIC_model <- function(df, method, sample_){
    fit_all <- regsubsets(as.formula(paste0('y_', sample_, '~.')), data=df, method=method, nvmax=10)
    results <- see_models(fit_all)
    cols <- results$Terms[which.min(results$AIC)]
    formula_str <- paste0(c('y_', sample_, '~', paste0( strsplit(cols, " ")[[1]][-1], collapse = '+ ')),collapse = '')
    fit <- lm(as.formula(formula_str), data=df)
    return(fit)
}

get_table <- function(fit,method, sample_){
    title_str = paste(method, sample_, sep='_')
    df <- data.frame(coef_name =names(fit$coef))
    df[paste('coef', title_str, sep='_')] <- unname(fit$coef)
    df[paste('ftest', title_str, sep='_')] <- pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],summary(fit)$fstatistic[3],lower.tail=FALSE)
    return(df)
}

get_analysis <- function() {
    list_of_df <- list()
    for(i in 1:2){
        df <- list(df_train, df_test)[[i]]
        sample_ <- c('train', 'test')[i]
        for(m in c('backward', 'exhaustive')){
            table <- get_table(
                    find_best_AIC_model(df, m, sample_),
                    m,
                    sample_
            )
            list_of_df <- c(list_of_df, 
                list(table)
            )
        }
    }
    return(list_of_df)
}

list_of_df <- get_analysis()
results_table<-reduce(list_of_df, full_join, by='coef_name')

### From the results table we can notice that backward-stepwise and best subset select the same 6 (constant excluded) parameters. Therefore, coefficient and F-tests are the same and F-test strongly rejects null.  The two methods select a common subset of 3 parameters also in the train sample. Again, F-test strongly rejects the null.

# sub task 9 -----------------------------------------
data("Wage", package = "ISLR")
library(olsrr)
### 1
names(Wage)
table(Wage$region)
### we can see that region is always the same
wage_filtered <- Wage |> select(-c(region, logwage))
str(Wage$education)
### education column is factor and will be handled by R as
contrasts(wage_filtered$education) <- contr.sum(levels(wage_filtered$education))
fit_9 <- lm(
    wage ~ . + I(age^2) + I(age^3),
    data = wage_filtered,
)

### 2
bs_9 <- ols_step_best_subset(fit_9)
bs_9_cols <- bs_9$metrics$predictors[which.min(bs_9$metrics$aic)]
bs_9_formula <-  paste0(
    c(
        'wage ~ ',
        paste0(strsplit(bs_9_cols, " ")[[1]], collapse = " + ")
    ),
    collapse=''
)

fit_9_bs <- lm(as.formula(bs_9_formula), data=wage_filtered)
## year: 1.2 capture the average abs. change in wage YoY (proxy for inflation)
## age: positive effect on wage
## In order: Married, Sperated, Widowed, Divorced, Never Married (baseline). In order. All the "sentimental" status diferent from Never Married bring (on avg) a positive effect on wage
## In order: White, Asian, Black, Other. Every race different from white predicts a lower avg wage.
## In order, (College Grad - Adv. Degree) (Some College - Adv. Degree), (HS - Adv. Degree), (< HS - Adv. Degree). Every effect is negative (more education -> greater salary), except for College Grad (get an Advanced degree is not convenient)
## In order, Inormation, Industrial. The first sector pays more on avg.ù
## In order, Very Good, Good. Good health predicts higher salary
## In order, Yes, No. Presence of Health Insurance predicts higher wage
## Considering that age >= 0, a neg coef on age^2 indicates that at some point an increase in age leads to a decrease in wage

### using poly(wage_filtere$age, 3, raw=T) would return a matrix in which columns are the powers from 1 to 3 of the original vector wage_filtere$age. Clearly, the regression will be equivalen to fit_9.

all(
    matrix(
    c(wage_filtered$age,
      wage_filtered$age^2,
      wage_filtered$age^3),
    byrow=F,
    ncol=3     
    ) == 
    poly(wage_filtered$age, 3, raw=T) 
)


### using poly(wage_filtere$age, 3) would return three a matrix of orthogonal vectors, thus avoiding all the problems of collinearity that my arise from fitting several powers of the same regressor. 

fit_9_ort <- lm(
    wage ~ . + poly(wage_filtered$age,3),
    data = wage_filtered|>select(-age),
)

### the estimate coefficient are different.

summary(fit_9)
summary(fit_9_ort)
### Also the significance is improved
