###################################################
### Example: South African Heart Disease
###################################################

data("SAheart", package = "ElemStatLearn")
SAheart <- SAheart |>
    dplyr::select(chd, sbp, tobacco, ldl, famhist, obesity, alcohol, age) |>
    dplyr::mutate(chd = factor(chd))

## Exploratory analysis
pairs(SAheart[, -1], pch = 19, col = SAheart$chd)
spineplot(chd ~ age, data = SAheart)

## GLM with a single covariate
model1 <- glm(chd ~ age, data = SAheart,
  family = binomial())
summary(model1)

library("effects")
Effect("age", model1)
eff1 <- allEffects(model1)
plot(eff1)

boxplot(fitted(model1, type = "response") ~ chd, data = SAheart,
  ylab = "Predicted probability", varwidth = TRUE, ylim = 0:1)

## GLM with two covariates
model2 <- glm(chd ~ age + famhist, data = SAheart, 
  family = binomial())
model2
summary(model2)

## Nested model comparison 
anova(model1, model2)

eff2 <- allEffects(model2)
eff2
plot(eff2)

boxplot(fitted(model2, type = "response") ~ chd, data = SAheart,
  ylab = "Predicted probability", varwidth = TRUE, ylim = 0:1)

## (Quasi-)complete separation
SAheart0 <- subset(SAheart, 
  (age <= 50 & chd == 0) | (age >= 50 & chd == 1))
model0 <- glm(chd ~ famhist + age, 
  data = SAheart0, family = binomial())
summary(model0)

## Full model
full.model <- glm(chd ~ ., data = SAheart,
  family = binomial())
printCoefmat(round(coef(summary(full.model)),
  digits = 2))

## Stepwise selected model
step.model <- step(full.model, trace = 0)
printCoefmat(round(coef(summary(step.model)),
  digits = 2))

## Regularized GLM 
mf <- model.frame(chd ~ ., data = SAheart)
X <- model.matrix(chd ~ ., data = mf)[, -1]
X <- scale(X)
y <- model.response(mf)

## LASSO
library("glmnet")
model <- glmnet(X, y, family = "binomial",
  nlambda = 500)

par(mfrow = c(1, 2))
plot(model, xlim = c(0, 2.8), col = 1:8)
text(2, coef(model, s = 0)[-1, "s1"],
  rownames(coef(model, s = 0))[-1], 
  col = 1:8, adj = 0)
plot(model, xvar = "lambda", xlim = c(-7.9, -1.5),
     col = 1:8)
text(-5.9, coef(model, s = 0)[-1, "s1"],
  rownames(coef(model, s = 0))[-1], 
  col = 1:8, adj = 1)

## Ridge
model <- glmnet(X, y, family = "binomial", 
  alpha = 0, nlambda = 500)

par(mfrow = c(1, 2))
plot(model, xlim = c(0, 2.5), col = 1:8)
text(1.95, coef(model, s = 0)[-1, "s1"],
  rownames(coef(model, s = 0))[-1], 
  col = 1:8, adj = 0)
plot(model, xvar = "lambda", xlim = c(-8, 5),
     col = 1:8)
text(-4.25, coef(model, s = 0)[-1, "s1"],
  rownames(coef(model, s = 0))[-1], 
  col = 1:8, adj = 1)


