###################################################
### Example: Prostate Cancer
###################################################
if (!requireNamespace("ElemStatLearn")) {
    URL <- "https://cran.r-project.org/src/contrib/Archive/ElemStatLearn"
    install.packages(file.path(URL, "ElemStatLearn_2015.6.26.2.tar.gz"))
}

data("prostate", package = "ElemStatLearn")
prostate.test <- prostate |>
  dplyr::filter(!train) |>
  dplyr::select(-train)
prostate <- prostate |> 
  dplyr::filter(train) |>
  dplyr::select(-train)


## Exploratory analysis
pairs(prostate, pch = 19)
as.dist(round(cor(prostate), digits = 3))

## Standardization based on the training data
prostate.scaled <-
  scale(prostate |> dplyr::select(-lpsa)) 
prostate.test.scaled <-
  scale(prostate.test |> dplyr::select(-lpsa), 
        center = attr(prostate.scaled, "scaled:center"),
        scale = attr(prostate.scaled, "scaled:scale"))
prostate.scaled <- prostate.scaled |> 
  cbind(data.frame(lpsa = prostate$lpsa))
prostate.test.scaled <- prostate.test.scaled |>
  cbind(data.frame(lpsa = prostate.test$lpsa))

## Linear model: full
linear.model <- lm(lpsa ~ ., data = prostate.scaled, y = TRUE)
printCoefmat(round(coef(summary(linear.model)), digits = 2))

## Linear model: subset
linear.model.sub <- lm(lpsa ~ . - age - lcp - gleason - pgg45, 
  data = prostate.scaled)
anova(linear.model.sub, linear.model)  

## Linear model: null model
null.model <- lm(lpsa ~ 1, data = prostate.scaled)

## Model comparison
prediction.error <- sapply(list("base model" = null.model, 
                                "full model" = linear.model, 
				"subset model" = linear.model.sub), 
  function(model) {
    yhat <- predict(model, newdata = prostate.test.scaled)
    y <- prostate.test$lpsa
    mean((y - yhat)^2)
})
round(prediction.error, digits = 3)

## Exhaustive best subset selection
par(mar = c(5.1, 4.1, 2.1, 0))
leaps.results <- summary(leaps::regsubsets(lpsa ~ ., data = prostate, 
     nbest = 70, nvmax = 9, really.big = TRUE))
plot(c(1, rowSums(leaps.results$which)),
     c(leaps.results$obj$nullrss, leaps.results$rss),
     pch = 19, col = "grey", ylim = c(-2, 102),
     main = "Exhaustive best subset selection",
     xlab = "Subset Size k", ylab = "Residual Sum-of-Squares")
lines(1:9,
      c(leaps.results$obj$nullrss,
        with(leaps.results, tapply(rss, rowSums(which), min))),
      pch = 19, type = "b", col = "red")      

## Ridge regression
par(mar = c(5.1, 4.1, 0.1, 0.1))
df <- seq(0, 8, length.out = 20)
X <- sqrt(nrow(prostate) / (nrow(prostate) - 1)) * model.matrix(linear.model)[, -1]
d <- svd(X)$d
lambda <- c(Inf, sapply(df[-1], function(y) 
  uniroot(function(x) sum(d^2 / (d^2 + x)) - y, c(0, 10^6),
          extendInt = "yes")$root))
ridge.models <- MASS::lm.ridge(lpsa ~ ., data = prostate.scaled,
  lambda = lambda)
matplot(df, coef(ridge.models)[, -1], 
  xlab = expression(paste("df(", lambda, ")")), 
  ylab = "Coefficients", type = "b", pch = 19, lty = 1, 
  xlim = c(0, 8.5))
text(8, coef(ridge.models)[nrow(coef(ridge.models)), -1],
    labels = colnames(coef(ridge.models))[-1], pos = 4)
abline(h = 0, lty = 2, col = "grey")

## Lasso fit
par(mar = c(5.1, 4.1, 2.1, 0.1))
lasso.fit <- lars::lars(model.matrix(linear.model)[, -1], linear.model$y, 
  normalize = FALSE)
plot(lasso.fit, lwd = 2)

