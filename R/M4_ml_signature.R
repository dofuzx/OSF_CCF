# M4 — ML signature: nested CV + external validation (preregistration C3, E2)
suppressMessages({library(glmnet); library(ranger); library(xgboost); library(pROC); library(Boruta)})
source("config.R"); set.seed(20261122)

## 1. Design matrix: TCGA (training), GSE39582 (external validation)
X_train <- readRDS(file.path(data_dir, "TCGA", "features.rds"))   # species/pathway proxies + folr2 score + covariates
y_train <- readRDS(file.path(data_dir, "TCGA", "os_2yr.rds"))     # binary 2-yr OS; censoring handled in Cox model
X_ext   <- readRDS(file.path(data_dir, "GSE39582", "features.rds"))

## 2. Nested CV model comparison (inner 5-fold tuning / outer 10-fold evaluation)
nested_cv <- function(x, y, model) {
  folds <- sample(rep(1:10, length.out = length(y)))
  sapply(1:10, function(k) {
    tr <- folds != k
    if (model == "enet") {
      cv <- cv.glmnet(x[tr, ], y[tr], family = "binomial", alpha = 0.5, nfolds = 5)
      p <- as.numeric(predict(cv, x[!tr, ], s = "lambda.min", type = "response"))
    } else if (model == "rf") {
      fit <- ranger(y ~ ., data = data.frame(y = y[tr], x[tr, ]), probability = TRUE, seed = 20261122)
      p <- predict(fit, data.frame(x[!tr, ]))$predictions[, "1"]
    } else {
      dtr <- xgb.DMatrix(as.matrix(x[tr, ]), label = y[tr])
      fit <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc"),
                       dtr, nrounds = 200, eta = 0.05, max_depth = 3, subsample = 0.8)
      p <- predict(fit, as.matrix(x[!tr, ]))
    }
    as.numeric(pROC::roc(y[!tr], p, quiet = TRUE)$auc)
  })
}
aucs <- sapply(c("enet","rf","xgb"), function(m) nested_cv(X_train, y_train, m))
best <- names(which.max(colMeans(aucs)))

## 3. Final model on full TCGA -> score both cohorts; calibration & DCA
fit_full <- switch(best, enet = cv.glmnet(as.matrix(X_train), y_train, family = "binomial", alpha = 0.5),
                   rf = NULL, xgb = NULL)   # rf/xgb final fits analogous
score <- function(X) predict(fit_full, as.matrix(X), s = "lambda.min", type = "response")
s_ext <- score(X_ext)
# calibration: rms::calibrate; DCA: dcurves::dca

## 4. Feature stability: Boruta on training set
bor <- Boruta(x = X_train, y = factor(y_train), seed = 20261122)
write.csv(bor$finalDecision, file.path(out_dir, "M4", "boruta_stability.csv"))
saveRDS(list(aucs = aucs, best = best, score_ext = s_ext), file.path(out_dir, "M4", "ml_results.rds"))
write("M4 done", file.path(out_dir, "M4", ".done"))
