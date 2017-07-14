#'---
#' title: Simple glmnet model
#' author: Žiga Avsec
#' wb:
#'   input: ["data/Concise/Splice_branchpoints/processed/branchpointer/train/filteredDescr.csv",
#'            "data/Concise/Splice_branchpoints/processed/branchpointer/test/filteredDescr.csv"]
#'   output: ["data/Concise/Splice_branchpoints/test_predictions/glmnet.csv"]
#'---
dt <- fread("data/Concise/Splice_branchpoints/processed/branchpointer/train/filteredDescr.csv")
dt_test <- fread("data/Concise/Splice_branchpoints/processed/branchpointer/train/filteredDescr.csv")

dt[, V1 := NULL]
dt_test[, V1 := NULL]

dt %>% head(3)

x_train <- model.matrix(set~. - 1, dt)
x_test <- model.matrix(set~.-1, dt_test)
y_train <- dt[, as.factor(set)]
y_test <- dt_test[, as.factor(set)]

## center and scale - shouldn't do for the other features
pp <- preProcess(x_train, method = c("center", "scale"))
x_train <- predict(pp, x_train)
x_test <- predict(pp, x_test)

x_train <- as(x_train, "sparseMatrix")
x_test <- as(x_test, "sparseMatrix")

library(doMC)
registerDoMC(cores=5)
fit <- cv.glmnet(x_train, as.integer(y_train == "HC"),
                 family = "binomial", alpha = .5, nfolds=5)
plot(fit)
## lambda=0 if for sure the right choice (large dataset)
fit <- glmnet(x_train, as.integer(y_train == "HC"), family = "binomial", alpha = .5, lambda = 0)
plot(fit)
y_test_pred <- predict(fit, x_test)[,1]
## AUC is .90
confusionMatrix(as.integer(y_test_pred>0), as.integer(y_test=="HC"))

## get the ID
dt_pred <- data.table(y_true = y_test, y_pred = y_test_pred)
write_csv(dt_pred, "/s/project/deepcis/Concise/Splice_branchpoints/test_predictions/glmnet.csv")

