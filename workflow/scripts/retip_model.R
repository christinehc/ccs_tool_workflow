library(Retip)
library(hash)
library(lightgbm)

#>Starts parallel computing
prep.wizard()

##############################################################################

#> Split in training and testing using caret::createDataPartition
# load pre-defined data
data <- readxl::read_excel(toString(snakemake@input[["data"]]), sheet = "Sheet1", col_types = c("text", "text", "text", "numeric"))

#> Calculate Chemical Descriptors from CDK
descs <- getCD(data)

#> Clean dataset from NA and low variance value
db <- proc.data(descs)

set.seed(101)
i_train <- caret::createDataPartition(db$XLogP, p = 0.8, list = FALSE)
train <- db[i_train, ]
test <- db[-i_train, ]
mat_train <- as.matrix(train)
mat_test <- as.matrix(test)
col_train <- ncol(mat_train)
col_test <- ncol(mat_test)
d_train <- lightgbm::lgb.Dataset(mat_train[, 2:col_train], label = mat_train[, 1])
lightgbm::lgb.Dataset.construct(d_train)
d_test <- lightgbm::lgb.Dataset.create.valid(d_train, mat_test[, 2:col_test], label = mat_test[, 1], max_bin = 100)
valid <- list(mat_test = d_test)
params <- list(objective = "regression", metric = "rmse", max_depth = 4, max_leaf = 20)

# note: throws errors. need to access lightgbm directly
# model <- fit.lightgbm(train, test)

# train model via multiple param sets, select best one, and save
model_cv <- lightgbm::lgb.cv(params, d_train,
    nrounds = 5000,
    nfold = 10, valid, verbose = 1, early_stopping_rounds = 1000,
    record = TRUE, eval_freq = 1L, stratified = TRUE
)
best.iter <- model_cv$best_iter
params <- list(objective = "regression_l2", metric = "rmse", max_depth = 4, max_leaf = 20)
model <- lightgbm::lgb.train(params, d_train,
    nrounds = best.iter, valid, verbose = 0, early_stopping_rounds = 1000, record = TRUE, eval_freq = 1L
)
print(paste0("End training"))
saveRDS.lgb.Booster(model, snakemake@output[["model"]])

xlsx::write.xlsx(train, snakemake@output[["training"]])
