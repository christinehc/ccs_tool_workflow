library(Retip)
library(lightgbm())
library(compiler)
getCD <- cmpfun(getCD)


# Define new chem desc function
getCD_memfix <- function(x) {
  gc(full = TRUE)

  print(paste0("Converting SMILES..."))

  for (i in 1:nrow(x)) {
    smi <- rcdk::parse.smiles(as.character(unlist(x[i, "SMILES"])))[[1]]
    smi1 <- rcdk::generate.2d.coordinates(smi)
    smi1 <- rcdk::get.smiles(smi, rcdk::smiles.flavors(c("CxSmiles")))
    x$SMILES[i] <- smi1
    print(paste0(i, " of ", nrow(x)))
  }
  smi <- NULL
  smi1 <- NULL
  gc()

  # select only one descriptors. This helps to remove compounds that makes errors
  descNames1 <- c("org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor")

  print(paste0("Checking for compound errors..."))

  # calculate only 1 descriptor for all the molecules
  mols_x <- rcdk::parse.smiles(as.character(unlist(x[1, "SMILES"])))
  descs1_x <- rcdk::eval.desc(mols_x, descNames1)
  mols_x <- NULL

  for (i in 2:nrow(x)) {
    mols1 <- rcdk::parse.smiles(as.character(unlist(x[i, "SMILES"])))
    descs1_x[i, ] <- rcdk::eval.desc(mols1, descNames1)
    print(paste0(i, " of ", nrow(x)))
  }
  mols1 <- NULL
  descNames1 <- NULL
  gc()

  # remove molecules that have NA values with only one descriptor
  x_na <- data.frame(descs1_x, x)
  descs1_x <- NULL

  x_na_rem <- x_na[stats::complete.cases(x_na), ]
  x_na <- NULL
  x_na_rem <- x_na_rem[, -c(1:6)]
  gc()

  # computing the whole descriptos on the good on the clean dataset
  print(paste0("Computing Chemical Descriptors 1 of ", nrow(x_na_rem), " ... Please wait"))

  # select all possible descriptors
  descNames <- rcdk::get.desc.names(type = "all")

  mols_x1 <- rcdk::parse.smiles(as.character(unlist(x_na_rem[1, "SMILES"])))[[1]]
  rcdk::convert.implicit.to.explicit(mols_x1)
  descs_x_loop <- rcdk::eval.desc(mols_x1, descNames)
  mols_x1 <- NULL
  gc()

  for (i in 2:nrow(x_na_rem)) {
    mols <- rcdk::parse.smiles(as.character(unlist(x_na_rem[i, "SMILES"])))[[1]]
    rcdk::convert.implicit.to.explicit(mols)
    descs_x_loop[i, ] <- rcdk::eval.desc(mols, descNames)
    print(paste0(i, " of ", nrow(x_na_rem)))
  }
  mols <- NULL
  descNames <- NULL
  gc()
  datadesc <- data.frame(x_na_rem, descs_x_loop)
  x_na_rem <- NULL
  descs_x_loop <- NULL
  gc()
  return(datadesc)
}



#>Starts parallel computing
prep.wizard()

data <- readxl::read_excel(snakemake@input[["data"]], col_types = c("text", "text", "text", "text"))
train <- readxl::read_excel(snakemake@input[["training"]])
model <- readRDS.lgb.Booster(toString(snakemake@input[["model"]]))

#> compute Chemical descriptors
desc <- getCD(data)

#> perform the RT spell
pred <- RT.spell(training = train, target = desc, model = model)

# save results
xlsx::write.xlsx(pred, snakemake@output[["result"]])
