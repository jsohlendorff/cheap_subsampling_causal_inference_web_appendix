library(data.table)
library(foreach)
library(lava)
library(targets)
source("get_lava_model.R")
source("bootstrap_ci.R")
source("~/cheap_subsampling_causal_inference_web_appendix/data_generating_mechanism.R")
tar_source("Ltmle") ## use targets to source all of the files
## see coefs
coefs

## make lava model from coefficients
lava_model <- get_lava_model(3, coefs)

## set up model object
mm <- list()
mm$outcome <- "heart_failure"
mm$censoring <- "Censored"
mm$comp.event <- "Dead"
mm$treatment <- "GS"
mm$time_covariates <- "insulin"
mm$time_horizon <- time_endpoint
mm$baseline <- c("sex","education","agegroups","tertile_income","index_heart_failure","diabetes_duration", "secondline_duration","first_2ndline")
mm$model <- lava_model

## simulate data

## get Ltmle data
set.seed(6)
data <- get_sim_data(mm, 1000)

## get data for use in Ltmle, i.., split up separately into the outcomes, regimen, baseline covariates and time-varying covariates
set.seed(6)
data <- get_sim_data(mm, 1000, TRUE)

x <- prepare_Ltmle(
  outcome_data = data$outcome,
  regimen_data = data$regimen,
  baseline_data = data$baseline_covariates,
  timevar_data = data$time_covariates,
  time_horizon = 3,
  censored_label = 0,
  name_outcome = "heart_failure",
  name_regimen = "GS",
  name_censoring = "Censored",
  name_competing_risk = "Dead",
  abar = list(treat=c(1,1,1), control = c(0,0,0)),
  SL.library = "glm",
  verbose = TRUE,
  gbounds = c(0,1)
)
f<-summary(do.call("Ltmle", x))
f

bs <- 5
## cheap_subsampling_ci
k_m <- 0.632
res_subsampling <- list()
for (b in seq_len(bs)) {
  ## subsample data of size m
  subsample <- sample(1:nrow(data$outcome), size = floor(k_m * nrow(data$outcome)), replace = FALSE)
  formatted_data_sub <-
    lapply(data, function(x) {
      x[subsample, ]
    })
  x <- prepare_Ltmle(
    outcome_data = formatted_data_sub$outcome,
    regimen_data = formatted_data_sub$regimen,
    baseline_data = formatted_data_sub$baseline_covariates,
    timevar_data = formatted_data_sub$time_covariates,
    time_horizon = 3,
    censored_label = 0,
    name_outcome = "heart_failure",
    name_regimen = "GS",
    name_censoring = "Censored",
    name_competing_risk = "Dead",
    abar = list(treat=c(1,1,1), control = c(0,0,0)),
    SL.library = "glm",
    verbose = TRUE,
    gbounds = c(0,1)
  )
  f_temp <- do.call("Ltmle", x)
  res_subsampling[[b]] <- summary(f_temp)[,c(1:3)]
}

res_subsampling <- rbindlist(res_subsampling)

## look at causal contrast
get_cheap_subsampling_ci(f[Target_parameter== "ATE", estimate], res_subsampling[Target_parameter == "ATE", estimate], floor(k_m * nrow(data$outcome)), nrow(data$outcome), 0.05)

res_non_parametric_bootstrap <- list()
for (b in seq_len(bs)) {
  ## subsample data of size m
  bootstrap_sample <- sample(1:nrow(data$outcome), size = nrow(data$outcome), replace = TRUE)
  formatted_data_boot <-
    lapply(data, function(x) {
      temp <- x[bootstrap_sample, ]
      temp[, pnr:= 1:.N]
      temp
    })
  x <- prepare_Ltmle(
    outcome_data = formatted_data_boot$outcome,
    regimen_data = formatted_data_boot$regimen,
    baseline_data = formatted_data_boot$baseline_covariates,
    timevar_data = formatted_data_boot$time_covariates,
    time_horizon = 3,
    censored_label = 0,
    name_outcome = "heart_failure",
    name_regimen = "GS",
    name_censoring = "Censored",
    name_competing_risk = "Dead",
    abar = list(treat=c(1,1,1), control = c(0,0,0)),
    SL.library = "glm",
    verbose = TRUE,
    gbounds = c(0,1)
  )
  f_temp <- do.call("Ltmle", x)
  res_non_parametric_bootstrap[[b]] <- summary(f_temp)[,c(1:3)]
}

## look at causal contrast
res_non_parametric_bootstrap <- rbindlist(res_non_parametric_bootstrap)
get_cheap_bootstrap_ci(f[Target_parameter== "ATE", estimate], res_non_parametric_bootstrap[Target_parameter == "ATE", estimate], nrow(outcome$data), 0.05)