library(lava)

get_lava_model <- function(time_horizon, coefs) {
  m <- lvm()
  recursive_function <- function(z, name, time_varying) {
    if (inherits(z, "interaction")) {
      y <- m
      message("adding: ", name)
      first_part <- z[[1]]
      second_part <- z[[2]]
      transform(y, formula = formula(paste0(name, "~", first_part, "+", second_part))) <- function(x) x[1] * x[2]
      m <<- y
    } else if (inherits(z, "normal")) {
      y <- m
      seqq <- seq_along(z)
      message("adding: ", name)
      intercept_loc <- grep("Intercept|mean", names(z))
      sd_loc <- grep("sd", names(z))
      distribution(y, name) <- normal.lvm(mean = z[intercept_loc], sd = z[sd_loc])
      index <- setdiff(seqq, c(intercept_loc, sd_loc))
      if (length(index) > 0) {
        regression(y,
          to = name,
          from = unlist(lapply(names(z)[index], function(y) {
            gsub("-", " - ", gsub(" ", "", y))
          }))
        ) <- z[index]
      }
      m <<- y
    } else if (is.numeric(z)) {
      y <- m
      message("adding: ", name)
      if (time_varying) {
        distribution(y, name) <- binomial.lvm(link = "logit")
        intercept(y, name) <- z[1]
        if (length(z) > 1) {
          regression(y,
            to = name,
            from = unlist(lapply(names(z)[-1], function(y) {
              gsub("-", " - ", gsub(" ", "", y))
            }))
          ) <- z[-1]
        }
      } else {
        ## this is slow when simulating
        levs <- names(z)
        message("adding: ", name)
        y <- categorical(y, name, labels = levs, p = z[-length(z)])
        eval(parse(text = paste0(
          "transform(y,formula(", paste0(name, gsub(" ", "", levs), "~", name),
          ")) <- function(x){1*(x[[1]]=='", levs, "')}"
        )))
      }
      m <<- y
    } else if (is.list(z)) {
      time_varying <- FALSE
      if (any(grepl("_\\d+$", names(z)))) { # time-varying
        time_varying <- TRUE
        z <- z[grepl(paste0("_", 0:time_horizon, "$", collapse = "|"), names(z))]
      }
      lapply(seq_along(z), function(i) recursive_function(z[[i]], names(z)[i], time_varying))
    } else {
      return(NULL) # Return NULL if it's not a numeric or list
    }
  }
  recursive_function(z = coefs, name = NULL, time_varying = FALSE)
  m
}

get_interventional_distribution <- function(x, name_regimen, abar) {
  if (is.list(abar)) {
    stop("abar must be a vector")
  }
  vars_to_use_reg <- paste0(name_regimen, "_", 0:(length(abar) - 1))
  for (i in seq_along(abar)) {
    distribution(x, vars_to_use_reg[i]) <- Binary.lvm(abar[i])
  }
  x
}

get_sim_data <- function(lava_model, sample_size, separate= FALSE, simdata=NULL){
  if (is.null(simdata)){
    simdata = as.data.table(sim(lava_model$model, sample_size))
  }
  simdata[, pnr := seq_len(.N)] ## Include pnr in simdata
  no_censor <- no_comp.risk <- no_time_covariates <- FALSE
  if (is.null(lava_model$censoring)){
    lava_model$censoring <- "C"
    ## add C_1, ..., C_{time_horizon} to simdata as ones
    simdata[, paste0(lava_model$censoring,"_",1:lava_model$time_horizon) := 1]
    no_censor <- TRUE
  }
  if (is.null(lava_model$comp.event)){
    lava_model$comp.event <- "D"
    ## add D_1, ..., D_{time_horizon-1} to simdata as zeros
    simdata[, paste0(lava_model$comp.event,"_",1:(lava_model$time_horizon-1)) := 0]
    no_comp.risk <- TRUE
  }
  if  (is.null(lava_model$time_covariates)){
    lava_model$time_covariates <- "L"
    ## add L_1, ..., L_{time_horizon-1} to simdata as zeros
    simdata[, paste0(lava_model$time_covariates,"_",0:(lava_model$time_horizon-1)) := 0]
    no_time_covariates <- TRUE
  }
  time_vars <- c(lava_model$censoring, lava_model$outcome, lava_model$comp.event, lava_model$time_covariates, lava_model$treatment)
  time_vars <- paste(rep(time_vars, lava_model$time_horizon), rep(1:lava_model$time_horizon, 1, each=length(time_vars)), sep = "_")  
  keep_names <- c("pnr",lava_model$baseline, paste0(c(lava_model$time_covariates,lava_model$treatment),"_0"), time_vars)
  ## remove Dead and time_covariates from last time point
  keep_names <- keep_names[!grepl(paste0(lava_model$comp.event,"_",lava_model$time_horizon,"|",lava_model$treatment,"_",lava_model$time_horizon,"|",paste0(paste0(lava_model$time_covariates,"_",lava_model$time_horizon),collapse="|")),keep_names)]
  simdata <- simdata[, ..keep_names]
  
  Y_nodes_position = match(paste0(lava_model$outcome,"_",1:lava_model$time_horizon), names(simdata))
  C_nodes_position = match(paste0(lava_model$censoring,"_",1:lava_model$time_horizon), names(simdata))
  if (lava_model$time_horizon>1){
    D_nodes_position = match(paste0(lava_model$comp.event,"_",1:(lava_model$time_horizon-1)), names(simdata))
  }
  else {
    D_nodes_position = integer(0)
  }
  L_nodes_position = list()
  for (ll in lava_model$time_covariates){
    if (lava_model$time_horizon>1){
      L_nodes_position[[ll]] = match(paste0(ll,"_",1:(lava_model$time_horizon-1)), names(simdata))
    }
    else {
      L_nodes_position[[ll]] = integer(0)
    } 
  }
  if (lava_model$time_horizon>1){
    A_nodes_position = match(paste0(lava_model$treatment,"_",1:(lava_model$time_horizon-1)), names(simdata))
  }
  else {
    A_nodes_position = integer(0)
  }
  # avoid repetition
  get_later_nodes <- function(D_nodes_position,C_nodes_position,Y_nodes_position,L_nodes_position,A_nodes_position,lava_model,k){
    increasing_seq <- function(from, to, ...) {
      if (length(from) == 0 || length(to) == 0 ||from > to) {
        return(NULL)
      } else {
        return(seq(from, to, ...))
      }
    }
    
    later_D_nodes=intersect(increasing_seq(k+1, D_nodes_position[length(D_nodes_position)]),D_nodes_position)
    later_C_nodes=intersect(increasing_seq(k+1, C_nodes_position[length(C_nodes_position)]),C_nodes_position)
    later_Y_nodes=intersect(increasing_seq(k+1, Y_nodes_position[length(Y_nodes_position)]),Y_nodes_position)
    later_A_nodes=intersect(increasing_seq(k+1, A_nodes_position[length(A_nodes_position)]),A_nodes_position)
    later_L_nodes=list()
    for (ll in lava_model$time_covariates){
      later_L_nodes[[ll]] = intersect(increasing_seq(k+1, L_nodes_position[[ll]][length(L_nodes_position[[ll]])]),L_nodes_position[[ll]])
    }
    list(later_D_nodes=later_D_nodes,later_C_nodes=later_C_nodes,later_Y_nodes=later_Y_nodes,later_L_nodes=later_L_nodes,later_A_nodes=later_A_nodes)
  }
  
  ## Make sure that no censoring or death happens after event
  for(k in Y_nodes_position){
    list2env(get_later_nodes(D_nodes_position,C_nodes_position,Y_nodes_position,L_nodes_position,A_nodes_position,lava_model,k),environment())
    if(any(has_outcome <- (simdata[[k]]%in%1))){
      for(l in later_D_nodes) {set(simdata,j=l,i=which(has_outcome),value=NA)}
      for(l in later_C_nodes) {set(simdata,j=l,i=which(has_outcome),value=1)} ## Because 1 means uncensored
      for(l in later_Y_nodes) {set(simdata,j=l,i=which(has_outcome),value=1)}
      for (ll in lava_model$time_covariates){
        for(l in later_L_nodes[[ll]]) {set(simdata,j=l,i=which(has_outcome),value=NA)}
      }
      for(l in later_A_nodes) {set(simdata,j=l,i=which(has_outcome),value=NA)}
    }
  }
  ## Make sure that no event or censoring happens after death
  if(length(lava_model$comp.event)>0){
    for(k in D_nodes_position){
      list2env(get_later_nodes(D_nodes_position,C_nodes_position,Y_nodes_position,L_nodes_position,A_nodes_position,lava_model,k),environment())
      if(any(has_died <- (simdata[[k]]%in%1))){
        for(l in later_Y_nodes) {set(simdata,j=l,i=which(has_died),value=0)}
        for(l in later_C_nodes) {set(simdata,j=l,i=which(has_died),value=1)} ## Because 1 means uncensored
        for(l in later_D_nodes) {set(simdata,j=l,i=which(has_died),value=1)}
        for (ll in lava_model$time_covariates){
          for(l in later_L_nodes[[ll]]) {set(simdata,j=l,i=which(has_died),value=NA)}
        }
        for(l in later_A_nodes) {set(simdata,j=l,i=which(has_died),value=NA)}
      }
    }
  }
  ## Make sure that no event or death happens after censoring
  if(length(lava_model$censoring)>0){ ## Make sure that no event happens after death
    for(k in C_nodes_position){
      list2env(get_later_nodes(D_nodes_position,C_nodes_position,Y_nodes_position,L_nodes_position,A_nodes_position,lava_model,k),environment())
      if(any(is_censored <- (simdata[[k]]%in%0))){
        for(l in later_Y_nodes) {set(simdata,j=l,i=which(is_censored),value=NA)}
        for(l in later_D_nodes) {set(simdata,j=l,i=which(is_censored),value=NA)}
        for(l in later_C_nodes) {set(simdata,j=l,i=which(is_censored),value=0)}
        for (ll in lava_model$time_covariates){
          for(l in later_L_nodes[[ll]]) {set(simdata,j=l,i=which(is_censored),value=NA)}
        }
        for(l in later_A_nodes) {set(simdata,j=l,i=which(is_censored),value=NA)}
      }
    }
  }
  if (separate) {
    outcome_vars_pos <- c(1,Y_nodes_position, D_nodes_position, C_nodes_position)
    outcome <- simdata[,..outcome_vars_pos]
    if (no_censor){
      outcome[,paste0(lava_model$censoring,"_",1:lava_model$time_horizon) := NULL]
    }
    if (no_comp.risk){
      outcome[,paste0(lava_model$comp.event,"_",1:(lava_model$time_horizon-1)) := NULL]
    }
    regimen_vars_pos_0 <- match(paste0(lava_model$treatment,"_0"), names(simdata))
    regimen_vars_pos <- c(1,regimen_vars_pos_0,A_nodes_position)
    regimen <- simdata[,..regimen_vars_pos]
    if (!no_time_covariates){
      time_covariates_vars_pos0 <- match(paste0(lava_model$time_covariates,"_0"), names(simdata))
      time_covariates_vars_pos <- c(1,time_covariates_vars_pos0,unlist(L_nodes_position))
      time_covariates <- simdata[,..time_covariates_vars_pos]
    }
    else {
      time_covariates <- NULL
    }
    baseline_vars_pos <- c(1,match(lava_model$baseline, names(simdata)))
    baseline_covariates <- simdata[,..baseline_vars_pos]
    return(list(outcome=outcome, regimen=regimen, time_covariates=time_covariates, baseline_covariates=baseline_covariates))
  } else {
    if (no_censor){
      simdata[, paste0(lava_model$censoring,"_",1:lava_model$time_horizon) := NULL]
    }
    if (no_comp.risk){
      simdata[, paste0(lava_model$comp.event,"_",1:(lava_model$time_horizon-1)) := NULL]
    }
    if  (no_time_covariates){
      simdata[, paste0(lava_model$time_covariates,"_",0:(lava_model$time_horizon-1)) := NULL]
    }
    return(simdata)
  }
}
