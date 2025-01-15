library(tidyverse)
library(activity)
library(camtraptor)
library(camtrapDensity)

# Von Mises mixture probability density

# INPUT
# x: vector of values at which to evaluate
# prm: named vector of distribution parameters;
#      names ki, mi, pi respectively k, m, p for the ith component,
#      p with 1 fewer than k/m, last is inferred since sum(pi) = 1

# OUTPUT
# A vector of probability densities
dvmsmix <- function(x, prm){
  k <- exp(prm[grep("k",names(prm))])
  m <- prm[grep("m",names(prm))]
  p <- prm[grep("p",names(prm))]
  if(length(k)==1) p  <- 1 else p <- c(1-sum(p), p)
  res <- rep(0, length(x))
  for(i in 1:length(k))
    res <- res + p[i]*activity::dvonm(x, m[i], k[i])
  res
}

# Generate random radian time observations from an empirical distribution function

# INPUT
# n: number of observations to generate
# fit: dataframe describing the empirical distribution, with columns:
#      x: a regular sequence of times from 0 to 2*pi
#      y: the (relative) probabilities of observing times

# OUTPUT
# A vector of radian times
rand_time <- function (n, fit){
  if (sum(c("x", "y") %in% names(fit)) != 2) 
    stop("fit must be a dataframe with (at least) columns named x and y")
  if (diff(range(diff(fit$x))) > 1e-04) 
    stop("x doesn't seem to be a regular sequence")
  cdf <- c(0, cumsum(head(fit$y, -1)) / sum(head(fit$y, -1)))
  rn <- stats::runif(n)
  res <- stats::approx(cdf, fit$x, rn)$y
  res <- res - mean(fit$x[1:2])
  res <- ifelse(res<0, res+2*pi, res)
  res
}

# Make a time x stratum matrix from a table of parameter values
#
# INPUT
#   table: a dataframe of parameter values with (at least) fields:
#     est: parameter estimates
#     timeID: (otional) integer time period IDs
#     stratumID: (optional) stratum IDs
#   str_data: dataframe of stratum data with (at least) a stratumID field
#   time_thresholds: vector of radian times demarcating time periods
#   nt: number of time points to generate for matrices
# OUTPUT
#   A matrix with dimensions nt rows and nrow(str_data) columns.
#
# DETAILS
# parameters$timeID values must be in 1:length(time_thresholds)
# parameters$stratumID values must match str_data$stratumID
# Each row of parameters must represent a unique combination of timeID 
# and stratumID (if both are present).
# time_thresholds must be strictly increasing, defining period 1 as
# tt[1]< <=tt[2] and so on.
pmat_from_table <- function(table, str_data, time_thresholds=NULL, nt=513){
  if(!"est" %in% names(table)) stop("table must contain est column")
  if("stratumID" %in% names(table) &
     !all(str_data$stratumID %in% unique(table$stratumID)))
    stop("Not all str_data$stratumID can be found in table$stratumID")
  if("timeID" %in% names(table)){
    if(is.null(time_thresholds)) 
      stop("time_thresholds must be provided if table contains timeID")
    if(min(time_thresholds)<0 | 
       max(time_thresholds)>2*pi | 
       any(diff(time_thresholds) <= 0))
      stop("time_thresholds must be an increasing sequence of radian times between 0 and 2*pi")
    tID <- unique(table$timeID)
    ntt <- length(time_thresholds)
    if(length(tID) != ntt | !all(range(tID)==c(1,ntt)))
      stop("table$timeID must contain integers from 1 to length(time_threshold)")
  }
  
  ni <- length(time_thresholds)
  ns <- nrow(str_data)
  sq <- seq(0, 2*pi, len=nt)
  
  if("timeID" %in% names(table)){
    i <- expand.grid(sq, time_thresholds) %>%
      apply(1, function(x) sign(diff(x))) %>%
      matrix(ncol=ni) %>%
      apply(1, sum) %>%
      (function(x) (ni-x)/2) %>%
      floor() %>%
      (function(x){
        x[x==0] <- ni
        x
      })
    if("stratumID" %in% names(table)){
      i <- i %>%
        expand.grid(str_data$stratumID) %>%
        apply(1, paste, collapse="") %>%
        match(paste0(table$timeID, table$stratumID))
    } else
      i <- match(i, table$timeID)
  } else{
    if("stratumID" %in% names(table)){
      i <- match(str_data$stratumID, table$stratumID) %>%
        rep(each=nt)
    } else
      i <- 1
  }
  res <- matrix(table$est[i], nrow=nt, ncol=ns)
  colnames(res) <- str_data$stratumID
  res
}

# Make speed/radius/angle parameter array from a function

# INPUT
# str_data: dataframe with a row per stratum and stratumID column
# func: a function with arguments time and/or stratum (or neither), used to generate time- and/or stratum-specific parameter values

# OUTPUT
# Returns a 513 times by strata matrix of parameter values
pmat_from_func <- function(func, str_data, nt=513){
  if(!any(class(func) %in% "function"))
    stop("func must be a function")
  args <- formalArgs(func)
  if(!all(args %in% c("time", "stratum")))
    stop("Function func cannot have arguments other than time and stratum")
  if(!"stratumID" %in% names(str_data)) 
    stop("str_data must contain a stratumID column")  
  
  narg <- length(args)
  nstr <- nrow(str_data)
  sq <- seq(0, 2*pi, len=nt)
  res <- if(narg==0) func() else
    if(narg == 1){
      if(args=="time") func(sq) else
        if(args=="stratum") rep(func(str_data$stratumID), each=nt)
    } else
      func(time = rep(sq, nstr),
           stratum = rep(str_data$stratumID, each=nt))
  res <- matrix(res, nrow=nt, ncol=nstr)
  colnames(res) <- str_data$stratumID
  res
}

# Generate randomised dataset

# INPUT
# dep_dat: dataframe of deployment data, columns deploymentID, stratumID, effort
# str_dat: dataframe of stratum data, columns stratumID, area
# N: scalar population size
# a: matrix of activity levels, eval times (wrapped) by strata
# p: matrix of proportion in stratum, eval times (wrapped) by strata
# sfunc/rfunc/afunc: NULL or a function with arguments time and stratum, passed to make_parmat
# size: negative binomial size parameter for random observations per deployment generation
# NB strata in above vectors/matrices must occur in a consistent order (there is no name matching)

# OUTPUT
# List holding input dep_dat and str_dat plus obs_dat, a dataframe of
# observations with columns deploymentID and time (radian time of day of 
# each observation)
data_gen <- function(dep_dat, str_dat, 
                     N, a, p, 
                     sfunc=NULL, rfunc=NULL, afunc=NULL,
                     size=10){
  nd <- nrow(dep_dat)
  ns <- nrow(str_dat) # number of strata
  
  spd <- pmat_from_func(sfunc, str_dat, 513)
  rad <- pmat_from_func(rfunc, str_dat, 513)
  ang <- pmat_from_func(afunc, str_dat, 513)
  q <- spd * rad * (2+ang)
  if(length(q) > 1) q <- t(q)  # parameter weight, stratum x time
  D <- t(p) * N / str_dat$area # density, stratum x time
  E <- dep_dat %>%
    group_by(stratumID) %>%
    summarise(E = sum(effort)) %>%
    pull(E) # effort, stratum
  P <- q * t(a) * D * E / (512 * pi) # expected obs, stratum x time
  ncams <- with(dep_dat, tapply(stratumID, stratumID, length)) # cameras, stratum
  Pdep <- (apply(P[,-1], 1, sum) / ncams)[dep_dat$stratumID] # expected observations, deployment
  Ydep <- rnbinom(nd, mu=Pdep, size=size) # realised observations, deployment
  Y <- tapply(Ydep, dep_dat$stratumID, sum) # realised observations, stratum
  
  sq <- (0:512)*2*pi/512 # generate random observation times
  f <- function(i){
    df <- data.frame(x=sq, y=P[i,])
    rand_time(Y[i], df)
  }
  obs_dat <- data.frame(deploymentID = rep(dep_dat$deploymentID, Ydep),
                        time = unlist(lapply(1:ns, f)))
  list(dep=dep_dat, str=str_dat, obs=obs_dat)
}

# Checks and modifies stratified activity data ready for input to sa function
# INPUT
# obs_data: dataframe of observations with fields deploymentID, time
# dep_data: dataframe of deployments with fields deploymentID, stratumID, effort
# str_data: dataframe of strata with fields stratumID, area
# OUTPUT
# list of dataframes:
#   obs: observations with stratumID field added
#   str: strata with area, effort and observations fields added
make_sa_data <- function(obs_data, dep_data, str_data){
  # Check for required fields and IDs)
  orf <- c("deploymentID", "time")
  drf <- c("deploymentID", "stratumID", "effort")
  srf <- c("stratumID", "area")
  msg <- function(df, ff) paste(df, "must contain columns:", paste(ff, collapse=", "))
  if(!all(orf %in% names(obs_data))) stop(msg("obs_data", orf))
  if(!all(drf %in% names(dep_data))) stop(msg("dep_data", drf))
  if(!all(srf %in% names(str_data))) stop(msg("str_data", srf))
  if(!all(obs_data$deploymentID %in% dep_data$deploymentID)) 
    stop("Cant find all observation deploymentIDs in dep_data")
  if(!all(dep_data$stratumID %in% str_data$stratumID)) 
    stop("Cant find all deployment stratumIDs in str_data")
  
  # add stratumID to obs_dat
  if(!"stratumID" %in% names(obs_data)) 
    obs_data <- dep_data %>%
    select(deploymentID, stratumID) %>%
    right_join(obs_data, by="deploymentID")
  
  # add n observations to str_dat, remove strata with no observations
  str_data <- obs_data %>%
    group_by(stratumID) %>%
    summarise(observations=n()) %>%
    right_join(str_data, by="stratumID") %>%
    mutate(observations = ifelse(is.na(observations), 0, observations))
  # add effort to str_dat
  str_data <- dep_data %>%
    group_by(stratumID) %>%
    summarise(effort=sum(effort)) %>%
    right_join(str_data, by="stratumID")
  
  list(obs=obs_data, str=str_data)
}

# function fits a stratified activity pattern  with or without bootstrapping
# INPUT
# dat: list of dataframes
#     obs with required fields: stratumID, time
#     str with required fields: stratumID, area, effort, observations
# wt: times by strata matrix of weights
sa <- function(dat, wt, adj=1, boot=FALSE){
  nt <- nrow(wt)
  odat <- lapply(dat$str$stratumID, function(s) subset(dat$obs, stratumID==s)$time)
  if(boot) 
    odat <- lapply(odat, function(x) sample(x, length(x), replace=TRUE))
  sq <- seq(0,2*pi,len=nt)
  f <- sapply(odat, function(d) 
    if(length(d)==0) rep(0,nt) else dvmkern(sq, d, adj=adj))
  colnames(f) <- dat$str$stratumID
  fw <- f * wt
  sum_fw <- apply(fw, 1, sum)
  pdf <- sum_fw / (2*pi*mean(head(sum_fw,-1))) #activity pdf over time
  act <- 1 / (2 * pi * max(pdf)) # overall activity level
  popdist <- fw / sum_fw # population distribution between strata
  pa <- popdist * pdf / max(pdf) # population distribution x activity level
  psum <- apply(head(popdist, -1), 2, sum)
  act_stratum <- apply(head(pa, -1), 2, sum) / ifelse(psum==0, 1, psum)
  list(act = act,
       act_stratum = act_stratum,
       pdf = pdf,
       popdist = popdist)
}

# gets summary error list from bootstrap results (boots)
get_errlist <- function(boots){
  lim <- c(0.025, 0.975)
  # overall activity level errors
  boot_act <- unlist(boots["act", ])
  se_act <- sd(boot_act)
  ci_act <- quantile(boot_act, lim)
  
  # stratum-specific activity level errors
  strata <- names(boots[[2]])
  ns <- length(strata)
  nt <- length(boots[[3]])
  reps <- ncol(boots)
  boot_act_stratum <- matrix(unlist(boots["act_stratum", ]), nrow=ns)
  se_act_stratum <- apply(boot_act_stratum, 1, sd)
  ci_act_stratum <- apply(boot_act_stratum, 1, quantile, lim)
  colnames(ci_act_stratum) <- names(se_act_stratum) <- strata
  
  # activity pattern pdf errors
  boot_pdf <- matrix(unlist(boots["pdf", ]), nrow=nt)
  se_pdf <- apply(boot_pdf, 1, sd)
  ci_pdf <- t(apply(boot_pdf, 1, quantile, lim))
  
  # population distribution errors
  boot_popdist <- array(unlist(boots["popdist", ]), dim=c(nt, ns, reps))
  se_popdist <- apply(boot_popdist, 1:2, sd)
  ci <- apply(boot_popdist, 1:2, quantile, lim)
  lcl_popdist <- ci[1,,]
  ucl_popdist <- ci[2,,]
  colnames(lcl_popdist) <- colnames(ucl_popdist) <- colnames(se_popdist) <- strata
  
  list(se = list(act=se_act, act_stratum=se_act_stratum, pdf=se_pdf, 
                 popdist=se_popdist),
       ci = list(act=ci_act, act_stratum=ci_act_stratum, pdf=ci_pdf,
                 popdist=list(lcl=lcl_popdist, ucl=ucl_popdist)))
}

# Fit a stratified activity distribution

# INPUT
# obs_dat: observations data with columns deploymentID and time
# dep_dat: deployment data with columns deploymentID and stratumID and effort
# str_dat: stratum data with columns stratumID and area

# OUTPUT
# A list with elements:
# data: input data, components dep_dat, str_dat, obs_dat
# actmods: stratum-specific activity models
# popdist: time by stratum matrix, proportion of population in each stratum time by time of day
# pdf: probability density function for overall activity pattern
# act: overall activity level
# act_stratum: stratum-specific activity levels
#obs_data <- dat$obs
#str_data <- dat$str
fitact_strat <- function(obs_data, dep_data, str_data, 
                         speed=NULL, radius=NULL, angle=NULL,
                         adj=1, reps=100, nt=513){
  
  check_param <- function(prm){
    strata <- as.character(str_data$stratumID)
    if(is.null(prm)) return(1) else
      if(class(prm)[1]=="matrix" & mode(prm)=="numeric" & !any(is.na(prm))){
        if(nrow(prm)==nt & ncol(prm)==nrow(str_data) & 
           identical(sort(colnames(prm)), sort(strata)))
          return(prm[, strata]) else
            stop("Parameter matrix must have dimensions nt by nrow(str_data), and have column names matching str_data$stratumID")
      } else
        stop("Parameter must be a matrix with no missing values")
  }
  
  dat <- make_sa_data(obs_data, dep_data, str_data)  
  
  speed <- check_param(speed)
  radius <- check_param(radius)
  angle <- check_param(angle)
  q <- speed * radius * (2+angle)
  A <- matrix(rep(dat$str$area, each=nt), nrow=nt)
  E <- rep(dat$str$effort, each=nt)
  Y <- rep(dat$str$observations, each=nt)
  wt <- Y * A / (q * E)
  
  est <- sa(dat, wt, adj)
  if(reps==0) err <- NULL else{
    boots <- pbapply::pbreplicate(reps, sa(dat, wt, adj, TRUE))
    err <- get_errlist(boots)
  }
  res <- list(data=list(obs=dat$obs, dep=dep_data, str=dat$str), 
              est=est, se=err$se, ci=err$ci)
}

time_of_day <- function(x, scale=c("radian", "hour", "proportion")){
  scale <- match.arg(scale)
  m <- switch(scale, 
              radian = 2*pi, 
              hour = 24, 
              proportion = 1)
  
  x <- lubridate::ymd_hms(x) %>% 
    lubridate::local_time() %>% 
    as.numeric()
  m * x  / (24 * 60^2) 
}

# Calculate SE of the mean of several estimates and their SEs
se_from_ses <- function(est, se){
  if(length(est) != length(se))
    stop("All input vectors must have the same length")
  sqrt(mean((mean(est) - est)^2 + se^2))
}

# generate activity pattern from parameter vector prm
pattern <- function(prm, x=seq(0, 2*pi, len=513)){
  if(length(prm)==1){
    1/(1+1/exp(prm*cos(x)))
  } else{
    if("max" %in% names(prm)){
      mx <- prm[names(prm)=="max"]
      prm <- prm[names(prm)!="max"]
    } else
      mx <- 1
    res <- dvmsmix(x, prm)
    mx * res / max(res)
  }
}

# plot patterns, ... are parameter vectors
plot_patterns <- function(prms){
  x <- seq(0, 2*pi, len=513)
  yy <- lapply(prms, pattern, x=x)
  plot(c(0,2*pi), 0:1, type="n", xlab="time", ylab="activity")
  for(i in 1:length(yy)) lines(x, yy[[i]], col=i)
}

# make crepuscular parameter vector
make_crep_prm <- function(kaps, pp, pks=pi*c(0.5, 1.5), max=NULL){
  kk <- rep(kaps, each=2)
  mm <- c(pks, mean(pks), wrap(mean(pks+c(2*pi,0))))
  names(kk) <- paste0("k", 1:4)
  names(mm) <- paste0("m", 1:4)
  names(pp) <- paste0("p", 2:4)
  prm <- c(kk, mm, pp, max=max)
}

# make diurnal parameter vector
make_diur_prm <- function(n, kap, on, off, max=NULL){
  kk <- rep(kap, n)
  mm <- seq(on, off, len=n)
  pp <- rep(1/n, n-1)
  names(kk) <- paste0("k", 1:n)
  names(mm) <- paste0("m", 1:n)
  names(pp) <- paste0("p", 2:n)
  c(kk, mm, pp, max=max)
}

# dep_data: deployment data with fields deploymentID, stratumID, effort
# str_data: stratum data with fields stratumID, area
# N: population size
# ppop: population distribution pattern parameter
# aprmH1, aprmH2: activity pattern parameters for habitats 1 and 2 respectively
# nobs: number of observations to generate
# reps: number of bootstrap replicates for resampling data in activity model fitting
run_scenario <- function(dep_data, str_data, N, 
                         pprm, aprmH1, aprmH2=aprmH1, 
                         nobs, reps=999){
  # generating activity & populations patterns
  a <- cbind(pattern(aprmH1), pattern(aprmH2))
  p <- cbind(1-pattern(pprm), pattern(pprm))
  # calculate true activity values and pdf
  apat <- p[,1]*a[,1] + p[,2]*a[,2]
  true <- list(act = mean(apat),
               act_stratum = c(sum(p[,1]*a[,1]) / sum(p[,1]),
                               sum(p[,2]*a[,2]) / sum(p[,2])),
               pdf = apat*512/(2*pi*sum(head(apat,-1))),
               popdist = p)
  
  # speed/radius/angle generation functions
  sfunc <- function() return(8) 
  rfunc <- function() return(0.005)
  afunc <- function() return(0.8)
  # generate synthetic data and fit stratified activity model
  dat <- data_gen(dep_data, str_data, N, a, p, sfunc, rfunc, afunc)
  dat$obs <- slice_sample(dat$obs, n=nobs)
  n <- dat$obs %>%
    left_join(dat$dep, by="deploymentID") %>%
    group_by(stratumID) %>%
    summarise(n=n())
  if(nrow(n)==1){
    n <- if(n$stratumID==1) c(n$n, 0) else c(0, n$n)
  } else
    n <- n$n
  strat <- fitact_strat(dat$obs, dat$dep, dat$str, adj=1.5, reps=reps)
  samp <- if(reps==0) "none" else "data"
  naive <- fitact(dat$obs$time, sample=samp, adj=1.5, reps=reps)
  list(data=dat, n=n, pprm=pprm, aprm=list(H1=aprmH1, H2=aprmH2),
       true=true, strat=strat$est, naive=list(pdf=naive@pdf[,2], act=naive@act))
}

# run scenario index s from scenarios table
# regenerates data reps times, but does not resample data within datasets
# stratification balances observation sample sizes, not deployments, between strata
run_scenario_s <- function(s, reps){
  print(paste("scenario", s))
  aln <- scenarios$alignment[s]
  pat <- scenarios$pattern[s]
  eff <- scenarios$effort[s]
  dep <- scenarios$depdat[s]
  str <- scenarios$strdat[s]
  aprmH1 <- aprm[[pat]][[1]] # parameters
  aprmH2 <- aprm[[pat]][[aln]]
  pp <- pprm[[pat]]
  # n deployments (to balance n observations between strata)
  A <- area[[str]]
  rel_traprate <- c(sum(pattern(aprmH1) * (1-pattern(pp))),
                    sum(pattern(aprmH2) * pattern(pp))) / A
  ndep <- switch(dep,
                 rep = round(deps * A / sum(A)),
                 strat = round(deps * (1 - rel_traprate / sum(rel_traprate))))
  # input data
  depdat <- data.frame(deploymentID = 1:sum(ndep),
                       stratumID = rep(1:nstr, ndep),
                       effort = 1000)
  strdat <- data.frame(stratumID = 1:nstr,
                       area = A)
  # replicated results
  replicate(reps, run_scenario(dep_data=depdat, str_data=strdat, N=N,
                               pprm=pp, aprmH1=aprmH1, aprmH2=aprmH2,
                               nobs=ssize[[eff]], reps=0), simplify=FALSE)
}

calc_errors <- function(allfits){
  calc <- function(fits){
    data.frame(
      n1 = sapply(fits, function(x) x$n[1]),
      n2 = sapply(fits, function(x) x$n[2]),
      act_strat = sapply(fits, function(x) x$strat$act - x$true$act),
      act_naive = sapply(fits, function(x) x$naive$act - x$true$act),
      act_str1 = sapply(fits, function(x) 
        x$strat$act_stratum[1] - x$true$act_stratum[1]),
      act_str2 = sapply(fits, function(x) 
        x$strat$act_stratum[2] - x$true$act_stratum[2]),
      popDist = sapply(fits, function(x)
        mean(abs(x$strat$popdist[,1] - x$true$popdist[,1]))),
      actPat_strat = sapply(fits, function(x)
        mean(abs(x$strat$pdf - x$true$pdf))),
      actPat_naive = sapply(fits, function(x)
        mean(abs(x$naive$pdf - x$true$pdf)))
    )
  }
  errors <- lapply(allfits, calc)
  i <- rep(1:length(errors), each=nrow(errors[[1]]))
  errors %>%
    bind_rows() %>%
    bind_cols(scenarios[i,])
}

# plot overall activity patterns
plot_activity_patterns <- function(fits, labels=TRUE){
  cf <- pi/12
  true <- fits[[1]]$true$pdf * cf
  strmat <- sapply(fits, function(x) x$strat$pdf) * cf
  naimat <- sapply(fits, function(x) x$naive$pdf) * cf
  dat <- data.frame(Time = seq(0, 24, len=513),
                    Key = rep(c("True", "Stratified", "Unstratified"), each=513),
                    pdf = c(true,
                            apply(strmat, 1, median),
                            apply(naimat, 1, median)),
                    lcl = c(true,
                            apply(strmat, 1, quantile, 0.025),
                            apply(naimat, 1, quantile, 0.025)),
                    ucl = c(true,
                            apply(strmat, 1, quantile, 0.975),
                            apply(naimat, 1, quantile, 0.975)))
  plt <- ggplot(dat, aes(x=Time, y=pdf, col=Key)) +
    geom_line() +
    geom_ribbon(aes(ymin=lcl, ymax=ucl, fill=Key), alpha=0.2, linetype=0) +
    theme_classic()
  if(!labels) plt <- plt +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none")
  plt
}

# plot population distribution pattern
plot_pop_distribution <- function(fits, labels=TRUE){
  ppat <- fits[[1]]$true$popdist[,1]
  pmat <- sapply(fits, function(x) x$strat$popdist[,1])
  dat <- data.frame(Time = seq(0, 24, len=513),
                    Key = rep(c("True", "Estimated"), each=513),
                    proportion = c(ppat, apply(pmat, 1, median)),
                    lcl = c(ppat, apply(pmat, 1, quantile, 0.025)),
                    ucl = c(ppat, apply(pmat, 1, quantile, 0.975)))
  plt <- ggplot(dat, aes(x=Time, y=proportion, col=Key)) +
    geom_line() +
    geom_ribbon(aes(ymin=lcl, ymax=ucl, fill=Key), alpha=0.2, linetype=0) +
    ylim(0,1) +
    theme_classic()
  if(!labels) plt <- plt +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none")
  plt
}

plot_naive_actmod <- function(fit, type=c("stratum", "overall", "both"), ...){
  type <- match.arg(type)
  m <- list()
  key <- character(0)
  if(type!="stratum"){
    key <- c(key, "overall")
    m <- c(m, fitact(fit$data$obs$time))
  }
  if(type!="overall"){
    obs <- fit$data$obs %>%
      left_join(fit$data$dep, by="deploymentID")
    for(s in fit$data$str$stratumID){
      key <- c(key, paste0("stratum", s))
      m <- c(m, fitact(subset(obs, stratumID==s)$time))
    }
  }
  getfrq <- function(m){
    n <- length(m@data)
    dat <- as.data.frame(m@pdf) %>%
      mutate(y = y * n*pi/12) %>%
      rename(time = x, frequency = y)
    if("lcl" %in% names(dat))
      dat <- dat %>%
      mutate(lcl = lcl*n*pi/12,
             ucl = ucl*n*pi/12)
    dat
  }
  dat <- lapply(m, getfrq) %>%
    bind_rows() %>%
    mutate(key = rep(key, each=513))
  p <- ggplot(dat, aes(time, frequency, col=key)) + 
    geom_line() +
    theme_classic()
  if("lcl" %in% names(dat))
    p <- p +
    geom_line(aes(time, lcl), linetype=2) +
    geom_line(aes(time, ucl), linetype=2)
  p
}

plot_scenario_patterns <- function(s){
  aln <- scenarios$alignment[s]
  pat <- scenarios$pattern[s]
  p <- pattern(pprm[[pat]])
  data.frame(Metric = fct_relevel(rep(c("Pop.distribution", "Rel.activity"), each=1026),
                                  c("Rel.activity", "Pop.distribution")),
             Stratum = rep(as.character(1:2), each=513),
             Time = seq(0, 24, len=513),
             Value = c(1-p, p, 
                       pattern(aprm[[pat]][[1]]),
                       pattern(aprm[[pat]][[aln]]))) %>%
    ggplot(aes(Time, Value, col=Stratum, lty=Metric)) + 
    geom_line() +
    scale_x_continuous(breaks = seq(0,24,6)) + 
    scale_y_continuous(breaks = 0:1) +
    theme_classic() +
    theme(axis.title = element_blank())
}

make_act_prms <- function(times, slopes, heights, nper, mx){
  mkfunc <- function(par){
    if(par[2]<par[1]) par[2] <- par[2]+2*pi
    tdif <- par[2] - par[1]
    kdif <- par[4] - par[3]
    bk <- kdif/tdif
    mm <- activity::wrap(seq(par[1], par[2]-tdif/par[5], len=par[5]) + tdif/(2*par[5]))
    names(mm) <- paste0("m", 1:length(mm))
    kk <- par[3] + bk * (mm - par[1])
    names(kk) <- paste0("k", 1:length(kk))
    c(mm, kk)
  }
  mk <- cbind(times, c(times[-1], times[1]),
              slopes, c(slopes[-1], slopes[1]),
              nper) %>%
    apply(1, mkfunc, simplify=FALSE) %>%
    unlist()
  mm <- mk[grepl("m", names(mk))]
  names(mm) <- paste0("m", 1:length(mm))
  kk <- mk[grepl("k", names(mk))]
  names(kk) <- paste0("k", 1:length(kk))
  hh <- rep(heights, nper)
  pp <- (hh/sum(hh))[-1]
  names(pp) <- paste0("p", 2:length(kk))
  c(mm, kk, pp, max=mx)
}
