# Three-arm confidential global permutation simulation
# Base R; reproducible, manuscript Section 7 implementation
set.seed(20261107)
OUT <- 'GRD_three_arm_results'
dir.create(OUT, showWarnings=FALSE)
B <- 5000L   # change to 100 for a quick smoke test
N_VALUES <- c(100L,300L,500L,1000L,3000L)
DELTAS <- c(0,0.2,0.5,0.8)
ALLOCATIONS <- list(equal=c(1/3,1/3,1/3), unequal=c(0.5,0.3,0.2))
LOOK_FRACTIONS <- c(.25,.50,.75,1)
THRESHOLD <- .80
SIGMA <- 1
# Rows are candidate mappings: true arm g -> anonymous label pi[g]
PERMS <- as.matrix(expand.grid(1:3,1:3,1:3))
PERMS <- PERMS[apply(PERMS,1,function(z) length(unique(z))==3),,drop=FALSE]
stopifnot(nrow(PERMS)==6L)
# E0: allocation only; E1: means; E2: means and variances;
# E3: means, variances, plus a simulated treatment-specific biomarker.
# Prior uncertainty is a modelling choice, not specified by the manuscript.
REGIMES <- c('E0','E1','E2','E3')
MODELS <- c('correct','misspecified')
# A historical means prior N(m_g, tau_mu^2). E2 assumes more precise
# historical mean information plus known common outcome variance.
TAU_MU <- .35
TAU_BIOM <- .30
BIOM_SD <- 1
BIOM_MEANS <- c(-.8,0,.8)
# Incorrect external knowledge: reverse named-treatment profiles.
external_profile <- function(mu, regime, model) {
  if(model=='misspecified') {
    mu <- rev(mu)
    bio <- rev(BIOM_MEANS)
    prob <- rev(c(.5,.3,.2))
  } else {
    bio <- BIOM_MEANS
    prob <- c(.5,.3,.2)
  }
  list(mu=mu, bio=bio, prob=prob)
}
# Posterior over all 6 mappings. For E0, known unequal allocation
# contributes multinomial likelihood; for equal allocation it is flat.
# For E1, historical means have predictive uncertainty tau_mu.
# E2 uses historical named-arm variances (here equal by design) and
# a tighter mean prior, so E2 is distinct from E1.
# E3 adds a biomarker with distinct named-treatment distributions.
posterior_mapping <- function(R,Y,Bio,alloc,mu,regime,model) {
  n_h <- tabulate(R,3)
  ybar <- vapply(1:3,function(h) mean(Y[R==h]),numeric(1))
  bbar <- vapply(1:3,function(h) mean(Bio[R==h]),numeric(1))
  prof <- external_profile(mu,regime,model)
  # Under equal allocation E0 remains exactly uniform. Under unequal
  # allocation E0 uses the publicly known allocation probabilities.
  known_p <- if(all(abs(alloc-1/3)<1e-12)) rep(1/3,3) else
    if(model=='correct') alloc else rev(alloc)
  logw <- numeric(6)
  for(k in 1:6) {
    pi <- PERMS[k,]
    logw[k] <- sum(n_h[pi]*log(known_p))
    if(regime!='E0') {
      tau <- if(regime=='E1') TAU_MU else TAU_MU/2
      # predictive variance of group mean, integrating prior uncertainty
      sd_pred <- sqrt(tau^2 + SIGMA^2/n_h[pi])
      logw[k] <- logw[k]+sum(dnorm(ybar[pi],prof$mu,sd_pred,log=TRUE))
    }
    if(regime=='E3') {
      sd_b <- sqrt(TAU_BIOM^2+BIOM_SD^2/n_h[pi])
      logw[k] <- logw[k]+sum(dnorm(bbar[pi],prof$bio,sd_b,log=TRUE))
    }
  }
  z <- exp(logw-max(logw)); z/sum(z)
}
metrics <- function(post, actual_pi) {
  H <- -sum(ifelse(post>0,post*log(pmax(post,1e-300)),0))
  C <- matrix(0,3,3,dimnames=list(paste0('anonymous_',1:3),paste0('treatment_',1:3)))
  for(k in 1:6) for(g in 1:3) C[PERMS[k,g],g] <- C[PERMS[k,g],g]+post[k]
  truth <- which(apply(PERMS,1,function(z) all(z==actual_pi)))
  # Ties resolved by random MAP selection (avoids artificial E0 advantage)
  winners <- which(abs(post-max(post))<1e-12)
  map_guess <- sample(winners,1)
  c(d_pi=max(post),entropy=H/log(6),n_eff=exp(H),leakage=1-H/log(6),
    d_max=max(C),map_correct=as.numeric(map_guess==truth),
    true_mapping_posterior=post[truth],
    d_group1=max(C[1,]),d_group2=max(C[2,]),d_group3=max(C[3,]))
}
# One-way ANOVA F and p without fitting a model repeatedly
anova_stats <- function(Y,group) {
  ns <- tabulate(group,3)
  means <- vapply(1:3,function(g) mean(Y[group==g]),numeric(1))
  grand <- mean(Y)
  ssb <- sum(ns*(means-grand)^2)
  ssw <- sum((Y-means[group])^2)
  f <- (ssb/2)/(ssw/(length(Y)-3))
  c(F=f,p=pf(f,2,length(Y)-3,lower.tail=FALSE))
}
# A single global permutation is reused at every cumulative look.
simulate_trial <- function(n,delta,alloc,regime,model,rep_id) {
  mu <- c(0,delta,2*delta)
  A <- sample(1:3,n,replace=TRUE,prob=alloc)
  # Ensure groups are nonempty at the earliest look; negligible at n>=100.
  while(length(unique(A[seq_len(max(3,floor(n*.25)))]))<3)
    A <- sample(1:3,n,replace=TRUE,prob=alloc)
  Y <- rnorm(n,mu[A],SIGMA)
  Bio <- rnorm(n,BIOM_MEANS[A],BIOM_SD)
  pi <- sample(1:3,3)
  R <- pi[A]
  looks <- unique(pmax(3L,as.integer(floor(n*LOOK_FRACTIONS))))
  output <- vector('list',length(looks))
  for(j in seq_along(looks)) {
    idx <- seq_len(looks[j]); a <- A[idx]; r <- R[idx]
    stA <- anova_stats(Y[idx],a)
    stR <- anova_stats(Y[idx],r)
    mA <- vapply(1:3,function(g) mean(Y[idx][a==g]),numeric(1))
    mR <- vapply(1:3,function(h) mean(Y[idx][r==h]),numeric(1))
    post <- posterior_mapping(r,Y[idx],Bio[idx],alloc,mu,regime,model)
    mt <- metrics(post,pi)
    output[[j]] <- data.frame(rep=rep_id,look=j,look_n=looks[j],
      t(mt),mean_diff=max(abs(mA-mR[pi])),
      F_diff=abs(stA['F']-stR['F']),p_diff=abs(stA['p']-stR['p']),
      stringsAsFactors=FALSE,check.names=FALSE)
  }
  do.call(rbind,output)
}
# Scenario-level summaries with Monte Carlo SEs.
summarize_one <- function(z) {
  measures <- c('d_pi','entropy','n_eff','leakage','d_max',
                'map_correct','true_mapping_posterior','d_group1',
                'd_group2','d_group3')
  result <- list()
  for(m in measures) {
    v <- z[[m]]
    result[[paste0(m,'_mean')]] <- mean(v)
    result[[paste0(m,'_mcse')]] <- sd(v)/sqrt(length(v))
  }
  result$max_mean_diff <- max(z$mean_diff)
  result$max_F_diff <- max(z$F_diff)
  result$max_p_diff <- max(z$p_diff)
  as.data.frame(result,check.names=FALSE)
}
SCENARIOS <- expand.grid(n=N_VALUES,delta=DELTAS,
  allocation=names(ALLOCATIONS),regime=REGIMES,model=MODELS,
  stringsAsFactors=FALSE)
# 5 * 4 * 2 * 4 * 2 = 320 scenarios, each B replications.
# For initial testing use B=100 and/or subset SCENARIOS.
summary_rows <- list(); monitoring_rows <- list()
for(s in seq_len(nrow(SCENARIOS))) {
  sc <- SCENARIOS[s,]
  cat(sprintf('Scenario %d/%d: n=%d delta=%.2f %s %s %s\n',
    s,nrow(SCENARIOS),sc$n,sc$delta,sc$allocation,sc$regime,sc$model))
  trials <- lapply(seq_len(B),function(b) simulate_trial(
    sc$n,sc$delta,ALLOCATIONS[[sc$allocation]],sc$regime,sc$model,b))
  z <- do.call(rbind,trials)
  bylook <- split(z,z$look)
  for(j in seq_along(bylook)) {
    zz <- bylook[[j]]
    summary_rows[[length(summary_rows)+1L]] <- cbind(sc,
      look=as.integer(names(bylook)[j]),look_n=unique(zz$look_n),
      summarize_one(zz))
  }
  bytrial <- split(z,z$rep)
  monitoring <- t(vapply(bytrial,function(zz) {
    crossing <- which(zz$d_pi>THRESHOLD)
    c(min_entropy=min(zz$entropy),max_d_pi=max(zz$d_pi),
      ever_cross=as.numeric(length(crossing)>0),
      first_cross_look=if(length(crossing)) zz$look[min(crossing)] else NA_real_)
  },numeric(4)))
  monitoring_rows[[s]] <- cbind(sc,
    min_entropy_mean=mean(monitoring[,'min_entropy']),
    max_d_pi_mean=mean(monitoring[,'max_d_pi']),
    threshold_cross_rate=mean(monitoring[,'ever_cross']),
    threshold_cross_mcse=sd(monitoring[,'ever_cross'])/sqrt(B),
    first_cross_mean=if(all(is.na(monitoring[,'first_cross_look']))) NA_real_ else
      mean(monitoring[,'first_cross_look'],na.rm=TRUE))
  if(s==1L) write.csv(z,file.path(OUT,'example_trial_level.csv'),row.names=FALSE)
  if(s%%10L==0L) {
    write.csv(do.call(rbind,summary_rows),file.path(OUT,'look_summary.csv'),row.names=FALSE)
    write.csv(do.call(rbind,monitoring_rows),file.path(OUT,'monitoring_summary.csv'),row.names=FALSE)
  }
}
summary_df <- do.call(rbind,summary_rows)
monitoring_df <- do.call(rbind,monitoring_rows)
write.csv(summary_df,file.path(OUT,'look_summary.csv'),row.names=FALSE)
write.csv(monitoring_df,file.path(OUT,'monitoring_summary.csv'),row.names=FALSE)
write.csv(SCENARIOS,file.path(OUT,'scenario_manifest.csv'),row.names=FALSE)
# Reference checks: E0 equal allocation gives uniform posterior for all data.
stopifnot(max(abs(posterior_mapping(rep(1:3,each=5),rnorm(15),rnorm(15),
  rep(1/3,3),c(0,.5,1),'E0','correct')-1/6))<1e-12)
# Invariance must hold to floating point tolerance.
stopifnot(max(summary_df$max_mean_diff)<1e-10,
          max(summary_df$max_F_diff)<1e-8,
          max(summary_df$max_p_diff)<1e-10)
cat('Completed. Outputs saved in ',OUT,'\n',sep='')
