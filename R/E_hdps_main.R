
#' E-value based collaborative targeted maximum likelihood estimation
#'
#' @param A A binary treatment variable
#' @param Y A outcome variable
#' @param W A matrix of covariates
#' @param family gaussian for a continuous outcome and binomial for a binary outcome
#' @param like_type RSS for a continuous outcome and loglike for a binary outcome
#' @param fluctuation linear for a continuous outcome and logistic for a binary outcome
#'
#' @return A list containing the estimated average treatment effect (ATE),its 95% confidence interval and standard error, and the selected covariates
#' @export
#'
#' @examples


library(broom)
#library(tiyesdyverse)
library(grf)
library(SuperLearner)
library(ctmle)
library(tmle)
library(dplyr)
library(EValue)
library(glmnet)
#source("ctmle_stage.R")


E_hdps<-function(A=A,Y=Y,W=W,family = "gaussian",like_type = "RSS",fluctuation = "linear") {

  data <- data.frame(A = A, Y = Y, W)
  stepwise_summary <- data.frame(variable = character(0), estimate = double(0), sd = double(0))
  E_value <- vector("list", length = p)

  initial_model <- lm(Y ~ ., data = data)

  for (i in 1:p) {
    covariate <- colnames(W)[i]
    formula <- as.formula(paste("Y ~ . -", covariate))
    model <- lm(formula, data = data)

    est <- coef(model)["A"]
    sd <- sd(data$Y)
    se <- sqrt(vcov(model)["A", "A"])
    stepwise_summary <- rbind(stepwise_summary, data.frame(variable = covariate, est = est, se = se))
    e_value <- evalues.OLS(est, sd = sd, se = se, delta = 0.9)
    E_value[[i]] <- e_value[[2, 1]]
  }
  E_value_df <- data.frame(variable = colnames(W), E_Value = unlist(E_value))
  sorted_indices <- order(E_value_df$E_Value, decreasing = TRUE)
  sorted_E_value_df <- E_value_df[sorted_indices, ]

  # CTMLE模型
  ctmleDiscrete <- function(Y=Y, A=A, W = W, Wg=W, Q=NULL,
                            preOrder, order , patience=patience,
                            Qbounds=NULL, cvQinit, Qform=NULL, SL.library=NULL,
                            alpha=.995, family= "gaussian",
                            gbound=.025,like_type = "RSS",
                            fluctuation="linear",
                            verbose=TRUE,  detailed=TRUE,
                            PEN=FALSE, V=5, folds = NULL,
                            stopFactor=10^6
  ) {
    n <- length(Y)
    maptoYstar <- fluctuation=="linear"
    if(is.null(colnames(W))){
      colnames(W) <-  paste("W",1:NCOL(W), sep="")
    }
    invalid.name <- nchar(colnames(W)) == 0
    if(any(invalid.name)){
      colnames(W)[invalid.name] <- paste(".internalW", which(invalid.name), sep="")
    }
    if(is.null(colnames(W))){
      colnames(W) <-  paste("W",1:NCOL(W), sep="")
    }
    invalid.name <- nchar(colnames(W)) == 0
    if(any(invalid.name)){
      colnames(W)[invalid.name] <- paste(".internalW", which(invalid.name), sep="")
    }
    if(!is.null(Q)){
      if(NCOL(Q)==2){
        Q <- cbind(QAW=(1-A)*Q[,1] + A*Q[,2], Q1W=Q[,2], Q0W=Q[,1])
      }}
    SL.library <- c("SL.glm", "SL.glmnet")
    Qinit <- stage1(Y, X = cbind(A, W), Qform = NULL, family = "gaussian", Qbounds = NULL, Qinit = Q,
                    alpha = 0.995, SL.library = SL.library, cvQinit = FALSE, maptoYstar = FALSE, verbose = TRUE)
    ##If preOrder is TRUE, then we run cv first to determine the optimal number of covariates to run CTMLE on. Saves Time and allows for stopping early by setting patience parameter
    if(preOrder==TRUE){
      if(is.null(order)){
        warning("The order is not specified. Naive order starting from W1 to Wp would be used")
      }else{
        Wg <- Wg[, order]
      }
      if(verbose) {cat("\t\t-----Cross-validation to select the number of variables to run CTMLE on-----\n")}
      cv.res <- cv(Y=Qinit$Y,X=data.frame(A,Wg),
                   Q=Qinit$Q, family=Qinit$family,
                   Qbounds=Qinit$Qbounds, ab=Qinit$ab,
                   like_type=like_type,gbound=gbound,
                   verbose=verbose, PEN=PEN ,
                   V=V, patience=patience, preOrder=preOrder,
                   folds = folds)
      best_k<- cv.res$best_k
      if(cv.res$best_k==1) best_k<- cv.res$best_k + 1 ##we add one because if best_k=1 then stage2 function gets error (requires at least one covariate). So we will construct one more estimator than needed but we'll end up taking the best one based on cv
      #note: we add one to best_k for X input to make sure we take at least one variable. This will calculate one more estimate than needed
      candidates <- stage2(Y=Qinit$Y, X=data.frame(A,Wg)[,1:best_k], Qinit$Q, family=Qinit$family, Qinit$Qbounds,ab=Qinit$ab,
                           ncandidates=best_k-1, training_set=rep(TRUE,length(Y)),
                           like_type="RSS",gbound=0.025, verbose=TRUE,
                           stopFactor=10^6, preOrder=preOrder)
    }
    if(preOrder==FALSE){
      candidates <- stage2(Y=Qinit$Y, X=data.frame(A,W), Qinit$Q,
                           family=Qinit$family, Qinit$Qbounds,ab=Qinit$ab,
                           ncandidates=NCOL(W), training_set=rep(TRUE,length(Y)),
                           like_type="RSS",gbound=0.025, verbose=TRUE,
                           stopFactor=10^6, preOrder=FALSE)
      if(verbose) {cat("\t\t-----Cross-validation to select the number of variables to run CTMLE on-----\n")}
      cv.res <- cv(Y=Qinit$Y,X=data.frame(A,W), Q=Qinit$Q,
                   family=Qinit$family, Qbounds=Qinit$Qbounds,
                   ab=Qinit$ab,  like_type="RSS",gbound=0.025,
                   verbose=TRUE, PEN=FALSE , V=5,
                   patience=10, preOrder=FALSE)
    }
    if(verbose){
      cat("\n\n")
      cat("\t terms included in final estimate: ", candidates$terms, "\n")
      cat("\t all estimates: ", candidates$results.all$est, "\n")
      cat("\t best_k: ", cv.res$best_k, "\n")
      cat("\t best estimate: ", candidates$results.all$est[cv.res$best_k], "\n\n")
    }
    # add varDstar and varIC on the original scale
    candidates$results.all$varIC_origscale <- candidates$results.all$varIC * ifelse(maptoYstar, (diff(Qinit$ab))^2, 1)
    candidates$results.all$varDstar_origscale <- candidates$results.all$varDstar * ifelse(maptoYstar, (diff(Qinit$ab))^2, 1)
    returnValc <- (list(candidates=candidates, best_k=cv.res$best_k,
                        est= candidates$results.all$est[cv.res$best_k],
                        CI = candidates$results.all$est[cv.res$best_k]+c(-1.96,1.96) * sqrt(candidates$results.all$varIC_origscale[cv.res$best_k]/n),
                        Qinit=Qinit,
                        ab = Qinit$ab,
                        likelihood= candidates$results.all$likelihood[cv.res$best_k],
                        varIC = candidates$results.all$varIC_origscale[cv.res$best_k],
                        varDstar = candidates$results.all$varDstar_origscale[cv.res$best_k],
                        varIC.cv=cv.res$varIC[cv.res$best_k],
                        var.psi=candidates$results.all$varIC_origscale[cv.res$best_k]/n,
                        penlikelihood.cv = cv.res$penlikelihood[cv.res$best_k],
                        cv.res=cv.res))
    returnValc$pvalue <- 2*pnorm(-abs(returnValc$est)/sqrt(returnValc$var.psi))
    just_enough <- c(2:4,6,8,9,11,14)

    returnValc$pvalue <- 2*pnorm(-abs(returnValc$est)/sqrt(returnValc$var.psi))
    if (!detailed) {
      returnValc <- cbind(returnValc[just_enough], pvalue=returnValc$pvalue)
    } else {
      returnValc$pvalue <- NULL
    }
    class(returnValc) <- "ctmle"
    return(returnValc)
  }
  # 构建CTMLE模型
  result<-ctmleDiscrete(Y, A, W = W, Wg = W, Q = NULL,
                        preOrder = TRUE, order = sorted_E_value_df$variable, patience = patience,
                        Qbounds = NULL, cvQinit = FALSE, Qform = NULL, SL.library = NULL,
                        alpha = 0.995, family = "gaussian",
                        gbound = 0.025, like_type = "RSS",
                        fluctuation = "linear",
                        verbose = FALSE, detailed = TRUE,
                        PEN = FALSE, V = 5, folds = NULL,
                        stopFactor = 10^6)
  selected_covariates<-result[["candidates"]][["terms"]]
  selected_covariates<-selected_covariates[-1]
  min_required_vars <- 2
  while (length(selected_covariates) < min_required_vars) {
    current_vars <- selected_covariates
    remaining_vars <- setdiff(sorted_E_value_df$variable, current_vars)

    result<-ctmleDiscrete(Y, A, W = W, Wg = W, Q = NULL,
                          preOrder = TRUE, order = remaining_vars, patience = patience,
                          Qbounds = NULL, cvQinit = FALSE, Qform = NULL, SL.library = NULL,
                          alpha = 0.995, family = "gaussian",
                          gbound = 0.025, like_type = "RSS",
                          fluctuation = "linear",
                          verbose = FALSE, detailed = TRUE,
                          PEN = FALSE, V = 5, folds = NULL,
                          stopFactor = 10^6)
    selected_covariates<-result[["candidates"]][["terms"]]
    selected_covariates<-selected_covariates[-1]
    selected_covariates <- union(current_vars,  selected_covariates)
  }

  returnVal <- list(est=result$est,
                    var.psi = result$var.psi,
                    sd = sqrt(result$var.psi),
                    se=sqrt(result$var.psi)/sqrt(N),
                    CI=result$CI,
                    selected_covariates=selected_covariates,
                    pvalue=result$pvalue)
  return(returnVal)


}
