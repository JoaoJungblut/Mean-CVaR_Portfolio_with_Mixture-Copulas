library(tidyverse)
library(tidyquant)
library(rugarch)
library(fGarch)
library(copula)
library(Rsolnp)

# TRANSFORMAR UNIF_DIST EM MATRIZ PARA RODAR NO FORLOOP INICIAL ASSIM COMO COP_PARAM

# importing data and calculating returns
returns <- tq_get(c("PETR4.SA", "VALE3.SA", "ITUB4.SA",
                    "BBAS3.SA", "ABEV3.SA"), 
                  from="1997-01-01") %>% 
  select(date, symbol, adjusted) %>% 
  group_by(symbol) %>% 
  mutate(return = log(adjusted) - log(dplyr::lag(adjusted))) %>% 
  ungroup() %>% 
  select(-adjusted) %>% 
  spread(symbol, return) %>% 
  na.omit()

# creating auxiliary matrices and list
N <- ncol(returns)-1
K <- 252
index_vector <- seq(1, nrow(returns), by=K)
names_vector <- names(returns)[-1]
cop_param <- unif_dist <- garch_pred <- garch_coef <- sigma <- residuals <- vector("list", length(index_vector)-1)
mod_garch <- try(ugarchspec(variance.model = list(model = "sGARCH",
                                                  garch.order = c(1, 1),
                                                  variance.targeting = TRUE),
                            mean.model = list(armaOrder = c(1, 0)),
                            distribution.model = "sstd"),
                 silent = TRUE)

#initializing copulas object
copt <- tCopula(param = 0.5, dim = N) 
copC <- claytonCopula(2, dim = N) # delta= 2
copG <- gumbelCopula(2, dim = N)  # theta= 2

#lower and upper bounds of the parameters and weights for bounded non linear opt.
lower <- c(0.1, 1, -0.9,(2+.Machine$double.eps), 0,0,0)     
upper <- c(copC@param.upbnd, copG@param.upbnd, 1,100, 1,1,1) #2+eps so that variance of t copula is defined

#initial guesses for weights = 1/3 each
par6 <- par5 <- par4 <- 1/3 

#negative log likelihood function for estimating copula weights and parameters
LLCG <- function(params, U, copC, copG, copt){ 
  slot(copC, "parameters") <- params[1]    #Initial Clayton parameter provided
  slot(copG, "parameters") <- params[2]    #Initial Gumbel parameter provided 
  slot(copt, "parameters") <- params[3:4]  #Initial t parameters provided (correlation and Degrees of Freedom)
  pi1 <- params[5] #weight of Clayton copula
  pi2 <- params[6] #weight of Gumbel copula
  pi3 <- params[7] #weight of t copula
  
  opt <- log(pi1 * dCopula(U, copC) + pi2 * dCopula(U, copG)
             + pi3 * dCopula(U, copt))  ##loglikelihood function to be optimized
  if(any(is.infinite(opt))){            ##esse IF está no Pfaff
    opt[which(is.infinite(opt))]<-0
  }
  -sum(opt)
}

#constrain function so that sum(weights)=1
eqfun <- function(params, U, copC, copG, copt){ 
  z <- params[5]+params[6]+params[7]
  return(z)
}


# iterating over n assets
for (i in 2:length(index_vector)) {
  
  # Matrix to save sigma forecasts
  unif_dist <- garch_pred <- sigma <- residuals <- matrix(nrow = K, ncol = N)
  garch_coef <- vector("list", length = N)
  
  for (j in 1:length(names_vector)) {
    t1 <- index_vector[i-1]
    t2 <- index_vector[i]-1
    
    x <- cbind(returns[t1:t2,(j+1)])
    
    garch_fit <- try(ugarchfit(mod_garch, data = x, solver = "hybrid"), 
                     silent = TRUE)
    
    residuals[,j] <- garch_fit@fit$residuals
    sigma[,j] <- garch_fit@fit$sigma
    garch_coef[[j]] <- garch_fit@fit$coef
    unif_dist[,j] <- psstd(q = residuals[,j]/sigma[,j],
                           nu = garch_coef[[j]][6],
                           xi = garch_coef[[j]][5])
  }
  
  ##Creating elliptical copula objects and estimating "initial guesses" for each copula parameter. 
  #Then, we maximize loglikelihood of the linear combination of the three copulas
  par1 <- fitCopula(copC, unif_dist, "itau", estimate.variance = T)@estimate #inversion of Kendall's tau for Clayton 
  par2 <- fitCopula(copG, unif_dist,"itau", estimate.variance = T)@estimate #inversion of Kendall's tau for Gumbel 
  par3 <- fitCopula(copt, unif_dist,"mpl", estimate.variance = F)@estimate ###mpl para poder estimar tambem DF. Na documentacao diz que nao pode usar 'itau' pois ele n estima DF.
  
  ##non linear constrained optimization (RSOLNP)
  opt <- solnp(pars = c(par1,par2,par3,par4,par5,par6), 
               fun = LLCG, LB = lower, UB = upper, 
               copt=copt,copC = copC, copG = copG, 
               U=unif_dist,eqfun = eqfun, eqB=c(1)) ####RSOLNP
  
  ##saving optimization parameters in a list
  cop_param <-opt$pars
  
  #clayton, t, gumbel and ctg variates matrix
  ctg <- Cc <- Cg <- Ct <- matrix(nrow = K, ncol = N)
  
  ##generating copula variates
  Cc[,]<- cop_param[5]*rCopula(n = K, 
                               copula = claytonCopula(param = cop_param[1], 
                                                      dim = N))   
  Cg[,]<- cop_param[6]*rCopula(n = K, 
                               copula = gumbelCopula(param = cop_param[2], 
                                                     dim = N))
  Ct[,]<- cop_param[7]*rCopula(n = K, 
                               copula = tCopula(param = cop_param[3], 
                                                df = cop_param[4], 
                                                dim = N))
  
  #linear combination of them
  ctg <- Cc + Ct + Cg 
  
  #for each asset, generate copula 'z' dependence structure 
  rtn_pred <- mean_pred <- sigma_pred <- zsim <- matrix(nrow = K, ncol = N) 
  for(j in 1:N){  
    zsim[,j] <- qsstd(ctg[,j], 
                      nu = garch_coef[[j]][6], 
                      xi = garch_coef[[j]][5]) / 
      sd(qsstd(ctg[,j], nu = garch_coef[[j]][6], 
               xi = garch_coef[[j]][5])) 
    
    rtn_t <- as.numeric(returns[(t2-1),(j+1)])
    sigma_t <- sigma[K,j]
    
    sigma_pred[1,j] <- sqrt(garch_coef[[j]][7] + #omega
                              garch_coef[[j]][3]*(rtn_t)^2 + #alpha1 
                              garch_coef[[j]][4]*(sigma_t)^2) #beta1
    mean_pred[1, j] <- garch_coef[[j]][1] + garch_coef[[j]][2]*rtn_t #mu e ar1
    rtn_pred[1, j] <- mean_pred[1, j] + sigma_pred[1, j]*zsim[1,j]
    
    for(t in 2:K){
      sigma_pred[t, j] <- sqrt(garch_coef[[j]][7] + #omega
                                 garch_coef[[j]][2]*(rtn_pred[(t-1), j])^2 + #alpha1
                                 garch_coef[[j]][3]*(sigma_pred[(t-1), j])^2) #beta1
      mean_pred[t, j] <- garch_coef[[j]][1] + garch_coef[[j]][2]*rtn_pred[(t-1), j] #mu e ar1
      rtn_pred[t, j] <- mean_pred[t, j] + sigma_pred[t, j]*zsim[t,j]
    }
    
  }
 
  
}



cop_portf_opt <- function(targetReturn, filename, nsim, type){
  cop_pars <- readRDS("copulaParams.Rds") #reading copula parameters
  garch_coefs <- readRDS("coef.Rds") #reading ArmaGarch parameters
  #arma_order <- readRDS("armaOrder.Rds")
  sigma_fit <- readRDS("sigma.Rds") #reading armaGarch sigma. We need this to estimate J one-steap ahead returns 
  residual_fit <- readRDS("residuos.Rds")  #reading armaGach fitted residuals 
  
  #####setting up a fGarch min CVaR optimization
  frontierSpec <- portfolioSpec() 
  setType(frontierSpec) <- "CVaR"  
  setSolver(frontierSpec) <- "solveRglpk.CVAR"  #solving as a linear program as Rockafellar & Uryasev (2000)
  setAlpha(frontierSpec) <- 0.05   #CVaR_0.95  
  setTargetReturn(frontierSpec) <- targetReturn  #daily target return constrain, we do this for 0, 0.00012, 0.00024 and 0.00036
  
  ####initializing returns and weights matrixes
  ret_sim <- matrix(0, nrow = nsim, ncol = 8) #initializaing sim.ret matrix
  cvar_opt <- matrix(0, nrow = 4550, ncol = 8) #initializing optimized portfolio weights matrix
  
  for(i in 1:length(sigma_fit)){  #we do everything for each optimization
    ##generating copula variates 
    if(type == 'mixture'){         #if we want mixture of copula, type should be 'mixture'
      cop_sim <- mix_cop_sim(cop_pars, i, nsim)
    } else if(type == 'gaussian'){ #else, 'gaussian'
      cop_sim <- gauss_cop_sim(cop_pars, i, nsim)
    }
    ##generating zsim 
    zsim <- z_cop_sim(garch_coefs, cop_sim, i, nsim)
    #sigma_per <- sigma_fit[[i]]  #getting fitted sigma and residuals we'll use in MA and Garch's forecasting
    #resid_per <- residual_fit[[i]] 
    
    ##K return scenarios generation, using z and armaGarch
    ret_sim <- ret_cop_sim(garch_coefs, zsim, sigma_fit[[i]], residual_fit[[i]], returns, i, nsim)
    
    ##optimizing portfolio using K simulated return for each assets, for optimization period i 
    retornofPort <- as.timeSeries(ret_sim[, 1:8])
    frontier1g <- fPortfolio::efficientPortfolio(data = retornofPort, spec = frontierSpec, constraints = "LongOnly")
    cvar_opt[i,1:8] <- fPortfolio::getWeights(frontier1g)   #storing resulting weights   
  }
  saveRDS(cvar_opt, file = filename) ##saving weights data, we repeat this for 0,3,6 and 9%
}
}