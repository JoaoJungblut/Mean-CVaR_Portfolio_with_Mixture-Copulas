library(tidyverse)
library(tidyquant)
library(rugarch)
library(fGarch)
library(copula)
library(Rsolnp)

# importing data and calculating returns
returns <- tq_get(c("PETR4.SA", "VALE3.SA", "ITUB4.SA",
                    "BBAS3.SA", "ABEV3.SA"), 
                  from="2020-01-01") %>% 
  select(date, symbol, adjusted) %>% 
  group_by(symbol) %>% 
  mutate(return = log(adjusted) - log(lag(adjusted))) %>% 
  ungroup() %>% 
  select(-adjusted) %>% 
  spread(symbol, return) %>% 
  na.omit()

# creating auxiliary matrices and list
n_assets <- ncol(returns)-1
data <- cbind(returns)
residuals <- matrix(nrow = nrow(data) - 250, ncol = ncol(data) - 1)
sigma <- matrix(nrow = nrow(data) - 250, ncol = ncol(data) - 1)
coef <- vector("list", 8)
unif_dist <- matrix(nrow = nrow(data) - 250, ncol = ncol(data) - 1)

# defining interval of estimation
t1 <- 1
t2 <- nrow(data)-250

# iterating over assets
for(i in 2:ncol(data)){
  print(i)
  x <- data[t1:t2, i]
  
  # estimating garch  
  mod_garch <- try(ugarchspec(variance.model = list(model="gjrGARCH",
                                                    garch.order=c(1,1),
                                                    variance.targeting=TRUE),
                              mean.model = list(arma.order=c(1,0)),
                              distribution.model = "sstd"),
                   silent = TRUE) 
  
  garch_fit <- try(ugarchfit(mod_garch, data=x, solver="hybrid"), silent=TRUE)
  
  # saving results  
  residuals[,i-1] <- garch_fit@fit$residuals
  sigma[,i-1] <- garch_fit@fit$sigma
  coef[[i-1]] <- garch_fit@fit$coef
  
  # generating CDF 
  unif_dist[,i-1] <- psstd(q = residuals[, i-1]/sigma[, i-1], 
                     nu = coef[[i-1]][8], 
                     xi = coef[[i-1]][7])
  #unif[,i-1] = pobs(residuals[,i-1]/sigma[,i-1])
}  

residuals <- vector("list", nrow(returns)-250)
sigma <- vector("list", nrow(returns)-250)
coef <- vector("list", nrow(returns)-250)
unif_dist <- vector("list", nrow(returns)-250)

for (i in 1:(nrow(returns) - 250)) {
  for (j in 2:ncol(returns)) {
    x <- cbind(returns[i:(i + 249), j])
    mod_garch <- try(ugarchspec(variance.model = list(model = "gjrGARCH",
                                                      garch.order = c(1, 1),
                                                      variance.targeting = TRUE),
                                mean.model = list(armaOrder = c(1, 0)),
                                distribution.model = "sstd"),
                     silent = TRUE)
    
    garch_fit <- try(ugarchfit(mod_garch, data = x, solver = "hybrid"), 
                     silent = TRUE)
    
    residuals[[i]][[j-1]] <- garch_fit@fit$residuals
    sigma[[i]][[j-1]] <- garch_fit@fit$sigma
    coef[[i]][[j-1]] <- garch_fit@fit$coef
    
    unif_dist[[i]][[j - 1]] <- psstd(q = residuals[[i]][[j-1]]/sigma[[i]][[j-1]],
                                     nu = coef[[i]][[j-1]][7],
                                     xi = coef[[i]][[j-1]][6])
    #unif[[i]][[j - 1]] = pobs(residuals[[i]][[j-1]]/sigma[[i]][[j-1]])
  }
}


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


#creating a list to store each period's parameters
cop_param <- vector('list', length(unif_dist)) 
  
#initializing copulas object
copt <- tCopula(param = 0.5, dim = n_assets) 
copC <- claytonCopula(2, dim = n_assets) # delta= 2
copG <- gumbelCopula(2, dim = n_assets)  # theta= 2
  
#lower and upper bounds of the parameters and weights for bounded non linear opt.
lower <- c(0.1, 1, -0.9,(2+.Machine$double.eps), 0,0,0)     
upper <- c(copC@param.upbnd, copG@param.upbnd, 1,100, 1,1,1) #2+eps so that variance of t copula is defined
  
  
for(i in 1:length(unif_dist)){ 
  ##pseudo-uniform [0,1] observations for each asset
  v<-as.matrix(do.call(cbind, unif_dist[[i]]))
  U<-v[,1:(n_assets)]  
    
  ##Creating elliptical copula objects and estimating "initial guesses" for each copula parameter. 
  #Then, we maximize loglikelihood of the linear combination of the three copulas
  par1 <- fitCopula(copC, U, "itau", estimate.variance = T)@estimate #inversion of Kendall's tau for Clayton 
  par2 <- fitCopula(copG, U,"itau", estimate.variance = T)@estimate #inversion of Kendall's tau for Gumbel 
  par3 <- fitCopula(copt, U,"mpl", estimate.variance = F)@estimate ###mpl para poder estimar tambem DF. Na documentacao diz que nao pode usar 'itau' pois ele n estima DF.
  par4 <- 1/3 #initial guesses for weights = 1/3 each
  par5 <- 1/3
  par6 <- 1/3
    
  ##non linear constrained optimization 
  opt <- solnp(pars = c(par1,par2,par3,par4,par5,par6), 
               fun = LLCG, LB = lower, UB = upper, 
               copt=copt,copC = copC, copG = copG, 
               U=U,eqfun = eqfun, eqB=c(1)) ####RSOLNP
  
  ##saving optimization parameters in a list
  cop_param[[i]]<-opt$pars  
}

# function to generate mixture copulas for potfolio optimization
mix_cop_sim <- function(cop_pars, i, nsim){
  #clayton, t, gumbel variates matrix
  Cc <- Cg <- Ct <- matrix(0,nrow = nsim, ncol = n_assets)   
  ctg <- matrix(0, nrow = nsim, ncol = n_assets)              
  ##generating copula variates
  Cc[,]<- cop_pars[[i]][[5]]*rCopula(n = nsim, 
                                     copula = claytonCopula(param = cop_pars[[i]][[1]], 
                                                            dim = n_assets))   
  Cg[,]<- cop_pars[[i]][[6]]*rCopula(n = nsim, 
                                     copula = gumbelCopula(param = cop_pars[[i]][[2]], 
                                                           dim = n_assets))
  Ct[,]<- cop_pars[[i]][[7]]*rCopula(n = nsim, 
                                     copula = tCopula(param = cop_pars[[i]][[3]], 
                                                      df = cop_pars[[i]][[4]], 
                                                      dim = n_assets))
  ctg <- Cc + Ct + Cg #linear combination of them 
  return(ctg)
}

gauss_cop_sim <- function(cop_pars, i, nsim){
  Gcop <- matrix(0, nrow = nsim, ncol = n_assets)
  Gcop[,] <- copula::rCopula(n = nsim, 
                             copula = normalCopula(param = cop_pars[[i]], 
                                                   dim = n_assets))
  return(Gcop)
}

# function for 'z' copula matrix, resulting from quantile of mixture
z_cop_sim <- function(garch_coefs, ctg, i, nsim){
  zsim <- matrix(0, nrow = nsim, ncol = 8) 
  #for each asset, generate copula 'z' dependence strucutre 
  for(j in 1:n_assets){  
    zsim[,j] <- qsstd(ctg[,j], 
                      nu = garch_coefs[[i]][[j]][[8]], 
                      xi = garch_coefs[[i]][[j]][[7]]) / # garch_coefs[[i]][[j]][[7]]
      sd(qsstd(ctg[,j], nu = garch_coefs[[i]][[j]][[8]], 
               xi = garch_coefs[[i]][[j]][[7]])) #nu = t's DF, xi = t's skew
  }
  return(zsim)
}

#simulated returns matrix
ret_cop_sim <- function(garch_coefs, zsim, sigma_per, resid_per, returns, i, nsim){
  ret_sim <- matrix(0, nrow = nsim, ncol = 8)    
  #getting real returns we'll use in one-step forward armaGarch's AR forecasting
  RZ <- returns[i:(1259+i),2:9] 
  for(j in 1:n_assets){   #K scenarios generation for each asset
    sigma_f_t1 <- tail(sigma_per[[j]],1) ##(t-1) sigma for GARCH forecasting
    e_f_t2_t1 <- tail(resid_per[[j]],2) ##(t-2, t-1) residuals for MA forecasting
    for(z in 1:nsim){
      ret_sim[z,j] = ((garch_coefs[[i]][[j]][[1]] * RZ[1260,j]) + (garch_coefs[[i]][[j]][[2]] * RZ[1259,j]) +  #AR1 * R_t-1, AR2 * R_t-2
                        (garch_coefs[[i]][[j]][[3]] * e_f_t2_t1[1]) + (garch_coefs[[i]][[j]][[4]] * e_f_t2_t1[2]) +  #MA1*e_t-1, MA2*e_t-2
                        (zsim[z,j] * (sqrt(garch_coefs[[i]][[j]][[9]]) +  ##alfa0
                                        sqrt(garch_coefs[[i]][[j]][[5]]) * e_f_t2_t1[2] + ##alfa1 * e_t-1
                                        sqrt(garch_coefs[[i]][[j]][[6]]) * sigma_f_t1))) ##beta1  * s_t-1
    }
  }
  return(ret_sim)
}
