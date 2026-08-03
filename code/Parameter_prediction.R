# Kinetic parameter prediction
mu_data  <- exp_data[1]
phi_data <- exp_data[-1]
cm0 <-  as.numeric(rho_cond[1]*M%*%q0_alt)[-p]
esat <- predict.parameters
rounding <- 1
kcat_ratio <- 3
default_kcat <- 10

# Estimates Km values assuming typical ratio Km = cm/esat for substrates and 

# forces some reactions to be irreversible if kcatb = 0 in the file
rev <- as.vector(1*(kcatb > 0))

K[(n_a+1):(p+n_a-1),] <- round(diag(cm0/esat)%*%(1*(M[-p,] < 0)) + is.reversible*diag(cm0)%*%(1*(M[-p,] > 0))%*%diag(rev),1)


# Estimates total saturation of each reaction assuming irrev. with cm = esat*Km
sat <- (esat/(esat + 1))^(colSums(1*(M[-p,] < 0)))


# assumes transporters are completely saturated
sat[1] <- 1

# Estimates now kcatf based on ...
b0 <- M%*%q0_alt
kcatf <- round(mu_data*q0_alt/(b0[p]*phi_data*sat), rounding)

# if some kcat was rounded to zero
kcatf[kcatf == 0] <- default_kcat
kcatf[kcatf == Inf] <- default_kcat
kcatf[is.nan(kcatf)] <- default_kcat


# estimating kcatb
kcatb <- round(kcatf/kcat_ratio, rounding)
kcatb[kcatb == 0] <- default_kcat/kcat_ratio
kcatb <- is.reversible*rev*kcatb
# ribosome is irreversible
kcatb[r] <- 0

# first separate Km matrices for substrates and for products
KS <- K*(Mtotal<0)
KP <- K*(Mtotal>0)

KS[KS == 0] <- Inf
if(is.reversible){
  KP[KP == 0] <- Inf
}else{
  KP[] <- Inf
}
KI[KI == 0] <- Inf


# scale kcats up or down until the predicted growth matches the experimental
if(rescale_kcats){
  tol <- 0.02
  n_conditions <- 1
  source(here("code", "GBA_solver.R"))

  print(paste0("Initial mu = ", mu_opt))
  n = 0
  if(abs(mu_opt - mu_data) > tol){
    print(paste0(" -> scaling kcats to mu = ", mu_data))
    
    while(abs(mu_opt - mu_data) > tol){
      n = n+1
      
      if(mu_opt < mu_data){
        kcatf <- round(kcatf*1.05, rounding)
      }else{
        kcatf <- round(kcatf*0.95, rounding)
      }
      
      kcatb <- round(kcatf/kcat_ratio, rounding)
      kcatb[kcatb == 0] <- 1
      kcatb <- is.reversible*rev*kcatb
      kcatb[r] <- 0 # ribosome is always irreversible

      source(here("code", "GBA_solver.R"))
      print(paste0("mu after ", n, " rounds = ", mu_opt))
      
      if(n > 40){
        break
      }
    }
  }
  
  source(here("code", "GBA_solver.R"))
  print(paste0("final mu after ", n, " rounds = ", mu_opt))
  
  n_conditions <- length(rho_cond)  
}



