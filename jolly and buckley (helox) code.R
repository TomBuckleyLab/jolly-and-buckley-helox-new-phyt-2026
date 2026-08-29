##### jolly and buckley (helox) code
#####
##### sample code to calculate intercellular relative humidity from
#####  gas exchange data collected using the helox compensation method,
#####  as described in Jolly and Buckley (2026) New Phytologist
#####
##### a sample data file (xxxx) is included for demonstration purposes
#####


### swvmf():
###   -function to calculate saturation water vapor mole fraction (mol/mol) from temperature (T, deg C)
swvmf <- function(T) {0.0061078*exp(17.27*T/(237.3 + T))}

### try_gs():
###
###  -takes input values of gasx variables
###    (En, Eh, gb1, gc1, stomatal ratio, wa and ws in helox and nitox)
###    (2nd arg is input vector iv = c(En, Eh, gb1, gc1, sr, wa, waprime, ws, wsprime))
###    as well as an 'input' estimate of stomatal conductance (gsin, first arg; nb this is one-sided abaxial gs)
###    and uses them to compute intercellular relative humidity (h),
###    then recalculates stomatal conductance ('output', gsout) from h,
###    and returns a vector containing 6 items:
###    (1) the squared difference between gsin and gsout,
###    (2-6) gsout itself, h, updated phi (gc1/gs), updated gamma (gs/gb1), and eta (gtwhelox/gtwnitox)
###  -optional args are the ratio of phi in helox and nitox (default = 1/2.33, assumes cuticular tpt is not gas diffusion)
###    and ratio of gamma in helox and nitox (default = 1.47, based on our chamber)
###  -calls calculate_gs() to get gs given gt, gb, gc
###  -is called by try_gs_err(), a dummy wrapper that returns only the squared error
###  -in the main calculation loop, optimize() minimizes the squared error in order to find the gs
###    that satisfies the system
###
try_gs <- function(gsin, iv, phiprime_phi = 1/2.33, gammaprime_gamma = 1.47) {
  ## extract inputs from input vector
  En = iv[1]
  Eh = iv[2]
  gb1 = iv[3]
  gc1 = iv[4]
  sr = iv[5]
  wa = iv[6]
  waprime = iv[7]
  ws = iv[8]
  wsprime = iv[9]
  
  ## estimate phi and gamma from gsin
  phi = gc1/gsin
  gamma = gsin/gb1
  
  ## calculate phi and gamma in helox
  phiprime = phi*phiprime_phi
  gammaprime = gamma*gammaprime_gamma
  
  ## calculate eta (ratio of gtw in helox/nitox)
  etanum = (1 + phiprime)/(gammaprime*(1 + phiprime) + 1) + (sr + phiprime)/(gammaprime*(sr + phiprime) + 1)
  etadenom = (1 + phi)/(gamma*(1 + phi) + 1) + (sr + phi)/(gamma*(sr + phi) + 1)
  eta = (En/Eh)*2.33*etanum/etadenom
  
  ## solve for h using quadratic
  q2 = (eta - 1)*wsprime*ws
  q1 = 2*(ws - eta*wsprime) + (eta + 1)*(wsprime*wa - ws*waprime)
  q0 = 2*(eta*waprime - wa) - (eta - 1)*wa*waprime
  h = (0.5/q2)*(-q1 - sqrt(q1*q1 - 4*q2*q0))
  
  ## recalculate gtw from h, then recalculate gs, phi and gamma
  gtout = En*(1 - 0.5*(h*ws + wa))/(h*ws - wa)
  gsout = calculate_gs(gtout, gb1, gc1, sr)
  phi = gc1/gsout
  gamma = gsout/gb1
  
  return(c((gsin - gsout)^2,
           gsout, h, phi, gamma, eta))
}

### try_gs_err():
###
###  -dummy wrapper for try_gs() that return only the first element of the return vector
###   (namely the squared difference between gsin and gsout)
###  -needed for optimize()
try_gs_err <- function(gsin, iv, phiprime_phi = 1/2.33, gammaprime_gamma = 1.47) {
  try_gs(gsin, iv, phiprime_phi, gammaprime_gamma)[1]
}  

### calculate_gs():
###
###  -invert eqn for gtw as f of gs, gb, gc, to calculate gs from gt, gb, gc
###  -gb, gc and gs here are one-sided values; gs specifically is abaxial
###  - sr is stomatal ratio (defined here as gs_adaxial/gs_abaxial)
### 
calculate_gs <- function(gt, gb1, gc1, sr) {
  y = gt/gb1
  x = gc1 + gb1
  q2 = sr*(2 - y)
  q1 = (sr + 1)*(gc1 + x - x*y)
  q0 = (2*gc1 - y*x)*x
  
  if(sr>0) {
    gs = (0.5/q2)*(-q1 + sqrt(q1*q1 - 4*q0*q2))
  } else {
    gs = -q0/q1    
  }
  return(gs)
}


############
############
### load in processed OPUS data
###

### set working directory
#setwd("...")
x <- read.csv("sample data.csv")

### input file has the following structure:
### data rows alternate between nitox and helox (e.g., row 1 = nitox, row 2 = helox)
### file has the following necessary columns:
###
### t - time of day (24-hour decimal hour) 
### E - leaf transpiration rate in mol m-2 s-1
### H2OS - sample stream water vapor mole fraction in mmol/mol
### TleafC - leaf temperature in degrees C
### gbw - two-sided (whole leaf) boundary layer conductance to water vapor, mol m-2 s-1
### gc - one-sided cuticular conductance (mol m-2 s-1)
### stomatal_ratio - ratio of adaxial to abaxial stomatal conductances


############
### make dataframe to contain processed results, with one row per N2-He switch
y <- data.frame(matrix(nrow=nrow(x)/2, ncol=24))
colnames(y) <- c("t", "gb2", "gc", "sr", 
                 "wsn", "wsh", "win", "wih", "wan", "wah",
                 "En", "Eh", "Tn", "Th", 
                 "gtn.est", "gsn.est", "gsn.corr", 
                 "phi", "gamma", "eta",
                 "hleaf.rough", "hleaf", "pct.dE", "deltaw")

## loop that runs through every pair of successive rows
for(p in 1:(nrow(x)/2)) {
  n = p*2-1 ## odd numbered rows starting from 1; this will be index for nitox data row
  h = n + 1 ## next row after n; this will be index for helox data row
  En = x$E[n]; Eh = x$E[h] ## transpiration rates in nitox and helox
  wan = 0.001*x$H2OS[n]; wah = 0.001*x$H2OS[h] ## sample (ambient) water vapor mole fractions 
  Tn = x$TleafC[n]; Th = x$TleafC[h] ## leaf temperatures
  wsn <- 0.0061078*exp(17.27*Tn/(237.3 + Tn)) ## saturated water vapor mole fractions
  wsh <- 0.0061078*exp(17.27*Th/(237.3 + Th))
  gb=x$gbw[n]
  gc=x$gc[n]
  sr=x$stomatal_ratio[n]
  gtn.est <- En*(1 - 0.5*(wsn + wan))/(wsn - wan) ## naive (saturated) estimation of total leaf conductance to water vapor in nitox
  gsn.est <- calculate_gs(gtn.est, gb/2, gc, sr) ## naive estimate of single-surface gsw in nitox

  iv=c(En, Eh, gb/2, gc, sr, wan, wah, wsn, wsh) ## vector of inputs to send to try_gs_err
  
  ### guts of the calculations:
  ###  - iteratively solve for value of gsw that is consistent with Eqns 5 and 8 in the main text, given the data
  ###  - unsaturation means gsn.est is an underestimate of true gsn, so lower-bound the search at gsn.est
  o <- optimize(try_gs_err, interval=c(gsn.est, 5*gsn.est), iv=iv) ## find solution
  rvec <- try_gs(o$minimum, iv) ## calculate resulting values of gsw in nitox, intercellular RH, phi, gamma and eta
  gsn.corr <- rvec[2] ## gsw in nitox corrected for unsaturation
  hleaf <- rvec[3] ## intercellular relative humidity
  phi <- rvec[4] ## ratio of gcw/gsw
  gamma <- rvec[5] ## ratio of gsw/gbw
  eta <- rvec[6] ## ratio of gtw in helox to that in nitox
  
  ## record everything in the dataframe y
  y$t[p] <- x$t[n]
  y$gb2[p] <- gb
  y$gc[p] <- gc
  y$sr[p] <- sr
  y$wsn[p] <- wsn
  y$wsh[p] <- wsh
  y$win[p] <- hleaf*wsn
  y$wih[p] <- hleaf*wsn
  y$wan[p] <- wan
  y$wah[p] <- wah
  y$En[p] <- En
  y$Eh[p] <- Eh
  y$Tn[p] <- Tn
  y$Th[p] <- Th
  y$gtn.est[p] <- gtn.est
  y$gsn.est[p] <- gsn.est
  y$gsn.corr[p] <- gsn.corr
  y$phi[p] <- phi
  y$gamma[p] <- gamma
  y$eta[p] <- eta
  y$hleaf.rough[p] <- (2.33*wah - wan)/(2.33*wsh - wsn) ## 'back of the envelope' estimate of hleaf (Eqn 9 in main text)
  y$hleaf[p] <- hleaf
  y$pct.dE[p] <- 100*(y$Eh[p] - y$En[p])/y$En[p] ## % change in E in helox vs nitox
  y$deltaw[p] <- 1000*(y$win[p] - y$wan[p]) ## deltaw in mmol/mol
}

