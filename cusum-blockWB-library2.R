
######################################################################
# CUSUM test with the given bandwidth q for Bartlett estimator
# return CUSUM statistic / p-value
######################################################################
cusum.bartlett <-function(data, q){
  n = length(data);
  if(missing(q)){ q = n^(1/3)};
  
  # if(missing(q)){ 
  #   rhohat =  ar.ols(data, FALSE, order=1);
  #   rhohat = rhohat$ar[1];
  #   q = 1.1447*(n*4*rhohat^2/(1-rhohat^2)^2)^(1/3);}
  #   q = floor(q);
  
  # Find break point
  gmean = mean(data);
  x1 = cumsum(data);
  j=seq(from=1, to=n, by=1);
  x2 = abs(x1 - j*gmean);
  khat = which.max(x2);
  kv = x2[khat];
  # Estimate long-run variance
  sn2 = longrun_bartlett(data, q);
  tn = kv/(sqrt(n*sn2));
  pval = sup_abs_bbridge(tn, 1);
  out = list();
  out$tn = tn;
  out$pval = pval;
  out$khat = khat;
  out$kv = kv;
  out$sn2 = sn2;
  out$q = q;
  return(out)
  
}

#################################################
# P-value calculation for Sup of Brownian Bridge
#################################################
sup_abs_bbridge <-function(tstat, nmax){
  # Calculate p-value for given test statistics
  # with nmax-number of max(sup|B_0(t)|) 
  # For one sup|B_0(t)| the formula is given by
  # 1 + 2*sum_{k=1}^{\infty} (-1}^k exp(-2*k^2*v^2)
  # Resnick, Adventures in Stochastic Processes p.533
  # default for nmax =1
  if(missing(nmax)){ nmax = 1};
  K = 100;
  k = seq(from=1, to=K, by=1);
  cumd=1+2*sum(((-1)^k)*exp(-2*tstat^2*k^2));
  pval = 1 - (cumd)^nmax;
  return(pval)
}


################################################### 
## Bartlett long-run variance calculation
###################################################
longrun_bartlett <-function(data, q){
  n = length(data);
  if(missing(q)){ q = n^(1/3)};
  # if(missing(q)){ 
  #   rhohat =  ar.ols(data, FALSE, order=1);
  #   rhohat = rhohat$ar[1];
  #   q = 1.1447*(n*4*rhohat^2/(1-rhohat^2)^2)^(1/3);
  # } # Data dependent bw
 
  if(q ==0){
    sn2 = var(data);
  } else{
  q = floor(q+1);
  wq = c(seq(from=1, to=q, by=1), q+1, seq(from=q, to=1, by=-1));
  xcov = acf(data, lag.max = q, type = "covariance", plot = FALSE, na.action = na.fail, demean = TRUE)
  sn2 = xcov$acf;
  sn2 = as.numeric(sn2);
  id2 = seq(from=q+1, to=2, by=-1)
  sn2 = c(sn2[id2], sn2)
  sn2 = sum(wq*sn2)/(q+1);
  }
  return(sn2);
}

#######################################
# Auxiliary function for block length
blocklength = function(data){
  n = length(data);
  a1 = acf(data, lag=n/4, plot = FALSE);
  a1 = as.numeric(a1$acf);
  a1 = a1[-1];
  m = which(abs(a1) <= 2/sqrt(n))[1];
  m =ifelse(is.na(m), 25, m)
  return(2*m+1)
}


######################################
# WB block length
# Using 2h rule in ACF
######################################
WBlength = function(rt){
  
  ## for mean chagne
  ## 2h rule for block length;
  t = length(rt);
  a1 = acf(rt, lag=t/4, plot = FALSE);
  a1 = as.numeric(a1$acf);
  a1 = a1[-1];
  m = which(abs(a1) <= 2/sqrt(t))[1];
  m1_m =ifelse(is.na(m), 25, m)
  
  ## for variance chagne
  
  khat = cusum.bartlett(rt)$khat
  
  rt_1 = c(rt[1:(khat)]-mean(rt[1:(khat)]),rt[(khat+1):t]-mean(rt[(khat+1):t]))	
  
  a1 = acf(abs(rt_1)^2, lag=t/4, plot = FALSE);
  a1 = as.numeric(a1$acf);
  a1 = a1[-1];
  m = which(abs(a1) <= 2/sqrt(t))[1];	
  m1_v = ifelse(is.na(m), 25, m)   
  
  
  n.b_m = floor(2*m1_m);  
  n.b_v = floor(2*m1_v);   
  
return(c(n.b_m, n.b_v))
}



cusum.bartlett.WB = function(data,boots,n.block) {
if(missing(boots)){ boots = 199 }
  t <- length(data)
  n.bl = n.block
  
  CUSUM.WB = numeric(boots)
  reject.WB = numeric(boots)
  
  ## CUSUM test 
  
  test <- cusum.bartlett(data);
  CUSUM = test$tn
  
  ## Wild bootstrap CUSUM test and its p-value
  w <- c(-1,1)  ## Rademacher distribution
  #w <- c(-(sqrt(5)-1)/2,(sqrt(5)+1)/2); p <- c((sqrt(5)+1)/(2*sqrt(5)),(sqrt(5)-1)/(2*sqrt(5)))
  
  CUSUM.WB <- numeric(boots)
  
  for(i in 1:boots){
    wt <- rep(sample(w,(t/n.bl+1),replace=TRUE),each=n.bl)[1:t]
    data.WB <- (data-mean(data))*wt
    CUSUM.WB[i] <- cusum.bartlett(data.WB)$tn
  }
  
  pval = sum(CUSUM.WB >= CUSUM)/boots
  try.reject <- try(ifelse(quantile(CUSUM.WB, prob=c(0.95))<CUSUM,1,0),TRUE)
  
  out = list();
  out$CUSUM.WB = CUSUM.WB
  out$cusumstat = CUSUM;
  out$cusumpval = test$pval;
  out$pval = pval;
  out$reject  =  try.reject;
  return(out)
}



########################################################################################
# Garch(1,1) with change parameters generation
# Generating GARCH(1,1) process with a parameter change. (cf. Kokoszka and Leipus (2000))
#########################################################################################

garchSimWithBreak <- function( n , omega1 , alpha1 , beta1 , omega2 , alpha2 , beta2 , k ) {
  # We generate two GARCH processes.
  # Note that the two processes are driven by a single sequence of i.i.d. innovations,
  # so, we can not use garchSim() directly.
  # The change of parameter occurs at time = k
  
  rt1 <- double(n)
  rt2 <- double(n)
  
  # Burning initial conditional variances
  h1 <- omega1 / ( 1 - alpha1 - beta1 )
  h2 <- omega2 / ( 1 - alpha2 - beta2 )
  xi <- rnorm( n + 200 )
  for( i in 1:200 ) {
    h1 <- omega1 + ( alpha1 * xi[i]^2 + beta1 ) * h1
    h2 <- omega2 + ( alpha2 * xi[i]^2 + beta2 ) * h2
  }
  
  
  for( i in 1:n ) {
    rt1[i] <- xi[i+200] * sqrt( h1 )
    rt2[i] <- xi[i+200] * sqrt( h2 )
    h1 <- omega1 + alpha1 * rt1[i]^2 + beta1 * h1
    h2 <- omega2 + alpha2 * rt2[i]^2 + beta2 * h2
  }
  
  rt <- c( rt1[1:k] , rt2[(k+1):n] )
  
  return(rt)
}


########################################################################################
# AR-Garch(1,1) with change parameters generation
#########################################################################################

ARgarchSimWithBreak <- function( n , phi, mu1, mu2, omega1 , alpha1 , beta1 , omega2 , alpha2 , beta2 , k ) {

  # Note that the two processes are driven by a single sequence of i.i.d. innovations,
  # so, we can not use garchSim() directly.
  # The change of parameter occurs at time = k

  innov = garchSimWithBreak(n+100, omega1, alpha1, beta1, omega2, alpha2, beta2, k);
  z = filter(innov, phi, method = "recursive");
  y = z[-(1:100)];
  y1 = y[1:k] + mu1;
  y2 = y[-(1:k)] + mu2;
  y = c(y1, y2);
return(y)
}



J.cusum = function(rt, boots){
  if(missing(boots)){ boots = 199}
  
  t = length(rt);
  khat = cusum.bartlett(rt)$khat
  rt_1 = c(rt[1:(khat)]-mean(rt[1:(khat)]),rt[(khat+1):t]-mean(rt[(khat+1):t]))	
  rt_2 = (rt_1)^2 - mean((rt_1)^2)
  fit1 = cusum.bartlett(rt);
  fit2 = cusum.bartlett(rt_2);
  J_CUSUM <- max(fit1$tn,fit2$tn);
  J_cusum.b <- min(fit1$pval, fit2$pval)
  
  ## 2*h rule for LRV/WB length
  mm = WBlength(rt);
  q1 = floor((mm[1]-1)/2);   q2 = floor((mm[2]-1)/2); 
  fit3 = cusum.bartlett(rt, q1);
  fit4 = cusum.bartlett(rt_2, q2);
  J_CUSUM.mm <- max(fit3$tn,fit4$tn);
  J_cusum.b.mm <- max(fit3$pval, fit4$pval)
  
  
  ## Andrew selection
  mAnd = 2*c(fit1$q, fit2$q)+1;
  out1 = J.cusum.boot(rt, rt_2, mAnd[1], mAnd[2], boots);
  pval1 = c(sum(out1$bstat1 > J_CUSUM))/boots; 

  ## 2*h rule for LRV/WB length
  out2 = J.cusum.boot(rt, rt_2, mm[1], mm[2], boots);
  pval2 = c(sum(out2$bstat1 > J_CUSUM.mm))/boots; 
  
#  ## Aggregation of several bandwidth
#  out3 = cusum.boot(rt, rt_2, M1, M2, boots);
#  pval5 = c(sum(out3$bstat1 > J_CUSUM),  sum(out3$bstat2 > J_CUSUM))/length(out3$bstat1); 
#  pval6 = c(sum(out3$bstat1 > J_CUSUM.mm),  sum(out3$bstat2 > J_CUSUM.mm))/length(out3$bstat1); 

  ## Fixed bandwidth as 2*sqrt(T)
  out3 = J.cusum.boot(rt, rt_2, floor(2*sqrt(t)), floor(2*sqrt(t)), boots);
  pval3 = c(sum(out1$bstat1 > J_CUSUM))/boots; 
  
  
  
  ## Fixed bandwidth as 2*log(T)
  out4 = J.cusum.boot(rt, rt_2, floor(2*log(t)), floor(2*log(t)), boots);
  pval4 = c(sum(out1$bstat1 > J_CUSUM))/boots; 
    
  out = list();
  
  out$J_CUSUM = c(J_CUSUM, J_CUSUM.mm);
  out$CUSUM.pval = c(J_cusum.b, J_cusum.b.mm);
  
  out$mAnd = mAnd;
  out$pval1 = pval1;

  out$mm = mm;
  out$pval2 = pval2;
  
  out$sqrt = floor(2*sqrt(t));
  out$pval3 = pval3;

  out$log = floor(2*log(t));
  out$pval4 = pval4;
  out$khat = ifelse(fit1$tn > fit2$tn, fit1$khat, fit2$khat);
  out$mean = ifelse(fit1$tn > fit2$tn, 1, 0);
  return(out)
}


J.cusum.boot = function(rt, rt_2, m1, m2, boots){
  
  if(missing(boots)){ boots = 199}
  #if(missing(J_CUSUM)){ J_CUSUM <- max(cusum.bartlett(rt)$tn,cusum.bartlett(rt_2)$tn)}
  m1 = floor(m1); m2 = floor(m2);
  t = length(rt);
  ## Apply BW for pvalue
  w <- c(-1,1)  ## Rademacher distribution
  #w <- c(-(sqrt(5)-1)/2,(sqrt(5)+1)/2); p <- c((sqrt(5)+1)/(2*sqrt(5)),(sqrt(5)-1)/(2*sqrt(5)))
  J_CUSUM.WB.var = J_CUSUM.WB = CUSUM.WB <- c()
  n.bl = m1; n.bl_2 = m2; q1 = floor((m1-1)/2); q2 = floor((m2-1)/2);
  
  for(i in 1:boots){
    
    #wt <- sample(w,t,replace=TRUE)
    #wt <- rep(sample(w,(t/n.bl+1),replace=TRUE,prob=p),each=n.bl)[1:t]
    wt <- rep(sample(w,(t/n.bl+1),replace=TRUE),each=n.bl)[1:t]
    data.WB <- rt*wt
    wt_2 <- rep(sample(w,(t/n.bl_2+1),replace=TRUE),each=n.bl_2)[1:t]
    data.WB_2 <- rt_2*wt_2   
    
    J_CUSUM.WB[i] <- max(cusum.bartlett(data.WB, q1)$tn,cusum.bartlett(data.WB_2, q2)$tn)
#    J_CUSUM.WB[i] <- max(cusum.bartlett(data.WB)$tn,cusum.bartlett(data.WB_2)$tn)
    J_CUSUM.WB.var[i] <- max(cusum.bartlett(data.WB, 0)$tn,cusum.bartlett(data.WB_2, 0)$tn)
  }
  
  return(list(bstat1 = J_CUSUM.WB, bstat2 = J_CUSUM.WB.var))
}




########################################################
# Binary segmentation to get changepoints
########################################################


JCUSUM.multi <- function(rt, boots, type){
  if(missing(type)){ type = c("log")}
  if(missing(boots)){ boots = 199}
  n = length(rt);
  
  Del =80; 
  out = list();
#  set.seed(13579)
  br = newbr = brhist = c(0, n); st=1; 
  TF = 1; pvalhist = tstahist = meanhist = NULL;
  nmax = 1; 
  while(st > 0 && nmax < 15){
    newTF = NULL;
    for(i in 1:(length(br)-1)){
      if(TF[i] > 0 && br[i+1] - br[i] > 2*Del){
        id = seq(from=br[i]+1, to = br[i+1], by=1);
        zz = rt[id];
        out = J.cusum(zz, boots)
        khat = out$khat;
        newk = khat + br[i];
        
        if(type == "andrew") { pval = out$pval1 };
        if(type == "politis") { pval = out$pval2 };
        if(type == "sqrt") { pval = out$pval3 };
        if(type == "log") { pval = out$pval4 };
        
        pvalhist = c(pvalhist, pval);
        tstahist = c(tstahist, out$J_CUSUM[1]);
        
        if( pval < .10 && khat > Del && (br[i+1] - newk) > Del){
          newbr = c(newbr, newk); 
          newTF = c(newTF, 1, 1); 
          brhist = c(brhist, newk);
#         pvalhist = rbind(pvalhist, pval);
#          tstahist = c(tstahist, out$J_CUSUM[1]);
          meanhist = c(meanhist, out$mean);
        } else{ newTF = c(newTF, 0); };
      } else{ newTF = c(newTF, 0);};
    }
    br = unique(newbr); br = sort(br);
    TF = newTF;
    st = sum(TF);
    nmax = nmax +1;
  }
  
  out = list();
  out$br = br;
  out$brhist = brhist;
  out$pvalhist = pvalhist;
  out$tstahist = tstahist;
  out$meanhist = meanhist;
  return(out)
}





########################################################
# BWB mean change - multi
########################################################


CUSUM.meanbw_multi <- function(rt, boots, type){
  if(missing(type)){ type = c("andrew")}
  if(missing(boots)){ boots = 199}
  n = length(rt);
  
  Del =200; 
  out = list();
#  set.seed(13579)
  br = newbr = brhist = c(0, n); st=1; 
  TF = 1; pvalhist = tstahist = meanhist = NULL;
  nmax = 1; 
  while(st > 0 && nmax < 15){
    newTF = NULL;
    for(i in 1:(length(br)-1)){
      if(TF[i] > 0 && br[i+1] - br[i] > 2*Del){
        id = seq(from=br[i]+1, to = br[i+1], by=1);
        zz = rt[id]; 
        t = length(zz);
        fit = cusum.bartlett(zz);
        khat = fit$khat
        #rt_1 = c(rt[1:(khat)]-mean(rt[1:(khat)]),rt[(khat+1):t]-mean(rt[(khat+1):t]))	
        # rt_2 = (rt_1)^2 - mean((rt_1)^2)
        # fit2 = cusum.bartlett(rt_2);
        #J_CUSUM <- max(fit1$tn,fit2$tn);
        #J_cusum.b <- min(fit1$pval, fit2$pval)
#        if(type == "andrew") { m1= fit$q };
        if(type == "andrew") { m1= floor(fit$q) };
        if(type == "politis") { m1 =  blocklength(zz); };
        if(type == "sqrt") { m1 = floor(1.5*sqrt(t)); };
        if(type == "log") { m1  = floor(1.5*log(t)); };
        
        pval = sum(cusum.mean.boot(zz, m1, boots)$bstat1 > fit$tn)/boots;
        newk = khat + br[i];
        
        pvalhist = c(pvalhist, pval);
        tstahist = c(tstahist, fit$tn);
        
        if( pval < .05 && khat > Del && (br[i+1] - newk) > Del){
          newbr = c(newbr, newk); 
          newTF = c(newTF, 1, 1); 
          brhist = c(brhist, newk);
        #  pvalhist = rbind(pvalhist, pval);
        #  tstahist = c(tstahist, fit$tn);
          meanhist = c(meanhist, out$mean);
        } else{ newTF = c(newTF, 0); };
      } else{ newTF = c(newTF, 0);};
    }
    br = unique(newbr); br = sort(br);
#    newbr = br;
    TF = newTF;
    st = sum(TF);
    nmax = nmax +1;
  }
  
  out = list();
  out$br = br;
  out$brhist = brhist;
  out$pvalhist = pvalhist;
  out$tstahist = tstahist;
  out$fitted = fit.breaks(rt, br);
  return(out)
}


fit.breaks = function(rt, br){
  fitted=0;
  for(i in 1:(length(br)-1)){
    id = seq(br[i]+1, br[i+1], 1);
    y1 = rep(mean(rt[id]), length(id)); 
    fitted = c(fitted, y1);
  }
  fitted = fitted[-1];
  return(fitted);
}



cusum.boot = function(rt, rt_2, M1, M2, boots){
  if(missing(boots)){ boots = 199}
  
  lm1 = length(M1); lm2 = length(M2);
  
  Bstat1 = Bstat2 = NULL;
  for(i in 1:lm1){
    for(j in 1:lm2){
      m1 = M1[i]; m2 = M2[j];
      b1 = cusum.mean.boot(rt, m1, boots);
      b2 = cusum.var.boot(rt_2, m2, boots);
      
      Bstat1= rbind(Bstat1, apply(cbind(b1$bstat1, b2$bstat1), 1, max));
      Bstat2= rbind(Bstat2, apply(cbind(b1$bstat2, b2$bstat2), 1, max));
      
    }
  }
  
  Bstat1.avg = as.vector(Bstat1);
  Bstat2.avg = as.vector(Bstat2);
  return(list(bstat1 = Bstat1.avg, bstat2 = Bstat2.avg))
}


cusum.mean.boot = function(rt, m1, boots){
  if(missing(boots)){ boots = 1999}
  t = length(rt);
  w <- c(-1,1)  ## Rademacher distribution
  #w <- c(-(sqrt(5)-1)/2,(sqrt(5)+1)/2); p <- c((sqrt(5)+1)/(2*sqrt(5)),(sqrt(5)-1)/(2*sqrt(5)))
  CUSUM.WB = CUSUM.WB.var = numeric(boots);
  n.bl = m1; 
  for(i in 1:boots){
    wt <- rep(sample(w,(t/n.bl+1),replace=TRUE),each=n.bl)[1:t]
    data.WB <- (rt-mean(rt))*wt
    CUSUM.WB[i] <- cusum.bartlett(data.WB)$tn;
#    CUSUM.WB[i] <- cusum.bartlett(data.WB)$tn;
    CUSUM.WB.var[i] <- cusum.bartlett(data.WB, 0)$tn;
  }
  
  #p1 = sum(J_CUSUM.WB >  J_CUSUM)/boots;
  #p2 = sum(J_CUSUM.WB.var >  J_CUSUM)/boots;
  #return(c(p1, p2))
  return(list(bstat1 = CUSUM.WB, bstat2 = CUSUM.WB.var))
  
}



cusum.var.boot = function(rt_2, m2, boots){
  if(missing(boots)){ boots = 199}
  w <- c(-1,1)  ## Rademacher distribution
  #w <- c(-(sqrt(5)-1)/2,(sqrt(5)+1)/2); p <- c((sqrt(5)+1)/(2*sqrt(5)),(sqrt(5)-1)/(2*sqrt(5)))
  CUSUM.WB = CUSUM.WB.var = numeric(boots);
  n.bl_2 = m2; q2 = floor((m2-1)/2);
  for(i in 1:boots){
    wt_2 <- rep(sample(w,(t/n.bl_2+1),replace=TRUE),each=n.bl_2)[1:t]
    data.WB_2 <- rt_2*wt_2   
    CUSUM.WB[i] <- cusum.bartlett(data.WB_2, q2)$tn;
#    CUSUM.WB[i] <- cusum.bartlett(data.WB_2)$tn;
    CUSUM.WB.var[i] <- cusum.bartlett(data.WB_2, 0)$tn;
  }
  
  return(list(bstat1 = CUSUM.WB, bstat2 = CUSUM.WB.var))
}




######################################################################
# CUSUM Multi with Bartlett estimator
######################################################################


cusum_multi <- function(data){
  T = length(data);
  result=list();

  ###############################################
  # Find the break point by CUSUM statistic
  ###############################################
  
  find_break <- function(data){
    n = length(data)
    gmean = mean(data);
    x1 = cumsum(data);
    j=seq(from=1, to=n, by=1);
    x2 = abs(x1 - j*gmean);
    khat = which(x2 == max(x2));
    value = x2[khat];
    out = cbind(khat, value);
    return(out);
  }
  
  
  # Step0: No break
  out = cusum.bartlett(data);
  result = out;
  result$type = NULL;
  result$kv = NULL;
  nmax = 1;
  # Step1 : 1 break
  if(out$pval < .1){
    nmax = nmax+1;
    khat = out$khat;
    dat1 = data[1:khat];
    id = seq(from=khat+1, to=T, by=1); 
    dat2 = data[id]; 
    out1 = cusum.bartlett(dat1);
    out2 = cusum.bartlett(dat2);
    tstat = c(out1$tn,out2$tn);
    br = c(0, khat, T);
    
    pval = sup_abs_bbridge(max(tstat), nmax)
    result$pval = c(result$pval, pval);
    result$tn = c(result$tn, max(tstat));
    result$khat = c(result$khat, br);
    result$sn2 = c(result$sn2, out1$sn2, out2$sn2);
    result$bw = c(result$bw, out1$bw, out2$bw);
    result$br = br;
    
    # Step1 : 2 breaks and more
    while (pval < .05){
      nmax = nmax+1;
      id = which(tstat == max(tstat));
      id2 = seq(from = br[id]+1, to = br[id+1], by =1);
      subdata = data[id2];
      # Apply BHKS again here
      br_out = find_break(subdata);
      khat = br_out[1];
      subdata1 = subdata[1:khat];
      id3 = seq(from=khat+1, to=length(subdata), by=1); 
      subdata2 = subdata[id3]; 
      out1 = cusum.bartlett(subdata1);
      out2 = cusum.bartlett(subdata2);
      tstat = append(tstat[-id], c(out1$tn, out2$tn), after=(id-1));
      br = append(br, br[id]+khat, after=id);
      pval = sup_abs_bbridge(max(tstat), nmax);
      
      result$pval = c(result$pval, pval);
      result$tn = c(result$tn, max(tstat));
      result$khat = c(result$khat, br);
      result$sn2 = c(result$sn2, out1$sn2, out2$sn2);
      result$bw = c(result$bw, out1$bw, out2$bw);
      result$br = br;
    }
    fitted=0;
    for(i in 1:(length(br)-1)){
      id = seq(br[i]+1, br[i+1], 1);
      y1 = rep(mean(data[id]), length(id)); 
      fitted = c(fitted, y1);
    }
    fitted = fitted[-1];
    result$fitted = fitted;
  }
  else{ result$fitted = rep(mean(data), T); }
  result$iter = nmax;
  return(result)
}



