BNLPWA <- function(data, func=NULL, method=c("mixture","mom","imom"), mixprob_dist=c("logit","genlogit","hypsec","gennormal"), scale_dist=c("polynom","doubleexp"), r=1, nu=1, wave.family="DaubLeAsymm", filter.number=6, bc="periodic")
{
	time = Sys.time()

  #...............Wavelet decomposition..................
  
  n <- length(data)
  wave_obj <- wavelet_trans(data,filter.number,wave.family,bc)
  dlj <- wave_obj$dlj
  sigmasq <- wave_obj$sigmasq


  #............. Hyperparameter specification ...........

  if(method=="mixture")
  {
		if(mixprob_dist=="logit")
		{
			gamma1.l_func <- function(C,l){ exp(C[1]-C[2]*l)/( 1+exp(C[1]-C[2]*l) ) }
			gamma2.l_func <- function(C,l){ exp(C[3]-C[4]*l)/( 1+exp(C[3]-C[4]*l) ) }
			switch(scale_dist,
				"polynom" = {
					tau1.l_func <- function(C,l){ C[5] * l^(-C[6]) }
			  	tau2.l_func <- function(C,l){ C[7] * l^(-C[8]) }
			  	npar_mixprob <- 4; npar_scale <- 4
				},
				"doubleexp" = {
					tau1.l_func <- function(C,l){ C[5] * exp(-C[6]*l) + C[7] * exp(-C[8]*l) }
				  tau2.l_func <- function(C,l){ C[9] * exp(-C[10]*l) + C[11] * exp(-C[12]*l) }
				  npar_mixprob <- 4; npar_scale <- 8
				}
			)
		} else
		if(mixprob_dist=="genlogit")
		{
			gamma1.l_func <- function(C,l){ (1+exp(-C[1]+C[2]*l))^{-C[3]} }
			gamma2.l_func <- function(C,l){ (1+exp(-C[4]+C[5]*l))^{-C[6]} }
			switch(scale_dist,
				"polynom" = {
					tau1.l_func <- function(C,l){ C[7] * l^(-C[8]) }
			  	tau2.l_func <- function(C,l){ C[9] * l^(-C[10]) }
			  	npar_mixprob <- 6; npar_scale <- 4
				},
				"doubleexp" = {
					tau1.l_func <- function(C,l){ C[7] * exp(-C[8]*l) + C[9] * exp(-C[10]*l) }
				  tau2.l_func <- function(C,l){ C[11] * exp(-C[12]*l) + C[13] * exp(-C[14]*l) }
				  npar_mixprob <- 6; npar_scale <- 8
				}
			)
		} else
		if(mixprob_dist=="hypsec")
		{
			gamma1.l_func <- function(C,l){ (2/pi) * atan(exp((pi/2)*(C[1]-C[2]*l))) } #(C[1]-l)/C[2]))
			gamma2.l_func <- function(C,l){ (2/pi) * atan(exp((pi/2)*(C[3]-C[4]*l))) } #(C[3]-l)/C[4]))
			switch(scale_dist,
				"polynom" = {
					tau1.l_func <- function(C,l){ C[5] * l^(-C[6]) }
  				tau2.l_func <- function(C,l){ C[7] * l^(-C[8]) }
  				npar_mixprob <- 4; npar_scale <- 4
				},
				"doubleexp" = {
					tau1.l_func <- function(C,l){ C[5] * exp(-C[6]*l) + C[7] * exp(-C[8]*l) }
				  tau2.l_func <- function(C,l){ C[9] * exp(-C[10]*l) + C[11] * exp(-C[12]*l) }
				  npar_mixprob <- 4; npar_scale <- 8
				}
			)
		} else
		if(mixprob_dist=="gennormal")
		{
			gamma1.l_func <- function(C,l){ 0.5 + sign(C[1]-l) * (1/(2*gamma(1/C[2]))) * gamma(1/C[2])*pgamma(abs((C[1]-l)/C[3])^C[2],shape=1/C[2],lower.tail=TRUE) }
			gamma2.l_func <- function(C,l){ 0.5 + sign(C[4]-l) * (1/(2*gamma(1/C[5]))) * gamma(1/C[5])*pgamma(abs((C[4]-l)/C[6])^C[5],shape=1/C[5],lower.tail=TRUE) }
			switch(scale_dist,
				"polynom" = {
					tau1.l_func <- function(C,l){ C[7] * l^(-C[8]) }
  				tau2.l_func <- function(C,l){ C[9] * l^(-C[10]) }
  				npar_mixprob <- 6; npar_scale <- 4
				},
				"doubleexp" = {
					tau1.l_func <- function(C,l){ C[7] * exp(-C[8]*l) + C[9] * exp(-C[10]*l) }
				  tau2.l_func <- function(C,l){ C[11] * exp(-C[12]*l) + C[13] * exp(-C[14]*l) }
				  npar_mixprob <- 6; npar_scale <- 8
				}
			)
		}
	} else
	if(method %in% c("mom","imom"))
	{
		if(mixprob_dist=="logit")
		{
			gamma.l_func <- function(C,l){ exp(C[1]-C[2]*l)/( 1+exp(C[1]-C[2]*l) ) }
			switch(scale_dist,
				"polynom" = {
					tau.l_func <- function(C,l){ C[3] * l^(-C[4]) }
	  			npar_mixprob <- 2; npar_scale <- 2
				},
				"doubleexp" = {
					tau.l_func <- function(C,l){ C[3] * exp(-C[4]*l) + C[5] * exp(-C[6]*l) }
	  			npar_mixprob <- 2; npar_scale <- 4
				}
			)
		} else
		if(mixprob_dist=="genlogit")
		{
			gamma.l_func <- function(C,l){ (1+exp(-C[1]+C[2]*l))^{-C[3]} }
			switch(scale_dist,
				"polynom" = {
					tau.l_func <- function(C,l){ C[4] * l^(-C[5]) }
	  			npar_mixprob <- 3; npar_scale <- 2
				},
				"doubleexp" = {
					tau.l_func <- function(C,l){ C[4] * exp(-C[5]*l) + C[6] * exp(-C[7]*l) }
	  			npar_mixprob <- 3; npar_scale <- 4
				}
			)
		} else
		if(mixprob_dist=="hypsec")
		{
			gamma.l_func <- function(C,l){ (2/pi) * atan(exp((pi/2)*(C[1]-C[2]*l))) } #(C[1]-l)/C[2]))
			switch(scale_dist,
				"polynom" = {
					tau.l_func <- function(C,l){ C[3] * l^(-C[4]) }
	  			npar_mixprob <- 2; npar_scale <- 2
				},
				"doubleexp" = {
					tau.l_func <- function(C,l){ C[3] * exp(-C[4]*l) + C[5] * exp(-C[6]*l) }
	  			npar_mixprob <- 2; npar_scale <- 4
				}
			)
		} else
		if(mixprob_dist=="gennormal")
		{
			gamma.l_func <- function(C,l){ 0.5 + sign(C[1]-l) * (1/(2*gamma(1/C[2]))) * gamma(1/C[2])*pgamma(abs((C[1]-l)/C[3])^C[2],shape=1/C[2],lower.tail=TRUE) }
			switch(scale_dist,
				"polynom" = {
					tau.l_func <- function(C,l){ C[4] * l^(-C[5]) }
	  			npar_mixprob <- 3; npar_scale <- 2
				},
				"doubleexp" = {
					tau.l_func <- function(C,l){ C[4] * exp(-C[5]*l) + C[6] * exp(-C[7]*l) }
	  			npar_mixprob <- 3; npar_scale <- 4
				}
			)
		}
	}

	#................. Log-likelihood function ..............

	minusll <- switch(
		method,
		"mixture" = 
		function(...,n,dlj,r,nu,sigmasq)
	  {
	    aux <- 0
	    C <- c(...)[order( as.numeric(gsub("\\D","",names(c(...)))) )] #extract digits from names as numeric, order columns by that
	    gamma1_l_list <- sapply(3:(log2(n)-1), function(l) gamma1.l_func(C, l))
			gamma2_l_list <- sapply(3:(log2(n)-1), function(l) gamma2.l_func(C, l))
			tau1_l_list <- sapply(3:(log2(n)-1), function(l) tau1.l_func(C, l))
			tau2_l_list <- sapply(3:(log2(n)-1), function(l) tau2.l_func(C, l))
	    for(l in 3:(log2(n)-1))
	    {
	      gamma1.l <- gamma1_l_list[l - 2]  # Adjust index due to starting l from 3
			  gamma2.l <- gamma2_l_list[l - 2]
			  tau1.l <- tau1_l_list[l - 2]
			  tau2.l <- tau2_l_list[l - 2]
			  two_power_l <- 2^l 
	      for(j in 1:(two_power_l))
	      {
	        d <- dlj[(two_power_l-1)+j]
	        M1 <- M1_func(r,d,tau1.l,sigmasq)
	        LA <- Lap_approx(d,nu,tau2.l,sigmasq)
	        d.star <- LA$d.max
	        sigma.star <- LA$sigma.d.max

	        aux <- aux + lhood_contrib(d,r,nu,M1,gamma1.l,gamma2.l,tau1.l,tau2.l,sigmasq,d.star,sigma.star)
	      }
	    }
	    -aux
	  },
	  "mom" =
	  function(...,n,dlj,r,nu,sigmasq)
    {
			aux <- 0
			C <- c(...)[order( as.numeric(gsub("\\D","",names(c(...)))) )] #extract digits from names as numeric, order columns by that
			gamma_l_list <- sapply(3:(log2(n)-1), function(l) gamma.l_func(C, l))
			tau_l_list <- sapply(3:(log2(n)-1), function(l) tau.l_func(C, l))
			
			for(l in 3:(log2(n)-1)){	
				gamma.l <- gamma_l_list[l - 2]
  			tau.l <- tau_l_list[l - 2]
  			two_power_l <- 2^l
				for(j in 1:(two_power_l))
				{
					d <- dlj[(two_power_l-1)+j]
					M1 <- M1_func(r,d,tau.l,sigmasq)
					aux1 <- lhood_contrib_indiv(method,d,r,gamma.l,tau.l,sigmasq,M1)
					aux <- aux + aux1
				}
			}
			-aux
    },
    "imom" =
    function(...,n,dlj,r,nu,sigmasq)
    {
      aux <- 0
      C <- c(...)[order( as.numeric(gsub("\\D","",names(c(...)))) )] #extract digits from names as numeric, order columns by that
      gamma_l_list <- sapply(3:(log2(n)-1), function(l) gamma.l_func(C, l))
			tau_l_list <- sapply(3:(log2(n)-1), function(l) tau.l_func(C, l))

      for(l in 3:(log2(n)-1)){
        gamma.l <- gamma_l_list[l - 2]
  			tau.l <- tau_l_list[l - 2]
        two_power_l <- 2^l
				for(j in 1:(two_power_l))
        {
          d <- dlj[(two_power_l-1)+j]
          LA <- Lap_approx(d,r,tau.l,sigmasq)
          d.star <- LA$d.max
          sigma.star <- LA$sigma.d.max

          aux1 <- lhood_contrib_indiv(method,d,r,gamma.l,tau.l,sigmasq,M1=NA,d.star,sigma.star)
          aux <- aux + aux1
        }
      }
      -aux
    }
	)

	#............. Hyperparameter estimation ...........

  npar <- npar_mixprob + npar_scale
  C <- rep(1, npar); names(C) <- paste0("C",1:npar)
  if(method=="mixture" & mixprob_dist=="logit") lower_mixprob <- c(-Inf,1e-50,-Inf,1e-50) else
  if(method=="mixture" & mixprob_dist=="genlogit") lower_mixprob <- c(-Inf,rep(1e-50,2),-Inf,rep(1e-50,2)) else
  if(method=="mixture" & mixprob_dist=="hypsec") lower_mixprob <- c(-Inf,1e-50,-Inf,1e-50) else
  if(method=="mixture" & mixprob_dist=="gennormal") lower_mixprob <- c(-Inf,rep(1e-50,2),-Inf,rep(1e-50,2)) else
  if(method %in% c("mom","imom") & mixprob_dist=="logit") lower_mixprob <- c(-Inf,1e-50) else
  if(method %in% c("mom","imom") & mixprob_dist=="genlogit") lower_mixprob <- c(-Inf,rep(1e-50,2)) else
  if(method %in% c("mom","imom") & mixprob_dist=="hypsec") lower_mixprob <- c(-Inf,1e-50) else
  if(method %in% c("mom","imom") & mixprob_dist=="gennormal") lower_mixprob <- c(-Inf,rep(1e-50,2))  
  lowerb <- c( lower_mixprob, rep(1e-50,npar_scale) )

  for(i in 1:50)
  {
  	ix <- 1:npar_mixprob; ix[c(1,2)] <- ix[c(2,1)]
  	for(j in c((npar_mixprob+1):npar,ix))
  	{
  		max_likelihood <- do.call(nlminb, c(list(start = C[j]),as.list(C[-j]),list(obj=minusll,lower=lowerb[j],upper=Inf,n=n,dlj=dlj,r=r,nu=nu,sigmasq=sigmasq)) )
		  C[j] <- max_likelihood$par		  		
  	}
  }
  hyparEst <- C

  #........ Posterior mean of the wavelet coef .......

  dlj.post.mean <- NULL
  for(l in 0:(log2(n)-1))
  for(j in 1:(2^l))
  {
    if(method=="mixture") aux <- dlj_postmean_mixture_func(n,r,nu,sigmasq,l,j,dlj,hyparEst,gamma1.l_func,gamma2.l_func,tau1.l_func,tau2.l_func) else aux <- dlj_postmean_indiv_func(n,r,nu,sigmasq,l,j,dlj,hyparEst,gamma.l_func,tau.l_func,method)
    dlj.post.mean <- c(dlj.post.mean,aux) 
  }

  #............ Posterior smoothed function ...........
  
  dlj.post.mean1 <- NULL
  for(l in 0:(log2(n)-1))
    dlj.post.mean1 <- c( dlj.post.mean[(2^l):(2^(l+1)-1)], dlj.post.mean1 )
  data.wavelet.post <- wave_obj$obj
  data.wavelet.post$D <- dlj.post.mean1
  data.post.mean <- wr(data.wavelet.post)
  
  if(!is.null(func)) {
  	MSE.mean <- mean((func - data.post.mean)^2)
  	return(list(data=data, func.post.mean=data.post.mean, wavelet.post.mean=dlj.post.mean, wavelet.empirical=dlj, hyperparam=hyparEst, MSE.mean=MSE.mean, sigma=sqrt(sigmasq), runtime=as.numeric(Sys.time()-time,units="secs") ) )
  } else{
  	return(list(data=data, func.post.mean=data.post.mean, wavelet.post.mean=dlj.post.mean, wavelet.empirical=dlj, hyperparam=hyparEst, sigma=sqrt(sigmasq), runtime=as.numeric(Sys.time()-time,units="secs") ) )
  }

}





wavelet_trans <- function(data,filter.number=6,wave.family="DaubLeAsymm",bc="periodic"){
  #...............Wavelet decomposition..................
  n <- length(data)
  data.wavelet <- wd(data, filter.number, wave.family, "wavelet", bc)

  dlj <- NULL
  for(l in 0:(log2(n)-1))
    dlj <- c( dlj , accessD(data.wavelet, level=l) )

  #...............Sigma squared Estimation................
  d_Lj <- accessD(data.wavelet, level=log2(n)-1)
  sigma <- mad(d_Lj, center = median(d_Lj), constant = 1/.6745, low = FALSE)
  sigmasq <- sigma^2

  return(list(dlj=dlj,sigmasq=sigmasq,obj=data.wavelet))
}

# doublefact <- function(x) {
#   if (x%%2 == 0) y <- x-1 else y <- x
#   factorial(y+1)/( 2^((y+1)/2) * factorial((y+1)/2) )
# }

# M1_func <- function(r,d,tau1.l,sigmasq){
#   # M1 <- 0
#   # for(i in 0:r) 
#   #   M1 <- M1 + factorial(2*r)/(factorial(2*i)*factorial(r-i)*2^(r-i)) * ( sqrt(tau1.l/(1+tau1.l)) * (d/sqrt(sigmasq)) )^(2*i)
#   # M1 <- 1/doublefact(2*r-1) * M1
#   # M1
  
#   #faster version
#   sqrt_term <- sqrt(tau1.l / (1 + tau1.l))
#   scaled_d <- d / sqrt(sigmasq)
#   factor <- sqrt_term * scaled_d
  
#   i <- 0:r
  
#   coeffs <- exp(lgamma(2 * r + 1) - lgamma(2 * i + 1) - lgamma(r - i + 1) - (r - i) * log(2))
#   powers <- factor^(2 * i)
  
#   M1 <- sum(coeffs * powers)  
#   M1 <- M1 / doublefact(2 * r - 1)
#   return(M1)
# }

# M2_func <- function(r,d,tau1.l,sigmasq){
#   # M2 <- 0
#   # for(i in 1:(r+1)) 
#   #   M2 <- M2 + factorial(2*r+1)/( factorial(2*i-1) * factorial(r+1-i) * 2^(r+1-i) ) * ( sqrt(tau1.l/(1+tau1.l)) * (d/sqrt(sigmasq)) )^(2*i-1)
#   # M2 <- 1/doublefact(2*r-1) * M2
#   # M2

#   # faster code
#   sqrt_term <- sqrt(tau1.l / (1 + tau1.l))
#   scaled_d <- d / sqrt(sigmasq)
#   factor <- sqrt_term * scaled_d
  
#   i <- 1:(r + 1)
  
#   coeffs <- exp(
#     lgamma(2 * r + 2) - lgamma(2 * i) - lgamma(r + 2 - i) - (r + 1 - i) * log(2)
#   )
#   powers <- factor^(2 * i - 1)
  
#   M2 <- sum(coeffs * powers)  
#   M2 <- M2 / doublefact(2 * r - 1)
#   return(M2)
# }

# h_func <- function(x,d,nu,tau2.l,sigmasq) { 
#   exp( min(709, log( (abs(x))^(-nu-1) ) -(2*sigmasq)^(-1) * (x^2-2*x*d) - tau2.l*sigmasq/x^2 ) ) 
# }

# L_func <- function(x,d,nu,tau2.l,sigmasq) { 
#   -(nu+1)*log(abs(x)) -(2*sigmasq)^(-1) * (x^2-2*x*d) - tau2.l*sigmasq/x^2 
# }

# L.dd_func <- function(x,nu,tau2.l,sigmasq) { 
#   (nu+1)/x^2 - 1/sigmasq - 6*tau2.l*sigmasq/x^4 
# }

# Lap_approx <- function(d,nu,tau2.l,sigmasq){
#   poly <- polynomial(c(-2*tau2.l*sigmasq^2, 0, (nu+1)*sigmasq, -d, 1))
#   d.solve <- Re(solve(poly))
#   L.values <- L_func(d.solve,d,nu,tau2.l,sigmasq)
#   L.dd.values <- L.dd_func(d.solve,nu,tau2.l,sigmasq)
#   id <- which(L.dd.values < 0)
#   id1 <- id[ which(L.values[id]==max(L.values[id])) ]
#   if(length(id1)==1){
#     d.max<- d.solve[id1]
#   } else{
#     if(length(id1)>1){
#       d.max<- d.solve[id1[which.max(L.dd.values[id1])]]
#     } else{
#       d.max<- 1e-50
#     }
#   }  
#   # d.max.set <- d.solve[ which(L.dd.values < 0) ]
#   # if(length(d.max.set)!=0) { 
#   # d.max <- d.max.set[ which.max(L_func(d.max.set,d=d,nu=nu,tau2.l=tau2.l,sigmasq=sigmasq)) ]
#   # L.values.d.max.set <- L_func(d.max.set,d=d,nu=nu,tau2.l=tau2.l,sigmasq=sigmasq)
#   # d.max.cand <- d.max.set[ L.values.d.max.set==max(L.values.d.max.set) ]
#   # d.max <- d.max.cand[ which.min( L.dd_func(d.max.cand,nu,tau2.l,sigmasq) ) ]
#   L.dd.dmax <- L.dd_func(d.max,nu,tau2.l,sigmasq)
#   sigma.d.max <- sqrt(-1/L.dd.dmax)
#   return(list(d.max=d.max,sigma.d.max=sigma.d.max))
# }

# lhood_contrib <- function(d,r,nu,M1,gamma1.l,gamma2.l,tau1.l,tau2.l,sigmasq,d.star,sigma.star){
#   log( exp( max(-745, log(gamma1.l) + (-r-0.5) * log(1+tau1.l) + log(M1) - d^2/(2*sigmasq*(1+tau1.l)) ) )  + exp( max(-745, log(1-gamma1.l) + log(gamma2.l) + nu/2*log(tau2.l*sigmasq) - log(gamma(nu/2)) - d^2/(2*sigmasq) + log(sqrt(2*pi)*sigma.star) + log(h_func(d.star,d,nu,tau2.l,sigmasq)) ) ) + exp( max(-745, log(1-gamma1.l) + log(1-gamma2.l) - d^2/(2*sigmasq) ) ) )
# }

# lhood_contrib_indiv <- function(method,d,r,M1,gamma.l,tau.l,sigmasq,d.star=NULL,sigma.star=NULL){
#   if(method=="mom"){
#   	out <- log( exp( max(-745, log(gamma.l) + (-r-.5)*log(1+tau.l) + log(M1) - d^2/(2*sigmasq*(1+tau.l))) ) + exp( max(-745, log(1-gamma.l) - d^2/(2*sigmasq)) )  )
#   } else{
#   	out <- log( exp( max(-745, log(gamma.l) + r/2*log(tau.l*sigmasq) - log(gamma(r/2)) - d^2/(2*sigmasq) + 
#              log(sqrt(2*pi)*sigma.star) + log(h_func(d.star,d,r,tau.l,sigmasq)) ) ) + exp( max(-745, log(1-gamma.l) - d^2/(2*sigmasq)) ) )
#   }
#   out
# }

# post_odds_func <- function(d,r,nu,M1,sigmasq,gamma1.l,gamma2.l,tau1.l,tau2.l,d.star,sigma.star){
#   if(gamma1.l!=0 & gamma2.l!=1){
#     O1 <-  gamma1.l/((1-gamma1.l) * (1-gamma2.l)) * (1+tau1.l)^(-r-0.5) * M1 * exp(1/(2*sigmasq) * tau1.l/(1+tau1.l) * d^2)
#     # O1_lj <- exp( min(709, log(gamma1.l)-log((1-gamma1.l))-log((1-gamma2.l))+(-r-0.5)*log((1+tau1.l))+log(M1)+1/(2*sigmasq) * tau1.l/(1+tau1.l) * d^2 ))
#   } else{
#     O1 <- 0
#   }
    
#   O2 <- gamma2.l/(1-gamma2.l) * (tau2.l*sigmasq)^(nu/2) * (gamma(nu/2))^(-1) * sqrt(2*pi) * sigma.star * h_func(d.star,d,nu,tau2.l,sigmasq)
#   # O2_lj <- exp( min(709, log(gamma2.l)-log(gamma2.l)+(nu/2)*log(tau2.l*sigma)-gamma(nu/2)+log(2*pi*sigma.star)+log(h(d.star,d,nu,tau2.l,sigmasq)) ))    

#   return(list(O1=O1,O2=O2))
# }

# post_mixprobs_func <- function(post_odds) {
#   O1 <- post_odds$O1
#   O2 <- post_odds$O2
#   if(O1!=Inf & O2!=Inf) 
#   {
#     p1 <- O1/(O1 + O2 + 1)
#     p2 <- O2/(O1 + O2 + 1)
#   } else
#   {
#     if(O1==Inf & O2==Inf)
#     {
#       p1 <- 0.5
#       p2 <- 0.5
#     } else
#     {
#       if(O1==Inf & O2!=Inf)
#       {
#         p1 <- 1
#         p2 <- 0
#       } else
#       {
#         p1 <- 0
#         p2 <- 1
#       }
#     }         
#   } 
#   return(list(p1=p1,p2=p2))
# }


dlj_postmean_mixture_func <- function(n,r,nu,sigmasq,l,j,dlj,hyperparam,gamma1.l_func,gamma2.l_func,tau1.l_func,tau2.l_func)
{
  postmean <- NULL
  d <- dlj[(2^l-1)+j]

  if((l==0)||(l==1)||(l==2)) postmean = d else
  {
    gamma1.l <- gamma1.l_func(hyperparam,l)
    gamma2.l <- gamma2.l_func(hyperparam,l)
    tau1.l <- tau1.l_func(hyperparam,l)
    tau2.l <- tau2.l_func(hyperparam,l)
    M1 <- M1_func(r,d,tau1.l,sigmasq)
    M2 <- M2_func(r,d,tau1.l,sigmasq)
    LA <- Lap_approx(d,nu,tau2.l,sigmasq)
    d.star <- LA$d.max
    sigma.star <- LA$sigma.d.max
   
    post_odds <- post_odds_func(d,r,nu,M1,sigmasq,gamma1.l,gamma2.l,tau1.l,tau2.l,d.star,sigma.star)
    O1_lj <- post_odds$O1
    O2_lj <- post_odds$O2
    post_mixprobs <- post_mixprobs_func(post_odds)
    p1_lj <- post_mixprobs$p1
    p2_lj <- post_mixprobs$p2
    
    aux = p1_lj * sqrt(tau1.l/(1+tau1.l)) * M2/M1 * sqrt(sigmasq) + p2_lj * d.star
    postmean <- c( postmean, aux )
  }
  postmean
}

# post_odds_func_indiv <- function(method,d,r,M1,sigmasq,gamma.l,tau.l,d.star=NULL,sigma.star=NULL){
# 	if(method=="mom"){
# 		O <- gamma.l/(1-gamma.l) * (1 + tau.l)^(-r-.5) * M1 * exp( 1/(2*sigmasq) * tau.l/(1+tau.l) * d^2 )
#   	if(is.nan(O)) O <- Inf
# 	} else{
# 		O <- gamma.l/(1-gamma.l) * (tau.l*sigmasq)^(r/2) * (gamma(r/2))^(-1) * sqrt(2*pi) * sigma.star * h_func(d.star,d,r,tau.l,sigmasq)
# 	}
# 	O
# }

dlj_postmean_indiv_func <- function(n,r,nu,sigmasq,l,j,dlj,hyperparam,gamma.l_func,tau.l_func,method)
{
  postmean <- NULL
  d <- dlj[(2^l-1)+j]

  if((l==0)||(l==1)||(l==2)) postmean = d else
  {
    gamma.l <- gamma.l_func(hyperparam,l)
    tau.l <- tau.l_func(hyperparam,l)
    if(method=="mom")
    {
      M1 <- M1_func(r,d,tau.l,sigmasq)
      M2 <- M2_func(r,d,tau.l,sigmasq)
      O_lj <- post_odds_func_indiv(method,d,r,sigmasq,gamma.l,tau.l,M1)
      # O_lj <- gamma.l/(1-gamma.l) * (1 + tau.l)^(-r-.5) * M1 * exp( 1/(2*sigmasq) * tau.l/(1+tau.l) * d^2 )
      # if(is.nan(O_lj)) O_lj <- Inf
      if(O_lj==Inf) p_lj <- 1  else p_lj <- O_lj/(1+O_lj)            
      aux <- p_lj * sqrt(sigmasq) * sqrt(tau.l/(1+tau.l)) * M2/M1 
    } else
    {
      LA <- Lap_approx(d,nu,tau.l,sigmasq)
      d.star <- LA$d.max
      sigma.star <- LA$sigma.d.max
      O_lj <- post_odds_func_indiv(method,d,r,sigmasq,gamma.l,tau.l,M1=NA,d.star,sigma.star)
      # O_lj <- gamma.l/(1-gamma.l) * (tau.l*sigmasq)^(r/2) * (gamma(r/2))^(-1) * sqrt(2*pi) * sigma.star * h_func(d.star,d,r,tau.l,sigmasq)
      if(O_lj==Inf) p_lj <- 1  else p_lj <- O_lj/(1+O_lj)      
      aux <- p_lj * d.star 
    }
    postmean <- c( postmean, aux )
  }
  postmean
}

























