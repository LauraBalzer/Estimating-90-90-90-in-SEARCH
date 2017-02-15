#=========================================
# R CODE TO EVALUATE THE UNAIDS 90/90/90 CASCADE COVERAGE
# IN THE SEARCH STUDY (INTERIM ANALYSIS - INTEVENTION ARM ONLY)

# Programmer: Laura Balzer
# lbbalzer@hsph.harvard.edu
#
# Requires the SuperLearner v2.0-2 
#	+ ltmle package v0.9-8-4
#
# Last update: February 14, 2017
#
# ******* THIS FILE IS HELPER FUNCTIONS
#==========================================


#===================================================#===================================================
# get.var - function to get inference via the delta method
# 		assumes inputed estimators are asymptotically linear
#		i.e. written in first order as an empircal mean of an influence curve (IC)
#	input:  point estimates (mu1, mu0), corresponding influence curves (IC1, IC0)
#		significance level
#	output: point estimate, var, wald-type CI 
#===================================================#===================================================
get.var<- function(mu1, mu0=NULL, IC1, IC0=NULL, alpha=0.05){
	
	mu1<- unlist(mu1)
		
	if(is.null(mu0)){ 
		# if single TMLE 
		psi<- mu1
		IC<- IC1
		log= F
	
	} else { 
		# if ratio of TMLEs (i.e. target = psi/psi0)
		mu0<- unlist(mu0)
		# get inference via the delta method on log scale
		psi<- log(mu1/mu0)
		IC <- 1/mu1*(IC1) - 1/mu0*IC0
		log=T
	}
	
	# variance of asy lin est is var(IC)/n
	var<- var(IC)/length(IC)
	# testing and CI	
	cutoff <- qnorm(alpha/2, lower.tail=F)
	se<- sqrt(var)
	CI.lo <- psi - cutoff*se
	CI.hi <- psi + cutoff*se

	if(log){
		est<- data.frame(psi=exp(psi), CI.lo=exp(CI.lo), CI.hi=exp(CI.hi) ) 
	}else{
		est<- data.frame(psi=psi, CI.lo=CI.lo, CI.hi=CI.hi) 
	}

	list(est=est, IC=IC)
}


#===================================================#===================================================
# SCREENING ALGORITHMS FOR SUPERLEARNER
# See SuperLearner help file for more info: ?SuperLearner
#===================================================#===================================================
screen.corRank10 <- function(Y, X, family, method = "pearson", rank = 10, ...) {
    listp <- apply(X, 2, function(x, Y, method) {
        ifelse(var(x) <= 0, 1, cor.test(x, y = Y, method = method)$p.value)
    }, Y = Y, method = method)
    whichVariable <- (rank(listp) <= rank)
    return(whichVariable)
}

screen.corRank20<- function(Y, X, family, method = "pearson", rank = 20, ...) {
    listp <- apply(X, 2, function(x, Y, method) {
        ifelse(var(x) <= 0, 1, cor.test(x, y = Y, method = method)$p.value)
    }, Y = Y, method = method)
    whichVariable <- (rank(listp) <= rank)
    return(whichVariable)
}


#===================================================#===================================================
# FUNCTIONS TO ENCODE OUR DETERMINISTIC KNOWLEDGE
# See ltmle help file for more info: ?ltmle
# Also see the Analysis Plan
#===================================================#===================================================

# deterministicQ_YES
# if detQ.variable==1, then outcome==1 with probability 1
deterministicQ_YES<- function(data, current.node, nodes, called.from.estimate.g) {
  L2.index <- which(names(data) == "detQ.variable")
  stopifnot(length(L2.index) == 1)
  L2.in.history <- L2.index < current.node
  if (! L2.in.history) return(NULL)
  
  is.deterministic <- data[,L2.index]==1
  return(list(is.deterministic=is.deterministic, Q.value=1))
}

# deterministicQ_NO
# if detQ.variable==0,  then outcome==0 with probability 1
deterministicQ_NO<- function(data, current.node, nodes, called.from.estimate.g) {
  L2.index <- which(names(data) == "detQ.variable")
  stopifnot(length(L2.index) == 1)
  L2.in.history <- L2.index < current.node
  if (! L2.in.history) return(NULL)
  
  is.deterministic <- data[,L2.index]== 0
  return(list(is.deterministic=is.deterministic, Q.value=0))
}

# deterministicQ_combo 
# cannot be suppressed if dead, outmigrated or not on ART
# cannot have Z*=1 if combo= (D=1 OR M=1 OR eART=0)
deterministicQ_combo<- function(data, current.node, nodes, called.from.estimate.g) {
  L2.index <- which(names(data) == "combo")
  stopifnot(length(L2.index) == 1)
  L2.in.history <- L2.index < current.node
  if (! L2.in.history) return(NULL)
  
  is.deterministic <- data[,L2.index]==1
  return(list(is.deterministic=is.deterministic, Q.value=0))
}


#===================================================#===================================================
#  FUNCTIONS TO MAKE OUTPUT PRETTY
#===================================================#===================================================

#-------------------------------------#---------------------------------------------------
# make.graph1: code to create Figure1
#		Runs on output from do.serial.analysis
#		See Cascade_RunSerialFig1
#-------------------------------------#---------------------------------------------------

make.graph1 <- function(out.t0, out.t1, out.t2, this.strata){
	
	Cascade<- NULL
	temp <- out.t0$Cascade
	colnames(temp)<- c('pt', 'CI.lo','CI.hi')
	rownames(temp)<- paste(this.strata, rownames(temp), sep=': ')
	rownames(temp)<- paste(rownames(temp), ' t', 0, sep='' )
	Cascade<- rbind(Cascade, round(temp[2:5, ],3), blank0=rep('',3))
	
	temp <- out.t1$Cascade
	colnames(temp)<- c('pt', 'CI.lo','CI.hi')
	rownames(temp)<- paste(this.strata, rownames(temp), sep=': ')
	rownames(temp)<- paste(rownames(temp), ' t', 1, sep='' )
	Cascade<- rbind(Cascade, round(temp[2:5, ],3), blank1=rep('',3))

	temp <- out.t2$Cascade
	colnames(temp)<- c('pt', 'CI.lo','CI.hi')
	rownames(temp)<- paste(this.strata, rownames(temp), sep=': ')
	rownames(temp)<- paste(rownames(temp), ' t', 2, sep='' )
	Cascade<- rbind(Cascade, round(temp[2:5, ],3) )

	
	#numerators/denominators
	temp0<- c(paste(out.t0$RAW['N.pDx',], out.t0$RAW['N.est.HIV.pos',], sep='/'),
					 paste(out.t0$RAW['N.eART',], out.t0$RAW['N.pDx',],  sep='/'),
					 paste( out.t0$RAW['N.est.suppression',], out.t0$RAW['N.eART',], sep='/'),
						 paste( out.t0$RAW['N.est.suppression',], out.t0$RAW['N.est.HIV.pos',], sep='/'))
	
	temp1<- c(paste(out.t1$RAW['N.pDx',], out.t1$RAW['N.est.HIV.pos',], sep='/'),
					 paste(out.t1$RAW['N.eART',], out.t1$RAW['N.pDx',],  sep='/'),
					 paste( out.t1$RAW['N.est.suppression',], out.t1$RAW['N.eART',], sep='/'),
					 paste( out.t1$RAW['N.est.suppression',], out.t1$RAW['N.est.HIV.pos',], sep='/'))

	temp2<- c(paste(out.t2$RAW['N.pDx',], out.t2$RAW['N.est.HIV.pos',], sep='/'),
					 paste(out.t2$RAW['N.eART',], out.t2$RAW['N.pDx',],  sep='/'),
					 paste( out.t2$RAW['N.est.suppression',], out.t2$RAW['N.eART',], sep='/'),
					 paste( out.t2$RAW['N.est.suppression',], out.t2$RAW['N.est.HIV.pos',], sep='/'))

	
	out = cbind(Cascade, Denominators=c(temp0, '', temp1, '', temp2))
	colnames(out)<- c(colnames(Cascade), 'Num/Den')
	out
}

#-------------------------------------#---------------------------------------------------
# make.eTable1: code to create eTable1
#		Runs on output from do.serial.analysis
#		See Cascade_RunSerialFig1
#-------------------------------------#---------------------------------------------------

make.eTable1 <- function(raw, Cascade.secondary, Cascade){
	
	CascadeEst<- Step1<- Step2<- Step3<- Step4 <- rep('', 3)
	OUT<- NULL
	OUT<- rbind(OUT, raw, CascadeEst=CascadeEst, Step1=Step1 )
	OUT<- rbind(OUT, Step1.secondary= get.fraction.inference( num=raw['N.Delta.pDx',] , den=raw['N.Delta.HIV.pos',],  
		inference=Cascade.secondary[ 'P(pDx=1 | Y*=1)', ]) )
	OUT<- rbind(OUT,  Step1.primary=get.fraction.inference( num=raw['N.pDx',], den=raw['N.est.HIV.pos',],  
		inference=Cascade[ 'P(pDx=1 | Y*=1)', ]) )	
	
	OUT<- rbind(OUT, Step2=Step2)
	OUT<- rbind(OUT, Step2.secondary=get.fraction.inference( num=raw['N.Delta.eART',], den=raw['N.Delta.pDx',],  
		inference=Cascade.secondary[ 'P(eART=1 | pDx=1)', ]))
	OUT<- rbind(OUT, Step2.primary=get.fraction.inference( num=raw['N.eART',], den=raw['N.pDx',],  
		inference=Cascade[ 'P(eART=1 | pDx=1)', ]))
		
	OUT<- rbind(OUT, Step3=Step3)
	OUT<- rbind(OUT, Step3.secondary= get.fraction.inference( num=raw['N.eART.and.obs.suppression',], 
		den=raw['N.eART.and.measuredVL',],  inference=Cascade.secondary['P(Supp*=1 | eART=1)', ]))
	OUT<- rbind(OUT,  Step3.primary=get.fraction.inference( num=raw['N.est.suppression',], den=raw['N.eART',],  
		inference=Cascade['P(Supp*=1 | eART=1)', ]))	
	
	OUT<- rbind(OUT, Step4=Step4)
	OUT<- rbind(OUT, Step4.secondary= get.fraction.inference( num=raw['N.eART.and.obs.suppression',], 
		den=raw['N.measured.VL',],  inference=Cascade.secondary[ 'P(Supp*=1 | Y*=1)', ]) )
	OUT<- rbind(OUT, Step4.primary=get.fraction.inference( num=raw['N.est.suppression',], 
		den=raw['N.est.HIV.pos',],  inference=Cascade[ 'P(Supp*=1 | Y*=1)', ]) )	
		
	OUT
}

#-------------------------------------#---------------------------------------------------
# Helper function for make.eTable1
#	uses paste function to get proper output
#-------------------------------------#---------------------------------------------------

get.fraction.inference <- function( num, den, inference){
	
	out<- rep(NA, 3)
	for(k in 1:3){
		ratio <- paste(num[k], den[k], sep='/')
		pt <-  inference[, grep(paste('pt',(k-1), sep=''), colnames(inference)) ]
		pt<- paste( ' ', round(pt*100,1), '% ', sep='' )
		CI.lo <- inference[, grep(paste('CI.lo',(k-1), sep=''), colnames(inference)) ] 
		CI.hi <- inference[, grep(paste('CI.hi',(k-1), sep=''), colnames(inference)) ] 

		CI<- paste( '(', round(CI.lo*100,1), '%, ', round(CI.hi*100,1), '%)', sep='')
		
		out[k] <- paste(ratio, pt, CI, sep='')
	}
	out	
}


#-------------------------------------#---------------------------------------------------
# make.graph2: code to create Figure2
#		Runs on output from do.long.no.int.censor
#		See Cascade_RunLongFig2
#-------------------------------------#---------------------------------------------------

make.graph2 <- function(out.t0, out.t1, out.t2, this.strata){
	
	these<- c(1:4, 8,9)
	OUT<- data.frame(c(out.t0$full[ these ,'pt0'], out.t1$full[these ,'pt1'], out.t2$full[these,'pt2']))
	rownames(OUT) <- c( paste(rownames(out.t0$full)[these], ' t', 0, sep='' ),
										paste(rownames(out.t1$full)[these], ' t', 1, sep='' ), 
										paste(rownames(out.t2$full)[these], ' t', 2, sep='' ))
	rownames(OUT)<- paste(this.strata, rownames(OUT), sep=': ')
	colnames(OUT)<- 'pt'
	OUT
}

#-------------------------------------#---------------------------------------------------
# make.eTable2: code to create eTable4
#		Runs on output from do.long.no.int.censor
#		See Cascade_RunLongFig2
#-------------------------------------#---------------------------------------------------
make.eTable2 <- function(raw, full){
	
	OUT<- NULL
	OUT<- rbind(OUT, Death=get.fraction.inference(num=raw['N.death',], den=raw['N.population',], 
		inference=full['p.died',] ))	
	OUT<- rbind(OUT, Outmigrate=get.fraction.inference(num=raw['N.migrate',], den=raw['N.population',], 
		inference=full['p.migrate',] ))	
	OUT<- rbind(OUT, New.Dx=get.fraction.inference(num=raw['N.new.dx',], den=raw['N.population',], 
		inference=full['p.new.dx',] ))	
	OUT<- rbind(OUT, No.ART=get.fraction.inference(num=raw['N.no.ART',], den=raw['N.population',], 
		inference=full['p.no.ART',] ))	
	
	OUT<- rbind(OUT, Supp.secondary=rep('',3))
	OUT<- rbind(OUT, VL.missed=get.fraction.inference(num=raw['N.VL.missed',], den=raw['N.population',], 
		inference=full['p.VL.missed',] ))	
	OUT<- rbind(OUT, Supp.failure.second=get.fraction.inference(num=raw['N.supp.failure.secondary',], den=raw['N.population',], 
		inference=full['p.supp.fail.secondary',] ))	
	OUT<- rbind(OUT, Supp.succ.second=get.fraction.inference(num=raw['N.supp.succ.secondary',], den=raw['N.population',], 
		inference=full['p.supp.succ.secondary',] ))	
	
	OUT<- rbind(OUT, Supp.primary=rep('',3))
	OUT<- rbind(OUT, Supp.failure.primary=get.fraction.inference(num=raw['N.est.supp.failure',], den=raw['N.population',], 
		inference=full['p.supp.fail',] ))	
	OUT<- rbind(OUT, Supp.succ.primary=get.fraction.inference(num=raw['N.est.supp.succ',], den=raw['N.population',], 
		inference=full['p.supp.succ',] ))	
	colnames(OUT)<- c('t0', 't1', 't2')
	OUT
}

