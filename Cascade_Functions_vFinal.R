#==========================================
# R CODE TO EVALUATE THE UNAIDS 90/90/90 CASCADE COVERAGE
# IN THE SEARCH STUDY (INTERIM ANALYSIS - INTEVENTION ARM ONLY)
#
# Programmer: Laura Balzer
# lbbalzer@hsph.harvard.edu
#
# Requires the SuperLearner v2.0-2 
#	+ ltmle package v0.9-8-4
#
# Last update: February 14, 2017
#
# MAIN FUNCTIONS FOR PRIMARY + SECONDARY ANALYSIS FOR
# 		FIGURES 1 AND 2
#==========================================================


#==========================================================#==========================================================
# do.serial.analysis - code to run primary and secondary analyses for 
#	serial cross-sectional analysis of prevalent HIV+ (Figure1)
#
#	input: data, baseline predictors, time t, indicator of inclusion in the analysis (these.obs),
#		variable names (prior Dx, ever ART use, VL testing, suppression)
#		SuperLearner library+ CV-scheme (SL.library, CV-size), 
#		unit of independence (id)
#
# output: raw output of observed and estimated numerators and denominators
#	primary and secondary analyses (pt estimates with 95% CI)
#==========================================================#==========================================================

do.serial.analysis<- function(data, baseline.pred, t, these.obs,
	pDx.variable, eART.variable, tstvl.variable, supp.variable,
	SL.library=NULL, CV.size=NULL, id=NULL){
	
	#--------------------------------------------------------------------------------------------------------
	# Define variables of interest at the relevant time
	pDx <- data[, paste(pDx.variable, t, sep='_')]
	eART <-  data[, paste(eART.variable, t, sep='_')]
	Delta <- data[, paste('delta', t, sep='_')]
	HIVpos <- data[, paste('y', t, sep='_')]
	TstVL <-  data[, paste(tstvl.variable, t, sep='_')]
	Supp <- data[, paste(supp.variable, t, sep='_')]
	
	#--------------------------------------------------------------------------------------------------------
	# Setup data frames for Cascade estimates
	Cascade<- data.frame(matrix(NA, nrow=5, ncol=3))
	rownames(Cascade)<- c('P(pDx=1 | Y*=1)', 'P(eART=1 | pDx=1)', 
		'P(Supp*=1 | eART=1)', 'P(Supp*=1 | Y*=1)', 'Product')
	colnames(Cascade) <- c('pt', 'CI.lo', 'CI.hi')
	Cascade.secondary<- Cascade
	
	#--------------------------------------------------------------------------------------------------------
	# Number in the population of interest
	N.population <-  sum(these.obs)
	
	#--------------------------------------------------------------------------------------------------------
	# Number with known HIV status & seen CHC/track = Number in population with Delta=1
	N.Delta <- sum(Delta & these.obs )
			
	#==========================================================
	# Prevalence of HIV at time t: P(Y*=1)  
	#==========================================================
	
	#************  Primary analysis  ****************  
	#outcome as observed HIV status
	Y<- HIVpos
	# intervention variable: seen at CHC/tracking with HIV status known
	A<- Delta
	
	# specify adjustment set other than baseline covariates
	# there is no variability in the rest of the variables given Y_{t-1}==0
	L <- c('chc', 'tr', 'tsthiv')
	adj <- get.adjustment.HIVpos(data, baseline.pred, L=L, t)
	
	# call ltmle using deterministic knowledge about known HIV status at prior time point (HIVpos_{t-1})
	Prob.HIVpos <- call.ltmle(data=data, t=t,  pred.A=adj$adj, A=A, Y= Y,
		id=id, these=these.obs, SL.library=SL.library, CV.size=CV.size, 
		deterministicQ= adj$detQ)		
	
	# # Side note: This is equivalent to the following 
	# # *** for t>0 only!! 
	# if(t>0){
		# pastHIV <- adj$adj$detQ.variable
		# temp <-  call.ltmle(data=data, t=t,  pred.A=adj$adj, A=A, Y= Y,
			# id=id , these=these.obs & !pastHIV, SL.library=SL.library, CV.size=CV.size)
			# c( Prob.HIVpos$est$pt, temp$e$pt*mean(!pastHIV[these.obs]) + mean(pastHIV[these.obs]))
	# }
	
	# Now obtain an adjusted estimate of the number HIV+ 
	# Number estimated to be HIV+ = (Estimated HIV prevalence) x (Population size) 
	N.est.HIV.pos<- round(Prob.HIVpos$e$pt*N.population, 0)


	#*************** Secondary analysis  ******************************
	# Number observed HIV+  = Number in population with Y=1  & Delta=1
	N.Delta.HIV.pos <- sum( Y &  A & these.obs )
	
	#Could also calculate using intervention=delta without adjustment
	# Latter is used to obtain inference
	Prob.HIVpos.secondary <- call.ltmle(data=data, t=t,  pred.A=NULL, A=A, Y= Y,
		id=id , these=these.obs, SL.library=NULL, CV.size=CV.size)
	c(Prob.HIVpos.secondary$est$pt, N.Delta.HIV.pos/N.Delta)
	
	rm(A, Y)


	#==========================================================	
	# Previous diagnosis at time t: P(pDx=1, Y*=1) 
	#==========================================================
	
	#************  Primary analysis  ****************  
	# Number with prior diagnosis in the population
	N.pDx<- sum( pDx & these.obs)
	
	#Could also calculate using a dummy intervention variable
	# Latter is used to obtain inference
	Prob.pDx <- call.ltmle(data=data, t=t,  pred.A= NULL, A=rep(1, nrow(data)), Y= pDx,
		id=id , these=these.obs, SL.library=NULL, CV.size=CV.size)
	 c(Prob.pDx$est$pt, N.pDx/N.population)
	
	
	#*************** Secondary analysis  ******************************
	# Number with prior diagnosis seen at CHC/track
	N.Delta.pDx <- sum( pDx &  Delta & these.obs )
	
	#Could also calculate using intervention=delta without adjustment
	# Latter is used to obtain inference
	Prob.pDx.secondary <- call.ltmle(data=data, t=t,  pred.A=NULL, A=Delta, Y= pDx,
		id=id , these=these.obs, SL.library=NULL, CV.size=CV.size)
	c(Prob.pDx.secondary$est$pt, N.Delta.pDx/N.Delta)
	
	
	#==========================================================		
	# Ever ART use at time t: P(eART=1, pDx=1, Y*=1) 
	#==========================================================
	
	#************  Primary analysis  ****************  
	# Number with eART in the population
	N.eART<- sum( eART & these.obs)
	
	#Could also calculate using a dummy intervention variable
	# Latter is used to obtain inference
	Prob.eART <- call.ltmle(data=data, t=t,  pred.A= NULL, A=rep(1, nrow(data)), Y= eART,
		id=id , these=these.obs, SL.library=NULL, CV.size=CV.size)
	c(Prob.eART$est$pt, N.eART/N.population)
	

	#*************** Secondary analysis  ******************************
	# Number with prior diagnosis seen at CHC/track
	N.Delta.eART <- sum( eART &  Delta & these.obs )
	
	#Could also calculate using intervention=delta without adjustment
	# Latter is used to obtain inference
	Prob.eART.secondary <- call.ltmle(data=data, t=t,  pred.A=NULL, A=Delta, Y= eART,
		id=id , these=these.obs, SL.library=NULL, CV.size=CV.size)
	c(Prob.eART.secondary$est$pt, N.Delta.eART/N.Delta)

	
	#==========================================================	
	# Suppression at t:  P(Supp*=1, eART=1, pDx=1, Y*=1) 
	#==========================================================	
	
	#****************Primary analysis   ******************************	
	# specify the outcome as observed suppression
	Y<- Supp
	
	# specify the intervention variable as having a measured VL
	# (TstVL=1 only if Delta=1)
	A<-  TstVL

	# specify history other than baseline covariates
	L<- c(pDx.variable, eART.variable, 'chc', 'tr', 'tsthiv', 'y', tstvl.variable, supp.variable, 'eversuppressed' )
	# no need to include (ever.test, delta) as these are a linear combination of above
	
	# get adjustment set
	adj <- get.adjustment.eART(data, baseline.pred,  L=L, t=t)
	
	Prob.supp.eART.pDx.HIVpos <- call.ltmle(data=data, t=t,  pred.A= adj, A=A, Y= Y, 
		id=id , these=these.obs, SL.library=SL.library, CV.size=CV.size, 
		 deterministicQ= deterministicQ_NO) 
	
	# #  Side note: This is equivalent to the following  
	# temp <-  call.ltmle(data=data, t=t,  pred.A= adj, A=A, 
		 # Y= Y, id=id , these=these.obs & eART, SL.library=SL.library, CV.size=CV.size)
	# c(Prob.supp.eART.pDx.HIVpos$est$pt, temp$e$pt*mean(eART[these.obs]))	
	
	# Now obtain an adjusted estimate of the number virally suppressed
	# Number estimated to be suppressed = (Estimated Suppression ) x (Population size)
	N.est.suppression <- round(Prob.supp.eART.pDx.HIVpos$e$pt*N.population, 0)

	rm(A,Y)
	
	
	#---------------------------------#---------------------------------	#----------------------------------------------------
	# USE THE ABOVE FOUR PROBABILITIES TO CALCULATE CASCADE COVERAGE 
	# 	USE THE DELTA.METHOD FOR INFERENCE IN PRIMARY ANALYSIS- See get.var() function
	#---------------------------------#---------------------------------	#----------------------------------------------------

	#==========================================================
	# Cascade step1: Proportion of HIV+ who prev dx:   P(pDx=1 | Y*=1)
	#==========================================================
	# ***** Primary analysis
	Cascade['P(pDx=1 | Y*=1)', ] <-  	get.var(mu1= Prob.pDx$est$pt, IC1= Prob.pDx$IC, 
																	mu0= Prob.HIVpos$est$pt, IC0=Prob.HIVpos$IC)$est
	
	#******* Secondary analysis 
	Cascade.secondary['P(pDx=1 | Y*=1)', ] <-  get.var(mu1= Prob.pDx.secondary$est$pt, IC1= Prob.pDx.secondary$IC, 
																		mu0= Prob.HIVpos.secondary$est$pt, IC0= Prob.HIVpos.secondary$IC)$est
	
	#==========================================================	
	# Cascade step2: Proportion of prev dx HIV+ who have ever been on ART:     P(eART=1 | pDx=1)
	#==========================================================
	# ***** Primary analysis
	Cascade['P(eART=1 | pDx=1)',] <- 	get.var(mu1= Prob.eART$est$pt, IC1= Prob.eART$IC, 
																		mu0= Prob.pDx$est$pt, IC0= Prob.pDx$IC)$est
	
	#******* Secondary analysis 
	Cascade.secondary['P(eART=1 | pDx=1)',] <- get.var(mu1= Prob.eART.secondary$est$pt, IC1= Prob.eART.secondary$IC, 
																				 mu0= Prob.pDx.secondary$est$pt, IC0= Prob.pDx.secondary$IC)$est
	
	
	#==========================================================	
	# Cascade step3: proportion of previously txt HIV+ who are suppressed:    P(Supp*=1 | eART=1)
	#==========================================================	
	
	#******* Primary analysis *****************
	Cascade['P(Supp*=1 | eART=1)', ] <- get.var(mu1= Prob.supp.eART.pDx.HIVpos$est$pt, IC1= Prob.supp.eART.pDx.HIVpos$IC,
					 mu0= Prob.eART$est$pt, IC0= Prob.eART$IC)$est
	
	
	#***********  Secondary analysis  **********************
	# Suppression among those with eART use 
	# Number with eART and observed suppression / Number with eART and measured viral load 
	
	# Number with eART and measured viral load
	N.eART.and.measuredVL <- sum( eART & TstVL &  these.obs)
	# Number with eART and observed suppression
	N.eART.and.obs.suppression <- sum(Supp & eART & TstVL & these.obs)
		
	#Could also calculate directly calculate suppression among eART=1
	# Latter is used to obtain inference
	Y<- Supp
	Cascade.secondary['P(Supp*=1 | eART=1)', ] <- call.ltmle(data=data, t=t,  pred.A= NULL, A= TstVL, Y= Y, 
		id=id , these=these.obs & eART, SL.library=NULL, CV.size=CV.size)$est 
	c( Cascade.secondary['P(Supp*=1 | eART=1)', 'pt'], N.eART.and.obs.suppression/N.eART.and.measuredVL)
	
	rm(Y)
	
	#==========================================================	
	# Cascade step4: proportion of  HIV+ who are suppressed: P(Supp*=1 | Y*=1)
	#==========================================================	
	#******* Primary analysis *****************
	Cascade['P(Supp*=1 | Y*=1)', ]<- get.var(mu1= Prob.supp.eART.pDx.HIVpos$est$pt, IC1= Prob.supp.eART.pDx.HIVpos$IC, 
										mu0= Prob.HIVpos$est$pt, IC0= Prob.HIVpos$IC)$est 
	
	#***********  Secondary analysis  **********************
	# Probability of suppression among all HIV+ 
	# Number with observed suppression / Number with measured viral load 
	
	# Number with measured viral load
	N.measured.VL <- sum(  TstVL &  these.obs)
	# Number  observed suppression
	N.obs.suppression <- sum(Supp & TstVL & these.obs)
		
	#Could also calculate directly calculate suppression among HIV+
	# Latter is used to obtain inference
	Cascade.secondary['P(Supp*=1 | Y*=1)', ] <- call.ltmle(data=data, t=t,  pred.A= NULL, A=TstVL,  Y= Supp, 
		id=id , these=these.obs, SL.library=NULL, CV.size=CV.size)$est
	c(Cascade.secondary['P(Supp*=1 | Y*=1)', 'pt'], N.obs.suppression/N.measured.VL)
	
		
	#==========================================================	
	# Cascade step4: via the product
	#==========================================================

	Cascade['Product', 1] <- Cascade['P(pDx=1 | Y*=1)',1]* Cascade['P(eART=1 | pDx=1)',1]*Cascade['P(Supp*=1 | eART=1)',1]
	Cascade.secondary['Product', 1] <- Cascade.secondary['P(pDx=1 | Y*=1)',1]*Cascade.secondary['P(eART=1 | pDx=1)',1]*Cascade.secondary['P(Supp*=1 | eART=1)',1]


	#==========================================================
	# Compiling the results
	#==========================================================	 
	
	# create a data frame for the raw & estimated numerators & denominators (eTable1)
	# N.eART.and.obs.suppression==N.obs.suppression; so drop latter
	RAW<- data.frame(rbind(N.population=N.population, N.Delta, N.Delta.HIV.pos, N.est.HIV.pos,
						N.Delta.pDx, N.pDx, N.Delta.eART, N.eART, N.eART.and.measuredVL, N.eART.and.obs.suppression, 
						N.measured.VL, N.est.suppression))
	colnames(RAW)<- paste('t',t, sep='')
	
	# Output HIV prevalence & Cascade Coverage
	Cascade <- rbind(HIVprev = Prob.HIVpos$e, Cascade)
	Cascade.secondary <-	rbind(HIVprev = Prob.HIVpos.secondary$e, Cascade.secondary)	
	
	colnames(Cascade)<- paste(colnames(Cascade), t,sep='')			
	colnames(Cascade.secondary)<- paste(colnames(Cascade.secondary), t,sep='')			
			
	list(RAW=RAW, Cascade=Cascade, Cascade.secondary=Cascade.secondary)
	
}

#===================================================#===================================================
#===================================================#===================================================
#===================================================#===================================================


#-------------------------------------------
# get.adjustment.HIVpos: function to get the adjustment set for the serial X analyses
#	when estimating prevalence
#--------------------------------------------

get.adjustment.HIVpos<- function(data, baseline.pred, L, t){
	
	if(t==0){
		# no deterministic knowledge at BL
		adj <- baseline.pred
		detQ <- NULL
	
	} else if(t==1){
		# if t=1, then deterministic knowledge that outcome (HIV+ or pDX)==1, if HIV0==1

		adj <- data.frame( baseline.pred, data[,paste(L, 0, sep='_')], data[,paste('y', 0, sep='_')] )
		# rename the final column as "detQ.variable"
		# this is the variable that we have knowledge about in deterministic Q
		colnames(adj)[ncol(adj)]<-  'detQ.variable'
		detQ <- deterministicQ_YES
		
	} else if(t==2) {
		# if t=2, then deterministic knowledge that outcome (HIV+ or pDX)==1, if HIV1==1

		adj <- data.frame( baseline.pred, data[,paste(L, 0, sep='_')],  data[,paste(L, 1, sep='_')], 
			data[,paste('y', 1, sep='_')] )
		# rename the final column as "detQ.variable"
		# this is the variable that we have knowledge about in deterministic Q
		colnames(adj)[ncol(adj)]<-  'detQ.variable'
		detQ <- deterministicQ_YES

	} 
	list(adj=adj, detQ=detQ)
}



#-------------------------------------------
# get.adjustment.eART: function to get the adjustment set for the serial X analyses
#		when estimating suppression
#--------------------------------------------
get.adjustment.eART <- function(data, baseline.pred, L, t){
	
	if(t==0){
		adj <- data.frame( baseline.pred, data[,paste(eART.variable, 0, sep='_')] )
	}else if(t==1){
		adj <- data.frame(baseline.pred, data[,paste(L, 0, sep='_')],   data[,paste(eART.variable, 1, sep='_')] )
	} else if (t==2){
		adj <- data.frame(baseline.pred, data[,paste(L, 0, sep='_')], data[,paste(L, 1, sep='_')],
						data[,paste(eART.variable, 2, sep='_')]  )
	} 
	# rename the final column as "detQ.variable"
	# this is the variable that we have knowledge about in deterministic Q
	colnames(adj)[ncol(adj)]<-  'detQ.variable'

	adj
}


#===================================================#===================================================
#===================================================#===================================================
#===================================================#===================================================


#-----------------------------------------------------------
# call.ltmle: function to call the ltmle package v0.9-8-4
#
# input: data, time t, adjustment set (pred.A), intervention variable (A), outcome Y, 
#		indicator of the indpt unit (id), indicator of inclusion in the analysis (these),
#		SuperLearner setup (SL library, CV.size),
#		whether ltmle should estimate time unitl completion (EST.TIME),
#		deterministicQ function
#
# output: list with 
#		estimates (pt, 95% confidence intervals)
#		influence curve at the household level (IC)
#--------------------------------------------------------------------

call.ltmle<- function(data, t, pred.A=NULL, A,  Y, 
	id=NULL, these, SL.library=NULL, CV.size=NULL,  
	EST.TIME=F, deterministicQ=NULL){
	

	# independent unit
	if(!is.null(id)){
		id<- id[these]
	}
	
	# create temporary data frame that is appropriately subsetted
	if(is.null(pred.A)){
		data.temp<- data.frame(A, Y)[these, ]
	} else{
		data.temp<- data.frame(pred.A, A, Y)[these, ] 
	}
				
	# Setup SuperLearner			
	if(!is.null(SL.library)){  
		SL.options<- list(SL.library=SL.library, verbose=T, cvControl=list(V= CV.size, stratifyCV=F))
	}else{
		SL.options<- NULL
	}
		
	est.temp<- ltmle(data=data.temp, Anodes='A', Lnodes=NULL, Ynodes='Y', abar=1,
		stratify = T, SL.options=SL.options, id=id,
		variance.method='ic', estimate.time=EST.TIME,
		deterministic.Q.function= deterministicQ)
	#print(est.temp$fit$g)
	#print(est.temp$fit$Q)
					
	IC<- est.temp$IC$tmle

	est<- data.frame(pt=est.temp$estimate["tmle"], 
		CI.lo=summary(est.temp)$treatment$CI[1], 
		CI.hi=summary(est.temp)$treatment$CI[2] 	)
	
	list(est=est,IC=IC)
}



#===================================================#===================================================
#===================================================#===================================================
#===================================================#===================================================


#==========================================================#==========================================================

# do.long.no.int.censor: code to do longitudinal analyses (primary + secondary) among known BL HIV+ 
#	 with no intervention on censoring (Figure2)
#
# input: data, time t, predictors (pred.A), population indicators (these.obs), 
#	variable names for prior Dx, ever ART use, VL measurement & Suppression, 
#	SL library, CV-scheme, id (indpt unit)
#
# output: raw output of observed and estimated numerators and denominators
#	primary and secondary analyses (pt estimates with 95% CI)
#	
#==========================================================#==========================================================


do.long.no.int.censor <- function(data, t, pred.A, these.obs,
	pDx.variable, eART.variable, tstvl.variable, supp.variable,
	SL.library, CV.size, id){


	# variable definition 
	death <- as.numeric(data[, paste('dead_by', t, sep='_')] )
	migrate <- as.numeric( data[, paste('move_by', t, sep='_')] )
	pDx <- as.numeric( data[, paste(pDx.variable, t, sep='_')] )
	eART <- as.numeric( data[, paste(eART.variable, t, sep='_')] )
	TstVL <- as.numeric( data[, paste(tstvl.variable, t, sep='_')] )
	supp <-  as.numeric( data[, paste(supp.variable, t, sep='_')] )
		
	# throughout we condition on being BL HIV+
	use.these<- these.obs	& data$y_0==1
	N.population <- sum(use.these)
	
	#------------------------------------------------------------------		
	# *********death*********** 
	Y= death
	tmle.died<- call.ltmle(data=data, t=t, pred.A=NULL, A=rep(1,nrow(data)), Y= Y, 
			id=id, these= use.these, SL.library=NULL)
	p.died<-  tmle.died$e
	
	N.death <- sum(Y & use.these)
	c(p.died$pt, N.death/N.population)
	rm(Y)
	
	#------------------------------------------------------------------		
	# *******out migrate & not dead****
	
	Y=as.numeric(migrate==1 & death==0)
	tmle.migrate <- call.ltmle(data=data, t=t, pred.A=NULL, A=rep(1,nrow(data)), Y=Y,
				id=id, these=use.these, SL.library=NULL)
	p.migrate <-   tmle.migrate$e
	
	N.migrate<- sum(Y & use.these)
	c(p.migrate$pt, N.migrate/N.population)
	rm(Y)
	
	#------------------------------------------------------------------		
	# ****new diagnosis & not outmigrated & not dead****
	
	Y=as.numeric(pDx ==0 & migrate==0 & death==0)
	tmle.new.dx <- call.ltmle(data=data, t=t, pred.A=NULL, A=rep(1,nrow(data)), Y=Y, 
				id=id, these=use.these, SL.library=NULL)
	p.new.dx <- tmle.new.dx$e
	
	N.new.dx<- sum(Y & use.these)
	c(p.new.dx$pt, N.new.dx/N.population)
	rm(Y)
	
	#------------------------------------------------------------------		
	#****prior diagnosis but no ART (& not outmigrated & not dead)****
	Y=as.numeric(eART==0 & pDx ==1 & migrate==0 & death==0)
	tmle.no.ART <- call.ltmle(data=data, t=t, pred.A=NULL, A=rep(1,nrow(data)), Y=Y, 
				id=id, these=use.these , SL.library=NULL)
	p.no.ART <-  tmle.no.ART$e
	
	N.no.ART <- sum(Y & use.these) 	
	c(p.no.ART$pt, N.no.ART/N.population)
	rm(Y)

	#---------------------------------------------------------------
	# ****** on ART but missed Viral Load measure *** Secondary Analysis **** 
	Y=as.numeric(TstVL==0 & eART==1 & pDx ==1 & migrate==0 & death==0)
	tmle.VL.missed <- call.ltmle(data=data, t=t, pred.A=NULL, A=rep(1,nrow(data)), Y=Y, 
				id=id, these=use.these , SL.library=NULL)
	p.VL.missed <-  tmle.VL.missed$e

	N.VL.missed <- sum(Y & use.these)  	
	c(p.VL.missed$pt, N.VL.missed/N.population)
	rm(Y)

	#---------------------------------------------------------------
	#****** Measured suppression failure *** Secondary Analysis **** 
	Y<- as.numeric( (supp==0) & (TstVL==1) & (eART==1) & (pDx==1) & (migrate==0) & (death==0) )	
	tmle.supp.failure.secondary <- call.ltmle(data=data, t=t, pred.A=NULL, A=rep(1,nrow(data)), Y=Y, 
				id=id, these=use.these , SL.library=NULL)
	p.supp.failure.secondary<-  tmle.supp.failure.secondary$e

	N.supp.failure.secondary <-  sum(Y & use.these)  	
	c(p.supp.failure.secondary$pt, N.supp.failure.secondary/N.population)
	rm(Y)
	
	#---------------------------------------------------------------
	# ****** Measured suppression success  *** Secondary Only****
	Y<- as.numeric( (supp==1) & (TstVL==1) & (eART==1) & (pDx==1) & (migrate==0) & (death==0) )	
	tmle.supp.succ.secondary <- call.ltmle(data=data, t=t, pred.A=NULL, A=rep(1,nrow(data)), Y=Y, 
				id=id, these=use.these , SL.library=NULL)
	p.supp.succ.secondary<-  tmle.supp.succ.secondary$e

	N.supp.succ.secondary <-  sum(Y & use.these)  	
	c(p.supp.succ.secondary$pt, N.supp.succ.secondary/N.population)
	rm(Y)
	
	#--------------------------------------------------------
	# SUPPRESSION SUCCESSES * PRIMARY ********
	#--------------------------------------------------------
	#***** suppressed & ever ART & previous diagnosis & not outmigrated & not dead *** Primary****
	
	# Prob(supp*=1, eART=1, pDx=1, move=0, dead=0 | Y0=1)  
	Y<- as.numeric( (supp==1) & (eART==1) & (pDx==1) & (migrate==0) & (death==0) )	
	
	# adjusted analysis using deterministic knowledge about suppression
	# cannot have Z*=1 if intervention= (D=1 OR M=1 OR eART=0)
	combo<- as.numeric( (eART==0) | (death==1) | (migrate==1) )
	# add this to the predictor set
	pred.A.temp<- data.frame(pred.A, combo)
	
	# call ltmle to calculate suppression success
	tmle.supp.succ <- call.ltmle(data=data, t=t,  pred.A=pred.A.temp, A= TstVL, 
		Y= Y,	id=id, these= use.these, SL.library=SL.library, CV.size=CV.size, 
		deterministicQ=deterministicQ_combo)
	p.supp.succ<- tmle.supp.succ$e
	
	# # Side note: This is equivalent to the following
	# mean(!combo[use.these])*call.ltmle(data=data, t=t,  pred.A=pred.A, A= TstVL, 
		# Y= Y,	id=id, these= use.these & !combo, SL.library=SL.library, CV.size=CV.size)$e$pt
			
	# Now obtain an adjusted estimate of the number virally suppressed
	# Number estimated to be suppressed = (Estimated prob suppressed ) x (Population size)
	N.est.supp.succ <- round(p.supp.succ$pt*N.population, 0)

	# ---------------------------------------------------------------
	# SUPPRESSION FAILURE
	#'****not suppressed & ever ART & previous diagnosis & not outmigrated & not dead *** Primary****
	# Prob(supp*=0, eART=1, pDx=1, move=0, dead=0 | Y0=1)
	
	# calculate using knowledge that the probabilities are mutually exclusive and exhaustive
	#	inference via the delta method
	pt <- p.died$pt + p.migrate$pt + p.new.dx$pt + p.no.ART$pt + p.supp.succ$pt
	IC <- tmle.died$IC + tmle.migrate$IC + tmle.new.dx$IC + tmle.no.ART$IC + tmle.supp.succ$IC
	temp<- get.var(mu1=pt, IC1=IC)$est
	p.supp.failure <- data.frame( pt=(1-temp$psi), CI.lo=(1-temp$CI.hi), CI.hi=(1-temp$CI.lo))
	
	# # Side note: This is equivalent to the following
	# Y<- as.numeric( (supp==0) & (eART==1) & (pDx==1) & (migrate==0) & (death==0) )	
	# mean(!combo[use.these])*call.ltmle(data=data, t=t,  pred.A=pred.A, A= TstVL, 
		# Y= Y,	id=id, these= use.these & !combo, SL.library=SL.library, CV.size=CV.size)$e$pt
	
	# Now obtain an adjusted estimate of the number virally unsuppressed
	# Number estimated to be unsuppressed = (Estimated prob unsuppressed ) x (Population size)
	N.est.supp.failure <- round(p.supp.failure$pt*N.population, 0)
	
	#-------
	# Reformat & Clean
	RAW<- data.frame(rbind(N.population=N.population, N.death, N.migrate, 
		N.new.dx, N.no.ART, N.VL.missed, N.supp.failure.secondary, N.supp.succ.secondary, 
		N.est.supp.failure, N.est.supp.succ )) 

		
	FULL<- data.frame(rbind( p.died=p.died, p.migrate=p.migrate, p.new.dx=p.new.dx, p.no.ART=p.no.ART, 
		p.VL.missed=p.VL.missed, p.supp.fail.secondary= p.supp.failure.secondary, 
		p.supp.succ.secondary= p.supp.succ.secondary, 
		p.supp.fail= p.supp.failure, p.supp.succ= p.supp.succ))	


	colnames(RAW)<- paste('t',t, sep='')
	colnames(FULL)<- paste(colnames(FULL),t, sep='')

	list(full=FULL, raw=RAW)
}

#-------------------------------------------
# get.adjustment.long.no.int.cens: function to get the adjustment set for longitudinal
#	analysis with no intervention on censoring (Figure2)
#
# input: data, baseline predictors, variables
#
# output: list with adjustment set at t={0,1,2}
#--------------------------------------------

get.adjustment.long.no.int.cens <- function(data, baseline.pred, pDx.variable, eART.variable, tstvl.variable, supp.variable){

	
	#---------------------------------------------------------------------
	# To avoid adjustment for collinear variables, we drop (ever.test, delta) from the adjustment set
	# We also do not need to include past death/migration history 
	# The population is defined for y_0=1; so we do not need to adjust for being HIV+ 
	#		We drop (y_t, tsthiv_t) for all t
	# For t>0, we also drop prior dx bc this will always be 1 by definition of our target population
	#
	# For simplicity, we also leave out (dead_by, move_by, eART) & account for these during estimation
	#		See do.long.no.int.censor() function
	#--------------------------------------------------------------------

	obs.nodes0 <- c(pDx.variable, eART.variable, 'chc', 'tr', tstvl.variable, supp.variable, 'eversuppressed' )
	obs.nodes  <- c(eART.variable, 'chc', 'tr',  tstvl.variable, supp.variable, 'eversuppressed' )
	
	adj0<- baseline.pred
	adj1<- data.frame(baseline.pred, data[,paste(obs.nodes0, 0, sep='_')] )
	adj2<- data.frame(baseline.pred, data[,paste(obs.nodes0, 0, sep='_')] ,  data[,paste(obs.nodes, 1, sep='_')] )

	RETURN<- list(adj0=adj0, adj1=adj1, adj2=adj2)
	RETURN	
}		

