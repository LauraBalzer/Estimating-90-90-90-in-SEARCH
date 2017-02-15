#=========================================
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
# SETUP & RUN PRIMARY & SECONDARY ANALYSES
# FOR SERIAL X-SECT ANALYSIS OF PREVALENT HIV+
#==========================================
rm(list=ls())

# Requires the SuperLearner v2.0-2 
library('SuperLearner')

# Also requires ltmle v0.9-8-4
# Available by 
# library('devtools')
# install_github("joshuaschwab/ltmle", ref="a6d907ec76c3f432de15ccd636066da19ec9931a") 
library('ltmle')

# load in necessary functions
source('Cascade_Functions_vFinal.R')
source('Cascade_Helper_vFinal.R')

# set the seed so that each time run all code; get same answer
set.seed(1)

# Can specify if want to run code on 'AllStrata' or  only for a given strata
# Strata names: 'Pooled' (all strata together),'Uganda', 'Kenya', 'Male','Female', 'Young', 'Old'
STRATA<- 'Old'

# Primary analysis is running with SuperLearner (DO.SL=T)
# *** If false, then runs main terms glm for adjustment *** use for debugging ONLY
DO.SL <-  T


#-------------------------------------#---------------------------------------------------
# SPECIFY PRIMARY + SENSITIVITY ANALYSES

# primary: restrict to baseline enumerated and stable only (ENUM=T)
ENUM <-   T
# primary: do Not include self-report of prior diagnosis (PRIOR.DX.SELF=F)
PRIOR.DX.SELF <-  F

#-------------------------------------#---------------------------------------------------
# LOAD IN THE DATA SET
load('Cascade_v29Sept.Rdata')


#-------------------------------------#---------------------------------------------------
#  SPECIFY THE TARGET POPULATION  & SUBSET THE DATA
#
# Primary analysis restricts to baseline enumerated, stable residents (ENUM=T)
# In secondary analysis, allow for inmigration (and other non-stable residents)

if(ENUM){	
	# Primary analysis is restricted to baseline enumerated + stable residents
	these <- inmigrant.0 <- inmigrant.1 <- inmigrant.2 <- (data$enumerated_bl==1 & data$newstable_bl=='stable')
	
	# drop from baseline adjustment set
	baseline.pred.temp<- subset(baseline.pred, select=-inmigrant )

} else{	
	# Secondary analysis for Serial X-sectional includes non-stable residents + in-migrants
	# t=0: include baseline enumerated; 
	# t=1: include BL enumerated + FUY1 inmigrants
	# t=2: include BL enumerated + FUY1 inmigrants + FUY2 inmigrants
	inmigrant.0 <- data$enumerated_postbl=='BL Enumerated'
	inmigrant.1 <- (data$enumerated_postbl=='BL Enumerated' | data$enumerated_postbl=='FUP1 Inmigrant' ) 	
	inmigrant.2 <-  these <- rep(T, nrow(data))
}

# Can subset now
print(sum(these))
data<- data[these, ]
baseline.pred<- baseline.pred[these,]
inmigrant.0 <- inmigrant.0[these]
inmigrant.1 <- inmigrant.1[these]
inmigrant.2 <- inmigrant.2[these]

# Children can age into the analysis
age.bound.0 <- data$age_bl>14			# aged 15yo+ at t=0
age.bound.1 <- data$age_bl>13			# aged 15yo+ at t=1
age.bound.2 <- data$age_bl>12			# aged 15yo+ at t=2
#-------------------------------------#---------------------------------------------------


#-------------------------------------#---------------------------------------------------
# WITH THE ABOVE INPUT, SPECIFY IF STRATA-SPECIFIC ANALYSIS
if(STRATA=='AllStrata'){
	strata <- c('Pooled','Uganda', 'Kenya', 'Male','Female', 'Young', 'Old')
} else{
	strata <- STRATA 
}
#-------------------------------------#---------------------------------------------------

#-------------------------------------#---------------------------------------------------
# WITH THE ABOVE INPUT, SPECIFY THE SUPERLEARNER LIBRARY FOR PRIMARY ANALYSIS
if(DO.SL){ 
	# SuperLearner library
	SL.library<- list( c('SL.glm', 'All'), c('SL.glm', 'screen.corRank20'), c('SL.glm', 'screen.corRank10'),
		c('SL.step', 'screen.corRank20'), c('SL.step', 'screen.corRank10'),
		 c('SL.gam', 'All'),		 c('SL.gam', 'screen.corRank20'), c('SL.gam', 'screen.corRank10'),
		 c('SL.mean', 'All'))
	# Doing 5-fold cross-validation
	CV.size <- 5
}else{
	# If no SuperLearner library, specified run glm with main terms for adjusted analysis
	SL.library<- CV.size<- NULL
}
print(SL.library)

#-------------------------------------#---------------------------------------------------
#-------------------------------------#---------------------------------------------------
# FINAL QCs TO THE DATA SET

# In primary and secondary analyses use a suppressed VL as indicator of prior diagnosis and ever ART use
pDx.variable<- 'pDxvl'
eART.variable<- 'e_artvl'

# Sensitivity analysis also use self-reported prior diagnosis
if(PRIOR.DX.SELF){	
	# Classify a subject as previously diagnosed based on prior record OR self-report 
	# (only among confirmed HIV+)
	temp<- which(data$pDxvl_0==0 & data$prdx_0==1 & data$y_0==1); print( length(temp) )
	data[temp, 'pDxvl_0']<- 1
	temp<- which(data$pDxvl_1==0 & data$prdx_0==1 & data$y_0==1); print( length(temp) )
	data[temp, 'pDxvl_1']<- 1
	temp<- which(data$pDxvl_2==0 & data$prdx_0==1 & data$y_0==1); print( length(temp) )
	data[temp, 'pDxvl_2']<- 1
}

# In primary and secondary analyses only using VL as measured at the CHC+tracking
tstvl.variable<- 'tst_vl_chc'
supp.variable<- 'supp_chc'

# In underlying dataset, missing HIV status is coded as 98
# As in the analysis plan, we define observed HIV status Y_t = TstHIV_t x Y*_t
data$y_0 <- data$tsthiv_0*data$y_0
data$y_1 <- data$tsthiv_1*data$y_1
data$y_2 <- data$tsthiv_2*data$y_2

# In underlying data set, suppression is already coded as Supp_t = TstVL_t x Supp*_t
		
#-------------------------------------#---------------------------------------------------
# CREATE THE "MEAT" OF THE FILE NAME FOR OUTPUT	
file.middle<- paste( 
	paste('enum', ENUM, sep='' ), 
	paste('pSelf', PRIOR.DX.SELF, sep=''),
	paste('SuperLearner', DO.SL, sep=''), 
	paste('v', format(Sys.time(), "%d%b%Y"),sep='') , sep='.' )



# -------------------------------------#---------------------------------------------------
# LOOP THROUGH STRATA AND RUN SERIAL X-SECTIONAL ANALYSIS
#-------------------------------------#---------------------------------------------------

for(k in 1:length(strata)){
	
	this.strata<- strata[k]
	save(file.middle, file=paste( 'Cascade.Serial', paste('Strata', this.strata,sep=''), file.middle, 'Rdata', sep='.'))

	if(this.strata =='Pooled'){
		strata.indicator <- rep(T, nrow(data))
		baseline.pred.temp <- baseline.pred
	
	} else if(this.strata =='Uganda'){
		strata.indicator <- data$country=='Uganda'
		baseline.pred.temp<- subset(baseline.pred, select=-c(Uganda))

	} else if(this.strata =='Kenya'){
		strata.indicator <- data$country=='Kenya'
		baseline.pred.temp<- subset(baseline.pred, select=-c(Uganda))
	
	} else if(this.strata =='Male'){
		strata.indicator <-  baseline.pred$male
		baseline.pred.temp<- subset(baseline.pred, select=-c(male))

	} else if(this.strata =='Female'){
		strata.indicator <-  !baseline.pred$male
		baseline.pred.temp<- subset(baseline.pred, select=-c(male))

	} else if(this.strata =='Young'){
		strata.indicator <- rep(T, nrow(data))
		# can age out
		age.bound.0<- age.bound.0 & data$age_bl<25  		# to be 'young' need to be <25
		age.bound.1<- age.bound.1 & data$age_bl<24
		age.bound.2<- age.bound.2 & data$age_bl<23
		baseline.pred.temp<- subset(baseline.pred, select=-c(age.25to34, age.35to44, age.45plus))
	
	} else if(this.strata =='Old'){
		strata.indicator <- rep(T, nrow(data))
		# can age in
		age.bound.0<- data$age_bl>24			# to be 'old' need to be > 24
		age.bound.1<- data$age_bl>23
		age.bound.2<- data$age_bl>22
		baseline.pred.temp <- baseline.pred
	
	} 	
	

	#-------------------------------------#---------------------------------------------------
	# Serial cross-sectional analyses are restricted to being alive, not outmigrated & relevant population (i.e strata, age, residency)
		
	these_0 <-  strata.indicator  & age.bound.0  & inmigrant.0    
	these_1 <-  !data$censor_by_1  & strata.indicator & age.bound.1  & inmigrant.1 
	these_2 <-  !data$censor_by_2  & strata.indicator  & age.bound.2  & inmigrant.2  

	print( colSums(cbind(these_0, these_1, these_2)) )
	head(baseline.pred.temp)
	
	#-------------------------------------#---------------------------------------------------
	# RUN ANALYSES FOR EACH TIME POINT
	# Throughout the household is the indpt unit
	id <- data$hhid	
	
	out.t0 <- do.serial.analysis(data=data, baseline.pred= baseline.pred.temp, 
			t=0, these.obs= these_0,
			pDx.variable= pDx.variable, eART.variable= eART.variable, 
			tstvl.variable= tstvl.variable, supp.variable= supp.variable,
			SL.library= SL.library, CV.size= CV.size, 
			id=id)
	
	out.t1 <- do.serial.analysis(data=data, baseline.pred= baseline.pred.temp, 
			t=1, these.obs= these_1,
			pDx.variable= pDx.variable, eART.variable= eART.variable, 
			tstvl.variable= tstvl.variable, supp.variable= supp.variable,
			SL.library= SL.library, CV.size= CV.size, 
			id=id)
	
	out.t2 <- do.serial.analysis(data=data, baseline.pred= baseline.pred.temp, 
			t=2, these.obs= these_2,
			pDx.variable= pDx.variable, eART.variable= eART.variable, 
			tstvl.variable= tstvl.variable, supp.variable= supp.variable,
			SL.library= SL.library, CV.size= CV.size, 
			id=id)
	
	#-------------------------------------#---------------------------------------------------
	# COMPILE THE RESULTS
	#-------------------------------------#---------------------------------------------------
	
	# make Figure1 and save the results
	Graph1 <- make.graph1(out.t0, out.t1, out.t2, this.strata)	
    write.csv(Graph1, file=paste( 'Cascade.Serial.Graph', paste('Strata', this.strata,sep=''), file.middle, 'csv', sep='.')) 

	# make corresponding eTable and save the results
	raw <-  cbind(out.t0$RAW,  out.t1$RAW, out.t2$RAW)	
	Cascade <- cbind(out.t0$Cascade,  out.t1$Cascade, out.t2$Cascade)
	Cascade.secondary <- cbind(out.t0$Cascade.secondary,  out.t1$Cascade.secondary, out.t2$Cascade.secondary)
	eTable1 <- make.eTable1(raw=raw, Cascade.secondary=Cascade.secondary, Cascade=Cascade)
	rownames(eTable1)<- paste(this.strata, rownames(eTable1), sep=':  ') 
	write.csv(eTable1, file=paste( 'Cascade.Serial.eTable', paste('Strata', this.strata,sep=''), file.middle, 'csv', sep='.'))
	
	save(out.t0, out.t1, out.t2, eTable1, file=paste( 'Cascade.Serial', paste('Strata', this.strata,sep=''), file.middle, 'Rdata', sep='.'))
		
}

