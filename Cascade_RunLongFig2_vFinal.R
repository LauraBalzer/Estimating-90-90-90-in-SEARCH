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
# FOR LONGITUDINAL ANALYSIS OF KNOWN BL HIV+ ADULTS
#		with no intervention on censoring
#	(i.e. deaths & outmigrations are treated as cascade failures)
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
STRATA<- 'AllStrata'

# Primary analysis is running with SuperLearner (DO.SL=T)
# 		If false, then runs main terms glm for adjustment *** use for debugging ONLY
DO.SL <-  T

#-------------------------------------#---------------------------------------------------
# LOAD IN THE DATA SET
load('Cascade_v29Sept.Rdata')


#-------------------------------------#---------------------------------------------------
#  SPECIFY THE TARGET POPULATION  & SUBSET THE DATA

# Analysis is resricted to baseline enumerated + stable residents + aged 15+ at BL
these <- (data$enumerated_bl==1 & data$newstable_bl=='stable' & data$age_bl>14)
	
# Can subset now
print(sum(these))
data<- data[these, ]
baseline.pred<- baseline.pred[these,]

# drop an indicator of being an inmigrant from baseline predictors
baseline.pred <- subset(baseline.pred, select=-inmigrant )

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


#------------------------------------------------------------------------------------------
# FINAL QCs TO THE DATA SET

# In primary and secondary analyses use a suppressed VL as indicator of prior diagnosis and ever ART use
pDx.variable<- 'pDxvl'
eART.variable<- 'e_artvl'

# In primary and secondary analyses only using VL as measured at the CHC+tracking
tstvl.variable<- 'tst_vl_chc'
supp.variable<- 'supp_chc'

# In underlying dataset, missing HIV status is coded as 98
data$y_0 <- data$tsthiv_0*data$y_0
data$y_1 <- data$tsthiv_1*data$y_1
data$y_2 <- data$tsthiv_2*data$y_2

# Underlying dataset already coded  Supp = TstVL x Supp*
		
#-------------------------------------#---------------------------------------------------
# CREATE THE "MEAT" OF THE FILE NAME FOR OUTPUT	
file.middle<- paste( 
	paste('SuperLearner', DO.SL, sep=''), 
	paste('v', format(Sys.time(), "%d%b%Y"),sep='') , sep='.' )

#----------------------------------------


#----------------------------------------
# Loop through strata
#----------------------------------------
for(k in 1:length(strata)){

	this.strata<- strata[k]
	save(file.middle, file=paste( 'Cascade.LongFig2', paste('Strata', this.strata,sep=''), file.middle, 'Rdata', sep='.'))
		
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
		strata.indicator <-  data$age_bl<25			# to be 'young' need to be <25
		baseline.pred.temp<- subset(baseline.pred, select=-c(age.25to34, age.35to44, age.45plus))
	
	} else if(this.strata =='Old'){
		strata.indicator <- data$age_bl>24			# to be 'old' need to be > 24
		baseline.pred.temp <- baseline.pred
	
	} 
	

	# Get adjustment sets
	ADJ<- get.adjustment.long.no.int.cens(data=data, baseline.pred=baseline.pred.temp, 
			pDx.variable, eART.variable, tstvl.variable, supp.variable)
	
	#-------------------------------------#---------------------------------------------------
	# RUN ANALYSES FOR EACH TIME POINT
	
	# Throughout the household is the indpt unit
	id <- data$hhid
		
	out.t0<- do.long.no.int.censor(data=data, t=0, 
		pred.A=ADJ$adj0, these.obs= strata.indicator,
		pDx.variable= pDx.variable, eART.variable= eART.variable, 
		tstvl.variable= tstvl.variable, supp.variable= supp.variable,
		SL.library= SL.library, CV.size = CV.size, 
		id=id)
		
	out.t1<- do.long.no.int.censor(data=data, t=1, 
		pred.A=ADJ$adj1, these.obs= strata.indicator,
		pDx.variable= pDx.variable, eART.variable= eART.variable, 
		tstvl.variable= tstvl.variable, supp.variable= supp.variable,
		SL.library= SL.library, CV.size=CV.size,
		id=id)
				
	out.t2<- do.long.no.int.censor(data=data, t=2, 
		pred.A=ADJ$adj2, these.obs= strata.indicator,
		pDx.variable= pDx.variable, eART.variable= eART.variable, 
		tstvl.variable= tstvl.variable, supp.variable= supp.variable,
		SL.library= SL.library, CV.size=CV.size,
		id=id)
				
	#-------------------------------------#---------------------------------------------------
	# COMPILE THE RESULTS
	#-------------------------------------#---------------------------------------------------
	
	# make Figure2 and save the results
	Graph2 <- make.graph2(out.t0, out.t1, out.t2, this.strata)
    write.csv(Graph2, file=paste( 'Cascade.Long.Graph2', paste('Strata', this.strata,sep=''), file.middle, 'csv', sep='.')) 

	# make corresponding eTable and save the results
	full <-cbind(out.t0$full, out.t1$full, out.t2$full)
	raw <- cbind(out.t0$raw, out.t1$raw, out.t2$raw)	
	eTable2 <- make.eTable2(raw=raw, full=full)
	rownames(eTable2)<- paste(this.strata, rownames(eTable2), sep=':  ') 
	write.csv(eTable2, file=paste( 'Cascade.Long.eTable', paste('Strata', this.strata,sep=''), file.middle, 'csv', sep='.'))
	
	save(out.t0, out.t1, out.t2, eTable2, Graph2,  file=paste( 'Cascade.LongFig2', paste('Strata', this.strata,sep=''), file.middle, 'Rdata', sep='.'))

}


