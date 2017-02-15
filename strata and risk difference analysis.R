library(devtools)
install_github("joshuaschwab/ltmle", ref="a6d907ec76c3f432de15ccd636066da19ec9931a")
library(ltmle)
library(matrixStats)
library(SuperLearner)
end.time.set <- 1:2
outcome.str.set <- c("ever_testvl", "e_artvl", "supp_chc")

screen.corRank10 <- function(Y, X, family, ...) screen.corRank(Y, X, family, rank = 10, ...)
screen.corRank20 <- function(Y, X, family, ...) screen.corRank(Y, X, family, rank = 20, ...)

SL.library <- list( c('SL.glm', 'All'), c('SL.glm', 'screen.corRank20'), c('SL.glm', 'screen.corRank10'),
                    c('SL.step', 'screen.corRank20'), c('SL.step', 'screen.corRank10'),
                    c('SL.gam', 'All'), c('SL.gam', 'screen.corRank20'), c('SL.gam', 'screen.corRank10'),
                    c('SL.mean', 'All'))
SL.options <- list(SL.library=SL.library, verbose=F, cvControl=list(V= 5, stratifyCV=F))

strata.set <- c("All", "pDxvl_0 == 0", "pDxvl_0 == 0 & cd4_0 >= 500", "pDxvl_0 == 0 & cd4_0 >= 0 & cd4_0 < 500", "pDxvl_0 == 0 & cd4_0 < 0", "pDxvl_0 == 1 & e_artvl_0 == 0", "pDxvl_0 == 1 & e_artvl_0 == 0 & cd4_0 >= 500", "pDxvl_0 == 1 & e_artvl_0 == 0 & cd4_0 >= 0 & cd4_0 < 500", "pDxvl_0 == 1 & e_artvl_0 == 0 & cd4_0 < 0", "e_artvl_0 == 1", "e_artvl_0 == 1 & tst_vl_chc_0 == 1 & supp_chc_0 == 1", "e_artvl_0 == 1 & tst_vl_chc_0 == 1 & supp_chc_0 == 0", "e_artvl_0 == 1 & tst_vl_chc_0 == 0") 

category.set <- list(c("male", ref="female"), c("age.15to24", "age.25to34", "age.35to44", ref="age.45plus"), c(ref="not.single", "single"),  c(ref="less.than.primary", "primary", "secondary.plus"), c(ref="formal", "informal.hi.occup", "informal.low.occup", "other.job", "jobless"), c("wealth1", "wealth2", "wealth3", "wealth4", ref="wealth5"), c(ref="moAway.less", "moAway1.plus"))

results.names <- c("adj.est", "adj.CI.low", "adj.CI.high", "N", "unadj.numer", "unadj.denom", "unadj.est", "unadj.CI.low", "unadj.CI.high")
results <- array(dim=c(length(strata.set), max(end.time.set), length(outcome.str.set), length(results.names)))
dimnames(results) <- list(strata.set, paste0("t=", 1:max(end.time.set)), outcome.str.set, results.names)
results.RD <- array(dim=c(length(unlist(category.set)), length(outcome.str.set), length(results.names)))
dimnames(results.RD) <- list(unlist(category.set), outcome.str.set, results.names)

AddResult <- function(r, r.unadjusted, N.baseline, unadj.numer, unadj.denom, name, is.risk.difference) {
  LtmleObj <- function(r) {
    if (identical(r, "reference")) {
      rep(NaN, 3)
    } else {
      if (class(r) == "ltmle") {
        c(r$estimates["tmle"], summary(r)$treatment$CI)
      } else if (class(r) == "ltmleEffectMeasures") {
        x <- summary(r)$effect.measures$ATE
        c(-x$estimate, -x$CI[2], -x$CI[1]) #outputs are in terms of "never test/suppressed" but x is in terms of "ever"
      } else {
        stop("bad class")
      }
    }
  }
  
  temp.result <- c(LtmleObj(r), N.baseline, unadj.numer, unadj.denom, LtmleObj(r.unadjusted))
  if (is.risk.difference) {
    results.RD[name, outcome.str, ] <<- temp.result
  } else {
    results[name, end.time, outcome.str, ] <<- temp.result
  }
}

GetDeterministicQFunction <- function(this.outcome.str, this.end.time) {
  if (this.outcome.str == "supp_chc") {
    #if eARTvl_t=0 then Supp.chc_t is determinsitically 0
    function(data, current.node, nodes, called.from.estimate.g) {
      index <- which(names(data) == paste0("e_artvl_", this.end.time))
      stopifnot(length(index) == 1)
      in.history <- index < current.node
      if (! in.history) return(NULL)
      return(list(is.deterministic=data[, index] == 0 & !is.na(data[, index]), Q.value=0))  
    }
  } else {
    NULL
  }
}

ltmle1 <- function(ddata, Anodes, Cnodes, Lnodes, Ynodes, abar, survivalOutcome, stratify, id, name, is.reference) {
  last.Ynode <- Ynodes[length(Ynodes)]
  
  if (is.reference) {
    r <- r.unadjusted <- "reference"
    is.risk.difference <- T
  } else {
    is.risk.difference <- is.list(abar)
    
    gform <- GetForm(ddata, Anodes, Cnodes, Lnodes, Ynodes, is.Qform = F, stratify, survivalOutcome)
    Qform <- GetForm(ddata, Anodes, Cnodes, Lnodes, Ynodes, is.Qform = T, stratify, survivalOutcome)
    r <- ltmle(data=ddata, Anodes=Anodes, Cnodes=Cnodes, Lnodes=Lnodes, Ynodes=Ynodes, abar=abar, estimate.time = F, survivalOutcome = survivalOutcome, variance.method = "ic", stratify=stratify, id=id, gform=gform, Qform=Qform, SL.options = SL.options, deterministic.Q.function = GetDeterministicQFunction(outcome.str, end.time))
    
    ACnodes <- names(ddata)[names(ddata) %in% c(Anodes, Cnodes)] #maintains order of A/C nodes
    data.unadjusted <- ddata[, c(ACnodes, last.Ynode)]
    Qform.intercept <- "Q.kplus1 ~ 1"
    names(Qform.intercept) <- last.Ynode
    gform.intercept <- paste(ACnodes, "~ 1")
    r.unadjusted <- ltmle(data=data.unadjusted, Anodes=Anodes, Cnodes=Cnodes, Lnodes=NULL, Ynodes=last.Ynode, abar=abar, estimate.time = F, variance.method = "ic", id=id, Qform=Qform.intercept, gform=gform.intercept)
    
    if (is.risk.difference) {
      unadj.est <- summary(r.unadjusted)$effect.measures$ATE$estimate
    } else {
      unadj.est <- r.unadjusted$estimates["tmle"]
    }
    
    if (is.risk.difference) {
      abar <- abar[[1]] #if this is a Risk Difference, abar is a list; if this is an Treatment Specific Mean, abar is a vector
    }
  }
  
  Y <- ddata[, last.Ynode]
  condition.index <- rowAlls(ddata[, Cnodes, drop=F] == "uncensored") 
  if (!is.null(Anodes)) {
    condition.index <- condition.index & rowAlls(ddata[, Anodes, drop=F] == matrix(abar, nrow=nrow(ddata), ncol=length(abar), byrow = T))
  } 
  
  unadj.numer <- sum(Y[condition.index])
  unadj.denom <- length(Y[condition.index])
  N.baseline <- length(Y) #doesn't consider censoring 
  
  AddResult(r, r.unadjusted, N.baseline, unadj.numer, unadj.denom, name, is.risk.difference) 
}

GetNonReference <- function(x) {
  x[names(x) != "ref"]
}

GetForm <- function(ddata, Anodes, Cnodes, Lnodes, Ynodes, is.Qform, stratify, survivalOutcome) {
  nodes <- ltmle:::CreateNodes(ddata, Anodes, Cnodes, Lnodes, Ynodes) #these are indicies (not names)
  if (is.Qform) {
    lhs <- rep("Q.kplus1", length(nodes$LY))
    node.set <- nodes$LY
  } else {
    lhs <- names(ddata)[nodes$AC]
    node.set <- nodes$AC
  }
  if (stratify) {
    stratify.names <- c(Cnodes, Anodes)
  } else {
    stratify.names <- Cnodes
  }
  if (survivalOutcome) {
    stratify.nodes <- c(stratify.names, Ynodes)
  }
  form <- NULL
  for (i in seq_along(node.set)) {
    cur.node <- node.set[i]
    tst_vl_chc_t.node <- which(names(ddata) == paste0("tst_vl_chc_", end.time)) 
    if (outcome.str == "supp_chc" && cur.node %in% c(tst_vl_chc_t.node, max(nodes$LY))) { #remove if this is the tst_vl_chc_t.node in gform or the last Y node in Qform
      remove.names <- paste0(c("chc", "tr"), "_", end.time)
    } else {
      remove.names <- NULL
    }
    parent.node.names <- setdiff(names(ddata)[ltmle:::sseq(1, cur.node - 1)], c(stratify.names, remove.names))
    if (length(parent.node.names) == 0) {
      form[i] <- paste(lhs[i], "~ 1")
    } else {
      form[i] <- paste(lhs[i], "~", paste(parent.node.names, collapse=" + "))     
    }
    names(form)[i] <- names(ddata)[cur.node]
  }
  return(form)
}

MakeCat <- function(cats) {
  in.cat <- rep(F, nrow(baseline.pred))
  for (cat1 in cats) {
    stopifnot(cat1 %in% names(baseline.pred))
    in.cat <- in.cat | baseline.pred[, cat1]
  }
  return(as.numeric(!in.cat))
}

load("~/Dropbox/Cascade_ForJosh_v7Oct/Cascade_v29Sept.Rdata")
for (tt in 0:2) {
  y.na <- data[, paste0("y_", tt)] == 98
  data[y.na, paste0("y_", tt)] <- 0 
}
#change reference job to formal.hi.occup
baseline.pred$other.job <- MakeCat(c("formal.hi.occup", "informal.hi.occup", "informal.low.occup", "jobless"))
baseline.pred$formal.hi.occup <- NULL

#change reference age to 45plus
baseline.pred$age.15to24 <- MakeCat(c("age.25to34", "age.35to44", "age.45plus"))
baseline.pred$age.45plus <- NULL

#drop inmigrant
stopifnot(all(baseline.pred$inmigrant[data$enumerated_bl == 1 & data$newstable_bl == "stable" & data$age_bl >= 15] == 0))
baseline.pred$inmigrant <- NULL

for (end.time in end.time.set) {
  for (outcome.str in outcome.str.set) {
    cat("\n\n----------- end.time = ", end.time, " outcome = ", outcome.str, " ---------- \n")
    
    if (outcome.str == "ever_testvl") {
      data.index <- data$pDxvl_0 == 0
    } else {
      data.index <- data$y_0 == 1
    }
    data.index <- data.index & data$enumerated_bl == 1 & data$newstable_bl == "stable" & data$age_bl >= 15
    
    ddata.without.baseline <- data.frame(matrix(nrow = nrow(data), ncol = 0))
    for (tt in 0:end.time) {
      C <- data.frame(BinaryToCensoring(is.censored = data[, paste0("censor_by_", tt)]))
      names(C) <- paste0("C", tt)
      ddata.without.baseline <- cbind(ddata.without.baseline, C, data[, paste0(c("pDxvl", "e_artvl", "chc", "tr", "tsthiv", "ever_testvl", "y", "tst_vl_chc", "supp_chc", "eversuppressed"), "_", tt)])
    }
    
    if (outcome.str == "supp_chc") {
      Anodes <- paste0("tst_vl_chc_", end.time)  #only the last tst_vl_chc
      abar <- 1
      survivalOutcome <- F
      stratify <- T
      this.strata.set <- strata.set
    } else {
      Anodes <- NULL
      abar <- NULL
      survivalOutcome <- T
      stratify <- F
      this.strata.set <- "All"
    } 
    
    Cnodes <- paste0("C", 0:end.time)
    Ynodes <- paste0(outcome.str, "_", 0:end.time)
    max.Ynode.index <- which(colnames(ddata.without.baseline) == Ynodes[length(Ynodes)])
    ddata.without.baseline <- ddata.without.baseline[, 1:max.Ynode.index] #cut columns after last Y 
    ddata <- cbind(baseline.pred, ddata.without.baseline)
    Lnodes <- setdiff(names(ddata.without.baseline), c(Anodes, Cnodes, Ynodes))
    
    #within strata
    for (strata.name in this.strata.set) {
      cat("strata.name = ", strata.name, "\n")
      if (strata.name == "All") {
        strata.index <- rep(T, nrow(data))
      } else {
        strata.index <- eval(parse(text=strata.name), data)
      }
      ltmle1(ddata = ddata[data.index & strata.index, ], Anodes = Anodes, Cnodes = Cnodes, Lnodes = Lnodes, Ynodes = Ynodes, abar = abar, survivalOutcome = survivalOutcome, stratify = stratify, id = data$hhid[data.index & strata.index], name=strata.name, is.reference=F)
    }
    
    #Risk Difference
    if (end.time == 2) { 
      for (category.index in seq_along(category.set)) {
        this.category <- category.set[[category.index]]
        other.categories <- unlist(category.set[-category.index])
        Anode.set <- GetNonReference(category.set[[category.index]])
        ddata.RD <- cbind(baseline.pred[, GetNonReference(other.categories)], baseline.pred[, GetNonReference(this.category), drop=F], ddata.without.baseline)
        
        cat("\n Aset = ", Anode.set, "\n")
        num.Anodes <- length(Anode.set)
        Anodes2 <- c(Anode.set, Anodes)
        abar.ref <- c(rep(0, num.Anodes), abar)
        for (j in seq_along(this.category)) {
          if (names(this.category)[j] == "ref") {
            ltmle1(ddata = ddata.RD[data.index, ], Anodes=Anodes2, Cnodes=Cnodes, Lnodes=NA, Ynodes=Ynodes, abar=abar.ref, survivalOutcome=NA, stratify = NA, id=NA, name=this.category[j], is.reference=T)
          } else {
            abar1 <- rep(0, num.Anodes)
            abar1[Anode.set == this.category[j]] <- 1
            abar.RD <- list(c(abar1, abar), abar.ref)
            ltmle1(ddata = ddata.RD[data.index, ], Anodes=Anodes2, Cnodes = Cnodes, Lnodes = Lnodes, Ynodes = Ynodes, abar=abar.RD, survivalOutcome = survivalOutcome, stratify=stratify, id=data$hhid[data.index], name=this.category[j], is.reference=F)
          }
        }
      }
    }
  }
}

# consistency checks
GetReference <- function(x) x[names(x)=="ref"]
stopifnot(isTRUE(all.equal(results[, , , "unadj.numer"] / results[, , , "unadj.denom"], results[, , , "unadj.est"], scale=1, tolerance=1e-5)))
for (k in seq_along(outcome.str.set)) {
  for (category.index in seq_along(category.set)) {
    this.category <- category.set[[category.index]]
    p <- results.RD[GetNonReference(this.category), k, "unadj.numer"] / results.RD[GetNonReference(this.category), k, "unadj.denom"]
    p.ref <- results.RD[GetReference(this.category), k, "unadj.numer"] / results.RD[GetReference(this.category), k, "unadj.denom"]
    stopifnot(isTRUE(all.equal(-(p - p.ref), results.RD[GetNonReference(this.category), k, "unadj.est"], scale=1, tolerance=1e-5)))
  }
}
