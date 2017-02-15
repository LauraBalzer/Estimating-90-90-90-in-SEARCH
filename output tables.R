### Table 1
Format1 <- function(x) {
  paste0(";", sum(x), " (", sprintf("%2.1f", 100 * mean(x)), "%)")
}

Add <- function(x=NULL) {
  row.cnt <<- row.cnt + 1
  if (is.null(x)) {
    val <- ""
  } else if (is.numeric(x)) {
    val <- Format1(x)
  } else if (is.character(x)) {
    rownames(output)[row.cnt] <<- x
    val <- Format1(baseline.pred[data.index, x])
  } else {
    stop("unexpected x")
  }
  output[row.cnt, region] <<- val
  return(NULL)
}

baseline.pred$age.45plus <- MakeCat(c("age.15to24", "age.25to34", "age.35to44"))
baseline.pred$no.or.missing.education <- MakeCat(c("primary", "secondary.plus"))
baseline.pred$formal <- MakeCat(c("informal.hi.occup", "informal.low.occup", "other.job", "jobless"))
baseline.pred$wealth5 <- MakeCat(paste0("wealth", 1:4))

region.names <- c("SW-Uganda", "E-Uganda", "Kenya", "all")
num.output.rows <- 26
output <- matrix("??", nrow=num.output.rows, ncol=4)
colnames(output) <- region.names
rownames(output) <- paste0("noname", 1:num.output.rows)

for (region in region.names) {
  row.cnt <- 0
  data.index <- data$enumerated_bl == 1 & data$newstable_bl == "stable" & data$age_bl >= 15
  if (region != "all") {
    data.index <- data.index & data$region_name == region
  }
  cat("N = ", sum(data.index), "\n")
  Add(data$tsthiv_0[data.index])
  Add(data$y_0[data.index & data$tsthiv_0 == 1])
  Add("male")
  Add()
  Add("age.15to24")
  Add("age.25to34")
  Add("age.35to44")
  Add("age.45plus")
  Add("single")
  Add()
  Add("no.or.missing.education")
  Add("primary")
  Add("secondary.plus")
  Add()
  Add("formal")
  Add("informal.hi.occup")
  Add("informal.low.occup")
  Add("other.job")
  Add("jobless")
  Add()
  Add("wealth1")
  Add("wealth2")
  Add("wealth3")
  Add("wealth4")
  Add("wealth5")
  Add("moAway1.plus")
}
print(data.frame(output))


### Table 2
RowText <- function(row) {
  if (any(is.nan(row))) {
    sprintf("%d/%d (%2.1f%%);Ref;Ref", row["unadj.numer"], row["unadj.denom"], 100 * row["unadj.numer"]/row["unadj.denom"])
  } else {
    sprintf("%d/%d (%2.1f%%);%2.1f%% (%2.1f%%, %2.1f%%);%2.1f%% (%2.1f%%, %2.1f%%);", row["unadj.numer"], row["unadj.denom"], 100 * row["unadj.numer"]/row["unadj.denom"], row["unadj.est"], row["unadj.CI.low"], row["unadj.CI.high"], row["adj.est"], row["adj.CI.low"], row["adj.CI.high"])
  }
}

results.p <- results
percentage <- c("adj.est", "adj.CI.low", "adj.CI.high", "unadj.est", "unadj.CI.low", "unadj.CI.high")
results.p[, , , percentage] <- results.p[, , , percentage] * 100

for (outcome.str in outcome.str.set) {
  cat("\n\n", outcome.str, "\n")
  max.row <- if (outcome.str == "supp_chc") nrow(results.p) else 1
  for (row.index in 1:max.row) {
    cat(rownames(results.p)[row.index], ";", results.p[row.index, "t=1", outcome.str, "N"], ";")
    cat(RowText(results.p[row.index, "t=1", outcome.str, ]), RowText(results.p[row.index, "t=2", outcome.str, ]), "\n")
  } 
}



### Table 3, eTable2, eTable3
results.RD.p <- results.RD
results.RD.p[, , percentage] <- results.RD.p[, , percentage] * 100
results.RD.p[, , "unadj.numer"] <- results.RD.p[, , "unadj.denom"] - results.RD.p[, , "unadj.numer"]

for (outcome.str in outcome.str.set) {
  cat("\n\n", outcome.str, "\n")
  for (row.index in 1:nrow(results.RD.p)) {
    if (any(sapply(category.set, function (x) x[[1]] == rownames(results.RD.p)[row.index]))) cat("\n") #insert blank line at start of category
    cat(rownames(results.RD.p)[row.index], ";")
    cat(RowText(results.RD.p[row.index, outcome.str, ]), "\n")
  } 
}



