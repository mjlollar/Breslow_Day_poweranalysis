#!/usr/bin/env Rscript

#### Running as bash script
# Run: Rscript bd_get_pvalues.r <your_input_file>
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
	stop("Zero arguments provided (one expected).n", call.=FALSE)
}

setwd("/home/matt/bd_scan_work_2023")
library("DescTools")
df <- read.csv('11_xa_null_bd_cells.csv', header=TRUE)
num_windows <- nrow(df) 

pvalues <- vector(mode="numeric", length=num_windows) #Initialize vector for pvalues
### Calculate BD test for each window comparison
for (i in 1:num_windows){
  ## Cell values with Fisher adjustment to prevent divide by zero occurrences
  bd_1 <- df$bd1[i] + 0.5
  bd_2 <- df$bd2[i] + 0.5
  bd_3 <- df$bd3[i] + 0.5
  bd_4 <- df$bd4[i] + 0.5
  bd_5 <- df$bd5[i] + 0.5
  bd_6 <- df$bd6[i] + 0.5
  bd_7 <- df$bd7[i] + 0.5
  bd_8 <- df$bd8[i] + 0.5
  ## Calculate Odds 
  odds_one <- ((bd_1) / (bd_1 + bd_2))
  odds_two <- ((bd_3) / (bd_3 + bd_4))
  odds_three <- ((bd_5) / (bd_5 + bd_6))
  odds_four <- ((bd_7) / (bd_7 + bd_8))
  max_odds <- max(c(odds_one, odds_two, odds_three, odds_four))
  if (max_odds %in% c(odds_two, odds_three, odds_four) == TRUE){
    #Skip calculation if max odds are not odds 1 or equal to odds one
    pvalues[i] <- 999
  } else {
    ## Calculate Breslow Day test
    bd_table <- xtabs(freq ~ ., cbind(expand.grid(phenotype=c("sterile", "fertile"),
                                                  window_one=c("focal", "non-focal"),
                                                  window_two=c("focal", "non-focal")),
                                      freq=c(bd_1, bd_2, bd_3, bd_4, bd_5, bd_6, bd_7, bd_8)))
    y <- BreslowDayTest(bd_table)$p.value
    pvalues[i] <- y
  }
}

### Print lowest pvalue to output
lowest_pvalue <- as.character(min(pvalues))
out_name <- paste(as.character(args[1]), '_null_lowest_pvalue.txt',sep='')
out_file <- file(out_name)
writeLines(lowest_pvalue, out_file)
#### Warning check (Should not occur)
if (any(is.logical(pvalues) == TRUE)){
  writeLines('     - Warning: The number of tests performed did not match the number of window comparisons')
}
close(out_file)
