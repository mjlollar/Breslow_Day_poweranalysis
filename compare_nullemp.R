#!/usr/bin/env Rscript

## Quick script to compared unordered null and emp pvalues from power analysis output

## Takes two inputs, csv files with one column listing pvalues for empirical (input 1) 
## and null (input 2) pvalues. 

#### Running as bash script
# Run: Rscript bd_get_pvalues.r <your_input_file>
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  stop("Zero arguments provided (two expected)).n", call.=FALSE)
}

##input
df.emp <- read.csv(args[1], header=FALSE)
df.null <- read.csv(args[2], header=FALSE)

#get 95%confidence
null.index <- order(df.null$V1)[50]
null.cutoff <- df.null$V1[null.index]
cat("Null 95% cutoff:",null.cutoff,"\n")

#get number of emp pavlues greater than null
emp.greater <- sum(df.emp$V1 < null.cutoff)

#calculate and output power percentage
total = (emp.greater / 10)
cat("Power: ",total,"%","\n",sep='')
