#	KS Drugbank Test: Program that performs our Kolmogorov-Smirnov
#		test using the sets of drugs corresponding to a target
#		as a query set and our ranked list of drugs.
#
#	Author: Yu, Alvin
#
#	Purpose: To analyze whether there is a correlation between the
#		top scoring hypothesized drugs to inhibit atherosclerosis
#		and a specific human target pathway
#
#	Usage: A full ranked list of drugs "Combined Results(Only unique drugs).txt"
#		will be needed. "TargetsforKStest Drugbank.txt" will also be needed, which
#		has a list of targets taken from Drugbank and their corresponding
#		drugs. Likewise, "TargetsforKStest TTD.txt" will be needed if
#		an analysis on the targets for TTD is done. The sROC package 
#		is necessary for pvalue calculation.
#
# 
# This file is part of the drugscan software package. drugscan is free
# software: you can redistribute it and/or modify it under the terms of
# Apache Software License version 2.0 (and incorporated into this package
# distribution in the file LICENSE).
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the Apache Software License version 2.0 for 
# details.
#
# Copyright Alvin Yu and Stephen Ramsey, Oregon State University
# 2016.02
#

#two sets of files needed for the execution of the test
results <- read.table("Combined Results(Only unique drugs).txt",sep="\t",quote="",header=TRUE,stringsAsFactors=FALSE)

### enter x=1 if doing Drugbank analysis ###
### enter x=2 if doing TTD analysis      ###
if(x==1){
targets <- read.table("TargetsforKStest Drugbank.txt",sep="\t",quote="",header=FALSE,stringsAsFactors=FALSE)
}

if(x==2){
targets <- read.table("TargetsforKStest TTD.txt",sep="\t",quote="",header=FALSE,stringsAsFactors=FALSE)
}

table <- matrix(nrow=nrow(targets),ncol=4)
table[,1] <- targets[,1]
randescores <- c(1)
scores <- c(1)
escore <- c(1)
pvalue <- c(1)
ngenes <- c(1)

#for loop runs one time for each target tested
for(a in 1:nrow(targets)) {
	b <- 2
	query <- c(1)

	#get the drugs associated with the target being tested as a query set
	while(nchar(targets[a,b]) > 4) { 
		query[b-1] <- targets[a,b]
		b <- b + 1
	}

#list of the drugs in order of their rankings starting at rank number 1
drugrank <- results[,1]
pelements <- vector(mode="numeric",length=length(drugrank))
scores <- 0

#Kolmogorov-Smirnov Test
for(c in 1:1000) {
	
	pelements[] <- 0
	querysetpos <- match(query,drugrank)
	querysetpos <- querysetpos[!is.na(querysetpos)]
		#PHIT
	pelements[querysetpos] <- 1 
	pelements[querysetpos] <- pelements[querysetpos] / sum(pelements[querysetpos])
		#PMISS
	pelements[-querysetpos] <- 1 
	pelements[-querysetpos] <- -pelements[-querysetpos] / sum(pelements[-querysetpos])

	scores <- cumsum(pelements)

#lists of random drugs for query set for pvalue calculation
	query <- drugrank[sample(1:length(drugrank),length(query),replace=FALSE)]

#max of the difference will be the escore
	randescores[c] <- scores[which.max(abs(scores))]
}

escore[a] <- randescores[c] #first is actual escore
randescores <- randescores[-1]

#calculate pvalue here
library(sROC)
pvalue[a] <- kCDF(randescores,xgrid=escore[a])$Fhat
if(pvalue[a] > 0.5) {
pvalue[a] <- kCDF(-randescores,xgrid=-escore[a])$Fhat
}
pvalue[a] <- 2*pvalue[a]

ngenes[a] <- length(query)
rm(query)


}

table[,2] <- escore
table[,3] <- pvalue
table[,4] <- ngenes
