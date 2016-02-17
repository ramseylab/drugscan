#	synthetic_data_analysis.R: Program that creates synthetic data sets to 
#		determine the optimal weighting value of the Kolomogorov-
#		Smirnov test
#
#	Author: Yu, Alvin
#
#	Purpose: To analyze whether a weight, unweighted, or mixed weighted
#		type of Kolmogorov-Smirnov test will produce the most
#		accurate results.
#
#	Usage: "HLdifvalues.txt" and "HTHLdifvalues.txt" are needed and produced 
#		from	"CompareGeneExpressionDrugVsVehicle.R". "cmap_instances_02_onlyHL60.txt"
#		is also needed and can be downloaded from the Broad Institue
#		Connectivity Map. The output of this script will be ten AUC
#		values for a certain type of weighted, unweighted, or mixed
#		weighted, with a specified q value. Package gplots, gtools, gdata,
#		caTools, bitop, AUC, and ROCR will need to be loaded.
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

pelements <- c(1,2,3)
randescores <- c(1,2,3)
#p is the weighting exponent that will weight a gene's differential expression value
p <- 1
deltaescorecaseone1 <- matrix(nrow=1000,ncol=10)
deltaescorecasetwo1 <- matrix(nrow=1000,ncol=10)
#q is the weighted constant for the probability equation to select genes for a query set
q <- 0.3

myhllogfoldchanges <- read.table("HLdifvalues.txt", sep="\t", row.names=1, header=TRUE, skip=0, stringsAsFactors=FALSE)
myhthllogfoldchanges <- read.table("HTHLdifvalues.txt", sep="\t", row.names=1, header=TRUE, skip=0, stringsAsFactors=FALSE)
cmapinstances <- read.table("cmap_instances_02_onlyHL60.txt", sep="\t", quote="\"", header=TRUE, stringsAsFactors=FALSE, fill=TRUE)

names(myhllogfoldchanges) <- substring(names(myhllogfoldchanges),2,nchar(myhllogfoldchanges))
names(myhthllogfoldchanges) <- substring(names(myhthllogfoldchanges),2,nchar(myhthllogfoldchanges))

for(col in 1:10) {
for(i in 1:1000){

#randomly choose a drug
randdrugrow <- sample(1:nrow(cmapinstances),size=1)
randdrugid <- cmapinstances[randdrugrow,1]

#check which gene chip the specified drug was scanned on
if(cmapinstances[randdrugrow,8] == "HG-U133A") {
	mylogfoldchanges <- myhllogfoldchanges
} else {
	mylogfoldchanges <- myhthllogfoldchanges
}

randdrugcol <- match(randdrugid,names(mylogfoldchanges))

#get and sort the fold changes of the drug
randdrugfoldchanges <- na.omit(mylogfoldchanges[,randdrugcol])
rankedfoldchanges <- sort(randdrugfoldchanges)

#create a random up query set of genes using a weighted probability analysis
upprobvec <- 1/(1:length(randdrugfoldchanges))^q
upprobvec <- upprobvec/sum(upprobvec)
upquery <- sample(1:length(randdrugfoldchanges),size=80,prob=upprobvec,replace=FALSE)

#same produce for the down query set, no replacement
down <- setdiff(1:length(randdrugfoldchanges),upquery)
downprobvec <- 1/(length(randdrugfoldchanges)+1-down)^q
downprobvec <- downprobvec/sum(downprobvec)
downquery <- sample(down,size=80,prob=downprobvec,replace=FALSE)

#this will be the summing vector of the KS test
pelements <- vector(mode="numeric",length=length(randdrugfoldchanges))

#Kolmogorov-Smirnov Test done once for each up and down query set

#up query
        pelements[] <- 0
		#PHIT  
 	  pelements[upquery] <- abs(randdrugfoldchanges[upquery]) ^ p
        pelements[upquery] <- pelements[upquery]/ sum(pelements[upquery])
		#PMISS 
	  pelements[-upquery] <- 1/(length(upquery)-length(randdrugfoldchanges))
        scores <- cumsum(pelements)
	  upescore <- scores[which.max(abs(scores))]

#down query
        pelements[] <- 0
		#PHIT  
 	  pelements[downquery] <- abs(randdrugfoldchanges[downquery]) ^ 0
        pelements[downquery] <- pelements[downquery]/ sum(pelements[downquery])
		#PMISS 
	  pelements[-downquery] <- 1/(length(downquery)-length(randdrugfoldchanges)) 	  
        scores <- cumsum(pelements)
	  downescore <- scores[which.max(abs(scores))]

#Kolmogorov Smirnov test for our random queries
for(j in 1:2) {
        pelements[] <- 0

#create a random query of length 80
        randquery <- sample(1:length(randdrugfoldchanges),size=80,replace=FALSE)

		#PHIT
 	  pelements[randquery] <- abs(randdrugfoldchanges[randquery]) ^ p
        pelements[randquery] <- pelements[randquery]/ sum(pelements[randquery])
		#PMISS 
 	  pelements[-randquery] <- 1/(length(randquery)-length(randdrugfoldchanges))

        scores <- cumsum(pelements)


#max of the difference is our final escore
	randescores[j] <- scores[which.max(abs(scores))]
}

#case one is the using the random queries
deltaescorecaseone1[i,col] <- randescores[1]-randescores[2]
#case two uses the biased query sets
deltaescorecasetwo1[i,col] <- upescore-downescore
}
}

predictions <- matrix(nrow=2000,ncol=10)
labels <- matrix(nrow=2000,ncol=10)

predictions[1:1000,1:10] <- deltaescorecasetwo1
predictions[1001:2000,] <- deltaescorecaseone1
labels[1:1000,] <- 1
labels[1001:2000,] <- 0

library(ROCR)
pred1 <- prediction(predictions[,1],labels[,1])
pred2 <- prediction(predictions[,2],labels[,2])
pred3 <- prediction(predictions[,3],labels[,3])
pred4 <- prediction(predictions[,4],labels[,4])
pred5 <- prediction(predictions[,5],labels[,5])
pred6 <- prediction(predictions[,6],labels[,6])
pred7 <- prediction(predictions[,7],labels[,7])
pred8 <- prediction(predictions[,8],labels[,8])
pred9 <- prediction(predictions[,9],labels[,9])
pred10 <- prediction(predictions[,10],labels[,10])

#AUC is our final statistic of area under the sensitivity vs false positive error rate
auc1 <- performance(pred1,"auc")
auc2 <- performance(pred2,"auc")
auc3 <- performance(pred3,"auc")
auc4 <- performance(pred4,"auc")
auc5 <- performance(pred5,"auc")
auc6 <- performance(pred6,"auc")
auc7 <- performance(pred7,"auc")
auc8 <- performance(pred8,"auc")
auc9 <- performance(pred9,"auc")
auc10 <- performance(pred10,"auc")
