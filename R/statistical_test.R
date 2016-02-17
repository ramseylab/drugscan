#	statistical_test.R: Program that performs the Kolmogorov-Smirnov test
#		on a list of gene differential expression values of HL60 cells of drug
#		treatment versus a control treatment using an outside query set
#		of Entrez gene IDs pertaining to a specific disease.
#
#	Author: Yu, Alvin
#
#	Purpose: To obtain enrichment scores for drugs to analyze their
#		correlation or anticorrelation with the query set provided
#
#	Usage: Manually input a query set as vector 'theirgenesetnames' 
#		before running. User also needs to provide a table of differential
#		values of drugs treated on a cell line in entrez gene ID form.
#		Output will be a data matrix with columns that show Drug Instance
#		ID, Drug Name, Concentration, Enrichment Score, and P-value.
#		sROC package is necessary for pvalue calculation.
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


###enter query set as "theirgenesetnames" in entrez ID###

#variables
drugid <- 0
p <- 0
drugcol <- 0
pvalue <- c(1,2,3)
escore <- c(1,2,3)
scores <- c(1,2,3)
randescores <- c(1,2,3)

#differential values of HL60 genes in response to various drugs
myhllogfoldchanges <- read.table("HLdifvalues.txt", sep="\t", row.names=1, header=TRUE, skip=0, stringsAsFactors=FALSE)
myhthllogfoldchanges <- read.table("HTHLdifvalues.txt", sep="\t", row.names=1, header=TRUE, skip=0, stringsAsFactors=FALSE)
#perturbagen information, downloaded from Broad Institute website
cmapinstances <- read.table("cmap_instances_02_onlyHL60.txt", sep="\t", quote="\"", header=TRUE, stringsAsFactors=FALSE, fill=TRUE)

data <- matrix(nrow = nrow(cmapinstances), ncol = 7)
colnames(data) <- c("Instance ID","Drug Name","Concentration","Escore","Neg log10 P-value","RandEscore Mean","RandEscore SD")
#transferring Instance ID, Drug Name, and Concentration to final matrix
data[,1] <- cmapinstances[,1]
data[,2] <- cmapinstances[,3]
data[,3] <- cmapinstances[,5]

names(myhllogfoldchanges) <- substring(names(myhllogfoldchanges),2,nchar(myhllogfoldchanges))
names(myhthllogfoldchanges) <- substring(names(myhthllogfoldchanges),2,nchar(myhthllogfoldchanges))

drugidlist <- cmapinstances[,1]

#for loop runs through one time for each drug
for(h in 1:length(drugidlist)) {

myidlist <- c("a","b","c")
myidrank <- c("a","b","c")
querymatch <- c(1,2,3)
foldranknum  <- c(1,2,3)
mydrugfoldchanges <- c(1,2,3)
myrankedfoldchanges <- c(1,2,3)

#what drug ID we're working with
drugid <- drugidlist[h]

#check if drug id is in HL or HTHL
if(cmapinstances[h,8] == "HG-U133A") {
	mylogfoldchanges <- myhllogfoldchanges
} else {
	mylogfoldchanges <- myhthllogfoldchanges
}

#find which column the drug id is in
drugcol <- match(drugid,names(mylogfoldchanges))

#get entrez number
entrez <- unlist(strsplit(row.names(mylogfoldchanges), "_at"))
#list of probe ids in numerical order
myidlist <- as.integer(entrez)
#find all genes which do not have a measured value for the drug
nafind <- which(is.na(mylogfoldchanges[,drugcol]))

myidlist <- myidlist[-nafind]
mydrugfoldchanges <- mylogfoldchanges[-nafind,drugcol]

#rank fold changes high to low
foldranknum <- rank(-mydrugfoldchanges, na="keep")

#this for loop puts all IDs into ranked order
for(i in 1:length(foldranknum)) {
	myidrank[foldranknum[i]] <- myidlist[i] #most up reg to down reg
	myrankedfoldchanges[foldranknum[i]] <- mydrugfoldchanges[i]
}

myidrank <- as.numeric(myidrank)

randquery <- as.integer(theirgenesetnames) #first run is actual escore

ngenes <- length(myidrank)
pelements <- vector(mode="numeric", length=ngenes)
scores <- 0

#PERFORM KOLMOGOROV SMIRNOV TEST
#first time through will calculate the "real" escore
#every score after that uses random query sets for p-value analysis
for(j in 1:1000) {

        pelements[] <- 0

        querysetpos <- match(randquery, myidrank)
        querysetpos <- querysetpos[!is.na(querysetpos)]
		#PHIT  
 	  pelements[querysetpos] <- abs(mydrugfoldchanges[querysetpos]) ^ p
        pelements[querysetpos] <- pelements[querysetpos]/ sum(pelements[querysetpos])
		#PMISS 
	  pelements[-querysetpos] <- abs(mydrugfoldchanges[-querysetpos]) ^ p
        pelements[-querysetpos] <- -pelements[-querysetpos]/ sum(pelements[-querysetpos])
        scores <- cumsum(pelements)

#lists of 50 random probes to be used next time through this for loop
	randquery <- myidlist[sample(1:length(myidlist),length(theirgenesetnames),replace=FALSE)]

#max of the difference, which will be the escore
	randescores[j] <- scores[which.max(abs(scores))]
}

escore[h] <- randescores[1] #first is actual escore
randescores <- randescores[-1]

#extra statistical information
data[h,6] <- mean(randescores)
data[h,7] <- sd(randescores)

#pvalue calculation done here
library(sROC)
pvalue[h] <- kCDF(randescores, xgrid=escore[h])$Fhat
if (pvalue[h] > 0.5) {
pvalue[h] <- kCDF(-randescores, xgrid=-escore[h])$Fhat
}
pvalue[h] <- 2*pvalue[h]
}

#turn pvalue into -log10 form(many p-values are very low, so this makes it easier to analyze)
pvalue <- -log10(pvalue)

#add results to data table
data[,4] <- escore
data[,5] <- pvalue

#perform write.table on 'data' to export data

