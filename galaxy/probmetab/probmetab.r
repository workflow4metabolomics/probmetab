#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file
# probmetab.r version="1.0.0"
# Author: Misharl Monsoor ABIMS TEAM mmonsoor@sb-roscoff.fr


# ----- LOG -----
log_file=file("probmetab.log", open = "wt")
sink(log_file)
sink(log_file, type = "out")

# ----- PACKAGE -----
cat("\tPACKAGE INFO\n")
pkgs=c("parallel","BiocGenerics", "Biobase", "Rcpp", "mzR", "igraph", "xcms","snow","CAMERA","batch","ProbMetab")
for(p in pkgs) {
	suppressWarnings( suppressPackageStartupMessages( stopifnot( library(p, quietly=TRUE, logical.return=TRUE, character.only=TRUE))))
	cat(p,"\t",as.character(packageVersion(p)),"\n",sep="")
}

source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
}
cat("\n\n")
# ----- ARGUMENTS -----
cat("\tARGUMENTS INFO\n") 
listArguments = parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects
write.table(as.matrix(listArguments), col.names=F, quote=F, sep='\t')

if (!is.null(listArguments[["zipfile"]])){
  zipfile= listArguments[["zipfile"]]; listArguments[["zipfile"]]=NULL
}

# ----- PROCESSING INFILE -----
cat("\tINFILE PROCESSING INFO\n")

# ----- INFILE PROCESSING -----

if(listArguments[["mode_acquisition"]]=="one") {
	load(listArguments[["xa"]])
	#Unzip the chromatograms file for plotting EIC pour the HTML file
	if(exists("zipfile"))
	{
		if (zipfile!=""){
			directory=unzip(zipfile)
		}
	}	
	if (!exists("xa")) {
		xa=xsAnnotate_object
	}
	source_local("lib.r")
	if (!exists("variableMetadata")) variableMetadata= getVariableMetadata(xa);
	
} else if(listArguments[["inputs_mode"]]=="two"){
	load(listArguments[["image_pos"]])
	#Unzip the chromatograms file for plotting EIC pour the HTML file
	if(exists("zipfile")) {
		if (zipfile!=""){
			directory=unzip(zipfile)
		}
	}
	if (!exists("xa")) {
		xa=xsAnnotate_object
	}
	xaP=xa
	source_local("lib.r")	
	if (!exists("variableMetadata")) variableMetadataP= getVariableMetadata(xa)
	else variableMetadataP=variableMetadata


	load(listArguments[["image_neg"]])
	#Unzip the chromatograms file for plotting EIC pour the HTML file
	if(exists("zipfile")) {
		
	 	if (zipfile!=""){
			directory=unzip(zipfile)
		}
	}
	if (!exists("xa")) {
		xa=xsAnnotate_object
	}
	xaN=xa
	source_local("lib.r")
	
	if (!exists("variableMetadata")) variableMetadataN= getVariableMetadata(xa)
	else variableMetadataN=variableMetadata
}

#Import the different functions
source_local("lib.r")
source_local("export.class.table-color-graph.R")

# ----- PROCESSING INFO -----
cat("\tMAIN PROCESSING INFO\n")

if(listArguments[["mode_acquisition"]]=="one") {
	results=probmetab(xa=xa, variableMetadata=variableMetadata,listArguments=listArguments)
} else if(listArguments[["inputs_mode"]]=="two"){
	results=probmetab(xaP=xaP, xaN=xaN,variableMetadataP=variableMetadataP, variableMetadataN=variableMetadataN, listArguments=listArguments)
}
#delete the parameters to avoid the passage to the next tool in .RData image
#rm(listArguments)
cat("\tDONE\n")
#saving R data in .Rdata file to save the variables used in the present tool
#save.image(paste("probmetab","RData",sep="."))

