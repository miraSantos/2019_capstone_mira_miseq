# dwnlding bioconductor ---------------------------------------------------------
source("http://bioconductor.org/biocLite.R")
biocLite(suppressUpdates = FALSE)
biocLite("ShortRead", suppressUpdates = FALSE)

# dwnlding  devtools and dada2 ---------------------------------------------------------
biocLite("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2") # use the ref="v1.10" (e.g.) argument to get specific versions

library("dada2")

help(package="dada2")
?derepFastq
?dada

#following this tutorial: https://benjjneb.github.io/dada2/tutorial.html
# set up basepath  ---------------------------------------------------------
base.path = "C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/2019_capstone_mira/"
list.files(base.path)
fq.base.path <- paste(base.path,"/fastq_files/",sep="")
setwd(fq.base.path)
subfile.names <- list.files(fq.base.path)
subfile.names

# set up subpaths ---------------------------------------------------------
sub.path <- c(vector(),1:length(subfile.names))
fq.file.names <- rep(list(list(list(),list())),length(sub.path))
names(fq.file.names) <- subfile.names
for(ii in 1:length(subfile.names)){
  sub.path[ii]<- paste(fq.base.path,subfile.names[ii],sep = "")
}

# set up fastqfilenames ---------------------------------------------------------
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
for (ii in 1: length(sub.path)){
  fnFs <- sort(list.files(sub.path[ii], pattern="R1.fastq", full.names = TRUE))
  fnRs <- sort(list.files(sub.path[ii], pattern="R2.fastq", full.names = TRUE))
  fq.file.names[[ii]][[1]] <- fnFs
  fq.file.names[[ii]][[2]] <- fnRs
  names(fq.file.names[[ii]]) <- c("fnFs","fnRs")  
}

# create_sample_names -----------------------------------------------------
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 6)

for(i in 1:length(sub.path)){
sample.names <- sapply(strsplit(basename(fq.file.names[[i]]$fnFs), "_"), `[`, 6)
}
#inspect quality
plotQualityProfile(fnFs[12])

#create new folder called filtered
dir.create(paste(base.path,"/filtered",sep=""))

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(base.path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(base.path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

