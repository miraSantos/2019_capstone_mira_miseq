library(dada2) #packageVersion("dada2")
path <- "C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/sequencing_data/fastq_files/Euk_1391f_EukBr"


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`,3)


# Plotting Quality Proflies -----------------------------------------------
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Filtering and Trimming --------------------------------------------------
filt.path <- "C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/sequencing_data/filtered"
filtFs <- file.path(filt.path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt.path, paste0(sample.names, "_R_filt.fastq.gz"))

FWD_PRIMER_LEN <- 16
REV_PRIMER_LEN <- 24
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=FALSE,
                     compress=TRUE, multithread=FALSE,trimLeft = c(FWD_PRIMER_LEN, REV_PRIMER_LEN)) # On Windows set multithread=FALSE


head(out)  


# Learning Error Rates ----------------------------------------------------

errF <- learnErrors(filtFs, multithread=FALSE) #Should I make multithread = FALSE here? I tried it with =TRUE and it took 4ever (>8min)
errR <- learnErrors(filtRs, multithread=FALSE) #okay I tried it with multithread =FALSE and it still took like >11 min :D
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplication -----------------------------------------------------------

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# Sample Inference --------------------------------------------------------

dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE)
dadaFs[[1]]


# Merge Paired Reads ------------------------------------------------------

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])



