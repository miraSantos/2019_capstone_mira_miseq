
# INSTALLATION ------------------------------------------------------------

#DADA2 Installation
biocLite(suppressUpdates = FALSE)
biocLite("ShortRead", suppressUpdates = FALSE)

biocLite("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2") # use the ref="v1.10" (e.g.) argument to get specific versions

#loading dada2 for sequence inference
library(dada2);packageVersion("dada2")

#installing Biostrings for finding reverse and complement of primers
library(Biostrings)

#installing cutadapt for trimming
system("pip3 install --user --upgrade cutadapt")
#when you install cutadapt the console should list where cutadapt.exe was stored. copy this filepath name below
cutadapt.path <- 'C:/Users/mps48/AppData/Roaming/Python/Python36/Scripts'

#installing phyloseq for visualization
install.packages("phyloseq")
library(phyloseq); packageVersion("phyloseq")

#installing ggplot2 for plotting (phyloseq needs this)
install.packages("ggplot2")
library(ggplot2); packageVersion("ggplot2")

#installing RColorBrewer to generate sets of distinct colors for plots
install.packages("RColorBrewer")
library(RColorBrewer); packageVersion("RColorBrewer")

#installing MSA for mutiple sequence alignment with phyloseq
BiocManager::install("msa")
library("msa")

# USER INPUT -------------------------------------------------------
#1. Create a folder called Seq_DADA2 and set the path to it below (make sure to use forward slashes)
base.path <-"D:/Seq_DADA2"

#2. place your folder of fastq files (make sure the folder is called fastq_files) inside the Seq_DADA2 folder. It's okay if they are zipped.
path <- paste0(base.path,"/fastq_files/")

#3. Download your preferred reference database (http://benjjneb.github.io/dada2/training.html) 
#save it to the Seq_DADA2 folder and write down the file name below
ref.path <-paste0(base.path,"/ref_databases/silva_132.18s.99_rep_set.dada2.fa.gz")

#4. setting up sample filenames
#find the pattern for forward and reverse fastq filenames and enter it into the 2nd argument of the list.files() function
# e.g. Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

#The first step of the pipeline will create unique sample names for each fastq file by splitting the full name by underscores into segments. You need to select the segment that is unique.
strsplit(basename(fnFs), "_")
#pick the segment that is unique and set it's index below
sample.name.num <- 3
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`,sample.name.num)

#PRIMER IN STRING FORM (5'-> 3')
forwardP <- "GTACACACCGCCCGTC"   # e.g. forwardP <-"ATCGATACT"
reverseP <- "TGATCCTTCTGCAGGTTCACCTAC"   # e.g. forwardP <-"ACCTATAGC"

#SET QUALITY THRESHOLDS see (https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming) for more details
q.score <- "10,15" # qscore <-(5'end score),(3'end score)

# Plotting Quality Profiles -----------------------------------------------
dir.create(paste0(base.path,"/quality.profile.plots/"))

for (ii in 1:length(fnFs)){
  pdf(paste0(base.path,"/quality.profile.plots/","Forward-",sample.names[ii],".pdf"))
  plotQualityProfile(fnFs[ii])
  dev.off()
  print(paste0(ii," of ", length(fnFs)))
}

for (jj in 1:length(fnRs)){
  pdf(paste0(base.path,"/quality.profile.plots/","Reverse-",sample.names[jj],".pdf"))
  plotQualityProfile(fnRs[jj])
  dev.off()
  print(paste0(jj," of ", length(fnRs)))
}

# Filtering and Trimming --------------------------------------------------
dir.create(paste0(base.path,"/filtered"))
filt.path <- paste0(base.path,"/filtered")

dir.create(paste0(base.path,"/filtered_CA"))
filt.CA.path <- paste0(base.path,"/filtered_CA")

filtFs_CA <- file.path(filt.CA.path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_CA <- file.path(filt.CA.path, paste0(sample.names, "_R_filt.fastq.gz"))

filtFs <- file.path(filt.path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt.path, paste0(sample.names, "_R_filt.fastq.gz"))

#Routing through cutadapt
setwd(cutadapt.path)
system("cutadapt --help")

#-a is for 3' trimming
#-g is for 5' trimming
for (ii in 1:length(fnFs)){
  text <- paste0("cutadapt -g ",forwardP," -a ",paste0(reverseComplement(DNAStringSet(reverseP)))," -q 2,2 -e 0.2 -l 150 -o ", filtFs_CA[ii]," ",fnFs[ii])
  system(text)
  print(paste(ii ," of ", length(fnFs)))
}

for (ii in 1:length(fnRs)){
  text <- paste0("cutadapt -g ",reverseP," -a ",paste0(reverseComplement(DNAStringSet(forwardP))),"  -q 2,2 -e0.2 -l 150  -o ", filtRs_CA[ii]," ",fnRs[ii])
  system(text)
  print(paste(ii ," of ", length(fnRs)))
}

out <- filterAndTrim(filtFs_CA, filtFs, filtRs_CA, filtRs, rm.phix=FALSE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)


# Learning Error Rates ----------------------------------------------------
errF <- learnErrors(filtFs, multithread=FALSE) #Took ~15 min!
errR <- learnErrors(filtRs, multithread=FALSE) #okay I tried it with multithread =FALSE and it still took like >11 min :D

dir.create(paste0(base.path,"/error_plots/"))

pdf(paste0(base.path,"/error_plots/","errF.pdf"))
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf(paste0(base.path,"/error_plots/","errR.pdf"))
plotErrors(errR, nominalQ=TRUE)
dev.off()


#Dereplication -----------------------------------------------------------

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# Sample Inference --------------------------------------------------------

dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE)
dadaFs[[1]]
dadaRs[[1]]


# Merge Paired Reads ------------------------------------------------------

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


# Construct Sequence Table ------------------------------------------------

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove Chimeras ---------------------------------------------------------

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)


sum(seqtab.nochim)/sum(seqtab)

# Track Reads Through Pipeline --------------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


# ASSIGN TAXONOMY ---------------------------------------------------------

taxa <- assignTaxonomy(seqtab.nochim,ref.path, multithread=FALSE)

#Inspect the Assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#check things

length(which(taxa.print == "Dinophyceae")) #dinoflagellates 9
length(which(taxa.print == "Bacillariophyceae")) #diatoms 63
length(which(taxa.print == "Ochrophyta")) 

seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
mult <- msa(seqs, method="ClustalW", type="dna", order="input")


length(which(taxa.print == "Dinoflagellata")) #dinoflagellates 544

