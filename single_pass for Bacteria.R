if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.8")


install.packages("C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/dada2-1.10",
                 repos = NULL,
                 type = "source",
                 dependencies = c("Depends", "Suggests","Imports"))


library("dada2") #packageVersion("dada2")
path <- "C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/sequencing_data/fastq_files/V4_515F_New_V4_806R_New"


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`,6)


# Plotting Quality Proflies -----------------------------------------------
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Filtering and Trimming --------------------------------------------------
filt.CA.path <- "C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/sequencing_data/filtered_CA"
dir.create(filt.CA.path)
filt.path <- "C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/sequencing_data/filtered"
filtFs_CA <- file.path(filt.CA.path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_CA <- file.path(filt.CA.path, paste0(sample.names, "_R_filt.fastq.gz"))

filtFs <- file.path(filt.path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt.path, paste0(sample.names, "_R_filt.fastq.gz"))

#RAN THROUGH CUTADAPT
#installing cutadapt
system("pip3 install --user --upgrade cutadapt")
setwd("C:/Users/mps48/AppData/Roaming/Python/Python36/Scripts")
system("cutadapt --help")

#check quality score encoding?
dir.create("C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/sequencing_data/filtered_cutadapt")
dir.create("C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/sequencing_data/filtered_cutadapt/Euk")
dir.create("C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/sequencing_data/filtered_cutadapt/Bact")

#-a is for 3' trimming
#-g is for 5' trimming
for (ii in 1:length(fnFs)){
  text <- paste("cutadapt -g GTGYCAGCMGCCGCGGTAA -a ATTAGAWACCCBNGTAGTCC  -q 2,2 -e 0.2 -l 150 -o ", filtFs_CA[ii]," ",fnFs[ii],sep="")
  system(text)
  print(paste(ii ," of ", length(fnFs)))
}

for (ii in 1:length(fnRs)){
  text <- paste("cutadapt -g GGACTACNVGGGTWTCTAAT -a TTACCGCGGCKGCTGRCAC  -q 2,2 -e0.2 -l 150  -o ", filtRs_CA[ii]," ",fnRs[ii],sep="")
  system(text)
  print(paste(ii ," of ", length(fnRs)))
}

# I go through dada2's filter because when I just fed the cutadapted files R crashed.
#notsure why this is happening
out <- filterAndTrim(filtFs_CA, filtFs, filtRs_CA, filtRs, rm.phix=FALSE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)


# Learning Error Rates ----------------------------------------------------

errF <- learnErrors(filtFs, multithread=FALSE) #Took ~15 min!
errR <- learnErrors(filtRs, multithread=FALSE) #okay I tried it with multithread =FALSE and it still took like >11 min :D
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

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


# Merge Paired Reads ------------------------------------------------------

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


# Construct Sequence Table ------------------------------------------------

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove Chimeras ---------------------------------------------------------

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)



# Track Reads Through Pipeline --------------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)



# Assign Taxonomy ---------------------------------------------------------

ref.path <- "C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/ref_databases/"
taxa <- assignTaxonomy(seqtab.nochim,paste0(ref.path,"silva_nr_v128_train_set.fa.gz"), multithread=FALSE)
taxa <- addSpecies(taxa, "silva_species_assignment_v128.fa.gz")



# Alternative Way w/ Decipher ---------------------------------------------

library(DECIPHER); packageVersion("DECIPHER")

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load(paste0(ref.path,"SILVA_SSU_r132_March2018.RData")) # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

# hand-off to phyloseq ----------------------------------------------------

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

samples.out <- rownames(seqtab.nochim)


#constructing phyloseq object

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),sample_data(samdf),tax_table(taxa))

