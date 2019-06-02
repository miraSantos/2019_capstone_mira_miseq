source("http://bioconductor.org/biocLite.R")
biocLite(suppressUpdates = FALSE)
biocLite("ShortRead", suppressUpdates = FALSE)

biocLite("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2") # use the ref="v1.10" (e.g.) argument to get specific versions

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.8")

library(dada2);packageVersion("dada2")
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
system("cd C:/Users/mps48/AppData/Roaming/Python/Python36/Scripts")
system("cutadapt --help")

#check quality score encoding?
dir.create("C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/sequencing_data/filtered_cutadapt")
dir.create("C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/sequencing_data/filtered_cutadapt/Euk")
dir.create("C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/sequencing_data/filtered_cutadapt/Bact")

#-a is for 3' trimming
#-g is for 5' trimming
for (ii in 1:length(fnFs)){
  text <- paste("cutadapt -a GTAGGTGAACCTGCAGAAGGATCA -g GTACACACCGCCCGTC -q 2,2 -e 0.2 -l 150 -o ", filtFs_CA[ii]," ",fnFs[ii],sep="")
  system(text)
  print(paste(ii ," of ", length(fnFs)))
}

for (ii in 1:length(fnRs)){
  text <- paste("cutadapt -a GACGGGCGGTGTGTAC  -g TGATCCTTCTGCAGGTTCACCTAC -q 2,2 -e0.2 -l 150  -o ", filtRs_CA[ii]," ",fnRs[ii],sep="")
  system(text)
  print(paste(ii ," of ", length(fnRs)))
}

# I go through dada2's filter because when I just fed the cutadapted files R crashed.
#notsure why this is happening
#i rerouted through dada2's filter. Sthing is happening that I am not understanding
#why can't i just go through filter and trims device.
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



# Assign Taxonomy ---------------------------------------------------------
set.seed(100)


setwd("C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/ref_databases/")
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=FALSE,tryRC = TRUE)

#OTHER DATABASES
taxa <- assignTaxonomy(seqtab.nochim, "gg_13_8_train_set_97.fa.gz", multithread=FALSE)
taxa <- assignTaxonomy(seqtab.nochim, "rdp_train_set_16.fa.gz", multithread=FALSE)

#add species (additonal step)
taxa <- addSpecies(taxa, "silva_species_assignment_v132.fa.gz")

#Inspect the Assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


#write data to a table and share
write.table(taxa.print, "F:/Sequencing_Data/tables_taxa/euks_taxa_print.txt", sep= "\t")
# Alternative Way w/ Decipher ---------------------------------------------
#I did not do this
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER", version = "3.8")

library(DECIPHER); packageVersion("DECIPHER")

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("~/tax/IDTaxa/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
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


# Evaluate accuracy -------------------------------------------------------

unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")


