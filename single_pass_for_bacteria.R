
# Installing Dada2 and setting up fastq filename paths --------------------


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

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.8")
library(dada2)

# ref.path <- "C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/ref_databases/"
ref.path <- "F:/Sequencing_Data/ref_databases/"

#note if you get error: cannot allocate vector of size 2.1 Gb. try gc() to cler up some memory without deleting anything
#in the global enviroment

#SILVA_v128
taxa <- assignTaxonomy(seqtab.nochim,paste0(ref.path,"silva_nr_v128_train_set.fa.gz"), multithread=FALSE)
taxa <- addSpecies(taxa,paste0(ref.path, "silva_species_assignment_v128.fa.gz") )

#SILVA_v132
taxa <- assignTaxonomy(seqtab.nochim,paste0(ref.path,"silva_nr_v132_train_set.fa.gz"), multithread=FALSE)
taxa <- addSpecies(taxa, paste0(ref.path,"silva_species_assignment_v132.fa.gz"))

#RDPa
taxa <- assignTaxonomy(seqtab.nochim,paste0(ref.path,"rdp_train_set_16.fa.gz"), multithread=FALSE)
taxa <- addSpecies(taxa,paste0(ref.path, "rdp_species_assignment_16.fa.gz"))

#GreenGenes
taxa <- assignTaxonomy(seqtab.nochim,paste0(ref.path,"gg_13_8_train_set_97.fa.gz"), multithread=FALSE)
taxa <- addSpecies(taxa,paste0(ref.path, "rdp_species_assignment_16.fa.gz"))

#Inspect the Assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

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
loading in phyloseq libraries -------------------------------------------


source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")



# Setting up variables ----------------------------------------------------


references <- read.csv("D:/Sequencing_Data/Station_Cast_Filter_References.csv",header = TRUE, sep = ",")


samples.out <- rownames(seqtab.nochim)
filter <- references$Filter 
site <- references$Station
cast <- references$Cast
day <- references$Day
samdf <- data.frame(filter=filter, site= site, cast=cast,day = day)
rownames(samdf) <- samples.out


# Constructing Phyloseq Object --------------------------------------------

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),sample_data(samdf),tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps

## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 232 taxa and 19 samples ]
## sample_data() Sample Data:       [ 19 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 232 taxa by 7 taxonomic ranks ]

# Bar_Plots ---------------------------------------------------------------

ps.din = subset_taxa(ps, Order == "")
plot_bar(ps.din)+geom_bar( stat="identity", position="stack")

ps.dia = subset_taxa(ps, Family == ")
plot_bar(ps.dia)
plot_bar(ps.dia)+geom_bar( stat="identity", position="stack")


# Plotting richness -------------------------------------------------------
plot_richness(ps, x="filter", measures=c("Shannon", "Simpson"))
plot_richness(ps, x="site", measures=c("Shannon", "Simpson"))
plot_richness(ps, x="cast", measures=c("Shannon", "Simpson"))
plot_richness(ps, x="day", measures=c("Shannon", "Simpson"))


# Ordination --------------------------------------------------------------

?ordinate()
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
ord.pcoa.bray <- ordinate(ps.prop, method="PCoA", distance="bray")
ord.cca.bray <- ordinate(ps.prop, method="CCA", distance="bray")
ord.dca.bray <- ordinate(ps.prop, method="DCA", distance="bray")
ord.rda.bray <- ordinate(ps.prop, method="RDA", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="day", title="Bray NMDS")
plot_ordination(ps.prop, ord.pcoa.bray, color="filter", title="Bray PCoA")+ facet_wrap(~day,3)
plot_ordination(ps.prop, ord.cca.bray, color="day", title="Bray CCA")
plot_ordination(ps.prop, ord.dca.bray, color="day", title="Bray DCA")
plot_ordination(ps.prop, ord.rda.bray, color="site", title="Bray DCA")

plot_ordination(ps.prop, ord.pcoa.bray, color="filter", title="Bray PCoA")+ facet_wrap(~day,3)
plot_ordination(ps.prop, ord.pcoa.bray, color="day", title="Bray PCoA")+ facet_wrap(~site,1)
plot_ordination(ps.prop, ord.pcoa.bray, color="filter", title="Bray PCoA")+ facet_wrap(~site,1)
plot_ordination(ps.prop, ord.pcoa.bray, color="", title="Bray PCoA")+ facet_wrap(~site,1)



# JUST OTU Plotting -------------------------------------------------------

p1 = plot_ordination(ps.prop,ord.pcoa.bray, type="taxa", color="Phylum", title="taxa")
print(p1)
print(p1+facet_wrap(~Phylum,3))

plankton = get_variable(ps, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")

p2 = plot_ordination(ps, ps.prop, type="samples", color="site", shape="filter") 
p2 + geom_polygon(aes(fill=SampleType)) + geom_point(size=5) + ggtitle("samples")

# Sorting Taxa ------------------------------------------------------------
if you don't sort and prune, it takes forever, data set is too much to plot'
plot_bar(ps.prop, x="day", fill="Family") + facet_wrap(~site, scales="free_x")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="day", fill="Family") + facet_wrap(~site, scales="free_x")

myPalette <- brewer.pal(12,"Paired")
top50 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:50]
ps.top50 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top50 <- prune_taxa(top50, ps.top50)
plot_bar(ps.top50, x="day", fill="Order") + facet_wrap(~site, scales="free_x")  + geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack") 

top100 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:100]
ps.top100 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top100 <- prune_taxa(top100, ps.top100)
plot_bar(ps.top100, x="day", fill="Family") + facet_wrap(~site, scales="free_x") + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")


# Transforming data -------------------------------------------------------



# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="day", title="Bray NMDS")

plot_heatmap(gpt, sample.label="SampleType")

