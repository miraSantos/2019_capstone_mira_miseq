
# INSTALLATION ------------------------------------------------------------
#CHECK TO MAKE SURE YOU HAVE LATEST VERSION OF R

#DADA2 Installation
chooseCRANmirror()
install.packages("BiocManager")

BiocManager::install("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2") 

#loading dada2 for sequence inference
library(dada2);packageVersion("dada2")

#installing Biostrings for finding reverse and complement of primers
library(Biostrings)

#installing phyloseq for visualization
BiocManager::install("phyloseq")
library(phyloseq); packageVersion("phyloseq")

#installing ggplot2 for plotting (phyloseq needs this)
install.packages("ggplot2")
library(ggplot2); packageVersion("ggplot2")

#installing RColorBrewer to generate sets of distinct colors for plots
install.packages("RColorBrewer")
library(RColorBrewer); packageVersion("RColorBrewer")

#installing MSA for mutiple sequence alignment with phyloseq
BiocManager:install("msa")
library("msa")

install.packages("phangorn")
library("phangorn")

#make sure you have latest version of python!! (python.org)

# USER INPUT -------------------------------------------------------
#make a folder called seq_dada2

#1. Create a folder called Seq_DADA2 and set the path to it below (make sure to use forward slashes)
base.path <-"D:/Seq_DADA2"

dir.create(base.path,"/cutadapt")
system(paste0("pip3 install --target=",base.path,"/cutadapt"," cutadapt"))


#2. place your folder of fastq files (make sure the folder is called fastq_files) inside the Seq_DADA2 folder. It's okay if they are zipped.
path <- paste0(base.path,"/fastq_files/")

#3. Download your preferred reference database (http://benjjneb.github.io/dada2/training.html) 
#save it to the Seq_DADA2 folder and write down the file name below
dir.create(paste0(base.path,"/ref_databases"))
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


#5. Setting up primers
#PRIMER IN STRING FORM (5'-> 3')
##18S PRIMERS
forwardP <- "GTACACACCGCCCGTC"   # e.g. forwardP <-"ATCGATACT"
reverseP <- "TGATCCTTCTGCAGGTTCACCTAC"   # e.g. forwardP <-"ACCTATAGC"

##16S PRIMERS
forwardP <- "GTACACACCGCCCGTC"   # e.g. forwardP <-"ATCGATACT"
reverseP <- "TGATCCTTCTGCAGGTTCACCTAC"   # e.g. forwardP <-"ACCTATAGC"


#6. Setting quality thresholds see (https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming) for more details
q.score <- "10,15" # qscore <-(5'end score),(3'end score)
e.score <- "0.2" #e.score e.g. 0.2
l.score <-"150" #l.score

#7. setting up metadata bout sampling
references <- read.csv("D:/Sequencing_Data/Station_Cast_Filter_References.csv",header = TRUE, sep = ",")

# Plotting Quality Profiles -----------------------------------------------
dir.create(paste0(base.path,"/quality.profile.plots/"))

for (ii in 1:length(fnFs)){
  plotQualityProfile(fnFs[ii])
  ggsave(paste0(base.path,'/quality.profile.plots/',"Forward-",sample.names[ii],'.pdf'))
  dev.cur()
  print(paste0(ii," of ", length(fnFs)))
}

for (jj in 1:length(fnRs)){
  plotQualityProfile(fnRs[jj])
  ggsave(paste0(base.path,'/quality.profile.plots/',"Reverse-",sample.names[jj],'.pdf'))
  dev.cur()
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
system(paste0(base.path,"/cutadapt/bin/cutadapt cutadapt --help"))

#-a is for 3' trimming
#-g is for 5' trimming
for (ii in 1:length(fnFs)){
  text <- paste0(base.path,"/cutadapt/bin/cutadapt -g ",forwardP," -a ",paste0(reverseComplement(DNAStringSet(reverseP)))," -q ",q.score," -e ",e.score," -l ",l.score," -o ",filtFs_CA[ii]," ",fnFs[ii])
  system(text)
  print(paste(ii ," of ", length(fnFs)))
}

for (ii in 1:length(fnRs)){
  text <- paste0(base.path,"/cutadapt/bin/cutadapt -g ",reverseP," -a ",paste0(reverseComplement(DNAStringSet(forwardP)))," -q ",q.score," -e 0.2 -l 150  -o ", filtRs_CA[ii]," ",fnRs[ii])
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

plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(base.path,"/error_plots/","errF.pdf"))
dev.off()

plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(base.path,"/error_plots/","errR.pdf"))
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
names(seqs) <- seqs



# Building the Phylogenetic Tree ------------------------------------------


mult <- msa(seqs, method="ClustalW", type="dna", order="input")

phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)


# Visualization with Phyloseq ---------------------------------------------

head(seqtab.nochim)

samples.out <- rownames(seqtab.nochim)
filter <- references$Filter 
site <- references$Station
cast <- references$Cast
day <- references$Day
samdf <- data.frame(filter=filter, site= site, cast=cast,day = day)
rownames(samdf) <- samples.out



ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),sample_data(samdf),tax_table(taxa),phy_tree(fitGTR$tree))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
taxa_names(ps) <- paste0("Seq", seq(ntaxa(ps)))

ntaxa(ps)
nsamples(ps)
sample_names(ps)[1:5]
rank_names(ps)
sample_variables(ps)
otu_table(ps)[1:5,1:5]
tax_table(ps)[1:5,1:5]
phy_tree(ps)
taxa_names(ps)[1:10]

# phy_tree ----------------------------------------------------------------

myTaxa = names(sort(taxa_sums(ps),decreasing=TRUE)[1:10])
ex1 = prune_taxa(myTaxa,ps)
plot(phy_tree(ex1),show.node.label = TRUE)

plot_tree(ex1, color = "site",label.tips="Phylum",ladderize = "left", justify="left",size="Abundance")

print(sums)
# Pre-Processing ----------------------------------------------------------

#Only OTUs wit ha mean greater than 10^-5 are kept
ps.1 = filter_taxa(ps, function(x) mean(x) > 1e-5, TRUE)

#removes samples with less than 20 reads
ps.2 = prune_samples(sample_sums(ps.1)>=20, ps.1)

#remove taxa not seen more than 3 times in the least 20% of the samples) protects gainst OTU with small mean and trivially large CV
ps.3 = filter_taxa(ps.2, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

#standardize abudance to median sequencing depth
total = median(sample_sums(ps.3))
standf = function(x, t=total) round(t * (x / sum(x)))
ps.4 = transform_sample_counts(ps.3, standf)
ps.4

ps.4.pr = subset_taxa(ps.4, Order !="Chloroplast")
ps.4.pr = subset_taxa(ps.4.pr, Class !="Chloroplastida")
ps.4.pr


#prints out percentages of taxa that were identified by the database
sums<-c(sum(is.na(taxa.print[,1]))/length(taxa.print[,1]),
        sum(is.na(taxa.print[,2]))/length(taxa.print[,2]),
        sum(is.na(taxa.print[,3]))/length(taxa.print[,3]),
        sum(is.na(taxa.print[,4]))/length(taxa.print[,4]),
        sum(is.na(taxa.print[,5]))/length(taxa.print[,5]),
        sum(is.na(taxa.print[,6]))/length(taxa.print[,6]),
        sum(is.na(taxa.print[,7]))/length(taxa.print[,7])
)

# Bar Plots ---------------------------------------------------------------


#PHYLUM.SITE.DAY
n.seq<-500
set.seed(1)
getPalette = colorRampPalette(brewer.pal(12, "Paired"),bias = 1,interpolate = "spline")
top <- names(sort(taxa_sums(ps.4.pr), decreasing=TRUE))[1:n.seq]
ps.top<- prune_taxa(top, ps.4.pr)
ps.top.glom <- tax_glom(ps.top,"Phylum")
ps.top.melted <- psmelt(ps.top.glom)
ps.top.melted$date.name <- factor(ps.top.melted$date.name, levels = c("Oct.4","Dec.4","Dec.6","Jan.16","Jan.17"))
ps.top.melted$filter.name <- factor(ps.top.melted$filter.name, levels = c("0.2","3","10"))
ggplot(ps.top.melted, aes(fill=Phylum,y=Abundance,x=factor(date.name)))+ 
  facet_wrap(~site.name.par,scales ="free_x")+
  geom_bar( stat="identity", position="fill")+
  scale_fill_manual(values = getPalette(8))+
  theme(axis.text.x = element_text(angle = 45,hjust=1),text=element_text(size=18))+
  xlab("Day")+
  ylab("Relative Abundance")

dir.create(paste0(base.path,"/plots"))
ggsave(paste0(base.path,"/plots/","phylum.site.day.pdf"))

#PHYLUM.SITE.FILTER
getPalette = colorRampPalette(brewer.pal(12, "Paired"),bias = 1,interpolate = "spline")
top <- names(sort(taxa_sums(ps.4.pr), decreasing=TRUE))[1:n.seq]
ps.top<- prune_taxa(top, ps.4.pr)
# ps.top.glom <- tax_glom(ps.top,"Phylum")
ps.top.melted <- psmelt(ps.top)
ps.top.melted$date.name <- factor(ps.top.melted$date.name, levels = c("Oct.4","Dec.4","Dec.6","Jan.16","Jan.17"))
ps.top.melted$filter.name <- factor(ps.top.melted$filter.name, levels = c("0.2","3","10"))
ps.top.melted$filter.name 
ggplot(ps.top.melted, aes(fill=Phylum,x=factor(filter.name),y=Abundance))+ 
  facet_wrap(~site.name.par,scales ="free_x")+
  geom_bar( stat="identity",position="fill")+
  scale_fill_manual(values = getPalette(8))+
  theme(text=element_text(size=18))+
  xlab(expression("Filter Size ("*mu*"m)"))+
  ylab("Relative Abundance")

ggsave(paste0(base.path,"/plots/","phylum.site.filter.pdf"))

# Ordination Plots --------------------------------------------------------
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.pca.bray <- ordinate(ps.prop, method="PCoA", distance="bray")

plot_ordination(ps.prop, ord.pcoa.bray, color="filter", title="Bray PCoA")+ facet_wrap(~day,3)
ggsave(paste0(base.path,"/plots/","pca.filter.day.pdf"))

plot_ordination(ps.prop, ord.pcoa.bray, color="day", title="Bray PCoA")+ facet_wrap(~site,1)
ggsave(paste0(base.path,"/plots/","pca.day.site.pdf"))

plot_ordination(ps.prop, ord.pcoa.bray, color="filter", title="Bray PCoA")+ facet_wrap(~site,1)
ggsave(paste0(base.path,"/plots/","pca.filter.site.pdf"))



