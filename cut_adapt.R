#install cutadapt with pip3 if you haven't installed it yet
system("pip3 install --user --upgrade cutadapt")
setwd("C:/Users/mps48/AppData/Roaming/Python/Python36/Scripts")
system("cutadapt --help")

#check quality score encoding?

dir.create("C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/sequencing_data/filtered/Euk")

dir.create("C:/Users/mps48/Documents/1111AAA_2019_Spring/Capstone/Sequencing/sequencing_data/filtered/Bact")

for (ii in 1:length(fnFs)){
  text <- paste("cutadapt -a GTAGGTGAACCTGCAGAAGGATCA -g GTACACACCGCCCGTC -q 2,2 -o ", filtFs[ii]," ",fnFs[ii],sep="")
system(text)
print(paste(ii ," of ", length(fnFs)))
}

for (ii in 1:length(fnRs)){
  text <- paste("cutadapt -a GTAGGTGAACCTGCAGAAGGATCA -g GTACACACCGCCCGTC -q 2,2 -o ", filtRs[ii]," ",fnRs[ii],sep="")
  system(text)
  print(paste(ii ," of ", length(fnRs)))
}

# for (ii in 1:length(fnRs)){
#   text <- paste("cutadapt -a AGACCAAGTCTCTGC -o ", filtFs1[ii]," ",fnFs[ii],sep="")
#   system(text)
# }
