##################################################################################################
##                                  Daniel Elias Martin Herranz                                 ##
##                                 MPhil in Computational Biology                               ##
##                                    University of Cambridge                                   ##
##                                          August 2015                                         ##
##################################################################################################

##### Create a FASTA file which contains the mitochondrial sequences with the filtered SNVs for 
#     each one of the samples. This FASTA file will be used afterwards to perform multiple sequence
#     alignment using Clustal Omega #####

setwd('/Users/dmh/Desktop/trial_caveman/mit_processing/extra_stuff/phyml/');
library('Biostrings');


## Read the reference sequence for the mitochondrial genome.

MT.genome <- as.character(unlist(strsplit(as.character(readDNAStringSet(filepath='/Users/dmh/Desktop/trial_caveman/mit_processing/extra_stuff/mutational_signature/Sarcophilus_harrisii.DEVIL7.0.70.dna.toplevel.fa',
                                                                        skip=27888, nrec=1)), '')));  # Each base is an element in a character vector.



## Read the filtered SNVs information and format it.

SNV.raw <- read.csv('/Users/dmh/Desktop/trial_caveman/mit_processing/baf_plots_filtered/filtered_variants.csv');
SNV.info <- SNV.raw[,2:ncol(SNV.raw)];
rownames(SNV.info) <- as.character(SNV.raw[,1]);
colnames(SNV.info) <- paste0(1:(length(SNV.raw)-1), ':', gsub('^X[0-9]+\\.', '', colnames(SNV.raw)[2:length(colnames(SNV.raw))]));

## Create a list for the different samples. Each elements contains the variants in that sample (e.g. 13759G>A, 11788G>A).
# The sample names are stored in the list names (e.g. 1:107T).

list.SNVs <- list();

for(i in 1:ncol(SNV.info)){
  
  sample <- colnames(SNV.info)[i];
  variants <- rownames(SNV.info)[which(SNV.info[,i] != '*')];
  list.SNVs[[i]]<- variants;
  names(list.SNVs)[i] <- sample;
   
}


## Create a list which contains the different mitochondrial sequences incorporating the SNVs.

list.mit.sequences <- list();

for(i in 1:length(list.SNVs)){
  
  sample <- names(list.SNVs)[i];
  variants <- list.SNVs[[i]];
  new.sequence <- MT.genome;
  
  for(j in 1:length(variants)){
  
    original.base <- substr(variants[j], nchar(variants[j])-2, nchar(variants[j])-2);
    new.base <- substr(variants[j], nchar(variants[j]), nchar(variants[j]));
    position <- as.numeric(substr(variants[j], 1, nchar(variants[j])-3));
    
    if(new.sequence[position] == original.base){
      
      new.sequence[position] <- new.base;
      
    }else{
      
      print('Error with original base in reference.')
    } 
  }
  
  list.mit.sequences[[i]] <- new.sequence;
  names(list.mit.sequences)[i] <- sample;
  
}


## Convert the list to a vector. Now the DNA sequences are strings and not vectors.

vector.mit.sequences <- sapply(1:length(list.mit.sequences), 
                               function(x){
                                 
                                 paste(list.mit.sequences[[x]], sep="", collapse="");
                                 
                               });


## Export the vector of sequences as FASTA file. The identifiers are the sample names.

set.sequences <- DNAStringSet(vector.mit.sequences);
names(set.sequences) <- names(list.mit.sequences);
writeXStringSet(set.sequences, file="/Users/dmh/Desktop/trial_caveman/mit_processing/extra_stuff/phyml/mit_sequences_with_SNVs.fa",format="fasta");

##### End of the script #####
