##################################################################################################
##                                  Daniel Elias Martin Herranz                                 ##
##                                 MPhil in Computational Biology                               ##
##                                    University of Cambridge                                   ##
##                                          August 2015                                         ##
##################################################################################################

##### Create several plots (mutational signature, VEP plot, rainfall plot) for the final 
#     filtered mitochondrial SNVs in Tasmanian devil samples #####

setwd('/Users/dmh/Desktop/trial_caveman/mit_processing/extra_stuff/list_of_variants/');

## Read the information from the variants.

raw.info <- read.csv('list_of_variants.csv');

# Exclude SNVs which are G/RE.

raw.info.no.re <- raw.info[-(which(raw.info[,3] =='G/RE')),];

## Extract the useful information in a matrix. One row per SNV, with the following columns:
#  SNV position, SNV mutation type, SNV class (G, S, HET), consequence (Non-coding transcript exon variant,
#  Upstream/downstream gene variant, Synonymous variant, Missense variant), samples. 

useful.info <- matrix(ncol=5, nrow=nrow(raw.info.no.re));
colnames(useful.info) <- c('SNV position', 'SNV mutation type', 
                           'SNV class', 'SNV consequence', 'Samples');

for (snv in 1:nrow(raw.info.no.re)){
  
  position <- substr(as.character(raw.info.no.re[snv,1]), 3, nchar(as.character(raw.info.no.re[snv,1]))-3);
  type <- substr(as.character(raw.info.no.re[snv,1]), 
                               nchar(as.character(raw.info.no.re[snv,1]))-2,
                               nchar(as.character(raw.info.no.re[snv,1])));
  class <- as.character(raw.info.no.re[snv,3]);
  consequence <- as.character(raw.info.no.re[snv,4]);
  sample <- as.character(raw.info.no.re[snv,2]);
  
  useful.info[snv,] <- rbind(position, type, class, consequence, sample);
  
}


## Convert the mutational types to the consensus nomenclature:
# C>A (=G>T), C>T (=G>A), C>G (=G>C), T>A (=A>T), T>C (=A>G), T>G (=A>C).

# Create an equivalence matrix (first column: consensus mut, second column: non-consensus mut)

equivalent.mut <- matrix(c('C>A', 'G>T', 
                           'C>T', 'G>A', 
                           'C>G', 'G>C', 
                           'T>A', 'A>T', 
                           'T>C', 'A>G', 
                           'T>G', 'A>C'),
                         nrow=6, ncol=2, byrow=TRUE);


# Sustitute the non-consensus mutations by the consensus mutations in useful.info.

useful.info.consensus <- useful.info;

for(mut in 1:nrow(equivalent.mut)){
  
  indices.to.change <- grep(equivalent.mut[mut,2], useful.info.consensus[,2]);
  useful.info.consensus[indices.to.change, 2] <- equivalent.mut[mut,1];
  
}



## 1. Mutational signatures. Create two plots: one showing the mutational signature only for
#  somatic SNVs and the other for the rest of SNVs.
#  Mutation types: C>A (=G>T), C>T (=G>A), C>G (=G>C), T>A (=A>T), T>C (=A>G), T>G (=A>C).

# Function which creates a mutational signature plot given a vector with the SNVs (snv.variants).
# The mutations in the vector have to be in the consensus form.
# e.g. c('C>T','C>T','T>A')
# main.title specifies the title in the plot.

mut.plot <- function(snv.variants, main.title){
  
  consensus.mut <- c('C>A', 'C>T', 'C>G', 'T>A', 'T>C', 'T>G');
  counts <- sapply(consensus.mut, function(x){length(grep(x, snv.variants))});
  barplot(counts, col=c('blue', 'red', 'black', 'grey', 'green', 'magenta'), ylim=c(0,20),
          main=main.title, cex.main = 0.9);
}


# Create two vectors: one containing all the SNVs and the other containing the somatic SNVs.

all.snvs <- useful.info.consensus[,2];
somatic.snvs <- useful.info.consensus[which(useful.info.consensus[,3] == 'S'),2];

# Create the two plots.

setwd('/Users/dmh/Desktop/trial_caveman/mit_processing/extra_stuff/mutational_signature/');

pdf('mutational_signatures_final.pdf', height=6, width=12);

par(mfrow=c(1,2));
mut.plot(all.snvs, 'Mutational signature for all the SNVs');
mut.plot(somatic.snvs, 'Mutational signature for the somatic SNVs');

dev.off();



## 2. Effect of SNVs.

library('ggplot2');

# Create a dataframe showing the information needed (value: % SNPs, sample:host/tumour, consequence). 
# Upstream/downstream SNVs are classified as 'Other'.

host.samples <- grep('host', useful.info.consensus[,5]);
tumour.samples <- grep('tumour', useful.info.consensus[,5]);
samples.data <- rep(NA, length(useful.info.consensus[,5]));
samples.data[host.samples] <- 'Host';
samples.data[tumour.samples] <- 'Tumour';

consequence.data <- gsub('Upstream, downstream gene variant', 'Other', useful.info.consensus[,4]);
value.data <- ifelse(samples.data=='Host', 100/length(host.samples), 100/length(tumour.samples));
data.effect <- data.frame(value=value.data, sample=samples.data, consequence=consequence.data);

# Order the data.

data.effect <- data.effect[order(data.effect[,3]),];

# Plot the results.

setwd('/Users/dmh/Desktop/trial_caveman/mit_processing/extra_stuff/VEP_plot/')

pdf('VEP_plot_final.pdf', height=6, width=6);

ggplot(data.effect, aes(x=sample, y=value, fill = consequence)) + geom_bar(stat ="identity") +
  xlab('Sample type') + ylab('% SNVs') +labs(fill='SNV consequence') +
  theme_bw();

dev.off();



## 3. Rainfall plot. Create a plot which represents genomic distance between two consecutive SNVs in
#  the Y-axis and the SNV order in the X-axis. Colour-code the SNPs by their consensus mutational type. 

# Calculate the genomic distance for the first SNV, taking into account that the mitochondrial genome 
# is circular: (length_mit_genome - position_last_SNV) + position_first_SNV.

distance.first.snv <- (16627 - as.numeric(useful.info.consensus[nrow(useful.info.consensus),1])) + 
  as.numeric(useful.info.consensus[1,1]);

# Calculate the genomic distance between the all SNVs.

genomic.distances <- c(distance.first.snv, diff(as.numeric(useful.info.consensus[,1])));

# Assign colours according to the mutational type to each SNV.

mut.types <- useful.info.consensus[,2];
colours.snv <- ifelse(mut.types=='C>A', 'blue', ifelse(mut.types=='C>T', 'red', 
                                                       ifelse(mut.types=='C>G', 'black', 
                                                              ifelse(mut.types=='T>A', 'grey',
                                                                     ifelse(mut.types=='T>C', 'green', 'magenta')))));

# Plot the results.

setwd('/Users/dmh/Desktop/trial_caveman/mit_processing/extra_stuff/rainfall_plot/');

pdf('rainfall_plot_final.pdf', height=6, width=8);

plot(1:length(genomic.distances), genomic.distances, col=colours.snv, pch=16, bty='l', cex=0.7,
     ylab='Genomic distance', xlab='', xaxt='n', cex.axis = 0.8, main='Rainfall plot');

legend('topleft', legend=c('C>A', 'C>T', 'C>G', 'T>A', 'T>C', 'T>G'), 
       col=c('blue', 'red', 'black', 'grey', 'green', 'magenta'),
       pch=16, ncol=3);

dev.off();

##### End of the script #####
