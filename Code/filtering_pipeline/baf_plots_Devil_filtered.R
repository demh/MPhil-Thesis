##################################################################################################
##                                  Daniel Elias Martin Herranz                                 ##
##                                 MPhil in Computational Biology                               ##
##                                    University of Cambridge                                   ##
##                                          August 2015                                         ##
##################################################################################################

##### Create BAF plots for final filtered mitochondrial SNVs in Tasmanian devil samples #####

setwd('/Users/dmh/Desktop/trial_caveman/mit_processing/baf_plots_filtered/');

## Read the information from the variants.

raw.data <- read.csv('filtered_variants.csv');
class.variants <- read.csv('list_of_variants.csv');

## Create a matrix which contains variant name (e.g.264G>A) and class (e.g.G).

variant.name.and.class <- matrix(nrow=nrow(class.variants), ncol=2);
variant.name.and.class[,1] <- substr(as.character(class.variants[,1]),
                                     3, nchar(as.character(class.variants[,1])));
variant.name.and.class[,2] <- as.character(class.variants[,3]);

## Format the information. Create a list of matrices, one matrix per sample.
#  each matrix containing the following columns: variant position, classification
#  of the variant (G, G/RE, S, HET), BAF.

final.list.per.sample <- list(); # Initialize list to contain all the information.

for(i in 2:ncol(raw.data)){
  
  sample <- i-1;
  sample.variant.name <- as.character(raw.data[which(raw.data[,i]!= '*'), 1]);
  sample.variant.position <- substr(sample.variant.name, 1, nchar(sample.variant.name)-3);
  sample.BAFs <- as.numeric(sapply(strsplit(as.character(raw.data[which(raw.data[,i]!= '*'), i]), 
                                            ':'), '[[', 2));
  indices.of.variants <- sapply(1:length(sample.variant.name), function(x) {
    grep(sample.variant.name[x], variant.name.and.class[,1])});
  sample.class <- variant.name.and.class[indices.of.variants,2];
  
  final.list.per.sample[[sample]] <- cbind(sample.variant.position, sample.class, sample.BAFs);
   
}


## For those variants that we added manually and which do not have a BAF value (ADD:NA), 
#  we add a BAF value based on the statistics from IGV.

# Create a matrix containing all the BAF values to be added. Each row has the following
# information: sample, variant position, BAF.

add.BAF.info <- rbind(c(23, 12864, 0.74),
                      c(9, 11788, 0.83),
                      c(6, 621, 0.68),
                      c(6, 12035, 0.60),
                      c(6, 13650, 0.70),
                      c(6, 16322, 0.73),
                      c(6, 12244, 0.77),
                      c(6, 13613, 0.70),
                      c(20, 621, 1.00));

# Add the BAF values.

for(i in 1:nrow(add.BAF.info)){
  
  sample <- add.BAF.info[i,1]; 
  variant <- add.BAF.info[i,2];
  BAF <- add.BAF.info[i,3];
  row.of.variant <- grep(as.character(variant), final.list.per.sample[[sample]][,1]);
  final.list.per.sample[[sample]][row.of.variant,3] <- BAF;
  
}


## Create the filtered BAF plots for the 40 samples. The SNVs will be colour-coded according
#  to the classification (G:blue, G/RE: green, S:red, HET: black).

for (i in 1:length(final.list.per.sample)){
  
  print(paste('Plotting sample', i, '...'));
  
  variants <- as.numeric(final.list.per.sample[[i]][,1]);
  classes <- final.list.per.sample[[i]][,2];
  BAFs <- as.numeric(final.list.per.sample[[i]][,3]);
  
  pdf(paste0('results/', i, '_BAF_plot_filtered.pdf'), height=6, width=12);
  
  colour.variants <- ifelse(classes == 'G', 'blue', ifelse(classes == 'G/RE', 'green', 
                                                           ifelse(classes == 'S', 'red', 'black')));
  
  plot(variants, BAFs, col='white', pch='.', xlim=c(1, 16627), ylim=c(0,1), 
       main=paste0('Sample ', i, ': filtered'), xlab='Position in mitochondrial genome', ylab='BAF', cex.axis=0.8);
  
  text(variants, BAFs, labels=as.character(variants), col=colour.variants, cex=0.7);
  
  legend('bottomleft', legend=c('Germline SNV', 'Germline SNV / Reference error', 'Somatic SNV', 'Heteroplasmic SNV'), 
         col=c('blue', 'green', 'red', 'black'), pch=16, cex=0.8);
  
  dev.off();
}

##### End of the script #####
