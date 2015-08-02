##################################################################################################
##                                  Daniel Elias Martin Herranz                                 ##
##                                 MPhil in Computational Biology                               ##
##                                    University of Cambridge                                   ##
##                                          August 2015                                         ##
##################################################################################################

##### Create BAF plots for unfiltered mitochondrial SNVs in Tasmanian devil samples #####

setwd('/Users/dmh/Desktop/trial_caveman/mit_processing/baf_plots/');


## Read argument from command line (number of sample to process). 

args <- commandArgs(TRUE);
sample.number <- args[1];


## Read the data from the VCF file. Keep only the lines with variants.

input <- readLines(con=paste0('/Users/dmh/Desktop/trial_caveman/mit_processing/flagged_files/', sample.number, '.vcf'));
input.variants <- input[grep('^#', input, invert=TRUE)];
all.info <- unlist(strsplit(input.variants, '\t')); 


## Extract the information for the variant position, variant flag and BAF.

info.variants <- matrix(nrow=length(input.variants), ncol=3);
colnames(info.variants) <- c('Position', 'Flag', 'BAF');

info.variants[,1] <- all.info[seq(from=2, to=length(all.info), by=11)];
info.variants[,2] <- all.info[seq(from=7, to=length(all.info), by=11)];
info.baf.all <- unlist(strsplit(all.info[seq(from=11, to=length(all.info), by=11)], ':'));
info.variants[,3] <- info.baf.all[seq(from=10, to=length(info.baf.all), by=10)];


## Create the BAF plots. 

# Those variants flagged as 'PASS' are depicted in blue. The rest of the variants are depicted in 'RED'.

pdf(paste0('results/', sample.number, '_BAF_plot.pdf'), height=6, width=12);

colour.variants <- ifelse(info.variants[,2] == 'PASS', 'blue', 'red');

plot(info.variants[,1], info.variants[,3], col='white', pch='.', xlim=c(1, 16627), ylim=c(0,1), 
     main=paste0('Sample ', sample.number), xlab='Position in mitochondrial genome', ylab='BAF', cex.axis=0.8);

text(as.numeric(info.variants[,1]), as.numeric(info.variants[,3]), labels=as.character(info.variants[,1]), col=colour.variants, cex=0.7);

legend('bottomleft', legend=c('Unflagged SNVs', 'Flagged SNVs (e.g. MQ, GI)'), col=c('blue', 'red'), pch=16, cex=0.5);

dev.off();

##### End of the script #####
