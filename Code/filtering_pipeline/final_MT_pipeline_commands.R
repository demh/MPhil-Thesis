##################################################################################################
##                                  Daniel Elias Martin Herranz                                 ##
##                                 MPhil in Computational Biology                               ##
##                                    University of Cambridge                                   ##
##                                          August 2015                                         ##
##################################################################################################

##### Final pipeline for mitochondrial SNVs filtering and annotation #####

#### Commands for the final results ####

source('/Users/dmh/Desktop/trial_caveman/mit_processing/final_MT_pipeline_functions.R');
setwd('/Users/dmh/Desktop/trial_caveman/mit_processing/');

# Read VCF files.

MT.variants <- read.vcf();

# Specify the tumour / host naming.

tumour.host.MT <- matrix(ncol=40, nrow=2);
tumour.host.MT[1,] <- 1:40;
tumour.host.MT[2,] <- c('107T', '107H',
                        '139T', '139H',
                        '140T', '140H',
                        '141T', 
                        '142T', '142H',
                        '143T', 
                        '144T', '144H',
                        '145T', '145H',
                        '146T', '146H',
                        '147T', '147H',
                        '148T', '148H',
                        '149T', '149H',
                        '150.2H', '150H', # Rename 150T as 150.2H
                        '151.2H', '151H', # Rename 151T as 151.2H
                        '152T', '152H',
                        '153T', '153H',
                        '154T', 
                        '155T', '155H',
                        '156T',
                        '157T', '157H',
                        '158T', '158H',
                        '159T', '159H');

# Construct table with variants and BAF.

MT.unfiltered.table <- table.variants(MT.variants, tumour.host.MT);

# Exclude variants flagged with flags other than PASS or MQ. NOT BINARY.

final.MT.only.MQ.PASS <- remove.variants.with.flag(MT.unfiltered.table, filters=c('PASS', 'MQ'), binary=F);

# Exclude variants flagged with flags other than PASS or MQ. BINARY.

final.MT.only.MQ.PASS.binary <- remove.variants.with.flag(MT.unfiltered.table, filters=c('PASS', 'MQ'), binary=T);

# Construct hierarchical clustering without taking into account possible host contamination.

variants.clustering(final.MT.only.MQ.PASS.binary, 'clustering_PASS_MQ_contains_HC.pdf');

# Flag variants which are host contamination in tumour samples.

final.MT.hc.table <- filter.hc(final.MT.only.MQ.PASS, 0.9, 0.5);

# Include results from the alleleCounter.pl script.

final.MT.with.AC <- add.alleleCounter.info(final.MT.hc.table);

# Keep variants flagged as PASS, MQ or AC (remove host contamination).

final.MT.filtered <- remove.variants.with.flag(final.MT.with.AC, filters=c('PASS', 'MQ', 'AC'), binary=F);

# Flag those variants which are putative HET variants.

final.MT.flagged.HET <- filter.het(final.MT.filtered, host.thr=0.8);

# Extract putative HET variants.

extract.HET.info <- apply(final.MT.flagged.HET, 1, function(x){x[grep('HET',x)]});
putative.HET.variants <- extract.HET.info[sapply(extract.HET.info, function(x){length(x) > 0})];

# Export final.MT.flagged.HET for manual checking of some of the variants.

# write.csv(final.MT.flagged.HET, file='/Users/dmh/Desktop/trial_caveman/mit_processing/all_variant_information_for_checking_second_version.csv');

# Flag with DISC (discard) those variants in hosts, with low BAF values and which are not HET variants. 
# In the case of variants present in tumours with low BAF values which are not HET variants (they are PASS or MQ),
# we keep them.

final.MT.flagged.DISC <- filter.disc(final.MT.flagged.HET, discard.thr=0.8);

# Keep those variants flagged as PASS, MQ or HET (remove DISC and AC). NOT BINARY.

final.MT.only.PASS.MQ.HET <- remove.variants.with.flag(final.MT.flagged.DISC, filters=c('PASS', 'MQ', 'HET'), binary=F);

# Keep those variants flagged as PASS, MQ or HET (remove DISC and AC). BINARY.

final.MT.only.PASS.MQ.HET.binary <- remove.variants.with.flag(final.MT.flagged.DISC, filters=c('PASS', 'MQ', 'HET'), binary=T);

# Create plot with hierarchical clustering.

variants.clustering(final.MT.only.PASS.MQ.HET.binary, 'changes_clustering.pdf');

# Create two matrices (include and exclude) which contains changes realted with the manual curation of the variants.

include.variants <- cbind(c('12864', '23'),
                          c('11788', '9'),
                          c('621', '6,20'),
                          c('12035', '6'),
                          c('13650', '6'),
                          c('16322', '6'),
                          c('12244', '6'),
                          c('13613', '6'));

exclude.variants <- cbind(c('867', '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40'),
                          c('15604', '3'),
                          c('254', '8,17'), 
                          c('6569', '7'),
                          c('2892', '5,17'),
                          c('2080', '35,19'),
                          c('15188', '35'),
                          c('1828', '39'),
                          c('4316', '32'),   # exclude?
                          c('8901', '27'),
                          c('8926', '27'),
                          c('13973', '15'),  # exclude?
                          c('8579', '15'),
                          c('3633', '15'),
                          c('5317', '29'),
                          c('4151', '29'),
                          c('4971', '29'));


# Flag the table with the information from manual curation.

final.curation.flagged <- manual.curation(final.MT.only.PASS.MQ.HET, include.variants, exclude.variants);

# Filter our table for the curated information (NOT BINARY).

final.after.curation.not.bin <- remove.variants.with.flag(final.curation.flagged, filters=c('PASS', 'MQ', 'HET', 'ADD'), binary=FALSE);

# Filter our table for the curated information (BINARY).

final.after.curation.bin <- remove.variants.with.flag(final.curation.flagged, filters=c('PASS', 'MQ', 'HET', 'ADD'), binary=TRUE);

# Create final hierarchical clustering.

variants.clustering(final.after.curation.bin , 'after_curation_clustering.pdf');

# Create text file with selected variants for using in the variant effect predictor and other extra analysis.

create.file.VEP(final.after.curation.not.bin, file.path='/Users/dmh/Desktop/trial_caveman/mit_processing/extra_stuff/mutational_signature/selected_variants.txt', chr='MT');


##### End of the script #####
