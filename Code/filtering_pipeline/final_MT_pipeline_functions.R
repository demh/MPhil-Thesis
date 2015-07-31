##################################################################################################
##                                  Daniel Elias Martin Herranz                                 ##
##                                 MPhil in Computational Biology                               ##
##                                    University of Cambridge                                   ##
##                                          August 2015                                         ##
##################################################################################################

##### Final pipeline for mitochondrial SNVs filtering and annotation #####

#### Functions ####

### 1. Read the VCF files and extract the variant information. Return a list with the variants information
#   for each smple. 

read.vcf <- function(){
  
  setwd('/Users/dmh/Desktop/trial_caveman/mit_processing/');
  
  variant.information <- vector("list", 40); 
  
  for (i in 1:length(variant.information)) {
    
    print(paste('Reading file', i, '...'));
    file <- readLines(con=paste0('flagged_files/', i, '.vcf'));
    file.variants <- file[grep('^#', file, invert=TRUE)]; # Delete comment lines, retain lines with variant information.
    variant.information[[i]] <- file.variants;
    rm(file); rm(file.variants);
  }
  
  # Return the results.
  
  return(variant.information);
  
}


### 2. Create a table (matrix) which displays interesting information from the VCF file, with variants as rows and samples as columns. 
#      For each variant present in a sample, specify the flags of the filtering and the BAF (PM: proportion of mutant allele). 
#      e.g. "PASS:9.5e-01".

# data.vcf: data as exported by read.vcf function (list).
# matched.info: matrix with two rows specifying the sample number and the sample ID, respectively.
#   e.g. Row_1: 1     2       3       4     5
#        Row_2: 133T  133H    140T    150T  150H

table.variants <- function(data.vcf, matched.info=NA){
  
  # Extract information for variant position, flags, substitution type (reference, alteration) and BAF.
  
  positions.per.sample <- vector('list', length(data.vcf));
  flags.per.sample <- vector('list', length(data.vcf));
  ref.per.sample <- vector('list', length(data.vcf));
  alt.per.sample <- vector('list', length(data.vcf));
  BAF.per.sample <- vector('list', length(data.vcf));
  
  for(i in 1:length(data.vcf)){
    
    print(paste('Extracting information for sample', i, '...'));
    separated.info <- unlist(strsplit(data.vcf[[i]], '\t'));
    positions.per.sample[[i]] <- separated.info[seq(from=2, to=length(separated.info), by=11)];
    flags.per.sample[[i]] <- separated.info[seq(from=7, to=length(separated.info), by=11)];
    ref.per.sample[[i]] <- separated.info[seq(from=4, to=length(separated.info), by=11)];
    alt.per.sample[[i]] <- separated.info[seq(from=5, to=length(separated.info), by=11)];
    BAF.unprocessed <- unlist(strsplit(separated.info[seq(from=11, to=length(separated.info), by=11)], ':'));
    BAF.per.sample[[i]] <- BAF.unprocessed[seq(from=10, to=length(BAF.unprocessed), by=10)];
  }
  
  # Initialize a matrix with all the variants as rows and all the samples as columns. 
  
  different.positions <- as.character(sort(as.numeric(unique(unlist(positions.per.sample)))));
  final.table <- matrix('*', nrow=(length(different.positions)), ncol=length(data.vcf));  # Asterisk represents absence of variant in that sample
  
  colnames(final.table) <- 1:40;
  
  if(is.matrix(matched.info)){ # Add tumour/host information.
    
    colnames(final.table) <- paste0(colnames(final.table), ':', matched.info[2,]);
    
  }
  
  # Obtain row names, including the systematic nomenclature for the substitution. 
  
  pos.ref.alt <- matrix(nrow=length(different.positions), ncol=3);
  pos.ref.alt[,1] <- different.positions;
  
  print('Processing substitution names ...');
  
  for(i in 1:length(different.positions)){
    
    for (j in 1:length(positions.per.sample)){
      
      if(different.positions[i] %in% positions.per.sample[[j]]) {
        
        pos.ref.alt[i,2] <- ref.per.sample[[j]][match(different.positions[i], positions.per.sample[[j]])];
        pos.ref.alt[i,3] <- alt.per.sample[[j]][match(different.positions[i], positions.per.sample[[j]])];
        break;
      }
    }
  }
  
  nomenclature <- paste0(pos.ref.alt[,1], pos.ref.alt[,2], '>', pos.ref.alt[,3]);
  rownames(final.table) <- nomenclature;
  
  # Fill in matrix with variant information for each sample (i.e. "flag:BAF")
  
  print('Filling in the matrix ...')
  
  for(i in 1:length(positions.per.sample)){
    
    indices.of.rows <- which(different.positions %in% positions.per.sample[[i]]);
    final.table[indices.of.rows,i] <- paste0(flags.per.sample[[i]], ':', BAF.per.sample[[i]]);
  }
  
  # Return the results.
  
  return(final.table);
  
}


### 3. Create a tab-separated file which contains chromosome and position of the variants called. This will be used by the alleleCounter.pl script
#      as the SNP file.

# data.table: data as exported by table.variants (matrix).
# chr: chromosome name.
# file.name: name of the SNP file.

create.snp.file <- function(data.table, chr='MT', file.name='SNP_file_MT.txt'){
  
  print(paste('Creating SNP file for chromosome', chr, '...'));
  
  positions.of.variants.unprocessed <- sapply(strsplit(rownames(data.table), '>'), '[[', 1);
  positions.of.variants <- as.numeric(sapply(positions.of.variants.unprocessed, function(x) {substr(x, 1, nchar(x)-1)}));
  chromosomes <- rep(chr, length(positions.of.variants));
  chr.and.pos <- cbind(chromosomes, positions.of.variants);
  write.table(chr.and.pos, file=file.name, sep='\t', quote=F, row.names=F, col.names=F);
  
}

### 4. Incorporate the information extracted with alleleCounter.pl into the unfiltered variant data. For those
#      variants already picked up by Caveman, check accross all samples whether there are some reads supporting 
#      the variant in a different sample. In this case, the BAF is calculated as the number of reads suporting 
#      the variant / number of total reads (as considered by alleleCounter.pl). The new variats are annotated
#      as AC:BAF.

# data.table: data as exported by remove.variants.with.flag, after filtering out for certain flags (e.g GI, RP, SE).

add.alleleCounter.info <- function(data.table){
  
  new.data.table <- data.table;  # Table to modify.
  
  for (sample in 1:ncol(data.table)){
    
    print(paste('Processing AC in sample', sample, '...'));
    
    # Read the alleleCounter data corresponding to the samples.
    
    AC.data <- read.table(file=paste0('/Users/dmh/Desktop/trial_caveman/mit_processing/allele_counter/alleleFrequencies_',
                                      sample, '_test.txt'), header=F, sep='\t');
    
    rownames(AC.data) <- AC.data[,2];
    AC.data <- AC.data[, 3:7];
    colnames(AC.data) <- c('A', 'C', 'G', 'T', 'Total');
    
    # Extract the read counts for those variants present in data.table.
    
    variants.positions.from.data.table <- as.numeric(substr(rownames(data.table), 1, nchar(rownames(data.table)) -3));
    read.counts.of.interest <- AC.data[as.numeric(rownames(AC.data)) %in% variants.positions.from.data.table,];
    
    # Select those positions which have not been called already by Caveman.
    
    positions.not.called <- as.numeric(which(data.table[,sample] == '*'));
    
    # Extract the expected substituted base from the positions not called and convert it to the
    # appropiate index (column in AC.data).
    
    selected.variant.names <- as.character(names(data.table[positions.not.called, sample]));
    subs.bases <- substr(selected.variant.names, nchar(selected.variant.names), nchar(selected.variant.names));
    subs.bases.indices <- as.numeric(factor(subs.bases, levels=c('A', 'C', 'G', 'T')));
    
    # Extract the number of reads supporting the 'mutant' variants and the total number of reads, for 
    # the positions not called by Caveman. Calculate the BAF for those positions.
    
    reads.mutant <- read.counts.of.interest[cbind(positions.not.called, subs.bases.indices)];
    reads.total <- read.counts.of.interest[cbind(positions.not.called, 5)];
    BAF <- reads.mutant/reads.total;
    
    # Fill in the new.data.table, for those cases where BAF > 0.
    
    BAF.formatted <- format(BAF, scientific=T, digits=2);
    new.data.table[positions.not.called[which(BAF > 0)],sample] <- paste0('AC:', BAF.formatted[which(BAF > 0)]);
  }
  
  return(new.data.table);
  
}


### 5. Flag those variants, present in tumours, that are probably host contamination/HC (BAF < threshold, present in matched host).

# data.table: data as exported by table.variants (matrix). It needs to include the tumour / host information in the col.names.
# thr.hc.matched: threshold value for host contamination when the matched host is present (BAF in tumour < thr.hc.matched).
# thr.hc.unmatched: threshold value for host contamination when the matched host is NOT present (BAF in tumour < thr.hc.unmatched).

filter.hc <- function(data.table, thr.hc.matched, thr.hc.unmatched){
  
  # Classify the tumours according to the presence or absence of matched host.
  
  all.tumours.info <- colnames(data.table)[grep('T', colnames(data.table))];
  all.tumours.ids <- sapply(strsplit(all.tumours.info, ':'), '[[', 2);
  
  tumours.matched <- c();
  tumours.unmatched <- c();
  
  for (tumour in all.tumours.ids){
    
    if(length(grep(gsub('T', 'H', tumour), colnames(data.table))) == 1){
      
      tumours.matched <- c(tumours.matched, tumour);
      
    }else{
      
      tumours.unmatched <- c(tumours.unmatched, tumour);
    }
  }
  
  print(paste('We found', length(tumours.matched), 'matched tumours and', length(tumours.unmatched), 'unmatched tumours.'));
  
  
  # Flagging the tumours with matched host.
  
  for(i in 1:length(tumours.matched)){
    
    print(paste('Flagging matched tumour', tumours.matched[i], '...'));
    
    index.tumour <- grep(tumours.matched[i], colnames(data.table));
    host.name <- gsub('T', 'H', tumours.matched[i]);
    index.host <- grep(host.name, colnames(data.table));
    
    BAF.values.tumour <- unlist(ifelse(data.table[,index.tumour] == '*', NA, strsplit(data.table[,index.tumour], ':')));
    BAF.values.tumour <- as.numeric(BAF.values.tumour[-grep('[A-Z;]+', BAF.values.tumour)]);
    BAF.condition <- ifelse(BAF.values.tumour < thr.hc.matched, TRUE, FALSE);
    BAF.condition[is.na(BAF.condition)] <- FALSE;
    
    conditions.matched <- (data.table[,index.host] != '*') & (BAF.condition);
    BAF.matched <- sapply(strsplit(data.table[conditions.matched,index.tumour], ':'), '[[', 2);
    data.table[conditions.matched,index.tumour] <- paste0('HC:', BAF.matched);
    
  }
  
  # Flagging the tumours without matched host.
  
  for(i in 1:length(tumours.unmatched)){
    
    print(paste('Flagging unmatched tumour', tumours.unmatched[i], '...'));
    
    index.tumour <- grep(tumours.unmatched[i], colnames(data.table));
    index.hosts <- grep('H', colnames(data.table));
    
    BAF.values.tumour <- unlist(ifelse(data.table[,index.tumour] == '*', NA, strsplit(data.table[,index.tumour], ':')));
    BAF.values.tumour <- as.numeric(BAF.values.tumour[-grep('[A-Z;]+', BAF.values.tumour)]);
    indexes.possible.HC <- ifelse(BAF.values.tumour < thr.hc.unmatched, TRUE, FALSE);
    indexes.possible.HC <- which(indexes.possible.HC == TRUE);
    
    variants.accross.hosts <- data.table[indexes.possible.HC, index.hosts] != '*';
    
    if(length(indexes.possible.HC) == 1){
      
      presence.in.any.host <- any(variants.accross.hosts);
      indexes.HC <- indexes.possible.HC[presence.in.any.host];
      
      if(length(indexes.HC) > 0){
        
        BAF.unmatched <- sapply(strsplit(data.table[indexes.HC,index.tumour], ':'), '[[', 2);
        data.table[indexes.HC,index.tumour] <- paste0('HC:', BAF.unmatched);
        
      }
    }
    
    if(length(indexes.possible.HC) > 1){
      
      presence.in.any.host <- apply(variants.accross.hosts, 1, any);
      indexes.HC <- indexes.possible.HC[presence.in.any.host];
      
      if(length(indexes.HC) > 0){
        
        BAF.unmatched <- sapply(strsplit(data.table[indexes.HC,index.tumour], ':'), '[[', 2);
        data.table[indexes.HC,index.tumour] <- paste0('HC:', BAF.unmatched);
        
      }
    }
  }
  
  # Return results.
  
  return(data.table);
  
}


### 6. Filter out those variants that are flagged with specific filters. Simplify the results (exclude variants not present in any sample
#      after filtering). Those variants that are flagged with more than one filter are also discarded.

# data.table: data as exported by table.variants or filter.hc (matrix).
# filters: vector specifying those flags that we want TO KEEP e.g. "c('PASS', 'MQ')".
# binary: convert to binary (1:present, 0:absent) output. Possible values: TRUE/FALSE.

remove.variants.with.flag <- function(data.table, filters=c('PASS', 'MQ'), binary=FALSE){
  
  # Extract the indices of the variants to keep.
  
  print('Extracting indices of variants to keep ...');
  
  indices.to.keep <- c();
  
  for (i in 1:length(filters)){
    
    # Discard those variants flagged with more than one filter.
    
    possible.to.keep <- grep(filters[i], data.table);
    discard <- which(nchar(sapply(strsplit(data.table[possible.to.keep], ':'), '[[', 1)) > nchar(filters[i]));
    
    if(length(discard) > 0){
      
      keep <- possible.to.keep[-discard];
      
    } else{
      
      keep <- possible.to.keep;
    }
    
    indices.to.keep <- c(indices.to.keep, keep);
  }
  
  indices.to.keep <- unique(sort(indices.to.keep));
  
  
  # Create a matrix with the same dimensions as data.table, which has TRUE for the indices.to.keep.
  
  boolean.filters <- matrix(nrow=nrow(data.table), ncol=ncol(data.table));
  boolean.filters[indices.to.keep] <- TRUE;
  boolean.filters[-indices.to.keep] <- FALSE;
  
  # Extract the variants selected. If not selected, output NA.
  
  print('Constructing matrix with variants to keep ...');
  
  for(i in 1:ncol(data.table)){
    
    data.table[,i] <- ifelse(boolean.filters[,i], data.table[,i], NA);
    
  }
  
  # Simplify the matrix (delete those variants not present in any sample).
  
  print('Performing matrix simplification ...');
  
  data.table <-  data.table[rowSums(is.na(data.table))!=ncol(data.table),];
  
  # Change the format of the output matrix to binary, in case it is requested. Return the results.
  
  if(binary){
    
    print('Binarizing matrix ...')
    data.table <- ifelse(is.na(data.table), 0, 1);
    return(data.table);
    
  } else{
    
    data.table[is.na(data.table)] <- '*'
    return(data.table);
  }
}


### 7. Flag those variants that are candidates for heteroplasmy / subclonality (HET). 
#      - HET is picked up in hosts when the variant is present only in one sample and BAF < host.thr
#      - HET is picked up in tumours when the variant is only present in tumours. In this case, these candidates 
#        need to be validated manually with BAF plots (determine level of host contamination).

# filtered.data.table: data filtered for unwanted flags (e.g. GI), including AC flagging and filtered for HC.
# host.thr: threshold for flagging a variant in a host as HET (BAF in host < host.thr). Default: 0.8.

filter.het <- function(filtered.data.table, host.thr=0.8){
  
  ## Initialize table to output.
  
  output.table <- filtered.data.table;
  
  
  ## Identify possible HET in hosts.
  
  # Extract indexes for host samples.
  
  hosts.indexes <- grep('H', colnames(filtered.data.table));
  
  # Find those variants which are present only in one sample and the samples that contain them.
  
  unique.variants <- which(apply(filtered.data.table, 1, function(x) sum(x == '*')) == ncol(filtered.data.table) - 1);
  samples.with.unique.variants <- sapply(unique.variants, function(x) {which(filtered.data.table[x,] != '*')});
  
  # Select those samples that are hosts and the corresponding unique variants.
  
  hosts.with.unique.variants <- samples.with.unique.variants[samples.with.unique.variants %in% hosts.indexes];
  unique.variants.in.hosts <- unique.variants[samples.with.unique.variants %in% hosts.indexes];
  
  # Extract BAF values for these unique variants in hosts and check whether they are < host.thr.
  
  unique.variants.in.hosts.info <- filtered.data.table[cbind(unique.variants.in.hosts, hosts.with.unique.variants)];
  BAF.unique.variants.hosts <- as.numeric(sapply(strsplit(unique.variants.in.hosts.info, ':'), '[[', 2));
  final.unique.hosts.selected.variants <- unique.variants.in.hosts[BAF.unique.variants.hosts < host.thr];
  final.selected.hosts.with.unique.variants <- hosts.with.unique.variants[BAF.unique.variants.hosts < host.thr];
  final.BAF.selected.unique.variants.hosts <- format(BAF.unique.variants.hosts[BAF.unique.variants.hosts < host.thr], scientific=T, digits=2);
  
  # Flag the final selected variants in hosts as possible HET.
  
  if(length(cbind(final.unique.hosts.selected.variants,final.selected.hosts.with.unique.variants)) > 0){
    
    output.table[cbind(final.unique.hosts.selected.variants,final.selected.hosts.with.unique.variants)] <- paste0('HET:', final.BAF.selected.unique.variants.hosts);
  }
  
  
  ## Identify possible HET in tumours.
  
  # Extract indexes for tumour samples.
  
  tumours.indexes <- grep('T', colnames(filtered.data.table));
  
  # For each variant, extract the samples in which it is present.
  
  samples.per.variant <- apply(filtered.data.table, 1, function(x) {which(x != '*')});
  
  # For each variant, check whether all the samples in which it is present are tumours.
  
  check.variants.for.tumours <- sapply(samples.per.variant, function(x) {length(grep('T',names(x))) == length(x)});
  
  # Flag the selected variants as possible HET.
  
  indexes.variants.unique <- which(check.variants.for.tumours == TRUE);
  indexes.variants.selected <- rep(indexes.variants.unique, as.numeric(sapply(samples.per.variant[check.variants.for.tumours], length)));
  indexes.samples.selected <- as.numeric(unlist(samples.per.variant[check.variants.for.tumours]));
  BAF.variants.in.tumours <- as.numeric(sapply(strsplit(output.table[cbind(indexes.variants.selected,indexes.samples.selected)], ':'), '[[', 2));
  final.BAF.variants.in.tumour <- format(BAF.variants.in.tumours, scientific=T, digits=2);
  
  if(length(cbind(indexes.variants.selected,indexes.samples.selected)) > 0){ # If there are any matches.
    output.table[cbind(indexes.variants.selected,indexes.samples.selected)] <- paste0('HET:', final.BAF.variants.in.tumour);
  }
  
  
  ## Return the results.
  
  return(output.table);
  
}


### 8. Flag with DISC (discard) those variants in hosts, with low BAF values and which are not HET variants. 

# data.table.het: data which contains HET flags for those variants that we believe are heteroplasmy candidates.
# discard.thr: threshold for flagging a variant as DISC (BAF in variant < discard.thr). It has to be adjusted 
#              taking into account the level of contamination expected in the samples. Default: 0.8.

filter.disc <- function(data.table.het, discard.thr=0.8){
  
  ## Initialize table to be returned.
  
  output.table <- data.table.het;
  
  ## Obtain the indexes of the host samples.
  
  indexes.hosts <- grep('H', colnames(data.table.het));
  
  ## Process the variants for each one of the hosts.
  
  for(host in indexes.hosts){
    
    # Extract information from variants found in host.
    
    indexes.variants.in.host <- which(data.table.het[,host] != '*');
    all.info.variants.in.host <- strsplit(data.table.het[indexes.variants.in.host,host], ':');
    flags.variants.in.host <- sapply(all.info.variants.in.host, '[[', 1);
    BAF.variants.in.host <- as.numeric(sapply(all.info.variants.in.host, '[[', 2));
    
    # Find those variants which are not flagged with HET and that have a BAF value < discard.thr.
    
    selected.variants.boolean <- (flags.variants.in.host != 'HET') & (BAF.variants.in.host < discard.thr);
    coordinates.selected.variants <- cbind(indexes.variants.in.host[selected.variants.boolean], rep(host, length(indexes.variants.in.host[selected.variants.boolean])));
    
    # Flag with DISC the variants selected.
    
    output.table[coordinates.selected.variants] <- gsub('^[A-Z]+', 'DISC', output.table[coordinates.selected.variants]);
  }
  
  ## Return the results.
  
  return(output.table);
}


### 9. Perform hierarchical clustering in a binary matrix. Export result as a heatmap in pdf.

# binary.matrix: same format as output from remove.variants.with.flag, using the binary=TRUE argument.
# plot.name: name of the plot.

variants.clustering <- function(binary.matrix, plot.name){
  
  library(cluster);
  
  class(binary.matrix) <- "numeric";
  palette <- colorRampPalette(c('#f0f3ff','#0033BB'))(256);
  
  pdf(plot.name, width=13, height=10);
  par(mfrow=c(1,1));
  heatmap(t(binary.matrix), Colv=F, cexRow=0.7, cexCol=0.7, scale='none',
          xlab="Variant Positions", ylab="Samples", col=palette);
  dev.off()
  
}


### 10. Introduce the changes in the final data table noticed after manual curation of the variants.
#      Specify the variants to be added and those to be excluded. The variants to be added are flagged as ADD:NA
#      and the ones to be deleted as DEL:NA. The function returns the final data table with the new flags. This
#      data table can then be processed with the remove.variants.with.flag function.


# final.table: filtered data to which we add the manual curation information. 
# include.matrix: matrix, with two rows, specifying the variant position to be added
#                 and the sample(s) (column number(s) in data table) that must contain it, respectively. 
#                 The variants listed here must be present in other samples. 
# exclude,matrix: matrix, with two rows, specifying the variant position to be deleted
#                 and the sample(s) (column number(s) in data table) that use to contain it, respectively.
#
#   e.g.  for matrices  Row_1 (variant position):      1786       10678          15678
#                       Row_2 (samples):               1,5,30      10             3,22    


manual.curation <- function(final.table, include.matrix, exclude.matrix){
  
  # Initialize the table to return.
  
  return.table <- final.table;
  
  # Flag variants to add.
  
  for (i in 1:ncol(include.matrix)){
    
    variant <- include.matrix[1,i];
    index.variant <- grep(variant, rownames(return.table));
    samples <- as.numeric(unlist(strsplit(include.matrix[2,i], ',')));
    coordinates <- cbind(rep(index.variant, length(samples)), samples); 
    return.table[coordinates] <- 'ADD:NA'; 
    
  }
  
  # Flag variants to remove.
  
  for (i in 1:ncol(exclude.matrix)){
    
    variant <- exclude.matrix[1,i];
    index.variant <- grep(variant, rownames(return.table));
    samples <- as.numeric(unlist(strsplit(exclude.matrix[2,i], ',')));
    coordinates <- cbind(rep(index.variant, length(samples)), samples);
    return.table[coordinates] <- 'DEL:NA'; 
    
  }
  
  # Return the modified table.
  
  return(return.table);
  
}


### 11. Create a text file with all the coordinates for the filtered variants, to be used
#      as the input in the Variant Effect Predictor of Ensembl.

# data.table: table containing the filtered variants (e.g. using the output from remove.variants.with.flag).
# file.path: string containing the path to save the file. Default: '/Users/dmh/Desktop/trial_caveman/mit_processing/VEP_input.txt'.
# chr: chromosome name. Default: 'MT'.

create.file.VEP <- function(data.table, file.path='/Users/dmh/Desktop/trial_caveman/mit_processing/VEP_input.txt', chr='MT'){
  
  # Extract the filtered variants.
  
  variant.names <- rownames(data.table);
  
  # Construct the matrix to output.
  
  final.matrix <- matrix(ncol=5, nrow=length(variant.names));
  
  final.matrix[,1] <- chr; # Chromosome name
  final.matrix[,2] <- substr(variant.names, 1, nchar(variant.names)-3); # First coordinate
  final.matrix[,3] <- final.matrix[,2];  # Second coordinate
  final.matrix[,4] <- gsub('>', '/', substr(variant.names, nchar(variant.names)-2, nchar(variant.names))); # Substitution type.
  final.matrix[,5] <- 1; # Length of substitution.
  
  # Export matrix as text file separated by white spaces.
  
  write.table(final.matrix, file=file.path, quote=F, row.names=F, col.names=F, sep=' ');
  
}


##### End of the script #####

