##################################################################################################
##                                  Daniel Elias Martin Herranz                                 ##
##                                 MPhil in Computational Biology                               ##
##                                    University of Cambridge                                   ##
##                                          August 2015                                         ##
##################################################################################################

##### Find possible driver genes contained in the amplified region of chromosome 1 (double-minute 
# chromosomes) of sample 88T #####

setwd('/Users/dmh/Desktop/Software/cn_calling/our_scripts/');


### Read all the data that we need.

data = read.table("88T.10kb_all.txt",sep="\t",header=F,row.names=NULL)
names(data) = c("chr","startpos","endpos","no.reads","bases.covered","segment.length","fractional.coverage")
data = data[order(data$chr,data$startpos),]
contig.lengths = read.table("chr_lengths.txt",sep="\t",header=T) 
contig.real.lengths = read.table('contig_and_length.txt', sep='\t', header=F)


### Create a function to extract the data from one chromosome, including the real genomic position.

# all.data: dataframe containing all data, as read previously.
# lengths: dataframe containing the names of the contigs (col 1) and their real lengths (col 2).
# chr: chromosome to extract the data from ('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'ChrU', 'Chrx').

extract.chr <- function(all.data, lengths, chr){
  
  ## Extract rows which contain data from our chromosome of choice.
  
  indices.chr.data <- which(substr(all.data$chr, 1, 4) == chr);
  chr.data <- all.data[indices.chr.data,];
  
  ## Calculate the real genomic coordinates for each window.
  
  # Create vectors which contain the numbers to add to different coordinates in contig.
  
  lengths.for.chr <- lengths[grep(chr, lengths[,1]),];
  
  cum.sum.contigs <- cumsum(as.numeric(lengths.for.chr[,2]));
  add.to.all.no.end.per.contig <- c(0, cum.sum.contigs[-length(cum.sum.contigs)]) + 1;
  add.to.end.per.contig <- cum.sum.contigs;
  
  # Transform to genomic coordinates.
  
  genomic.coord <- chr.data$startpos;
  contigs.numbers <- as.integer(unlist(strsplit(as.character(chr.data$chr), '_'))[seq(from=3, to=length(unlist(strsplit(as.character(chr.data$chr), '_'))), by=3)]);
  
  sum.to.each.position <- rle(NA); # Construct vector to sum to each position.
  sum.to.each.position$lengths <- rle(contigs.numbers)$lengths;     
  sum.to.each.position$values <-  add.to.all.no.end.per.contig; 
  sum.to.each.position <- inverse.rle(sum.to.each.position);
  
  genomic.coord <- genomic.coord + sum.to.each.position # Inside contigs.
  
  # Append to our data per chromosome.
  
  chr.data$genome.pos <- genomic.coord;
  
  
  ## Return data for chromosome.
  
  return(chr.data);
  
}


### Create a function to find the intersection between the genes in our amplified contigs and the COSMIC database.

# contigs.of.interest: specify the contigs we are interested in as a vector. e.g. c('Chr1_supercontig_000000182', 'Chr1_supercontig_000000183').
# synonyms: include gene name synonyms in our search in COSMIC database (TRUE/FALSE).
# orthologs: for those genes that do not have an external gene name in the Tasmanian Devil, try to find its external gene name through the human
#            ortholog (TRUE/FALSE).

overlap.with.cosmic <- function(contigs.of.interest, synonyms=TRUE, orthologs=TRUE){
  
  library("biomaRt");
  
  ## Extract gene names of genes in selected contigs.
  
  # Connect to Tasmanian Devil data in ENSEMBL database.
  
  ensembl <- useMart("ensembl",dataset="sharrisii_gene_ensembl");
  
  # Convert from our IDs to ENSEMBL contig IDs.
  
  equivalence.contigs <- read.table("/Users/dmh/Desktop/utilities/scaffold_localID2acc_for_7_1",sep="\t",header=F,row.names=NULL);
  contigs.selected <- as.character(equivalence.contigs[which(as.character(equivalence.contigs[,1]) %in% contigs.of.interest),2]); # Ensembl IDs.
  
  # Process taking into account human orthologs or not.
  
  if(orthologs){
    
    # Extract all gene names in the contigs selected, including information regarding the human orthologs.
    
    filters <- listFilters(ensembl)[1,1]; # Filter by contig name
    attributes <- list('external_gene_name', 'ensembl_gene_id', 'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_orthology_confidence'); # Output
    pre.gene.names.extracted <- getBM(attributes=attributes, filters=filters, values=contigs.selected, mart=ensembl);
    
    # Extract external gene names from the human dataset, for the human orthologs.
    
    ensembl.human <- useMart("ensembl",dataset="hsapiens_gene_ensembl");
    
    human.orthologs <- pre.gene.names.extracted[,3];
    filters.human <- 'ensembl_gene_id'; # Filter by Ensembl Gene ID from human.
    attributes.human <- list('external_gene_name', 'ensembl_gene_id');
    genes.orth <- getBM(attributes=attributes.human, filters=filters.human, values=human.orthologs, mart=ensembl.human);
    
    # Rename the columns for both dataframes.
    
    names(pre.gene.names.extracted) <- c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_orth_human", "hsapiens_homolog_orthology_confidence");
    names(genes.orth) <- c("external_gene_name", "ensembl_gene_id_orth_human");
    
    # Merge both dataframes, by human ortholog Ensembl ID.
    
    all.information.genes <- merge(pre.gene.names.extracted, genes.orth, by='ensembl_gene_id_orth_human', all=TRUE, sort=FALSE);
    
    # Reorder the information in all.information.genes. There will be 5 columns: devil_ensembl_id, devil_external_name, human_orth_ensembl_ids, human_orth_external_names 
    # and orth_confidence.
    
    col1 <- unique(all.information.genes[,3]);
    col2 <- c();
    col3 <- c();
    col4 <- c();
    col5 <- c();
    
    for (gene in 1:length(col1)){
      
      indices <- which(all.information.genes[,3] == col1[gene]);
      
      col2 <- c(col2, all.information.genes[indices[1], 2]);
      col3 <- c(col3, paste(all.information.genes[indices, 1], collapse=';'));
      col4 <- c(col4, paste(all.information.genes[indices, 5], collapse=';'));
      col5 <- c(col5, paste(all.information.genes[indices, 4], collapse=';'));
      
    }
    
    
    organized.all <- data.frame(devil_ensembl_id=col1, devil_external_name=col2, human_orth_ensembl_ids=col3, human_orth_external_names=col4,
                                orth_confidence=col5);
    
    # Create a data.frame with devil_ensembl_id as the first column and external_name as the second column. If devil_external_name is present, then
    # we take this as the external_name. Otherwise, we take the human_orth_external_names. The third column specifies where the external name
    # comes from (DEVIL, 0, 1).
    
    col1 <- as.character(organized.all[,1]);
    indices.with.devil <- which(organized.all[,2] != '');
    indices.with.orth <- which(organized.all[,2] == '');
    col2 <- rep(NA, nrow(organized.all));
    col2[indices.with.devil] <- as.character(organized.all[indices.with.devil,2]);
    col2[indices.with.orth] <- as.character(organized.all[indices.with.orth, 4]);
    col3 <- rep(NA, nrow(organized.all));
    col3[indices.with.devil] <- 'DEVIL';
    col3[indices.with.orth] <- as.character(organized.all[indices.with.orth, 5]);
    
    selected.names <- data.frame(devil_ensembl_id=col1, external_name=col2, origin=col3);
    
    # Extract gene names to check in COSMIC.
    
    final.names <- unlist(strsplit(as.character(selected.names[,2]), ';'));
    final.names <- final.names[which(final.names != 'NA')];
    
  }else{
    
    # Extract all gene names in the contigs selected.
    
    filters <- listFilters(ensembl)[1,1]; # Filter by contig name
    attributes <- list('external_gene_name', 'ensembl_gene_id'); # Output
    
    gene.names.extracted <- getBM(attributes=attributes, filters=filters, values=contigs.selected, mart=ensembl);
    
    final.names <- gene.names.extracted[,1];
    
  }
  
  
  ## Check whether these genes (final.names) are indexed in the COSMIC database.
  
  # Load the gene census from COSMIC.
  
  cosmic.gene.census <- read.csv('/Users/dmh/Desktop/utilities/cancer_gene_census.csv');
  
  
  # Find those genes present in the COSMIC gene census.
  
  if(synonyms){ # Including gene name synonyms.
    
    gene.synonyms <- unlist(strsplit(as.character(cosmic.gene.census$Synonyms), ','));
    all.gene.names.with.synonyms <- c(as.character(cosmic.gene.census[,1]), gene.synonyms);
    
    matches.synonyms <- final.names %in% all.gene.names.with.synonyms;
    genes.in.census.with.synonyms <- final.names[matches.synonyms];
    
    return(genes.in.census.with.synonyms);
    
  } else{ # Without including gene name synonyms.
    
    matches <- final.names %in% as.character(cosmic.gene.census[,1]);
    genes.in.census <- final.names[matches];
    
    return(genes.in.census);
    
  }
}


### Find genes in amplified region in Chr 1.

## Extract data from chromosome 1.

chr1.data <- extract.chr(data, contig.real.lengths, 'Chr1');
plot(chr1.data$genome.pos, chr1.data$no.reads,pch=20,col="red",ylim=c(0,10000),cex=0.5,xlab="", ylab="Depth/10kb")

## Plot vertical lines to identify region.

for (i in seq(from=1, to=700000000, by=10000000)){
  
  abline(v=i,lty=1,col="blue",lwd=2)
}

abline(v=250000000, lty=1,col="green",lwd=2)
abline(v=270000000, lty=1,col="green",lwd=2)


## Extract data from that region.

indices.interest <- (chr1.data$genome.pos >= 250000000 &  chr1.data$genome.pos <= 270000000);
data.of.interest <- chr1.data[indices.interest,];

plot(data.of.interest$genome.pos, data.of.interest$no.reads,pch=20,col="red",ylim=c(0,10000),cex=0.5,xlab="", ylab="Depth/10kb")

for (i in seq(from=250000000, to=270000000, by=1000000)){
  
  abline(v=i,lty=1,col="blue",lwd=2)
}

## Zoom in.

indices.1 <- (chr1.data$genome.pos >= 251000000 &  chr1.data$genome.pos <= 252500000);
indices.2 <- (chr1.data$genome.pos >= 253000000 &  chr1.data$genome.pos <= 258500000);
indices.3 <- (chr1.data$genome.pos >= 260800000 &  chr1.data$genome.pos <= 267500000);

data.of.interest.1 <- chr1.data[indices.1,];
data.of.interest.2 <- chr1.data[indices.2,];
data.of.interest.3 <- chr1.data[indices.3,];

## Find candidate driver genes in COSMIC database.

contigs.of.interest.chr1 <- c('Chr1_supercontig_000000181', 
                              'Chr1_supercontig_000000182', 'Chr1_supercontig_000000183', 'Chr1_supercontig_000000184', 'Chr1_supercontig_000000185',
                              'Chr1_supercontig_000000190', 'Chr1_supercontig_000000191', 'Chr1_supercontig_000000192', 'Chr1_supercontig_000000193', 
                              'Chr1_supercontig_000000194', 'Chr1_supercontig_000000195');

chr1.no.syn.no.orth.genes <- overlap.with.cosmic(contigs.of.interest.chr1, F, F);
ch1.syn.no.orth.genes <- overlap.with.cosmic(contigs.of.interest.chr1, T, F);
chr1.no.syn.orth.genes <- overlap.with.cosmic(contigs.of.interest.chr1, F, T);
chr1.syn.orth.genes <- overlap.with.cosmic(contigs.of.interest.chr1, T, T);

##### End of the script #####
