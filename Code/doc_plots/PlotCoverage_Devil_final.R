##################################################################################################
##                                  Daniel Elias Martin Herranz                                 ##
##                                 MPhil in Computational Biology                               ##
##                                    University of Cambridge                                   ##
##                                          August 2015                                         ##
##################################################################################################

##### Create the DOC (depth of coverage) plots for the Tasmanian devil samples #####

## We use as input the results from the search of coverageBed.

## The input file with the coverage has the following columns:
# 1.chromosome/contig     2.window_start    3.window_end      4.number_features_in_BAM_overlapping_with_window (in our case number of reads in the window)
# 5.number_bases_window_no_zero_coverage    6.window_length   7.fraction_bases_window_with_non-zero_coverage 


## Arguments from the command line.

args <- commandArgs(trailingOnly = TRUE); 

# 1st arg: absolute path where the coverageBed output file is located.

input.path <- args[1];

# 2nd arg: absolute path where we want to output the plots.

output.path <- args[2];

# 3rd arg: absolute path where the chr_lengths.txt file is located. This file contains the length for the
#          different contigs in a tab-separated manner, including a header. e.g.
#
#          chr  length
#          Chr1_supercontig_000000000	  2190000
#          Chr1_supercontig_000000001	  380000
#          Chr1_supercontig_000000002	  1080000

chr.path <- args[3];

# rest of the args: integer numbers specifying the maximum value in the Y-axis of the plot (zooming).

zoom.levels <- args[-(1:3)];


## Read all the data that we need.

print('Reading data ...')

data = read.table(input.path,sep="\t",header=F,row.names=NULL);
names(data) = c("chr","startpos","endpos","no.reads","bases.covered","segment.length","fractional.coverage")
data = data[order(data$chr,data$startpos),]
contig.lengths = read.table(chr.path, sep="\t",header=T)


## Calculate cumulative starting positions. If the genome.pos object has been saved previosuly,
# this saves a lot of time from the loop. It has to be placed in the same folder as the script.

cum.lengths = cumsum(contig.lengths$length)
data$genome.pos = data$startpos

# for(c in 2:nrow(contig.lengths)){
#   print(c);
#   data$genome.pos[data$chr==contig.lengths$chr[c]] = data$genome.pos[data$chr==contig.lengths$chr[c]] + cum.lengths[c-1]
# }

load("genome_pos.RData");
data$genome.pos <- genome.pos; # After loading genome_pos, don't redo the loop.


## Create a function to make the DOC plot.

# plot.name: name of the plot (end in .png)
# output.path: path to output the plot.
# zoom: level of zooming in the plot (indicates the maximum value in the Y-axis). It must be 'max' or an integer. 

plot.cn <- function(plot.name, output.path, zoom) {
  
  setwd(output.path);
  
  png(plot.name, width=2000,height=400)
  
  # Plot information per chromosome.
  
  indexes.change.of.chr <- cumsum((rle(as.integer(as.factor(substr(contig.lengths[,1], 1, 4)))))$lengths);
  cum.length.per.chr <- cum.lengths[indexes.change.of.chr];
  #chr.names <- levels(as.factor(substr(contig.lengths[,1], 1, 4)));
  chr.names <- c('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'ChrX', 'ChrU');
  
  if(zoom == 'max'){
    
    plot(data$genome.pos,data$no.reads,pch=20,col="red",ylim=c(0,max(data$no.reads)),cex=0.5,xlab="",ylab="Depth/10kb", xaxt="n")
    abline(v=0,lty=1,col="blue",lwd=4)
    for(c in 1:(length(cum.length.per.chr)-1)){ # Do not plot the title ChrU
      abline(v=cum.length.per.chr[c],lty=1,col="blue",lwd=4)
      if(c==1){
        text(cum.length.per.chr[c]/2, max(data$no.reads),chr.names[c],pos=1,cex=2)
      }else{
        text((cum.length.per.chr[c-1]+cum.length.per.chr[c])/2,max(data$no.reads),chr.names[c],pos=1,cex=2)
      }
    }
  }
  
  if(is.numeric(zoom)){
    
    plot(data$genome.pos,data$no.reads,pch=20,col="red",ylim=c(0,zoom),cex=0.5,xlab="",ylab="Depth/10kb", xaxt="n")
    abline(v=0,lty=1,col="blue",lwd=4)
    for(c in 1:(length(cum.length.per.chr)-1)){ # Do not plot the title ChrU
      abline(v=cum.length.per.chr[c],lty=1,col="blue",lwd=4)
      if(c==1){
        text(cum.length.per.chr[c]/2, zoom,chr.names[c],pos=1,cex=2)
      }else{
        text((cum.length.per.chr[c-1]+cum.length.per.chr[c])/2, zoom,chr.names[c],pos=1,cex=2)
      }
    }
  }
  
  dev.off();
  
}


## Plot the results.

name.of.input <- unlist(strsplit(input.path, '/'))[length(unlist(strsplit(input.path, '/')))];
name.for.plot <- substr(name.of.input, 1, nchar(name.of.input)-4); # Remove .txt from the name

for(zoom in zoom.levels){
  
  print(paste('Plotting', zoom, '...'));
  plot.cn(paste0(name.of.input, '_', zoom, '.png'), output.path, as.numeric(zoom));
  
}

plot.cn(paste0(name.of.input, '_max', '.png'), output.path, 'max');

##### End of the script #####
