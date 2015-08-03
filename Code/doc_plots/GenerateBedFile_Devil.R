##################################################################################################
##                                  Daniel Elias Martin Herranz                                 ##
##                                 MPhil in Computational Biology                               ##
##                                    University of Cambridge                                   ##
##                                          August 2015                                         ##
##################################################################################################

# We are thankful to David C. Wedge, from the  Wellcome Trust Sanger Institute, who kindly provided the 
# original code and pipeline outline, which we adapted for our own purposes.

##### Generate Bed file with genomic windows of a certain size for the Tasmanian devil samples #####

setwd("/lustre/scratch104/sanger/dmh/cn_analysis");

# contig_and_length.txt is a tab-separated file which contains the name of each contig and its length.
# It does not have a header. e.g.
#
#         Chr1_supercontig_000000000      2180112
#         Chr1_supercontig_000000001	    373546
#         Chr1_supercontig_000000002	    1075293
#         Chr1_supercontig_000000003	    174856
#         Chr1_supercontig_000000004	    817733

contig.and.length = read.table("contig_and_length.txt",sep="\t",header=F,stringsAsFactors=F,row.names=NULL);
names(contig.and.length)[1:2] = c("chr","pos");

chrs = unique(contig.and.length$chr);
maxpos = vector(mode="integer",length=length(chrs));
bed.table = array(NA,c(0,3));

windowsize = 10000; # Adjust window size

for(c in 1:length(chrs)){  
  print(c);
  chr=chrs[c]
  maxpos[c] = windowsize*ceiling(max(contig.and.length$pos[contig.and.length$chr==chr])/windowsize)
  bed.table = rbind(bed.table,cbind(chr,format(seq(0,maxpos[c]-windowsize,windowsize),scientific=F,trim=T),format(seq(windowsize,maxpos[c],windowsize),scientific=F,trim=T)))
}

write.table(cbind(chrs,maxpos),"chr_lengths.txt",col.names = c("chr","length"),row.names=F,quote=F,sep="\t")
write.table(bed.table,paste("windows_",windowsize/1000,"kb.bed",sep=""),sep="\t",quote=F,row.names=F,col.names=F)

##### End of the script #####
