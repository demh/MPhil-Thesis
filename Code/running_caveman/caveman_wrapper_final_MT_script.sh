#!/bin/bash

FILE=$1

mkdir $FILE

caveman.pl \
-outdir /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/$FILE \
-reference /lustre/scratch104/sanger/dmh/references/7.0/Sarcophilus_harrisii.DEVIL7.0.70.dna.toplevel.fa.fai \
-tumour-bam /lustre/scratch104/sanger/dmh/trial_caveman/bam_files/$FILE \
-normal-bam /lustre/scratch104/sanger/dmh/trial_caveman/simulated_reads/alignment/aln-pe_MT_picard.bam \
-ignore-file /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/ignore_all_but_MT.tsv \
-tumour-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_tumour.txt \
-normal-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_normal.txt \
-tum-cn-default 5 \
-norm-cn-default 1 \
-species 'DEVIL' \
-species-assembly 'devil7.0' \
-flag-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files \
-germline-indel /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files/germline_indel.bed \
-unmatched-vcf /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/unmatched_vcf \
-seqType 'genomic' \
-normal-contamination 0 \
-flagConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.vcf.config.ini \
-flagToVcfConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.to.vcf.convert.ini \
-process 'setup' \
-annot-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/annotation_bed_files


caveman.pl \
-outdir /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/$FILE \
-reference /lustre/scratch104/sanger/dmh/references/7.0/Sarcophilus_harrisii.DEVIL7.0.70.dna.toplevel.fa.fai \
-tumour-bam /lustre/scratch104/sanger/dmh/trial_caveman/bam_files/$FILE \
-normal-bam /lustre/scratch104/sanger/dmh/trial_caveman/simulated_reads/alignment/aln-pe_MT_picard.bam \
-ignore-file /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/ignore_all_but_MT.tsv \
-tumour-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_tumour.txt \
-normal-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_normal.txt \
-tum-cn-default 5 \
-norm-cn-default 1 \
-species 'DEVIL' \
-species-assembly 'devil7.0' \
-flag-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files \
-germline-indel /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files/germline_indel.bed \
-unmatched-vcf /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/unmatched_vcf \
-seqType 'genomic' \
-normal-contamination 0 \
-flagConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.vcf.config.ini \
-flagToVcfConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.to.vcf.convert.ini \
-process 'split' \
-index 27889 \
-annot-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/annotation_bed_files



caveman.pl \
-outdir /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/$FILE \
-reference /lustre/scratch104/sanger/dmh/references/7.0/Sarcophilus_harrisii.DEVIL7.0.70.dna.toplevel.fa.fai \
-tumour-bam /lustre/scratch104/sanger/dmh/trial_caveman/bam_files/$FILE \
-normal-bam /lustre/scratch104/sanger/dmh/trial_caveman/simulated_reads/alignment/aln-pe_MT_picard.bam \
-ignore-file /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/ignore_all_but_MT.tsv \
-tumour-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_tumour.txt \
-normal-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_normal.txt \
-tum-cn-default 5 \
-norm-cn-default 1 \
-species 'DEVIL' \
-species-assembly 'devil7.0' \
-flag-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files \
-germline-indel /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files/germline_indel.bed \
-unmatched-vcf /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/unmatched_vcf \
-seqType 'genomic' \
-normal-contamination 0 \
-flagConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.vcf.config.ini \
-flagToVcfConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.to.vcf.convert.ini \
-process 'split_concat' \
-annot-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/annotation_bed_files



caveman.pl \
-outdir /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/$FILE \
-reference /lustre/scratch104/sanger/dmh/references/7.0/Sarcophilus_harrisii.DEVIL7.0.70.dna.toplevel.fa.fai \
-tumour-bam /lustre/scratch104/sanger/dmh/trial_caveman/bam_files/$FILE \
-normal-bam /lustre/scratch104/sanger/dmh/trial_caveman/simulated_reads/alignment/aln-pe_MT_picard.bam \
-ignore-file /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/ignore_all_but_MT.tsv \
-tumour-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_tumour.txt \
-normal-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_normal.txt \
-tum-cn-default 5 \
-norm-cn-default 1 \
-species 'DEVIL' \
-species-assembly 'devil7.0' \
-flag-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files \
-germline-indel /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files/germline_indel.bed \
-unmatched-vcf /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/unmatched_vcf \
-seqType 'genomic' \
-normal-contamination 0 \
-flagConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.vcf.config.ini \
-flagToVcfConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.to.vcf.convert.ini \
-process 'mstep' \
-annot-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/annotation_bed_files



caveman.pl \
-outdir /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/$FILE \
-reference /lustre/scratch104/sanger/dmh/references/7.0/Sarcophilus_harrisii.DEVIL7.0.70.dna.toplevel.fa.fai \
-tumour-bam /lustre/scratch104/sanger/dmh/trial_caveman/bam_files/$FILE \
-normal-bam /lustre/scratch104/sanger/dmh/trial_caveman/simulated_reads/alignment/aln-pe_MT_picard.bam \
-ignore-file /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/ignore_all_but_MT.tsv \
-tumour-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_tumour.txt \
-normal-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_normal.txt \
-tum-cn-default 5 \
-norm-cn-default 1 \
-species 'DEVIL' \
-species-assembly 'devil7.0' \
-flag-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files \
-germline-indel /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files/germline_indel.bed \
-unmatched-vcf /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/unmatched_vcf \
-seqType 'genomic' \
-normal-contamination 0 \
-flagConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.vcf.config.ini \
-flagToVcfConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.to.vcf.convert.ini \
-process 'merge' \
-annot-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/annotation_bed_files



caveman.pl \
-outdir /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/$FILE \
-reference /lustre/scratch104/sanger/dmh/references/7.0/Sarcophilus_harrisii.DEVIL7.0.70.dna.toplevel.fa.fai \
-tumour-bam /lustre/scratch104/sanger/dmh/trial_caveman/bam_files/$FILE \
-normal-bam /lustre/scratch104/sanger/dmh/trial_caveman/simulated_reads/alignment/aln-pe_MT_picard.bam \
-ignore-file /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/ignore_all_but_MT.tsv \
-tumour-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_tumour.txt \
-normal-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_normal.txt \
-tum-cn-default 5 \
-norm-cn-default 1 \
-species 'DEVIL' \
-species-assembly 'devil7.0' \
-flag-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files \
-germline-indel /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files/germline_indel.bed \
-unmatched-vcf /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/unmatched_vcf \
-seqType 'genomic' \
-normal-contamination 0 \
-flagConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.vcf.config.ini \
-flagToVcfConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.to.vcf.convert.ini \
-process 'estep' \
-annot-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/annotation_bed_files



caveman.pl \
-outdir /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/$FILE \
-reference /lustre/scratch104/sanger/dmh/references/7.0/Sarcophilus_harrisii.DEVIL7.0.70.dna.toplevel.fa.fai \
-tumour-bam /lustre/scratch104/sanger/dmh/trial_caveman/bam_files/$FILE \
-normal-bam /lustre/scratch104/sanger/dmh/trial_caveman/simulated_reads/alignment/aln-pe_MT_picard.bam \
-ignore-file /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/ignore_all_but_MT.tsv \
-tumour-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_tumour.txt \
-normal-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_normal.txt \
-tum-cn-default 5 \
-norm-cn-default 1 \
-species 'DEVIL' \
-species-assembly 'devil7.0' \
-flag-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files \
-germline-indel /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files/germline_indel.bed \
-unmatched-vcf /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/unmatched_vcf \
-seqType 'genomic' \
-normal-contamination 0 \
-flagConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.vcf.config.ini \
-flagToVcfConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.to.vcf.convert.ini \
-process 'merge_results' \
-annot-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/annotation_bed_files



caveman.pl \
-outdir /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/$FILE \
-reference /lustre/scratch104/sanger/dmh/references/7.0/Sarcophilus_harrisii.DEVIL7.0.70.dna.toplevel.fa.fai \
-tumour-bam /lustre/scratch104/sanger/dmh/trial_caveman/bam_files/$FILE \
-normal-bam /lustre/scratch104/sanger/dmh/trial_caveman/simulated_reads/alignment/aln-pe_MT_picard.bam \
-ignore-file /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/ignore_all_but_MT.tsv \
-tumour-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_tumour.txt \
-normal-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_normal.txt \
-tum-cn-default 5 \
-norm-cn-default 1 \
-species 'DEVIL' \
-species-assembly 'devil7.0' \
-flag-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files \
-germline-indel /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files/germline_indel.bed \
-unmatched-vcf /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/unmatched_vcf \
-seqType 'genomic' \
-normal-contamination 0 \
-flagConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.vcf.config.ini \
-flagToVcfConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.to.vcf.convert.ini \
-process 'add_ids' \
-annot-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/annotation_bed_files



caveman.pl \
-outdir /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/$FILE \
-reference /lustre/scratch104/sanger/dmh/references/7.0/Sarcophilus_harrisii.DEVIL7.0.70.dna.toplevel.fa.fai \
-tumour-bam /lustre/scratch104/sanger/dmh/trial_caveman/bam_files/$FILE \
-normal-bam /lustre/scratch104/sanger/dmh/trial_caveman/simulated_reads/alignment/aln-pe_MT_picard.bam \
-ignore-file /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/ignore_all_but_MT.tsv \
-tumour-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_tumour.txt \
-normal-cn /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/CN_normal.txt \
-tum-cn-default 5 \
-norm-cn-default 1 \
-species 'DEVIL' \
-species-assembly 'devil7.0' \
-flag-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files \
-germline-indel /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag_bed_files/germline_indel.bed \
-unmatched-vcf /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/unmatched_vcf \
-seqType 'genomic' \
-normal-contamination 0 \
-flagConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.vcf.config.ini \
-flagToVcfConfig /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/flag.to.vcf.convert.ini \
-process 'flag' \
-annot-bed-files /lustre/scratch104/sanger/dmh/trial_caveman/run_caveman_wrapper/post_parameters/annotation_bed_files



