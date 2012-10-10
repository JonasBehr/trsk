#!/bin/bash

#bamfile1=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans/polyA_left_sam_stranded.mapped.2.bam
#bamfile2=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans/polyA_right_sam_stranded.mapped.2.bam

#bamfile1=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans.dummy/polyA_left_trim_new.bam
#bamfile2=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans.dummy/polyA_right_trim_new.bam


#genome_config_dir=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/elegans/genome_dir/
#genome_config=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/elegans/genome_dir/genome.config



#output_dir=/fml/ag-raetsch/home/jonas/tmp/RNA_seq_splice_label_unbiased
output_dir=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/human/reads_jin/splice_predictions
genome_config_dir=/fml/ag-raetsch/nobackup/projects/rgasp/genomes/hg19_14/hg19.gio/
genome_config=${genome_config_dir}/genome.config



time ./splice_labels_from_RNA_seq $genome_config $output_dir /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.*seq.5_me4_mm0_mi200000.sorted.bam


matlab -r "addpath ~/svn/projects/mGene_core/; paths; addpath ~/svn/projects/mGene_core/RNA_seq_label; train_splice_signals($genome_config_dir, $output_dir)"
