#!/bin/bash

bamfile=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans/polyA_left_sam.filtered.2.bam
#bamfile=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/tracks/elegans/polyA_left/polyA_left_I+_el15_mm1_spliced.bam

genome_config=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/elegans/genome_dir/genome.config

time ./spliced_nucleotide_db $bamfile  $genome_config
