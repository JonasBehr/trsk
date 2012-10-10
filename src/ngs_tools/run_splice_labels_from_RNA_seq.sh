
#organism=elegans
#organism=pacificus
#organism=japonica
organism=Col_0
#organism=drosophila
#fn_genome=/fml/ag-raetsch/nobackup/projects/sequencing_runs/mouse/cufflinks/gio/genome.config
#fn_genome=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/genomes/$organism/${organism}_condensed.gio/genome.config
#dir_genome=/fml/ag-raetsch/share/projects/rgasp/genomes/elegans/elegans.gio/
#dir_genome=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/genomes/japonica2/japonica2_condensed.gio/
dir_genome=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/genomes/A_thaliana_magic/Col_0/Col_0.gio/
#dir_genome=/fml/ag-raetsch/share/projects/rgasp/genomes/drosophila/drosophila.gio
fn_genome=$dir_genome/genome.config
#fn_out=/tmp/
fn_out=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/$organism/layer_1_label_gen_alt_in_positive_class/
#fn_out=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/$organism/layer_1_label_gen/
#fn_bam=/fml/ag-raetsch/nobackup/projects/sequencing_runs/mouse/reads/cufflinks_SRX017795.sanitized.sorted.chr.bam
#fn_bam=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/worms/pac.bam
#fn_bam=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/worms/pacificus.bam
#fn_bam=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/worms/japonica.bam
fn_bam=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/reads/Col_0_Chr_full_ns.bam
#fn_bam=/fml/ag-raetsch/nobackup/projects/sequencing_runs/D_melanogaster/reads/D_melanogaster_L3.5.bam
#fn_bam=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans/polyA_left_sam_stranded.mapped.2.bam
#fn_bam2=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans/polyA_right_sam_stranded.mapped.2.bam

#echo "./splice_labels_from_RNA_seq $fn_genome $fn_out $fn_bam"

## run 
#time ./splice_labels_from_RNA_seq $fn_genome $fn_out $fn_bam $fn_bam2
time ./splice_labels_from_RNA_seq $fn_genome $fn_out $fn_bam $fn_bam2


mat_src=~jonas/svn/projects/mGene_core
#echo "cd $mat_src; dbstop error; paths; train_splice_signals('$dir_genome', '$fn_out')"
matlab -nojvm -nosplash -r "cd $mat_src; dbstop error; paths; train_splice_signals('$dir_genome', '$fn_out')"
