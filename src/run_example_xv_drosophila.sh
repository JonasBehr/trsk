base_dir=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/drosophila/reads
gio_file=/fml/ag-raetsch/share/projects/rgasp/genomes/drosophila/drosophila.gio/genome.config
bam_file=/fml/ag-raetsch/nobackup/projects/sequencing_runs/D_melanogaster/reads/D_melanogaster_L3.5.best.bam

mkdir $base_dir

#opt="-maxin 20000 -maxel 8000 -reglen 0.66"
#genes-inscf1-exm3-excut3-indt50-exd10-tf1.0-toff50-incut3.gff3
#elegans_options="-maxel 8000 -reglen 0.66 -maxic 10000 -minic 20 -maxin 50000 -mm 0 -exm 3 -indt 100 -exd 10 -tf 1.0 -inscf 3 -excut 3 -toff 50 -el 15"
elegans_options="-maxel 8000 -reglen 0.99 -maxic 10000 -minic 20 -maxin 50000 -mm 0 -exm 3 -indt 100 -exd 10 -tf 1.0 -inscf 3 -excut 3 -toff 50 -el 15"
gff_file=${base_dir}/genes-elegans-options.gff3
reg_file=${base_dir}/genes-elegans-options-regions-ss-reglen0.99.txt
log_file=${base_dir}/genes-elegans-options.log
time ./infer_genes -gio $gio_file -bam $bam_file -gff $gff_file -reg $reg_file $elegans_options -ss | tee $log_file

