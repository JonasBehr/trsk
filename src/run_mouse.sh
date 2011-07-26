base_dir=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/mouse/reads
gio_file=/fml/ag-raetsch/nobackup/projects/sequencing_runs/mouse/cufflinks/gio/genome.config
bam_file=/fml/ag-raetsch/nobackup/projects/sequencing_runs/mouse/cufflinks/GSM521257_60hr.bam

mkdir $base_dir

for reglen in  0.5 #0.99 1.5 2
do
	options="-maxel 8000 -reglen $reglen -maxic 10000 -minic 20 -maxin 50000 -mm 0 -exm 3 -indt 100 -exd 10 -tf 1.0 -inscf 3 -excut 3 -toff 200 -el 15"
	gff_file=${base_dir}/genes_toff200.gff3
	reg_file=${base_dir}/genes-reglen$reglen_toff200.txt
	log_file=${base_dir}/genes_toff200.log
	time ./infer_genes -gio $gio_file -bam $bam_file -gff $gff_file -reg $reg_file $options -ss | tee $log_file
done
