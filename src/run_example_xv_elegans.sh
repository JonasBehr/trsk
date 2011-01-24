

base_dir=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/elegans/reads
gio_file=/fml/ag-raetsch/share/projects/rgasp/genomes/elegans/elegans.gio/genome.config

mkdir $base_dir

opt="-maxin 20000 -maxel 8000 -inscf 1 -excut 1  "

#for bam in ".best." ".";
for mm in 1 2 3
do
	for id in 25 50 100;
	do
		for ed in 5 10 20;
		do
			for exm in 1 3 5;
			do 
				#bam_file=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/worms/elegans.5${bam}bam # mem consumption >60G
				bam1=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans/polyA_left_sam_stranded.mapped.2.bam
				bam2=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans/polyA_right_sam_stranded.mapped.2.bam
				#gff_file=${base_dir}/genes-exm$exm-indt$id-exd$ed${bam}gff3
				gff_file=${base_dir}/genes-mm$mm-exm$exm-indt$id-exd$ed.gff3
				echo $gff_file
				options="$opt -mm $mm -exm $exm -indt $id -exd $ed"
				#time ./infer_genes -gio $gio_file -bam $bam_file -gff $gff_file $options
				time ./infer_genes -gio $gio_file -bam $bam1 -bam $bam2 -gff $gff_file $options
			done
		done
	done
done

