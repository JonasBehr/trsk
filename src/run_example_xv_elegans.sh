

base_dir=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/elegans/reads
gio_file=/fml/ag-raetsch/share/projects/rgasp/genomes/elegans/elegans.gio/genome.config

mkdir $base_dir

opt="-maxin 20000 -maxel 8000 -ss"

#for bam in ".best." ".";
mm=1
for inscf in 3
do
	for id in 100 #25 50 100 150
	do
		for ed in 10 #5 10 20
		do
			exm=3 #1 3 5 # 1 was also not bad a bit less specific
			for excut in 3
			do 
				tf=1.0 #for tf in 1.0 # 2.0 #0.0 0.25 0.5
				for toff in 50 #0 10 50 100
				do 
					#bam_file=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/worms/elegans.5${bam}bam # mem consumption >60G
					bam1=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans/polyA_left_sam_stranded.mapped.2.bam
					bam2=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/elegans/polyA_right_sam_stranded.mapped.2.bam
					#gff_file=${base_dir}/genes-exm$exm-indt$id-exd$ed${bam}gff3
					gff_file=${base_dir}/genes-inscf$inscf-exm$exm-excut$excut-indt$id-exd$ed-tf${tf}-toff${toff}_short_intron_reject.gff3
					gff_file=/tmp/genes-inscf$inscf-exm$exm-excut$excut-indt$id-exd$ed-tf${tf}-toff${toff}_short_intron_reject.gff3
					reg_file=${base_dir}/genes-inscf$inscf-exm$exm-excut$excut-indt$id-exd$ed-tf${tf}-toff${toff}_short_intron_reject_region_ss.txt
					echo $gff_file
					options="$opt -mm $mm -exm $exm -indt $id -exd $ed -tf $tf -inscf $inscf -excut $excut -toff $toff"
					#time ./infer_genes -gio $gio_file -bam $bam_file -gff $gff_file $options
					echo $options
					#time ./infer_genes -gio $gio_file -bam $bam1 -bam $bam2 -gff $gff_file -reg $reg_file $options 
				done
			done
		done
	done
done

