

base_dir=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/elegans/reads
gio_file=/fml/ag-raetsch/share/projects/rgasp/genomes/elegans/elegans.gio/genome.config

mkdir $base_dir

#opt="-maxin 20000 -maxel 8000 -ss"
opt="-maxel 8000 -ss -reglen 0.66 -maxic 10000 -minic 20" # get more close to the old matlab region creation skript

#for bam in ".best." ".";
maxin=50000
mm=0
el=15
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
					gff_file=${base_dir}/genes-inscf$inscf-exm$exm-excut$excut-indt$id-exd$ed-tf${tf}-toff${toff}-el$el-mm$mm-maxin${maxin}.gff3
					reg_file=${base_dir}/genes-inscf$inscf-exm$exm-excut$excut-indt$id-exd$ed-tf${tf}-toff${toff}-el$el-mm$mm-maxin${maxin}_regions.txt
					echo $gff_file
					options="$opt -maxin $maxin -mm $mm -exm $exm -indt $id -exd $ed -tf $tf -inscf $inscf -excut $excut -toff $toff -el $el"
					#time ./infer_genes -gio $gio_file -bam $bam_file -gff $gff_file $options
					#echo $options
					time ./infer_genes -gio $gio_file -bam $bam1 -bam $bam2 -gff $gff_file -reg $reg_file $options 
				done
			done
		done
	done
done

