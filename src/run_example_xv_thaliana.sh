base_dir=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/thaliana/reads
#gio_file=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/genomes/A_thaliana_magic/Col_0/Col_0.gio/genome.config
gio_file=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/genomes/A_thaliana_magic/Col_0/genome.gio/genome.config
#bam_file=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/reads/Col_0_Chr_nss.bam

bam_file1=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana/reads/SRX006192.new.bam
bam_file2=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana/reads/SRX006681.new.bam
bam_file3=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana/reads/SRX006682.new.bam
bam_file4=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana/reads/SRX006688.new.bam
bam_file5=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana/reads/SRX006690.new.bam
bam_file6=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana/reads/SRX006692.new.bam
bam_file7=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana/reads/SRX006704.new.bam

mkdir $base_dir

elegans_options="-maxel 8000 -reglen 1.5 -maxic 10000 -minic 20 -maxin 50000 -mm 0 -exm 3 -indt 100 -exd 10 -tf 1.0 -inscf 3 -excut 3 -toff 50 -el 10"
gff_file=${base_dir}/genes-elegans-options.gff3
reg_file=${base_dir}/genes-elegans-options-regions-ss-reglen1.5.txt
log_file=${base_dir}/genes-elegans-options.log
echo "infer_genes -gio $gio_file -bam $bam_file1 -bam $bam_file2 -bam $bam_file3 -bam $bam_file4 -bam $bam_file5 -bam $bam_file6 -bam $bam_file7 -gff $gff_file -reg $reg_file $elegans_options -ss "
time ./infer_genes -gio $gio_file -bam $bam_file1 -bam $bam_file2 -bam $bam_file3 -bam $bam_file4 -bam $bam_file5 -bam $bam_file6 -bam $bam_file7 -gff $gff_file -reg $reg_file $elegans_options -ss | tee $log_file

#mm=0
#for inscf in 1 3
#do
#	for id in 50 100 150
#	do
#		for ed in 10 #5 10 20
#		do
#			exm=1 #1 3 5 # 1 was also not bad a bit less specific
#			excut=3
#			for maxic in 5000 #20000
#			do 
#				tf=1.0 #for tf in 1.0 # 2.0 #0.0 0.25 0.5
#				for toff in 10 #0 10 50 100
#				do 
#					gff_file=${base_dir}/genes-inscf$inscf-exm$exm-excut$excut-indt$id-exd$ed-tf${tf}-toff${toff}-mm0-el9.gff3
#					reg_file=${base_dir}/genes-inscf$inscf-exm$exm-excut$excut-indt$id-exd$ed-tf${tf}-toff${toff}-mm0-el9_regions-maxic$maxic.txt
#					echo $gff_file
#					options="$opt -mm $mm -exm $exm -indt $id -exd $ed -tf $tf -inscf $inscf -excut $excut -toff $toff -maxic $maxic"
#					time ./infer_genes -gio $gio_file -bam $bam_file1 -bam $bam_file2 -bam $bam_file3 -bam $bam_file4 -bam $bam_file5 -bam $bam_file6 -bam $bam_file7 -gff $gff_file -reg $reg_file $options
#				done
#			done
#		done
#	done
#done

