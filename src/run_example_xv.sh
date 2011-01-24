

base_dir=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/human/reads
gio_file=/fml/ag-raetsch/nobackup/projects/rgasp.2/genomes/human/hg19/hg19.gio/genome.config
bam_file=${base_dir}/sanitized.bam

opt="-maxin 100000 -maxel 8000 "

for i in 3; # 1 2 3 5 9; 
do
	for e in 1; # 2 3;
	do
		for em in 3; #1 2 3 4 5;
		do
			for id in 25 50 100;
			do
				for ed in 3 5 10 20;
				do
					gff_file=${base_dir}/genes-inscf$i-excut$e-exm$em-indt$id-exd$ed.gff3
					echo $gff_file
					options="$opt -inscf $i -excut $e -exm $em -indt $id -exd $ed"
					time ./infer_genes -gio $gio_file -bam $bam_file -gff $gff_file $options
				done
			done
		done
	done
done

