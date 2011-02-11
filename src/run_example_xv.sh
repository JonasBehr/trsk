

base_dir=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/human/reads
gio_file=/fml/ag-raetsch/nobackup/projects/rgasp.2/genomes/human/hg19/hg19.gio/genome.config
bam_file=${base_dir}/sanitized.bam

opt="-maxin 100000 -maxel 8000 "

for i in 3; # 1 2 3 5 9; 
do
	e="3" #for e in 1; # 2 3;
	for tf in 1 #0.0 0.25 0.75 1.0 1.5 2.0;
	do
		for em in 3; #1 2 3 4 5;
		do
			for id in 100; #25 50 100;
			do
				for ed in 10; #3 5 10 20;
				do
					for incut in 3 #0.0 0.25 0.5 0.75 1
					do
						#gff_file=${base_dir}/genes-inscf$i-excut$e-exm$em-indt$id-exd$ed-tf$tf-incut${incut}_term_off.gff3
						gff_file=/tmp/genes-inscf$i-excut$e-exm$em-indt$id-exd$ed-tf$tf-incut${incut}_term_off.gff3
						echo $gff_file
						options="$opt -inscf $i -excut $e -exm $em -indt $id -exd $ed -tf $tf -incut $incut"
						time ./infer_genes -gio $gio_file -bam $bam_file -gff $gff_file $options
					done
				done
			done
		done
	done
done

