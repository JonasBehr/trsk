base_dir=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/human/reads_jin
gio_dir=/fml/ag-raetsch/nobackup/projects/rgasp/genomes/hg19_14/hg19.gio/
gio_file=${gio_dir}/genome.config
#bam_file="/fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.cb_f_seq.5.sorted.bam \
#     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.ht_m_seq.5.sorted.bam \
#     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.lv_f_seq.5.sorted.bam \
#     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.cb_m_seq.5.sorted.bam \
#     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.kd_f_seq.5.sorted.bam \
#     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.lv_m_seq.5.sorted.bam \
#     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.ht_f_seq.5.sorted.bam \
#     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.kd_m_seq.5.sorted.bam \
#     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.tt_m_seq.5.sorted.bam"
bam_file="/fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.cb_f_seq.5_me4_mm0_mi200000.sorted.bam \
     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.ht_m_seq.5_me4_mm0_mi200000.sorted.bam \
     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.lv_f_seq.5_me4_mm0_mi200000.sorted.bam \
     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.cb_m_seq.5_me4_mm0_mi200000.sorted.bam \
     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.kd_f_seq.5_me4_mm0_mi200000.sorted.bam \
     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.lv_m_seq.5_me4_mm0_mi200000.sorted.bam \
     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.ht_f_seq.5_me4_mm0_mi200000.sorted.bam \
     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.kd_m_seq.5_me4_mm0_mi200000.sorted.bam \
     -bam /fml/ag-raetsch/nobackup2/projects/bonobo/alignments/jin.human.tt_m_seq.5_me4_mm0_mi200000.sorted.bam"

mkdir -p $base_dir

#best results for human rgasp data: trsk genes-inscf3-excut3-exm3-indt100-exd10-tf0.25-incut3

for reglen in  0.5 #0.99 1.5 2
do
	options="-maxel 8000 -reglen $reglen -maxic 10000 -minic 20 -maxin 50000 -mm 0 -exm 3 -indt 100 -exd 10 -tf 1.0 -inscf 3 -excut 3 -toff 50 -el 15"
	gff_file=${base_dir}/genes.gff3
	reg_file=${base_dir}/genes-reglen$reglen.txt
	log_file=${base_dir}/genes.log
	time ./infer_genes -gio $gio_file -bam $bam_file -gff $gff_file -reg $reg_file $options -ss #| tee $log_file
done


# parse gff file
~/svn/releases/mGeneToolbox-0.2.0/release/src/parsegff/GFFParser.sh $gff_file ${base_dir}/genes.mat ${gio_dir}

pred_dir=/fml/ag-raetsch/nobackup/projects/rgasp/mgene_predictions/human/lsl/trsk_jin_unfiltered/output/genome_wide_predictions/
mkdir -p $pred_dir
ln -s ${base_dir}/genes.mat $pred_dir/genes.mat

matlab -r "addpath ~/svn/projects/mGene_core; paths; eval_conf('human', 'trsk_jin_unfiltered', 0, 0); exit"
