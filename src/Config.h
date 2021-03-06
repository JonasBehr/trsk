#ifndef _CONFIG_H__
#define _CONFIG_H__
#include <stdio.h>
#include <stdlib.h>
#include <vector>
	using std::vector;

class Config 
{
	public: 	
		Config();
		Config(int argc, char *argv[])
		{
			default_values();
			parseCommandLine(argc, argv);
		};
		void default_values();
		int parseCommandLine(int argc, char *argv[]);
		void print(_IO_FILE*& fd)
		{
			fprintf(fd,"\nRead filter:\n");
			fprintf(fd,"\tmax_intron_len (-maxin):\t%i\n", max_intron_len);
			fprintf(fd,"\tmm_filter (-mm):\t\t%i\n", mm_filter);
			fprintf(fd,"\tel_filter (-el):\t\t%i\n", el_filter);
			fprintf(fd,"\nTranscript filter:\n");
			fprintf(fd,"\tmax_exon_len (-maxel):\t\t%i\n", max_exon_len);
			fprintf(fd,"\tmin_exon_len (-minel):\t\t%i\n", min_exon_len);
			fprintf(fd,"\texon_mean (-exm):\t\t%i\n", exon_mean);
			fprintf(fd,"\texon_term_thresh (-exts):\t%i\n", exon_term_thresh);
			fprintf(fd,"\texon_drop (-exd):\t\t%i\n", exon_drop);
			fprintf(fd,"\texon_cut (-excut):\t\t%i\n", exon_cut);
			fprintf(fd,"\tintron_cut (-incut):\t\t%.2f\n", intron_cut);
			fprintf(fd,"\tintron_conf (-incf):\t\t%i\n", intron_conf);
			fprintf(fd,"\tintron_dist (-indt):\t\t%i\n", intron_dist);
			fprintf(fd,"\tintron_seed_conf (-inscf):\t%i\n", intron_seed_conf);
			fprintf(fd,"\treject_retained_introns (-ri):\t%i\n", reject_retained_introns);
			fprintf(fd,"\tterm_filter (-tf):\t\t%.2f\n", term_filter);
			fprintf(fd,"\nORF:\n");
			fprintf(fd,"\tfind_orf (-orf/-noorf):\t\t%i\n", find_orf);
			fprintf(fd,"\tmin_orf_len (-orflen):\t\t%i\n", min_orf_len);
			fprintf(fd,"\tmin_orf_sep (-orfsep):\t\t%.2f\n", min_orf_sep);
			fprintf(fd,"\tterm_offset (-toff):\t\t%i\n", term_offset);

			fprintf(fd,"\nRegions:\n");
			if (!have_reg_file)
				fprintf(fd,"(only effective if reg file (-reg) is given)\n");
			fprintf(fd,"\tregion_rel_length (-reglen):\t%.2f\n", region_rel_length);
			fprintf(fd,"\tmin_intergenic_len (-minic):\t%i\n", min_intergenic_len);
			fprintf(fd,"\tmax_intergenic_len (-maxic):\t%i\n", max_intergenic_len);
			fprintf(fd,"\tintergenic_win (-iwin):\t\t%i\n", intergenic_win);
			fprintf(fd,"\tstrand_specific (-ss/-nss):\t%i\n", strand_specific);

			fprintf(fd,"\nSignals:\n");
			if (have_tss_seq_file)
				fprintf(fd,"\ttss sequence file: %s\n", tss_seq_file);
			if (have_tis_seq_file)
				fprintf(fd,"\ttis sequence file: %s\n", tis_seq_file);
			fprintf(fd,"\n");
		};
		int max_exon_len;
		int min_exon_len;
		int max_intron_len;
		int min_intergenic_len;
		int intergenic_win; 
		int max_intergenic_len; 
		int intron_conf;
		int intron_dist;
		int intron_seed_conf;
		bool reject_retained_introns;
		float term_filter;
		int mm_filter; 
		int el_filter;

		/** number of positions with coverage lower than this 
		 * are counted to decide if exon is rejected or not
		 */
		int exon_mean;
		int exon_term_thresh;
		int exon_drop;
		int exon_cut;
		float intron_cut;

		/** determines how long intergenic regions are 
		 *	relative to the next observed island of coverage
		 */
		float region_rel_length; 

		vector<char*> bam_files;
		char* gio_file;
		char* gff_file;
		char* reg_file; 
		char* tis_seq_file; 
		char* tss_seq_file;

		bool have_bam_file;
		bool have_gio_file;
		bool have_reg_file; 
		bool have_gff_file;
		bool have_tis_seq_file;
		bool have_tss_seq_file;
		
		bool strand_specific; 

		/** orf stuff */
		bool find_orf;
		int min_orf_len; 
		float min_orf_sep;
		int term_offset;// find more ORFs by setting the boundaries of the genes a bit wider

		bool find_single_exon_genes_orf;
};
#endif
