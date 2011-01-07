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
			fprintf(fd,"\nFilter:\n");
			fprintf(fd,"\tmax_exon_len (-maxel):\t\t%i\n", max_exon_len);
			fprintf(fd,"\tmin_exon_len (-minel):\t\t%i\n", min_exon_len);
			fprintf(fd,"\tmax_intron_len (-maxin):\t%i\n", max_intron_len);
			fprintf(fd,"\texon_cut (-excut):\t\t%i\n", exon_cut);
			fprintf(fd,"\tintron_conf (-incf):\t\t%i\n", intron_conf);
			fprintf(fd,"\tintron_seed_conf (-inscf):\t%i\n", intron_seed_conf);
			fprintf(fd,"\nORF:\n");
			fprintf(fd,"\tfind_orf (-orf/-noorf):\t\t%i\n", find_orf);
			fprintf(fd,"\tmin_orf_len (-orflen):\t\t%i\n", min_orf_len);
			fprintf(fd,"\tmin_orf_sep (-orfsep):\t\t%.2f\n", min_orf_sep);
			fprintf(fd,"\nRegions:\n");
			fprintf(fd,"\tregion_rel_length (-reglen):\t%.2f\n", region_rel_length);
			fprintf(fd,"\tmin_intergenic_len (-minic):\t%i\n", min_intergenic_len);
			fprintf(fd,"\tmax_intergenic_len (-maxic):\t%i\n", max_intergenic_len);
			fprintf(fd,"\tintergenic_win (-iwin):\t\t%i\n", intergenic_win);
			fprintf(fd,"\tstrand_specific (-ss/-nss):\t%i\n", strand_specific);
			fprintf(fd,"\n");
		};
		int max_exon_len;
		int min_exon_len;
		int max_intron_len;
		int min_intergenic_len;
		int intergenic_win; 
		int max_intergenic_len; 
		int intron_conf;
		int intron_seed_conf;

		/** number of positions with coverage lower than this 
		 * are counted to decide if exon is rejected or not
		 */
		int exon_cut;

		/** determines how long intergenic regions are 
		 *	relative to the next observed island of coverage
		 */
		float region_rel_length; 

		vector<char*> bam_files;
		char* gio_file;
		char* gff_file;
		char* reg_file; 
		bool have_bam_file;
		bool have_gio_file;
		bool have_reg_file; 
		bool have_gff_file;
		
		bool strand_specific; 

		/** orf stuff */
		bool find_orf;
		int min_orf_len; 
		float min_orf_sep;

};
#endif
