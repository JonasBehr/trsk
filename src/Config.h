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
			fprintf(fd,"max_exon_len:\t%i\n", max_exon_len);
			fprintf(fd,"min_exon_len:\t%i\n", min_exon_len);
			fprintf(fd,"max_intron_len:\t%i\n", max_intron_len);
			fprintf(fd,"min_intergenic_len:\t%i\n", min_intergenic_len);
			fprintf(fd,"intergenic_win:\t%i\n", intergenic_win);
			fprintf(fd,"exon_cut:\t\t%i\n", exon_cut);
			fprintf(fd,"region_rel_length:\t%.2f\n", region_rel_length);
//			fprintf(fd,"strand_specific\t%b\n", strand_specific);
//			fprintf(fd,"find_orf:\t%b\n", find_orf);
		};
		int max_exon_len;
		int min_exon_len;
		int max_intron_len;
		int min_intergenic_len;
		int intergenic_win; 

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
		bool find_orf;
};
#endif
