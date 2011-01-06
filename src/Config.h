#ifndef _CONFIG_H__
#define _CONFIG_H__
#include <stdio.h>
#include <stdlib.h>

class Config 
{
	public: 	
		Config();
		int parseCommandLine(int argc, char *argv[]);
		void print(_IO_FILE*& fd)
		{
			fprintf(fd,"max_exon_len:\t%i\n", max_exon_len);
		};
		int max_exon_len;
		int min_exon_len;
		int max_intron_len;
		int min_intergenic_len;
		int intergenic_win; 
		bool strand_specific; 
};
#endif
