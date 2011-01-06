#include "Config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Config::Config()
{
		max_exon_len = 8000;
		max_intron_len = 30000;
		min_exon_len = 10;
		min_intergenic_len = 50;
		intergenic_win = 100;
		strand_specific = false;
}
int Config::parseCommandLine(int argc, char *argv[])
{
    int i;
    char not_defined;
    char has_index = 0;
    char has_genome = 0;

    for (i = 1; i < argc; i++) 
	{
        not_defined = 1;

        //genome file
        if (strcmp(argv[i], "-i") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -i\n") ;
               //usage();
                exit(1);
            }
		}
	}
}
