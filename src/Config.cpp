#include "Config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Config::Config()
{
	default_values();
}
void Config::default_values()
{
	max_exon_len = 8000;
	max_intron_len = 30000;
	min_exon_len = 10;
	min_intergenic_len = 50;
	intergenic_win = 100;
	exon_cut = 3;
	strand_specific = false;
	region_rel_length = 0.25;
	have_bam_file = false;
	have_gio_file = false;
	have_reg_file = false; 
	have_gff_file = false;

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

        if (strcmp(argv[i], "-maxel") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -maxel\n") ;
               //usage();
                exit(1);
            }
			i++; 
			max_exon_len = atoi(argv[i]);
		}
        else if (strcmp(argv[i], "-minel") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -minel\n") ;
               //usage();
                exit(1);
            }
			i++; 
			min_exon_len = atoi(argv[i]);
		}
	    else if (strcmp(argv[i], "-reglen") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -reglen\n") ;
               //usage();
                exit(1);
            }
			i++; 
			region_rel_length = atof(argv[i]);
		}
		else if (strcmp(argv[i], "-ss") == 0) 
		{

            not_defined = 0;
			strand_specific = true;
		}
		else if (strcmp(argv[i], "-bam") == 0)
		{
			not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -bam\n") ;
               //usage();
                exit(1);
            }
			i++; 
			bam_files.push_back(argv[i]);
			have_bam_file = true;
		}
		else if (strcmp(argv[i], "-gio") == 0)
		{
			not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -gio\n") ;
               //usage();
                exit(1);
            }
			i++; 
			gio_file = argv[i];
			have_gio_file = true;
		}
		else if (strcmp(argv[i], "-reg") == 0)
		{
			not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -reg\n") ;
               //usage();
                exit(1);
            }
			i++; 
			reg_file = argv[i];
			have_reg_file = true;
		}
		else if (strcmp(argv[i], "-gff") == 0)
		{
			not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -gff\n") ;
               //usage();
                exit(1);
            }
			i++; 
			gff_file = argv[i];
			have_gff_file = true;
		}
	}
}
