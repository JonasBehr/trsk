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
	max_intergenic_len = 20000;
	intergenic_win = 100;
	exon_cut = 3;
	strand_specific = false;
	region_rel_length = 0.25;
	have_bam_file = false;
	have_gio_file = false;
	have_reg_file = false; 
	have_gff_file = false;
	find_orf = true;
	min_orf_len = 300;
	min_orf_sep = 0.7;
	intron_conf = 1;
	intron_seed_conf = 3;
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
        else if (strcmp(argv[i], "-maxin") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -maxin\n") ;
               //usage();
                exit(1);
            }
			i++; 
			max_intron_len = atoi(argv[i]);
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
        else if (strcmp(argv[i], "-incf") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -incf\n") ;
               //usage();
                exit(1);
            }
			i++; 
			intron_conf = atoi(argv[i]);
		}
        else if (strcmp(argv[i], "-inscf") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -inscf\n") ;
               //usage();
                exit(1);
            }
			i++; 
			intron_seed_conf = atoi(argv[i]);
		}
        else if (strcmp(argv[i], "-excut") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -excut\n") ;
               //usage();
                exit(1);
            }
			i++; 
			exon_cut = atoi(argv[i]);
		}
        else if (strcmp(argv[i], "-iw") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -iw\n") ;
               //usage();
                exit(1);
            }
			i++; 
			intergenic_win = atoi(argv[i]);
		}
        else if (strcmp(argv[i], "-minic") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -minic\n") ;
               //usage();
                exit(1);
            }
			i++; 
			min_intergenic_len = atoi(argv[i]);
		}
        else if (strcmp(argv[i], "-maxic") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -maxic\n") ;
               //usage();
                exit(1);
            }
			i++; 
			max_intergenic_len = atoi(argv[i]);
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
		else if (strcmp(argv[i], "-nss") == 0) 
		{

            not_defined = 0;
			strand_specific = false;
		}
		else if (strcmp(argv[i], "-orf") == 0) 
		{

            not_defined = 0;
			find_orf = true;
		}
		else if (strcmp(argv[i], "-noorf") == 0) 
		{

            not_defined = 0;
			find_orf = false;
		}
		else if (strcmp(argv[i], "-orflen") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -orflen\n") ;
               //usage();
                exit(1);
            }
			i++; 
			min_orf_len = atoi(argv[i]);
		}
		    else if (strcmp(argv[i], "-orfsep") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -orfsep\n") ;
               //usage();
                exit(1);
            }
			i++; 
			min_orf_sep = atof(argv[i]);
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
