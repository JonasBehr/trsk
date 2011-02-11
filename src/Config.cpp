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
	max_intron_len = 20000;
	mm_filter = 1;
	el_filter = 8;

	max_exon_len = 8000;
	min_exon_len = 10;
	exon_mean = 5;
	exon_term_thresh = 3;
	exon_drop = 5;
	exon_cut = 3;
	intron_conf = 1;
	intron_dist = 0;
	intron_seed_conf = 3;
	reject_retained_introns = false;
	term_filter = 2.0;

	find_orf = true;
	min_orf_len = 300;
	min_orf_sep = 0.7;

	region_rel_length = 0.25;
	min_intergenic_len = 50;
	max_intergenic_len = 20000;
	intergenic_win = 100;
	strand_specific = false;

	have_bam_file = false;
	have_gio_file = false;
	have_reg_file = false; 
	have_gff_file = false;

}
int Config::parseCommandLine(int argc, char *argv[])
{
    bool not_defined;

    for (int i = 1; i < argc; i++) 
	{
        not_defined = true;

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
        else if (strcmp(argv[i], "-mm") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -mm\n") ;
               //usage();
                exit(1);
            }
			i++; 
			mm_filter = atoi(argv[i]);
		}
        else if (strcmp(argv[i], "-el") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -el\n") ;
               //usage();
                exit(1);
            }
			i++; 
			el_filter = atoi(argv[i]);
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
        else if (strcmp(argv[i], "-indt") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -indt\n") ;
               //usage();
                exit(1);
            }
			i++; 
			intron_dist = atoi(argv[i]);
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
        else if (strcmp(argv[i], "-exm") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -exm\n") ;
               //usage();
                exit(1);
            }
			i++; 
			exon_mean = atoi(argv[i]);
		}
        else if (strcmp(argv[i], "-exts") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -exts\n") ;
               //usage();
                exit(1);
            }
			i++; 
			exon_term_thresh = atoi(argv[i]);
		}
        else if (strcmp(argv[i], "-exd") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -exd\n") ;
               //usage();
                exit(1);
            }
			i++; 
			exon_drop = atoi(argv[i]);
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
		else if (strcmp(argv[i], "-incut") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -incut\n") ;
               //usage();
                exit(1);
            }
			i++; 
			intron_cut = atof(argv[i]);
		}
		else if (strcmp(argv[i], "-ri") == 0) 
		{

            not_defined = 0;
			reject_retained_introns = true;
		}
		else if (strcmp(argv[i], "-tf") == 0) 
		{

            not_defined = 0;
            if (i + 1 > argc - 1) 
			{
                fprintf(stderr, "ERROR: Argument missing for option -tf\n") ;
               //usage();
                exit(1);
            }
			i++; 
			term_filter = atof(argv[i]);
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

		if (not_defined)
		{
			fprintf(stderr, "ERROR: unknown argument %s\n", argv[i]) ;
			print(stdout);
			exit(1);
		}
	}
}
