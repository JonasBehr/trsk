#include <string>
	using std::string;
#include <vector>
	using std::vector;
#include <sys/stat.h>
#include <algorithm>
#include "region.h"
#include "genome.h"
#include "gene.h"
#include "gene_tools.h"
#include "splice_labels_from_RNA_seq.h"
//extern void save_score_pos(char* path_prefix, char** score_names, float** scores, int num_scores, int* pos, int num_pos);

bool GeneTools::check_consensus(int pos, const char* seq, int len, vector<string*> motifs)
{
	
    bool match = false;
    for (int j=0; j<motifs.size(); j++)
    {   
		if (match == true)// last motif matched
			return true;
        for (int k=0; k<motifs[j]->size(); k++)
        {   
            if (pos+k>=len)
            {   
				match = false;
				break;
            }
            if (seq[pos+k]==motifs[j]->at(k))
                match = true;
            else
            {   
                match = false;
                break;
            }
        }
    }
	return match;
}

void GeneTools::find_max_orf(char* seq, int len, int* tis, int* stop, int* second_best)
{
	vector<string*> tis_cons; 
	tis_cons.push_back(new string("atg", 3));
	vector<string*> stop_cons; 
	stop_cons.push_back(new string("tag", 3));
	stop_cons.push_back(new string("taa", 3));
	stop_cons.push_back(new string("tga", 3));

	int max_len = 0;
	*second_best = 0;

	for (int frame=0; frame<3; frame++)
	{
		int first_tis = -1; 
		int first_stop = -1;
		for (int i=frame; i<len; i+=3)
		{
			if (first_tis==-1 && check_consensus(i, seq, len, tis_cons))
			{
				first_tis = i;
			}
			if (first_tis!=-1 && check_consensus(i, seq, len, stop_cons))
			{
				first_stop = i;
				int len = first_stop-first_tis; 
				if (len>max_len)
				{
					*second_best = max_len;
					max_len = len; 
					*tis = first_tis;
					*stop = first_stop; 
				}
				else if (len>*second_best)
				{
					*second_best = len;
				}
				first_tis = -1;
			}
		}
	}
	delete tis_cons[0];
	delete stop_cons[0];
	delete stop_cons[1];
	delete stop_cons[2];
}
void GeneTools::find_all_orfs(char* seq, int len, vector<int>* tis, vector<int>* stop, int min_len)
{
	vector<string*> tis_cons; 
	tis_cons.push_back(new string("atg", 3));
	vector<string*> stop_cons; 
	stop_cons.push_back(new string("tag", 3));
	stop_cons.push_back(new string("taa", 3));
	stop_cons.push_back(new string("tga", 3));

	for (int frame=0; frame<3; frame++)
	{
		int first_tis = -1; 
		int first_stop = -1;
		for (int i=frame; i<len; i+=3)
		{
			if (first_tis==-1 && check_consensus(i, seq, len, tis_cons))
			{
				first_tis = i;
			}
			if (first_tis!=-1 && check_consensus(i, seq, len, stop_cons))
			{
				first_stop = i;
				int len = first_stop-first_tis; 
				if (len>min_len)
				{
					tis->push_back(first_tis);
					stop->push_back(first_stop); 
				}
				first_tis = -1;
			}
		}
	}
	delete tis_cons[0];
	delete stop_cons[0];
	delete stop_cons[1];
	delete stop_cons[2];
}

vector<Region*> GeneTools::init_regions(const char* gio_fname)
{
	Genome* gio = new Genome(); 
	gio->init_genome((char*) gio_fname);
	const char* strands = "+-";

	vector<Region*> regions; 
	for (int i=0; i<gio->num_contigs; i++)
	{
		for (int s=0; s<2; s++)
		{
			int start = 1; 
			int stop = gio->contig_len(i)-1;
			Region* reg = new Region(start, stop, i, strands[s]);
			reg->set_gio(gio);
			regions.push_back(reg);
		}
	}
	return regions;
}

bool GeneTools::write_tss_labels(vector<Gene*>* genes, Region* region, char* dirname)
{
	vector<example_t> examples;
	for (int r=0; r<genes->size(); r++)
	{
		if (genes->at(r)->chr_num==region->chr_num&&genes->at(r)->strand==region->strand)
			genes->at(r)->generate_tss_labels(&examples);
	}
	printf("got %i tss labels\n", (int) examples.size());

	// sort training examples by position
	sort(examples.begin(), examples.end(), exp_compare);

	char basename[1000];

	sprintf(basename, "%s/contig_%i%c", dirname, region->chr_num+1, region->strand);
	
	// create directory
	struct stat st;
    if(stat(dirname,&st) != 0)
    {   
        if (mkdir(dirname, 01777))
        {   
            fprintf(stderr, "cannot create dir: %s\n", dirname);
            exit(-1);
        }
    }

	int* pos = new int[examples.size()];
	float* label = new float[examples.size()];
    for (int i=0; i<examples.size(); i++)    
    {
       label[i] = examples[i].label;
       pos[i] = examples[i].pos;
    }


	printf("writing to file %s\n", basename);
	char* score_name = (char*) "label";
	save_score_pos(basename, &score_name, &label, 1, pos, examples.size());
	return true;
}
