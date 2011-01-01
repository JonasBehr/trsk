#include <string>
	using std::string;
#include <vector>
	using std::vector;
#include "region.h"
#include "genome.h"
#include "gene.h"
#include "gene_tools.h"

bool GeneTools::check_consensus(int pos, char* seq, int len, vector<string*> motifs)
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
			int stop = gio->contig_len(i);
			Region* reg = new Region(start, stop, i, strands[s]);
			reg->set_gio(gio);
			regions.push_back(reg);
		}
	}
	return regions;
}
