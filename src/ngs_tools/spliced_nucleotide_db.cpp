#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <algorithm>
#include "get_reads_direct.h"
#include "genome.h"

#include <vector>
  using std::vector;
#include <string>
  using std::string;
#define MAXLINE 1000
typedef struct {
	char** consensus;	 
	int cons_len; 
	int offset; 
	int cons_num; 
	char* sig_name; 
} signal_t;
typedef struct {
	int pos;
	int label;
} example_t;

void vec_unique(vector<int>* vec);
void find_all_occurences(vector<int>* pos, char* seq, int seq_len, signal_t sig, int start_pos);
void efficient_clear(vector<int>* vec, int pos);
void remove_intersection(vector<int>* vec1, vector<int>* vec2);
bool exp_compare(example_t exp1, example_t exp2){ return (exp1.pos<exp2.pos);};
void combine_with_label(vector<example_t>* all_exp, vector<int> pos, int label);
void create_dir(char* path_prefix, char* dn_output, int chr_num, signal_t sig, char strand);
bool check_cons(int pos, signal_t sig, char* seq);
void cons_hist(vector<int> pos, signal_t sig, char* DNA_seq);
extern void save_score_pos(char* path_prefix, char** score_names, float** scores, int num_scores, int* pos, int num_pos);

int main(int argc, char* argv[])
{

	char* fn_bam = argv[1];
	char* fn_genome_config = argv[2];
	char* dn_output = argv[3];
	int subsample = 1000;

	Genome* g = new Genome();
	int ret = g->init_genome(fn_genome_config);	

	if (ret<0)
	{
		fprintf(stderr, "error reading genome info file: %s\n", fn_genome_config);
		return -1;
	}

	fprintf(stdout, "number of contigs: %d\n", g->num_contigs);

	signal_t acc;
	acc.sig_name = (char*) "acc"; 
	acc.cons_num = 1; 
	acc.consensus = new char*[acc.cons_num]; 
	acc.consensus[0] = (char*) "ag"; 
	acc.cons_len = 2; 
	acc.offset = 2; 

	signal_t don;
	don.sig_name = (char*) "don"; 
	don.cons_num = 2; 
	don.consensus = new char*[don.cons_num]; 
	don.consensus[0] = (char*) "gc"; 
	don.consensus[1] = (char*) "gt"; 
	don.cons_len = 2; 
	don.offset = -1; 

	//for (int chr_num=0; chr_num<g->num_contigs; chr_num++)
	for (int chr_num=0; chr_num<1; chr_num++)
	{
		for (int strand=0; strand<2; strand++)
		{
			char* DNA_seq = g->read_flat_file(chr_num);
			int seq_len = g->contig_len(chr_num);
			char s;
			if (strand==0)
				s = '+';
			else
			{
				s = '-';
				g->complement(DNA_seq, seq_len);
				//reverse(DNA_seq[0], DNA_seq[seq_len-1]);
			}

			char region[MAXLINE];
			sprintf(region, "%s:%i-%i", g->get_contig_name(chr_num), 1, seq_len/10);
			
			vector<CRead*> reads;
  			get_reads_from_bam(fn_bam, region, &reads, s, subsample);

			fprintf(stdout, "number of reads: %d\n", (int) reads.size());


			/* retrieve positive examples from rna_seq data
			 * *********************************************/
			vector<int> acc_pos; 
			vector<int> don_pos;

			for (int i=0; i<reads.size(); i++)
			{
				reads[i]->get_acc_splice_sites(&acc_pos);
				reads[i]->get_don_splice_sites(&don_pos);
			}

			if (strand==1)
			{
				don.consensus[0] = (char*) "tg";
				don.consensus[1] = (char*) "cg";
				don.offset = 2;
				acc.consensus[0] = (char*) "ga";
				acc.offset = -1;
			}
			fprintf(stdout, "number of acc sites: %d\n", (int) acc_pos.size());
			fprintf(stdout, "number of don sites: %d\n", (int) don_pos.size());
		
			vector<example_t> introns;

			for (int i=0; i<acc_pos.size(); i++)
			{
				example_t exp;
				exp.pos = acc_pos[i];
				exp.label = don_pos[i];
				introns.push_back(exp);
			}
			sort(introns.begin(), introns.end(), exp_compare);

			
			vector<example_t> unique_introns;
			for (int i=1; i<introns.size(); i++)
			{
				if (introns[i].label!=introns[i-1].label||introns[i].pos!=introns[i-1].pos)
					unique_introns.push_back(introns[i]);
			}

			cons_hist(acc_pos, acc, DNA_seq);
			cons_hist(don_pos, don, DNA_seq);

			for (int i=0; i<0; i++)
			{
				if (!check_cons(unique_introns[i].pos+1, acc, DNA_seq)||!check_cons(unique_introns[i].label, don, DNA_seq))
					continue;

				char str_buff[101];
				for (int k=0; k<50; k++)
					str_buff[k] = DNA_seq[unique_introns[i].label-47+k];
				for (int k=0; k<50; k++)
					str_buff[k+50] = DNA_seq[unique_introns[i].pos-2+k];
				str_buff[100]=0;
				fprintf(stdout, "%s\n", str_buff);
			}

			// clean up 
			for (int i=0; i<reads.size(); i++)
				delete reads[i];
			reads.clear();
			acc_pos.clear();
			don_pos.clear();
		}
	}
	

	delete g;
}

void cons_hist(vector<int> pos, signal_t sig, char* DNA_seq)
{
			int num = 5;
			int offsets[2*num+1];
			for (int x=-num; x<num+1; x++)
				offsets[x+num] = 0;
			int all = 0;
			for (int i=0; i<pos.size(); i++)
			{
				all++;
				for (int x=-num; x<num+1; x++)
					if (check_cons(pos[i]+x, sig, DNA_seq))
						offsets[x+num]++;

			}
			fprintf(stdout, "\n%s, all: %i\n", sig.sig_name, all);
			for (int x=-num; x<num+1; x++)
				fprintf(stdout, "%i:%i\t", x, offsets[x+num]);
			fprintf(stdout, "\n");

}

bool check_cons(int pos, signal_t sig, char* seq)
{
	bool match = false; 
	for (int c=0; c<sig.cons_num; c++)
	{   
	    bool cmatch = true;
	    for (int j=0; j<sig.cons_len; j++)
	        cmatch = cmatch && seq[pos-sig.offset+j]==sig.consensus[c][j];
	    match = match || cmatch;
	}
	if (match)
		return true;
	return false;
}
void combine_with_label(vector<example_t>* all_exp, vector<int> pos, int label)
{
	for (int i=0; i<pos.size(); i++)
	{
		example_t exp;
		exp.pos = pos[i];
		exp.label = label;
		all_exp->push_back(exp);
	}

}

void create_dir(char* path_prefix, char* dn_output, int chr_num, signal_t sig, char strand)
{
	sprintf(path_prefix, "%s/%s/", dn_output, sig.sig_name);
	struct stat st;
	if(stat(path_prefix,&st) != 0)
	{
		if (mkdir(path_prefix, 01777))
		{
			fprintf(stderr, "cannot create dir: %s\n", path_prefix );
			exit(-1);
		}
	}
	sprintf(path_prefix, "%scontig_%d%c", path_prefix, chr_num+1, strand);
}

void remove_intersection(vector<int>* vec1, vector<int>* vec2)
{
	int k = 0;
	// assuming entries are sorted
	for (int i=0; i<vec1->size(); i++)
	{
		while (k<vec2->size() && vec2->at(k)<vec1->at(i))
			k++;
		
		if (k>=vec2->size())
			break;

		if (vec2->at(k)==vec1->at(i))
		{
			efficient_clear(vec2, k);
			efficient_clear(vec1, i);
		}
	}
	vec_unique(vec1);
	vec_unique(vec2);
}

void efficient_clear(vector<int>* vec, int pos)
{
	// erase takes time to copy everything
	// therefore the element is overwritten with 
	// the preceding element and unique is later 
	// used to delete all at once
	if (pos<vec->size()-1)
		vec->at(pos) = vec->at(pos+1);
	else if (pos<vec->size())
		vec->pop_back();
}

void vec_unique(vector<int>* vec)
{
	vector<int>::iterator it;
	sort(vec->begin(), vec->end());
	it = unique(vec->begin(), vec->end());
	vec->resize(it-vec->begin());
}

void find_all_occurences(vector<int>* pos, char* seq, int seq_len, signal_t sig, int start_pos)
{   
    for (int i=0; i<seq_len-sig.cons_len+1; i++)
    {
        bool match = false; 
        for (int c=0; c<sig.cons_num; c++)
        {   
            bool cmatch = true;
            for (int j=0; j<sig.cons_len; j++)
                cmatch = cmatch && seq[i+j]==sig.consensus[c][j];
            match = match || cmatch;
        }
        if (match)
        {
            pos->push_back(i+sig.offset+start_pos);
        }
    }
}

