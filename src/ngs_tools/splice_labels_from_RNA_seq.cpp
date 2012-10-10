#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <stdint.h>
#include <algorithm>
#include "get_reads_direct.h"
#include "genome.h"

#include <vector>
  using std::vector;
#include <string>
  using std::string;
#include "splice_labels_from_RNA_seq.h"

int get_splice_labels(int argc, char* argv[])
{

	bool write_examples = false;
	bool write_spf = false;
	FILE* write_fd = stdout ;
	
	int argn = 1;
	if ( strcmp(argv[1], "-v")==0 )
	{
		// write training exampltes to std out
		write_examples = true;
		argn++;
	}
	if ( strcmp(argv[1], "-s")==0 )
	{
		// write spf files to std out
		write_spf = true;
		argn++;
		char* write_fn = argv[argn++];
		write_fd = fopen(write_fn, "w+") ;
	}
		
	char* fn_genome_config = argv[argn++];
	char* dn_output = argv[argn++];//directory name
	char* fn_bam = argv[argn];
	int subsample = 1000;

	fprintf(stdout, "genome_config: %s\n", fn_genome_config);
	fprintf(stdout, "dn_output: %s\n", dn_output);
	fprintf(stdout, "fn_bam: %s\n", fn_bam);

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

	//for (int chr_num=0; chr_num<1; chr_num++)
	for (int chr_num=0; chr_num<g->num_contigs; chr_num++)
	{
		for (int strand=0; strand<2; strand++)
		{
			char* DNA_seq = g->read_flat_file(chr_num);
			int seq_len = g->contig_len(chr_num);
			char s;
			if (strand==0)
			{
				s = '+';
				acc.consensus[0] = (char*) "ag"; 
				acc.offset = 2; 
				don.consensus[0] = (char*) "gc"; 
				don.consensus[1] = (char*) "gt"; 
				don.offset = -1; 
			}
			else
			{
				s = '-';
				g->complement(DNA_seq, seq_len);
				acc.consensus[0] = (char*) "ga";
				acc.offset = -1;
				don.consensus[0] = (char*) "tg";
				don.consensus[1] = (char*) "cg";
				don.offset = 2;
			}

			char region[MAXLINE];
			sprintf(region, "%s:%i-%i", g->get_contig_name(chr_num), 1, seq_len);
			//sprintf(region, "%s:%i-%i", g->get_contig_name(chr_num), 1, 1000000);
			
			vector<CRead*> all_reads;
			for (int i=argn; i<argc; i++)
			{
				fn_bam = argv[i];
				fprintf(stdout, "getting reads from file: %s\n", fn_bam);
  				get_reads_from_bam(fn_bam, region, &all_reads, s, subsample);
			}

			fprintf(stdout, "number of reads (not filtered): %d\n", (int) all_reads.size());
			/* filter reads
			* **************/	
			int exon_len_filter = 15;
			int filter_mismatch = 0;
			int intron_len_filter = 20000;
			vector<CRead*> reads;
			for (uint32_t i=0; i<all_reads.size(); i++) 
			{
				if (all_reads[i]->max_intron_len()<intron_len_filter && all_reads[i]->min_exon_len()>exon_len_filter && all_reads[i]->get_mismatches()<=filter_mismatch && all_reads[i]->multiple_alignment_index==0)
					reads.push_back(all_reads[i]);
			}


			fprintf(stdout, "number of reads: %d\n", (int) reads.size());

			/* retrieve all positive examples from rna_seq data
			 * these will be removed from the negative class later
			 * *********************************************/
			vector<int> all_acc_pos; 
			vector<int> all_don_pos;

			for (uint32_t i=0; i<all_reads.size(); i++)
			{
				vector<int> tmp;
				all_reads[i]->get_acc_splice_sites(&tmp);
				for (uint32_t j=0; j<tmp.size(); j++)
					all_acc_pos.push_back(tmp[j]);
				tmp.clear();
				all_reads[i]->get_don_splice_sites(&tmp);
				for (uint32_t j=0; j<tmp.size(); j++)
					all_don_pos.push_back(tmp[j]);
			}
			vec_unique(&all_acc_pos);
			vec_unique(&all_don_pos);

			/* retrieve filtered positive examples from rna_seq data
			 * *********************************************/
			vector<int> acc_pos; 
			vector<int> acc_pos_noncon; 
			vector<int> don_pos;
			vector<int> don_pos_noncon;

			for (uint32_t i=0; i<reads.size(); i++)
			{
				vector<int> tmp;
				reads[i]->get_acc_splice_sites(&tmp);
				for (uint32_t j=0; j<tmp.size(); j++)
					if (check_consensus(tmp[j], acc, DNA_seq))
						acc_pos.push_back(tmp[j]);
					else
						acc_pos_noncon.push_back(tmp[j]);
				tmp.clear();
				reads[i]->get_don_splice_sites(&tmp);
				for (uint32_t j=0; j<tmp.size(); j++)
					if (check_consensus(tmp[j], don, DNA_seq))
						don_pos.push_back(tmp[j]);
					else
						don_pos_noncon.push_back(tmp[j]);
			}

			cons_hist(acc_pos_noncon, acc, DNA_seq);
			cons_hist(don_pos_noncon, don, DNA_seq);

			// make unique
			vec_unique(&acc_pos);
			vec_unique(&don_pos);
			fprintf(stdout, "number of acc sites: %d\n", (int) acc_pos.size());
			fprintf(stdout, "number of don sites: %d\n", (int) don_pos.size());
		

			/* retrieve negative examples from rna_seq data
			 * *********************************************/
			vector<int> acc_neg; 
			vector<int> don_neg;
		
			int mindist = 10;
			
			for (uint32_t i=0; i<reads.size(); i++)
			{
				for (uint32_t j=0; j<reads[i]->block_starts.size(); j++)
				{
					if (reads[i]->block_lengths[j]<2*mindist)
						continue;
					int start_pos = reads[i]->start_pos + reads[i]->block_starts[j] + mindist;
					int seq_len = reads[i]->block_lengths[j]-2*mindist;
					find_all_occurences(&acc_neg, &DNA_seq[start_pos], seq_len, acc, start_pos);
					find_all_occurences(&don_neg, &DNA_seq[start_pos], seq_len, don, start_pos);
				}
			}

			vec_unique(&acc_neg);
			vec_unique(&don_neg);
			fprintf(stdout, "number of decoy acc sites: %d\n", (int) acc_neg.size());
			fprintf(stdout, "number of decoy don sites: %d\n", (int) don_neg.size());

			//fprintf(stdout, "remove intersection between positive and negative examples\n");
			//remove_intersection(&acc_pos, &acc_neg);
			//remove_intersection(&don_pos, &don_neg);
			fprintf(stdout, "remove positive examples form negative class\n");
			remove_intersection(&all_acc_pos, &acc_neg);
			remove_intersection(&all_don_pos, &don_neg);

			fprintf(stdout, "number of acc sites: %d\n", (int) acc_pos.size());
			fprintf(stdout, "number of don sites: %d\n", (int) don_pos.size());
			fprintf(stdout, "number of decoy acc sites: %d\n", (int) acc_neg.size());
			fprintf(stdout, "number of decoy don sites: %d\n", (int) don_neg.size());

			fprintf(stdout, "\n");

			/* combine positive and negative examples, create label
			 * and sort
			 * *********************************************/
			vector<example_t> acc_exp;
			combine_with_label(&acc_exp, acc_pos, 1);
			combine_with_label(&acc_exp, acc_neg, -1);
			sort(acc_exp.begin(), acc_exp.end(), exp_compare);

			vector<example_t> don_exp;
			combine_with_label(&don_exp, don_pos, 1);
			combine_with_label(&don_exp, don_neg, -1);
			sort(don_exp.begin(), don_exp.end(), exp_compare);


			acc_pos.reserve(acc_exp.size());
			float* acc_label = new float[acc_exp.size()];
			don_pos.reserve(don_exp.size());
			float* don_label = new float[don_exp.size()];
			if (strand==0)
			{
				for (uint32_t i=0; i<acc_exp.size(); i++)	
				{
					acc_label[i] = acc_exp[i].label;
					acc_pos[i] = acc_exp[i].pos-1;
				}
				for (uint32_t i=0; i<don_exp.size(); i++)	
				{
					don_label[i] = don_exp[i].label;
					don_pos[i] = don_exp[i].pos+2;
				}
			}
			else
			{
				for (uint32_t i=0; i<acc_exp.size(); i++)	
				{
					acc_label[i] = acc_exp[i].label;
					acc_pos[i] = acc_exp[i].pos +1;
				}
				for (uint32_t i=0; i<don_exp.size(); i++)	
				{
					don_label[i] = don_exp[i].label;
					don_pos[i] = don_exp[i].pos;
				}
			}
	
			/* save to file
			 * *********************************************/
			char score_name[] = "label";
			int num_scores = 1;
			char path_prefix[MAXLINE];

			create_dir(path_prefix, dn_output, chr_num, acc, s);
			save_score_pos(path_prefix, (char**) &score_name, &acc_label, num_scores, &acc_pos[0], acc_exp.size());

			create_dir(path_prefix, dn_output, chr_num, don, s);
			save_score_pos(path_prefix, (char**) &score_name, &don_label, num_scores, &don_pos[0], don_exp.size());
			
			if (write_examples)
			{
				for (uint32_t i=0; i<acc_exp.size(); i++)
				{
					char str_buff[142];	
					for (int k=0; k<141; k++)
						str_buff[k] = DNA_seq[acc_pos[i]-70+k];
					str_buff[141] = 0;
					fprintf(stdout, "#%1.f %s\n", acc_label[i], str_buff);
				}
			}
			if (write_spf)
			{
				for (uint32_t i=0; i<acc_exp.size(); i++)
				{
					fprintf(write_fd, "%s\tacc\tlabel\t%i\t%c\t%i\n", g->get_contig_name(chr_num), acc_pos[i], (strand==0||strand=='+')  ? '+' : '-', (int)acc_label[i]);
				}
				for (uint32_t i=0; i<don_exp.size(); i++)
				{
					fprintf(write_fd, "%s\tdon\tlabel\t%i\t%c\t%i\n", g->get_contig_name(chr_num), don_pos[i], (strand==0||strand=='+') ? '+' : '-', (int)don_label[i]);
				}
			}
			
			// clean up 
			for (uint32_t i=0; i<all_reads.size(); i++)
				delete all_reads[i];
			all_reads.clear();
			reads.clear();

			acc_pos.clear();
			don_pos.clear();
			acc_neg.clear();
			don_neg.clear();
			delete[] acc_label;
			delete[] don_label;
		}
	}
	if (write_fd!=stdout)
		fclose(write_fd) ;
	delete g;
	return 0;
}

bool check_consensus(int pos, signal_t sig, char* seq)
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

void cons_hist(vector<int> pos, signal_t sig, char* DNA_seq)
{
			int num = 5;
			int offsets[2*num+1];
			for (int x=-num; x<num+1; x++)
				offsets[x+num] = 0;
			int all = 0;
			for (uint32_t i=0; i<pos.size(); i++)
			{
				all++;
				for (int x=-num; x<num+1; x++)
					if (check_consensus(pos[i]+x, sig, DNA_seq))
						offsets[x+num]++;

			}
			fprintf(stdout, "\n%s, all: %i\n", sig.sig_name, all);
			for (int x=-num; x<num+1; x++)
				fprintf(stdout, "%i:%i\t", x, offsets[x+num]);
			fprintf(stdout, "\n");

}

void combine_with_label(vector<example_t>* all_exp, vector<int> pos, int label)
{
	for (uint32_t i=0; i<pos.size(); i++)
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
	char tmp[MAXLINE];
	sprintf(tmp, "%s", path_prefix);
	sprintf(path_prefix, "%scontig_%d%c", tmp, chr_num+1, strand);
}

void remove_intersection(vector<int>* vec1, vector<int>* vec2)
{
	uint32_t k = 0;
	// assuming entries are sorted
	for (uint32_t i=0; i<vec1->size(); i++)
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

void efficient_clear(vector<int>* vec, uint32_t pos)
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
bool exp_compare(example_t exp1, example_t exp2)
{ 
	return (exp1.pos<exp2.pos);
}

