#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <algorithm>
	using std::reverse;
#include "genome.h"
#include "region.h"
#include "read.h"
#include "splice_labels_from_RNA_seq.h"
#include <string>
	using std::string;
#include <vector>
	using std::vector;
#include "gene.h"
#include "gene_tools.h"

void Gene::find_orf(int min_len, float separation)
{
	char* mRNA_seq; 
	int len; 
	get_mRNA_seq(&mRNA_seq, &len); 

	//printf("strand: %c, len: %i \n%s\n", strand, len, mRNA_seq);

	int tis = 0;
	int stop = 0;
	int second_best = 0;
	GeneTools::find_max_orf(mRNA_seq, len, &tis, &stop, &second_best);

	//printf("tis:%i, stop:%i, len:%i second:%i\n", tis, stop, stop-tis, second_best);

	int orflen = stop-tis;
	if (orflen<min_len || ((double) second_best)/orflen>separation)
	{
		//printf("rejected: len:%i second:%i, frac: %1.2f\n", stop-tis, second_best, second_best/((double) orflen));
		return;
	}
	string tis_cons(&mRNA_seq[tis], 3);
	//printf("tis_cons: %s\n", tis_cons.c_str()); 
	string stop_cons(&mRNA_seq[stop], 3);
	//printf("stop_cons: %s\n", stop_cons.c_str()); 


	int dna_tis;
	int dna_stop;
	if (strand=='-')
	{
		dna_tis = map_rna_to_dna(len-tis+1);
		dna_stop = map_rna_to_dna(len-stop-2);
	}
	else
	{
		dna_tis = map_rna_to_dna(tis+1);
		dna_stop = map_rna_to_dna(stop+4);
	}

	delete[] mRNA_seq;
	//printf("dna_tis: %i, dna_stop:%i\n", dna_tis-start, dna_stop-start);

	//string dna_tis_cons(&seq[dna_tis-start], 3);
	//printf("dna tis_cons: %s\n", dna_tis_cons.c_str()); 
	//string dna_stop_cons(&seq[dna_stop-start], 3);
	//printf("dna stop_cons: %s\n", dna_stop_cons.c_str()); 
	if (dna_tis>-1 && dna_stop>-1)
		split_exons(dna_tis, dna_stop);

}
void Gene::split_exons(int tis, int stop)
{
	utr5exons.clear();
	cds_exons.clear();
	utr3exons.clear();
	if (strand=='+')
	{
		for (int i=0; i<exons.size(); i++)
		{
			if (exons[i].second<tis)
			{
				utr5exons.push_back(exons[i]);
			}
			else if (exons[i].first<tis && exons[i].second<stop)
			{
				// exon contains tis, but not stop codon
				segment utr_part(exons[i].first, tis-1);
				utr5exons.push_back(utr_part);
				segment cds_part(tis, exons[i].second);//might be of length 1 if tis==exons[i].second
				cds_exons.push_back(cds_part);
			}
			else if (exons[i].first<tis && exons[i].second>=stop)
			{
				// exon contains tis and stop codon
				segment utr5part(exons[i].first, tis-1);
				utr5exons.push_back(utr5part);
				segment cds_part(tis, stop-1);
				cds_exons.push_back(cds_part);
				segment utr3part(stop, exons[i].second);
				utr3exons.push_back(utr3part);
			}
			else if (exons[i].first>=tis && exons[i].second<stop)
			{
				cds_exons.push_back(exons[i]);
			}
			else if (exons[i].first>=stop)
			{
				utr3exons.push_back(exons[i]);
			}
			else if (exons[i].first>=tis && exons[i].second>=stop)
			{
				segment cds_part(exons[i].first, stop-1);
				cds_exons.push_back(cds_part);
				segment utr3part(stop, exons[i].second);
				utr3exons.push_back(utr3part);
			}
			else
			{
				fprintf(stderr, "Error, this case is not captured\n");
			}
		}
	}
	else if (strand=='-')
	{
		for (int i=0; i<exons.size(); i++)
		{
			if (exons[i].second<stop)
			{
				utr3exons.push_back(exons[i]);
			}
			else if (exons[i].first<stop && exons[i].second<tis)
			{
				segment utr_part(exons[i].first, stop-1);
				utr3exons.push_back(utr_part);
				segment cds_part(stop, exons[i].second);
				cds_exons.push_back(cds_part);
			}
			else if (exons[i].first<stop && exons[i].second>=tis)
			{
				segment utr3part(exons[i].first, stop-1);
				utr3exons.push_back(utr3part);
				segment cds_part(stop, tis-1);
				cds_exons.push_back(cds_part);
				segment utr5part(tis, exons[i].second);
				utr5exons.push_back(utr5part);
			}
			else if (exons[i].first>=stop && exons[i].second<tis)
			{
				cds_exons.push_back(exons[i]);
			}
			else if (exons[i].first>=tis)
			{
				utr5exons.push_back(exons[i]);
			}
			else if (exons[i].first>=stop && exons[i].second>=tis)
			{
				segment cds_part(exons[i].first, tis-1);
				cds_exons.push_back(cds_part);
				segment utr5part(tis, exons[i].second);
				utr5exons.push_back(utr5part);
			}
			else
			{
				fprintf(stderr, "Error, this case is not captured c-strand\n");
			}
		}
	}
}

int Gene::map_rna_to_dna(int rna_pos)
{
	if (exons.size()==0)
	{
		fprintf(stderr, "no exons specified\n");
		return -1;
	}
	int cur_offset = exons[0].first;
	int mrna_len = 0;

	int exon_cnt = 0;
	int prev_mrna_len=0; 
	while (rna_pos>mrna_len && exon_cnt<exons.size())
	{
		prev_mrna_len = mrna_len;
		mrna_len += exons[exon_cnt].second-exons[exon_cnt].first+1;
		exon_cnt++;
	}
	if (exon_cnt>0)
		exon_cnt--;
	if (rna_pos<=mrna_len)
		return rna_pos-prev_mrna_len+exons[exon_cnt].first;
	else
		return -1;
}

void Gene::get_mRNA_seq(char** mRNA_seq, int* len)
{
	if (start>=stop)
	{
		printf("gene start >= gene stop\n");
		print(stdout);
	}
	char* genome_seq = this->get_sequence();
	
	string* sseq = new string("");
	for (int i=0; i<exons.size();i++)
	{
		int exon_len = exons[i].second-exons[i].first+1;
		int exon_start_local = exons[i].first-start;
		if (strand=='-'&&i>0)
			exon_start_local = exons[i].first-start-1;
		else if(strand=='-'&&i<exons.size()-1)
			exon_len = exon_len-1;
		string* exon_seq = new string(&genome_seq[exon_start_local], exon_len);
		*sseq += *exon_seq;
		delete exon_seq;
	}
	*len = (int) sseq->size();
	if (strand=='-')
	{
		reverse(sseq->begin(),(sseq->end()));
		gio->complement((char*) sseq->c_str(), *len);
	}
	*mRNA_seq = new char[sseq->size()+1];
	strcpy(*mRNA_seq, sseq->c_str());

	delete sseq;
}


bool Gene::write_window(_IO_FILE*& fd, string tmp_seq, int center_pos, int left_offset, int right_offset, int label)
{

    // calculate window coordinates
    int window_start = center_pos - left_offset;
    int window_stop = center_pos + right_offset;

    // check bounds
    if (window_start >= 0 && window_stop < get_length()) 
    {
        // extract string
        string window = tmp_seq.substr(window_start, window_stop - window_start);

        // write to file
        fprintf(fd, "%s %i %c\n", window.c_str(), label, strand);

        return true;
    }

    return false;

}


bool Gene::generate_tis_labels(vector<example_t>* examples)
{

    if (!is_coding()) 
    {
        return false;
    }
	char* tmp_seq = this->get_sequence();

    // find actual TIS
    int tis_pos = 0;
    int local_tis_pos = 0;

    // deal with strand weirdness
    if (strand=='+') 
    {
        tis_pos = cds_exons[0].first;
        local_tis_pos = cds_exons[0].first - start - 1;
    } 
    else
    {
        tis_pos = cds_exons.back().second;
        local_tis_pos = stop - cds_exons.back().second + 1;

		reverse(tmp_seq, tmp_seq+get_length());
		gio->complement(tmp_seq, get_length());
    }

	example_t exm;
	exm.pos = tis_pos;
	exm.label = 1;
	examples->push_back(exm);


    // create motif    
    vector<string*> tis_cons; 
    tis_cons.push_back(new string("atg", 3));

	assert(GeneTools::check_consensus(local_tis_pos, tmp_seq, get_length()-1, tis_cons));
	if (false)
	{
		printf("tis:");
		for (int i=-3;i<4;i++)
			printf(" %i:%i",i, GeneTools::check_consensus(local_tis_pos+i, tmp_seq, get_length()-1, tis_cons));
		printf("\n");
	}

    // find negative examples
    for (int i=0; i!=get_length()-5; i++)
    {
        if (GeneTools::check_consensus(i, tmp_seq, get_length()-1, tis_cons))
        {
            // accept, if it is not the positive example
            if (i != tis_pos)
            {
				if (strand=='+')
					exm.pos = i+start+1;
				else
					exm.pos = stop-i+1;
				exm.label = -1;
				examples->push_back(exm);
            }
        }
    }
    return true;
}


bool Gene::generate_tss_labels(vector<example_t>* examples)
{
	int win = 300;

    // init coordinates
    int tss_pos = 0;
	int left = 0;
	int right = 0;

    // deal with negative strand
    if (strand=='+') 
    {
        tss_pos = exons[0].first - 1;
		right = exons.back().second; 
		left = std::min(tss_pos+win, right);// some distance to the true tss site
    } 
    else
    {
		tss_pos = exons.back().second + 1;
		left = exons[0].first;
		right = std::max(tss_pos-win, left);
	}

	example_t exm;
	exm.pos = tss_pos;
	exm.label = 1;
	examples->push_back(exm);

	// number of negative examples
	int take=std::min(100, (right-left)/100); 

	
	srand(time(NULL));
    // find negative examples
    for (int i=0; i<take; i++)
    {
  		int neg_pos = left+rand()%(right-left);
		
		assert(neg_pos>=left);	
		assert(neg_pos<=right);	

		exm.pos = neg_pos;
		exm.label = -1;
		examples->push_back(exm);
    }
    return true;
}

