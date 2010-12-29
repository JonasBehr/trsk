#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <algorithm>
	using std::reverse;
#include "genome.h"
#include "region.h"
#include "read.h"
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

	printf("strand: %c, len: %i \n%s\n", strand, len, mRNA_seq);


	int tis = 0;
	int stop = 0;
	int second_best = 0;
	GeneTools::find_max_orf(mRNA_seq, len, &tis, &stop, &second_best);

	printf("tis:%i, stop:%i, len:%i second:%i\n", tis, stop, stop-tis, second_best);

	int orflen = stop-tis;
	if (orflen<min_len || ((double) second_best)/orflen>separation)
	{
		printf("rejected: len:%i second:%i, frac: %1.2f\n", stop-tis, second_best, second_best/((double) orflen));
		return;
	}
	string tis_cons(&mRNA_seq[tis], 3);
	printf("tis_cons: %s\n", tis_cons.c_str()); 
	string stop_cons(&mRNA_seq[stop], 3);
	printf("stop_cons: %s\n", stop_cons.c_str()); 

	int dna_tis = map_rna_to_dna(tis);
	int dna_stop = map_rna_to_dna(stop);
	if (strand=='-')
	{
		dna_tis = map_rna_to_dna(len-tis);
		dna_stop = map_rna_to_dna(len-stop);
		dna_tis-=3;//move to start of consensus
		dna_stop-=3;
	}

	printf("dna_tis: %i, dna_stop:%i\n", dna_tis-start, dna_stop-start);

	string dna_tis_cons(&seq[dna_tis-start], 3);
	printf("dna tis_cons: %s\n", dna_tis_cons.c_str()); 
	string dna_stop_cons(&seq[dna_stop-start], 3);
	printf("dna stop_cons: %s\n", dna_stop_cons.c_str()); 

	split_exons(dna_tis, dna_stop);
}
void Gene::split_exons(int tis, int stop)
{
	if (strand=='+')
	{
		stop +=4; 
		tis +=1;
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
		stop +=1; 
		tis +=4;
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
	exon_cnt--;
	if (rna_pos<=mrna_len)
		return rna_pos-prev_mrna_len+exons[exon_cnt].first;
	else
		return -1;
}

void Gene::get_mRNA_seq(char** mRNA_seq, int* len)
{
	if (!seq)
		load_genomic_sequence();

	string* sseq = new string("");
	for (int i=0; i<exons.size();i++)
	{
		int exon_len = exons[i].second-exons[i].first+1;
		int exon_start_local = exons[i].first-start;
		if (strand=='-'&&i>0)
			exon_start_local = exons[i].first-start-1;
		else if(strand=='-'&&i<exons.size()-1)
			exon_len = exon_len-1;
		string* exon_seq = new string(&seq[exon_start_local], exon_len);
		*sseq += *exon_seq;
	}

	*mRNA_seq = (char*) sseq->c_str(); 
	*len = (int) sseq->size();
	if (strand=='-')
	{
		reverse(sseq->begin(),(sseq->end()));
		gio->complement((char*) sseq->c_str(), *len);
	}
}
