#include <algorithm>
#include "assert.h"
#include "region.h"
#include "genome.h"
#include "get_reads_direct.h"
#include "read.h"
#include "math_tools.h"
#include "graph_tools.h"
#include <vector>
  using std::vector;
#include <fstream>

bool compare_second(segment intr1, segment intr2);

/** default constructor*/
Region::Region()
{
	start = NULL; 
	stop = NULL;
	strand = NULL;
	chr_num = NULL;
	chr = NULL;
	seq = NULL;
	coverage = NULL;
	intron_coverage = NULL;
	gio = NULL;
	
}

/** constructor*/
Region::Region(int pstart, int pstop, int pchr_num, char pstrand, const char* gio_fname)
{
	start = pstart; 
	stop = pstop;
	strand = pstrand;
	chr_num = pchr_num;
	chr = NULL;
	seq = NULL;
	coverage = NULL;
	intron_coverage = NULL;
	fd_out = stdout;

	// initialize genome information object
	gio = new Genome(); 
	int ret = gio->init_genome((char*) gio_fname);
	if (ret<0)
	{   
		fprintf(stderr, "error reading genome info file: %s\n", gio_fname);
		return;
	}
}
/** constructor*/
Region::Region(int pstart, int pstop, int pchr_num, char pstrand)
{
	start = pstart; 
	stop = pstop;
	strand = pstrand;
	chr_num = pchr_num;
	chr = NULL;
	seq = NULL;
	coverage = NULL;
	intron_coverage = NULL;
	fd_out = stdout;
	gio = NULL;
}
/** constructor*/
Region::Region(int pstart, int pstop, char* pchr, char pstrand)
{
	start = pstart; 
	stop = pstop;
	strand = pstrand;
	chr_num = NULL;
	chr = pchr;
	seq = NULL;
	coverage = NULL;
	intron_coverage = NULL;
	fd_out = stdout;
	gio = NULL;
}
/** constructor*/
Region::Region(Region* reg)
{
	start = reg->start; 
	stop = reg->stop;
	strand = reg->strand;
	chr_num = reg->chr_num;
	if (reg->chr)
	{
		chr = new char[100];
		strcpy(chr, reg->chr);
	}
	else
	{
		chr = NULL;
	}
	seq = NULL;
	coverage = NULL;
	intron_coverage = NULL;
	fd_out = reg->fd_out;
	gio = NULL;
}

/** destructor*/
Region::~Region()
{
	delete[] seq;	
	
	clear_reads();

	delete[] chr;
	delete[] coverage;
	delete[] intron_coverage;
	intron_list.clear();
	unique_introns.clear();
	intron_counts.clear();
}

void Region::clear_reads()
{
	for (int i=0; i<all_reads.size(); i++)
		delete all_reads[i];
	all_reads.clear();
	// reads is a subset of all_reads, therefore the 
	// destructor for each read has already been called
	// only the pointers have to be removed
	reads.clear(); 
}

void Region::set_gio(Genome* pgio)
{
	gio = pgio;
}

void Region::load_genomic_sequence()
{
	if (!check_region())
	{
		printf("load_genomic_sequence: check_region failed\n");
		print(stderr);
		exit(-1);
	}
	seq = gio->read_flat_file(chr_num, start, stop);
}

void Region::print(_IO_FILE*& fd)
{
	fprintf(fd, "region %s\n", get_region_str());
	fprintf(fd, "region start:\t%i\n", start);
	fprintf(fd, "region stop:\t%i\n", stop);
	fprintf(fd, "region strand:\t%c\n", strand);
	if (chr_num)
		fprintf(fd, "region chr_num:\t%i\n", chr_num);
	if (gio && chr_num)
		fprintf(fd, "region chr:\t%s\n", gio->get_contig_name(chr_num));
	else if (chr)
		fprintf(fd, "region chr:\t%s\n", chr);
}

char* Region::get_region_str()
{
	char* reg_str = new char[1000];
	if (chr)
		sprintf(reg_str, "%s:%i-%i", chr, start, stop);
	else if (gio && chr_num)
		sprintf(reg_str, "%s:%i-%i", gio->get_contig_name(chr_num), start, stop);
	else
	{
		fprintf(stderr, "genome information object not set");
		exit(-1);
	}
	return reg_str;
}

void Region::get_reads(char** bam_files, int num_bam_files, int intron_len_filter, int filter_mismatch, int exon_len_filter, bool mm_filter)
{
	char* reg_str = get_region_str();
	int subsample = 1000;
	for (int i=0; i<num_bam_files; i++)
	{
    	char* fn_bam = bam_files[i];
	    fprintf(fd_out, "getting reads from file: %s\n", fn_bam);
		get_reads_from_bam(fn_bam, reg_str, &all_reads, strand, subsample);
	}
	delete[] reg_str; 
	fprintf(fd_out, "number of reads (not filtered): %d\n", (int) all_reads.size());
	/* filter reads
	* **************/
	//int exon_len_filter = 8;
	//int filter_mismatch = 1;
	//int intron_len_filter = 100000;
	int intron_filter = 0;
	int exon_filter = 0;
	int mismatch_filter = 0;
	int multi_filter= 0;
	for (int i=0; i<all_reads.size(); i++)
	{
		if (all_reads[i]->max_intron_len()>=intron_len_filter)
			intron_filter++;
		if (all_reads[i]->min_exon_len()<=exon_len_filter)
			exon_filter++;
		if (all_reads[i]->get_mismatches()>filter_mismatch)
			mismatch_filter++;
		if (all_reads[i]->multiple_alignment_index!=0 && mm_filter)
			multi_filter++;

	    if (all_reads[i]->max_intron_len()<intron_len_filter && all_reads[i]->min_exon_len()>exon_len_filter && all_reads[i]->get_mismatches()<=filter_mismatch && (all_reads[i]->multiple_alignment_index==0 || !mm_filter))
	        reads.push_back(all_reads[i]);
	}
	fprintf(fd_out, "number of reads: %d\n", (int) reads.size());
	fprintf(fd_out, "Filter: intron_len:%i, min_exon_len:%i, mismatch:%i, multi_mapper:%i\n", intron_filter, exon_filter, mismatch_filter, multi_filter);
}

void Region::compute_coverage() 
{
	int num_pos = stop-start+1;
	if (!coverage)
		coverage = new uint32_t[num_pos];

	memset(coverage, 0, num_pos*sizeof(uint32_t));
	//for (int i=0; i<num_pos; i++)
	//	coverage[i] = 0;
	for (int i=0; i<reads.size(); i++)
	{
		// add read contribution to coverage
		reads[i]->get_coverage(start, stop, coverage);
	}
}

float Region::get_coverage_global(int pstart, int pstop)
{
	if (!coverage)
		compute_coverage();
	
	assert(pstart<=pstop);

	int sum = 0;
	for (int i = pstart; i<=pstop; i++)
	{
		if (!(i>=start && i<=stop))
		{
			fprintf(stdout, "start%i pstart%i stop%i pstop%i\n",start, pstart, stop, pstop);	
		}
		assert(i>=start && i<=stop);
		sum+=coverage[i-start];
	}
	return ((float) sum)/(pstop-pstart+1);
}

float Region::get_coverage_seg(int i)
{
	assert(i<segments.size());
	if (seg_cov.size()==0)
	{
		compute_seg_cov();
	}
	return seg_cov[i];
}

void Region::compute_seg_cov()
{
	seg_cov.clear();
	for (int s=0; s<segments.size(); s++)
		seg_cov.push_back(get_coverage_global(segments[s].first, segments[s].second));
}

void Region::compute_intron_coverage() 
{
	int num_pos = stop-start+1;
	if (!intron_coverage)
		intron_coverage = new uint32_t[num_pos];

	for (int i=0; i<num_pos; i++)
		intron_coverage[i] = 0;
	for (int i=1; i<unique_introns.size(); i++)
	{
		int from_pos = std::max(start, unique_introns[i].first);
        int to_pos = std::min(stop, unique_introns[i].second);
		for (int j=from_pos; j<to_pos; j++)
		{
			intron_coverage[j] += intron_counts[i];
		}
	}
}

int Region::get_intron_conf(int intron_start, int intron_stop)
{
	//printf("get_intron_conf for %i -> %i\n", intron_start, intron_stop);
	for (int i=1; i<unique_introns.size(); i++)
	{
		//if (unique_introns[i].first==intron_start || unique_introns[i].second==intron_stop)
		//printf("found intron %i -> %i\n", unique_introns[i].first, unique_introns[i].second);
		if (unique_introns[i].first==intron_start && unique_introns[i].second==intron_stop)
			return intron_counts[i];
	}
	return 0;
}

int Region::get_read_starts(int from, int to)
{
	int ret = 0;
	for (int i=0; i<reads.size(); i++)
		if (reads[i]->start_pos>=from && reads[i]->start_pos<=to)
			ret++;
	return ret;
}
int Region::get_read_ends(int from, int to)
{
	int ret = 0;
	for (int i=0; i<reads.size(); i++)
		if (reads[i]->get_last_position()>=from && reads[i]->get_last_position()<=to)
			ret++;
	return ret;
}
void Region::find_tss_and_cleave(vector<int>* pos, vector<int>* starts, vector<int>* stops, float pval)
{
	// process read starts and ends
	vector<int> rstarts;
	vector<int> rends;
	for (int i=0; i<reads.size(); i++)
	{
		rstarts.push_back(reads[i]->start_pos);
		rends.push_back(reads[i]->get_last_position());
	}
	sort(rstarts.begin(), rstarts.end());
	sort(rends.begin(), rends.end());

	int* rs = &rstarts[0];
	int* re = &rends[0];

	int num_pos = pos->size();
	for (int i=0; i<num_pos-1; i++)
	{
		if (pos->at(i)==pos->at(i+1))
			continue;
		
		int win=30;
		vector<int> bins;
		vector<int> start_cnts;
		vector<int> end_cnts;
		for (int j=pos->at(i); j<pos->at(i+1)-win; j+=win)
		{
			bins.push_back(j);
			start_cnts.push_back(0);
			end_cnts.push_back(0);
			while (rs<&rstarts.back() && *rs<j+win)
			{
				rs++;
				if (*rs>=j)
					start_cnts.back()++;
			}
			//printf("%i, %i\n", start_cnts.back(), get_read_starts(j+1, j+win));
			while (re<&rends.back() && *re<j+win)
			{
				re++;
				if (*re>=j)
					end_cnts.back()++;
			}
			//printf("%i, %i\n", end_cnts.back(), get_read_ends(j+1, j+win));
		}
		for (int b=0; b<((int) bins.size())-1; b++)
		{
			//int cnt1 = get_read_starts(j, j+win-1);
			//int cnt2 = get_read_starts(j+win, j+2*win-1);

			//int cnt3 = get_read_ends(j, j+win-1);
			//int cnt4 = get_read_ends(j+win, j+2*win-1);

			int cnt1 = start_cnts[b];
			int cnt2 = start_cnts[b+1];

			int cnt3 = end_cnts[b];
			int cnt4 = end_cnts[b+1];

			float p_val1 = bino_cfd(cnt1, cnt1+cnt2, 0.5);
			float p_val2 = bino_cfd(cnt4, cnt3+cnt4, 0.5);
			if (p_val1<pval)
			{
				//printf("add tss%i: cnt1:%i, cnt2:%i, pval:%.6f\n", j+win, cnt1, cnt2, p_val1);
				pos->push_back(bins[b]+win);
				starts->push_back(bins[b]+win);
			}
			else if (p_val2<pval)
			{
				//printf("add cleave%i: cnt1:%i, cnt2:%i, pval:%.6f\n", j+win, cnt1, cnt2, p_val2);
				pos->push_back(bins[b]+win);
				stops->push_back(bins[b]+win);
			}
		}
	}
}

void Region::generate_segment_graph(float seg_filter, float tss_pval)
{
	compute_intron_list();
	vector<int> pos;
	vector<int> starts;
	vector<int> stops;

	// add splice sites, tss, and tts from transcripts
	for (int i=0; i<transcripts.size(); i++)
	{
		starts.push_back(transcripts[i].front().first-1);
		stops.push_back(transcripts[i].back().second);
		for (int j=0; j<transcripts[i].size(); j++)
		{
			pos.push_back(transcripts[i][j].first-1);
			pos.push_back(transcripts[i][j].second);
		}
	}
	//printf("pos from anno:");
	//for (int i=0; i<pos.size(); i++) printf("%i, ", pos[i]);printf("\n");

	if (reads.size()==0 && pos.size()==0)
	{
		segment seg(start, stop);
		segments.push_back(seg);
	}
	else
	{
		pos.push_back(start-1);
		for (int i=0; i<unique_introns.size(); i++)
		{
			// this is how introns are defined:
			//       ^
			//     /   \
			//    |     |
			//NNNNGC...AGNNN
			if (unique_introns[i].first>start && unique_introns[i].first<stop)
				pos.push_back(unique_introns[i].first-1);// last exonic position
			if (unique_introns[i].second>start && unique_introns[i].second<stop)
				pos.push_back(unique_introns[i].second);// the following segment starts with the first exonic position
		}

		//printf("pos from spliced reads: "); 
		//for (int i=0; i<pos.size(); i++) printf("%i, ", pos[i]); printf("\n");

		// add segment boundaries where coverage drops to zero
		compute_coverage();
		int num_pos = stop-start+1;
		for (int i=0; i<num_pos-1; i++)
		{
			if (coverage[i]==0 && coverage[i+1]>0)
				pos.push_back(i+1+start-1);// start next segment with the first exonic position
			else if (coverage[i]>0 && coverage[i+1]==0)
				pos.push_back(i+start);// this is the last exonic position
		}

		sort(pos.begin(), pos.end());
		pos.push_back(stop);
		
		//for (int i=0; i<pos.size(); i++) printf("%i, ", pos[i]);
		//printf("\n");

		find_tss_and_cleave(&pos, &starts, &stops, tss_pval);

		sort(pos.begin(), pos.end());

		//for (int i=0; i<pos.size(); i++) printf("%i, ", pos[i]);
		//printf("\n");
		
		//for (int i=0; i<starts.size(); i++) printf("%i, ", starts[i]);
		//printf("\n");
		//for (int i=0; i<stops.size(); i++) printf("%i, ", stops[i]);
		//printf("\n");

		
		// create segments
		for (int i=0; i<((int)pos.size())-1; i++)
		{
			if (pos[i]==pos[i+1])
				continue;
			segment seg(pos[i]+1, pos[i+1]);
			segments.push_back(seg);
			//printf("%i->%i\n", pos[i]+1, pos[i+1]);
		}

		// filter segments without coverage
		vector<segment> seg_filtered;
		for (int i=0; i<segments.size(); i++)
		{
			if (!(segments[i].first<=segments[i].second))
			{
				
				print_segments_and_coverage(stdout); 
				assert(false);
			}
			float cov = 0;
			for (int j=segments[i].first; j<segments[i].second; j++)
				if (coverage[j-start]>0)
					cov++;

			cov = cov/(segments[i].second-segments[i].first+1);
			//float cov = get_coverage_global(segments[i].first, segments[i].second);
			//printf("%i->%i  %.4f, %.4f\n", segments[i].first, segments[i].second, cov, seg_filter);
			
			// check if segment is part of any transcript
			bool keep = is_annotated(i);


			// check if removal of this segment cuts a path between to 
			// paired end reads
			if (pair_mat.size()>0)
			{
				compute_admat(starts, stops);

				// remove all connections to the current segment
				vector<vector<float> > admat_wo = admat;
				for (int j=0; j<admat.size(); j++)
				{
					admat_wo[i][j] = NO_CONNECTION;
					admat_wo[j][i] = NO_CONNECTION;
				}


				for (int j=0; j<pair_mat.size(); j++)
				{
					for (int k=0; k<pair_mat.size(); k++)
					{
						if (pair_mat[j][k]>0)
						{
							int np1 = count_num_paths(admat, j+1, k+1);
							int np2 = count_num_paths(admat_wo, j+1, k+1);
							if (np2==0 && np1>0)
							{
								// had only one connection via the current segment
								keep = true;
							}
						}
					}
				}
			}


			if (cov>seg_filter || keep) 
			{
				seg_filtered.push_back(segments[i]);	
			}
		}
		//printf("\t%i (%i) segments passed filter\n", (int)seg_filtered.size(), (int) segments.size());
		if (seg_filtered.size()>0)
		{
			segments.clear();
			segments = seg_filtered;
		}
	}
	compute_admat(starts, stops);


	// remove short tss and tts segment next to splice sites
	// those originate mose likely from wrong alignments of 
	// reads that should be spliced
	for (int i=0; i<segments.size(); i++)
	{
		bool discard = false;
		if (is_annotated(i))
			continue;

		if (segments[i].second-segments[i].first>20)
			continue;

		if (is_initial(i+1))
		{
			if (i+1<segments.size() && is_acceptor_ss(i+1+1))// +1 because admat is 1 based
			{
				discard = true;
			}
		}
		else if(is_terminal(i+1))
		{
			if (i>0 && is_donor_ss(i-1+1))// +1 because admat is 1 based 
			{
				discard = true;
			}
		}

		if (discard)
		{
			vector<segment> seg_filtered;
			for (int j=0; j<segments.size(); j++)
			{
				if (j!=i)
					seg_filtered.push_back(segments[j]);
			}
			segments = seg_filtered;
			compute_admat(starts, stops);
		}
	}

	// compute transcript paths
	for (int i=0; i<transcripts.size(); i++)
	{
		vector<int> path; 
		int seg = 0; 
		while (segments[seg].first!=transcripts[i][0].first)
		{
			seg++;
			assert(seg<segments.size());
		}
		assert(seg<segments.size());

		path.push_back(seg);
		int exon = 0;
		while (exon<transcripts[i].size())
		{
			assert(segments[seg].first==transcripts[i][exon].first);
			while (segments[seg].second != transcripts[i][exon].second)
			{
				assert(seg+1<segments.size());
				// make sure segments are connected
				assert(segments[seg].second+1==segments[seg+1].first);
				seg++;
				//add segment to the path
				path.push_back(seg);
			}

			exon++;
			if (exon<transcripts[i].size())
			{
				// find segment
				while (segments[seg].first!=transcripts[i][exon].first)
				{
					seg++; 
					assert(seg<segments.size());
				}
				assert(segments[seg].first==transcripts[i][exon].first);
				path.push_back(seg);
				// assert that intron connection exists
			}
		}
		transcript_paths.push_back(path);
		// print
		//for (int exon=0; exon<transcripts[i].size(); exon++)
		//	printf("exon: %i->%i\n", transcripts[i][exon].first, transcripts[i][exon].second );

		//for (int seg=0; seg<segments.size(); seg++)
		//	printf("segment: %i->%i\n", segments[seg].first, segments[seg].second);

		//printf("path:  ");
		//for (int p=0; p<path.size(); p++)
		//	printf("%i ", path[p]);
		//printf("\n");
			
	}
}
bool Region::is_annotated(int i)
{
	bool ret = false;
	for (int k=0; k<transcripts.size(); k++)
		for (int j=0; j<transcripts[k].size(); j++)
			if (transcripts[k][j].first<=segments[i].first && transcripts[k][j].second>=segments[i].second)
				ret = true;
	return ret;
}
bool Region::is_acceptor_ss(int i)
{
	bool ret = false;
	for (int j=i; j<segments.size(); j++)
	{
		if (admat[i][j]>NO_CONNECTION && segments[i-1].second<segments[j-1].first-1)
			ret = true;
	}
	return ret;
}
bool Region::is_donor_ss(int i)
{
	bool ret = false;
	for (int j=1; j<i; j++)
	{
		if (admat[j][i]>NO_CONNECTION && segments[j-1].second<segments[i-1].first-1)
			ret = true;
	}
	return ret;
}
bool Region::is_initial(int i)
{
	assert(i<admat[0].size());
	return admat[0][i]>NO_CONNECTION || i==0;
}
bool Region::is_terminal(int i)
{
	int num_seg = segments.size();
	assert(i<admat[num_seg].size());
	return admat[num_seg+1][i]>NO_CONNECTION || i==num_seg;
}
vector<int> Region::find_max_path()
{
	// find path with highest expression
	int max_seg = -1;
	float max_cov = -1;
	for (int i=0; i<segments.size(); i++)
	{
		float cov = get_coverage_global(segments[i].first, segments[i].second);
		if (cov>max_cov)
		{
			max_cov=cov;
			max_seg=i;
		}
	}
	vector<int> path;
	path.push_back(max_seg);


	// search upstream
	int i=max_seg;
	while (!is_initial(i+1))
	{
		vector<int> parents = get_parents(i);	
		assert(parents.size()>0);
		float max_cov = -3;
		int next = -1;
		for (int j=0; j<parents.size(); j++)
		{
			int k = parents[j];
			float cov=-1;
			if (admat[i][k]==NEIGHBOR)
			{
				float ci = get_coverage_seg(i);
				float ck = get_coverage_seg(k);
				cov = std::min(ci, ck);
			}
			else if (admat[i][k]==NO_CONNECTION)
			{
				fprintf(stderr, "This is a bug: no connection from %i->%i\n", i, k);
			}
			else
			{
				cov = admat[i][k];
			}
			if (cov>max_cov)
			{
				max_cov = cov;
				next = j;
			}
		}
		assert(next>-1);
		i = next;
		fprintf(stdout, "%i  ", next);
		path.push_back(next);
	}
	// search downstream
	i=max_seg;
	while (!is_terminal(i))
	{
		vector<int> children = get_children(i);	
		float max_cov = -1;
		int next = -1;
		for (int j=0; j<children.size(); j++)
		{
			int k = children[j];
			float cov=-1;
			if (admat[i][k]==NEIGHBOR)
			{
				float ci = get_coverage_seg(i);
				float ck = get_coverage_seg(k);
				cov = std::min(ci, ck);
			}
			else
			{
				cov = admat[i][k];
			}
			if (cov>max_cov)
			{
				max_cov = cov;
				next = j;
			}
		}
		i = next;
		fprintf(stdout, "%i  ", next);
		path.push_back(next);
	}
	// sort path
	sort(path.begin(), path.end());

	return path;
}

void Region::add_bias_counts(vector<int>* vec)
{
	int num_bins = vec->size();
	vector<float> exonic_cov;

	vector<int> path = find_max_path();
	//printf("\n");
	for (int i=0; i<path.size(); i++)
	{
		int j=path[i];
		for (int p=segments[j].first; p<segments[j].second; p++)
		{
			float cov = get_coverage_global(p, p);
			exonic_cov.push_back(p);
	//		printf("%.1f ", cov);
		}
	}
	//printf("\n");

	int exonic_len = exonic_cov.size();
	float step = exonic_len/num_bins;

}

void Region::print_segments_and_coverage(_IO_FILE*& fd)
{
	for (int i=0; i<segments.size(); i++) 
	{
		float cov = get_coverage_global(segments[i].first, segments[i].second);
		fprintf(fd, "%i\t%i\t%.4f\n", segments[i].first, segments[i].second, cov);
	}
}

void Region::compute_admat(vector<int> starts, vector<int> stops)
{
	// check that segments are sorted
	for (int i=1; i<segments.size(); i++)
		assert(segments[i].first>=segments[i-1].second);

	// allocate memory and set all entries to NO_CONNECTION
	int num_seg = segments.size();
	admat.resize(num_seg+2); // source and sink node
	for (int i=0; i<admat.size(); i++)
	{
		admat[i].resize(num_seg+2);
		for (int j=0; j<admat.size(); j++)
			admat[i][j] = NO_CONNECTION;
	}

	// add connections from RNA-seq introns
	update_coverage_information();

	// connect neighboring segments
	for (int j=1; j<segments.size(); j++)
	{
		if (segments[j-1].second+1==segments[j].first)
		{
			admat[j][j+1] = NEIGHBOR;
			admat[j+1][j] = NEIGHBOR;
			//printf("connecting neighbor (%i->%i) %i, %i\n", segments[j-1].second, segments[j].first, j, j+1);
		}
		//else
			//printf("not connecting neighbor (%i->%i) %i, %i\n", segments[j-1].second, segments[j].first, j, j+1);
	}
	// connect nodes to source and sink
	for (int j=0; j<segments.size(); j++)
	{
		vector<int> parents = get_parents(j+1);	
		vector<int> children = get_children(j+1);	
		bool is_start=false;
		for (int k=0; k<starts.size(); k++)
			if(segments[j].first==starts[k]+1)
				is_start = true;
		bool is_stop=false;
		for (int k=0; k<stops.size(); k++)
			if(segments[j].second==stops[k])
				is_stop = true;

		//printf("%i->%i, start%i, stop%i parents:%i, children:%i\n", segments[j].first, segments[j].second, is_start, is_stop, (int)parents.size(), (int)children.size());

		if (parents.size()==0 || is_start)
		{
			admat[0][j+1] = CONNECTION;
			admat[j+1][0] = CONNECTION;
		}
		if (children.size()==0 || is_stop)
		{
			admat[num_seg+1][j+1] = CONNECTION;
			admat[j+1][num_seg+1] = CONNECTION;
		}
	}
}

void Region::update_coverage_information()
{
	compute_intron_list();

	// add connections from transcripts
	vector<segment> trans_introns;
	for (int i=0; i<transcripts.size(); i++)
	{
		int num_exons = transcripts[i].size();
		for (int j=0; j<num_exons-1; j++)
		{
			segment* seg = new segment();
			// intron coordinates should be: 
			//       |-------|
			//XXXXXXX---------XXXXXXX
			// exon1   intron  exon2
			seg->first = transcripts[i][j].second+1;
			seg->second = transcripts[i][j+1].first-1;
			trans_introns.push_back(*seg);
		}
	}
	// match transcript introns with RNA-seq introns
	int match_cnt = 0;
	for (int i=0; i<trans_introns.size(); i++)
	{
		bool matched = false;
		for (int j=0; j<unique_introns.size(); j++)
		{
			if (trans_introns[i].first==unique_introns[j].first && trans_introns[i].second==unique_introns[j].second)
				matched = true;
		}
		if (!matched)
		{
			unique_introns.push_back(trans_introns[i]);
			intron_counts.push_back(0);
		}
		else
			match_cnt++;
	}
	if (trans_introns.size()>5 && unique_introns.size()>5 && match_cnt<2)
	{
		printf("Warning: annotated introns do not match RNA-seq introns\n");
		for (int i=0; i<trans_introns.size(); i++)
		{
			printf("%i->%i\n", trans_introns[i].first, trans_introns[i].second);
		}
		for (int i=0; i<unique_introns.size(); i++)
		{
			printf("%i->%i (%i)\n", unique_introns[i].first, unique_introns[i].second, intron_counts[i]);
		}
	}

	// remove previous coverage intormation
	for (int i=0; i<admat.size(); i++)
		for (int j=0; j<admat.size(); j++)
			if (admat[i][j]>=CONNECTION)
				admat[i][j] = CONNECTION;

	// add counts for introns
	for (int i=0; i<unique_introns.size(); i++)
	{
		int start = unique_introns[i].first-1;
		int stop = unique_introns[i].second+1;
		int seg1=-1;
		int seg2=-1;
		for (int j=0; j<segments.size(); j++)
		{
			if (segments[j].second==start)
			{
				seg1 = j+1;// one based because of source node
			}
			if (segments[j].first==stop)
			{
				seg2 = j+1;// one based because of source node
			}
		}
		//assert(seg1>=0 && seg2>=0);
		if (!(seg1>=0 && seg2>=0))
			continue;
		//{
		//	printf("intron missing: %i->%i\n", start, stop);
		//}
		//printf("intron: %i->%i %i->%i (%i)\n", start, stop, seg1, seg2, intron_counts[i]);
		if (intron_counts[i]==0)
		{
			assert(admat[seg1][seg2]<=0);
			assert(admat[seg2][seg1]<=0);
		}
		admat[seg2][seg1] = intron_counts[i];
		admat[seg1][seg2] = intron_counts[i];
	}
}
int Region::compute_pair_mat()
{
	// get read pairs
	sort(reads.begin(), reads.end(), CRead::compare_by_read_id);

	vector<int> pair_ids;
	int num_pos = stop-start+1;
	bool no_pair_filter = true;
	
	int take_cnt = 0;
	int discard_cnt = 0; 
	int discard_strand_cnt = 0; 
	//printf("reads.size(): %i\n", reads.size());
	// find consecutive reads with the same id
	for (int i=0; i<((int) reads.size())-1; i++) 
	{
		//printf("reads[%i]->read_id: %s\n", i, reads[i]->read_id);
		int j = i+1;
		while(j<reads.size() && strcmp(reads[i]->read_id, reads[j]->read_id) == 0) 
		{
			if ((reads[i]->left && reads[j]->right) || (reads[j]->left && reads[i]->right) && (reads[i]->reverse != reads[j]->reverse)) 
			{
				if (reads[i]->get_last_position()==-1 || reads[j]->get_last_position()==-1)
					break;
				if (false)//(reads[i]->strand[0]=='0' && reads[j]->strand[0]=='0' )
				{ 
					// discard pairs without strand information
					discard_strand_cnt++;
				}
				else if (reads[i]->get_last_position()<reads[j]->start_pos && reads[j]->start_pos-reads[i]->get_last_position()<60000) 
				{
					pair_ids.push_back(i);
					pair_ids.push_back(j);
					take_cnt++;
				}
				else if (reads[i]->start_pos>reads[j]->get_last_position() && reads[j]->get_last_position()-reads[i]->start_pos<60000)
				{
					pair_ids.push_back(i);
					pair_ids.push_back(j);
					take_cnt++;
				}
				else
				{
					if (no_pair_filter>0 && reads[i]->start_pos<reads[j]->start_pos)
					{
						pair_ids.push_back(i);
						pair_ids.push_back(j);
						take_cnt++;
					}
					else if (no_pair_filter>0)
					{
						pair_ids.push_back(j);
						pair_ids.push_back(i);
						take_cnt++;
					}
					else
						discard_cnt++;
					//printf("istart:%i, ilast:%i  jstart:%i, jlast: %i\n", reads[i]->start_pos, reads[i]->get_last_position(), reads[j]->start_pos, reads[j]->get_last_position());
				}
			}
			else
				discard_cnt++;
			j++;
		}
	}

	//pair_mat = zeros(size(segments, 1));
	pair_mat.clear();
	for (int i=0; i<segments.size(); i++)
	{
		vector<int> row(segments.size(), 0);
		pair_mat.push_back(row);
	}

	//for p = 1:size(pairs, 2)
	for (int i=0; i<pair_ids.size(); i+=2)
	{
	//	if read_stops(pairs(1, p)+1)<read_stops(pairs(2, p)+1)
	//		r1 = pairs(1, p)+1;
	//		r2 = pairs(2, p)+1;
	//	else
	//		r1 = pairs(2, p)+1;
	//		r2 = pairs(1, p)+1;
	//	end
		int rid1 = pair_ids[i];
		int rid2 = pair_ids[i];
		if (reads[rid1]->get_last_position()>reads[rid2]->get_last_position())
		{
			// switch reads
			int tmp = rid2;
			rid2 = rid1;
			rid1 = tmp;
		}
		int st1 = -1;
		int en1 = -1;
		int st2 = -1;
		int en2 = -1;
		for (int j=0; j<segments.size(); j++)
		{
			if (reads[rid1]->start_pos>=segments[j].first && reads[rid1]->start_pos<=segments[j].second)
				st1 = j;
			if (reads[rid1]->get_last_position()>=segments[j].first && reads[rid1]->get_last_position()<=segments[j].second)
				en1 = j;
			if (reads[rid2]->start_pos>=segments[j].first && reads[rid2]->start_pos<=segments[j].second)
				st2 = j;
			if (reads[rid2]->get_last_position()>=segments[j].first && reads[rid2]->get_last_position()<=segments[j].second)
				en2 = j;
		}
		if (st1>-1 && st2>-1 && en1>-1 && en2>-1)
		{
			pair_mat[st1][st2]++;
			pair_mat[st2][st1]++;

			pair_mat[en1][en2]++;
			pair_mat[en2][en1]++;

			pair_mat[st1][en2]++;
			pair_mat[en2][st1]++;

			pair_mat[en1][st2]++;
			pair_mat[st2][en1]++;
		}
	//	% find the segments where the left read starts and right read ends
	//	% and increase pair counter between these segments
	//	for x = 1:4
	//		if x==1 % start <-> stop
	//			s1 = find(segments(:,1)<=read_starts(r1)&segments(:,2)>=read_starts(r1));
	//			s2 = find(segments(:,1)<=read_stops(r2)&segments(:,2)>=read_stops(r2));
	//		elseif x==2 % start <-> start
	//			s1 = find(segments(:,1)<=read_starts(r1)&segments(:,2)>=read_starts(r1));
	//			s2 = find(segments(:,1)<=read_starts(r2)&segments(:,2)>=read_starts(r2));
	//		elseif x==3 % stop <-> stop
	//			s1 = find(segments(:,1)<=read_stops(r1)&segments(:,2)>=read_stops(r1));
	//			s2 = find(segments(:,1)<=read_stops(r2)&segments(:,2)>=read_stops(r2));
	//		elseif x==4 % stop <-> start
	//			s1 = find(segments(:,1)<=read_stops(r1)&segments(:,2)>=read_stops(r1));
	//			s2 = find(segments(:,1)<=read_starts(r2)&segments(:,2)>=read_starts(r2));
	//		end

	//		if isempty(s1) || isempty(s2)
	//			continue
	//		end

	//		assert(length(s1)==1)
	//		assert(length(s2)==1)
	//		
	//		pair_mat(s1, s2) = pair_mat(s1, s2)+1;
	//		if s1~=s2
	//			pair_mat(s2, s1) = pair_mat(s2, s1)+1;
	//		end
	//	end
	//end
	}
}


vector<int> Region::get_parents(int node)
{
	vector<int> parents;
	for (int j=1; j<node; j++)
	{
		if (admat[j][node]>=NEIGHBOR)
			parents.push_back(j);
	}
	return parents;
}
vector<int> Region::get_children(int node)
{
	vector<int> children;
	for (int j=node+1; j<=segments.size(); j++)
	{
		if (admat[j][node]>=NEIGHBOR)
			children.push_back(j);
	}
	return children;
}

void Region::write_segment_graph(_IO_FILE*& fd)
{
	fprintf(fd, "region\t%s\t%c\t%i\t%i\n", chr, strand, start, stop);

	int num_seg = segments.size();
	fprintf(fd, "segments\t%i\n", num_seg);
	print_segments_and_coverage(fd);
	for (int i=0; i<admat.size(); i++)
	{
		for (int j=0; j<admat.size(); j++)
		{
			if (admat[i][j]>NEIGHBOR)
				fprintf(fd, "%i\t%i\n", i+1, j+1);
		}
	}
	fprintf(fd, "end\n");
}


int Region::write_binary(std::ofstream* ofs)
{
	// region meta info
	ofs->write((char *) &start, sizeof(int));
	ofs->write((char *) &stop, sizeof(int));
	int len = strlen(chr);
	assert(len>0);
	ofs->write((char *) &len, sizeof(int));
	ofs->write((char *) chr, len);
	ofs->write((char *) &strand, sizeof(char));

	// splice graph
	len = segments.size();
	ofs->write((char *) &len, sizeof(int));
	ofs->write((char *) &segments[0], len*sizeof(segment));
	if (seg_cov.size()==0)
		compute_seg_cov();
	//if (len==1 && seg_cov[0]<1e-3)
	//{
	//	printf("seg_cov for single segment locus: %.2f\n", seg_cov[0]);
	//	printf("number of reads (filtered): %i\n", (int)reads.size());
	//}
	ofs->write((char *) &seg_cov[0], len*sizeof(float));
	
	//splice graph admat
	len = admat.size();
	ofs->write((char *) &len, sizeof(int));
	vector<int> sp;
	vector<float> val;
	for (int i=0; i<len; i++)
		for (int j=i+1; j<len; j++)
			if (admat[i][j]>NEIGHBOR)
			{
				sp.push_back(i);
				sp.push_back(j);
				val.push_back(admat[i][j]);
			}
	len = sp.size();
	ofs->write((char *) &len, sizeof(int));
	ofs->write((char *) &sp[0], len*sizeof(int));
	len = val.size();
	ofs->write((char *) &len, sizeof(int));
	ofs->write((char *) &val[0], len*sizeof(float));


	// transcript paths
	int num_trans = transcript_paths.size();
	ofs->write((char *) &num_trans, sizeof(int));
	for (int i=0; i<num_trans; i++)
	{
		len = transcript_paths[i].size();
		ofs->write((char *) &len, sizeof(int));
		if (len>0)
		{
			ofs->write((char *) &transcript_paths[i][0], len*sizeof(int));
		}
	}

	// pair matrix
	len = pair_mat.size();
	ofs->write((char *) &len, sizeof(int));
	if (len==0)
		return -1;
	sp.clear();
	vector<int> pair_cnt;
	for (int i=0; i<len; i++)
		for (int j=i+1; j<len; j++)
			if (pair_mat[i][j]>0)
			{
				sp.push_back(i);
				sp.push_back(j);
				pair_cnt.push_back(pair_mat[i][j]);
			}
	len = sp.size();
	ofs->write((char *) &len, sizeof(int));
	ofs->write((char *) &sp[0], len*sizeof(int));
	len = pair_cnt.size();
	ofs->write((char *) &len, sizeof(int));
	ofs->write((char *) &pair_cnt[0], len*sizeof(int));

	return 1;
}
int Region::read_binary(std::ifstream* ifs)
{
	// region meta info
	ifs->read((char *) &start, sizeof(int));
	ifs->read((char *) &stop, sizeof(int));
	int len = 0;
	ifs->read((char *) &len, sizeof(int));
	if (len<=0)
	{
		fprintf(stderr, "read_binary: could not read from file: len:%i\n", len);
		return -1;
	}
	chr = new char[len+1];
	ifs->read((char *) chr, len);
	chr[len] = '\0';
	ifs->read((char *) &strand, sizeof(char));

	//printf("%s:%i..%i %c\n", chr, start, stop, strand);
	// splice graph
	ifs->read((char *) &len, sizeof(int));
	assert(len>0);
	segment seg[len];
	ifs->read((char *) seg, len*sizeof(segment));
	for (int i=0; i<len; i++)
	{
		segments.push_back(seg[i]);
	}
	seg_cov.resize(len);
	ifs->read((char *) &seg_cov[0], len*sizeof(float));

	//splice graph admat
	ifs->read((char *) &len, sizeof(int));
	assert(len>0);
	assert(admat.size()==0);
	for (int i=0; i<len; i++)
	{
		vector<float> row(len, -2);
		admat.push_back(row);
	}
	ifs->read((char *) &len, sizeof(int));
	assert(len>0);
	int sp[len];
	ifs->read((char *) sp, len*sizeof(int));
	int len2 = 0;
	ifs->read((char *) &len2, sizeof(int));
	assert(len2*2==len);
	float val[len2];
	ifs->read((char *) val, len2*sizeof(float));
	for (int i=0; i<len2; i++)
	{
		int j = sp[2*i];
		int k = sp[2*i+1];
		admat[j][k] = val[i];
		admat[k][j] = val[i];
	}
	
	fd_out = stdout;

	// transcript paths
	int num_trans = 0;
	ifs->read((char *) &num_trans, sizeof(int));
	for (int i=0; i<num_trans; i++)
	{
		ifs->read((char *) &len, sizeof(int));
		if (len>0)
		{
			vector<int> path(len, 0);
			ifs->read((char *) &path[0], len*sizeof(int));

			for (int j=0; j<path.size(); j++)
				assert(path[j]<segments.size());
			transcript_paths.push_back(path);
		}
	}

	//pair_mat
	ifs->read((char *) &len, sizeof(int));
	if (len==0)
		return 0;

	pair_mat.clear();
	for (int i=0; i<segments.size(); i++)
	{
		vector<int> row(segments.size(), 0);
		pair_mat.push_back(row);
	}
	ifs->read((char *) &len, sizeof(int));
	int pairs[len];
	ifs->read((char *) pairs, len*sizeof(int));
	ifs->read((char *) &len, sizeof(int));
	int pair_cnt[len];
	ifs->read((char *) pair_cnt, len*sizeof(int));

	for (int i=0; i<len; i+=2)
	{
		int j = pairs[i]; 
		int k = pairs[i+1];
		pair_mat[j][k] = pair_cnt[i/2];
	}
	return 0;
}
void Region::compute_intron_list()
{
	intron_list.clear();
	unique_introns.clear();
   	intron_counts.clear();

	vector<int> introns;
	for (int i=0; i<reads.size(); i++)
	{
		reads[i]->get_introns(&introns);
    }
	for (int i=0; i<introns.size(); i+=2)
	{
		segment intr(introns[i], introns[i+1]);
		intron_list.push_back(intr);
	}
	fprintf(fd_out, "found %i introns\n", (int) intron_list.size());
	
	// sort by intron start
	sort(intron_list.begin(), intron_list.end(), compare_second);

	if (intron_list.size()>0)
	{
    	unique_introns.push_back(intron_list[0]);
		intron_counts.push_back(1);
    	for (int i=1; i<intron_list.size(); i++)
    	{
			//printf("intron: %i->%i\n", intron_list[i].first, intron_list[i].second);
    	    if (intron_list[i].first!=intron_list[i-1].first || intron_list[i].second!=intron_list[i-1].second)
    	    {
				//printf("uintron: %i->%i\n", unique_introns.back().first, unique_introns.back().second);
    	       	unique_introns.push_back(intron_list[i]);
				intron_counts.push_back(1);
    	    }
			else
			{
    	    	intron_counts.back()++;
			}
    	}
	}

	//for (int i=1; i<unique_introns.size(); i++)
	//	assert(unique_introns[i].first!=unique_introns[i-1].first || unique_introns[i].second!=unique_introns[i-1].second);

	//vector<segment*> same_start;
	//for (int i=0; i<intron_list.size(); i++)
	//{
	//	if (i<intron_list.size()-1 && intron_list[i].first==intron_list[i+1].first)
	//	{
	//		// collect introns with same start
	//		same_start.push_back(&intron_list[i]);
	//	}
	//	else
	//	{
	//		same_start.push_back(&intron_list[i]);
	//		sort(same_start.begin(), same_start.end(), compare_second);
	//		unique_introns.push_back(*(same_start[0]));
	//		intron_counts.push_back(1);
	//		for (int j=1; j<same_start.size(); j++)
	//		{
	//			if (same_start[j]->second!=same_start[j-1]->second)
	//			{
	//				unique_introns.push_back(*(same_start[j]));
	//				intron_counts.push_back(1);
	//			}
	//			else
	//			{
	//				//printf("intron_counts before: %i ", (*(intron_counts.back())));
	//				(intron_counts.back())++; 
	//				//printf("after: %i \n", (*intron_counts.back()));
	//			}
	//		}
	//		same_start.clear();
	//	}
	//}
	fprintf(fd_out, "found %i unique introns\n", (int) unique_introns.size());

#ifdef READ_DEBUG
	//check unique list
	printf("DEBUG: Check 'unique introns'-list contains all introns\n");
	int idx = 0;
	for (int i=0; i<intron_list.size(); i++)
	{
		bool match = false;
		while (unique_introns[idx].first<intron_list[i].first && idx<unique_introns.size())
			idx++;
		int idx2 = idx;
		while (unique_introns[idx2].first==intron_list[i].first && idx2<unique_introns.size())
		{
			if (unique_introns[idx2].second==intron_list[i].second)
				match = true;
			idx2++;
		}
		assert(match);
	}
#endif

	int sum = 0;
	for (int i=0; i<intron_counts.size(); i++)
		sum+=intron_counts[i]; 
	//if (sum!=intron_list.size())
	//{
	//	printf("Error %i!=%i\n", sum, (int)intron_list.size());
	//	for (int i=0; i<intron_list.size(); i++)
	//	{
	//		printf("intron: %i->%i\n", intron_list[i].first, intron_list[i].second);
	//	}
	//	for (int i=0; i<unique_introns.size(); i++)
	//	{
	//		printf("uintron: %i->%i\t%i\n", unique_introns[i].first, unique_introns[i].second, intron_counts[i]);
	//	}
	//}
	assert(sum==intron_list.size());
}



bool compare_second(segment intr1, segment intr2)
{
	if (intr1.first<intr2.first)
		return true;
	if (intr1.first>intr2.first)
		return false;
	return (intr1.second<intr2.second);
}


