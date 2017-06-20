/*
Copyright (c) 2011 Expression Analysis / Gunjan Hariani, Erik Aronesty

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

$Id$
*/
const char * VERSION = "1.01 $Id$";

#include <ctype.h>
#include <stdio.h>

void usage( FILE * f ) {
  fprintf( f,
	   "\nUsage: fastq-stats [options] <fastq-file>\n\n"
	   "Version: %s\n" 
	   "\n"
	   "Produces lots of easily digested statistics for the files listed\n" 
	   "\n"
	   "Options\n"
	   "\n"
	   "-c     cyclemax: max cycles for which following quality stats are produced [35]\n"
	   "-w INT window: max window size for generating duplicate read statistics [2000000]\n"
	   "-d     debug: prints out debug statements\n"
	   "-D     don't do duplicate read statistics\n"
	   "-s INT number of top duplicate reads to display\n"
	   "-x FIL output fastx statistics (requires an output filename)\n"
	   "-b FIL output base breakdown by per phred quality at every cycle.\n"
	   "       It sets cylemax to longest read length\n"
	   "-L FIL Output length counts \n\n"
	   
	   "\n" 
	   "The following data are printed to stdout:\n" "\n"
	   "  reads			: #reads in the fastq file\n"
	   "  len 	                : read length. mean and stdev are provided for variable read lengths\n"
	   "  phred			: phred scale used\n"
	   "  window-size		: Number of reads used to generate duplicate read statistics\n"
	   "  cycle-max		: Number of bases to assess for duplicity\n"
	   "  dups			: Number of reads that are duplicates\n"
	   "  %%dup			: Pct reads that are duplcate\n"
	   "  unique-dup seq	: Number sequences that are duplicated\n"
	   "  min dup count		: Smallest duplicate tally for any duplicate sequence\n"
	   "  dup seq <rank> <count> <sequence> \n"
	   "  			: Lists top 10 most frequent duplicate reads along with count mean and stdev\n"
	   "  qual			: Base Quality min, max and mean\n"
	   "  %%A,%%T,%%C,%%G		: base percentages\n" 
	   "  total bases		: total number of bases\n" 
	   "\n"
	   ,VERSION);
  
} //end usage function

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <string>
#include <sparsehash/sparse_hash_map> // or sparse_hash_set, dense_hash_map, ...
#include <iostream>
#include "fastq-lib.h"
#include "gcModel.h"

using namespace std;

#define T_A 0
#define T_C 2
#define T_G 6
#define T_T 19
#define roundgt0(x) (long)(x+0.5)

class ent {
public:
    std::string seq;
    int cnt;

    ent(const std::string &s, int c) { seq=s; cnt=c; };
    static bool comp_cnt (const ent &a, const ent &b) {
        return a.cnt > b.cnt;
    };

};

class countPerCycle {
	public:
	int basecount[26];
	int qc;
	double qsum;
	//vector<int> qual; 
	int counts_by_qual[127];
	int qmin;
	int qmax;
	
	countPerCycle() {
		//26 english alphabets for A/C/G/T char
		for(int i=0; i<26; i++) {
			basecount[i]=0;
		}
		for(int i=0; i<127; i++) {
			counts_by_qual[i] = 0;
		}
		qc = 0;
		qsum = 0;
		qmin = 10000;
		qmax = 0;
	};
};

class count_perCycle_perQual {
	public:
	int counts_by_qual[127];

	count_perCycle_perQual() {
		for(int i=0; i<127; i++) {
			counts_by_qual[i] = 0;
		}
	};
};


void usage( FILE * f );
double std_dev( double count , double total, double sqsum );
double quantile( const std::vector <int> & vec, double p );
double quantiles_with_counts(int* v, int start, int end, double p, bool dbug);
std::string string_format( const std::string &fmt, ... );

extern int optind;
bool nodup = 0;
google::sparse_hash_map <std::string, int> dups;

vector <std::string> dup_reads; // do i need this

int window = 2000000;
int cyclemax = 35; 
int gcCyclemax = 100; // to compare with fastqc, seq is rounded to nearest 100 to reduce # of gc models; for < 200 length, this is teh same as max=100
float gcSum;
uint64_t gcTotal;

int show_max = 10;
bool debug = 0;
bool fastx = 0;
char *fastx_outfile = NULL;
bool brkdown = 0;
char *brkdown_outfile = NULL;
bool len_hist = 0;
vector<int> vlen; //all read lengths
char *lenhist_outfile = NULL;
bool gc = 0;
char *gc_outfile = NULL;

int main( int argc, char**argv ) {

	int index;
	char c;
	optind = 0;
	char *filename = NULL;

// bad change to working syntax... breaks things!
//	if(argc < 2) {usage(stdout); return 0;}

	while ( (c = getopt (argc, argv, "?DdL:g:x:b:c:w:s:h")) != -1) {
		switch (c) {
			case 'c': cyclemax = atoi(optarg); break;
			case 'D': ++nodup; break;
			case 'd': ++debug; break;
			case 'w': window = atoi(optarg); break;
			case 's': show_max = atoi(optarg); break;
			case 'x': fastx_outfile = optarg; ++fastx; break;
			case 'b': brkdown_outfile = optarg; ++brkdown; break;
			case 'L': ++len_hist; lenhist_outfile = optarg; break;
			case 'g': gc_outfile = optarg; ++gc; break;
			case 'h': usage(stdout); return 0;
			case '?':
					  if (!optopt) {
						  usage(stdout); return 0;
					  } else if(optopt && strchr("gbxcws", optopt)) {
					 // 		fprintf(stderr, "Option -%c requires an argument.\n", optopt);
					  } else {
					//	  fprintf (stderr, "Unknown option \n", optopt);
					  }
					  usage(stderr);
					  return 1;
		}
	}
	
	filename = argv[optind];

	int lenmax = 0;
	int lenmin = 100000000;
	double lensum = 0;
	double lenssq = 0;
	double nbase = 0;
	int qualmax = 0;
	int qualmin = 100000;
	double qualsum = 0;
	double qualssq = 0;
	int errs = 0;
	long long nreads = 0;
	int ndups = 0;
	double dupss = 0;
	bool fixlen = 0; //is fixed length
	FILE *file;
	struct fq newFq; meminit(newFq);
	bool isgz;
	vector<countPerCycle> qcStats (1);
	vector<count_perCycle_perQual> qcStats_by_qual (1);
	int phred = 64;
	double ACGTN_count[26];
	double total_bases = 0;


	for(int i=0; i<26; i++) {
		ACGTN_count[i] = 0;
	}
	dups.set_deleted_key("<>");

	if(debug) {
		cout << endl;
		cout << "Parameters: " << endl;
		printf("cyclemax: %d, window: %d, nodup: %d, debug: %d, showmax: %d, fastx: %d, outfile: %s, breakdown: %s, gc: %s\n",
		       cyclemax, window, nodup, debug, show_max, fastx, fastx_outfile, brkdown_outfile, gc_outfile);
		cout << endl;
	}

	if(gc) {
	  gcInit(gcCyclemax);
	}

	//read file
	file = filename ? gzopen(filename,"r",&isgz) : stdin;
	while(read_fq(file,nreads++,&newFq)) {

		if(newFq.seq.n != newFq.qual.n) {
			errs++;
		}
		
		if(nreads == 10000) {
			if(!std_dev((double)nreads,lensum,lenssq)) {
				fixlen = 1;
			}
		}
		
		total_bases += newFq.seq.n;
		if(len_hist) {
			if(newFq.seq.n > vlen.size()) 
				vlen.resize(newFq.seq.n+1);
			++vlen[newFq.seq.n];
		}

		if(!fixlen) {
			if(newFq.seq.n > lenmax) {
				lenmax = newFq.seq.n;
			}
			if(newFq.seq.n < lenmin) {
				lenmin = newFq.seq.n;
			}
			lensum += newFq.seq.n;
			lenssq += newFq.seq.n*newFq.seq.n;
		}
	
		
		if((newFq.seq.n > qcStats.size()) && (fastx)) {
			qcStats.resize(newFq.seq.n,countPerCycle());

		} 
		
		if((newFq.seq.n > qcStats_by_qual.size()) && (brkdown) && (!fastx)) {
			qcStats_by_qual.resize(newFq.seq.n,count_perCycle_perQual());
		}

		int gcTally = 0;
		//compute quality stats for the first cyclemax bases
		for(int i=0; i < newFq.seq.n; i++) {
			int ascii_val = (int) newFq.qual.s[i];
			if(fastx && ((nreads < window) || (nreads%10 == 0))) {
				qcStats[i].qc++;
				qcStats[i].counts_by_qual[ascii_val]++;
				qcStats[i].qsum += ascii_val;
				qcStats[i].basecount[(toupper(newFq.seq.s[i])-65)]++;
				if(ascii_val < qcStats[i].qmin) {
					qcStats[i].qmin = ascii_val;
				}
				if(ascii_val > qcStats[i].qmax) {
					qcStats[i].qmax = ascii_val;
				}
			}
			if(brkdown && (!fastx) && ((nreads < window) || (nreads%10 == 0))) {
				qcStats_by_qual[i].counts_by_qual[ascii_val]++;
			}

			if (i < cyclemax) {
				nbase++;

				if(ascii_val > qualmax) {
					qualmax = ascii_val;
				}
				if(ascii_val < qualmin) {
					qualmin = ascii_val;
				}
				qualsum += ascii_val;
				qualssq += ascii_val*ascii_val;

				ACGTN_count[(toupper(newFq.seq.s[i])-65)]++;
			}
			
			if (gc && i < gcCyclemax) {
			  if(toupper(newFq.seq.s[i]) == 'G' || toupper(newFq.seq.s[i]) == 'C') {
			    gcTally++;
			  }
			}
		}
		if(gc) {
		  int gcReadLength = newFq.seq.n > gcCyclemax? gcCyclemax : newFq.seq.n;
		  gcProcessSequence(gcReadLength, gcTally);
		  gcSum += (float)( gcTally )/gcReadLength;
		  gcTotal++;
		}

		if(!nodup) {//if you want to look at duplicate counts
			if(newFq.seq.n > cyclemax) {
				newFq.seq.s[cyclemax] = '\0';
				newFq.seq.n = cyclemax;
			}

			if(nreads < window) {
				dups[newFq.seq.s]++;
			} else {
				if(dups.find(newFq.seq.s) != dups.end()) {
					dups[newFq.seq.s]++;
				}//make sure the element already exists in the key
	
				if(nreads==window) {
					google::sparse_hash_map<string,int>::iterator it = dups.begin();
					while(it != dups.end()) {
						if((*it).second <= 1) {
							dups.erase(it);
						}
						it++;
					} //end while loop
				}
			}//if nreads > window
		} //end if you want to look for dups
	
	} //end reading all fastq reads

	nreads--;

	int inputReadError = gzclose(file, isgz);


	if(gc) {
	  FILE *myfile;
	  myfile = fopen(gc_outfile, "w");
	  gcPrintDistribution(myfile);
	  gcClose();
	}

	std::vector<ent> dup_sort;
	google::sparse_hash_map<string,int>::iterator it = dups.begin();
	while(it != dups.end()) {
		if((*it).second > 1) {
			ent e((*it).first,(*it).second);
			//printf("seq: %s dups:%d\n", e.seq.c_str(), e.cnt);
			dup_sort.push_back(e);
			ndups += (*it).second;
			dupss += (*it).second*(*it).second;
		} 
		it++;
	} //end while loop
	dups.clear();

	std::sort(dup_sort.begin(),dup_sort.end(),ent::comp_cnt);
	
	if(nreads < window) {
		window = nreads;
	}

	if(nreads < 1) {
		cout << "No reads in " << filename << ", not generating output" << endl;
		return 0;
	}
	//autodetect phred
	if(qualmin < 64) {
		phred = 33;
	}
	printf("reads\t%lld\n",nreads);

	if(!fixlen) {
		printf("len\t%d\n", lenmax);
		printf("len mean\t%.4f\n", (double)lensum/nreads);
		if(nreads > 1) {
			printf("len stdev\t%.4f\n", std_dev((double)nreads,lensum,lenssq));
		}
		printf("len min\t%d\n", lenmin);
	} else {
		printf("len\t%d\n",lenmax);
	}
	
	printf("phred\t%d\n", phred);
	if(errs > 0) {
		printf("errors\t%d\n", errs);
	}


	printf("window-size\t%d\n", window);
	printf("cycle-max\t%d\n", cyclemax);

	if(fastx) {

		if(brkdown) {
			FILE *myfile;
			myfile = fopen(brkdown_outfile,"wd");
			fprintf(myfile,"Cycle\tQuality\tCount\n");
		
			for(int i=0; i<qcStats.size(); i++) {
				for(int j=qcStats[i].qmin; j<=qcStats[i].qmax; j++) {
					fprintf(myfile,"%d\t%d\t%d\n",(i+1),(j-phred),qcStats[i].counts_by_qual[j]);
				}
			}
			fclose(myfile);
		}


		FILE *myfile;
		myfile = fopen(fastx_outfile,"wd");
		fprintf(myfile,"column\tcount\tmin\tmax\tsum\tmean\tQ1\tmed\tQ3\tIQR\tlW\trW\tA_count\tC_count\tG_count\tT_count\tN_count\tMax_count\n");
		for(int i=0; i<qcStats.size(); i++) {
			int A_tot = 0;
			int C_tot = 0;
			int G_tot = 0;
			int T_tot = 0;
			int N_tot = 0;
			for(int j=0; j<26; j++) {
				if(j==T_A) {
					A_tot += qcStats[i].basecount[j];
				} else if(j==T_C) {
					C_tot += qcStats[i].basecount[j];
				} else if(j==T_G) {
					G_tot += qcStats[i].basecount[j];
				} else if(j==T_T) {
					T_tot += qcStats[i].basecount[j];
				} else {
					N_tot += qcStats[i].basecount[j];
				}
			}

			double q1 = quantiles_with_counts(qcStats[i].counts_by_qual,qcStats[i].qmin,qcStats[i].qmax,.25,0)-phred;
			double med = quantiles_with_counts(qcStats[i].counts_by_qual,qcStats[i].qmin,qcStats[i].qmax,.5,0)-phred;
			double q3 = quantiles_with_counts(qcStats[i].counts_by_qual,qcStats[i].qmin,qcStats[i].qmax,.75,0)-phred;
			
			double iqr = q3-q1;
			int lW = 0;
			int rW = 0;
			
			int low_bound = round(q1-iqr*1.5);
			if(low_bound <= (qcStats[i].qmin-phred)) {
				lW = qcStats[i].qmin-phred;
			} else {
				for(int low=(low_bound+phred);low<=qcStats[i].qmax;low++) {
					if(qcStats[i].counts_by_qual[low] > 0) {
						lW = low-phred;
						low = qcStats[i].qmax+1;
					}
				}
			}
	
			int up_bound = round(q3+iqr*1.5);
			if(up_bound >= (qcStats[i].qmax-phred)) {
				rW = qcStats[i].qmax-phred;
			} else {
				for(int up=(up_bound+phred);up>=qualmin;up--) {
					if(qcStats[i].counts_by_qual[up] > 0) {
						rW = up-phred;
						up = qcStats[i].qmin-1;
					}
				}
			}

			fprintf(myfile,"%d\t%d\t%d\t%d\t%.0f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t", (i+1), qcStats[i].qc, (qcStats[i].qmin-phred),
			                        (qcStats[i].qmax-phred), (qcStats[i].qsum-qcStats[i].qc*phred), 
									(qcStats[i].qsum/qcStats[i].qc-phred),
									 q1, med, q3,iqr, lW, rW);
			fprintf(myfile,"%d\t%d\t%d\t%d\t%d\t%lld\n", A_tot, C_tot, G_tot, T_tot, N_tot,nreads);
		
		}
		fclose(myfile);
	}

	if(brkdown && (!fastx)) {
		FILE *myfile;
		myfile = fopen(brkdown_outfile,"wd");
		fprintf(myfile,"Cycle\tQuality\tCount\n");
		for(int i=0; i<qcStats_by_qual.size(); i++) {
			for(int j=qualmin; j<=qualmax; j++) {
				fprintf(myfile,"%d\t%d\t%d\n",(i+1),(j-phred),qcStats_by_qual[i].counts_by_qual[j]);
			}
		}
		fclose(myfile);
	}

	if(len_hist) {
		FILE *myfile;
		myfile = fopen(lenhist_outfile,"wd");
		fprintf(myfile,"Length\tCount\n");
		for(int len_i=0; len_i<=vlen.size(); len_i++) {
			if(vlen[len_i]) {
				fprintf(myfile,"%d\t%d\n", len_i,vlen[len_i]);
			}
		}
		fclose(myfile);
	}

	int uniq_dup = (int)dup_sort.size();
	if(debug) {
		cout << endl;
		cout << "unique duplicates\t" << uniq_dup << endl;
		cout << "total duplicates\t" << ndups << endl; 
		cout << endl;
	}
	if (uniq_dup && !nodup) {
		printf("dups\t%d\n",ndups-uniq_dup);
		printf("%%dup\t%.4f\n", ((double)(ndups-uniq_dup)/nreads)*100);
	    int uniq_dup = (int)dup_sort.size();
	    printf("unique-dup seq\t%d\n", uniq_dup);
		printf("min dup count\t%d\n", dup_sort.back().cnt);


		for(int i=0; i<show_max; i++) {
			if(i < dup_sort.size()) {
				if(dup_sort.at(i).cnt != 0) {
					cout << "dup seq \t" << (i+1) << "\t" <<  (dup_sort.at(i).cnt-1) << "\t" << dup_sort.at(i).seq << endl;
				}
			} else { i = show_max; }
		}

		if(uniq_dup > 1) {
			printf("dup mean\t%.4f\n", (double)ndups/uniq_dup);
			printf("dup stddev\t%.4f\n", (std_dev((double)uniq_dup, ndups, dupss)));
		}
	}
	printf("qual min\t%d\n", qualmin-phred);
	printf("qual max\t%d\n", qualmax-phred);
	printf("qual mean\t%.4f\n", ((double)qualsum/nbase)-phred);
	printf("qual stdev\t%.4f\n", std_dev((double)nbase,qualsum,qualssq));

	
	if(gc) {
        // put these where they belong
		if (debug)
			printf("gcTotal\t%lu\tgcSum\t%f\n\n", gcTotal, gcSum);
        printf("pct-gc cycle-max\t%d\n", gcCyclemax);
        printf("pct-gc mean\t%.2f\n", 100.0 * gcSum / gcTotal);
    }

	printf("%%A\t%.4f\n", ((double)ACGTN_count[T_A]/nbase*100));
	printf("%%C\t%.4f\n", ((double)ACGTN_count[T_C]/nbase*100));
	printf("%%G\t%.4f\n", ((double)ACGTN_count[T_G]/nbase*100));
	printf("%%T\t%.4f\n", ((double)ACGTN_count[T_T]/nbase*100));
	double ACGT_total  = ACGTN_count[T_A] + ACGTN_count[T_C] + ACGTN_count[T_G] + ACGTN_count[T_T];
	printf("%%N\t%.4f\n", ((double)(nbase-ACGT_total)/nbase*100));
	printf("total bases\t%.0f\n",total_bases);

    if (inputReadError) {   
        printf("error\t%s\n", "error during close, output may be invalid");
    }

    // fail if input read failed....  even if we don't know why and reported all the stats
    return inputReadError;

} //end main method

double quantile( const std::vector <int> & vec, double p ) {
	int l = vec . size();
	double t = ( (double) l- 1 ) * p;
	int it = (int) t;
	int v = vec [it];
	if ( t > (double) it ) {
		return ( v + (t-it) * ( vec [ it + 1 ] - v ) );
	}
	else {
		return v;
	}
} //end quantile function

std::string string_format( const std::string &fmt, ... ) {
	int n, size = 100;
	std::string str;
	va_list ap;
	while (1) {
		str . resize(size);
		va_start( ap, fmt );
		int n =
			vsnprintf( ( char * ) str . c_str(), size, fmt . c_str(), ap );
		va_end(ap);
		if ( n > -1 && n < size ) return str;
		if ( n > -1 ) size = n + 1;
		else size *= 2;
	}
} //end string_format function

double std_dev(double count, double total, double sqsum) {
	if(debug) {
		cout << endl;
		cout << "count " << count << " total " << total << " sqsum " << sqsum << endl;
		cout << endl;
	}
	return sqrt(sqsum/(count-1)-(total/count *total/(count-1)));
}

double quantiles_with_counts(int *v, int start, int end, double p, bool dbug) {
	int v_size = 0;
	for(int i=start; i<=end; i++) {
		if(dbug) 
			cout << "i: " << i << " v[i]: " << v[i] << endl;
		v_size += v[i];
	}
	
	double q = p*(v_size-1);
	int count_skip = (int) q;
	double val = -1;
	bool v_fill = 0;
	int v_next = -1;

	if(dbug) {
		cout << "p : " << p << endl;
		cout << "v-size: " << v_size << endl;
		cout << "q : " << q << endl;
		cout << "count-skip: " << count_skip << endl;
	}
	int tot=0;
	for(int i=start; i<=end; i++) {
		tot += v[i];
		if(tot>count_skip && !v_fill) {
			val = i;
			if(dbug)
				cout << "val : " << val << " val-count: " << v[i] << endl;
			v_fill = 1;
		}
		if(tot>(count_skip+1)) {
			v_next = i;
			if(dbug)
				cout << "val_next : " << v_next << " val-count: " << v[i] << endl;
			i = end+1;			
		}
	}
	
	if(q > count_skip) {
		if(dbug)
			cout << "v_next - val " << (v_next-val) << endl;
		return (val + (q-count_skip)*(v_next-val));
	} else {
		return val;
	}
}

