#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <sys/stat.h>

#include <string>
#include <queue>
#include <google/sparse_hash_map> // or sparse_hash_set, dense_hash_map, ...
#include <google/dense_hash_map> // or sparse_hash_set, dense_hash_map, ...

#include <api/BamReader.h>
#include <api/BamWriter.h>

#include "bamtools_pileup_engine.h"
#include "fastq-lib.h"

const char * VERSION = "1.0";

using namespace BamTools;
using namespace std;

void usage(FILE *f);

#define meminit(l) (memset(&l,0,sizeof(l)))
#define debug(s,...) if (debug) fprintf(stderr,s,##__VA_ARGS__)
#define warn(s,...) ++errs; fprintf(stderr,s,##__VA_ARGS__)
#define stat_out(s,...) fprintf(stdout,s,##__VA_ARGS__)
#define stdev(cnt, sum, ssq) sqrt((((double)cnt)*ssq-pow((double)sum,2)) / ((double)cnt*((double)cnt-1)))

double quantile(const std::vector<int> &vec, double p);
double quantile(const std::vector<double> &vec, double p);
double pnorm(double x);
double qnorm(double x);

std::string string_format(const std::string &fmt, ...);

int debug=0;
int errs=0;
extern int optind;

class xAlign : public BamAlignment {
public:
	int Index;
};

// entry in faidx table
#define FAIBUFN 4096
class Faient {
	FILE *f_ref;
	char buf[FAIBUFN+1];
	long long buf_p;
public:
	long long len;
	long long offset;
	int dataw;
	int linew;
public:
	Faient() {
		len=offset=dataw=linew=buf_p=-1;f_ref=0;
	}
	Faient(const Faient &e) {
		// copy constructor... so vectors work
		len=e.len; offset=e.offset; dataw=e.dataw; linew=e.linew; if (e.buf_p>=0) memcpy(buf, e.buf, FAIBUFN); buf_p=e.buf_p; f_ref=e.f_ref;
	}
	Faient(FILE *f, long long l, long long o, int dw, int lw) {
		// init with info read from file
		f_ref=f; len=l; offset=o; dataw=dw; linew=lw; buf_p=-1;
	}
	char getat(long long pos);
};

class Faidx {
	// hash of faidx entries by name
	google::dense_hash_map<string, Faient> ents;
public:
	Faidx(const string &fa);
	Faient &getent(const string &name) {return ents[name];}
};

class Noise {
public:
	Noise() {noise=0;depth=0;};
	Noise(int d, double n, double q, double mq) {depth=d; noise=n;qnoise=q;mnqual=mq;};
	double noise;
	double qnoise;
	int depth;
	double mnqual;
};

double quantile_depth(const vector<Noise> &vec, double p);

bool noisebydepth (const Noise &a, const Noise &b) { return (a.depth>b.depth);}

class RefVisitor : public PileupVisitor {
    public:
        RefVisitor() : PileupVisitor() {}
        RefVector Refs;
        vector<Faient *> Fais;
};


class VarStatVisitor : public RefVisitor {
    public:
        VarStatVisitor() : RefVisitor() {tot_len=0; num_reads=0;};

    // PileupVisitor interface implementation
    public:
        // prints coverage results ( tab-delimited )
        void Visit(const PileupPosition& pileupData);

    public:
	double tot_len;
	int num_reads;
	vector<Noise> stats;
};

class VarCallVisitor : public RefVisitor {

    public:
        VarCallVisitor() : RefVisitor() {};

    // PileupVisitor interface implementation
    public:
        // prints coverage results ( tab-delimited )
        void Visit(const PileupPosition& pileupData);

    public:
        const RefVector *refs;

};

bool alLessThan (const BamAlignment &a, const BamAlignment &b) {
	if (a.RefID==b.RefID) {
		return a.Position < b.Position;
	} else {
		return a.RefID < b.RefID;
	}
}

bool hasdata(const string &file) {
	struct stat st;
	if (!stat(file.c_str(), &st)) {
		return false;
	}
	return st.st_size > 0;
}

// sort order
bool operator < (const class xAlign &a, const class xAlign &b) {return (!alLessThan(a,b));}

int minsampdepth=10;
double pct_depth=0;
double pct_qdepth=0;
double pct_balance=.10;
int min_depth=0;
int min_adepth=0;
int max_dupreads=5;		// skip duplicate reads greater than this maximum
double zygosity=.5;		// set to .1 for 1 10% admixture, or even .5 

int parse_bams(RefVisitor &v, int in_n, char **in, Faidx *pfai);
FILE *noise_f=NULL;
double alpha=.05;
int phred=33;
double phi(double x);
int main(int argc, char **argv) {
	char c;
	const char *noiseout=NULL;
	const char *ref=NULL;
	optind = 0;
	int umindepth=0;
	int uminadepth=0;
	int upctqdepth=0;
	int do_stats=0;
	int do_varcall=0;
	while ( (c = getopt_long(argc, argv, "?dsvhe:m:N:r:p:a:",NULL,NULL)) != -1) {
		switch (c) {
			case 'd': ++debug; break;
			case 'h': usage(stdout); return 0;
			case 'm': umindepth=atoi(optarg); break;
			case 'a': uminadepth=atoi(optarg);break;
			case 'b': pct_balance=atof(optarg)/100.0; break;
			case 'p': upctqdepth=atoi(optarg); break;
			case 'e': alpha=atof(optarg); break;
			case 'r': ref=optarg; break;
			case 'N': noiseout=optarg; break;
			case 's': do_stats=1; break;
			case 'v': do_varcall=1; break;
			case '?':
					  if (!optopt) {
						  usage(stdout); return 0;
					  } else if (optopt && strchr("ox", optopt))
						  fprintf (stderr, "Option -%c requires an argument.\n", optopt);
					  else if (isprint(optopt))
						  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
					  else
						  fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
					  usage(stderr);
					  return 1;
		}
	}

	if (umindepth > minsampdepth) {
		minsampdepth=umindepth;
	}

	if (!do_stats && !do_varcall) {
		warn("Specify -s for stats only, or -v to do variant calling\n");
		usage(stderr);
		return 1;
	}

	if (!ref) {
		warn("Need a reference sequence (-r) parameter, try -h for help\n");
		return 1;
	}

	if (!hasdata(string(ref)+".fai")) {
		system(string_format("samtools faidx '%s'", ref).c_str());
	}

	Faidx fai(ref);

	if (noiseout) {
		noise_f = fopen(noiseout, "w");
		if (!noise_f) {
			warn("Can't write %s: %s\n", noiseout, strerror(errno));
			exit(1);
		}
	}

	// set argv to '-' if stdin
	const char *stdv[3] = {argv[0],"-",NULL};
	if (!argv[optind]) {
		argc=2;
		argv = (char **) stdv;
		optind=1;
	}

	char **in=&argv[optind];
	int in_n = argc-optind;

	if (do_stats) {
		VarStatVisitor vstat;

		parse_bams(vstat, in_n, in, &fai);

		// sort by depth descending
		sort(vstat.stats.begin(), vstat.stats.end(), noisebydepth);

		// flip 3 and 1 because sorted in descending order for sampling (above)
		double depth_q3=quantile_depth(vstat.stats, .25);
		double depth_q2=quantile_depth(vstat.stats, .50);
		double depth_q1=quantile_depth(vstat.stats, .75);

		// number of locii to compute error rate
		int ncnt=min(10000,vstat.stats.size());

		int i;
		double nsum=0, nssq=0, dsum=0, dmin=vstat.stats[0].depth, qnsum=0, qnssq=0, qualsum=0;
		for (i=0;i<ncnt;++i) {
			nsum+=vstat.stats[i].noise;
			nssq+=vstat.stats[i].noise*vstat.stats[i].noise;
			dsum+=vstat.stats[i].depth;
			qnsum+=vstat.stats[i].qnoise;
			qnssq+=vstat.stats[i].qnoise*vstat.stats[i].qnoise;
			qualsum+=vstat.stats[i].mnqual;
			if (vstat.stats[i].depth < dmin) dmin = vstat.stats[i].depth;
		}

		double noise_mean =nsum/ncnt;
		double noise_dev = stdev(ncnt, nsum, nssq);
		double qnoise_mean =qnsum/ncnt;
		double qnoise_dev = stdev(ncnt, qnsum, qnssq);
		double qual_mean = qualsum/ncnt;

		stat_out("min sampling depth\t%d\n", minsampdepth);
		stat_out("alpha\t%f\n", alpha);

		stat_out("qual mean\t%.4f\n", qual_mean);
		stat_out("noise mean\t%.4f\n", noise_mean);
		stat_out("noise dev\t%.4f\n", noise_dev);
		stat_out("qnoise mean\t%.4f\n", qnoise_mean);
		stat_out("qnoise dev\t%.4f\n", qnoise_dev);

		if (qnoise_mean >= noise_mean ) {
			stat_out("error\tpoor quality estimates\n");
		}

		stat_out("noise depth mean\t%.4f\n", dsum/ncnt);
		stat_out("noise depth min\t%.4f\n", dmin);
		stat_out("noise cnt\t%d\n", ncnt);

		stat_out("depth q1\t%.4f\n", depth_q1);
		stat_out("depth median\t%.4f\n", depth_q2);
		stat_out("depth q3\t%.4f\n", depth_q3);

		dsum=0;
		for (i=0;i<vstat.stats.size();++i) {
			dsum+=vstat.stats[i].depth;
		}

		// pick a min depth that makes sense based on available data
		min_depth=depth_q1 > 10 ? depth_q1 : 10;
		min_depth=max(dmin,min_depth);

		int locii_gtmin=0;
		for (i=0;i<vstat.stats.size();++i) {
			if (vstat.stats[i].depth > min_depth) {
				++locii_gtmin;
			}
		}
		stat_out("locii gt min depth\t%d\n", locii_gtmin);

		double stdevfrommean=-qnorm((alpha/sqrt(locii_gtmin))/2);
		stat_out("qnorm adj\t%f\n", stdevfrommean);

		pct_qdepth=qnoise_mean+qnoise_dev*stdevfrommean;
		stat_out("minority qpct\t%.4f\n", 100*pct_qdepth);

		max_dupreads=depth_q2/(vstat.tot_len/vstat.num_reads);
	//  qpct has better sens/spec, discourage use

	//	pct_depth=noise_mean+noise_dev*stdevfrommean;
	//	stat_out("minority pct\t%.4f\n", 100*pct_depth);

		// min depth required to make a het call, based on the pct_depth and a p-value based on the # of locii tested
		double logz = log(1/zygosity);
		min_adepth = -log(alpha/sqrt(locii_gtmin))/logz;

		stat_out("min depth\t%d\n", min_depth);
		stat_out("min call depth\t%dn", min_adepth);
	}

	if (do_varcall) {
		if (umindepth) min_depth=umindepth;
		if (upctqdepth) pct_qdepth=(double)upctqdepth/100;
		if (uminadepth) min_adepth=uminadepth;

		if (min_depth || (!pct_depth  && !pct_qdepth)) {
			warn("Outputting all variations, no minimum depths specified\n");
		}
		warn("min depth\t%d\n", min_depth);
		warn("min call depth\t%d\n", min_adepth);
		warn("minority qpct\t%d\n", (int)(100*pct_qdepth));
		warn("min balance\t%d\n", (int)(100*pct_balance));

		VarCallVisitor vcall;
		parse_bams(vcall, in_n, in, &fai);
	}
}

double qnorm(double q) {
     if(q == .5)
          return 0;

     q = 1.0 - q;

     double p = (q > 0.0 && q < 0.5) ? q : (1.0 - q);
     double t = sqrt(log(1.0 / pow(p, 2.0)));

     double c0 = 2.515517;
     double c1 = 0.802853;
     double c2 = 0.010328;

     double d1 = 1.432788;
     double d2 = 0.189269;
     double d3 = 0.001308;

     double x = t - (c0 + c1 * t + c2 * pow(t, 2.0)) /
                    (1.0 + d1 * t + d2 * pow(t, 2.0) + d3 * pow(t, 3.0));
    
     if(q > .5)
          x *= -1.0;

     return x;
}

double pnorm(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

int parse_bams(RefVisitor &v, int in_n, char **in, Faidx *pfai) {
	int i;
	RefVector mref;
	BamReader inbam[in_n];
	const char *p;
	for (i=0;i<in_n;++i) {
		if (!inbam[i].Open(in[i]) ) {
			warn("Error reading '%s': %s\n", in[i], strerror(errno));
			return false;
		}
		SamHeader theader = inbam[i].GetHeader();
		RefVector refs = inbam[i].GetReferenceData();
		int i;
		for (i = 0; i < refs.size(); ++i) {
			if (i >= mref.size()) {
				mref.resize(i+1);
				mref[i] = refs[i];
			} else {
				if (refs[i].RefLength != mref[i].RefLength) {
					warn("All files must have the same bam header order, please resort\n");
					exit(1);
				}
			}
		}

	}

	// multiple file pileup	
    priority_queue<xAlign> q;				// ordered queue
	xAlign al;
	for (i=0;i<in_n;++i) {					// read one from each
		inbam[i].GetNextAlignment(al);
		al.Index=i;
		q.push(al);					// sort them in order
	}
	PileupEngine pileup;					

	v.Refs=mref;

	v.Fais.resize(v.Refs.size());

	if (pfai) {
		// check to see if all refereence names are referred to in the FAI specified
		int i;
		for (i=0;i<v.Refs.size();++i) {
			v.Fais[i]=&(pfai->getent(v.Refs[i].RefName));
			if (!(v.Fais[i])->len) {
				warn("No reference data for '%s' in fa", v.Refs[i].RefName.c_str());
				exit(1);
			}
		}
	} else {
		for (i=0;i<v.Refs.size();++i) {
			v.Fais[i]=NULL;
		}
	}

	pileup.AddVisitor(&v);					// visitor watcher, gets called at each position
	int n=0;
	while(!q.empty()) {
		++n;
		const xAlign &tal = q.top();			// next in order
		pileup.AddAlignment(tal); 			// add top align
		if (inbam[tal.Index].GetNextAlignment(al)) {	// read from same file
			al.Index=tal.Index;
			q.pop();				// pop off top
			q.push(al);				// push new
		} else {
			q.pop();				// just pop off, no more in that file
		}
	}
	return n;
}

#define b2i(c) ((c)=='A'?0:(c)=='a'?0:(c)=='C'?1:(c)=='c'?1:(c)=='G'?2:(c)=='g'?2:(c)=='T'?3:(c)=='t'?3:(c)=='-'?4:(c)=='+'?5:6)
#define i2b(i) (i==0?'A':i==1?'C':i==2?'G':i==3?'T':i==4?'-':i==5?'+':'?')

bool hitoloint (int i,int j) { return (i>j);}

class vcall {
public:
	vcall() {base='\0'; fwd=rev=qual=0;};
	char base;
	int qual, fwd, rev;
	int depth() const {return fwd+rev;};
};

bool hitolocall (const vcall &i,const vcall &j) {return ((i.depth())>(j.depth()));}

class PileupSummary {
public:
	PileupSummary(const PileupPosition& pileupData);
	int Depth;
	int TotQual;
	int TotLen;
	int NumReads;
	vector<vcall> Calls;
	int SkipDupReads;
	int MaxDepthByPos;
};

PileupSummary::PileupSummary(const PileupPosition& pileupData) {
	Depth = pileupData.PileupAlignments.size();
	SkipDupReads = 0;
	MaxDepthByPos = 0;
	TotLen = 0;
	NumReads = 0;
	int i;
	vector<int> depthbypos;
	for (i=0;i<Depth;++i) {
		const PileupAlignment &p = pileupData.PileupAlignments[i];
		if (p.PositionInAlignment >= depthbypos.size()) {
			depthbypos.resize(p.PositionInAlignment+1);
		}
		depthbypos[p.PositionInAlignment]++;

		if (depthbypos[p.PositionInAlignment] > max_dupreads) {
			++SkipDupReads;
			continue;
		}

		if (p.PositionInQuery == 0) {
			++NumReads;
			TotLen+=p.Alignment.Qualities.size();
		}

		if (p.PositionInQuery >= p.Alignment.Qualities.size()) {
			if (p.PositionInQuery < 0) {
				// ignore deletions at the front of a query, makes no sense
				continue;
			}
			warn("Bad query position %d in %s\n", p.PositionInQuery, p.Alignment.Name.c_str());
			continue;
		}

		if (p.PositionInAlignment >= p.Alignment.AlignedBases.size()) {
			// why does this happen?
			continue;
		}

		char q = p.Alignment.Qualities[p.PositionInQuery];
		char c = p.Alignment.AlignedBases[p.PositionInAlignment];

		int j = b2i(c);
		if (j >= Calls.size()) {
			int was = Calls.size();
			Calls.resize(j+1);
			int t; for (t=was;t<=j;++t) {
				Calls[t].base=i2b(t);;
			}
		}

		if (p.Alignment.IsReverseStrand()) {
			++Calls[j].rev;
		} else {
			++Calls[j].fwd;
		}

		Calls[j].qual+=q-phred;

		if (p.IsCurrentDeletion) {
		} else {
			if (p.IsNextInsertion) {
				c = '+';
				j = b2i(c);
				if (j >= Calls.size()) {
					int was = Calls.size();
					Calls.resize(j+1);
					int t; for (t=was;t<=j;++t) {
						Calls[t].base=i2b(t);;
					}
				}
				if (p.Alignment.IsReverseStrand()) {
					++Calls[j].rev;
				} else {
					++Calls[j].fwd;
				}
			} else {
			}
		}
	}

	for (i=0;i<depthbypos.size();++i) {
		if (depthbypos[i] > MaxDepthByPos) {
			MaxDepthByPos = depthbypos[i];
		}
	}

	Depth=0;
	for (i=0;i<5 && i < Calls.size();++i) {		// total depth (exclude inserts for tot depth, otherwise they are double-counted)
		Depth+=Calls[i].depth();
	}

	TotQual=0;
	for (i=0;i<5 && i < Calls.size();++i) {		// total depth (exclude inserts for tot depth, otherwise they are double-counted)
		TotQual+=Calls[i].qual;
	}
}

void VarCallVisitor::Visit(const PileupPosition& pileupData) {
	int depth = pileupData.PileupAlignments.size();
	if (depth < min_depth) 
		return;

	PileupSummary p(pileupData);	
	if (p.Depth < min_depth)
		return;

	int i;
	if (p.Calls.size() > 6) 
		p.Calls.resize(6);	// toss N's before sort
	sort(p.Calls.begin(), p.Calls.end(), hitolocall);

	string pil;
	for (i=0;i<p.Calls.size();++i) {		// all calls
		double pct = (double) p.Calls[i].depth()/p.Depth;
		double qpct = (double) p.Calls[i].qual/p.TotQual;
		if (!p.Calls[i].base)
			continue;
		if (pct > pct_depth && qpct > pct_qdepth && (p.Calls[i].depth() >= min_adepth)) {
			double bpct = (double) min(p.Calls[i].fwd,p.Calls[i].rev)/p.Calls[i].depth();
			if (bpct > pct_balance) {
				pil += string_format("\t%c:%d,%d", p.Calls[i].base,p.Calls[i].fwd,p.Calls[i].rev);	
			} else {
//				++skipped_balance;
			}
		} else
			break;
	}
	if (pil.size() > 0) {
		string ref;
		char base = 'N';
		if (pileupData.RefId < Refs.size()) { 
			ref=Refs[pileupData.RefId].RefName;
			base=Fais[pileupData.RefId]->getat(pileupData.Position+1);
		}
		printf("%s\t%d\t%c\t%d\t%d%s\n",ref.c_str(), pileupData.Position, base, p.Depth, p.SkipDupReads,pil.c_str());
	}
}


void VarStatVisitor::Visit(const PileupPosition& pileupData) {
	int depth = pileupData.PileupAlignments.size();
	if (depth < minsampdepth)
		return;

	PileupSummary p(pileupData);	
	if (p.Depth < minsampdepth)
		return;

	int i;
	if (p.Calls.size() > 5) 
		p.Calls.resize(5);		// toss N's and '+' before sort

	sort(p.Calls.begin(), p.Calls.end(), hitolocall);

	double noise = p.Calls.size() > 2 ? (double) p.Calls[2].depth()/p.Depth : 0;
	double qnoise = p.Calls.size() > 2 ? (double) p.Calls[2].qual/p.TotQual : 0;
	double mnqual = (double)p.TotQual/p.Depth;
	char pbase = p.Calls.size() > 2 ? p.Calls[2].base : '.';

	if (noise_f) {
		fprintf(noise_f,"%d\t%c\t%f\t%f\n", p.Depth, pbase, noise, qnoise, mnqual);
	}

	tot_len += p.TotLen;
	num_reads += p.NumReads;
	stats.push_back(Noise(p.Depth, noise, qnoise, mnqual));
}


void usage(FILE *f) {
        fprintf(f,
"Usage: varcall <-s|-v> [options] file1 [file2...]\n"
"Version: %s\n"
"\n"
"Either outputs stats for the list of files, or performs variant calling\n"
"\n"
"Options:\n"
"\n"
"-s             Calculate statistics\n"
"-v             Calculate variants bases on supplied parameters\n"
"-S FIL         Read statistics from previous run\n"
"-A ANNOT       Calculate in-target stats using the annotation file\n"
"\n"
"Stats Output:\n"
"\n"
"Contains mean, median, quartile information for depth, base quality, read len,\n"
"mapping quality, indel levels. Also estimates parameters suitable for\n"
"variant calls, and can be passed directly to this program for variant calls\n"
"\n"
"Variant Output:\n"
"\n"
        ,VERSION);
}

std::string string_format(const std::string &fmt, ...) {
       int n, size=100;
       std::string str;
       va_list ap;
       while (1) {
       str.resize(size);
       va_start(ap, fmt);
       int n = vsnprintf((char *)str.c_str(), size, fmt.c_str(), ap);
       va_end(ap);
       if (n > -1 && n < size)
           return str;
       if (n > -1)
           size=n+1;
       else
           size*=2;
       }
}

// returns quantile depth 
double quantile_depth(const std::vector<Noise> &vec, double p) {
        int l = vec.size();
        double t = ((double)l-1)*p;
        int it = (int) t;
        int v=vec[it].depth;
        if (t > (double)it) {
                return (v + (t-it) * (vec[it+1].depth - v));
        } else {
                return v;
        }
}

double quantile(const std::vector<int> &vec, double p) {
        int l = vec.size();
        double t = ((double)l-1)*p;
        int it = (int) t;
        int v=vec[it];
        if (t > (double)it) {
                return (v + (t-it) * (vec[it+1] - v));
        } else {
                return v;
        }
}

double quantile(const std::vector<double> &vec, double p) {
        int l = vec.size();
        double t = ((double)l-1)*p;
        int it = (int) t;
        double v=vec[it];
        if (t > (double)it) {
                return (v + p * (vec[it+1] - v));
        } else {
                return v;
        }
}

Faidx::Faidx(const string &fa) {
	ents.set_empty_key("n/a");

	FILE *ref=fopen(fa.c_str(),"r");
	if (!ref) {
			warn("Can't open %s: %s\n", fa.c_str(), strerror(errno));
			exit(1);
	}
	FILE *f = fopen((fa+".fai").c_str(),"r");
	if (!f) {
			warn("Can't open %s.fai: %s\n", fa.c_str(), strerror(errno));
			exit(1);
	}
	line l; meminit(l);
	while(read_line(f, l)>0) {
		char *cols[7];
		char *sp;
		char *t=strtok_r(l.s, " \t", &sp);
		int i=0;
		while(t) {
			cols[i++]=t;
			t=strtok_r(NULL, " \t", &sp);
			if (i > 5) break;
		};
		if (i < 5 || (atoll(cols[4]) == 0)) {
			warn("Invalid FAI file '%s', requires 5 columns, all nonzero\n", (fa+".fai").c_str());
			continue;
		}
		// add reference
		ents[cols[0]]=Faient(ref, atoll(cols[1]), atoll(cols[2]), atoi(cols[3]), atoi(cols[4]));
	}
}

char Faient::getat(long long pos) {
	// taken from faidx.c (samtools)
	long long fpos=offset + ((pos / dataw) * linew) + (pos % dataw);

	// readahead buffer (possibly useless?  benchmark)
	if (buf_p == -1 || pos < buf_p || pos >= (buf_p+FAIBUFN) ) {
		 if (fseek(f_ref, fpos, SEEK_SET)) {
				return '\0';
		 }
		 int n = fread(buf, 1, FAIBUFN, f_ref);
		 if (n == 0) {
				return '\0';
		 }
		 if (n < FAIBUFN) {
			memset(buf+n, 0, FAIBUFN-n);
		 }
		 buf_p=fpos;
	}
	return buf[fpos-buf_p];
}


