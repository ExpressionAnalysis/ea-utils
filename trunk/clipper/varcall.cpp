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
#include <list>

#include <google/sparse_hash_map> // or sparse_hash_set, dense_hash_map, ...
#include <google/dense_hash_map> // or sparse_hash_set, dense_hash_map, ...

#include "fastq-lib.h"

const char * VERSION = "0.8.6";

#define MIN_READ_LEN 20

using namespace std;
using namespace google;

void usage(FILE *f);

#define meminit(l) (memset(&l,0,sizeof(l)))
#define debug(s,...) if (debug) fprintf(stderr,s,##__VA_ARGS__)
#define warn(s,...) ++errs; fprintf(stderr,s,##__VA_ARGS__)
#define die(s,...) (fprintf(stderr,s,##__VA_ARGS__), exit(1))
#define stat_out(s,...) fprintf(stat_fout,s,##__VA_ARGS__)
#define stdev(cnt, sum, ssq) sqrt((((double)cnt)*ssq-pow((double)sum,2)) / ((double)cnt*((double)cnt-1)))

double quantile(const std::vector<int> &vec, double p);
double quantile(const std::vector<double> &vec, double p);
double pnorm(double x);
double qnorm(double x);

// basic utils
std::vector<char *> split(char* str, const char* delim);
std::string string_format(const std::string &fmt, ...);
void to_upper(const std::string str);

int debug=0;
int errs=0;
extern int optind;

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

class PileupEnt {
public:
	bool is_rev;
	bool is_start;
	bool f;
	const char *b;
	int q;
	int m;
	int p;
};

class vcall {
public:
    vcall() {base='\0'; fwd=rev=qual=is_ref=0;};
    char base;
	bool is_ref;
    int qual, fwd, rev;
	vector <string> seqs;
    int depth() const {return fwd+rev;};
};

bool hitolocall (const vcall &i,const vcall &j) {return ((i.depth())>(j.depth()));}

class Read {
public:
    int MapQ;
    string Seq;
    Read() {MapQ=0;};
};

class PileupReads {
public:
    double MeanReadLen() {return ReadBin.size() ? TotReadLen/ReadBin.size() : MIN_READ_LEN;}
    int TotReadLen;
    deque<Read> ReadBin;
    list<Read> ReadList;
    PileupReads() {TotReadLen=0;}
};

class PileupSummary {
public:
    string Chr;
    int Pos;
    char Base;
    int Depth;
    int TotQual;
    int NumReads;
    vector<vcall> Calls;

    int SkipDupReads;
    int SkipMinMapq;
    int SkipMinQual;
    int MaxDepthByPos;
    int RepeatCount;
    char RepeatBase;

	PileupSummary(char *line, PileupReads &reads);
    PileupSummary() { Base = '\0'; Pos=-1; };
};

class PileupVisitor {
    public:
        char InputType;
        PileupReads Reads;
        PileupVisitor() {InputType ='\0';}
		void Parse(char *dat) {PileupSummary p(dat, Reads); Visit(p);};
		virtual void Visit(PileupSummary &dat)=0;
		virtual void Finish()=0;
};

class VarStatVisitor : public PileupVisitor {
    public:
    VarStatVisitor() : PileupVisitor() {tot_depth=0; num_reads=0;};

    void Visit(PileupSummary &dat);
    void Finish() {};

    public:
	double tot_depth;
	int num_reads;
	vector<Noise> stats;
};

class VarCallVisitor : public PileupVisitor {

    deque<PileupSummary> Win;
    void VisitX(PileupSummary &dat);

    public:
    int WinMax;
    VarCallVisitor() : PileupVisitor() {SkippedDepth=0;WinMax=0;};

    void Visit(PileupSummary &dat);
    void Finish();

	int SkippedDepth;
};

bool hasdata(const string &file) {
	struct stat st;
	if (stat(file.c_str(), &st)) {
		return false;
	}
	return st.st_size > 0;
}


int minsampdepth=10;
double pct_depth=0;
double pct_qdepth=0;
double pct_balance=0;
char *debug_xchr=NULL;
int debug_xpos=0;
int min_depth=0;
int min_mapq=1;
int min_qual=3;
int repeat_filter=8;
double artifact_filter=1;
int min_adepth=0;
int min_idepth=0;
int no_baq=0;
double zygosity=.5;		        // set to .1 for 1 10% admixture, or even .05 for het/admix

void parse_bams(PileupVisitor &v, int in_n, char **in, const char *ref);

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
	int uminidepth=0;
	double upctqdepth=0;
	int do_stats=0;
	int do_varcall=0;
	while ( (c = getopt_long(argc, argv, "?dsvBhe:m:N:x:f:p:a:q:Q:i:D:R:b:",NULL,NULL)) != -1) {
		switch (c) {
			case 'd': ++debug; break;
			case 'h': usage(stdout); return 0;
			case 'm': umindepth=atoi(optarg); break;
			case 'q': min_qual=atoi(optarg); break;
			case 'Q': min_mapq=atoi(optarg); break;
			case 'R': repeat_filter=atoi(optarg); break;
			case 'a': uminadepth=atoi(optarg);break;
			case 'D': artifact_filter=atof(optarg);break;
			case 'i': uminidepth=atoi(optarg);break;
			case 'x': {
					debug_xchr=optarg;
					char *p=strchr(debug_xchr, ':');
					if (!p) die("Invalid param for -x");
					*p='\0';
					debug_xpos=atoi(++p);
					if (!p) die("Invalid param for -x, need pos");
					break;
				}
			case 'b': pct_balance=atof(optarg)/100.0; break;
			case 'B': no_baq=1; break;
			case 'p': upctqdepth=atof(optarg); break;
			case 'e': alpha=atof(optarg); break;
			case 'f': ref=optarg; break;
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

	if (!do_stats && !do_varcall || do_stats && do_varcall) {
		warn("Specify -s for stats only, or -v to do variant calling\n\n");
		usage(stderr);
		return 1;
	}

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
		FILE *stat_fout=stdout;			// stats to stdout

		if (do_varcall) 				// unless varcalling at the same time
			stat_fout=stderr;

		VarStatVisitor vstat;

		parse_bams(vstat, in_n, in, ref);

		// sort by depth descending
		sort(vstat.stats.begin(), vstat.stats.end(), noisebydepth);

		// flip 3 and 1 because sorted in descending order for sampling (above)
		double depth_q3=quantile_depth(vstat.stats, .25);
		double depth_q2=quantile_depth(vstat.stats, .50);
		double depth_q1=quantile_depth(vstat.stats, .75);
		double depth_qx=quantile_depth(vstat.stats, .90);

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
		min_depth=depth_qx > 10 ? (depth_qx < 100 ? depth_qx : 100) : 10;
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
		stat_out("min pct qual\t%.4f\n", 100*pct_qdepth);

	//  qpct has better sens/spec, discourage use

	//	pct_depth=noise_mean+noise_dev*stdevfrommean;
	//	stat_out("minority pct\t%.4f\n", 100*pct_depth);

		// min depth required to make a het call, based on the pct_depth and a p-value based on the # of locii tested
		double logz = log(1/zygosity);
		min_adepth = -log(alpha/sqrt(locii_gtmin))/logz;

		stat_out("min depth\t%d\n", min_depth);
		stat_out("min call depth\t%d\n", min_adepth);
	}

	if (do_varcall) {
		if (umindepth) min_depth=umindepth;
		if (upctqdepth > 0) pct_qdepth=(double)upctqdepth/100;
		if (uminadepth) min_adepth=uminadepth;
		if (uminidepth) min_idepth=uminidepth;

		if (min_depth || (!pct_depth  && !pct_qdepth)) {
			warn("warning\toutputting all variations, no minimum depths specified\n");
		}

		warn("version\tvarcall-%s\n", VERSION);
		warn("min depth\t%d\n", min_depth);
		warn("min call depth\t%d\n", min_adepth);
		warn("min pct qual\t%d\n", (int)(100*pct_qdepth));

		warn("min balance\t%d\n", (int)(100*pct_balance));
		warn("min qual\t%d\n", min_qual);
		warn("min map qual\t%d\n", min_mapq);

		VarCallVisitor vcall;

        if (repeat_filter > 0) {
		    warn("homopolymer filter\t%d\n", repeat_filter);
            vcall.WinMax=repeat_filter+repeat_filter+1;
        }

		parse_bams(vcall, in_n, in, ref);

        if (vcall.InputType == 'B') {
        	warn("baq correct\t%s\n", (no_baq?"no":"yes"));
        }
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

void parse_bams(PileupVisitor &v, int in_n, char **in, const char *ref) {

	if (!in_n) {
		warn("No input files, quitting\n");
		exit(1);
	}

	int i, bam_n=0;
	for (i=0;i<in_n;++i) {
		if (!strcmp(fext(in[i]), ".bam")) {
			++bam_n;
		}
	}

	if (bam_n != in_n) {
		if (bam_n > 0) {
			warn("Can't mix bams and other input files\n");
			exit(1);
		} else {
			if (in_n > 1) {
				warn("Can't handle multiple pileups... TODO\n");
				exit(1);
			} else {
				warn("input\t%d pileup\n", in_n);
                v.InputType='P';
			}
		}
	} else {
		warn("input\t%d bam\n", bam_n);
        v.InputType='B';
	}

	int is_popen = 0;
	FILE *fin;

	if (bam_n) {
        if (!ref) {
            warn("Need a reference file (-f) parameter, try -h for help\n");
            exit(1);
        }

        if (!hasdata(string(ref)+".fai")) {
            int ret=system(string_format("samtools faidx '%s'", ref).c_str());
            if (ret) {
                warn("Need a %s.fai file, run samtools faidx\n", ref);
                exit(1);
            }
        }


		const char *nobaq = no_baq ? "-B" : "";

		string mpil_cmd = string_format("samtools mpileup -d 100000 %s -f '%s'", nobaq, ref);

		int i;
		for (i=0;i<in_n;++i) {
			mpil_cmd += " '";
			mpil_cmd += in[i];
			mpil_cmd += "' ";
		}

		warn("command\t%s\n", mpil_cmd.c_str());

		fin = popen(mpil_cmd.c_str(), "r");
		if (!fin) 
			exit(1);

		is_popen = 1;
	} else {
		fin = fopen(in[0], "r");
		if (!fin) {
			warn("%s: %s", in[0], strerror(errno));
		}
	}

    line l; meminit(l);
	int cnt=0;
    while(read_line(fin, l)>0) {
	//	chr      2       G       6       ^9,^+.^*,^2,^&.^&,      &.'&*-  9+*2&&  166,552,643,201,299,321
		v.Parse(l.s);
		++cnt;
	}

	if (is_popen) pclose(fin); else fclose(fin);

	if (cnt == 0) {
		warn("No data in pileup, quitting\n");
		exit(1);
	}
}

#define b2i(c) ((c)=='A'?0:(c)=='a'?0:(c)=='C'?1:(c)=='c'?1:(c)=='G'?2:(c)=='g'?2:(c)=='T'?3:(c)=='t'?3:(c)=='*'?4:(c)=='-'?5:(c)=='+'?6:7)
#define i2b(i) (i==0?'A':i==1?'C':i==2?'G':i==3?'T':i==4?'*':i==5?'-':i==6?'+':'?')

bool hitoloint (int i,int j) { return (i>j);}

int track_readlen[10000];


PileupSummary::PileupSummary(char *line, PileupReads &rds) {

	vector<char *> d=split(line, "\t");

	if (d.size() < 6) {
		warn("Can't read pileup : %d fields, need 6 columns\n", (int) d.size());
		exit(1);
	}

	const char * p_qual=d[5];

	Chr=d[0];
	Pos=atoi(d[1]);
	Base=*(d[2]);
	Depth = atoi(d[3]);
	SkipDupReads = 0;
	SkipMinQual = 0;
	SkipMinMapq = 0;
	MaxDepthByPos = 0;
	RepeatCount = 0;
	RepeatBase = '\0';
	NumReads = 0;

	int i;
	vector<int> depthbypos;

	const char *cur_p = d[4];

    list<Read>::iterator read_i = rds.ReadList.begin();

    int eor=0;
	for (i=0;i<Depth;++i,++read_i) {
		bool sor=0;
		
		if (*cur_p == '^') {
			sor=1;
			++cur_p;
            Read x;
            x.MapQ = *cur_p-phred;
            ++cur_p;
            if (read_i != rds.ReadList.end()) {
                ++read_i;
            }
            read_i=rds.ReadList.insert(read_i,x);
		}

        if (read_i == rds.ReadList.end()) {
            warn("warning\tread start without '^', partial pileup: '%s'\n", cur_p);
            Read x;
            x.MapQ = -1;
            read_i=rds.ReadList.insert(read_i,x);
        }

        int pia = read_i->Seq.length()+1;
		if (pia >= depthbypos.size()) {
			depthbypos.resize(pia+1);
		}
		depthbypos[pia]++;


		if (sor) 
			++NumReads;

		char q = p_qual[i]-phred;				// qual char
		char mq = read_i->MapQ;
		char o = *cur_p;				// orig call
		char c = toupper(o);			// uppercase/ref 
		bool is_ref = 0;

		if (o == '.' || o == ',') {	
			c = Base;					// ref instead
			is_ref = 1;
		}

		bool skip = 0;

		if (artifact_filter > 0 && (depthbypos[pia] > artifact_filter * (1+(Depth/rds.MeanReadLen())))) {
			++SkipDupReads;
			skip=1;
		} else if (mq < min_mapq) {
			++SkipMinMapq;
			skip=1;
		} else if (q < min_qual) {
			++SkipMinQual;
			skip=1;
		} else {
			int j = b2i(c);
			if (j >= Calls.size()) {
				int was = Calls.size();
				Calls.resize(j+1);
				int t; for (t=was;t<=j;++t) {
					Calls[t].base=i2b(t);
				}
			}
			if (is_ref) 
				Calls[j].is_ref = 1;

			if ( o == ',' || o == 'a' || o == 'c' || o == 't' || o == 'g' ) {
				++Calls[j].rev;
			} else {
				++Calls[j].fwd;
			}

			Calls[j].qual+=q;
		}

		if (c == '-' || c == '+') {
            warn("invalid pileup, at '%s', indel not attached to read?\n", cur_p);
		} else {
		    if (c != '*') 
                read_i->Seq += c;
			++cur_p;
		}

        if (*cur_p == '+' || *cur_p == '-') {
            c = *cur_p;
            char *end_p;
            int len = strtol(++cur_p, &end_p, 10);
            string ins_seq(end_p, len);
            to_upper(ins_seq);
            read_i->Seq += ins_seq;
            if (!skip) {
                int j = b2i(c);
                if (j >= Calls.size()) {
                    int was = Calls.size();
                    Calls.resize(j+1);
                    int t; for (t=was;t<=j;++t) {
                        Calls[t].base=i2b(t);
                    }
                }
                if ( o == ',' || o == 'a' || o == 'c' || o == 't' || o == 'g' ) {
                    ++Calls[j].rev;
                } else {
                    ++Calls[j].fwd;
                }
                Calls[j].qual+=q;
                Calls[j].seqs.push_back(ins_seq);
            }
            cur_p=end_p+len;
        }

        if (*cur_p == '$') {
            if (read_i->MapQ > -1) {
                rds.TotReadLen+=read_i->Seq.size();
                rds.ReadBin.push_back(*read_i);
                if (rds.ReadBin.size() > min(1000,Depth*2)) {
                    rds.ReadBin.pop_front();
                    rds.TotReadLen-=rds.ReadBin.front().Seq.size();
                }
            }
//            printf("%d\t%s\n", read_i->MapQ, read_i->Seq.c_str());
            read_i=rds.ReadList.erase(read_i);
            --read_i;
            ++cur_p;
            ++eor;
        }
	}

    if ((Depth-eor) != rds.ReadList.size()) {
        warn("warning\tdepth is %d, but read list is: %d\n", Depth, (int) rds.ReadList.size());
    }

	if (*cur_p == '-' || *cur_p == '+') {
		char *end_p;
		int len = strtol(++cur_p, &end_p, 10);
		// keep this
		string idl(end_p, len);
		cur_p=end_p+len;
	}

	if (*cur_p) {
		warn("Failed to parse pileup %s\n", d[4]);
		exit(1);
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

PileupSummary JunkSummary;

void VarCallVisitor::Visit(PileupSummary &p) {
    if (WinMax < 5) {
        VisitX(p);
        return;
    }

    if (Win.size() && (Win.back().Pos != (p.Pos - 1) )) {
        if (p.Base != '-') {
            if (Win.back().Pos < p.Pos && ((p.Pos - Win.back().Pos) <= (WinMax/2))) {
                while (Win.back().Pos < (p.Pos - 1)) {
                    // visit/pop, add a placeholder
                    JunkSummary.Base = '-';
                    JunkSummary.Pos = Win.back().Pos + 1;
                    Visit(JunkSummary);
                }
            } else {
                while (Win.size()) {
                    // visit/pop, but don't add anything, until it's empty
                    JunkSummary.Base = '@';
                    JunkSummary.Pos = 0;
                    Visit(JunkSummary);
                }
            }
        } 
    }

    if (p.Base != '@')              // see above, false-add to clear queue
        Win.push_back(p);

    if (Win.size() > WinMax)        // queue too big?  pop
        Win.pop_front();
    
    int i;
    int lrc=0,rrc=0;                // left repeat count, right repeat count
    char lrb, rrb;                  // left repeat base...
    int vx;

    if (Win.size() < WinMax/2) {    // small window?  look at leading edge only
        vx = Win.size()-1;
    } else {
        vx = WinMax/2;              // larger window? look at midpoint
    }

    if (Win[vx].Base == '-') 
        return;

    if (vx > 1) {                   // look left
        lrb = Win[vx-1].Base;
        for (i=vx-2; i >= 0; --i) { // increment repeat count
            if (Win[i].Base == lrb) 
                ++lrc;
            else 
                break;
        }
    }
    if (vx < (Win.size()-2)) {
        rrb = Win[vx+1].Base;
        for (i=vx+1; i < Win.size(); ++i) {
            if (Win[i].Base == rrb)
                ++rrc;
            else
                break;
        }
    }

    // maximum repeat count and associated base
    if (lrb == rrb ) {
        Win[vx].RepeatCount = lrc+rrc;
        Win[vx].RepeatBase = lrb;
    } else if (lrb > rrb) {
        Win[vx].RepeatCount = lrc;
        Win[vx].RepeatBase = lrb;
    } else {
        Win[vx].RepeatCount = rrc;
        Win[vx].RepeatBase = rrb;
    }

	if (debug_xpos) {
        if (Win[vx].Pos == debug_xpos && !strcmp(debug_xchr,Win[vx].Chr.data())) {
            fprintf(stderr,"xpos-window\t");
            for (i=0;i<Win.size();++i) {
                fprintf(stderr,"%c", Win[i].Base);
            }
            fprintf(stderr,"\n");
        }
    }

    VisitX(Win[vx]);
}

void VarCallVisitor::Finish() {
    int vx = WinMax/2;
    while (vx < Win.size()) {
        VisitX(Win[vx++]);
    }
}

void VarCallVisitor::VisitX(PileupSummary &p) {
	if (p.Depth < min_depth) {
		++SkippedDepth;
		return;
	}

	if (debug_xpos) {
		if (p.Pos != debug_xpos)
			return;
		if (strcmp(debug_xchr,p.Chr.data())) 
			return;
	}

	int i;
	if (p.Calls.size() > 6) 
		p.Calls.resize(7);	// toss N's before sort
	sort(p.Calls.begin(), p.Calls.end(), hitolocall);

	int need_out = -1;
	int skipped_balance=0;
	int skipped_indel=0;
	int skipped_depth=0;
	int skipped_repeat=0;
	string pil;
	for (i=0;i<p.Calls.size();++i) {		// all calls
//		printf("d:%d b: %c, pd: %d\n", (int) p.Calls[i].depth(), p.Calls[i].base, p.Depth);
	
		double pct = (double) p.Calls[i].depth()/p.Depth;
		double qpct = (double) p.Calls[i].qual/p.TotQual;

		if (!p.Calls[i].base)
			continue;

		if (pct > pct_depth && qpct >= pct_qdepth && (p.Calls[i].depth() >= min_adepth)) {
			double bpct = (double) min(p.Calls[i].fwd,p.Calls[i].rev)/p.Calls[i].depth();
			// balance is meaningless at low depths
			if ((bpct >= pct_balance) || (p.Calls[i].depth()<4)) {
				if (p.Calls[i].base == '+' || p.Calls[i].base == '-') {
                    // yuk ... time to think about a possible indel call
					if (p.Calls[i].depth() >= min_idepth) {
						// should really pick more than 1
						// but need to allow "similar" indels to pile up
                        // should group into distinct bins, using some homology thing
						sort(p.Calls[i].seqs.begin(), p.Calls[i].seqs.end());
						string prev, maxs;
						int pcnt=0, maxc=0, j;
						for (j=0;j<p.Calls[i].seqs.size();++j) {
							if (prev == p.Calls[i].seqs[j]) {
								++pcnt;
							} else {
								if (pcnt > maxc) {
									maxs=prev;
									maxc=pcnt;
								}
								prev=p.Calls[i].seqs[j];
								pcnt=1;
							}
						}
						if (pcnt > maxc) {
							maxs=prev;
							maxc=pcnt;
						}
						if (maxc >= min_idepth && maxc >= min_adepth) {
                            // only calls 1 indel at a given position
                            if (p.RepeatCount < repeat_filter) {
                                if (need_out == -1) 
                                    need_out = i;
                                pil += string_format("\t%c%s:%d,%d", p.Calls[i].base, maxs.c_str(),maxc,p.Calls[i].qual/p.Calls[i].depth());	
                                skipped_indel+=p.Calls[i].depth()-maxc;
                            } else {
                                skipped_repeat+=p.Calls[i].depth();
                            }
						} else {
							skipped_indel+=p.Calls[i].depth();
						}
					} else {
						skipped_indel+=p.Calls[i].depth();
					}
				} else {
					if (!p.Calls[i].is_ref || debug_xpos) { 
                        if (need_out == -1)
                            need_out = i;
                    }
					pil += string_format("\t%c:%d,%d", p.Calls[i].base,p.Calls[i].depth(),p.Calls[i].qual/p.Calls[i].depth());	
				}
			} else {
				skipped_balance+=p.Calls[i].depth();
			}
		} else {
			skipped_depth+=p.Calls[i].depth();
			break;
		}
	}
	if (need_out>=0||debug_xpos) {
        double pct_allele = 100.0 * p.Calls[need_out].depth() / (double) p.Depth;
		printf("%s\t%d\t%c\t%d\t%d\t%2.2f%s\n",p.Chr.c_str(), p.Pos, p.Base, p.Depth, skipped_depth+skipped_balance+p.SkipDupReads+p.SkipMinMapq+p.SkipMinQual, pct_allele, pil.c_str());
		if (debug_xpos) {
		    fprintf(stderr,"xpos-skip-dup\t%d\n",p.SkipDupReads);
		    fprintf(stderr,"xpos-skip-mapq\t%d\n",p.SkipMinMapq);
		    fprintf(stderr,"xpos-skip-qual\t%d\n",p.SkipMinQual);
		    fprintf(stderr,"xpos-skip-bal\t%d\n",skipped_balance);
		    fprintf(stderr,"xpos-skip-depth\t%d\n",skipped_depth);
		    fprintf(stderr,"xpos-skip-indel\t%d\n",skipped_indel);
		    fprintf(stderr,"xpos-skip-repeat\t%d\n",skipped_repeat);
            if (repeat_filter > 0) {
                fprintf(stderr,"repeat-count\t%d\n",p.RepeatCount);
                fprintf(stderr,"repeat-base\t%c\n",p.RepeatBase);
            }
			exit(0);
		}
	}
}


void VarStatVisitor::Visit(PileupSummary &p) {
	if (p.Depth < minsampdepth)
		return;

    // insert and deletions have their own, separate noise levels

	int ins_depth = p.Calls.size() > 6 ? p.Calls[6].depth() : 0;
	int ins_qual = p.Calls.size() > 6 ? p.Calls[6].qual : 0;
	double ins_noise = 0;
	double ins_qnoise = 0;
	if (p.Calls.size() > 1 && p.Calls[1].depth() > ins_depth && ins_depth > 0) {
		ins_noise = (double) ins_depth/p.Depth;
		ins_qnoise = (double) ins_qual/p.TotQual;
	}

	int del_depth = p.Calls.size() > 5 ? p.Calls[5].depth() : 0;
	int del_qual = p.Calls.size() > 5 ? p.Calls[5].qual : 0;
	double del_noise = 0;
	double del_qnoise = 0;
	if (p.Calls.size() > 1 && p.Calls[1].depth() > del_depth && del_depth > 0) {
		del_noise = (double) del_depth/p.Depth;
		del_qnoise = (double) del_qual/p.TotQual;
	}

    // snp's are "noise" if there are 3 alleles at a given position
	int i;
	if (p.Calls.size() > 5) 
		p.Calls.resize(5);		// toss N's and inserts before sort

	sort(p.Calls.begin(), p.Calls.end(), hitolocall);

	double noise = p.Calls.size() > 2 ? (double) p.Calls[2].depth()/p.Depth : 0;
	double qnoise = p.Calls.size() > 2 ? (double) p.Calls[2].qual/p.TotQual : 0;

	double mnqual = (double)p.TotQual/p.Depth;

	char pbase = p.Calls.size() > 2 ? p.Calls[2].base : '.';

	if (noise_f) {
		fprintf(noise_f,"%d\t%c\t%f\t%f\n", p.Depth, pbase, noise, qnoise, mnqual);
/*
        if (ins_noise > 0) {
		    fprintf(noise_f,"%d\t%c\t%f\t%f\n", p.Depth, '+', ins_noise, ins_qnoise, mnqual);
        }
        if (del_noise > 0) {
		    fprintf(noise_f,"%d\t%c\t%f\t%f\n", p.Depth, '-', del_noise, del_qnoise, mnqual);
        }
*/
	}

	tot_depth += p.Depth;
	num_reads += p.NumReads;
	stats.push_back(Noise(p.Depth, noise, qnoise, mnqual));
}


void usage(FILE *f) {
        fprintf(f,
"Usage: varcall <-s|-v> <-f REF> [options] bam1 [bam2...]\n"
"Version: %s (BETA)\n"
"\n"
"Either outputs summry stats for the list of files, or performs variant calling\n"
"\n"
"Options:\n"
"\n"
"-s          Calculate statistics\n"
"-v          Calculate variants bases on supplied parameters\n"
"-f          Reference fasta (required if using bams)\n"
"-m          Min locii depth (0)\n"
"-a          Min allele depth (0)\n"
"-p          Min allele pct by quality (0)\n"
"-q          Min qual (3)\n"
"-C          Min mapping quality (1)\n"
"-b          Min pct balance (strand/total) (0)\n"
"-D FLOAT    Max duplicate read fraction (depth/length per position) (1)\n"
"-B          Turn off BAQ correction (false)\n"
"-R          Homopolymer repeat indel filtering (8)\n"
"-V          Output vcf format\n"
"-x CHR:POS  Output this pos only, then quit\n"
"-N FIL      Output noise stats to file\n"
"-S FIL      Read statistics and params from previous run with -s\n"
"-A ANNOT    Calculate in-target stats using the annotation file (requires -o)\n"
"-o PREFIX   Output prefix\n"
"\n"
"Input files\n"
"\n"
"Files must be sorted bam files with bai index files available.  Alternatively,\n"
"a single pileup file can be supplied.\n"
"\n"
"Output files\n"
"\n"
"Varcalls go to stdout.  Stats go to stdout, or stderr if varcalling too\n"
"\n"
"If an output prefix is used, files are created as follows:\n"
"   PREFIX.vstats    Pileup stats, and param hints\n"
"   PREFIX.var       Variant calls in tab delimited format\n"
"   PREFIX.vcf       Variant calls, in vcf format\n"
"   PREFIX.varsum    Summary of variant calls\n"
"   PREFIX.tgt       On-target stats detail\n"
"   PREFIX.tgtsum    Summary of on-target stats\n"
"   PREFIX.noise     Noise stats detail\n"
"\n"
"Stats Output:\n"
"\n"
"Contains mean, median, quartile information for depth, base quality, read len,\n"
"mapping quality, indel levels. Also estimates parameters suitable for\n"
"variant calls, and can be passed directly to this program for variant calls\n"
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
		   if (n > -1 && n < size) {
			   str.resize(n);
			   return str;
		   }
		   if (n > -1)
			   size=n+1;
		   else
			   size*=2;
       }
}

void to_upper(const std::string str) {
	std::string::iterator it;
	int i;
	for ( i=0;i<str.size();++i ) {
		((char *)(void *)str.data())[i]=toupper(((char *)str.data())[i]);
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

std::vector<char *> split(char* str,const char* delim)
{
    char* token = strtok(str,delim);
    std::vector<char *> result;
    while(token != NULL)
    {
        result.push_back(token);
        token = strtok(NULL,delim);
    }
    return result;
}

