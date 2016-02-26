/*
Copyright (c) 2012 Erik Aronesty

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
*/

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>

#include <gsl/gsl_randist.h>

#include <sys/stat.h>

#include <string>
#include <queue>
#include <list>

#include <sparsehash/sparse_hash_map> // or sparse_hash_set, dense_hash_map, ...
#include <sparsehash/dense_hash_map> // or sparse_hash_set, dense_hash_map, ...
#include "tidx/tidx.h"

#include "fastq-lib.h"

const char * VERSION = "0.96.819";

#define MIN_READ_LEN 20
#define DEFAULT_LOCII 1000000


using namespace std;
using namespace google;

void usage(FILE *f);


FILE *stat_fout=NULL;

// #define DEBUG 1

#define meminit(l) (memset(&l,0,sizeof(l)))
#ifdef DEBUG
    #define debug(s,...) fprintf(stderr,s,##__VA_ARGS__)
#else
    #define debug(s,...)
#endif
#undef warn
#define warn(s,...) (++errs, fprintf(stderr,s,##__VA_ARGS__))
#define die(s,...) (fprintf(stderr,s,##__VA_ARGS__), exit(1))
#define stat_out(s,...) fprintf(stat_fout,s,##__VA_ARGS__)
#define stdev(cnt, sum, ssq) sqrt((((double)cnt)*ssq-pow((double)sum,2)) / ((double)cnt*((double)cnt-1)))
#define log10(x) (log(x)/log(10))

//////////// BASE 2 INTEGER

#define T_A 0
#define T_C 1
#define T_G 2
#define T_T 3
#define T_SDEL 4
#define T_NDEL 5
#define T_INS 6
#define T_N 7
#define T_CNT 8
#define b2i(c) ((c)=='A'?0:(c)=='a'?0:(c)=='C'?1:(c)=='c'?1:(c)=='G'?2:(c)=='g'?2:(c)=='T'?3:(c)=='t'?3:(c)=='*'?4:(c)=='-'?5:(c)=='+'?6:7)
#define i2b(i) (i==0?'A':i==1?'C':i==2?'G':i==3?'T':i==4?'*':i==5?'-':i==6?'+':'?')
#define toupper(c) ((c)>='a' && (c)<='z' ? ((c)-32) : (c))

// die unless it's a number
int ok_atoi(const char *s) {
    if (!s || (!isdigit(*s)&& !(*s=='-'))) {
        die("%s is not a number\n", s);
    }
    return atoi(s);
}

double quantile(const std::vector<int> &vec, double p);
double quantile(const std::vector<double> &vec, double p);
double pnorm(double x);
double qnorm(double x);
int rand_round(double x);

// basic utils
std::vector<char *> split(char* str, char delim);
int split(char **buf, char* str, char delim);
std::string string_format(const std::string &fmt, ...);
void to_upper(const std::string str);
void rename_tmp(std::string f);

int errs=0;
extern int optind;
int g_lineno=0;
double vse_rate[T_CNT][T_CNT];

class Faidx {
public:
    typedef struct {
        int len;
        long long offset;
        int line_blen;
        int line_len;
    } Faient;

    sparse_hash_map<string, Faient> faimap;

    string fa_n;
    FILE *fa_f;

    Faidx() {fa_f=NULL;};
    void Load(const char *path);                    // open file, read fai

    // read into buffer
    bool Fetch(char *buf, const string &chr, int pos_from, int pos_to) {
        Fetch(buf, Chrdex(chr), pos_from, pos_to);
    };

    // read into buffer, with cached Chrdex
    bool Fetch(char *buf, const Faient *ent, int pos_from, int pos_to);
    const Faient * Chrdex(const string &chr) {
        return &(faimap[chr]);
    }
};

class Noise {
public:
	Noise() {noise=0;depth=0;};
	Noise(char r, char v, int d, double n, double q, double mq) {ref=r, var=v, depth=d; noise=n;qnoise=q;mnqual=mq;};
    char ref;
    char var;
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
    vcall() {base='\0'; mn_qual=mq0=fwd=rev=qual=is_ref=qual_ssq=mq_sum=mq_ssq=tail_rev=tail_fwd=fwd_q=rev_q=0; agreement=diversity=0.0;}
    char base;
	bool is_ref;
    int qual, fwd, rev, mq0, mn_qual, qual_ssq, mq_sum, mq_ssq, tail_rev, tail_fwd, fwd_q, rev_q;
    double diversity, agreement;
	vector <string> seqs;
    int depth() const {return fwd+rev;}
    int mq_rms() const {return sqrt(mq_ssq/depth());}
    int qual_rms() const {return sqrt(qual_ssq/depth());}
};

class vfinal {
public:
    vfinal(vcall &c) {max_idl_cnt=0; padj=1; pcall = &c;};
    vfinal & operator=(vfinal const&x) {max_idl_seq=x.max_idl_seq; max_idl_cnt=x.max_idl_cnt; padj=x.padj; pcall=x.pcall;}
    vcall *pcall;
    string max_idl_seq;
    int max_idl_cnt;
    double padj;
    bool is_indel() {return max_idl_cnt > 0;};
};

bool hitolocall (const vcall &i,const vcall &j) {return ((i.depth())>(j.depth()));}
bool sortreffirst (const vfinal &i,const vfinal &j) {return (i.pcall->is_ref&&!j.pcall->is_ref)||((i.pcall->is_ref==j.pcall->is_ref) && ((i.pcall->depth())>(j.pcall->depth())));}

class Read {
public:
    int MapQ;
    int Pos;
    string Seq;
    Read() {MapQ=Pos=0;};
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
    bool InTarget;
    int Regions;

    int SkipN;
    int SkipAmp;
    int SkipDupReads;
    int SkipMinMapq;
    int SkipMinQual;

    int RepeatCount;
    char RepeatBase;

	void Parse(char *line, PileupReads &reads, tidx *annot=NULL, char annot_type='\0');
    PileupSummary() { Base = '\0'; Pos=-1; };
};

class PileupManager;

class PileupSubscriber {
public:
    PileupManager *Manager;
    virtual void Visit(PileupSummary &dat) = 0;
    virtual void Finish() {};
    PileupSubscriber(PileupManager &man);
    PileupSubscriber() {Manager = NULL;}
    void SetManager(PileupManager &man);
};

class PileupManager  {
friend class PileupSubscriber;

private:
    void Visit(PileupSummary &dat);
    void VisitX(PileupSummary &dat, int windex);
    PileupSummary Pileup;

protected:
    vector<PileupSubscriber *> Kids;

public:

    string Reference;
    char InputType;
    int WinMax;             // flanking window size
    int WinDex;             // current index into the window (ususally midpoint)

    deque<PileupSummary> Win;

    int UseAnnot;
    tidx AnnotDex;          // start/stop index file
    char AnnotType;         // b (bed) or g (gtf - preferred)

    PileupReads Reads;
 
    PileupManager() {InputType ='\0'; WinMax=0; WinDex=0; UseAnnot=0; AnnotType='\0';}

    void Finish();

    void Parse(char *dat);

    void LoadAnnot(const char *annot_file);
    void FillReference(int refSize);
};

class VarStatVisitor : public PileupSubscriber {
    public:
    VarStatVisitor() : PileupSubscriber() {tot_locii=0; tot_depth=0; num_reads=0; stats.reserve(1000000); ins_stats.reserve(1000000); del_stats.reserve(1000000);};
    VarStatVisitor(PileupManager &man) : PileupSubscriber(man) {tot_locii=0; tot_depth=0; num_reads=0;};

    void Visit(PileupSummary &dat);
    void Finish() {};

	double tot_depth;
	int tot_locii;
	int num_reads;
	vector<Noise> stats;
	vector<Noise> ins_stats;
	vector<Noise> del_stats;
};


class VarCallVisitor : public PileupSubscriber {
    public:

    VarCallVisitor(PileupManager &man) : PileupSubscriber(man) {
        SkippedAnnot=0;
        SkippedDepth=0;
        Hets=0;
        Homs=0;
        Locii=0;
    };

    void Visit(PileupSummary &dat);

	int SkippedDepth;
	int SkippedAnnot;
	int Locii;
	int Hets;
	int Homs;
};

bool hasdata(const string &file) {
	struct stat st;
	if (stat(file.c_str(), &st)) {
		return false;
	}
	return st.st_size > 0;
}


int minsampdepth=20;
double pct_depth=0;
double pct_qdepth=0;
double global_error_rate=0;
double max_phred;
double vse_max_phred[T_CNT][T_CNT];
int total_locii=-1;
double pct_balance=0.05;       // at least 1 reverse read for every 20
char *debug_xchr=NULL;
int debug_xpos=0;
int debug_level=0;
int min_depth=1;
int min_mapq=0;
double min_diversity=0.15;     // only skip huge piles in one spot... even a little diversity is OK
double min_agreement=0.15;     // only really relevent when depths are high ... this is a bare minimum score
int min_qual=3;
int repeat_filter=7;
double artifact_filter=1;
int min_adepth=2;
int read_tail_pct=.6;
int read_tail_len=4;
int min_idepth=3;
int no_baq=0;
double zygosity=.5;		        // set to .1 for 1 10% admixture, or even .05 for het/admix
bool output_ref=0;              // set to 1 if you want to output reference-only positions
bool no_indels=0;

void parse_bams(PileupManager &v, int in_n, char **in, const char *ref);
void check_ref_fai(const char * ref);

FILE *noise_f=NULL, *var_f = NULL, *varsum_f = NULL, *tgt_var_f = NULL, *tgt_cse_f = NULL, *vcf_f = NULL, *eav_f=NULL, *cse_f=NULL;

double alpha=.05;
int phred=33;
double phi(double x);

FILE *openordie(const char *path, const char *mode) {
    FILE *f=fopen(path, mode);
    if (!f) {
        warn("Can't open-%s %s: %s\n", mode, path, strerror(errno));
        exit(1);
    }
    return f;
}

int str_in(const char *needle, const char **haystack) {
    int i=-1;
    while (*haystack) {
        ++i;
        if (!strcasecmp(needle, *haystack)) {
            return i;
        }
        ++haystack;
    }
    return -1;
}

void output_stats(VarStatVisitor &vstat);

Faidx faidx;
bool pcr_annot = false;

int main(int argc, char **argv) {
	char c;
	const char *ref=NULL;
	optind = 0;
	int umindepth=0;
	int uminadepth=-1;
	int uminidepth=0;
	double upctqdepth=0;
	int do_stats=0;
	int do_varcall=0;

    char *out_prefix = NULL;
    char *target_annot = NULL;
    const char *read_stats = NULL;


// list of default output formats used when -o is specified
#define MAX_F 20
    const char *format_list[MAX_F]={"var", "eav", "noise", "varsum", NULL};

// option characters... use \1, \2... for long-only

    #define OPT_PCR_ANNOT '\1'
    #define OPT_DEBUG_LEVEL '\2'
    #define OPT_NO_INDELS '\3'
    #define OPT_FILTER_ANNOT 'A'

// long options
    static struct option long_options[] = {
       {"pcr-annot", 1, 0, OPT_PCR_ANNOT},
       {"filter-annot", 1, 0, OPT_FILTER_ANNOT},
       {"repeat-filter", 1, 0, 'R'},
       {"no-indels", 0, 0, OPT_NO_INDELS},
       {"agreement", 1, 0, 'G'},
       {"diversity", 1, 0, 'd'},
       {"version", 0, 0, 'V'},
       {"debug", 1, 0, OPT_DEBUG_LEVEL},
       {0, 0, 0, 0}
    };

	while ( (c = getopt_long(argc, argv, "?sv0VBhe:m:x:f:p:a:g:q:Q:i:o:D:R:b:L:S:F:A:G:d:",long_options,NULL)) != -1) {
		switch (c) {
			case OPT_PCR_ANNOT: target_annot=optarg; pcr_annot=true; break;
			case OPT_FILTER_ANNOT: target_annot=optarg; pcr_annot=false; break;
			case OPT_NO_INDELS: no_indels=true; break;
			case 'h': usage(stdout); return 0;
			case 'm': umindepth=ok_atoi(optarg); break;
			case 'q': min_qual=ok_atoi(optarg); break;
			case 'o': out_prefix=optarg; break;
			case 'Q': min_mapq=ok_atoi(optarg); break;
			case 'V': printf("Version: %s\n", VERSION); exit(0); break;
			case 'R': repeat_filter=ok_atoi(optarg); break;
			case 'a': uminadepth=ok_atoi(optarg);break;
			case 'D': artifact_filter=atof(optarg);break;
			case 'i': uminidepth=ok_atoi(optarg);break;
			case 'd': min_diversity=atof(optarg); break;
			case 'G': min_agreement=atof(optarg); break;
			case '0': min_qual=0; umindepth=0; min_mapq=0; repeat_filter=0; uminadepth=0; artifact_filter=0; uminidepth=0; min_diversity=0; alpha=1; pct_balance=0; upctqdepth=0; min_agreement=0; break;
			case 'x': {
					debug_xchr=optarg;
					char *p=strrchr(debug_xchr, ':');
					if (!p) die("Invalid param for -x");
					*p='\0';
					debug_xpos=atoi(++p);
					if (!p) die("Invalid param for -x, need pos");
                    ++debug_level;
					break;
				}
			case OPT_DEBUG_LEVEL: 
                    debug_level = atoi(optarg);
			case 'b': pct_balance=atof(optarg)/100.0; break;
			case 'B': no_baq=1; break;
			case 'p': upctqdepth=atof(optarg); break;
			case 'e': alpha=atof(optarg); break;
			case 'g': global_error_rate=atof(optarg); break;
			case 'L': total_locii=ok_atoi(optarg); break;
			case 'f': ref=optarg; break;
			case 's': do_stats=1; break;
			case 'S': read_stats=optarg; break;
			case 'v': do_varcall=1; break;
			case 'F': {
                char *tok, *saved; int i=0;
                for (tok = strtok_r(optarg, "%", &saved); tok && i < MAX_F; tok = strtok_r(NULL, " ,", &saved)) {
                    format_list[i++]=tok;
                }
                format_list[i]=NULL;
                break;
            }
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


	if (!do_stats && !do_varcall) {
		warn("Specify -s for stats only, or -v to do variant calling\n\n");
		usage(stderr);
		return 1;
	}

    if (out_prefix && do_varcall) {
        var_f = openordie(string_format("%s.var.tmp", out_prefix).c_str(), "w");

        fprintf(var_f,"%s\t%s\t%s\t%s\t%s\t%s\t%s%s\n","chr", "pos", "ref", "depth", "skip", "pct", (target_annot&&!pcr_annot) ? "target\t" : pcr_annot ? "regions\t" : "", "...");

        varsum_f = openordie(string_format("%s.varsum.tmp", out_prefix).c_str(), "w");

        if (target_annot && !pcr_annot) {
            // targted only output
            tgt_var_f = openordie(string_format("%s.tgt.var.tmp", out_prefix).c_str(), "w");
            fprintf(tgt_var_f,"%s\t%s\t%s\t%s\t%s\t%s\t%s\n","chr", "pos", "ref", "depth", "skip", "pct", "...");
        }
 
        if (str_in("vcf", format_list)>=0) {
            vcf_f = openordie(string_format("%s.vcf.tmp", out_prefix).c_str(), "w");
        }
        if (str_in("eav", format_list)>=0) {
            eav_f = openordie(string_format("%s.eav.tmp", out_prefix).c_str(), "w");
        }

        if (str_in("cse", format_list)>=0) {
            check_ref_fai(ref);
            cse_f = openordie(string_format("%s.cse.tmp", out_prefix).c_str(), "w");
            faidx.Load(ref);
            // targted only output
            if (target_annot && ! pcr_annot) 
                tgt_cse_f = openordie(string_format("%s.tgt.cse.tmp", out_prefix).c_str(), "w");
        }
    } else {
        var_f = stdout;
        varsum_f = stderr;
    } 

	if (umindepth > minsampdepth) {
		minsampdepth=umindepth;
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

    // not really random
    srand(1);

    meminit(vse_rate);

	if (do_stats) {
        if (out_prefix) {
            stat_fout = openordie(string_format("%s.stats", out_prefix).c_str(), "w");
            noise_f = openordie(string_format("%s.noise", out_prefix).c_str(), "w");
		    fprintf(noise_f,"%s\t%s\t%s\t%s\t%s\t%s\n", "depth", "ref", "var", "noise", "qnoise", "qmean");
        }
        if (!stat_fout)  {
            if (do_varcall)
    		    stat_fout=stderr;			// stats to stderr
            else
    		    stat_fout=stdout;			// stats to stdout
        }
        // do stats by myself
        PileupManager pman;
		VarStatVisitor vstat(pman);
		parse_bams(pman, in_n, in, ref);
        output_stats(vstat);
        if (out_prefix) {
            fclose(stat_fout);
            fclose(noise_f);
            noise_f=NULL;
            stat_fout=NULL;
        }
    }
 
    if (do_varcall && out_prefix) {
        // run stats again... with new levels if any
        stat_fout = openordie(string_format("%s.vstats", out_prefix).c_str(), "w");
    }

    if (read_stats){
        FILE * f = fopen(read_stats, "r");
        if (!f) {
            warn("File %s does not exist, quitting\n", read_stats);
            exit(1);
        }
        line l; meminit(l);
        char *val;

        double noise_mean=0;
        double noise_dev=0;

        double qnoise_mean=0;
        double qnoise_dev=0;

        double vse_mean[T_CNT][T_CNT];
        double vse_dev[T_CNT][T_CNT];

        meminit(vse_mean);
        meminit(vse_dev);

        while(read_line(f, l)>0) {
            if (val=strchr(l.s, '\t')) {
                *val='\0'; ++val;
                if (!strcasecmp(l.s, "min depth")) {
                    if (umindepth && umindepth > atoi(val)) {
			            fprintf(varsum_f,"warning\tsampling depth was less than variation depth\n");
                    }
                    if (!umindepth) umindepth=atoi(val); 
                } else if (!strcasecmp(l.s, "min pct qual")) {
                    if (upctqdepth<=0) upctqdepth=atof(val); 
                } else if (!strcasecmp(l.s, "noise mean")) {
                    if (global_error_rate<=0) noise_mean=atof(val); 
                } else if (!strcasecmp(l.s, "noise dev")) {
                    if (global_error_rate<=0) noise_dev=atof(val); 
                } else if (!strcasecmp(l.s, "qnoise mean")) {
                    if (global_error_rate<=0) qnoise_mean=atof(val); 
                } else if (!strcasecmp(l.s, "qnoise dev")) {
                    if (global_error_rate<=0) qnoise_dev=atof(val); 
                } else if (!strcasecmp(l.s, "locii >= min depth")) {
                    if (total_locii<0) total_locii=atoi(val); 
                } else if (!strcasecmp(l.s, "alpha")) {
                    if (alpha<=0) alpha=atof(val); 
                } else if (!strncasecmp(l.s, "vnoise", 6)) {
                    char ref, var; char typ[10];
                    if ( sscanf(l.s, "vnoise %s %c:%c", typ, &ref, &var) == 3) {
                        if (*typ == 'm') {
                            // warn("vse_mean  %d, %d : %f\n", b2i(ref),b2i(var), atof(val));
                            vse_mean[b2i(ref)][b2i(var)]=atof(val); 
                        } else if (*typ == 'd') {
                            vse_dev[b2i(ref)][b2i(var)]=atof(val);
                        } else {
                            die("Invalid stats format : %s\n", l.s);
                        }
                    }
                }
            }
        }

        int i, j;
        for (i=0;i<T_CNT;++i) {
            for (j=0;j<T_CNT;++j) {
                if (vse_mean[i][j]>0) {
                    vse_rate[i][j] = vse_mean[i][j] + vse_dev[i][j];
                } else {
                    vse_rate[i][j] = noise_mean+noise_dev;
                }
            }
        }

        if (noise_mean > 0) {
            global_error_rate=noise_mean+noise_dev;
        }
    }

    // for speed, do this once...
    max_phred = -log10(global_error_rate)*10;
    meminit(vse_max_phred);

    // init table
    int i, j;
    for (i=0;i<T_CNT;++i) {
        for (j=0;j<T_CNT;++j) {
            if (vse_rate[i][j]>0) {
                vse_max_phred[i][j]=-log10(vse_rate[i][j])*10;
            } else {
                vse_max_phred[i][j]=max_phred;
            }
            // warn("%d %d %f\n", i, j, vse_max_phred[i][j]);
        }
    }
 
    if (total_locii<0) total_locii=DEFAULT_LOCII;
    if (total_locii==0) total_locii=1;          // no adjustment

    if (eav_f) {
        fprintf(eav_f,"chr\tpos\tref\tdepth\tnum_states\ttop_consensus\ttop_freq\tvar_base\tvar_depth\tvar_qual\tvar_strands\tforward_strands\treverse_strands\t%cval\tdiversity\tagreement\t%s\n", (total_locii>1?'e':'p'), (target_annot&&!pcr_annot) ? "in_target\t" : pcr_annot ? "regions\t" : "");
    }

	if (do_varcall) {
		if (umindepth) min_depth=umindepth;
		if (upctqdepth > 0) pct_qdepth=(double)upctqdepth/100;
		if (uminadepth>=0) min_adepth=uminadepth;
		if (uminidepth) min_idepth=uminidepth;

		if (!min_depth || (!pct_depth  && !pct_qdepth)) {
			fprintf(varsum_f,"warning\toutputting all variations, no minimum depths specified\n");
		}

        if (pct_qdepth==0.0 && !min_adepth) {
            output_ref=1;
        }

		fprintf(varsum_f,"version\tvarcall-%s\n", VERSION);
		fprintf(varsum_f,"min depth\t%d\n", min_depth);
		fprintf(varsum_f,"min call depth\t%d\n", min_adepth);
		fprintf(varsum_f,"alpha\t%f\n", alpha);
		fprintf(varsum_f,"min pct qual\t%d\n", (int)(100*pct_qdepth));

		fprintf(varsum_f,"min balance\t%d\n", (int)(100*pct_balance));
		fprintf(varsum_f,"artifact filter\t%f\n", artifact_filter);
		fprintf(varsum_f,"min qual\t%d\n", min_qual);
		fprintf(varsum_f,"min map qual\t%d\n", min_mapq);
		fprintf(varsum_f,"error rate\t%f\n", global_error_rate);
		fprintf(varsum_f,"locii used for adjustment\t%d\n", total_locii);

        PileupManager pman;

		VarCallVisitor vcall(pman);
		VarStatVisitor vstat;

		if (stat_fout) {
            vstat.SetManager(pman);
        }
       
        if (target_annot) {
            pman.LoadAnnot(target_annot);
       }

        if (cse_f) {
            pman.WinMax=21;
        } else if (repeat_filter > 0) {
		    fprintf(varsum_f,"homopolymer filter\t%d\n", repeat_filter);
            pman.WinMax=repeat_filter+repeat_filter+3;
        } else {
            pman.WinMax=5;
        }

        if (vcf_f) {
            // print VCF header
            fprintf(vcf_f, "%s\n", "##fileformat=VCFv4.1");
        }
        if (cse_f) {
            fprintf(cse_f, "Chr\tPos\tRef\tA\tC\tG\tT\ta\tc\tg\tt\tAq\tCq\tGq\tTq\taq\tcq\tgq\ttq\tRefAllele\tAd\tCd\tGd\tTd\tAg\tCg\tGg\tTg%s\n", pcr_annot ? "\tRegions" : "");
        }

		parse_bams(pman, in_n, in, ref);

        if (pman.InputType == 'B') {
        	fprintf(varsum_f,"baq correct\t%s\n", (no_baq?"no":"yes"));
        }
        fprintf(varsum_f,"locii\t%d\n", vcall.Locii);
        fprintf(varsum_f,"hom calls\t%d\n", vcall.Homs);
        fprintf(varsum_f,"het calls\t%d\n", vcall.Hets);
        fprintf(varsum_f,"locii below depth\t%d\n", vcall.SkippedDepth);
        fprintf(varsum_f,"locii outside annot\t%d\n", vcall.SkippedAnnot);

        if (out_prefix) {
            // close it all
            fclose(var_f);
            fclose(varsum_f);
            if (vcf_f) fclose(vcf_f);
            if (eav_f) fclose(eav_f);
            if (noise_f) fclose(noise_f);
            if (cse_f) fclose(cse_f);
            if (tgt_var_f) fclose(tgt_var_f);
            if (tgt_cse_f) fclose(tgt_cse_f);

            rename_tmp(string_format("%s.var.tmp", out_prefix));
            rename_tmp(string_format("%s.varsum.tmp", out_prefix));

            if (vcf_f) rename_tmp(string_format("%s.vcf.tmp", out_prefix));
            if (eav_f) rename_tmp(string_format("%s.eav.tmp", out_prefix));
            if (cse_f) rename_tmp(string_format("%s.cse.tmp", out_prefix));

            if (tgt_var_f) rename_tmp(string_format("%s.tgt.var.tmp", out_prefix));
            if (tgt_cse_f) rename_tmp(string_format("%s.tgt.cse.tmp", out_prefix));

            if (stat_fout) 
                output_stats(vstat);
        }
	}
}

void rename_tmp(std::string f) {
    std::string notmp = f;
    size_t pos = notmp.find(".tmp");
    if (pos >= 0) {
        notmp.replace(notmp.find(".tmp"),4,""); 
        rename(f.c_str(),notmp.c_str());
    }
}

// normal distribution
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

void parse_bams(PileupManager &v, int in_n, char **in, const char *ref) {

	if (!in_n) {
		die("No input files, quitting\n");
	}

	int i, bam_n=0;
	for (i=0;i<in_n;++i) {
		if (!strcmp(fext(in[i]), ".bam")) {
			++bam_n;
		}
	}

	if (bam_n != in_n) {
		if (bam_n > 0) {
			die("Can't mix bams and other input files\n");
		} else {
			if (in_n > 1) {
				die("Can't handle multiple pileups... TODO\n");
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
        check_ref_fai(ref);

		const char *nobaq = no_baq ? "-B" : "";

		string mpil_cmd = string_format("samtools mpileup -Q 0 -d 100000 %s -f '%s'", nobaq, ref);

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
        if (!strcmp(in[0], "-")) {
            fin=stdin;
        } else {
            if (!strcmp(fext(in[0]), ".gz")) {
                string gunz = string_format("gunzip -c '%s'", in[0]);
                fin = popen(gunz.c_str(), "r");
                is_popen = 1;
            } else {
                fin = fopen(in[0], "r");
            }
            if (!fin) {
                warn("%s: %s", in[0], strerror(errno));
                exit(1);
            }
        }
	}

    line l; meminit(l);
    g_lineno=0;
    if (fin) {
        while(read_line(fin, l)>0) {
            ++g_lineno;
        //	chr      2       G       6       ^9,^+.^*,^2,^&.^&,      &.'&*-  9+*2&&  166,552,643,201,299,321
            v.Parse(l.s);
        }
        v.Finish();

        if (is_popen) pclose(fin); else fclose(fin);
    }

	if (g_lineno == 0) {
		warn("No data in pileup, quitting\n");
		exit(1);
	}
}


bool hitoloint (int i,int j) { return (i>j);}

class q_calls {public: q_calls() {meminit(call);} int call[8];};
vector<int> depthbypos;
vector<q_calls> depthbyposbycall;
typedef struct  {
    string Chr;
    int Beg;
    int End;
} ChrRange;


char *_dat[256];
inline void PileupSummary::Parse(char *line, PileupReads &rds, tidx *adex, char atype) {

	int dsize=split(_dat, line, '\t');

	if (dsize < 6) {
		warn("Can't read pileup : %d fields, need 6 columns, line %d\n", (int) dsize, g_lineno);
		exit(1);
	}

	const char * p_qual=_dat[5];

	Chr=_dat[0];
	Pos=atoi(_dat[1]);
	Base=*(_dat[2]);
	Depth = atoi(_dat[3]);
	SkipDupReads = 0;
	SkipN = 0;
	SkipAmp = 0;
	SkipMinQual = 0;
	SkipMinMapq = 0;
	RepeatCount = 0;
	RepeatBase = '\0';
	NumReads = 0;
    InTarget = 0;
    Regions = 0;

	int i;

	const char *cur_p = _dat[4];

    list<Read>::iterator read_i = rds.ReadList.begin();
    
    memset(depthbypos.data(),0,depthbypos.size()*sizeof(depthbypos[0]));
    memset(depthbyposbycall.data(),0,depthbyposbycall.size()*sizeof(depthbyposbycall[0]));

    int eor=0;

    // list of amplicon range objects
    vector<ChrRange> amps;

    if (pcr_annot && adex) {
        string s = adex->lookup(Chr.data(), Pos + (atype=='b' ? -1 : 0), "^");
        if (s.length()) {
            vector<char *> a=split((char *)s.data(), '^');
            Regions=a.size()-1;
            // skip leading entry...
            for(i=1;i<a.size();++i) {
                vector<char *> f=split(a[i], '\t');
                // create new range object
                if (f.size() >= 3) {
                    ChrRange amp;
                    amp.Chr=f[0];
                    amp.Beg=atoi(f[1]);
                    amp.End=atoi(f[2]);
                    if (atype=='b') {
                        ++amp.Beg;
                    } 
                    if ((amp.End < amp.Beg) || !amp.Beg) {
                        die("Annotation file must be in bed or gtf format, or at least a 1-based inclusive set of ranges\n"); 
                    }
//                    warn("AMP: %s:%d-%d\n",amp.Chr.data(),amp.Beg,amp.End);
                    amps.push_back(amp);
                }
            }
        }
    }

    if (debug_xpos) {
        if (Pos == debug_xpos && !strcmp(debug_xchr,Chr.data())) {
            fprintf(stderr,"depth: %d, meanreadlen: %f\n", Depth, rds.MeanReadLen());
        }
    }

    int meanreadlen = rds.MeanReadLen();
    int maxdepthbypos = meanreadlen <= 0 ? 10 : max(10, round(10.0 * artifact_filter * (Depth/(double)meanreadlen)));

    Calls.clear();

    int j;
    int pia_len=0;
	for (i=0;i<Depth;++i,++read_i) {
		bool sor=0;

		if (*cur_p == '^') {
			sor=1;
			++cur_p;
            Read x;
            x.MapQ = *cur_p-phred;
            x.Pos = Pos;
            ++cur_p;
            if (read_i != rds.ReadList.end()) {
                ++read_i;
            }
            read_i=rds.ReadList.insert(read_i,x);
            meanreadlen = rds.MeanReadLen();
		}

        if (read_i == rds.ReadList.end()) {
            warn("warning\tread start without '^', partial pileup: '%s'\n", cur_p);
            Read x;
            x.MapQ = 0;
            x.Pos = -1;
            read_i=rds.ReadList.insert(read_i,x);
            meanreadlen = rds.MeanReadLen();
        }


        // position of read relative to my position
        int pia = read_i->Pos >= 0 ? Pos-read_i->Pos : 0;

        pia = pia % (meanreadlen*2);

		if (pia >= depthbypos.size()) {
			depthbypos.resize(pia+1);
			depthbyposbycall.resize(pia+1);
		}
        if (pia >= pia_len)
            pia_len=pia+1;

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

		if (o == '>' || o == '<') {
            /// splice!  don't add to SEQ or depth, etc!	
			c = 'N';					// no call
			is_ref = 1;
		}

		bool skip = 0;
        bool ampok = !pcr_annot || !adex;

        if (!ampok) {
            for (j=0;j<amps.size();++j) {
                int apos = read_i->Pos + meanreadlen + 1;
                int bpos = read_i->Pos + meanreadlen;
                int cpos = read_i->Pos + meanreadlen - 1;
                if (apos == amps[j].End || bpos == amps[j].End || cpos == amps[j].End) {
                    ampok=1;
                }
                if (read_i->Pos == amps[j].Beg) {
                    ampok=1;
                    break;
                }
            }
        }


        if (debug_xpos) {
            if (debug_level >= 3) {
                if (Pos == debug_xpos && !strcmp(debug_xchr,Chr.data())) {
                    fprintf(stderr, "DEBUG: PIA: %d, DBP: %d, rrou: %d, nonRR: %f, maxdbp: %d, filt: %d, f1: %d, f2: %d\n", pia,
                        depthbypos[pia], max(1,rand_round(0.5 + artifact_filter * (Depth/rds.MeanReadLen()))), artifact_filter * (Depth/rds.MeanReadLen()),
                        maxdepthbypos, (10*(depthbypos[pia]-1))+(i%10), ((10*(depthbypos[pia]-1))+(i%10)) > maxdepthbypos, depthbypos[pia] > max(1,rand_round(0.5+artifact_filter * (Depth/rds.MeanReadLen()))) );
                }
            }
        }


// before dup filtering

        int call_index = b2i(c);
        depthbyposbycall[pia].call[call_index]++;

		if (!ampok) {
//            warn("SKIP: %d-%d, %c\n", read_i->Pos, (int)(read_i->Pos+rds.MeanReadLen()), o);
            ++SkipAmp;
            skip=1;
		} else if (c == 'N') {
			++SkipN;
			skip=1;
            // ok, instead of rolling a random number, even things out
		} else if (mq < min_mapq) {
			++SkipMinMapq;
			skip=1;
		} else if (q < min_qual) {
			++SkipMinQual;
			skip=1;
		} else if (artifact_filter > 0 && (((10*(depthbypos[pia]-1))+(i%10)) > maxdepthbypos )) {
			++SkipDupReads;
			skip=1;
		} else {
			int j = call_index;

			if (call_index >= Calls.size()) {
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
				Calls[j].rev_q += q;
			} else if ( c != 'N' ) {
				++Calls[j].fwd;
				Calls[j].fwd_q += q;
			}

			Calls[j].qual+=q;
            Calls[j].mn_qual+=min(mq,q);
            Calls[j].mq_ssq+=mq*mq;
            Calls[j].mq_sum+=mq;
            Calls[j].qual_ssq+=q*q;
/*
            if (pia <= read_tail_len || (rds.MeanReadLen()-pia) <= read_tail_len) {
                if ( o == ',' || o == 'a' || o == 'c' || o == 't' || o == 'g' ) {
                    ++Calls[j].tail_rev;
                } else {
                    ++Calls[j].tail_fwd;
                }
            }
*/

            if (vcf_f) {
                if (mq == 0) 
                    Calls[j].mq0++; 
            }
		}

		if (c == '-' || c == '+') {
            warn("invalid pileup, at '%s', indel not attached to read?\n", cur_p);
		} else {
		    if (c != '*' && c != 'N') 
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
                    Calls[j].rev_q+=q;
                } else {
                    ++Calls[j].fwd;
                    Calls[j].fwd_q+=q;
                }
                Calls[j].qual+=q;
                Calls[j].mn_qual+=min(q, mq);
                Calls[j].qual_ssq+=q*q;
                Calls[j].mq_ssq+=mq*mq;
                Calls[j].mq_sum+=mq;
                Calls[j].seqs.push_back(ins_seq);
            }
            cur_p=end_p+len;
        }

        if (*cur_p == '$') {
            if (read_i->MapQ > -1) {
                rds.TotReadLen+=read_i->Seq.size();
                rds.ReadBin.push_back(*read_i);
                if (rds.ReadBin.size() > min(1000,Depth*2)) {
                    rds.TotReadLen-=rds.ReadBin.front().Seq.size();
                    rds.ReadBin.pop_front();
                }
            }
//            printf("%d\t%s\n", read_i->MapQ, read_i->Seq.c_str());
            read_i=rds.ReadList.erase(read_i);
            meanreadlen = rds.MeanReadLen();
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
		warn("Failed to parse pileup %s\n", _dat[4]);
		exit(1);
	}

	Depth=0;
	for (i=0;i<5 && i < Calls.size();++i) {		// total depth (exclude inserts for tot depth, otherwise they are double-counted)
		Depth+=Calls[i].depth();
	}

    if (debug_xpos) {
        if (Pos == debug_xpos && !strcmp(debug_xchr,Chr.data())) {
            fprintf(stderr,"xpos-max-per-pos\t%f\n", maxdepthbypos/10.0);
            fprintf(stderr,"xpos-mean-readlen\t%f\n", rds.MeanReadLen());
            fprintf(stderr,"xpos-depth-list\t");
            for (i=0;i<pia_len;++i) {
                fprintf(stderr,"%d,", depthbypos[i]);
            }
            fprintf(stderr,"\n");

            if (debug_level >= 2) {
                for (i=0;i < Calls.size();++i) {		// total depth (exclude inserts for tot depth, otherwise they are double-counted)
                    if (Calls[i].depth()>0) {
                        fprintf(stderr,"xpos-depth-list-%c\t",Calls[i].base);
                        for (j=0;j<depthbyposbycall.size();++j) {
                            fprintf(stderr,"%d,", depthbyposbycall[j].call[i]);
                        }
                        fprintf(stderr,"\n");
                    }
                } 
            }
        }
    }

//    bool quit=0;
	for (i=0;i < Calls.size();++i) {		// total depth (exclude inserts for tot depth, otherwise they are double-counted)
        if (Calls[i].depth()>0) {
            double expected=Calls[i].depth()/meanreadlen;
            double all_pct=Calls[i].depth()/(double)Depth;
            double p2=(Calls[i].depth()+2*all_pct*pia_len)/(double)(Depth+2*pia_len);
            double total_v=0;
            double diff;
            double moment_den=0;
            double p1;
            double pdiff;
            double cube_v=0;
            int j;

/*
            if(Depth>100 && Calls[i].depth() < 100) {
                printf("exp: %g, depth: %d, cdep: %d, p2: %g\nv=c(", expected, Depth, Calls[i].depth(), p2);
                for (j=0;j<depthbyposbycall.size();++j) {
                    printf("%d,", depthbyposbycall[j].call[i]);
                }
                printf(")\nd=c(");
                for (j=0;j<depthbyposbycall.size();++j) {
                    printf("%d,", depthbypos[j]);
                }
                printf(")\n");
            }

*/
            int poscnt=0;
            for (j=0;j<pia_len;++j) {
                // diversity calc
                diff=floor(fabs(depthbyposbycall[j].call[i]-expected));
                total_v+=diff*diff;
                if (depthbyposbycall[j].call[i] > 0) {
                    ++poscnt;
                }

                p1=((depthbyposbycall[j].call[i]+all_pct*2)/((double)depthbypos[j]+2));
                pdiff=fabs(p1-p2);
/*
                if (Pos == debug_xpos && !strcmp(debug_xchr,Chr.data())) {
                    warn("base:%c, depth:%d, dbc:%d, p1: %g, pdiff: %g, maxc: %g, all_pct: %g\n", Calls[i].base, depthbypos[j], depthbyposbycall[j].call[i], p1, pdiff,  max(depthbypos[j]*pdiff*pdiff*pdiff-all_pct,0), all_pct);
                }
*/
                cube_v += max(depthbypos[j]*pdiff*pdiff*pdiff-all_pct,0);
            }
            double shift_v = max(0, total_v-2*Calls[i].depth());
            Calls[i].diversity = max(0,1-shift_v/max(1,(pow(Calls[i].depth()-expected,2)-2*Calls[i].depth())));
            if (poscnt==1) Calls[i].diversity = 0;

            double wt4_od=pow(cube_v/((double)(Depth+2*pia_len)),1.0/3.0)/all_pct;
            Calls[i].agreement = max(0,1-wt4_od);

            if (debug_xpos) {
                if (Pos == debug_xpos && !strcmp(debug_xchr,Chr.data())) {
                    if (debug_level > 2) warn("base:%c, cube_v:%g, deno:%g\n", Calls[i].base, cube_v, ((double)(Depth+2*pia_len)));
                    fprintf(stderr,"xpos-agree-%c\t%g\n",Calls[i].base, Calls[i].agreement);
                    fprintf(stderr,"xpos-diver-%c\t%g\n",Calls[i].base, Calls[i].diversity);
                }
            }

/*
            if(Depth>100 && Calls[i].depth() < 100) {
                printf("num %g, den %g\n", (p2*(1-p2)*(Depth+2*pia_len)), moment_den);
                printf(", Agree: %g", Calls[i].agreement);
                printf(", Divers: %g\n", Calls[i].diversity);
                quit=1;
            }
*/
        }
	}
//    if(quit) die("\nDIED\n");

	TotQual=0;
	for (i=0;i<5 && i < Calls.size();++i) {		// total depth (exclude inserts for tot depth, otherwise they are double-counted)
		TotQual+=Calls[i].qual;
	}
}

PileupSummary JunkSummary;

inline void PileupManager::Parse(char *dat) {
    Pileup.Parse(dat, Reads, UseAnnot ? &AnnotDex : NULL, AnnotType);
    Visit(Pileup);
}

void PileupManager::Visit(PileupSummary &p) {

    if (WinMax < 3) {
        // no real window ... just go straight
        VisitX(p, -1);
        return;
    }

    if (p.Base != '-' && p.Base != '@') {
        if (Win.size() && (Win.back().Pos != (p.Pos - 1) )) {
            if (Win.back().Pos < p.Pos && ((p.Pos - Win.back().Pos) <= (WinMax/2))) {
                while (Win.back().Pos < (p.Pos - 1)) {
                    // visit/pop, add a placeholder
                    JunkSummary.Base = '-';
                    JunkSummary.Pos = Win.back().Pos + 1;
                    Visit(JunkSummary);
                }
            } else {
                while (Win.size() && Win[WinMax/2].Base != '@') {
                    // visit/pop, but don't add anything, until it's empty
                    JunkSummary.Base = '@';
                    JunkSummary.Pos = 0;
                    Visit(JunkSummary);
                }
            }
        }
    }

    // initialize the window with nothing, if it's not full
    while (Win.size() < WinMax) {
        JunkSummary.Base = '@';
        JunkSummary.Pos = 0;
        Win.push_back(JunkSummary);
    }

    Win.push_back(p);

    //debug("Visit: %d\n", p.Pos);

    if (Win.size() > WinMax)        // queue too big?  pop
        Win.pop_front();
    
    int i;
    int lrc=0,rrc=0;                // left repeat count, right repeat count
    char lrb, rrb;                  // left repeat base...
    int vx;

    if (Win.size() < WinMax) {    // small window?  look at leading edge only
        return;
    } else {
        vx = WinMax/2;              // larger window? look at midpoint
    }

    if (Win[vx].Base == '-' || Win[vx].Base == '@') 
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
        for (i=vx+2; i < Win.size(); ++i) {
            if (Win[i].Base == rrb)
                ++rrc;
            else
                break;
        }
    }

    // repeat counts are now 1-based, not 0-based
    ++lrc;
    ++rrc;

    // maximum repeat count and associated base
    if (lrb == rrb ) {
        Win[vx].RepeatCount = lrc+rrc;
        Win[vx].RepeatBase = lrb;
    } else if (lrc > rrc) {
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

    double drms = 0; 
    if (vx < Win.size()-1) {
        int i;
		int dminus = b2i('-');
		int dstar = b2i('*');

        if (Win[vx].Calls.size() > dminus && Win[vx].Calls[dminus].depth() > 0) {
            if (Win[vx+1].Calls.size() > dstar && Win[vx+1].Calls[dstar].depth() > 0) {
                // baq adjustment works at the 'star' not at the 'indel', so adjust qual using the next locus
               double adj=Win[vx+1].Calls[dstar].qual_rms()/(double)Win[vx].Calls[dminus].qual_rms();
               if (debug_xpos) {
                    if (Win[vx].Pos == debug_xpos && !strcmp(debug_xchr,Win[vx].Chr.data())) {
                        fprintf(stderr,"xpos-adj-qual\t%d to %d (%f)\n", Win[vx].Calls[dminus].qual_rms(),Win[vx+1].Calls[dstar].qual_rms(), adj);
                    }
               }
               Win[vx].Calls[dminus].qual *= adj; 
               Win[vx].Calls[dminus].qual_ssq *= adj;
            } else {    
                vcall none;
                if (debug_xpos) {
                    if (Win[vx].Pos == debug_xpos && !strcmp(debug_xchr,Win[vx].Chr.data())) {
                        fprintf(stderr,"xpos-skip-del-qual\t%d\n", Win[vx].Calls[dminus].depth());
                    }
                }
                Win[vx].Calls[dminus] = none;
            }
        }
    }

    VisitX(Win[vx], vx);
}

void PileupManager::Finish() {
    // finish out the rest of the pileup, with the existing window
    int vx = WinMax/2+1;
    while (vx < Win.size()) {
        ///debug("Finish: %d\n", Win[vx].Pos);
        VisitX(Win[vx++], vx);
    }
    int i;
    for (i=0;i<Kids.size();++i) {
        Kids[i]->Finish();
    }
}

void PileupManager::VisitX(PileupSummary &p, int windex) {

    WinDex=windex;

    if (UseAnnot) {
        // index lookup only.... not string lookup
        const std::vector <long int> * v = &(AnnotDex.lookup(p.Chr.data(), p.Pos + (AnnotType=='b' ? -1 : 0)));
        if (v && v->size()) {
            p.InTarget=1;
        }
    }

    if (no_indels) {
        if (p.Calls.size() > 4) {
            p.Calls.resize(4);
        }
    }

    int i;
    for (i=0;i<Kids.size();++i) {
        Kids[i]->Visit(p);
    }
}

void VarCallVisitor::Visit(PileupSummary &p) {
    //debug("VisitX: %d\n", p.Pos);

	if (debug_xpos) {
		if (p.Pos != debug_xpos)
			return;
		if (strcmp(debug_xchr,p.Chr.data())) 
			return;
	}

    if (pcr_annot) {
        if (!p.InTarget) {
            if (debug_xpos) {
                fprintf(stderr,"xpos-skip-annot\t1\n");
            }
            ++SkippedAnnot;
            return;
        }
    }

	if (p.Depth < min_depth) {
        if (debug_xpos) {
            fprintf(stderr,"xpos-skip-depth\t%d < %d\n",p.Depth, min_depth);
		    fprintf(stderr,"xpos-skip-dup\t%d\n",p.SkipDupReads);
		    fprintf(stderr,"xpos-skip-n\t%d\n",p.SkipN);
		    fprintf(stderr,"xpos-skip-amp\t%d\n",p.SkipAmp);
		    fprintf(stderr,"xpos-skip-mapq\t%d\n",p.SkipMinMapq);
		    fprintf(stderr,"xpos-skip-qual\t%d\n",p.SkipMinQual);
			exit(0);
        }
		++SkippedDepth;
		return;
	}

	int ins_fwd = p.Calls.size() > 6 ? p.Calls[6].fwd : 0;
	int ins_rev = p.Calls.size() > 6 ? p.Calls[6].rev : 0;

	int i;
	if (p.Calls.size() > 6) 
		p.Calls.resize(7);	// toss N's before sort

    static char regions[64] = "";
    if (pcr_annot) {
        sprintf(regions, "\t%d", p.Regions); 
    } 

    // OUTPUT CSE BEFORE REORDRED BASES!
    if (cse_f) {
        if (p.Calls.size() < 4) 
            p.Calls.resize(4);	// cse needs 4 calls

        Manager->FillReference(21);
    
        // cse format... no need to sort or call anything
        if (p.Calls[T_A].depth()||p.Calls[T_C].depth()|| p.Calls[T_G].depth()|| p.Calls[T_T].depth()) {
            // silly 15 decimals to match R's default output ... better off with the C default
            static char cse_buf[8192]; 
            #define MEANQ(base,dir) (p.Calls[base].dir?(p.Calls[base].dir##_q/(double)p.Calls[base].dir):0)
           sprintf(cse_buf,"%s\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g%s\n",p.Chr.c_str(), p.Pos, toupper(p.Base)
                    , p.Calls[T_A].fwd, p.Calls[T_C].fwd, p.Calls[T_G].fwd, p.Calls[T_T].fwd
                    , p.Calls[T_A].rev, p.Calls[T_C].rev, p.Calls[T_G].rev, p.Calls[T_T].rev
                    , MEANQ(T_A,fwd), MEANQ(T_C,fwd), MEANQ(T_G,fwd), MEANQ(T_T,fwd)
                    , MEANQ(T_A,rev), MEANQ(T_C,rev), MEANQ(T_G,rev), MEANQ(T_T,rev)
                    , Manager->Reference.data()
                    , p.Calls[T_A].diversity
                    , p.Calls[T_C].diversity
                    , p.Calls[T_G].diversity
                    , p.Calls[T_T].diversity
                    , p.Calls[T_A].agreement
                    , p.Calls[T_C].agreement
                    , p.Calls[T_G].agreement
                    , p.Calls[T_T].agreement
                    , regions
            );
            fputs(cse_buf, cse_f);
            // cse requires separate output for on-target (instead of another column)
            if (tgt_cse_f && p.InTarget) {
                fputs(cse_buf, tgt_cse_f);
            }
        }
    }

	sort(p.Calls.begin(), p.Calls.end(), hitolocall);

	int need_out = -1;
	int skipped_balance=0;
	int skipped_alpha=0;
	int skipped_indel=0;
	int skipped_tail_hom=0;
	int skipped_depth=0;
	int skipped_repeat=0;
	int skipped_diversity=0;
	int skipped_agreement=0;

    vector<vfinal> final_calls;
    for (i=0;i<p.Calls.size();++i) {		// all calls
        //        printf("CALL TOP: depth:%d base: %c, pd: %d, calls: %d\n", (int) p.Calls[i].depth(), p.Calls[i].base, p.Depth, (int) p.Calls.size());

        double pct = (double) p.Calls[i].depth()/p.Depth;
        double qpct = (double) p.Calls[i].qual/p.TotQual;

        if (!p.Calls[i].base)
            continue;

        if (!p.Calls[i].depth())
            continue;

        double bpct = (double) min(p.Calls[i].fwd,p.Calls[i].rev)/p.Calls[i].depth();

        // REBALANCE READS.... CUTTING OFF HIGH COLUMNS
        if (pct >= pct_depth && qpct >= pct_qdepth && (p.Calls[i].depth() >= min_adepth)) {
            if (bpct < pct_balance) {
                int fwd_adj=0, rev_adj=0;
                // f=b*(f+r); r=f/b-f; adj=r-(f/b-f)
                if (p.Calls[i].fwd < p.Calls[i].rev) {
                    rev_adj = (int) p.Calls[i].rev - ( p.Calls[i].fwd/pct_balance  - p.Calls[i].fwd );
                } else {
                    fwd_adj = (int) p.Calls[i].fwd - ( p.Calls[i].rev/pct_balance  - p.Calls[i].rev );
                }
                if (fwd_adj + rev_adj > 1 && bpct > 0) {
                    // adjust call down
                    p.Calls[i].qual -= (rev_adj+fwd_adj)*(p.Calls[i].qual/p.Calls[i].depth()); 
                    p.Calls[i].mq_sum -= (rev_adj+fwd_adj)*(p.Calls[i].mq_sum/p.Calls[i].depth());
                    p.Calls[i].qual_ssq -= (rev_adj+fwd_adj)*(p.Calls[i].qual_ssq/p.Calls[i].depth());
                    p.Calls[i].mq_ssq -= (rev_adj+fwd_adj)*(p.Calls[i].mq_ssq/p.Calls[i].depth());
                    p.Calls[i].rev -= rev_adj;
                    p.Calls[i].fwd -= fwd_adj;
                    skipped_balance+=rev_adj+fwd_adj;

                    // fixed bpct
                    bpct = (double) min(p.Calls[i].fwd,p.Calls[i].rev)/p.Calls[i].depth();
                } else {
                    // it's junk anyway
                }

                // fix depths after adjustment!
                pct = (double) p.Calls[i].depth()/p.Depth;
                qpct = (double) p.Calls[i].qual/p.TotQual;
            }
        }

        if (pct >= pct_depth && qpct >= pct_qdepth && (p.Calls[i].depth() >= min_adepth)) {
            // balance is meaningless at low depths
            if ((bpct >= pct_balance) || (p.Calls[i].depth()<4)) {
                // reads come from diverse positions
                if (p.Calls[i].diversity >= min_diversity) {
                    if (p.Calls[i].agreement >= min_agreement) {
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
                                    if ((repeat_filter == 0) || (p.RepeatCount < repeat_filter)) {
                                        // maybe use rms here... see if it helps
                                        double mean_qual = p.Calls[i].qual/(double)p.Calls[i].depth();
                                        double err_rate = mean_qual < max_phred ? pow(10,-mean_qual/10.0) : global_error_rate;
                                        // expected number of non-reference = error_rate*depth
                                        double pval=(p.Depth*err_rate==0)?0:gsl_ran_poisson_pdf(p.Calls[i].depth(), p.Depth*err_rate);
                                        double padj=total_locii ? pval*total_locii : pval;           // multiple-testing adjustment

                                        if (alpha>=1 || padj <= alpha) {
                                            vfinal final(p.Calls[i]);

                                            double mq_padj=max(total_locii*pow(10,-p.Calls[i].mq_sum/10.0),padj);      // never report pval as better than the total mapping quality
                                            if (debug_xpos) fprintf(stderr,"xpos-debug-pval\tbase:%c, err:%g, pval:%g, padj:%g, mq_padj:%g, mq_sum:%d\n", p.Calls[i].base, err_rate, pval, padj, mq_padj, p.Calls[i].mq_sum);

                                            if (mq_padj > 1) mq_padj=1;

                                            if (need_out == -1) 
                                                need_out = i;

                                            //                                    printf("FINAL: depth:%d base: %s\n", (int) maxc, maxs.c_str());
                                            final.padj=mq_padj;
                                            final.max_idl_cnt=maxc;
                                            final.max_idl_seq=maxs;
                                            final_calls.push_back(final);
                                        } else {
                                            skipped_alpha+=p.Calls[i].depth();
                                        }
                                        // implicitly skip all the ohter indel calls at the same locus
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
                            if (p.Calls[i].base == '*' && (
                                        ((repeat_filter > 0) && (p.RepeatCount >= repeat_filter)) || 
                                        (p.Calls[i].depth() < min_idepth)
                                        )) {
                                skipped_indel+=p.Calls[i].depth();
                            } else {
                                // subtract inserts from reference .. perhaps > 0 is correct here....
                                if (p.Calls[i].is_ref && (ins_rev+ins_fwd) > max(min_idepth,min_adepth)) {
                                    p.Calls[i].fwd-=ins_fwd;
                                    p.Calls[i].rev-=ins_rev;
                                }

                                double mean_qual = p.Calls[i].qual/(double)p.Calls[i].depth();

                                /*
                                   if ( (repeat_filter > 0) && (p.RepeatCount >= repeat_filter) ) {
                                   p.Calls[i].fwd-=p.Calls[i].tail_fwd; 
                                   p.Calls[i].rev-=p.Calls[i].tail_rev;
                                   skipped_tail_hom+=p.Calls[i].tail_fwd+p.Calls[i].tail_rev;
                                   }
                                 */
                                if (p.Calls[i].depth() >= min_adepth && p.Calls[i].depth() > 0) {
                                    double err_rate = mean_qual < vse_max_phred[b2i(p.Base)][b2i(p.Calls[i].base)] ? pow(10,-mean_qual/10.0) : vse_rate[b2i(p.Base)][b2i(p.Calls[i].base)];
                                    // expected number of non-reference bases at this position is error_rate*depth
                                    double pval=(p.Depth*err_rate==0)?0:gsl_ran_poisson_pdf(p.Calls[i].depth(), p.Depth*err_rate);
                                    double padj=total_locii ? pval*total_locii : pval;           // multiple-testing adjustment

                                    if (alpha>=1 || padj <= alpha) {
                                        double mq_padj=max(total_locii*pow(10,-p.Calls[i].mq_sum/10.0),padj);      // never report as better than the mapping quality

                                        if (mq_padj > 1) mq_padj=1;

                                        if (debug_xpos) fprintf(stderr,"xpos-debug-pval\tbase:%c, err:%g, pval:%g, padj:%g, mq_padj:%g, mq_sum:%d\n", p.Calls[i].base, err_rate, pval, padj, mq_padj, p.Calls[i].mq_sum);

                                        if (!p.Calls[i].is_ref || debug_xpos || output_ref) {
                                            if (need_out == -1)
                                                need_out = i;
                                        }
                                        vfinal final(p.Calls[i]);
                                        final.padj=mq_padj;
                                        final_calls.push_back(final);
                                    } else {
                                        skipped_alpha+=p.Calls[i].depth();
                                    }
                                }
                            }
                        }
                    } else {
		                if (debug_xpos) {
                            warn("xpos-skipped-agree-%c\t%g\n", p.Calls[i].base, p.Calls[i].agreement);
                        }
                        skipped_agreement+=p.Calls[i].depth();
                    } 
                } else {
                    skipped_diversity+=p.Calls[i].depth();
                } 
            } else {
                skipped_balance+=p.Calls[i].depth();
            }
        } else {
            // depth is too low now.... technically you can just add all the rest of the calls to skipped_depth without checking
            skipped_depth+=p.Calls[i].depth();
        }
    }

    ++Locii;

	if (need_out>=0||debug_xpos) {

        if (final_calls.size() > 1) {
//            printf("HERE1 %c/%c\n", final_calls[0].pcall->base, final_calls[1].pcall->base);
            if(final_calls[1].pcall->is_ref) {
                vfinal tmp=final_calls[1];
                final_calls[1]=final_calls[0];
                final_calls[0]=tmp;
//                printf("HERE2 %c/%c\n", final_calls[0].pcall->base, final_calls[1].pcall->base);
            }
        }

//        printf("allele_count: %d\n", (int) final_calls.size());

        
        int total_call_depth=0;
        int i;
        for (i=0;i<final_calls.size();++i) {
            total_call_depth+=final_calls[i].pcall->depth();
        }
        double pct_allele = 0;
        if (need_out >=0) {
            // more than 1 call at this position = Het
            if (final_calls.size() > 1) {
                if (final_calls[0].pcall->is_ref) {
                    pct_allele = 100.0 * final_calls[1].pcall->depth() / (double) total_call_depth;
                } else {
                    // no reference seen... but still het?
                    pct_allele = 100.0 * final_calls[0].pcall->depth() / (double) total_call_depth;
                }
                ++Hets;
            } else {
                pct_allele = 100.0 * final_calls[0].pcall->depth() / (double) total_call_depth;
                ++Homs;
            }
        }

        /// INTERNAL VAR FILE
        if (var_f) {
            int i;
            string pil;
            for (i=0;i<final_calls.size();++i) {
               vfinal &f=final_calls[i];
               if (f.is_indel()) {
                    pil += string_format("\t%c%s:%d,%d,%.1e,%.2g,%.2g",f.pcall->base,f.max_idl_seq.c_str(),f.max_idl_cnt,f.pcall->qual/f.pcall->depth(),f.padj, f.pcall->diversity, f.pcall->agreement);
               } else {
                    pil += string_format("\t%c:%d,%d,%.1e,%.2g,%.2g",f.pcall->base,f.pcall->depth(),f.pcall->qual/f.pcall->depth(),f.padj, f.pcall->diversity, f.pcall->agreement);
               }
            }
            fprintf(var_f,"%s\t%d\t%c\t%d\t%d\t%2.2f%s%s\n",p.Chr.c_str(), p.Pos, p.Base, p.Depth, skipped_diversity+skipped_agreement+skipped_alpha+skipped_depth+skipped_balance+p.SkipAmp+p.SkipN+p.SkipDupReads+p.SkipMinMapq+p.SkipMinQual, pct_allele, Manager->UseAnnot==1?(p.InTarget ? "\t1" : "\t0"):"", pil.c_str());

            if (tgt_var_f) {
                if (p.InTarget) {
                    fprintf(tgt_var_f,"%s\t%d\t%c\t%d\t%d\t%2.2f%s\n",p.Chr.c_str(), p.Pos, p.Base, p.Depth, skipped_diversity+skipped_agreement+skipped_alpha+skipped_depth+skipped_balance+p.SkipAmp+p.SkipN+p.SkipDupReads+p.SkipMinMapq+p.SkipMinQual, pct_allele, pil.c_str());
                }
            }
        }

       if (vcf_f) {
            for (i=0;i<final_calls.size();++i) {
               vfinal &f=final_calls[i];
               int qual = f.padj>0?min(40,10*(-log10(f.padj))):40;

               if (f.is_indel()) {
                    string base;
                    string alt;
                    if (f.pcall->base =='-') {
                        base = p.Base + f.max_idl_seq;
                        alt = p.Base;
                    } else {
                        base = p.Base;
                        alt = p.Base + f.max_idl_seq;
                    }
                    double freq_allele = f.max_idl_cnt / (double) p.Depth;
                    fprintf(vcf_f,"%s\t%d\t.\t%s\t%s\t%2d\tPASS\tMQ=%d;BQ=%d;DP=%d;AF=%2.2f\n", 
                        p.Chr.c_str(), p.Pos, base.c_str(), alt.c_str(), qual, 
                        (int) f.pcall->mq_rms(),
                        (int) f.pcall->qual_rms(),
                        total_call_depth,
                        freq_allele);
                } else {
                    char alt = f.pcall->base;
                    if (f.pcall->is_ref) 
                        alt = '.';
                    double freq_allele = f.pcall->depth() / (double) p.Depth;
                    fprintf(vcf_f,"%s\t%d\t.\t%c\t%c\t%d\tPASS\tMQ=%d;BQ=%d;DP=%d;AF=%2.2f\n",
                        p.Chr.c_str(), p.Pos, p.Base, alt, qual,
                        (int) f.pcall->mq_rms(),
                        (int) f.pcall->qual_rms(),
                        total_call_depth,
                        freq_allele);
                }
           }
        }

        if (eav_f) {
//            printf(eav_f,"chr\tpos\tref\tdepth\tnum_states\ttop_consensus\ttop_freq\tvar_base\tvar_depth\tvar_qual\tvar_strands\tforward_strands\treverse_strands\n");
            string top_cons, var_base, var_depth, var_qual, var_strands, forward, reverse, diversity, agreement;
           
            float padj=final_calls[0].padj;
            if (final_calls[0].pcall->is_ref && final_calls.size() > 1) {
                padj=final_calls[1].padj;
            }
            for (i=0;i<final_calls.size();++i) {
                vfinal &f=final_calls[i];
                if (i < 2) {
                    if (i > 0) top_cons += "/";
                    top_cons += f.pcall->base;
                }
                if (i > 0) var_base += "/";
                var_base += f.pcall->base;
                if (f.is_indel()) {
                    if (i < 2) {
                        top_cons += f.max_idl_seq;
                    }
                    var_base += f.max_idl_seq;
                }
                if (i > 0) var_depth+= ";";
                var_depth+= string_format("%d",f.pcall->depth());
                if (i > 0) var_qual+= ";";
                var_qual+= string_format("%d",f.pcall->qual_rms());
                if (i > 0) var_strands+= ";";
                var_strands+= string_format("%d",(f.pcall->fwd>0)+(f.pcall->rev>0));
                if (i > 0) forward += ";";
                forward+= string_format("%d",f.pcall->fwd);
                if (i > 0) reverse += ";";
                reverse+= string_format("%d",f.pcall->rev);
                if (i > 0) agreement += ";";
                agreement+= string_format("%g",f.pcall->agreement);
                if (i > 0) diversity += ";";
                diversity+= string_format("%g",f.pcall->diversity);
            }
            fprintf(eav_f,"%s\t%d\t%c\t%d\t%d\t%s\t%2.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%.1e\t%s\t%s%s\n",p.Chr.c_str(), p.Pos, p.Base, p.Depth, (int) final_calls.size(),top_cons.c_str(), pct_allele, var_base.c_str(), var_depth.c_str(), var_qual.c_str(), var_strands.c_str(), forward.c_str(), reverse.c_str(), padj, diversity.c_str(), agreement.c_str(), Manager->UseAnnot==1?(p.InTarget?"\t1":"\t0"):regions);
        }

		if (debug_xpos) {
		    fprintf(stderr,"xpos-skip-dup\t%d\n",p.SkipDupReads);
		    fprintf(stderr,"xpos-skip-mapq\t%d\n",p.SkipMinMapq);
		    fprintf(stderr,"xpos-skip-qual\t%d\n",p.SkipMinQual);
		    fprintf(stderr,"xpos-skip-bal\t%d\n",skipped_balance);
		    fprintf(stderr,"xpos-skip-depth\t%d\n",skipped_depth);
		    fprintf(stderr,"xpos-skip-indel\t%d\n",skipped_indel);
		    fprintf(stderr,"xpos-skip-repeat\t%d\n",skipped_repeat);
		    fprintf(stderr,"xpos-skip-diversity\t%d\n",skipped_diversity);
		    fprintf(stderr,"xpos-skip-agreement\t%d\n",skipped_agreement);
		    fprintf(stderr,"xpos-skip-alpha\t%d\n",skipped_alpha);
            if (repeat_filter > 0) {
                fprintf(stderr,"repeat-count\t%d\n",p.RepeatCount);
                fprintf(stderr,"repeat-filter\t%d\n",repeat_filter);
                fprintf(stderr,"repeat-base\t%c\n",p.RepeatBase);
            }
			exit(0);
		}
	}
}

void PileupManager::FillReference(int refSize) {
    int flank=(refSize-1)/2;
    Reference.resize(refSize);

    // if you're in the middle of a window
    if (WinDex==flank && Win.size()==refSize) {
        bool needfai=0;
        int i;
        for (i=WinDex-flank;i<refSize;++i) {
            if (!isalpha(Win[i].Base)) {
                needfai=1;
                break;
            } else {
                Reference[i-(WinDex-flank)]=Win[i].Base;
            }
        }
        if (needfai) {
            faidx.Fetch((char *)Reference.data(), Pileup.Chr, Pileup.Pos-flank-1, Pileup.Pos+flank-1);
        }
    } else {
        faidx.Fetch((char *)Reference.data(), Pileup.Chr, Pileup.Pos-flank-1, Pileup.Pos+flank-1);
    }
}

void PileupManager::LoadAnnot(const char *path) {
    FILE *f = fopen(path,"r");
    if (!f) {
        warn("Can't open %s : %s\n", path, strerror(errno));
        exit(1);
    }

    AnnotType = '\0';

    if (!strcmp(fext(path), ".bed"))  
        AnnotType='b';

    if (!strcmp(fext(path), ".gtf"))  
        AnnotType='g';


    if (!AnnotType) { 
        // try to detect it?
        line l; meminit(l);
        int cnt=0;
        while(read_line(f, l)>0) {
            vector<char *> d=split(l.s, '\t');
            if (d.size() >= 7) {
                AnnotType = (*d[5]=='+' || *d[5] == '-') ? 'b' : AnnotType;

                if (!AnnotType)  
                    AnnotType = (*d[6]=='+' || *d[6] == '-') ? 'g' : AnnotType;
            }
            break;
        }
        warn("detect-annot\t%s\n", AnnotType == 'g' ? "gtf" : AnnotType == 'b' ? "bed" : "unknown");
    }

    fclose(f);

    if (!AnnotDex.read(path)) {
        // try building if we know the type
        if (AnnotType) {
        //    void build(const char *path, const char *sep, int nchr, int nbeg, int nend, int skip_i, char skip_c, int sub_e);
            AnnotDex.build(path, "\t", 0, 1, 2, 0, '#', AnnotType=='b' ? 1 : 0);
        }
        die("Either %s.tidx must be a valid tidx indexed file, or %s must be a BED or GTF file\n", path, path);
    }

    UseAnnot=1;
}

void VarStatVisitor::Visit(PileupSummary &p) {
	tot_locii += 1;

//    if (tot_locii % 10000 == 0) {
//        fprintf(stderr,"\r%.4f          ", 100*(float)tot_locii/777385781.0);
//    }

	if (p.Depth < minsampdepth)
		return;

    // insert and deletions have their own, separate noise levels

	int ins_depth = p.Calls.size() > 6 ? p.Calls[6].depth() : 0;
	int ins_qual = p.Calls.size() > 6 ? p.Calls[6].qual : 0;
	double ins_noise = 0;
	double ins_qnoise = 0;
	if (p.Calls.size() > 1 && p.Depth > 0 && ins_depth > 0) {
		ins_noise = (double) ins_depth/p.Depth;
		ins_qnoise = (double) ins_qual/p.TotQual;
	}

	int del_depth = p.Calls.size() > 5 ? p.Calls[5].depth() : 0;
	int del_qual = p.Calls.size() > 5 ? p.Calls[5].qual : 0;
	double del_noise = 0;
	double del_qnoise = 0;
	if (p.Calls.size() > 1 && p.Depth > 0 && del_depth > 0) {
		del_noise = (double) del_depth/p.Depth;
		del_qnoise = (double) del_qual/p.TotQual;
	}

    // snp's are "noise" ... since they are supposedly rare.
	int i;
	if (p.Calls.size() > 5) 
		p.Calls.resize(5);		// toss N's and inserts before sort

	sort(p.Calls.begin(), p.Calls.end(), hitolocall);

	double noise;
	double qnoise;
    if (p.Calls.size() > 1) {
        // assume non-reference is noise
        noise = (double) p.Calls[1].depth()/p.Depth;
        qnoise = (double) p.Calls[1].qual/p.TotQual;
        if (noise > .25) {
            // unless maybe that was a het or something....
            if (p.Calls.size() > 2) {
                // but 3rd allele is always noise
                noise = (double) p.Calls[2].depth()/p.Depth;
                qnoise = (double) p.Calls[2].qual/p.TotQual;
            } else {
                // this is weird... but ok
                noise = 0;
                qnoise = 0;
            }
        }
    } else {
        noise = qnoise = 0.0;
    }

	double mnqual = (double)p.TotQual/p.Depth;

	char pbase = p.Calls.size() > 1 ? p.Calls[1].base : '.';

	if (noise_f) {
		fprintf(noise_f,"%d\t%c\t%c\t%f\t%f\t%f\n", p.Depth, p.Base, pbase, noise, qnoise, mnqual);
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
	stats.push_back(Noise(p.Base, p.Calls[1].base, p.Depth, noise, qnoise, mnqual));
	ins_stats.push_back(Noise(p.Base, p.Calls[1].base, p.Depth, ins_noise, ins_qnoise, mnqual));
	del_stats.push_back(Noise(p.Base, p.Calls[1].base, p.Depth, del_noise, del_qnoise, mnqual));
}


void usage(FILE *f) {
        fprintf(f,
"Usage: varcall <-s|-v> <-f REF> [options] bam1 [bam2...]\n"
"Version: %s (BETA)\n"
"\n"
"Either outputs summry stats for the list of files, or performs variant calling\n"
"\n"
"Options (later options override earlier):\n"
"\n"
"-s          Calculate statistics\n"
"-v|version  Calculate variants bases on supplied parameters (see -S)\n"
"-f          Reference fasta (required if using bams, ignored otherwise)\n"
"-m          Min locii depth (1)\n"
"-a          Min allele depth (2)\n"
"-p          Min allele pct by quality (0)\n"
"-q          Min qual (3)\n"
"-Q          Min mapping quality (0)\n"
"-b          Min pct balance (strand/total) (0)\n"
"-D FLOAT    Max duplicate read fraction (depth/length per position) (1)\n"
"-d FLOAT    Minimum diversity (CV from optimal depth) (0.25)\n"
"-G FLOAT    Minimum agreement (Weighted CV of positional variation) (0.25)\n"
"-0          Zero out all filters, set e-value filter to 1, report everything\n"
"-B          If running from a BAM, turn off BAQ correction (false)\n"
"-R          Homopolymer repeat indel filtering (8)\n"
"-e FLOAT    Alpha filter to use, requires -l or -S (.05)\n"
"-g FLOAT    Global minimum error rate (default: assume phred is ok)\n"
"-l INT      Number of locii in total pileup used for bonferroni (1 mil)\n"
"-x CHR:POS  Output this pos only, then quit\n"
"-S FILE     Read in statistics and params from a previous run with -s (do this!)\n"
"-A ANNOT    Calculate in-target stats using the annotation file (requires -o)\n"
"-o PREFIX   Output prefix (works with -s or -v)\n"
"-F files    List of file types to output (var, varsum, eav, vcf)\n"
"\n"
"Extended Options\n"
"\n"
"--pcr-annot   BED      Only include reads adhering to the expected amplicons\n"
"--stranded    TYPE     Can be FR (the default), FF, FR.  Used with pcr-annot\n"
"--diversity|d FLOAT    Alias for -d\n"
"--agreement|G FLOAT    Alias for -G\n"
"--no-indels            Ignore all indels\n"
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
"   PREFIX.var         Variant calls in tab delimited 'varcall' format\n"
"   PREFIX.eav         Variant calls in tab delimited 'ea-var' format\n"
"   PREFIX.cse         Variant calls in tab delimited 'varprowl' format\n"
"   PREFIX.vcf         Variant calls, in vcf format\n"
"   PREFIX.varsum      Summary of variant calls\n"
"   PREFIX.tgt.var     On-target version of .var\n"
"   PREFIX.tgt.cse     On-target version of .cse\n"
"   PREFIX.tgt.varsum  On-target version of .varsum\n"
"\n"
"Stats Output:\n"
"\n"
"Contains mean, median, quartile information for depth, base quality, read len,\n"
"mapping quality, indel levels. Also estimates parameters suitable for\n"
"variant calls, and can be passed directly to this program for variant calls\n"
"\n"
"If an output prefix is used, files are created as follows:\n"
"\n"
"   PREFIX.stats       Stats output\n"
"   PREFIX.noise       Non-reference, non-homozygous allele summary\n"
"   PREFIX.xnoise      Like noise, but with context-specific rates\n"
"\n"
"Filtering Details:\n"
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
	static int i;
    static char c;
	for ( i=0;i<str.size();++i ) {
        c=((char *)(void *)str.data())[i];
		if (c >= 'a' && c <= 'z') {
            ((char *)(void *)str.data())[i]=c-32;
        }
	}
}

// returns quantile depth 
double quantile_depth(const std::vector<Noise> &vec, double p) {
        int l = vec.size();
        assert(l > 0);
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
    char *sav;
    char* token = strtok_r(str,delim, &sav);
    std::vector<char *> result;
    while(token != NULL)
    {
        result.push_back(token);
        token = strtok_r(NULL,delim, &sav);
    }
    return result;
}

std::vector<char *> split(char* str, char delim)
{
    char*p=strchr(str,delim);
    std::vector<char *> result;
    while(p != NULL)
    {
        *p='\0';
        result.push_back(str);
        str=p+1;
        p=strchr(str,delim);
    }
    if (*str) 
        result.push_back(str);
    return result;
}

int split(char **buf, char* str, char delim)
{
    char **bb=buf;
    char *p=strchr(str,delim);
    std::vector<char *> result;
    while(p != NULL)
    {
        *p='\0';
        *bb=str;
        ++bb;
        str=p+1;
        p=strchr(str,delim);
    }
    if (*str) {
        *bb=str;
        bb++;
    }
    *bb=NULL;
    return bb-buf;
}

int rand_round(double x) {
    return floor(x)+((rand()>(x-int(x))) ? 1 : 0);
//warn("rr:%f=%d\n",x);
}


void check_ref_fai(const char * ref) {
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
}

void Faidx::Load(const char *path) {
    fa_f = openordie(path, "r");
    FILE *fp = openordie(string_format("%s.fai", path).c_str(), "r");

    char *buf = (char*)calloc(0x10000, 1);
    char *p;
    while (!feof(fp) && fgets(buf, 0x10000, fp)) {
        for (p = buf; *p && isgraph(*p); ++p);
        *p = 0; ++p;
        Faient e;
        sscanf(p, "%d%lld%d%d", &e.len, &e.offset, &e.line_blen, &e.line_len);
        string chr = buf;
        faimap[chr]=e;
    }
    free(buf);
}

bool Faidx::Fetch(char *buf, const Faient *ent, int pos_from, int pos_to) {
    int len = (pos_to-pos_from+1);

    if (fseek(fa_f, ent->offset + pos_from / ent->line_blen * ent->line_len + pos_from % ent->line_blen, SEEK_SET) == -1) {
        return false;
    }
    int l = 0;
    char c;
    while (((c=fgetc(fa_f))!= EOF) && (l < len)) {
        if (isgraph(c)) buf[l++] = c;
    }
    return l==len;
}


PileupSubscriber::PileupSubscriber(PileupManager &man) {
    Manager = NULL; 
    SetManager(man);
};

void PileupSubscriber::SetManager(PileupManager &man) { 
    assert(!Manager);
    Manager = &man; 
    man.Kids.push_back(this); 
};



void output_stats(VarStatVisitor &vstat) {
    stat_out("version\tvarcall-%s\n", VERSION);
    stat_out("min depth\t%d\n", minsampdepth);
    stat_out("alpha\t%f\n", alpha);

    if (vstat.stats.size()) {
        // sort by depth descending
        sort(vstat.stats.begin(), vstat.stats.end(), noisebydepth);

        // flip 3 and 1 because sorted in descending order for sampling (above)
        double depth_q3=quantile_depth(vstat.stats, .25);
        double depth_q2=quantile_depth(vstat.stats, .50);
        double depth_q1=quantile_depth(vstat.stats, .75);
        double depth_qx=quantile_depth(vstat.stats, .95);

        // number of locii to compute error rate
        int ncnt=min(100000,vstat.stats.size());

        int i;
        double nsum=0, nssq=0, dsum=0, dmin=vstat.stats[0].depth, qnsum=0, qnssq=0, qualsum=0;

        double ins_nsum=0, ins_nssq=0, del_nsum=0, del_nssq=0;

        double qvsum[T_CNT][T_CNT], qvssq[T_CNT][T_CNT]; int qvcnt[T_CNT][T_CNT];
        meminit(qvsum);
        meminit(qvssq);
        meminit(qvcnt);

        for (i=0;i<ncnt;++i) {
            if (vstat.stats[i].depth < depth_q1) {
                continue;
            }

            int ref_i, var_i;
            ref_i=b2i(vstat.stats[i].ref);
            var_i=b2i(vstat.stats[i].var);

            if (ref_i < T_N && var_i < T_N) {
                qvsum[ref_i][var_i]+=vstat.stats[i].qnoise;
                qvssq[ref_i][var_i]+=vstat.stats[i].qnoise*vstat.stats[i].qnoise;
                qvcnt[ref_i][var_i]+=1;
            }

            nsum+=vstat.stats[i].noise;
            nssq+=vstat.stats[i].noise*vstat.stats[i].noise;
            dsum+=vstat.stats[i].depth;
            qnsum+=vstat.stats[i].qnoise;
            qnssq+=vstat.stats[i].qnoise*vstat.stats[i].qnoise;
            qualsum+=vstat.stats[i].mnqual;
            if (vstat.stats[i].depth < dmin) dmin = vstat.stats[i].depth;
            ins_nsum+=vstat.ins_stats[i].noise;
            ins_nssq+=vstat.ins_stats[i].noise*vstat.ins_stats[i].noise;
            del_nsum+=vstat.del_stats[i].noise;
            del_nssq+=vstat.del_stats[i].noise*vstat.del_stats[i].noise;
        }

        double noise_mean =nsum/ncnt;
        double noise_dev = stdev(ncnt, nsum, nssq);
        double qnoise_mean =qnsum/ncnt;
        double qnoise_dev = stdev(ncnt, qnsum, qnssq);
        double qual_mean = qualsum/ncnt;
        double ins_noise_mean =ins_nsum/ncnt;
        double ins_noise_dev = stdev(ncnt, ins_nsum, ins_nssq);
        double del_noise_mean =del_nsum/ncnt;
        double del_noise_dev = stdev(ncnt, del_nsum, del_nssq);

        stat_out("qual mean\t%.4f\n", qual_mean);
        stat_out("noise mean\t%.6f\n", noise_mean);
        stat_out("noise dev\t%.6f\n", noise_dev);
        stat_out("qnoise mean\t%.6f\n", qnoise_mean);
        stat_out("qnoise dev\t%.6f\n", qnoise_dev);
        stat_out("ins freq\t%.6f\n", ins_noise_mean);
        stat_out("ins freq dev\t%.6f\n", ins_noise_dev);
        stat_out("del freq\t%.6f\n", del_noise_mean);
        stat_out("del freq dev\t%.6f\n", del_noise_dev);

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

        int locii_gtmin=0;
        for (i=0;i<vstat.stats.size();++i) {
            if (vstat.stats[i].depth >= min_depth) {
                ++locii_gtmin;
            }
        }
        stat_out("locii >= min depth\t%d\n", locii_gtmin);
        stat_out("locii\t%d\n", vstat.tot_locii);

        double stdevfrommean=-qnorm((alpha/locii_gtmin)/2);
        stat_out("qnorm adj\t%f\n", stdevfrommean);


        pct_qdepth=qnoise_mean+qnoise_dev*stdevfrommean;
        stat_out("min pct qual\t%.4f\n", 100*pct_qdepth);

        meminit(vse_rate);

        // variation-specific error rates
        int j;
        for (i=0;i<T_CNT;++i) {
        for (j=0;j<T_CNT;++j) {
            if (i!=j && (qvcnt[i][j] > 20)) {
                double qn_mean =qvsum[i][j]/qvcnt[i][j];
                double qn_dev = stdev(qvcnt[i][j], qvsum[i][j], qvssq[i][j]);
                stat_out("vnoise mean %c:%c\t%.6f\n", i2b(i), i2b(j), qn_mean);
                stat_out("vnoise dev %c:%c\t%.6f\n", i2b(i), i2b(j), qn_dev);
                vse_rate[i][j]=qn_mean+qn_dev;          // vse rate
            } else {
                vse_rate[i][j]=noise_mean+noise_dev;              // global rate
            }
        }
        }

        
        // now set params... as if you just read them in
        // this should mirror "read stats"
        min_depth = minsampdepth;
        global_error_rate = noise_mean+noise_dev;
        total_locii = locii_gtmin;
    }
}


