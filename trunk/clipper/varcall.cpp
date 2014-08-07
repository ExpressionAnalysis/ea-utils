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

#define SVNREV atoi(strchr("$Revision$", ':')+1)
const char * VERSION = "0.9";

#define MIN_READ_LEN 20
#define DEFAULT_LOCII 1000000


using namespace std;
using namespace google;

void usage(FILE *f);

// #define DEBUG 1

#define meminit(l) (memset(&l,0,sizeof(l)))
#ifdef DEBUG
    #define debug(s,...) fprintf(stderr,s,##__VA_ARGS__)
#else
    #define debug(s,...)
#endif
#undef warn
#define warn(s,...) ++errs; fprintf(stderr,s,##__VA_ARGS__)
#define die(s,...) (fprintf(stderr,s,##__VA_ARGS__), exit(1))
#define stat_out(s,...) fprintf(stat_fout,s,##__VA_ARGS__)
#define stdev(cnt, sum, ssq) sqrt((((double)cnt)*ssq-pow((double)sum,2)) / ((double)cnt*((double)cnt-1)))
#define log10(x) (log(x)/log(10))

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
std::string string_format(const std::string &fmt, ...);
void to_upper(const std::string str);
void rename_tmp(std::string f);

int errs=0;
extern int optind;
int g_lineno=0;

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
    vcall() {base='\0'; mn_qual=mq0=fwd=rev=qual=is_ref=qual_ssq=mq_sum=mq_ssq=tail_rev=tail_fwd=fwd_q=rev_q=0;}
    char base;
	bool is_ref;
    int qual, fwd, rev, mq0, mn_qual, qual_ssq, mq_sum, mq_ssq, tail_rev, tail_fwd, fwd_q, rev_q;
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
    bool InTarget;

    int SkipN;
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

        bool UseAnnot;
        tidx AnnotDex;          // start/stop index file
        char AnnotType;         // b (bed) or g (gtf - preferred)
        
        PileupReads Reads;
        PileupVisitor() {InputType ='\0';}
		void Parse(char *dat) {PileupSummary p(dat, Reads); Visit(p);};
        void LoadAnnot(const char *annot_file);
		virtual void Visit(PileupSummary &dat)=0;
		virtual void Finish()=0;
};

class VarStatVisitor : public PileupVisitor {
    public:
    VarStatVisitor() : PileupVisitor() {tot_locii=0; tot_depth=0; num_reads=0;};

    void Visit(PileupSummary &dat);
    void Finish() {};

    public:
	double tot_depth;
	int tot_locii;
	int num_reads;
	vector<Noise> stats;
	vector<Noise> ins_stats;
	vector<Noise> del_stats;
};

class VarCallVisitor : public PileupVisitor {

    deque<PileupSummary> Win;
    void VisitX(PileupSummary &dat, int windex);

    public:
    int WinMax;
    VarCallVisitor() : PileupVisitor() {SkippedDepth=0;WinMax=0;Hets=0;Homs=0;Locii=0;};

    void Visit(PileupSummary &dat);
    void Finish();

	int SkippedDepth;
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
int total_locii=-1;
double pct_balance=0;
char *debug_xchr=NULL;
int debug_xpos=0;
int min_depth=1;
int min_mapq=0;
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

void parse_bams(PileupVisitor &v, int in_n, char **in, const char *ref);
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


Faidx faidx;

int main(int argc, char **argv) {
	char c;
	const char *noiseout=NULL;
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
    char *read_stats = NULL;


#define MAX_F 20
    const char *format_list[MAX_F]={"var", "eav", "noise", "varsum", NULL};

	while ( (c = getopt_long(argc, argv, "?svVBhe:m:N:x:f:p:a:g:q:Q:i:o:D:R:b:L:S:F:A:",NULL,NULL)) != -1) {
		switch (c) {
			case 'h': usage(stdout); return 0;
			case 'm': umindepth=ok_atoi(optarg); break;
			case 'q': min_qual=ok_atoi(optarg); break;
			case 'o': out_prefix=optarg; break;
			case 'Q': min_mapq=ok_atoi(optarg); break;
			case 'V': printf("Version: %s.%d\n", VERSION, SVNREV); exit(0); break;
			case 'R': repeat_filter=ok_atoi(optarg); break;
			case 'A': target_annot=optarg; break;
			case 'a': uminadepth=ok_atoi(optarg);break;
			case 'D': artifact_filter=atof(optarg);break;
			case 'i': uminidepth=ok_atoi(optarg);break;
			case 'x': {
					debug_xchr=optarg;
					char *p=strrchr(debug_xchr, ':');
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
			case 'g': global_error_rate=atof(optarg); break;
			case 'L': total_locii=ok_atoi(optarg); break;
			case 'f': ref=optarg; break;
			case 'N': noiseout=optarg; break;
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


	if (!do_stats && !do_varcall || do_stats && do_varcall) {
		warn("Specify -s for stats only, or -v to do variant calling\n\n");
		usage(stderr);
		return 1;
	}

    if (out_prefix) {
        if (!do_varcall) {
            warn("Specify -o with -v only\n\n");
            usage(stderr);
            return 1;
        }

        var_f = openordie(string_format("%s.var.tmp", out_prefix).c_str(), "w");

        fprintf(var_f,"%s\t%s\t%s\t%s\t%s\t%s\t%s%s\n","chr", "pos", "ref", "depth", "skip", "pct", target_annot ? "target\t" : "", "...");

        varsum_f = openordie(string_format("%s.varsum.tmp", out_prefix).c_str(), "w");

        if (target_annot) {
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
            if (target_annot) 
                tgt_cse_f = openordie(string_format("%s.tgt.cse.tmp", out_prefix).c_str(), "w");
        }
    } else {
        var_f = stdout;
        varsum_f = stderr;
    } 

	if (umindepth > minsampdepth) {
		minsampdepth=umindepth;
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

    // not really random
    srand(1);

    max_phred = -log10(global_error_rate)*10;

	if (do_stats) {
		FILE *stat_fout=stdout;			// stats to stdout

		if (do_varcall) 				// unless varcalling at the same time
			stat_fout=stderr;

		VarStatVisitor vstat;

		parse_bams(vstat, in_n, in, ref);

		stat_out("version\tvarcall-%s.%d\n", VERSION, SVNREV);
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
            for (i=0;i<ncnt;++i) {
                if (vstat.stats[i].depth < depth_q1) {
                    continue;
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
        }
	}

    if (read_stats){
        FILE * f = fopen(read_stats, "r");
        if (!f) {
            warn("File %s does not exist, quitting\n", read_stats);
            exit(1);
        }
        line l; meminit(l);
        char *val;
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
                    if (global_error_rate<=0) global_error_rate=atof(val); 
                } else if (!strcasecmp(l.s, "locii >= min depth")) {
                    if (total_locii<0) total_locii=atoi(val); 
                } else if (!strcasecmp(l.s, "alpha")) {
                    if (alpha<=0) alpha=atof(val); 
                }
            }
        }
    }

    if (total_locii<0) total_locii=DEFAULT_LOCII;
    if (total_locii==0) total_locii=1;          // no adjustment

    if (eav_f) {
        fprintf(eav_f,"chr\tpos\tref\tdepth\tnum_states\ttop_consensus\ttop_freq\tvar_base\tvar_depth\tvar_qual\tvar_strands\tforward_strands\treverse_strands\t%cval\n",total_locii>1?'e':'p');
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

		fprintf(varsum_f,"version\tvarcall-%s.%d\n", VERSION, SVNREV);
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

		VarCallVisitor vcall;

        if (target_annot) {
            vcall.LoadAnnot(target_annot);
       }

        if (cse_f) {
            vcall.WinMax=21;
        } else if (repeat_filter > 0) {
		    fprintf(varsum_f,"homopolymer filter\t%d\n", repeat_filter);
            vcall.WinMax=repeat_filter+repeat_filter+3;
        } else {
            vcall.WinMax=5;
        }

        if (vcf_f) {
            // print VCF header
            fprintf(vcf_f, "%s\n", "##fileformat=VCFv4.1");
        }
        if (cse_f) {
            fprintf(cse_f, "Chr\tPos\tRef\tA\tC\tG\tT\ta\tc\tg\tt\tAq\tCq\tGq\tTq\taq\tcq\tgq\ttq\tRefAllele\n");
        }

		parse_bams(vcall, in_n, in, ref);

        if (vcall.InputType == 'B') {
        	fprintf(varsum_f,"baq correct\t%s\n", (no_baq?"no":"yes"));
        }
        fprintf(varsum_f,"locii\t%d\n", vcall.Locii);
        fprintf(varsum_f,"hom calls\t%d\n", vcall.Homs);
        fprintf(varsum_f,"het calls\t%d\n", vcall.Hets);
        fprintf(varsum_f,"locii below depth\t%d\n", vcall.SkippedDepth);

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
            if (noise_f) rename_tmp(string_format("%s.noise.tmp", out_prefix));

            if (tgt_var_f) rename_tmp(string_format("%s.tgt.var.tmp", out_prefix));
            if (tgt_cse_f) rename_tmp(string_format("%s.tgt.cse.tmp", out_prefix));
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

void parse_bams(PileupVisitor &v, int in_n, char **in, const char *ref) {

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

#define T_A 0
#define T_C 1
#define T_G 2
#define T_T 3
#define T_SDEL 4
#define T_NDEL 5
#define T_INS 6
#define T_N 7
#define b2i(c) ((c)=='A'?0:(c)=='a'?0:(c)=='C'?1:(c)=='c'?1:(c)=='G'?2:(c)=='g'?2:(c)=='T'?3:(c)=='t'?3:(c)=='*'?4:(c)=='-'?5:(c)=='+'?6:7)
#define i2b(i) (i==0?'A':i==1?'C':i==2?'G':i==3?'T':i==4?'*':i==5?'-':i==6?'+':'?')

bool hitoloint (int i,int j) { return (i>j);}

int track_readlen[10000];


PileupSummary::PileupSummary(char *line, PileupReads &rds) {

	vector<char *> d=split(line, '\t');

	if (d.size() < 6) {
		warn("Can't read pileup : %d fields, need 6 columns, line %d\n", (int) d.size(), g_lineno);
		exit(1);
	}

	const char * p_qual=d[5];

	Chr=d[0];
	Pos=atoi(d[1]);
	Base=*(d[2]);
	Depth = atoi(d[3]);
	SkipDupReads = 0;
	SkipN = 0;
	SkipMinQual = 0;
	SkipMinMapq = 0;
	MaxDepthByPos = 0;
	RepeatCount = 0;
	RepeatBase = '\0';
	NumReads = 0;
    InTarget = 0;

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

		if (o == '>' || o == '<') {	
			c = 'N';					// no call
			is_ref = 1;
		}

		bool skip = 0;

        // probably should not be adding anything here... but the old code added 1 and floored... new code adds .5 and rounds... which is comparable
        // really.. should just be adding zero, the reason the old code had it was because of a lack of max()
		if (c == 'N') {
			++SkipN;
			skip=1;
		} else if (artifact_filter > 0 && (depthbypos[pia] > max(1,rand_round(0.5+artifact_filter * (Depth/rds.MeanReadLen()))))) {
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

void VarCallVisitor::Finish() {
    // finish out the rest of the pileup, with the existing window
    int vx = WinMax/2+1;
    while (vx < Win.size()) {
        ///debug("Finish: %d\n", Win[vx].Pos);
        VisitX(Win[vx++], vx);
    }
}

void VarCallVisitor::VisitX(PileupSummary &p, int windex) {
    //debug("VisitX: %d\n", p.Pos);

	if (debug_xpos) {
		if (p.Pos != debug_xpos)
			return;
		if (strcmp(debug_xchr,p.Chr.data())) 
			return;
	}

	if (p.Depth < min_depth) {
        if (debug_xpos) {
            fprintf(stderr,"xpos-skip-depth\t%d < %d\n",p.Depth, min_depth);
		    fprintf(stderr,"xpos-skip-dup\t%d\n",p.SkipDupReads);
		    fprintf(stderr,"xpos-skip-n\t%d\n",p.SkipN);
		    fprintf(stderr,"xpos-skip-mapq\t%d\n",p.SkipMinMapq);
		    fprintf(stderr,"xpos-skip-qual\t%d\n",p.SkipMinQual);
        }
		++SkippedDepth;
		return;
	}

    if (UseAnnot) {
        // index lookup only.... not string lookup
        const std::vector <long int> * v = &(AnnotDex.lookup(p.Chr.data(), p.Pos + (AnnotType=='b' ? -1 : 0)));
        if (v && v->size()) {
            p.InTarget=1;
        }
    }

	int ins_fwd = p.Calls.size() > 6 ? p.Calls[6].fwd : 0;
	int ins_rev = p.Calls.size() > 6 ? p.Calls[6].rev : 0;

	int i;
	if (p.Calls.size() > 6) 
		p.Calls.resize(7);	// toss N's before sort

    // OUTPUT CSE BEFORE REORDRED BASES!
    if (cse_f) {
        if (p.Calls.size() < 4) 
            p.Calls.resize(4);	// cse needs 4 calls
    
        char ref21[22];
        if (windex==10 && Win.size()==21) {
            bool needfai=0;
            for (i=windex-10;i<21;++i) {
                if (!isalpha(Win[i].Base)) {
                    needfai=1;
                    break;
                } else {
                    ref21[i-(windex-10)]=Win[i].Base;
                }
            }
            if (needfai) {
                faidx.Fetch(ref21, p.Chr, p.Pos-10, p.Pos+10);
            }
        }
        ref21[21]='\0';

        // cse format... no need to sort or call anything
        if (p.Calls[T_A].depth()||p.Calls[T_C].depth()|| p.Calls[T_G].depth()|| p.Calls[T_T].depth()) {
            // silly 15 decimals to match R's default output ... better off with the C default
            static char cse_buf[8192]; 
            sprintf(cse_buf,"%s\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%s\n",p.Chr.c_str(), p.Pos, toupper(p.Base)
                    , p.Calls[T_A].fwd, p.Calls[T_C].fwd, p.Calls[T_G].fwd, p.Calls[T_T].fwd
                    , p.Calls[T_A].rev, p.Calls[T_C].rev, p.Calls[T_G].rev, p.Calls[T_T].rev
                    , p.Calls[T_A].fwd_q/(double)p.Calls[T_A].fwd, p.Calls[T_C].fwd_q/(double)p.Calls[T_C].fwd, p.Calls[T_G].fwd_q/(double)p.Calls[T_G].fwd, p.Calls[T_T].fwd_q/(double)p.Calls[T_T].fwd
                    , p.Calls[T_A].rev_q/(double)p.Calls[T_A].rev, p.Calls[T_C].rev_q/(double)p.Calls[T_C].rev, p.Calls[T_G].rev_q/(double)p.Calls[T_G].rev, p.Calls[T_T].rev_q/(double)p.Calls[T_T].rev
                    , ref21
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
                            double err_rate = mean_qual < max_phred ? pow(10,-mean_qual/10.0) : global_error_rate;
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
                    pil += string_format("\t%c%s:%d,%d,%.1e",f.pcall->base,f.max_idl_seq.c_str(),f.max_idl_cnt,f.pcall->qual/f.pcall->depth(),f.padj);
               } else {
                    pil += string_format("\t%c:%d,%d,%.1e",f.pcall->base,f.pcall->depth(),f.pcall->qual/f.pcall->depth(),f.padj);
               }
            }
            fprintf(var_f,"%s\t%d\t%c\t%d\t%d\t%2.2f%s%s\n",p.Chr.c_str(), p.Pos, p.Base, p.Depth, skipped_alpha+skipped_depth+skipped_balance+p.SkipN+p.SkipDupReads+p.SkipMinMapq+p.SkipMinQual, pct_allele, UseAnnot?(p.InTarget ? "\t1" : "\t0"):"", pil.c_str());

            if (tgt_var_f) {
                if (p.InTarget) {
                    fprintf(tgt_var_f,"%s\t%d\t%c\t%d\t%d\t%2.2f%s\n",p.Chr.c_str(), p.Pos, p.Base, p.Depth, skipped_alpha+skipped_depth+skipped_balance+p.SkipN+p.SkipDupReads+p.SkipMinMapq+p.SkipMinQual, pct_allele, pil.c_str());
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
            string top_cons, var_base, var_depth, var_qual, var_strands, forward, reverse;
           
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
            }
            fprintf(eav_f,"%s\t%d\t%c\t%d\t%d\t%s\t%2.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%.1e\n",p.Chr.c_str(), p.Pos, p.Base, p.Depth, (int) final_calls.size(),top_cons.c_str(), pct_allele, var_base.c_str(), var_depth.c_str(), var_qual.c_str(), var_strands.c_str(), forward.c_str(), reverse.c_str(), padj);
        }

		if (debug_xpos) {
		    fprintf(stderr,"xpos-skip-dup\t%d\n",p.SkipDupReads);
		    fprintf(stderr,"xpos-skip-mapq\t%d\n",p.SkipMinMapq);
		    fprintf(stderr,"xpos-skip-qual\t%d\n",p.SkipMinQual);
		    fprintf(stderr,"xpos-skip-bal\t%d\n",skipped_balance);
		    fprintf(stderr,"xpos-skip-depth\t%d\n",skipped_depth);
		    fprintf(stderr,"xpos-skip-indel\t%d\n",skipped_indel);
//		    fprintf(stderr,"xpos-skip-tail-imbalance\t%d\n",skipped_tail_hom);
		    fprintf(stderr,"xpos-skip-repeat\t%d\n",skipped_repeat);
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


void PileupVisitor::LoadAnnot(const char *path) {
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
                AnnotType = (*d[5]=='+' || *d[5] == '-') ? 'b' : '\0';  
                AnnotType = (*d[6]=='+' || *d[6] == '-') ? 'g' : AnnotType;
            }
            break;
        }
        if (AnnotType)
            warn("detect annot\t%c\n", AnnotType);
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
	ins_stats.push_back(Noise(p.Depth, ins_noise, ins_qnoise, mnqual));
	del_stats.push_back(Noise(p.Depth, del_noise, del_qnoise, mnqual));
}


void usage(FILE *f) {
        fprintf(f,
"Usage: varcall <-s|-v> <-f REF> [options] bam1 [bam2...]\n"
"Version: %s.%d (BETA)\n"
"\n"
"Either outputs summry stats for the list of files, or performs variant calling\n"
"\n"
"Options (later options override earlier):\n"
"\n"
"-s          Calculate statistics\n"
"-v          Calculate variants bases on supplied parameters (see -S)\n"
"-f          Reference fasta (required if using bams, ignored otherwise)\n"
"-m          Min locii depth (0)\n"
"-a          Min allele depth (0)\n"
"-p          Min allele pct by quality (0)\n"
"-q          Min qual (3)\n"
"-Q          Min mapping quality (0)\n"
"-b          Min pct balance (strand/total) (0)\n"
"-D FLOAT    Max duplicate read fraction (depth/length per position) (1)\n"
"-B          Turn off BAQ correction (false)\n"
"-R          Homopolymer repeat indel filtering (8)\n"
"-e FLOAT    Alpha filter to use, requires -l or -S (.05)\n"
"-g FLOAT    Global minimum error rate (default: assume phred is ok)\n"
"-l INT      Number of locii in total pileup used for bonferroni (1 mil)\n"
"-x CHR:POS  Output this pos only, then quit\n"
"-N FIL      Output noise stats to FIL\n"
"-S FIL      Read in statistics and params from a previous run with -s (do this!)\n"
"-A ANNOT    Calculate in-target stats using the annotation file (requires -o)\n"
"-o PREFIX   Output prefix (note: overlaps with -N)\n"
"-F files    List of file types to output (var, varsum, eav, vcf)\n"
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
"   PREFIX.var       Variant calls in tab delimited 'varcall' format\n"
"   PREFIX.eav       Variant calls in tab delimited 'ea-var' format\n"
"   PREFIX.cse       Variant calls in tab delimited 'varprowl' format\n"
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
"Filtering Details:\n"
"\n"
        ,VERSION, SVNREV);
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
    char* token = strtok(str,delim);
    std::vector<char *> result;
    while(token != NULL)
    {
        result.push_back(token);
        token = strtok(NULL,delim);
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

