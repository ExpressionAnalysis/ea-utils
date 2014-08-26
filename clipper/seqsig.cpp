
#include <iostream>
#include <sparsehash/sparse_hash_map>
#include "fastq-lib.h"

#include <time.h>

#include <string>
#include <vector>

#include <stdio.h>
#include <getopt.h>
#include <stdarg.h>

#define meminit(l) (memset(&l,0,sizeof(l)))
#define VERSION "0.5"

FILE *zopenordie(const char *path, const char *mode);
FILE *openordie(const char *path, const char *mode);
FILE *popenordie(const char *path, const char *mode);
std::string string_format(const std::string &fmt, ...);

using namespace std;
void seq2char(uint32_t s, int k, char *p);
// true if different
void compare_files(const char *f1, const char *f2);
template <class vtype> double quantile(const vtype &vec, double p);

#define AS_UINT(x) ((int *)(void *)(x))

// main

unsigned long long basemap[256];
unsigned long long basemaprc[256];
bool lowcomplex(char *s, int n, int max);
void usage(FILE *f);

// coefficients used to produce the blended similarity score (range 0 to 1, where 1 is identical)

double scoef=.75;
double qcoef=350;
int main(int argc, char **argv) {
    int i;
    clock_t clock1;

    clock1=clock();

    meminit(basemap);
    for (i=0;i<4;++i) {
        basemap[("acgt")[i]]=i;
        basemap[("ACGT")[i]]=i;
        basemaprc[("tgca")[i]]=i;
        basemaprc[("TGCA")[i]]=i;
    }

    // max 16... lower values result in better sensitivity, lower specificity expecially with blended/mixed samples

    int k=12;
    int x=200;
    int sampsize=500000;
    float lowcom_pct=0.40;
    double ftop=.99;
    double fdel=.015;
    double fmin=.50;

    // run comparison on args?
    bool docompare=0;


    char c;
    char *endp;
    char **inf;
    int inc=0;

    char *contampath=NULL;

    while ( (c = getopt_long(argc, argv, "chk:x:n:l:t:d:m:r:S:Q:",NULL,NULL)) != -1) {
        switch (c) {
            case 'h': usage(stdout); return 0;
            case 'c': docompare=1; break;
            case 'x': x=atoi(optarg); break;
            case 'k': k=atoi(optarg); if(k>16) fail("Max k is 16\n"); break;
            case 'n': sampsize=atoi(optarg); break;
            case 'l': lowcom_pct=strtod(optarg,&endp); break;
            case 't': ftop=strtod(optarg,&endp); break;
            case 'S': scoef=strtod(optarg,&endp); break;
            case 'Q': qcoef=strtod(optarg,&endp); break;
            case 'r': contampath=optarg; break;
            case 'd': fdel=strtod(optarg,&endp); break;
            case 'm': fmin=strtod(optarg,&endp); break;
            case '\1': inf[++inc]=optarg; break;
            case '?':
                      if (!optopt) {
                          usage(stdout); return 0;
                      } else if (optopt && strchr("kxnltdmrSQ", optopt))
                          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                      else if (isprint(optopt))
                          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                      else
                          fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                      usage(stderr);
                      return 1;
        }
    }
    
    inf=&(argv[optind]);
    inc=argc-optind;

    if (!inc) {
        usage(stderr);
        return 1;
    }

    if (!x || !k || !sampsize || ftop<.8 || ftop > 1 || fdel > ftop || fmin > ftop || inc < 1) {
        fail("Bad params, quitting\n");
    }

    if (docompare) {
        compare_files(inf[0],inf[1]);
        exit(0);
    }

    char*in=inf[0];

    FILE *fin;
    if (contampath) {
        if (!strcmp(fext(in),".dz")||!strcmp(fext(in),".dsrc")) {
            fin=popenordie(string_format("dsrc d -s %s | bowtie -p 3 -S %s -q - | samtools view -S -f 4 - 2> /dev/null", in, contampath).c_str(), "r");
        } else {
            fin=popenordie(string_format("bowtie -p 3 -S %s -q %s | samtools view -S -f 4 - 2> /dev/null", contampath, in).c_str(), "r");
        }
    } else {
        fin=zopenordie(in,"r");
    }

    struct fq fq;

    meminit(fq);
    int n=0;

    google::sparse_hash_map<uint32_t, int> kh;

    std::string mer, scnt;
    mer.resize(k);
    int *pdat;
    int kcnt=0;
    unsigned long long lm; //, lmrc;
    unsigned long long kmask=~(((unsigned long long) ~0) << (k*2));
//    printf("kmask:%llx\n",kmask);

   
// minimum # of low-complexity positions in 32 bases to be called "low complexity"
    int lowcom_count = 32 * lowcom_pct;
    int skiplow=0;
    while (contampath?read_fq_sam(fin,n,&fq):read_fq(fin, n, &fq)) {

        if (fq.seq.n > 32)          // only first 32 bases considered
            fq.seq.n = 32;

        // sequence as 64-bit value
        lm=0;
//        lmrc=0;

        if (lowcomplex(fq.seq.s, fq.seq.n, lowcom_count)) {
            ++skiplow;
            continue;
        }

//        printf("s:%s\n", fq.seq.s);

// skip too short
        if (fq.seq.n < k) {
            ++skiplow;
            continue;
        }

        ++n;

        for(i=0;i<fq.seq.n;++i) {   // map to 64 bit sequence
            lm=lm|(basemap[fq.seq.s[i]]<<(i*2));
            // lmrc=(lmrc<<2)|basemaprc[fq.seq.s[i]];
        }
//        printf("orig: %s\n", string(fq.seq.s,32).c_str());
//        printf("s:%s,lm:%llx\n",fq.seq.s,lm);
        for (i=0;i<fq.seq.n-k+1;++i) {  // create unint_16_t kmers
            uint32_t v=lm&kmask;
            if (!kh[v]) ++kcnt;
            kh[v]++;
//            printf("i:%d,v:%hx\n",i,v);
            lm=lm>>2;               // remove first base

            // reverse

//            char t[33];
//            seq2char(v, k, t);

//            printf("orig: %s, lm: %llx, asv: %hx, mapped: %s\n", string(fq.seq.s+i,k).c_str(), lm, v, t);
        }
        if (n>=sampsize) break;
    }
    if (n==0) {
        fail("Insufficient sequence content\n");
    }

//    printf("entries: %d\n", kcnt);


    struct {uint32_t s; int v;} topx[x]; meminit(topx);

    google::sparse_hash_map<uint32_t,int>::iterator it = kh.begin();
    vector<int> dist;
    while(it!=kh.end()) {
        // *not* making the distinction between 2 singleton mappings and 1 paired here
        int val = it->second;
        dist.push_back(val);
        if (val > topx[x-1].v) {
            for (i=0;i<x;++i) {
                if (val > topx[i].v) {
                    break;
                }
            }
            int found=i;
            // shift remaining down (max 8 shifts)
            for (i=x-1;i>found;--i) {
                topx[i]=topx[i-1];
            }
            topx[found].v=val;
            topx[found].s=it->first;
        }
        ++it;
    }
    sort(dist.begin(),dist.end());

    double qtop=quantile(dist,ftop);


    printf("k\t%d\n", k);
    printf("top-x\t%d\n", x);
    printf("q-start\t%.3f\n", ftop);
    printf("q-delta\t%.3f\n", fdel);
    printf("q-min\t%.3f\n", fmin);
    printf("reads-used\t%d\n", n);
    printf("lowcom-pct\t%.3f\n", lowcom_pct*100);
    printf("quantile-top\t%f\n", qtop);
    printf("skip-low\t%d\n", skiplow);
    printf("run-time\t%f\n", (clock()-clock1)/(double)CLOCKS_PER_SEC);

    printf("--------\n");

    int j;
    char buf[k+1];
    for (i=0;i<x;++i) {
        seq2char(topx[i].s, k, buf);
        printf("%s %.3f\n", buf, topx[i].v/qtop);
    }

    printf("--------\n");

    double f;
    vector<double>dsig;
    for (f=ftop-fdel;f>=fmin;f-=fdel) {
        double r=quantile(dist, f)/qtop;
        dsig.push_back(r);
        printf("%.6f\n",r);
    }

}

FILE *zopenordie(const char *path, const char *mode) {
    bool isgz;
    FILE *f=gzopen(path, mode, &isgz);
    if (!f) {
        warn("Can't open-%s '%s': %s\n", mode, path, strerror(errno));
        exit(1);
    }
    return f;
}

FILE *openordie(const char *path, const char *mode) {
    bool isgz;
    FILE *f=fopen(path, mode);
    if (!f) {
        warn("Can't open-%s '%s': %s\n", mode, path, strerror(errno));
        exit(1);
    }
    return f;
}

FILE *popenordie(const char *cmd, const char *mode) {
    FILE *f=popen(cmd, mode);
    if (!f) {
        warn("Can't open-%s '%s': %s\n", mode, cmd, strerror(errno));
        exit(1);
    }
    return f;
}



void seq2char(uint32_t s, int k, char *p) {
    int i;
    for(i=0;i<k;++i) {
        *p++=("ACGT")[(s>>(i<<1))&0x3];
    }
    *p='\0';
}

bool lowcomplex(char *s, int n, int lowcom_thresh) {
    int lc=0;
    int i;

    for (i=1; i<n; ++i) {
        // N's always match everything
        if (s[i] == 'N' || (s[i] == s[i-1])) {
            lc+=1;
        } else if (i >= 3 && (s[i] == s[i-2] && s[i-1] == s[i-3])) {
            lc+=1;
        } else if (i >= 5 && (s[i] == s[i-3] && s[i-1] == s[i-4] && s[i-2] == s[i-5])) {
            lc+=1;
        } else if (i >= 7 && (s[i] == s[i-4] && s[i-1] == s[i-5] && s[i-2] == s[i-6] && s[i-3] == s[i-7])) {
            lc+=1;
        }
    }

// enough lc to pass?
    return lc > lowcom_thresh;
}

uint32_t char2seq(char *p, int k) {
    int i;
    uint32_t r=0;
    for(i=0;i<k;++i) {   // map to 64 bit sequence
        r=r|(basemap[p[i]]<<(i*2));
    }
    return r;
}

uint32_t char2seqrc(char *p, int k) {
    int i;
    uint32_t r=0;
    for(i=0;i<k;++i) {   // map to 64 bit sequence
        r=r<<2|(basemaprc[p[i]]);
    }
    return r;
}



class sig {
public:
    int k;                  // kmer size
    int svx;                // number of entries in sequence vector
    int n;                  // sample size
    double ftop;            // top quantile
    double fdel;            // delta
    double fmin;            // min quality
    double lowcom_pct;      // complexity filter used
    string fn;              // file read from (if any)
    int qvx;                // number of entries in quantile vector

    void read(const char *file);

    typedef struct {
        uint32_t seq;       // sequence
        uint32_t rcseq;     // reverse complement of same
        int lev;
    } ent;

    typedef struct {
        float sdist;        // s distance
        float qdist;        // q distance
        float score;      // blended score
    } comp;

    ent *svec;
    double *qvec;

    sig() {init();};
    ~sig() {clear();};

    void init() {svec=NULL;svx=qvx=n=k=0;qvec=NULL;};
    void clear() {if(svec) free(svec); if (qvec) free(qvec); init();}

    sig::comp compare(sig &other);
};

void sig::read(const char *file) {
    clear();

    FILE *f;
    bool isop=0;
    if (!strcmp(file,"-")) {
        f=stdin;
    } else {
        f=openordie(file,"r");
        isop=1;
    }

    fn=file;

    char *s=NULL; size_t a; int l;
    char *p, *e;
    double v;
    while((l=getline(&s, &a, f))>1) {
        if(p=strchr(s,'\t')) {
            *p++='\0';
            v=strtod(p, &e);
        } else {
            break;
        }

        if(p=strchr(s,':'))         // colon ignored 
            *p='\0';

        if (!strcmp(s, "k")) {
            k=v;
        } else if (!strcmp(s, "top-x")) {
            svx=v;
        } else if (!strcmp(s, "q-start")) {
            ftop=v;
        } else if (!strcmp(s, "q-delta")) {
            fdel=v;
        } else if (!strcmp(s, "q-min")) {
            fmin=v;
        } else if (!strcmp(s, "lowcom-pct")) {
            lowcom_pct=v/100;
        } else if (!strcmp(s, "reads-used")) {
            n=v;
        }
    }

    if (!k || !svx || ftop==0 || fdel==0 || fmin == 0 || n ==0) {
        fail("File '%s': Invalid header.\n", file);
    }

    svec = (sig::ent *)malloc(sizeof(*svec)*svx);
    int i = 0;
    while(i<svx && ( (l=getline(&s, &a, f)) > 1 )) {
        s[k]='\0';
        svec[i].seq=char2seq(s,k); 
        svec[i].rcseq=char2seqrc(s,k); 
        svec[i].lev=atoi(s+k+1);
        ++i;
    }

    // number of quantile vector entries
    qvx = (ftop-fmin-fdel)/fdel;

    // discard separator
    l=getline(&s, &a, f);

    qvec = (double *)malloc(sizeof(*qvec)*qvx);
    i = 0;
    while(i<qvx && ( (l=getline(&s, &a, f)) > 1 )) {
        qvec[i]=strtod(s,&e);
        ++i;
    }
    if (isop) fclose(f);
}

sig::comp sig::compare(sig &o) {
    if (svx!=o.svx || k!=o.k || fmin != o.fmin || fdel != o.fdel || ftop != o.ftop || lowcom_pct != o.lowcom_pct) {
        fail("Can't compare sigs with different parameters\n");
    }
    // influences certainty

    double diff=(double)abs(n-o.n)/(double(n+o.n)/2);

    // top x sequence vector comparison
    int i, j;
    int matches=0;
    double sdistance=0;
    for (i=0;i<svx;++i) {
        bool found=0;
        for (j=0;j<svx;++j) {
            // cap distance at 1000
            if (svec[i].seq==o.svec[j].seq || svec[i].seq==o.svec[j].rcseq) {
                found=1;
                matches+=1;
                // the max range on these is about 1000
                sdistance+=min(1000,fabs(svec[i].lev-o.svec[j].lev));
                break;
            }
        }
        if (!found) {
            // 10x worse than finding it at all...
            sdistance+=10000;
        }
    }
    sdistance=(sdistance/(double)svx)/10000;

    // top x quantile vector comparison
    double qdistance=0;
    for (i=1;i<qvx;++i) {
        // the biggest this can be is 1, but realistically, it will be very, very small
        qdistance+=fabs(fabs(qvec[i]-qvec[i-1])-fabs(o.qvec[i]-o.qvec[i-1]));
    }
    qdistance = qdistance/qvx;

    sig::comp c;
    c.sdist=sdistance;
    c.qdist=qdistance;

    // blended score - based on a fit to a model of expected outputs which included some mixtures
    c.score=max(0,min(1,1-scoef*sdistance-qcoef*qdistance));

    return c;
}


void compare_files(const char *i1, const char *i2) {
    sig s1;
    sig s2;

    s1.read(i1);
    s2.read(i2);

    sig::comp c = s1.compare(s2);

    printf("sequence distance\t%f\n", c.sdist);
    printf("quantile distance\t%f\n", c.qdist);
    printf("similarity\t%f\n", c.score);
}

template <class vtype>
double quantile(const vtype &vec, double p) {
        int l = vec.size();
        if (!l) return 0;
        double t = ((double)l-1)*p;
        int it = (int) t;
        int v=vec[it];
        if (t > (double)it) {
                return (v + (t-it) * (vec[it+1] - v));
        } else {
                return v;
        }
}


void usage(FILE *f) {
    fprintf(f,

"Usage: seqsig [options] file1\n"
"   or: seqsig -c file1 file2\n"
"Version: %s\n"
"\n"
"Generate qsig files for fastqs, or compare them\n"
"\n"
"Sig generation options:\n"
"    -k INT    kmer size (12)\n"
"    -x INT    number of top kmers to output (200)\n"
"    -n INT    number of reads to use (500000)\n"
"    -l REAL   low complexity filter (0.40), where 1 is no filter\n"
"    -t REAL   top quantile to use (0.90)\n"
"    -d REAL   quantile delta to iterate (0.15)\n"
"    -m REAL   minimum quantile to output/test (0.50)\n"
"    -r FASTA  bowtie indexed file of contamination/spike-ins to ignore\n"
"\n"
"Comparison options:\n"
"    -Q REAL   quantile coefficient for score output (350)\n"
"    -S REAL   sequence coefficient for score output (0.75)\n"
"\n"
"Misc options:\n"
"    -h        this help\n"
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


