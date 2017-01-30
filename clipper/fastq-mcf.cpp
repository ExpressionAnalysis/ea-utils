/*
Copyright (c) 2011 Expression Analysis / Erik Aronesty

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

/* 

See "void usage" below for usage.

*/

#include <sparsehash/sparse_hash_map> // or sparse_hash_set, dense_hash_map, ...
#include <string>

#include "fastq-lib.h"

#define VERSION "1.05"

#define MAX_ADAPTER_NUM 1000
#define SCANLEN 15
#define SCANMIDP ((int) SCANLEN/2)
#define MAX_FILES 5
#define MAX_REF 10
#define B_A     0
#define B_C     1
#define B_G     2
#define B_T     3
#define B_N     4
#define B_CNT   5
#define MAXWARN 10
#define MAX_PHRED 100

struct ad {
	char *id;  int nid;  size_t naid; 
	char *seq; int nseq; size_t naseq;
	char escan[SCANLEN+1]; 			// scan sequence
	int bcnt[MAX_FILES];			// number found at beginning
	int bcntz[MAX_FILES];			// number found at beginning
	int ecnt[MAX_FILES];			// number found at end
	int ecntz[MAX_FILES];			// number found at end

	char end[MAX_FILES];			// 'b' or 'e'
	int thr[MAX_FILES];			// min-length for clip
};

int read_fa(FILE *in, int rno, struct ad *ad);		// 0=done, 1=ok, -1=err+continue
int meanqwin(const char *q, int qn, int i, int w);     // mean quality within window win, at position i
bool evalqual(struct fq &fq, int file_num);

int char2bp(char c);
char bp2char(int b);
void saveskip(FILE **fout, int fo_n, struct fq *fq);

void valid_arg(char c, const char *a);

void usage(FILE *f, const char *msg=NULL);
int debug=0;
int warncount = 0;

// used to filter out other genomes, spike in controls, etc

const char *cmd_align_se = "bowtie -S %i -f %1";
const char *cmd_align_pe = "bowtie -S %i -1 %1 -2 %2";

// quality filter args
int qf_mean=0, qf_max_ns=-1, qf_xgt_num=0, qf_xgt_min=0, qf_max_n_pct=-1;
int qf2_mean=0, qf2_max_ns=-1, qf2_xgt_num=0, qf2_xgt_min=0, qf2_max_n_pct=-1;

// qual adjust
class adjustment {
public:
    int pos;
    int adj;
    adjustment() {pos=adj=0;}
};
std::vector<adjustment> cycle_adjust;
int phred_adjust[MAX_PHRED];
int phred_adjust_max=0;
bool have_phred_adjust=false;

std::string arg2cmdstr(int argc, char** argv);

// phred used
char phred = 0;

google::sparse_hash_map <std::string, int> dupset;
int dupmax = 40000000;              // this should be configurable, but right now it isn't
int max_in_buffer = 2400000;

class inbuffer {
    int max_buf;
public:
    inbuffer() {fin=0; gz=0; bp=0; max_buf=max_in_buffer;};
    ~inbuffer() {close();};

    FILE *fin;      
    bool gz;
    int bp;
    std::vector<std::string> buf;

    ssize_t getline(char **lineptr, size_t *n) {
        if (bp < buf.size()) {
            // return bufffered
            int l=buf[bp].length();         // length without null char
            if (!*lineptr || *n < (l+1)) {
                // alloc with room for null
                *lineptr=(char*)realloc(*lineptr,*n=(l+1));
            }
            memcpy(*lineptr,buf[bp].data(),l);
            (*lineptr)[l]='\0';
            ++bp;
            return l;
        } else {
            int l=::getline(lineptr, n, fin);
            if (max_buf > 0) {
                if (buf.size() > max_buf) {
					if (debug) fprintf(stderr, "Clearing buffer at %d lines\n", (int) buf.size());
                    buf.resize(0);
                    bp=0;
                    max_buf = 0;
                } else {
                    if (l > 0) {
                        buf.push_back(std::string(*lineptr));
                        ++bp;
                    }
                }
            }
            return l;
        }
    }

    int read_fq(int rno, struct fq *fq, const char *name=NULL) {
        if (bp < buf.size()) {
            fq->id.n=getline(&fq->id.s, &fq->id.a);
            fq->seq.n=getline(&fq->seq.s, &fq->seq.a);
            fq->com.n=getline(&fq->com.s, &fq->com.a);
            fq->qual.n=getline(&fq->qual.s, &fq->qual.a);
            if (fq->qual.n <= 0)
                    return 0;

            if (fq->id.s[0] != '@' || fq->com.s[0] != '+' || fq->seq.n != fq->qual.n) {
                    const char *errtyp = (fq->seq.n != fq->qual.n) ?  "length mismatch" : fq->id.s[0] != '@' ? "no '@' for id" : "no '+' for comment";
                    if (name) {
                        fprintf(stderr, "Malformed fastq record (%s) in file '%s', line %d\n", errtyp, name, rno*2+1);
                    } else {
                        fprintf(stderr, "Malformed fastq record (%s) at line %d\n", errtyp, rno*2+1);
                    }
                    return -1;
            }

            fq->seq.s[--fq->seq.n] = '\0';
            if (fq->seq.s[fq->seq.n-1] == '\r') {
                fq->seq.s[--fq->seq.n] = '\0';
            }
            fq->qual.s[--fq->qual.n] = '\0';
            if (fq->qual.s[fq->qual.n-1] == '\r') {
                fq->qual.s[--fq->qual.n] = '\0';
            }
 
            return fq->qual.n >= 0;  // github issue 46, 53
        } else {
            return ::read_fq(fin, rno, fq, name);
        }
    }

    void reset() {
        assert(max_buf > 0);
        bp=0;
    }

    bool full() {
        return buf.size()>=max_buf;
    }

    int close() {
       int ret=true;
       if (fin) {
            ret = gz ? pclose(fin) : fclose(fin);
            fin=NULL;
       }
        return ret;
    }
};

int main (int argc, char **argv) {
	char c;
	bool eol;
	int nmin = 1, nkeep = 19, nmax=0;
    int qf2_min_len=0;
	float minpct = 0.25;
	int pctdiff = 10;
	int sampcnt = 300000;			// # of reads to sample to determine adapter profile, and base skewing
	int xmax = -1;
	float scale = 2.2;
	int noclip=0;
        int nreadsout=0;         // max # of reads to output, all by default
	char end[MAX_FILES]; meminit(end);
	float skewpct = 2; 			// any base at any position is less than skewpct of reads
	float pctns = 20;			// any base that is more than 20% n's
	bool rmns = 1;				// remove n's at the end of the read
	int qthr = 7;				// remove end of-read with quality < qthr
	int qwin = 1;				// remove end of read with mean quality < qthr
	int ilv3 = -1;
	int duplen = 0;
	int dupskip = 0;
    int min_start_trim = 0;
    int min_end_trim = 0;
    bool noexec = 0;
    bool hompol_filter = 0;
    bool lowcom_filter = 0;
    float hompol_pct = .92;
    float lowcom_pct = .90;
    bool keeponlyclip=0;

    dupset.set_deleted_key("<>");

	int i;

	char *afil = NULL;
	char *ifil[MAX_FILES]; meminit(ifil);
	const char *ofil[MAX_FILES]; meminit(ofil);
	int i_n = 0;
	int o_n = 0;
	int e_n = 0;
	bool skipb = 0;
	char *fref[MAX_REF]; meminit(fref); 
	int fref_n = 0;
    char *qspec = NULL;

    static struct option long_options[] = {
       {"keep-clipped", 0, 0, 0},
       {"qual-mean", 1, 0, 0},
       {"max-ns", 1, 0, 0},
       {"max-output-reads", 1, 0, 'O'},
       {"qual-gt", 1, 0, 0},
       {"min-len", 1, 0, 'l'},
       {"cycle-adjust", 1, 0, 0},
       {"phred-adjust", 1, 0, 0},
       {"phred-adjust-max", 1, 0, 0},
       {"mate-qual-mean", 1, 0, 0},
       {"mate-max-ns", 1, 0, 0},
       {"mate-qual-gt", 1, 0, 0},
       {"mate-min-len", 1, 0, 0},
       {"homopolymer-pct", 1, 0, 0},
       {"lowcomplex-pct", 1, 0, 0},
       {"min-start-trim", 1, 0, 0},
       {"min-end-trim", 1, 0, 0},
       {0, 0, 0, 0}
    };

    meminit(phred_adjust);

    int option_index = 0;
    while (	(c = getopt_long(argc, argv, "-nf0uXUVHKSRdbehp:o:O:l:s:m:t:k:x:P:q:L:C:w:F:D:",long_options,&option_index)) != -1) {
		switch (c) {
			case '\0':
                { 
                    const char *oname=long_options[option_index].name;
                    if(!strcmp(oname,        "qual-mean")) {
                        qf_mean=qf2_mean=atoi(optarg);
                    } else if(!strcmp(oname,        "keep-clipped")) {
                        keeponlyclip=1;
                    } else if(!strcmp(oname, "mate-qual-mean")) {
                        qf2_mean=atoi(optarg);
                    } else if (!strcmp(oname, "min-start-trim")) {
        		min_start_trim = atoi(optarg);
                    } else if (!strcmp(oname, "min-end-trim")) {
        		min_end_trim = atoi(optarg);
                    } else if(!strcmp(oname, "homopolymer-pct")) {
                        hompol_pct=atof(optarg)/100.0;
                        hompol_filter=1;
                    } else if(!strcmp(oname, "lowcomplex-pct")) {
                        lowcom_pct=atof(optarg)/100.0;
                        lowcom_filter=1;
                    } else if(!strcmp(oname, "qual-gt")) {
                        if (!strchr(optarg, ',')) {
                            fprintf(stderr, "Error, %s requires NUM,THR as argument\n", oname);
                            exit(1);
                        }
                        qf_xgt_num=qf2_xgt_num=atoi(optarg);
                        qf_xgt_min=qf2_xgt_min=atoi(strchr(optarg, ',')+1);
                    } else if(!strcmp(oname, "mate-qual-gt")) {
                        if (!strchr(optarg, ',')) { 
                            fprintf(stderr, "Error, %s requires NUM,THR as argument\n", oname);
                            exit(1);
                        }
                        qf2_xgt_num=atoi(optarg);
                        qf2_xgt_min=atoi(strchr(optarg, ',')+1);
                    } else if(!strcmp(oname, "cycle-adjust")) {
                        if (!strchr(optarg, ',')) {
                            fprintf(stderr, "Error, %s requires CYC,ADJ as argument\n", oname);
                            exit(1);
                        }
                        adjustment a;
                        a.pos=atoi(optarg);
                        a.adj=atoi(strchr(optarg, ',')+1);
                        cycle_adjust.push_back(a);
                    } else if(!strcmp(oname, "phred-adjust-max")) {
                        phred_adjust_max=atoi(optarg);
                    } else if(!strcmp(oname, "phred-adjust")) {
                        if (!strchr(optarg, ',')) {
                            fprintf(stderr, "Error, %s requires CYC,ADJ as argument\n", oname);
                            exit(1);
                        }
                        int phred=atoi(optarg);
                        int adj=atoi(strchr(optarg, ',')+1);
                        assert(phred<MAX_PHRED && phred >= 0);
                        if (adj)
                            have_phred_adjust=true;
                        phred_adjust[phred]=adj;
                    } else if(!strcmp(oname, "max-ns")) {
                        if (strchr(optarg,'%')) {
                            qf_max_n_pct=atoi(optarg);
                            qf2_max_n_pct=atoi(optarg);
                        } else {
                            qf_max_ns=atoi(optarg);
                            qf2_max_ns=atoi(optarg);
                        }

                    } else if(!strcmp(oname, "mate-max-ns")) {
                        if (strchr(optarg,'%')) {
                            qf2_max_n_pct=atoi(optarg);
                        } else {
                            qf2_max_ns=atoi(optarg);
                        }
                    } else if(!strcmp(oname, "mate-min-len")) {
                        qf2_min_len=atoi(optarg);
                    }
                    break;
                }
			case '\1': 
				if (!afil) 
					afil = optarg; 
				else if (i_n<MAX_FILES) 
					ifil[i_n++] = optarg; 
				else {
					usage(stderr, "Too many input files."); return 1;
				}
				break;
			case 't': minpct = atof(optarg); break;
			case 'm': nmin = atoi(optarg); break;
			case 'l': nkeep = atoi(optarg); break;
			case 'L': nmax = atoi(optarg); break;
			case 'K': keeponlyclip=1; break;
			case '0': nmax=0; skewpct=0; pctns=0; rmns=0; qthr=0; nkeep=0; ilv3=-1;  break;
			case 'u': ilv3=1; break;
			case 'U': ilv3=0; break;
			case 'H': hompol_filter=1; break;
			case 'X': lowcom_filter=1; break;
			case 'k': skewpct = atof(optarg); break;
			case 'q': qthr = atoi(optarg); valid_arg(c,optarg); break;
			case 'Q': qspec = optarg; break;
			case 'w': qwin = atoi(optarg); break;
			case 'C': sampcnt = atoi(optarg); if (sampcnt*8 > max_in_buffer) max_in_buffer = sampcnt * 8; break;
			case 'F': fref[fref_n++] = optarg; break;
			case 'x': pctns = atof(optarg); break;
			case 'R': rmns = false; break;
			case 'V': printf("Version: %s\n", VERSION); return 0; break;
			case 'p': pctdiff = atoi(optarg); break;
			case 'P': phred = (char) atoi(optarg); break;
			case 'D': duplen = atoi(optarg); break;
			case 'h': usage(stdout); return 1; 
			case 'o': if (!o_n < MAX_FILES) 
						  ofil[o_n++] = optarg;
					  break;
               		case 'O': nreadsout = atoi(optarg); break;
			case 's': scale = atof(optarg); break;
			case 'S': skipb = 1; break;
			case 'i': if (i_n<MAX_FILES)
						  ifil[i_n++] = optarg; 
					  else
						  return usage(stderr, "Too many input files."), 1;
					  break;
			case 'n': noclip = 1; break;
			case 'd': ++debug; break;
			case 'b': end[e_n++] = 'b'; break;
			case 'e': end[e_n++] = 'e'; break;
			case '?': 
					  if (strchr("polsmtkx", optopt))
						  fprintf (stderr, "Option -%c requires an argument.\n", optopt);
					  else if (isprint(optopt))
						  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
					  else
						  fprintf (stderr,
								  "Unknown option character `\\x%x'.\n",
								  optopt);
					  usage(stderr);
					  return 1;
		}
	}

    if (duplen > 75) {
		fprintf(stderr, "WARNING: duplen of %d is probably too long, do you really need it?\n", duplen);
    }

	if (i_n == 1 && o_n == 0) {
		ofil[o_n++]="-";
	}

	if (!noclip && o_n != i_n) {
		fprintf(stderr, "Error: number of input files must match number of '-o' output files.\n");
		return 1;
	}

	if (argc < 3 || !afil || !i_n) {
		usage(stderr);
		return 1;
	}

	FILE *ain = NULL;
	if (strcasecmp(afil, "n/a") && strcasecmp(afil, "/dev/null") && strcasecmp(afil, "NUL")) {
		ain = fopen(afil, "r");
		if (!ain) {
			fprintf(stderr, "Error opening adapter file '%s': %s\n",afil, strerror(errno));
			return 1;
		}
	}

	FILE *fstat = stderr;
	if (!noclip && strcmp(ofil[0], "-")) {
		fstat = stdout;
	}
	if (noclip) {
		fstat = stdout;
	}

	fprintf(fstat, "Command Line: %s\n", arg2cmdstr(argc, argv).c_str());

	FILE *fout[MAX_FILES]; meminit(fout);
	bool gzout[MAX_FILES]; meminit(gzout);
    inbuffer fin[MAX_FILES];

   // if (debug) fprintf(stderr,"i_n:%d, ifil[0]:%s\n",i_n, ifil[0]);

	for (i=0;i<i_n;++i) {
	    if ((i_n==1) && !strcmp(ifil[0], "-")) {
            fin[i].fin=stdin;
            fin[i].gz=0;
        } else {
            fin[i].fin=gzopen(ifil[i], "r", &fin[i].gz);
        }
	}

	struct ad ad[MAX_ADAPTER_NUM+1];
	memset(ad, 0, sizeof(ad));

	int acnt=0, ok=0, rno=0;	// adapter count, ok flag, record number

	if (ain) {
		while (acnt < MAX_ADAPTER_NUM && (ok = read_fa(ain, rno, &ad[acnt]))) {
			++rno;
			if (ok < 0)
				break;
			// copy in truncated to max scan length
			strncpy(ad[acnt].escan, ad[acnt].seq, SCANLEN);
			ad[acnt].escan[SCANLEN] = '\0';
			//fprintf(stderr, "escan: %s, %s\n", ad[acnt].id, ad[acnt].escan);
			++acnt;
		}

		if (acnt == 0) {
			fprintf(stderr, "No adapters in file '%s'\n",afil);
		}
	}

	fprintf(fstat, "Scale used: %g\n", scale);
	int maxns = 0;						// max sequence length
	int avgns[MAX_FILES]; meminit(avgns);			// average sequence length per file
	// read length
	for (i=0;i<i_n;++i) {

		char *s = NULL; size_t na = 0; int nr = 0, ns = 0, totn[MAX_FILES]; meminit(totn);
		char *q = NULL; size_t naq = 0; int nq =0;
		int j;
		int ilv3det=2;
        int skipped = 0;

        struct stat st;
        stat(ifil[i], &st);

		while (fin[i].getline(&s, &na) > 0) {
			if (*s == '@')  {
				// look for illumina purity filtering flags
				if (ilv3det==2) {
					ilv3det=0;
					const char *p=strchr(s, ':');
					if (p) {
						++p;
						if (isdigit(*p)) {
							p=strchr(s, ' ');
							if (p) {
								++p;
								if (isdigit(*p)) {
									++p;
									if (*p ==':') {
										++p;
										if (*p =='Y') {
											// filtering found
											ilv3det=1;
										} else if (*p =='N') {
											// still illumina
											ilv3det=2;
										}
									}
								}
							}
						}
					}
				}

				if ((ns=fin[i].getline(&s, &na)) <=0) {
					// reached EOF
					if (debug) fprintf(stderr, "Dropping out of sampling loop\n");
					break;
				}

				nq=fin[i].getline(&q, &naq);
				nq=fin[i].getline(&q, &naq);		// qual is 2 lines down

				// skip poor quals/lots of N's when doing sampling
				if (st.st_size > (sampcnt * 500) && (skipped < sampcnt) && poorqual(i, ns, s, q)) {
					if (debug) fprintf(stderr, "Skip poorqual\n");
                    ++skipped;
					continue;
                }

				if (phred == 0) {
					--nq;
					for (j=0;j<nq;++j) {
						if (q[j] < 64) {
							if (debug) fprintf(stderr, "Using phred 33, because saw: %c\n", q[j]);
							// default to sanger 33, if you see a qual < 64
							phred = 33;
							break;
						}
					}
				}
				--ns;                                   // don't count newline for read len
				++nr;
				avgns[i] += ns;
				if (ns > maxns) maxns = ns;

				// just 10000 reads for readlength sampling
				if (nr >= 10000) {	
					if (debug) fprintf(stderr, "Read 10000\n");
					break;
				}
			} else {
				fprintf(stderr, "Invalid FASTQ format : %s\n", ifil[i]);
				break;
			}
		}
		if (ilv3det == 1 && (ilv3 == -1)) {
			ilv3=1;
		}
		if (debug) fprintf(stderr,"Ilv3det: %d\n", ilv3det);
		if (s) free(s);
		if (q) free(q);
		if (nr)
			avgns[i] = avgns[i]/nr;
	}

	if (ilv3 == -1) {
		ilv3 = 0;
	}

	if (ilv3) {
		fprintf(fstat, "Filtering Illumina reads on purity field\n");
	}

	// default to illumina 64 if you never saw a qual < 33
	if (phred == 0) phred = 64;
	fprintf(fstat, "Phred: %d\n", phred);

	for (i=0;i<i_n;++i) {
		if (avgns[i] == 0) {
			fprintf(stderr, "No records in file %s\n", ifil[i]);
			exit(1);
		}
	}

	for (i=0;i<i_n;++i) {
        fin[i].reset();
	}

	if (debug) fprintf(stderr,"Max ns: %d, Avg[0]: %d\n", maxns, avgns[0]);

	// total base count per read position in sample
	int balloc = maxns;
	bool dobcnt = 1;
	if (maxns > 500) {
		dobcnt = 0;
		balloc = 1;
	}

	int bcnt[MAX_FILES][2][balloc][6]; meminit(bcnt);
	int qcnt[MAX_FILES][2]; meminit(qcnt);
	char qmin=127, qmax=0;
	int nsampcnt = 0;
    double stat_lowcom_total=0, stat_lowcom_ssq=0, stat_lowcom_b4_total=0, stat_lowcom_b4_ssq=0;
    long stat_lowcom_cnt=0, stat_lowcom_b4_cnt=0;
    int skipunclip=0;

	for (i=0;i<i_n;++i) {

		struct stat st;
		stat(ifil[i], &st);

		// todo, use readfq
		char *s = NULL; size_t na = 0; int ns = 0, nr = 0;
		char *q = NULL; size_t naq = 0; int nq =0;
		char *d = NULL; size_t nad = 0; int nd =0;

        int skipped = 0;
		while ((nd=fin[i].getline(&d, &nad)) > 0) {
			if (*d == '@')  {
				if ((ns=fin[i].getline(&s, &na)) <=0) 
					break;
				nq=fin[i].getline(&q, &naq);
				nq=fin[i].getline(&q, &naq);		// qual is 2 lines down

				--nq; --ns;				// don't count newline for read len

				// skip poor quals/lots of N's when doing sampling (otherwise you'll miss some)
                if (ns == 0) {  // github issue 46, 53
                    ++skipped;
                    continue;
                }
				if ((st.st_size > (sampcnt * 500)) && (skipped < sampcnt) && poorqual(i, ns, s, q)) {
                    ++skipped;
					continue;
                }

				if (nq != ns) {
					if (warncount < MAXWARN) {
						fprintf(stderr, "Warning, corrupt quality for sequence: %s", s);
						++warncount;
					}
					continue;
				}

				if (i > 0 && avgns[i] < 11) 			// reads of avg length < 11 ? barcode lane, skip it
					continue;

				if (ilv3) {					// illumina purity filtering
					char * p = strchr(d, ' ');
					if (p) {
						p+=2;
						if (*p==':') {
							++p;
							if (*p == 'Y') {
								continue;
							}
						}
					}
				}

				++nr;

				// to be safe, we don't assume reads are fixed-length, not any slower, just a little more code
				if (dobcnt) {
					int b;
					for (b = 0; b < ns/2 && b < maxns; ++b) {
						++bcnt[i][0][b][char2bp(s[b])];		// count from begin
						++bcnt[i][0][b][B_CNT];			// count of samples at position
						++bcnt[i][1][b][char2bp(s[ns-b-1])];	// count from end
						++bcnt[i][1][b][B_CNT];			// count of samples at offset-from-end position
					}
				}
				qcnt[i][0]+=((q[0]-phred)<qthr);		// count of q<thr for last (first trimmable) base
				qcnt[i][1]+=((q[ns-1]-phred)<qthr);	
				//fprintf(stderr,"qcnt i%d e0=%d, e1=%d\n", i, qcnt[i][0], qcnt[i][1]);

				// BUF conains only the first 15 characters of the sequence
				int a;
				char buf[SCANLEN+1];
				strncpy(buf, s, SCANLEN);
				buf[SCANLEN]='\0';
				for(a=0;a<acnt;++a) {
					char *p;
					// search whole seq for 15 char "end" of adap string
					if (p = strstr(s+1, ad[a].escan)) { 
						if (debug > 1) fprintf(stderr, "  END S: %s A: %s (%s), P: %d, SL: %d, Z:%d\n", s, ad[a].id, ad[a].escan, (int) (p-s), ns, (p-s) == ns-SCANLEN);
                        // found at the very end
						if ((p-s) == ns-SCANLEN) 
							++ad[a].ecntz[i];
						++ad[a].ecnt[i];
					}
					// search 15 char begin of seq in longer adap string
					int slen;
					if (debug > 1) fprintf(stderr, "COMPARE: %d <= %d, adseq: %s, buf: %s\n", SCANLEN, ad[a].nseq, ad[a].seq, buf);
					// if the 15bp sequence is smaller than the adapter size
					if (SCANLEN <= ad[a].nseq) {
						slen = SCANLEN;
						// search for the truncated buffer in the ADAPTER
						p = strstr(ad[a].seq, buf);
					} else {
						// search for the adapter at the beginning of the buffer only... if it's short
						slen = ad[a].nseq;
						if (!strncmp(ad[a].seq,buf,ad[a].nseq)) 
							p=ad[a].seq;
						else
							p=NULL;
					}
					if (p) { 
						if (debug > 1) fprintf(stderr, "BEGIN S: %s A: %s (%s), P: %d, SL: %d, Z:%d\n", buf, ad[a].id, ad[a].seq, (int) (p-ad[a].seq), ns, (p-ad[a].seq )  == ad[a].nseq-slen);
                        // found the end of the adapter
						if (p-ad[a].seq == ad[a].nseq-slen) 
							++ad[a].bcntz[i];
						++ad[a].bcnt[i];
					}
				}
			}
			if (fin[i].full() || nr >= sampcnt)		// enough samples 
				break;
		}
		if (s) free(s);
		if (d) free(d);
		if (q) free(q);
		if (i == 0 || avgns[i] >= 11) {
			if (nsampcnt == 0 || nr < nsampcnt)			// fewer than max, set for thresholds
				nsampcnt=nr;
		}
	}

	if (nsampcnt == 0) {
		fprintf(stderr, "ERROR: Unable to read file for subsampling\n");
		exit(1);
	}

	sampcnt = nsampcnt;
	int sktrim[i_n][2]; meminit(sktrim);

	// look for severe base skew, and auto-trim ends based on it
	int needqtrim=0;
	if (dobcnt) {
	if (sampcnt > 0 && skewpct > 0) {
		for (i=0;i<i_n;++i) {
			if (avgns[i] < 11) 			// reads of avg length < 11 ? barcode lane, skip it
				continue;
			int e;
			for (e = 0; e < 2; ++e) {
				// 5% qual less than low-threshold?  need qualtrim
				if (qthr > 0 && (100.0*qcnt[i][e])/sampcnt > 5) {
					needqtrim = 1;
				}

				int p;
				for (p = 0; p < maxns/2; ++p) {
					int b;

					int skth = (int) ( (float) bcnt[i][e][p][B_CNT] * ( skewpct / 100.0 ) ) ;	// skew threshold
					int thr_n = (int) ( (float) bcnt[i][e][p][B_CNT] * ( pctns / 100.0 ) );		// n-threshold

					if (debug > 1) 
						fprintf(stderr,"Sk Prof [%d, %d]: skth=%d, bcnt=%d, ncnt=%d, a=%d, c=%d, g=%d, t=%d\n", e, p, skth, 
								bcnt[i][e][p][B_CNT], bcnt[i][e][p][B_N], bcnt[i][e][p][B_A], 
								bcnt[i][e][p][B_C], bcnt[i][e][p][B_G], bcnt[i][e][p][B_T]);

					if (skth < 10)						// too few samples to detect skew
						continue;

					int tr = 0;
					for (b = 0; b < 4; ++b) {
						if (bcnt[i][e][p][b] < skth) {			// too few bases of this type
							tr=1;
							if (debug > 1) 
								fprintf(stderr, "Skew at i:%d e:%d p:%d b:%d\n", i, e, p, b);
							break;
						}
					}
					if (bcnt[i][e][p][B_N] > thr_n) {			// too many n's
						if (debug > 1) 
							fprintf(stderr, "Too many N's at i:%d e:%d p:%d b:%d ( %d > %d )\n", i, e, p, b, bcnt[i][e][p][B_N], thr_n);
						tr=1;
					}

					if (tr) {
						if (p == sktrim[i][e]) {				// adjacent, so increase trim
							++sktrim[i][e];
						} else {
							fprintf(fstat, "Within-read Skew: Position %d from the %s of reads is skewed!\n", p, e==0?"start":"end");
						}
					}
				}
			}
		}
	}

	}

	int e;
	bool someskew = false;
	for (i=0;i<i_n;++i) {
		int totskew = sktrim[i][0] + sktrim[i][1];
		if ((maxns - totskew) < nkeep) {
			if (totskew > 0) {
				fprintf(fstat, "Warning: Too much skewing found (%d), disabling skew clipping\n", totskew);
			}
			meminit(sktrim);
			break;
		}
	}

	for (i=0;i<i_n;++i) {
		for (e=0;e<2;++e) {
			if (sktrim[i][e] > 0) {
				fprintf(fstat, "Trim '%s': %d from %s\n",  e==0?"start":"end", sktrim[i][e], ifil[i]);
				someskew=true;
			}
		}
	}

	int athr = (int) ((float)sampcnt * minpct) / 100;
	fprintf(fstat, "Threshold used: %d out of %d\n", athr+1, sampcnt);

	int a;
	int newc=0;
	for(a=0;a<acnt;++a) {
		int any=0;
		for (i=0;i<i_n;++i) {
			if (debug) fprintf(stderr, "ad:%s, EC:%d, BC:%d, ECZ: %d, BCZ: %d\n", ad[a].id, ad[a].ecnt[i], ad[a].bcnt[i], ad[a].ecntz[i], ad[a].bcntz[i]);
			if (ad[a].ecnt[i] > athr || ad[a].bcnt[i] > athr) {
				int cnt;
				// heavily weighted toward start/end maches
				if ((ad[a].ecnt[i] + 10*ad[a].ecntz[i]) >= (ad[a].bcnt[i] + 10*ad[a].bcntz[i])) {
					ad[a].end[i]='e';
					cnt = ad[a].ecnt[i];
				} else {
					ad[a].end[i]='b';
					cnt = ad[a].bcnt[i];
				}

                char *p;
                if (p=strstr(ad[a].id, "_3p")) {
                    if (p[3] == '\0' || p[3] == '_') {
                        ad[a].end[i]='e'; 
					    cnt = ad[a].ecnt[i];
                    }
                } else if (p=strstr(ad[a].id, "_5p")) {
                    if (p[3] == '\0' || p[3] == '_') {
                        ad[a].end[i]='b';
                        cnt = ad[a].bcnt[i];
                    }
                }

				// user supplied end.... don't clip elsewhere
				if (end[i] && ad[a].end[i] != end[i])
					continue;

				if (scale >= 100) 
					ad[a].thr[i] = ad[a].nseq;
				else
					ad[a].thr[i] = min(ad[a].nseq,max(nmin,(int) (-log(cnt / (float) sampcnt)/log(scale))));

				fprintf(fstat, "Adapter %s (%s): counted %d at the '%s' of '%s', clip set to %d", ad[a].id, ad[a].seq, cnt, ad[a].end[i] == 'e' ? "end" : "start", ifil[i], ad[a].thr[i]);
				if (abs((ad[a].bcnt[i]-ad[a].ecnt[i])) < athr/4) {
					fprintf(fstat, ", warning end was not reliable: %s/%s\n", ad[a].id, ad[a].seq);
				} else {
					fputc('\n', fstat);
				}
				++any;
			}
		}
		if (!any) 
			continue;
		ad[newc++]=ad[a];
	}

	acnt=newc;

	if (acnt == 0 && !someskew && !needqtrim && !ilv3) {
		fprintf(fstat, "No adapters found");
		if (skewpct > 0) fprintf(fstat, ", no skewing detected"); 
		if (qthr > 0) fprintf(fstat, ", and no trimming needed");
		fprintf(fstat, ".\n");
		if (noclip) exit (1);			// for including in a test
	} else {
		if (debug) fprintf(stderr, "acnt: %d, ssk: %d, needq: %d\n", acnt, someskew, needqtrim);
		if (noclip) {
			if (acnt == 0) fprintf(fstat, "No adapters found. ");
			if (someskew) fprintf(fstat, "Skewing detected. "); 
			if (needqtrim) fprintf(fstat, "Quality trimming is needed. ");
			fprintf(fstat, "\n");
		}
	}

	if (noclip)
		exit(0);

	for (i=0;i<o_n;++i) {
		if (!strcmp(ofil[i],"-")) {
			fout[i]=stdout;
		} else {
			fout[i]=gzopen(ofil[i], "w", &gzout[i]);
		}
	}

	FILE *fskip[MAX_FILES]; meminit(fskip);
	bool gzskip[MAX_FILES]; meminit(gzskip);

	if (skipb) {
		for (i=0;i<o_n;++i) {
			if (!strcmp(ofil[i],"-")) {
				fskip[i]=stderr;
			} else {
				char *skipfil = (char *) malloc(strlen(ofil[i])+10);
				if (!strcmp(fext(ofil[i]),".gz")) {
					char *p=(char *)strrchr(ofil[i],'.');
					*p='\0';
					sprintf(skipfil, "%s.skip.gz", ofil[i]);
					*p='.';
				} else {
					sprintf(skipfil, "%s.skip", ofil[i]);
				}
				if (!(fskip[i]=gzopen(skipfil, "w", &gzskip[i]))) {
					fprintf(stderr, "Error opening skip file '%s': %s\n",skipfil, strerror(errno));
					return 1;
				}
				free(skipfil);
			}
		}
	}

	struct fq fq[MAX_FILES];	
	memset(&fq, 0, sizeof(fq));

	int nrec=0;
	int wrec=0;
	int nerr=0;
	int nok=0;
	int ntooshort=0;
	int ntoohompol=0;
	int ntoolowcom=0;
	int nfiltered=0;
	// total per read
	int ntrim[MAX_FILES]; meminit(ntrim);

	// total per end
	int cnttrim[MAX_FILES][2]; meminit(cnttrim);
	double tottrim[MAX_FILES][2]; meminit(tottrim);
	double ssqtrim[MAX_FILES][2]; meminit(ssqtrim);
	int trimql[MAX_FILES]; meminit(trimql);
	int trimqb[MAX_FILES]; meminit(trimqb);
	int nilv3pf=0;	// number of illumina version 3 purity filitered
	int read_ok;

	if (i_n > 0)
		fprintf(fstat, "Files: %d\n", i_n);

	for (i=0;i<i_n;++i) {
        fin[i].reset();
	}

    google::sparse_hash_map <std::string, int>::const_iterator lookup_it;

    bool io_ok = true;
    while (read_ok=fin[0].read_fq(nrec, &fq[0])) {
                if (nreadsout && (wrec == nreadsout)) break;
		for (i=1;i<i_n;++i) {
			int mok=fin[1].read_fq(nrec, &fq[i]);
			if (mok != read_ok) {
				fprintf(stderr, "# of rows in mate file '%s' doesn't match, quitting!\n", ifil[i]);
				return 1;
			}
		}
		++nrec;
        if (fq[0].qual.n == 0) { // github issue 46, 53
            ++nfiltered;
            continue;
        } else if (i_n > 1) {
            if (fq[1].qual.n == 0) {
                ++nfiltered;
                continue;
            }
        }
		if (read_ok < 0) {
			++nerr;
			continue;
		}

		if (ilv3) {
			char * p = strchr(fq[0].id.s, ' ');
			if (p) {
				p+=2;
				if (*p==':') {
					++p;
					if (*p == 'Y') {
						++nilv3pf;
						if (skipb) saveskip(fskip, i_n, fq);
						continue;
					}
				}
			}
		}

		// chomp

		int dotrim[MAX_FILES][2];
		int skip = 0;							// skip whole record?
		int hompol_seq=0;
		int hompol_cnt=0;
		int lowcom_seq=0;
		int lowcom_cnt=0;
		int f;	
        bool didclip=0;
		for (f=0;f<i_n;++f) {
			dotrim[f][0] = sktrim[f][0];					// default, trim to detected skew levels
			dotrim[f][1] = sktrim[f][1];
			// trim to minimum, if specified
			dotrim[f][0] = max(dotrim[f][0], min_start_trim);
			dotrim[f][1] = max(dotrim[f][1], min_end_trim);
			if (avgns[f] < 11)  
				// reads of avg length < 11 ? barcode lane, skip it
				continue;


            if (have_phred_adjust) {
                for (i=0;i<fq[f].qual.n;++i) {
                   if (phred_adjust[fq[f].qual.s[i]-phred]) {
                        fq[f].qual.s[i]+=phred_adjust[fq[f].qual.s[i]-phred];
                   } 
                }
            }

            if (phred_adjust_max) {
                for (i=0;i<fq[f].qual.n;++i) {
                   if ((fq[f].qual.s[i]-phred)>phred_adjust_max) {
                        fq[f].qual.s[i]=phred_adjust_max+phred;
                   } 
                }
            }


            for (i=0;i<cycle_adjust.size();++i) {
                if (abs(cycle_adjust[i].pos) < fq[f].qual.n) {
                    if (cycle_adjust[i].pos>0) {
                        fq[f].qual.s[cycle_adjust[i].pos-1]+=cycle_adjust[i].adj;
                    } else {
                        fq[f].qual.s[fq[f].qual.n+cycle_adjust[i].pos]+=cycle_adjust[i].adj;
                    }
                }
            }


			if (rmns) {
				for (i=dotrim[f][0];i<(fq[f].seq.n);++i) {
					// trim N's from the front
					if (fq[f].seq.s[i] == 'N') 
						dotrim[f][0] = i + 1;
					else
						break;
				}
				for (i=dotrim[f][1];i<(fq[f].seq.n);++i) {
					// trim N's from the end
					if (fq[f].seq.s[fq[f].seq.n-i-1] == 'N')
						dotrim[f][1] = i + 1;
					else 
						break;
				}
			}

            if (hompol_filter) {
                char p; int h = 0;
                for (i = dotrim[f][0]+1;i<fq[f].seq.n;++i) {
                    // N's always match everything
                    if (fq[f].seq.s[i] == 'N' || (fq[f].seq.s[i] == fq[f].seq.s[i-1])) {
                        ++hompol_seq;
                    }
                    ++hompol_cnt;
                }
            }

            if (lowcom_filter) {
                char p; int h = 0;
                for (i = dotrim[f][0]+1;i<fq[f].seq.n;++i) {
                    // N's always match everything
                    if (fq[f].seq.s[i] == 'N' || (fq[f].seq.s[i] == fq[f].seq.s[i-1])) {
                        ++lowcom_seq;
                    } else if (i >= dotrim[f][0]+3 && (fq[f].seq.s[i] == fq[f].seq.s[i-2] && fq[f].seq.s[i-1] == fq[f].seq.s[i-3])) {
                        ++lowcom_seq;
                    } else if (i >= dotrim[f][0]+5 && (fq[f].seq.s[i] == fq[f].seq.s[i-3] && fq[f].seq.s[i-1] == fq[f].seq.s[i-4] && fq[f].seq.s[i-2] == fq[f].seq.s[i-5])) {
                        ++lowcom_seq;
                    } else if (i >= dotrim[f][0]+7 && (fq[f].seq.s[i] == fq[f].seq.s[i-4] && fq[f].seq.s[i-1] == fq[f].seq.s[i-5] && fq[f].seq.s[i-2] == fq[f].seq.s[i-6] && fq[f].seq.s[i-3] == fq[f].seq.s[i-7])) {
                        ++lowcom_seq;
                    } 
                    ++lowcom_cnt;
                }
            }

			if (qthr > 0) {
				bool istrimq = false;

				// trim qual from the begin
				for (i=dotrim[f][0];i<(fq[f].seq.n);++i) {
					if (qwin > 1 && (meanqwin(fq[f].qual.s,fq[f].seq.n,i,qwin)-phred) < qthr) {
						++trimqb[f];
						istrimq = true;
						dotrim[f][0] = i + 1;
					} else if ((fq[f].qual.s[i]-phred) < qthr) {
						++trimqb[f];
						istrimq = true;
						dotrim[f][0] = i + 1;
					} else
						break;
				}


                // trim qual from the end ... stop at what you trimmed from the front!
				for (i=dotrim[f][1];i<(fq[f].seq.n-dotrim[f][0]);++i) {
					if (qwin > 1 && (meanqwin(fq[f].qual.s,fq[f].seq.n,fq[f].seq.n-i-1,qwin)-phred) < qthr) {
						++trimqb[f];
						istrimq = true;
						dotrim[f][1] = i + 1;
					} else if ((fq[f].qual.s[fq[f].seq.n-i-1]-phred) < qthr) {
						++trimqb[f];
						istrimq = true;
						dotrim[f][1] = i + 1;
					} else 
						break;
				}

                // denominator
				if (istrimq) trimql[f]+=1;
			}

			int bestscore_e = INT_MAX, bestoff_e = 0, bestlen_e = 0; 
			int bestscore_b = INT_MAX, bestoff_b = 0, bestlen_b = 0; 

			for (i =0; i < acnt; ++i) {
				if (debug) fprintf(stderr, "seq[%d]: %s %d\n", f, fq[f].seq.s, fq[f].seq.n);

				if (!ad[i].end[f])
					continue;

				int nmatch = ad[i].thr[f];
				if (!nmatch) nmatch = ad[i].nseq;			// full match required if nmin == 0

				// how far in to search for a match?
				int mx = ad[i].nseq;
				if (xmax) {
					mx = fq[f].seq.n;
					if (xmax > 0 && (xmax+ad[i].nseq) < mx)
						mx = xmax+ad[i].nseq;			// xmax is added to adapter length
				}

				if (debug)
					fprintf(stderr, "adapter: %s, adlen: %d, nmatch: %d, mx: %d\n", ad[i].seq, ad[i].nseq, nmatch, mx);

				if (ad[i].end[f] == 'e') {
					int off;
					for (off = nmatch; off <= mx; ++off) {		// off is distance from tail of sequence
						char *seqtail = fq[f].seq.s+fq[f].seq.n-off; 	// search at tail
						int ncmp = off<ad[i].nseq ? off : ad[i].nseq;
						int mind = (pctdiff * ncmp) / 100;
						int d = hd(ad[i].seq,seqtail,ncmp);		// # differences
						if (debug>1)
							fprintf(stderr, "tail: %s, bestoff: %d, off: %d, ncmp: %d, mind: %d, hd %d\n", seqtail, bestoff_e, off, ncmp, mind, d);
						if (d <= mind) {
							// squared-distance over length
							int score = (1000*(d*d+1))/ncmp;
							if (score <= bestscore_e) {			// better score?
								bestscore_e = score;			// save max score
								bestoff_e = off;			// offset at max
								bestlen_e = ncmp;			// cmp length at max
							}
							if (d == 0 && (ncmp == ad[i].nseq)) {
								break;
							}
						}
					}
				} else {
					int off;
					for (off = nmatch; off <= mx; ++off) {              // off is distance from start of sequence
						int ncmp = off<ad[i].nseq ? off : ad[i].nseq;	// number we are comparing
						char *matchtail = ad[i].seq+ad[i].nseq-ncmp;    // tail of adapter
						char *seqstart = fq[f].seq.s+off-ncmp;		// offset into sequence (if any)
						int mind = (pctdiff * ncmp) / 100;
						int d = hd(matchtail,seqstart,ncmp);            // # differences
						if (debug>1)
							fprintf(stderr, "bestoff: %d, off: %d, ncmp: %d, mind: %d, hd %d\n", bestoff_e, off, ncmp, mind, d);

						if (d <= mind) {
							int score = (1000*(d*d+1))/ncmp;
							if (score <= bestscore_b) {                       // better score?
								bestscore_b = score;                      // save max score
								bestoff_b = off;                          // offset at max
								bestlen_b = ncmp;                         // cmp length at max
							}
							if (d == 0 && (ncmp == ad[i].nseq)) {
								break;
							}
						}
					}
				}
		    }
	
			int adapcliplen = bestoff_b ? bestoff_b : bestoff_e;

			// lengthen trim based on best level
			if (bestoff_b > dotrim[f][0])
				dotrim[f][0]=bestoff_b;

			if (bestoff_e > dotrim[f][1])
				dotrim[f][1]=bestoff_e;

			int totclip = min(fq[f].seq.n,dotrim[f][0] + dotrim[f][1]);

//			if (debug > 1) fprintf(stderr,"totclip %d\n", totclip);

            if (totclip > 0) {
                // keep length > X, X based on mate
                int tkeep = f == 0 ? nkeep : qf2_min_len > 0 ? qf2_min_len : nkeep;

				if ( (fq[f].seq.n-totclip) < tkeep) {
					// skip all reads if one is severely truncated ??
					// maybe not... ?
					skip = 1;
					break;
				}

				// count number of adapters clipped, not the number of rows trimmed
				if ( adapcliplen > 0 ) {
					++ntrim[f];
                    didclip=1;
                }

				// save some stats
				if (bestoff_b > 0) {
					cnttrim[f][0]++;
					tottrim[f][0]+=bestoff_b;
					ssqtrim[f][0]+=bestoff_b*bestoff_b;
				} else if (bestoff_e > 0) {
					cnttrim[f][1]++;
					tottrim[f][1]+=bestoff_e;
					ssqtrim[f][1]+=bestoff_e*bestoff_e;
				}

			} else {
				// skip even if the original was too short
				if (fq[f].seq.n < nkeep) 
					skip = 1;
			}
		}

        if (keeponlyclip && !didclip) {
            ++skipunclip;
            skip=1;
        }

        int hompol_skip=0;
        if (hompol_filter) {
            int hompol_max = hompol_pct * hompol_cnt;
            if (debug>0) printf("%s: hompol cnt:%d, max:%d, seq:%d\n", fq[0].id.s, hompol_cnt, hompol_max, hompol_seq);
            if (hompol_seq>=hompol_max) hompol_skip = skip = true;
        }

        int lowcom_skip=0;
        if (!hompol_skip && lowcom_filter) {
            int lowcom_max = lowcom_pct * lowcom_cnt;
            if (debug>0) printf("%s: lowcom cnt:%d, max:%d, seq:%d\n", fq[0].id.s, lowcom_cnt, lowcom_max, lowcom_seq);
            if (lowcom_seq>=lowcom_max) lowcom_skip = skip = true;
            if (!lowcom_skip) { 
                stat_lowcom_total+=((double)lowcom_seq/(double)lowcom_cnt);
                stat_lowcom_ssq+=pow(((double)lowcom_seq/(double)lowcom_cnt),2);
                stat_lowcom_cnt+=1;
            }
            stat_lowcom_b4_total+=((double)lowcom_seq/(double)lowcom_cnt);
            stat_lowcom_b4_ssq+=pow(((double)lowcom_seq/(double)lowcom_cnt),2);
            stat_lowcom_b4_cnt+=1;
        }


		if (!skip) {
			int f;
			for (f=0;f<o_n;++f) {
                if (dotrim[f][1] >= strlen(fq[f].seq.s)) {
					if (debug) fprintf(stderr,"trimmming full sequence from end (%d), %s", dotrim[f][1], fq[f].id.s);
                    skip=1;
                    continue;
                }
				if (dotrim[f][1] > 0) {
					if (debug) fprintf(stderr,"trimming %d from end, %s", dotrim[f][1], fq[f].id.s);
					fq[f].seq.s[fq[f].seq.n -=dotrim[f][1]]='\0';
					fq[f].qual.s[fq[f].qual.n-=dotrim[f][1]]='\0';
				}
				if (dotrim[f][0] > 0) {
					if (debug) fprintf(stderr,"trimming %d from begin, %s", dotrim[f][0], fq[f].id.s);
					fq[f].seq.n -= dotrim[f][0];
					fq[f].qual.n -= dotrim[f][0];
                    if (fq[f].seq.n < 0) {
                        fq[f].seq.n = 0;
                        fq[f].qual.n = 0;
                    }
					memmove(fq[f].seq.s ,fq[f].seq.s +dotrim[f][0],fq[f].seq.n );
					memmove(fq[f].qual.s,fq[f].qual.s+dotrim[f][0],fq[f].qual.n);
					fq[f].seq.s[fq[f].seq.n]='\0';
					fq[f].qual.s[fq[f].qual.n]='\0';
				}
				if (nmax > 0) {
					if (fq[f].seq.n >= nmax ) {
						fq[f].seq.s[nmax]='\0';
						fq[f].qual.s[nmax]='\0';
					}
				}
                if (avgns[f]>=11 && !evalqual(fq[f],f)) {
                    skip = 2;                       // 2==qual
                }
			}

            if (duplen > 0 && !skip) {
                // lookup dupset
                for (f=0;!skip&&f<o_n;++f) {
                    if (avgns[f]>=11) {
                        char t;
                        if (fq[f].seq.a > duplen) {
                            // truncate if needed
                            t = fq[f].seq.s[duplen];
                            fq[f].seq.s[duplen] = '\0';
                        }
                        lookup_it = dupset.find(fq[f].seq.s);
                        if (lookup_it != dupset.end()) {
                            skip=1;                 // 1==dup
                        } else {
                            if (dupset.size() < dupmax) {
                                dupset[fq[f].seq.s]=1;
                            }
                        }
                        if (fq[f].seq.a > duplen) {
                            // restore full length
                            fq[f].seq.s[duplen] = t;
                        }
                    }
                }
            }
            if (!skip) {
               for (f=0;f<o_n;++f) {
                    io_ok=io_ok&&(fputs(fq[f].id.s,fout[f])>=0);
                    io_ok=io_ok&&(fputs(fq[f].seq.s,fout[f])>=0);
                    io_ok=io_ok&&(fputc('\n',fout[f])>=0);
                    io_ok=io_ok&&(fputs(fq[f].com.s,fout[f])>=0);
                    io_ok=io_ok&&(fputs(fq[f].qual.s,fout[f])>=0);
                    io_ok=io_ok&&(fputc('\n',fout[f])>=0);
                }
	       wrec++;
            } else {
                if (skipb) saveskip(fskip, i_n, fq);
                if (skip==2) ++nfiltered;
                if (skip==1) ++dupskip;
            }
		} else {
			if (skipb) saveskip(fskip, i_n, fq);
            if (hompol_skip) {
    			++ntoohompol;
            } else if (lowcom_skip) {
    			++ntoolowcom;
            } else {
    			++ntooshort;
            }
		}
	}

	for (i=0;i<i_n;++i) {
		if (fout[i])  { io_ok = io_ok && ( gzout[i] ? !pclose(fout[i]) : !fclose(fout[i]) ); }
        fin[i].close();
		if (fskip[i]) { if (gzskip[i]) pclose(fskip[i]); else fclose(fskip[i]); }
	}

    if (!io_ok) {
	    fprintf(fstat, "Error during file close, possible partial write, failing\n");
    }

	fprintf(fstat, "Total reads: %d\n", nrec);
	fprintf(fstat, "Too short after clip: %d\n", ntooshort);
    if (nfiltered)
	fprintf(fstat, "Filtered on quality: %d\n", nfiltered);
    if (dupskip)
	fprintf(fstat, "Filtered on duplicates: %d\n", dupskip);
    if (ntoohompol)
	fprintf(fstat, "Filtered on hompolymer: %d\n", ntoohompol);
    if (ntoolowcom)
	fprintf(fstat, "Filtered on low complexity: %d\n", ntoolowcom);
    if (stat_lowcom_b4_total > 0) {
    	fprintf(fstat, "Mean lowcom score: %2.2f(%2.2f), %2.2f(%2.2f) after\n", 
                100*(stat_lowcom_b4_total/(double)stat_lowcom_b4_cnt), 100*stdev(stat_lowcom_b4_cnt,stat_lowcom_b4_total,stat_lowcom_b4_ssq), 
                100*(stat_lowcom_total/(double)stat_lowcom_cnt), 100*stdev(stat_lowcom_cnt,stat_lowcom_total,stat_lowcom_ssq)
        );
    }

	int f;
	if (i_n == 1) {
		f=0;
		for (e=0;e<2;++e) {
			if (cnttrim[f][e]>0) {
				fprintf(fstat, "Clipped '%s' reads: Count: %d, Mean: %.2f, Sd: %.2f\n", e==0?"start":"end", cnttrim[f][e], (double) tottrim[f][e] / cnttrim[f][e], stdev(cnttrim[f][e], tottrim[f][e], ssqtrim[f][e]));
			}
		}
		if (trimql[f] > 0) {
			fprintf(fstat, "Trimmed %d reads by an average of %.2f bases on quality < %d\n", trimql[f], (float) trimqb[f]/trimql[f], qthr);
		}
	} else
		for (f=0;f<i_n;++f) {
			for (e=0;e<2;++e) {
				if (cnttrim[f][e]>0) {
					fprintf(fstat, "Clipped '%s' reads (%s): Count %d, Mean: %.2f, Sd: %.2f\n", e==0?"start":"end", ifil[f], cnttrim[f][e], (double) tottrim[f][e] / cnttrim[f][e], stdev(cnttrim[f][e], tottrim[f][e], ssqtrim[f][e]));
				}
			}
			if (trimql[f] > 0) {
				fprintf(fstat, "Trimmed %d reads (%s) by an average of %.2f bases on quality < %d\n", trimql[f], ifil[f], (float) trimqb[f]/trimql[f], qthr);
			}
		}
	if (nilv3pf > 0) {
		fprintf(fstat, "Filtered %d reads on purity flag\n", nilv3pf);
	}
	if (skipunclip > 0) {
		fprintf(fstat, "Skipped %d unclipped reads\n", skipunclip);
	}
	if (nerr > 0) {
		fprintf(fstat, "Errors (%s): %d\n", ifil[f], nerr);
		return 2;
	}
    if (!io_ok) {
        return 3;
    }
	return 0;
}

int read_fa(FILE *in, int rno, struct ad *fa) {
// note: this only reads one line of sequence!
	fa->nid = getline(&fa->id, &fa->naid, in);
	fa->nseq = getline(&fa->seq, &fa->naseq, in);
	if (fa->nseq <= 0)
		return 0;
	if (fa->id[0] != '>') {
		fprintf(stderr, "Malformed adapter fasta record at line %d\n", rno*2+1);
		return -1;
	}
	// chomp
	fa->seq[--fa->nseq] = '\0';
	fa->id[--fa->nid] = '\0';
	char *p = fa->id+1;
	while (*p == ' ') {
		++p;
	}
	memmove(fa->id, p, strlen(p)+1);
	fa->nid=strlen(fa->id);

	// rna 2 dna
	int i;
	for (i=0;i<fa->nseq;++i) {
		if (fa->seq[i]=='U') fa->seq[i] = 'T';
	}
	return 1;
}

void usage(FILE *f, const char *msg) {
	if(msg)
		fprintf(f, "%s\n", msg);

	fprintf(f, 
"Usage: fastq-mcf [options] <adapters.fa> <reads.fq> [mates1.fq ...] \n"
"Version: %s\n"
"\n"
"Detects levels of adapter presence, computes likelihoods and\n"
"locations (start, end) of the adapters.   Removes the adapter\n"
"sequences from the fastq file(s).\n"
"\n"
"Stats go to stderr, unless -o is specified.\n"
"\n"
"Specify -0 to turn off all default settings\n"
"\n"
"If you specify multiple 'paired-end' inputs, then a -o option is\n" 
"required for each.  IE: -o read1.clip.q -o read2.clip.fq\n"
"\n"
"Options:\n"
"    -h       This help\n"
"    -o FIL   Output file (stats to stdout)\n"
"    -O N     Only output the first N records (all)\n"
"    -s N.N   Log scale for adapter minimum-length-match (2.2)\n"
"    -t N     %% occurance threshold before adapter clipping (0.25)\n"
"    -m N     Minimum clip length, overrides scaled auto (1)\n"
"    -p N     Maximum adapter difference percentage (10)\n"
"    -l N     Minimum remaining sequence length (19)\n"
"    -L N     Maximum remaining sequence length (none)\n"
"    -D N     Remove duplicate reads : Read_1 has an identical N bases (0)\n"
"    -k N     sKew percentage-less-than causing cycle removal (2)\n"
"    -x N     'N' (Bad read) percentage causing cycle removal (20)\n"
"    -q N     quality threshold causing base removal (10)\n"
"    -w N     window-size for quality trimming (1)\n"
"    -H       remove >95%% homopolymer reads (no)\n"
"    -X       remove low complexity reads (no)\n"
//"    -F FIL  remove sequences that align to FIL\n"
"    -0       Set all default parameters to zero/do nothing\n"
"    -U|u     Force disable/enable Illumina PF filtering (auto)\n"
"    -P N     Phred-scale (auto)\n"
"    -R       Don't remove N's from the fronts/ends of reads\n"
"    -n       Don't clip, just output what would be done\n"
"    -K       Only keep clipped reads\n"
"    -S       Save all discarded reads to '.skip' files\n"
"    -C N     Number of reads to use for subsampling (300k)\n"
"    -d       Output lots of random debugging stuff\n"
"\n"
"Minimum trimming options:\n"
"    --min-start-trim    NUM       Always trim at least NUM bases from start\n"
"    --min-end-trim      NUM       Always trim at least NUM bases from end\n"
"\n"
"Quality adjustment options:\n"
"    --cycle-adjust      CYC,AMT   Adjust cycle CYC (negative = offset from end) by amount AMT\n"
"    --phred-adjust      SCORE,AMT Adjust score SCORE by amount AMT\n"
"    --phred-adjust-max  SCORE     Adjust scores > SCORE to SCOTE\n"
"\n"
"Filtering options*:\n"
"    --[mate-]qual-mean  NUM       Minimum mean quality score\n"
"    --[mate-]qual-gt    NUM,THR   At least NUM quals > THR\n" 
"    --[mate-]max-ns     NUM       Maxmium N-calls in a read (can be a %%)\n"
"    --[mate-]min-len    NUM       Minimum remaining length (same as -l)\n"
"    --homopolymer-pct   PCT       Homopolymer filter percent (95)\n"
"    --lowcomplex-pct    PCT       Complexity filter percent (95)\n"
"    --keep-clipped                Only keep clipped (same as -K)\n"
"    --max-output-reads   N        Only output first N records (same as -O)\n"
"\n"
"If mate- prefix is used, then applies to second non-barcode read only\n"
/*
"Config:\n"
"\n"
"Some options are best set globally, such as the aligner to use\n"
"for filtering, these can be ENV vars, or in /etc/ea-utils.conf:\n"
"\n"
"Command line options, if specified, always override config vars.\n"
"\n"
"When uses as environment vars, they are all caps, and with\n"
"EAUTILS_ as the prefix (IE: EAUTILS_PHRED=33)\n"
"\n"
"    phred           (auto)\n"
"    sample_reads    100000\n"
"    scale_clip_len  2.2\n"
"    trim_skew       2\n"
"    trim_quality    10\n"
"    min_clip_len    0\n"
"    min_seq_remain  15\n"
"    max_adap_diff   20\n"
//"    cmd_align_se    bowtie -S %%i -f %%1\n"
//"    cmd_align_pe    bowtie -S %%i -1 %%1 -2 %%2\n"
//"\n"
//"Command lines must return SAM formatted lines, %%i is the filter FIL,\n"
//"%%1 and %%2 are the first and second fastq's\n"
*/
"\n"
"Adapter files are 'fasta' formatted:\n"
"\n"
"Specify n/a to turn off adapter clipping, and just use filters\n"
"\n"
"Increasing the scale makes recognition-lengths longer, a scale\n"
"of 100 will force full-length recognition of adapters.\n"
"\n"
"Adapter sequences with _5p in their label will match 'end's,\n"
"and sequences with _3p in their label will match 'start's,\n"
"otherwise the 'end' is auto-determined.\n"
"\n"
"Skew is when one cycle is poor, 'skewed' toward a particular base.\n"
"If any nucleotide is less than the skew percentage, then the\n"
"whole cycle is removed.  Disable for methyl-seq, etc.\n"
"\n"
"Set the skew (-k) or N-pct (-x) to 0 to turn it off (should be done\n"
"for miRNA, amplicon and other low-complexity situations!)\n"
"\n"
"Duplicate read filtering is appropriate for assembly tasks, and\n"
"never when read length < expected coverage.  -D 50 will use\n"
"4.5GB RAM on 100m DNA reads - be careful. Great for RNA assembly.\n"
"\n"
"*Quality filters are evaluated after clipping/trimming\n"
"\n"
"Homopolymer filtering is a subset of low-complexity, but will not\n"
"be separately tracked unless both are turned on.\n"
	,VERSION);
}

inline int char2bp(char c) {
        if (c == 'A' || c == 'a') return B_A;
        if (c == 'C' || c == 'c') return B_C;
        if (c == 'G' || c == 'g') return B_G;
        if (c == 'T' || c == 't') return B_T;
        return B_N;
}

inline int char2bp_rc(char c) {
        if (c == 'A' || c == 'a') return B_T;
        if (c == 'C' || c == 'c') return B_G;
        if (c == 'G' || c == 'g') return B_C;
        if (c == 'T' || c == 't') return B_A;
        return B_N;
}


inline char bp2char(int b) {
        if (b == B_A) return 'A';
        if (b == B_C) return 'C';
        if (b == B_G) return 'G';
        if (b == B_T) return 'T';
        return 'N';
}

void saveskip(FILE **fout, int fo_n, struct fq *fq)  {
	int f;
	for (f=0;f<fo_n;++f) {
		fputs(fq[f].id.s,fout[f]);
		fputs(fq[f].seq.s,fout[f]);
		fputc('\n',fout[f]);
		fputs(fq[f].com.s,fout[f]);
		fputs(fq[f].qual.s,fout[f]);
		fputc('\n',fout[f]);
	}
}

int meanqwin(const char *q, int qn, int i, int w) {
	if (w > qn) w=qn/4;                         // maximum window is length/4
	int s = i-w/2;                              // start/end window
	int e = i+w/2;
	if (s < 0) {e-=s;s=0;}                      // shift window over if you're past the start
	if (e >= qn) {s-=((e-qn)+1);e=qn-1;}        // shift window over if you're past the end
	int t = 0;
	for (i=s;i<=e;++i) {
		t+=q[i];	
	}
	return t / (e-s+1);                         // mean quality within the window at that position
}

bool evalqual(struct fq &fq, int file_num) {
    int t_mean, t_max_ns, t_xgt_num, t_xgt_min, t_max_n_pct;


    if (file_num <= 0) {
        // applies to file 1
        t_mean=qf_mean;
        t_max_ns=qf_max_ns;
        t_max_n_pct=qf_max_n_pct;
        t_xgt_num=qf_xgt_num;
        t_xgt_min=qf_xgt_min;
    } else {
        // applies to file 2 or greater, only if they are set
        t_mean=qf2_mean > 0 ? qf2_mean : qf_mean;
        t_max_ns=qf_max_ns > -1 ? qf2_max_ns : qf_max_ns;
        t_max_n_pct=qf_max_n_pct > -1 ? qf2_max_n_pct : qf_max_n_pct;
        t_xgt_num=qf_xgt_num > 0 ? qf2_xgt_num : qf_xgt_num;
        t_xgt_min=qf_xgt_min > 0 ? qf2_xgt_min : qf_xgt_min;
    }

    if (t_max_n_pct>=0) {
        t_max_ns=max(t_max_ns,(fq.qual.n*100)/t_max_n_pct);
    }

    if (t_mean > 0) {
        int t = 0;
        int i;
        for (i=0;i<=fq.qual.n;++i) {
            t+=fq.qual.s[i];
        }
        if ((t/fq.qual.n-phred) < t_mean) {
            return false;
        }
    }
    if (t_max_ns >= 0) {
        int t = 0;
        int i;
        for (i=0;i<=fq.seq.n;++i) {
            t+=(fq.seq.s[i]=='N');
        }

//        if (debug > 2) fprintf(stderr,"maxn: max:%d,t:%d,i:%d,id:%s", t_max_ns, t, i, fq.id.s);

        if (t > t_max_ns) {
            return false;
        }
    }

    if (t_xgt_num > 0) {
        int t = 0;
        int i;
        int h = t_xgt_min+phred;
        for (i=0;i<=fq.qual.n;++i) {
            t+=(fq.qual.s[i]>=h);
        }
        if (t < t_xgt_num) {
            return false;
        }
    }
    return true;
}

void valid_arg(char opt, const char *arg) {
    if (!arg || !*arg || *arg == '-') {
        fprintf(stderr,"Option '%c' requires an argument.\n\n", opt);
        usage(stderr); 
        exit(1);
    }
}

bool  arg_int_pair(const char *optarg, int &a, int&b) {
    if (!strchr(optarg, ',')) {
        return false;
    }
    a=atoi(optarg);
    b=atoi(strchr(optarg, ',')+1);
}

char *arg2cmd(int argc, char** argv) {
    char *buf=NULL;
    int n = 0;
    int k, i;
    for (i=1; i <argc;++i) {
        int k=strlen(argv[i]);
        buf=( char *)realloc(buf,n+k+4);
        char *p=buf+n;
        char endq=0;
        // this is a poor mans quoting, which is good enough for anything that's not rediculous
        if (strchr(argv[i], ' ')) {
            if (!strchr(argv[i], '\'')) {
                *p++='\'';
                endq='\'';
            } else {
                *p++='\"';
                endq='\"';
            }
        }
        memcpy(p, argv[i], k);
        p+=k;
        if (i < (argc-1)) *p++=' ';
        if (endq) *p++=endq;
        *p='\0';
        n = p-buf;
    }
    return buf;
}

std::string arg2cmdstr(int argc, char **argv) {
    char *tmp=arg2cmd(argc, argv);
    std::string ret=tmp;
    free(tmp);
    return ret;
}


/* vim: set noai ts=4 sw=4: */
