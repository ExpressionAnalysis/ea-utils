#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <sys/stat.h>

/*

Currently only works with adapter sequences that are at the END of read lines.
See usage below.

*/

#define MAX_ADAPTER_NUM 1000
#define max(a,b) (a>b?a:b)
#define min(a,b) (a<b?a:b)
#define SCANLEN 15
#define SCANMIDP ((int) SCANLEN/2)
#define MAX_FILES 3
#define meminit(l) (memset(&l,0,sizeof(l)))
#define stdev(cnt, sum, ssq) sqrt((((double)cnt)*ssq-pow((double)sum,2)) / ((double)cnt*((double)cnt-1)))
#define B_A     0
#define B_C     1
#define B_G     2
#define B_T     3
#define B_N     4
#define B_CNT   5

struct fq {
	char *id;   int nid;   size_t naid;
	char *seq;  int nseq;  size_t naseq;
	char *com;  int ncom;  size_t nacom;
	char *qual; int nqual; size_t naqual;
};

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
int read_fq(FILE *in, int rno, struct fq *fq);		// 0=done, 1=ok, -1=err+continue

int char2bp(char c);
char bp2char(int b);

void usage(FILE *f, char *msg=NULL);
int hd(char *a, char *b, int n);
int debug=0;
int main (int argc, char **argv) {
	char c;
	bool eol;
	int nmin = 1, nkeep = 15;
	float minpct = 0.25;
	int pctdiff = 20;
	int sampcnt = 40000;				// # of reads to sample to determine adapter profile, and base skewing
	int xmax = -1;
	float scale = 2.2;
	int noclip=0;
	char end[MAX_FILES]; meminit(end);
	float skewpct = 2; 			// any base at any position is less than skewpct of reads
	float pctns = 5;				// any base that is more than 10% n's
	bool rmns = 1;				// any base that is more than 10% n's
	int qthr = 10;				// remove end of-read with quality < qthr
	char phred = 64;
	bool force = 0;

	int i;
	
	char *afil = NULL;
	char *ifil[MAX_FILES]; meminit(ifil);
	char *ofil[MAX_FILES]; meminit(ofil);
	int i_n = 0;
	int o_n = 0;
	int e_n = 0;

	while (	(c = getopt (argc, argv, "-nfRdbehp:o:l:s:m:t:k:x:P:q:")) != -1) {
		switch (c) {
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
		case 'f': force = true; break;
		case 'k': skewpct = atof(optarg); break;
		case 'q': qthr = atoi(optarg); break;
		case 'x': pctns = atof(optarg); break;
		case 'R': rmns = false; break;
		case 'p': pctdiff = atoi(optarg); break;
		case 'P': phred = (char) atoi(optarg); break;
		case 'h': usage(stdout); return 1; 
		case 'o': if (!o_n < MAX_FILES) 
				ofil[o_n++] = optarg;
			  break;
		case 's': scale = atof(optarg); break;
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

        FILE *ain = fopen(afil, "r");
        if (!ain) {
                fprintf(stderr, "Error opening file '%s': %s\n",afil, strerror(errno));
                return 1;
        }

	FILE *fstat = stderr;
	if (strcmp(ofil[0], "-") && !noclip) {
		fstat = stdout;
	}
	if (noclip) {
		fstat = stdout;
	}

	FILE *fin[MAX_FILES];  meminit(fin);
	FILE *fout[MAX_FILES]; meminit(fout);

	for (i=0;i<i_n;++i) {
		if (!(fin[i]=fopen(ifil[i], "r"))) {
			fprintf(stderr, "Error opening file '%s': %s\n",ifil[i], strerror(errno));
			return 1;
		}
	}


	struct ad ad[MAX_ADAPTER_NUM+1];
	memset(ad, 0, sizeof(ad));

	int acnt=0, ok=0, rno=0;	// adapter count, ok flag, record number
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

	if (acnt == 0) 
		exit(1);

	fprintf(fstat, "Scale used: %g\n", scale);

	int maxns = 0;						// max sequence length
	int avgns[MAX_FILES]; meminit(avgns);			// average sequence length per file
	// read length
        for (i=0;i<i_n;++i) {
            struct stat st;
            stat(ifil[i], &st);
            fseek(fin[i], st.st_size > sampcnt ? (st.st_size-sampcnt)/3 : 0, 0);

            char *s = NULL; size_t na = 0; int nr = 0, ns = 0, totn[MAX_FILES]; meminit(totn);
            while (getline(&s, &na, fin[i]) > 0) {
                if (*s == '@')  {
                        if ((ns=getline(&s, &na, fin[i])) <=0)
                                break;
                        --ns;                                   // don't count newline for read len
                        ++nr;
	                avgns[i] += ns;
			if (ns > maxns) maxns = ns;

			if (nr >= 1000)	
				break;
		}
	    }

	    if (nr)
		    avgns[i] = avgns[i]/nr;
	}

	for (i=0;i<i_n;++i) {
		if (avgns[i] == 0) {
			fprintf(stderr, "No records in file %s\n", ifil[i]);
			exit(1);
		}
	}

	if (debug) printf("Max ns: %d, Avg[0]: %d\n", maxns, avgns[0]);

	// total base count per read position in sample
	int bcnt[MAX_FILES][2][maxns][6]; meminit(bcnt);
	int nr;
	for (i=0;i<i_n;++i) {

	    struct stat st;
	    stat(ifil[i], &st);
	    fseek(fin[i], st.st_size > sampcnt ? (st.st_size-sampcnt)/3 : 0, 0);

	    char *s = NULL; size_t na = 0; int nr = 0, ns = 0;
            while (getline(&s, &na, fin[i]) > 0) {
		if (*s == '@')  {
			if ((ns=getline(&s, &na, fin[i])) <=0) 
				break;
			--ns;					// don't count newline for read len
			++nr;

			if (avgns[i] < 11) 			// reads of avg length < 11 ? barcode lane, skip it
				continue;

			// to be safe, we don't assume reads are fixed-length, not any slower, just a little more code
			int b;
			for (b = 0; b < ns/2 && b < maxns; ++b) {
				++bcnt[i][0][b][char2bp(s[b])];		// count from begin
				++bcnt[i][0][b][B_CNT];			// count of samples at position
				++bcnt[i][1][b][char2bp(s[ns-b-1])];	// count from end
				++bcnt[i][1][b][B_CNT];			// count of samples at offset-from-end position
			}
			
			int a;
			char buf[SCANLEN+1];
			strncpy(buf, s, SCANLEN);
			for(a=0;a<acnt;++a) {
				char *p;
				if (p = strstr(s, ad[a].escan)) { 
			//		if (debug) fprintf(stderr, "END  : A: %s (%s), P: %d, SL: %d, Z:%d\n", ad[a].id, ad[a].escan, p-s, ns, (p-s) == ns-SCANLEN);
					if ((p-s) == ns-SCANLEN) 
						++ad[a].ecntz[i];
					++ad[a].ecnt[i];
				}
					
				if ((p = strstr(ad[a].seq, buf))) { 
			//		if (debug) fprintf(stderr, "BEGIN: A: %s (%s), P: %d, SL: %d, Z:%d\n", ad[a].id, ad[a].seq, p-ad[a].seq, ns, (p-ad[a].seq )  == ad[a].nseq-SCANLEN);
					if (p-ad[a].seq == ad[a].nseq-SCANLEN) 
						++ad[a].bcntz[i];
					++ad[a].bcnt[i];
				}
			}
		}
		if (nr >= sampcnt)		// enough samples 
			break;
            }
	    if (nr < sampcnt)			// fewer than max, set for thresholds
		sampcnt=nr;
	}

	// look for severe base skew, and auto-trim ends based on it
	int sktrim[i_n][2]; meminit(sktrim);
	for (i=0;i<i_n;++i) {
		if (avgns[i] < 11) 			// reads of avg length < 11 ? barcode lane, skip it
			continue;
		int e;
		for (e = 0; e < 2; ++e) {
			int p;
			for (p = 0; p < maxns/2; ++p) {
				int b;

				int skth = (int) ( (float) bcnt[i][e][p][B_CNT] * ( skewpct / 100.0 ) ) ;	// skew threshold
				int thr_n = (int) ( (float) bcnt[i][e][p][B_CNT] * ( pctns / 100.0 ) );		// n-threshold
	
				if (debug > 1) 
					printf("Sk Prof [%d, %d]: skth=%d, bcnt=%d, ncnt=%d, a=%d, c=%d, g=%d, t=%d\n", e, p, skth, 
						bcnt[i][e][p][B_CNT], bcnt[i][e][p][B_N], bcnt[i][e][p][B_A], 
						bcnt[i][e][p][B_C], bcnt[i][e][p][B_G], bcnt[i][e][p][B_T]);

				if (skth < 10)						// too few samples to detect skew
					continue;

				int tr = 0;
				for (b = 0; b < 4; ++b) {
					if (bcnt[i][e][p][b] < skth) {			// too few bases of this type
						tr=1;
						break;
					}
				}
				if (bcnt[i][e][p][B_N] > thr_n) 			// too many n's
					tr=1;

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

	int e;
	bool someskew = false;
	for (i=0;i<i_n;++i) {
		int totskew = sktrim[i][0] + sktrim[i][1];
		if (maxns - totskew < nkeep) {
			fprintf(fstat, "Warning: Too much skewing found, disabling skew clipping\n");
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
		    if (ad[a].ecnt[i] > athr || ad[a].bcnt[i] > athr) {
			int cnt;
			if (debug) fprintf(stderr, "EC: %d, BC:%d, ECZ: %d, BCZ: %d\n", ad[a].ecnt[i], ad[a].bcnt[i], ad[a].ecntz[i], ad[a].bcntz[i]);
			// heavily weighted toward start/end maches
			if ((ad[a].ecnt[i] + 10*ad[a].ecntz[i]) >= (ad[a].bcnt[i] + 10*ad[a].bcntz[i])) {
				ad[a].end[i]='e';
				cnt = ad[a].ecnt[i];
			} else {
				ad[a].end[i]='b';
				cnt = ad[a].bcnt[i];
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
				fprintf(fstat, ", warning end was not reliable\n", ad[a].id, ad[a].seq);
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

	if (acnt == 0 && !someskew && !force) {
		fprintf(fstat, "No adapters found and no skewing detected\n");
		if (noclip) exit (1);			// for including in a test
		exit(0);				// not really an error, check size of output files
	}

	if (noclip)
		exit(0);

	for (i=0;i<o_n;++i) {
		if (!strcmp(ofil[i],"-")) {
			fout[i]=stdout;
		} else if (!(fout[i]=fopen(ofil[i], "w"))) {
                        fprintf(stderr, "Error opening output file '%s': %s\n",ofil[i], strerror(errno));
                        return 1;
                }
	}

	struct fq fq[MAX_FILES];	
        memset(&fq, 0, sizeof(fq));

	int nrec=0;
	int nerr=0;
	int nok=0;
	int ntooshort=0;
	// total per read
	int ntrim[MAX_FILES]; meminit(ntrim);

	// total per end
	int cnttrim[MAX_FILES][2]; meminit(cnttrim);
	double tottrim[MAX_FILES][2]; meminit(tottrim);
	double ssqtrim[MAX_FILES][2]; meminit(ssqtrim);
	int trimql[MAX_FILES]; meminit(trimql);
	int trimqb[MAX_FILES]; meminit(trimqb);

	int nbtrim=0;
	int read_ok;

	if (i_n > 0)
		fprintf(fstat, "Files: %d\n", i_n);

	for (i=0;i<i_n;++i)
		fseek(fin[i], 0, 0);

	while (read_ok=read_fq(fin[0], nrec, &fq[0])) {
		for (i=1;i<i_n;++i) {
			int mok=read_fq(fin[i], nrec, &fq[i]);
			if (mok != read_ok) {
				fprintf(stderr, "# of rows in mate file '%s' doesn't match, quitting!\n", ifil[i]);
				return 1;
			}
		}
		++nrec;
		if (read_ok < 0) {
			++nerr;
			continue;
		}

		// chomp

		int dotrim[MAX_FILES][2];
		bool skip = 0;							// skip whole record?
		int f;	
		for (f=0;f<i_n;++f) {
		    dotrim[f][0] = sktrim[f][0];					// default, trim to detected skew levels
		    dotrim[f][1] = sktrim[f][1];

		    if (rmns) {
			    for (i=dotrim[f][0];i<fq[f].nseq/3;++i) {
			    		// trim N's from the front
					if (fq[f].seq[i] == 'N') 
						dotrim[f][0] = i + 1;
			    }
			    for (i=dotrim[f][1];i<fq[f].nseq/3;++i) {
			    		// trim N's from the end
					if (fq[f].seq[fq[f].nseq-i-1] == 'N')
						dotrim[f][1] = i + 1;
			    }
		    }

                    if (qthr > 0) {
			    bool istrimq = false;
                            for (i=dotrim[f][1];i<fq[f].nseq/2;++i) {
                                        // trim N's from the end
                                        if ((fq[f].qual[fq[f].nseq-i-1]-phred) < qthr) {
						++trimqb[f];
						istrimq = true;
                                                dotrim[f][1] = i + 1;
					}
                            }
			    if (istrimq) ++trimql[f];
                    }

		    int bestscore_e = INT_MAX, bestoff_e = 0, bestlen_e = 0; 
		    int bestscore_b = INT_MAX, bestoff_b = 0, bestlen_b = 0; 

		    for (i =0; i < acnt; ++i) {
		    	if (debug) fprintf(stderr, "seq[%d]: %s %d\n", f, fq[f].seq, fq[f].nseq);

			if (!ad[i].end[f])
				continue;

			int nmatch = ad[i].thr[f];
			if (!nmatch) nmatch = ad[i].nseq;			// full match required if nmin == 0
	
			// how far in to search for a match?
			int mx = ad[i].nseq;
			if (xmax) {
				 mx = fq[f].nseq;
				 if (xmax > 0 && (xmax+ad[i].nseq) < mx)
					mx = xmax+ad[i].nseq;			// xmax is added to adapter length
			}

			if (debug)
				fprintf(stderr, "adapter: %s, adlen: %d, nmatch: %d, mx: %d\n", ad[i].seq, ad[i].nseq, nmatch, mx);

			if (ad[i].end[f] == 'e') {
			    int off;
			    for (off = nmatch; off <= mx; ++off) {		// off is distance from tail of sequence
				char *seqtail = fq[f].seq+fq[f].nseq-off; 	// search at tail
				int ncmp = off<ad[i].nseq ? off : ad[i].nseq;
				int mind = (pctdiff * ncmp) / 100;
				int d = hd(ad[i].seq,seqtail,ncmp);		// # differences
				if (debug>1)
					fprintf(stderr, "tail: %s, bestoff: %d, off: %d, ncmp: %d, mind: %d, hd %d\n", seqtail, bestoff_e, off, ncmp, mind, d);
				if (d <= mind) {
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
                                char *seqstart = fq[f].seq+off-ncmp;		// offset into sequence (if any)
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
		    // lengthen trim based on best level
		    if (bestoff_b > dotrim[f][0])
			dotrim[f][0]=bestoff_b;

		    if (bestoff_e > dotrim[f][1])
			dotrim[f][1]=bestoff_e;

		    int totclip = dotrim[f][0] + dotrim[f][1];

		    if (totclip > 0) {
		    	if ( (fq[f].nseq-totclip) < nkeep) {
					// skip all reads if one is severely truncated ??
					// maybe not... ?
					skip = 1;
					break;
			}

			// count number of adapters clipped, not the number of rows trimmed
			if (bestoff_b > 0 || bestoff_e > 0) 
				++ntrim[f];

			// save some stats
			if (bestoff_b > 0) {
				cnttrim[f][0]++;
				tottrim[f][0]+=bestoff_b;
				ssqtrim[f][0]+=bestoff_b*bestoff_b;
			}
			if (bestoff_e > 0) {
				cnttrim[f][1]++;
				tottrim[f][1]+=bestoff_e;
				ssqtrim[f][1]+=bestoff_e*bestoff_e;
			}

		    }
		}

		if (!skip) {
			int f;
			for (f=0;f<o_n;++f) {
                                if (dotrim[f][0] > 0) {
                                        if (debug) printf("trimming %d from begin\n", dotrim[f][0]);
                                        memmove(fq[f].seq ,fq[f].seq +dotrim[f][0],fq[f].nseq -=dotrim[f][0]);
                                        memmove(fq[f].qual,fq[f].qual+dotrim[f][0],fq[f].nqual-=dotrim[f][0]);
                                        fq[f].seq[fq[f].nseq]='\0';
                                        fq[f].qual[fq[f].nqual]='\0';
                                }
				if (dotrim[f][1] > 0) {
					if (debug) printf("trimming %d from end\n", dotrim[f][1]);
					fq[f].seq [fq[f].nseq -dotrim[f][1]]='\0';
					fq[f].qual[fq[f].nqual-dotrim[f][1]]='\0';
				}
				fputs(fq[f].id,fout[f]);
				fputs(fq[f].seq,fout[f]);
				fputc('\n',fout[f]);
				fputs(fq[f].com,fout[f]);
				fputs(fq[f].qual,fout[f]);
				fputc('\n',fout[f]);
			}
		} else 
			++ntooshort;
	}
	fprintf(fstat, "Total reads: %d\n", nrec);
	fprintf(fstat, "Too short after clip: %d\n", ntooshort);

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
	if (nerr > 0) {
		fprintf(fstat, "Errors (%s): %d\n", ifil[f], nerr);
		return 2;
	}
	return 0;
}

// returns number of differences
inline int hd(char *a, char *b, int n) {
	int d=0;
	//if (debug) fprintf(stderr, "hd: %s,%s ", a, b);
	while (*a && *b && n > 0) {
		if (*a != *b) ++d;
		--n;
		++a;
		++b;
	}
	//if (debug) fprintf(stderr, ", %d/%d\n", d, n);
	return d+n;
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

int read_fq(FILE *in, int rno, struct fq *fq) {
        fq->nid = getline(&fq->id, &fq->naid, in);
        fq->nseq = getline(&fq->seq, &fq->naseq, in);
        fq->ncom = getline(&fq->com, &fq->nacom, in);
        fq->nqual = getline(&fq->qual, &fq->naqual, in);
        if (fq->nqual <= 0)
                return 0;
        if (fq->id[0] != '@' || fq->com[0] != '+' || fq->nseq != fq->nqual) {
                fprintf(stderr, "Malformed fastq record at line %d\n", rno*2+1);
                return -1;
        }

	// chomp
	fq->seq[--fq->nseq] = '\0';
	fq->qual[--fq->nqual] = '\0';
        return 1;
}

void usage(FILE *f, char *msg) {
	if(msg)
		fprintf(f, "%s\n", msg);

	fprintf(f, 
"usage: fastq-mcf [options] <adapters.fa> <reads.fq> [mates1.fq ...] \n"
"\n"
"Detects levels of adapter presence, computes likelihoods and\n"
"locations (start, end) of the adapters.   Removes the adapter\n"
"sequences from the fastq file, and sends it to stdout.\n\n"
"Stats go to stderr, unless -o is specified.\n"
"\n"
"If you specify multiple 'paired-end' inputs, then a -o option is\n" 
"required for each.  IE: -o read1.clip.q -o read2.clip.fq\n"
"\n"
"Options:\n"
"	-h	This help\n"
"	-o FIL	Output file (stats to stdout)\n"
"	-s N.N	Log scale for clip pct to threshold (2.5)\n"
"	-t N	%% occurance threshold before clipping (0.25)\n"
"	-m N	Minimum clip length, overrides scaled auto (1)\n"
"	-p N	Maximum adapter difference percentage (20)\n"
"	-l N	Minimum remaining sequence length (15)\n"
"	-k N	sKew percentage causing trimming (2)\n"
"	-q N	quality threshold causing trimming (10)\n"
"	-f	force output, even if not much will be done\n"
"	-P N	phred-scale (64)\n"
"	-x N	'N' (Bad read) percentage causing trimming (10)\n"
"	-R      Don't remove N's from the fronts/ends of reads\n"
"	-n	Don't clip, just output what would be done\ni"
"\n"
"Increasing the scale makes recognition-lengths longer, a scale\n"
"of 100 will force full-length recognition.\n"
"\n"
"Set the skew (-k) or N-pct (-x) to 100 or 0 to turn it off\n"
	);
}

inline int char2bp(char c) {
        if (c == 'A' || c == 'a') return B_A;
        if (c == 'C' || c == 'c') return B_C;
        if (c == 'G' || c == 'g') return B_G;
        if (c == 'T' || c == 't') return B_T;
        return B_N;
}

inline char bp2char(int b) {
        if (b == B_A) return 'A';
        if (b == B_C) return 'C';
        if (b == B_G) return 'G';
        if (b == B_T) return 'T';
        return 'N';
}

