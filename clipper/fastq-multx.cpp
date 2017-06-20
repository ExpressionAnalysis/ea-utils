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

#include "fastq-lib.h"

#define MAX_BARCODE_NUM 6000
#define MAX_GROUP_NUM 500
// factor to divide max by
#define THFIXFACTOR 20
#define endstr(e) (e=='e'?"end":e=='b'?"start":"n/a")

const char * VERSION = "1.03";

// barcode
struct bc {
	line id;
	line seq;
	char *out[6];			// one output per input
	FILE *fout[6];
	bool gzout[6];
	uint64_t cnt;			// count found
	bool shifted;			// count found in 1-shifted position
	char * dual;			// is this a dual-indexed barcode?  if so, this points to the second index.
	int dual_n;			// length of dual
};

// group of barcodes
struct group {
	char *id;
	int tcnt;			// number of codes past thresh
	int i;				// my index
};

// barcode group
struct bcg {
	struct bc b;			// barcode
        line group;			// group (fluidigm, truseq, etc)
        int bcnt[6];			// matched begin of file n
        int ecnt[6];			// matched end of file n
        int bscnt[6];			// matched begin of file n, shifted by 1
        int escnt[6];			// matched end of file n, shifted by 1
        int dbcnt[6];                   // dual matched begin of file n
        int decnt[6];                   // dual matched end of file n
	struct group *gptr;		
};

struct group* getgroup(char *s);

void usage(FILE *f);
static int debug=0;
// it's times like this when i think a class might be handy, but nah, not worth it
typedef struct bnode {
	char *seq;
	int cnt;
} bnode;

struct group grs[MAX_GROUP_NUM];
static int grcnt=0;

struct bc bc[MAX_BARCODE_NUM+1];
static int bcnt=0;

static int pickmax=0;
static int pickmax2=0;
static void *picktab=NULL;
void pickbest(const void *nodep, const VISIT which, const int depth);
int bnodecomp(const void *a, const void *b) {return strcmp(((bnode*)a)->seq,((bnode*)b)->seq);};
static float pickmaxpct=0.10;
void getbcfromheader(struct fq *fqin, struct fq *bc, char **s2=NULL, int *ns2=NULL);
void getbcfromheader(char *s, int *ns, char **q=NULL, char **s2=NULL, int *ns2=NULL);

int ignore;
size_t ignore_st;


int main (int argc, char **argv) {
	char c;
	bool trim = true;
	int mismatch = 1;
	int distance = 2;
	int poor_distance = 0;       // count of skipped reads on distance only
	int quality = 0;
	char end = '\0';
	char dend = '\0';
	bool dual = false;
	char *in[6];
	const char *out[6];
	int f_n=0;
	int f_oarg=0;
	const char* guide=NULL;		// use an indexed-read
	const char* list=NULL;		// use a barcode master list
	char verify='\0';
	bool noexec = false;
    bool seqnames = false;
	const char *group = NULL;
    bool usefile1 = false;
    int phred = 33;
    double threshfactor = 1;
    int bcinheader = 0;

	int i;
	bool omode = false;	
	char *bfil = NULL;
	while (	(c = getopt (argc, argv, "-DzxnHhbeosv:m:B:g:L:l:G:q:d:t:")) != -1) {
		switch (c) t:{
		case '\1': 
                       	if (omode) {
				if (f_oarg<5)
					out[f_oarg++] = optarg;
				else {
					usage(stderr); return 1;
				}
			} else if (!bfil && !guide && !list) 
				bfil = optarg; 
			else if (f_n<5) {
				in[f_n++] = optarg; 
			} else {
				usage(stderr); return 1;
			}
			break;
                case 'o': omode=true; break;
                case 's': seqnames=true; break;
                case 'v': 
			if (strlen(optarg)>1) {
				fprintf(stderr, "Option -v requires a single character argument");
				exit(1);
			}
			verify = *optarg; break;
		case 'b': end = 'b'; break;
		case 'h': usage(stdout); exit(0); break;
		case 'H': bcinheader = 1; usefile1=1; break;
		case 'e': end = 'e'; break;
		case 'G': group = optarg; break;
		case 'g': 
			guide = optarg;
			in[f_n++] = optarg;
			out[f_oarg++] = "n/a";
			break;
		case 'l': list = optarg; usefile1=0; break;
		case 'L': list = optarg; usefile1=1; break;
		case 'B': bfil = optarg; list = NULL; break;
		case 'x': trim = false; break;
		case 'n': noexec = true; break;
		case 't': threshfactor = atof(optarg); break;
		case 'm': mismatch = atoi(optarg); break;
		case 'd': distance = atoi(optarg); break;
		case 'q': quality = atoi(optarg); break;
		case 'D': ++debug; break;
		case '?': 
		     if (strchr("vmBglG", optopt))
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

	if (group && !list) {
		fprintf(stderr, "Error: -G only works with -l\n");
		return 1;
	}

    if ((list && guide) || (list && bfil) || (guide && bfil)) {
            fprintf(stderr, "Error: Only one of -B -g or -l\n");
            return 1;
    }

	if (f_n != f_oarg) {
		fprintf(stderr, "Error: number of input files (%d) must match number of output files following '-o'.\n", f_n);
		return 1;
	}

	if (argc < 3 || !f_n || (!bfil && !guide && !list)) {
		usage(stderr);
		return 1;
	}

    quality+=phred;

	FILE *fin[6];
	bool gzin[6]; meminit(gzin);
	for (i = 0; i < f_n; ++i) {
		fin[i]=gzopen(in[i],"r",&gzin[i]);
	}

	// set all to null, zero
	meminit(bc);


	// 3 ways to get barcodes
	if (list) {
		// use a list of barcode groups... determine the best set, then use the determined set 
		struct bcg *bcg = (struct bcg *) malloc(sizeof(*bcg) * MAX_GROUP_NUM * MAX_BARCODE_NUM);
		if (!bcg) {
                        fprintf(stderr, "Out of memory\n");
                        return 1;
		}
		memset(bcg, 0, sizeof(*bcg) * MAX_GROUP_NUM * MAX_BARCODE_NUM);
		int bgcnt=0;
		int b;
        FILE *lin = fopen(list, "r");
        if (!lin) {
                fprintf(stderr, "Error opening file '%s': %s\n",list, strerror(errno));
                return 1;
        }
        // read barcode groups
        int ok;

        while (bgcnt < (MAX_GROUP_NUM * MAX_BARCODE_NUM) && (ok = read_line(lin, bcg[bgcnt].b.id))) {
            if (ok <= 0) break;
            if (bcg[bgcnt].b.id.s[0]=='#') continue;
            bcg[bgcnt].b.id.s=strtok(bcg[bgcnt].b.id.s, "\t\n\r ");
            bcg[bgcnt].b.seq.s=strtok(NULL, "\t\n\r ");
            char *g=strtok(NULL, "\n\r");
			if (!g) {
				if (bgcnt==0){
					fprintf(stderr,"Barcode guide list needs to be ID<whitespace>SEQUENCE<whitespace>GROUP");
					return 1;
				} else {
					continue;
				}
			}
			if (group) {
				if (strcasecmp(group, g)) {
					continue;
				}
			}
            if (!strcmp(bcg[bgcnt].b.seq.s,"seq")) continue;

            // dual indexed indicated by a dash in the sequence...
			if (bcg[bgcnt].b.dual=strchr(bcg[bgcnt].b.seq.s,'-')) {
				*bcg[bgcnt].b.dual = '\0';
				++bcg[bgcnt].b.dual;
				bcg[bgcnt].b.dual_n = strlen(bcg[bgcnt].b.dual);
			}
            // group pointer for this group
			bcg[bgcnt].gptr = getgroup(g);
            bcg[bgcnt].b.id.n=strlen(bcg[bgcnt].b.id.s);
            bcg[bgcnt].b.seq.n=strlen(bcg[bgcnt].b.seq.s);

            if (debug) fprintf(stderr, "BCG: %d bc:%s n:%d\n", bgcnt, bcg[bgcnt].b.seq.s, bcg[bgcnt].b.seq.n);
            ++bgcnt;
        }

		if (!bgcnt) {
			fprintf(stderr,"No barcodes %s from guide list %s.\n", group ? "matched" : "read", list);
			return 1;
		}

        int sampcnt = 200000;
        struct stat st;
		int fsum[f_n], fmax[f_n]; int bestcnt=0, besti=-1, bestdual=0;
		int dfsum[f_n], dfmax[f_n]; int dbestcnt=0, dbesti=-1;
		meminit(fsum); meminit(fmax); meminit(dfsum); meminit(dfmax);

        // subsample to determine group to use
		for (i=0;i<(usefile1?1:f_n);++i) {
            char *s = NULL; size_t na = 0; int nr = 0, ns = 0;
            char *q = NULL; size_t nq = 0;
			double tots=0, totsq=0;
		    char *s2 = NULL; int ns2=0;
            char *ignore_s=NULL;
	
			stat(in[i], &st);

			while ((ns=getline(&s, &na, fin[i])) > 0) {
				if (*s != '@')  {
					fprintf(stderr,"Invalid fastq file: %s.\n", in[i]);
					exit(1);
				}

                if (bcinheader) {
                    // read in 3 more lines (seq, comment, qual) and ignore them
                    ignore=getline(&ignore_s, &ignore_st, fin[i]);
                    ignore=getline(&ignore_s, &ignore_st, fin[i]);
                    ignore=getline(&ignore_s, &ignore_st, fin[i]);
                    getbcfromheader(s, &ns, &q, &s2, &ns2);
                    nq=ns;
               }  else {
                    // read in 3 more lines (seq, comment, qual)
                    if ((ns=getline(&s, &na, fin[i])) <=0)
                        break;

                    ignore=getline(&q, &ignore_st, fin[i]);
                    nq=getline(&q, &ignore_st, fin[i]);

    				s[--ns]='\0'; q[--nq]='\0';
                }

// skip if quality is below average
				if (st.st_size > (sampcnt * 500) && poorqual(i, ns, s, q)) 
					continue;
	
				for (b=0;b<bgcnt;++b) {
                    // matches front of read?
					if (!strncasecmp(s, bcg[b].b.seq.s, bcg[b].b.seq.n)) {
						++bcg[b].bcnt[i];
					} else if (!strncasecmp(s+1, bcg[b].b.seq.s, bcg[b].b.seq.n)) {
                        // shifted read?
						++bcg[b].bscnt[i];
                    }

					if (ns >= bcg[b].b.seq.n && !strcasecmp(s+ns-bcg[b].b.seq.n, bcg[b].b.seq.s)) {
						++bcg[b].ecnt[i]; 
					} else if (ns > bcg[b].b.seq.n && !strncasecmp(s+ns-bcg[b].b.seq.n-1, bcg[b].b.seq.s, bcg[b].b.seq.n)) {
						++bcg[b].escnt[i]; 
					}

					if (bcg[b].b.dual) {
                        const char * t=s;
                        int nt=ns;
                        if (bcinheader) {             // barcode in header? use stuff after '+' sign
                            t=s2;
                            nt=ns2;
                        }
						if (!strncasecmp(t, bcg[b].b.dual, bcg[b].b.dual_n)) {
							++bcg[b].dbcnt[i];
                        }
						if (ns >= bcg[b].b.dual_n && !strcasecmp(t+nt-bcg[b].b.dual_n, bcg[b].b.dual)) {
							++bcg[b].decnt[i];
                        }
					}
				}	
				
				++nr;
                // got enough reads?
				if (nr >= sampcnt) 
                    break;
			}

			for (b=0;b<bgcnt;++b) {
				// highest count
				int hcnt = (int) (max(bcg[b].bcnt[i],bcg[b].ecnt[i]) * log(bcg[b].b.seq.n));
				fsum[i]+=hcnt;
				if (hcnt > fmax[i])
					fmax[i]=hcnt;

				if (fsum[i] > bestcnt)  {
                    if (debug > 1) 
                        fprintf(stderr,"file %d(%s), bcg: %s, file-sum: %d, bestsum: %d\n", i, in[i], bcg[b].gptr->id, fsum[i], bestcnt);

					bestcnt=fsum[i];
					besti=i;
					bestdual=(bcg[b].b.dual!=NULL);
				}

                if (debug > 1) 
                    fprintf(stderr,"dual %d(%s), bcg: %s, file-sum: %d, bestsum: %d\n", i, in[i], bcg[b].gptr->id, dfsum[i], dbestcnt);

				if (bcg[b].b.dual) {
					// highest count
					int dcnt = (int) (max(bcg[b].dbcnt[i],bcg[b].decnt[i]) * log(bcg[b].b.dual_n));
					dfsum[i]+=dcnt;
					if (dcnt > dfmax[i])
						dfmax[i]=dcnt;
					if (dfsum[i] > dbestcnt)  {
                        if (debug > 1) 
                            fprintf(stderr,"dual %d(%s), bcg: %s, file-sum: %d, bestsum: %d\n", i, in[i], bcg[b].gptr->id, dfsum[i], dbestcnt);
						dbestcnt=dfsum[i];
						dbesti=i;
					}
				}
                        }
			if (debug > 0) fprintf(stderr,"file-best %d sum:%d, max:%d\n", besti, fsum[besti], fmax[besti]);
			if (debug > 0 && bestdual) fprintf(stderr,"dual file-best %d sum:%d, max:%d\n", dbesti, dfsum[dbesti], dfmax[dbesti]);
		}

        // chosen file is "besti"
		i=usefile1?0:besti;

		int gmax=0, gindex=-1, scnt = 0, ecnt=0, dscnt = 0, decnt = 0;
		int thresh = (int) (pickmaxpct*fmax[i]); 

		if (debug > 0) fprintf(stderr,"besti: %d thresh: %d, dual: %d\n", besti, thresh, bestdual);
		for (b=0;b<bgcnt;++b) {
			int hcnt = (int) (max(bcg[b].bcnt[i],bcg[b].ecnt[i]) * log(bcg[b].b.seq.n));
			if (debug > 1) fprintf(stderr,"cnt: %s %s hc:%d bc:%d ec: %d\n", bcg[b].b.id.s, bcg[b].b.seq.s, hcnt, bcg[b].bcnt[i], bcg[b].ecnt[i]);
			if (hcnt >= thresh) {
				// increase group count	
				bcg[b].gptr->tcnt += hcnt;
				if (bcg[b].gptr->tcnt > gmax) {
					gindex=bcg[b].gptr->i;
					gmax=bcg[b].gptr->tcnt;
				}
			}
		}
		if (gindex == -1) {
            fprintf(stderr, "Unable to determine barcode group\n");
			exit(1);
		}
//		printf("gmax: %d, gindex %d, %s, thresh: %d\n", gmax, gindex, grs[gindex].id, thresh);

        for (b=0;b<bgcnt;++b) {
			if (bcg[b].gptr->i == gindex) {
				if (bcg[b].bcnt[i] > bcg[b].ecnt[i]) {
					scnt+=bcg[b].dbcnt[i];
				} else if (bcg[b].bcnt[i] < bcg[b].ecnt[i]) {
					ecnt+=bcg[b].decnt[i];
				}
				if (bcg[b].dbcnt[dbesti] > bcg[b].decnt[dbesti]) {
					dscnt+=bcg[b].dbcnt[dbesti];
				} else if (bcg[b].dbcnt[dbesti] < bcg[b].decnt[dbesti]) {
					decnt+=bcg[b].decnt[dbesti];
				}
			}
		};
		end = scnt >= ecnt ? 'b' : 'e';

		if (debug) fprintf(stderr,"scnt: %d, ecnt, %d, end: %c\n", scnt, ecnt, end);

		thresh/=threshfactor;
		if (bestdual) 
		    thresh/=5;

		// since this is a known good set, use a very low threshold, just to catch them all
        fprintf(stderr, "Using Barcode Group: %s on File: %s (%s), Threshold %2.2f%%\n", 
        grs[gindex].id, in[i], endstr(end), 100.0 * (float) ((float)thresh/THFIXFACTOR)/sampcnt);

		if (bestdual) {
			dend = dscnt >= decnt ? 'b' : 'e';
			fprintf(stderr, "Dual index on File: %s (%s)\n", in[dbesti], endstr(dend));
			dual=true;
			for (b=0;b<bgcnt;++b) {
				// trim down a bit, but later should trim down to "both-match"
				if (bcg[b].gptr->i == gindex) {
					if (bcg[b].decnt[dbesti] < bcg[b].ecnt[i]) 
						bcg[b].ecnt[i] = bcg[b].decnt[dbesti];
					if (bcg[b].dbcnt[dbesti] < bcg[b].bcnt[i]) 
						bcg[b].bcnt[i] = bcg[b].dbcnt[dbesti];
				}
			}
		}

        for (b=0;b<bgcnt;++b) {
			if (bcg[b].gptr->i == gindex) {
				int cnt = (end == 'e' ? (bcg[b].ecnt[i]+bcg[b].escnt[i]) : ( bcg[b].bcnt[i] + bcg[b].bscnt[i] ));
				if (cnt > thresh/THFIXFACTOR) {
					// count exceeds threshold... use it
					bc[bcnt]=bcg[b].b;
					if ((end == 'e' && (bcg[b].escnt[i] < 1.2*bcg[b].ecnt[i])) ||
					    (end == 'b' && (bcg[b].bscnt[i] < 1.2*bcg[b].bcnt[i]))
					  ) {
						if (!dual)
							fprintf(stderr, "Using Barcode %s (%s)\n", bcg[b].b.id.s, bcg[b].b.seq.s);

						if (debug) fprintf(stderr, "Debug Barcode %s (%s-%s) ... ecnt:%d, escnt:%d,bcnt:%d, bscnt:%d\n", bcg[b].b.id.s, bcg[b].b.seq.s, bcg[b].b.dual, bcg[b].ecnt[i], bcg[b].escnt[i], bcg[b].bcnt[i], bcg[b].bscnt[i]);

					} else {
						bc[bcnt].shifted=1;

						if (!dual)
							fprintf(stderr, "Using Barcode %s (%s) shifted\n", bcg[b].b.id.s, bcg[b].b.seq.s);

						if (debug) printf("Debug Barcode %s (%s-%s) shifted ... ecnt:%d, escnt:%d,bcnt:%d, bscnt:%d\n", bcg[b].b.id.s, bcg[b].b.seq.s, bcg[b].b.dual, bcg[b].ecnt[i], bcg[b].escnt[i], bcg[b].bcnt[i], bcg[b].bscnt[i]);
					}
					++bcnt;
				}
			}
		}

		if (i != 0) {
			// in[0] needs to be the guide file
			FILE *f = fin[0];
			char *n = in[0];
			const char *o = out[0];
			bool gzi = gzin[0];
			fin[0]=fin[i];
			in[0]=in[i];
			out[0]=out[i];
			gzin[0]=gzin[i];
			fin[i]=f;
			in[i]=n;
			out[i]=o;
			gzin[i]=gzi;
			// swap file in to position 1 if dual
			if (dual && dbesti != 1) {
				FILE *f = fin[1];
				char *n = in[1];
				const char *o = out[1];
				bool gzi = gzin[1];
				fin[1]=fin[dbesti];
				in[1]=in[dbesti];
				out[1]=out[dbesti];
				gzin[1]=gzin[dbesti];
				fin[dbesti]=f;
				in[dbesti]=n;
				out[dbesti]=o;
				gzin[dbesti]=gzi;
			}
		}
		if (bcg) free(bcg);
	} else if (guide) {
		// use the first file as a "guide file" ... and select a set of codes from that
		FILE *gin = fin[0];

		int blen = 0;
	
		int sampcnt = 100000;
		struct stat st;
		stat(guide, &st);

		char *s = NULL; size_t na = 0; int nr = 0, ns = 0;
		char *q = NULL; size_t nq = 0;

// small sample to get lengths
		double tots=0, totsq=0;
		while ((ns=getline(&s, &na, gin)) > 0) {
			if (*s != '@')  {
				fprintf(stderr,"Invalid fastq file: %s.\n", in[0]);
				exit(1);
			}

            if (bcinheader) {
                    ignore=getline(&q, &ignore_st, fin[i]);
                    ignore=getline(&q, &ignore_st, fin[i]);
                    ignore=getline(&q, &ignore_st, fin[i]);
                    /// no dual barcode detection allowed
                    getbcfromheader(s, &ns);
                    printf("bc is %s\n", s);
            } else {
                if ((ns=getline(&s, &na, gin)) <=0)
                    break;
                ignore=getline(&q, &ignore_st, gin);
                ignore=getline(&q, &ignore_st, gin);
            }

			--ns;
			tots+=ns;
			totsq+=ns*ns;
			++nr;
			if (nr >= 200) break;
		}
		double dev = stdev(nr, tots, totsq);

		// short, and nonvarying (by much, depends on the tech used)
		if (dev < .25 && roundl(tots/nr) < 12) {
			// most probably a barcode-only read
			blen = (int) round(tots/nr);
			end = 'b';
		} else if (round(tots/nr) < 12) {
			fprintf(stderr, "File %s looks to be barcode-only, but it's length deviation is too high (%.4g)\n", in[0], dev);
			return 1;
		} else {
			fprintf(stderr, "File %s isn't a barcode-only file, try using -l instead\n", in[0]);
			return 1;
		}

		fprintf(stderr, "Barcode length used: %d (%s)\n", blen, endstr(end));

		// load a table of possble codes
		pickmax=0;
		picktab=NULL;
		bnode * ent = NULL;
        while ((ns=getline(&s, &na, gin)) > 0) {
			if (*s != '@')  {
				fprintf(stderr,"Invalid fastq file: %s.\n", in[i]);
				exit(1);
			}

            if (bcinheader) {
                ignore=getline(&q, &ignore_st, fin[i]);
                ignore=getline(&q, &ignore_st, fin[i]);
                ignore=getline(&q, &ignore_st, fin[i]);
                getbcfromheader(s, &ns, &q);
                printf("bc is %s\n", s);
            } else {
                if ((ns=getline(&s, &na, gin)) <=0)
                    break;

                ignore=getline(&q, &ignore_st, gin);
                if (getline(&q, &ignore_st, gin) != ns)
                    break;
			    s[--ns]='\0'; q[ns]='\0';
            }

			if (st.st_size > (sampcnt * 500) && poorqual(i, ns, s, q)) 
				continue;

            ++nr;

			char *p;
			if (end == 'b') {
				p=s;
			} else {
				p=s+nr-blen;
			}
			p[blen]='\0';
			if (!ent) {		// make a new ent 
				ent = (bnode *) malloc(sizeof(*ent));
			    ent->seq=(char*)malloc(blen+1);;
            }

			if (strchr(p, 'N')||strchr(p, 'n'))
				continue;

			ent->cnt=0;
			strcpy(ent->seq, p);

			bnode *fent = * (bnode**)  tsearch(ent, &picktab, bnodecomp);

			if (fent == ent)	// used the ent, added to tree
				ent = NULL;	// need a new one

			++fent->cnt;

			if (fent->cnt > pickmax) 
                pickmax=fent->cnt;
			else if (fent->cnt > pickmax2) 
                pickmax2=fent->cnt;
			
			if (nr > sampcnt)
				break;
		}
		pickmax=max(1,(int)(pickmaxpct*pickmax2));
		fprintf(stderr, "Threshold used: %d\n", pickmax);
		twalk(picktab, pickbest);
	} else {
		// user specifies a list of barcodes, indexed read is f[0] and f[1] if dual
		FILE *bin = fopen(bfil, "r");
		if (!bin) {
			fprintf(stderr, "Error opening file '%s': %s\n",bfil, strerror(errno));
			return 1;
		}

		bcnt = 0;
		int ok;
		while (bcnt < MAX_BARCODE_NUM && (ok = read_line(bin, bc[bcnt].id))) {
			if (ok <= 0) break;
			if (bc[bcnt].id.s[0]=='#') continue;
			bc[bcnt].id.s=strtok(bc[bcnt].id.s, "\t\n\r ");
			bc[bcnt].seq.s=strtok(NULL, "\t\n\r ");
			if (!bc[bcnt].seq.s) {
				fprintf(stderr, "Barcode file '%s' required format is 'ID SEQ'\n",bfil);
				return 1;
			}
            if (bc[bcnt].dual=strchr(bc[bcnt].seq.s,'-')) {
                *bc[bcnt].dual = '\0';
                ++bc[bcnt].dual;
				bc[bcnt].dual_n = strlen(bc[bcnt].dual);
				dual=true;
            }
			bc[bcnt].id.n=strlen(bc[bcnt].id.s);
			bc[bcnt].seq.n=strlen(bc[bcnt].seq.s);
			if (debug) fprintf(stderr, "BC: %d bc:%s n:%d\n", bcnt, bc[bcnt].seq.s, bc[bcnt].seq.n);
			++bcnt; 
		}

        fprintf(stderr, "Using Barcode File: %s\n", bfil);
	}

	if (noexec) {
		int b;
        	for (b=0;b<bcnt;++b) {
			fprintf(stdout, "%s %s\n", bc[b].id.s, bc[b].seq.s);
		}
		exit(0);
	}

	// for whatever reason, the end is not supplied... easy enough to determine accurately
	// or it's dual... which means we need to resample stuff
    if (bcinheader && !end) {
        end = 'b';
    }

    if (end == '\0' || dual) {
        for (i=0;i<f_n;++i) {
            if (!gzin[i])
                fseek(fin[i],0,0);
            else {
                pclose(fin[i]);
                fin[i]=gzopen(in[i],"r",&gzin[i]);
            }
        }

        int sampcnt = dual ? 200000 : 10000;
        struct stat st;
        stat(in[0], &st);
        char *s = NULL; size_t na = 0; int nr = 0, ns = 0;
        char *q = NULL; size_t nq = 0;
        int ne=0, nb=0, dne=0, dnb=0, tcount=0, read_ok=0;

        int *recount = dual ? ((int *) malloc(sizeof(int)*bcnt)) : NULL;
        if (dual) memset(recount, 0, sizeof(int)*bcnt);

        struct fq fq[2]; meminit(fq);

        while (read_ok=read_fq(fin[0], nr, &fq[0])) {
            fq[0].id.s[--fq[0].id.n]='\0';

            if (dual)
                read_fq(fin[1], nr, &fq[1]);
            ++nr;

            if (st.st_size > (sampcnt * 500) && poorqual(0, fq[0].seq.n, fq[0].seq.s, fq[0].qual.s)) 
                continue;

            if (dual)
                if (st.st_size > (sampcnt * 500) && poorqual(1, fq[1].seq.n, fq[1].seq.s, fq[1].qual.s))
                    continue;

            for (i=0;i<bcnt;++i) {
                int dok = 0;
                if (debug > 1) fprintf(stderr, "check %s vs %s: %s vs %s", fq[0].id.s, bc[i].id.s, fq[0].seq.s, bc[i].seq.s);
                if (!strncmp(fq[0].seq.s, bc[i].seq.s, bc[i].seq.n)) {
                    ++nb;
                    ++dok;
                } else if (!strncmp(fq[0].seq.s+fq[0].seq.n-bc[i].seq.n, bc[i].seq.s, bc[i].seq.n)) {
                    ++ne;
                    ++dok;
                }
                if (dual) {
                    if (debug > 1) fprintf(stderr, ", dual: %s vs %s, ", fq[1].seq.s, bc[i].dual);
                    if (!strncmp(fq[1].seq.s, bc[i].dual, bc[i].dual_n)) {
                        ++dnb;
                        ++dok;
                    } else if (!strncmp(fq[1].seq.s+fq[1].seq.n-bc[i].dual_n, bc[i].dual, bc[i].dual_n)) {
                        ++dne;
                        ++dok;
                    }
                }
                if (debug > 1) fprintf(stderr, ", dok:%d, i:%d\n", dok, i);
                if (dok == 2) {
                    ++recount[i];
                    ++tcount;
                    break;
                }
            }

            if (nr >= sampcnt) 
                break;
        }

        end = (ne > nb) ? 'e' : 'b';
        fprintf(stderr, "End used: %s\n", endstr(end));

        if (dual && list) {
            // trim down possiblities to reduce number of open files, and small stub files
            dend = (dne > dnb) ? 'e' : 'b';
            fprintf(stderr, "Dual-end used: %s\n", endstr(dend));
            int ocnt = bcnt;
            // this should allow up to a 300 plex
            int thresh = max(1,tcount/1000);
            thresh /= threshfactor;
            bcnt=0;
            if (debug)
                fprintf(stderr, "dual resample threshold: %d out of %d\n", thresh, tcount);
            for (i=0;i<ocnt;++i) {
                if (recount[i] >= thresh) {
                    fprintf(stderr, "Using Barcode %s (%s-%s)\n", bc[i].id.s, bc[i].seq.s, bc[i].dual);
                    if (debug)
                        fprintf(stderr, "%d >= %d\n", recount[i], thresh);
                    bc[bcnt].seq=bc[i].seq;
                    bc[bcnt].id=bc[i].id;
                    bc[bcnt].dual=bc[i].dual;
                    bc[bcnt].dual_n=bc[i].dual_n;
                    ++bcnt;
                } else {
                    if (debug)
                        fprintf(stderr, "skipping barcode %s (%s-%s), %d < %d\n", bc[i].id.s, bc[i].seq.s, bc[i].dual, recount[i], thresh);
                }
            }
        }
    }

	if (bcnt == 0) { 
		fprintf(stderr, "No barcodes defined, quitting.\n");
		exit(1);
	}

	// one beyond barcode count is unmatched
	bc[bcnt].id.s=(char *)"unmatched";

	// TODO: output barcode read ...but only for unmatched?
	int b;
    for (b=0;b<=bcnt;++b) {
		size_t nameseq_len = strlen(bc[b].id.s);
		if ((b < bcnt) && seqnames) {
                  nameseq_len = strlen(bc[b].seq.s);
                  if (bc[b].dual)
                    nameseq_len += bc[b].dual_n + 1;
		}

		for (i=0;i<f_n;++i) {
			if (!strcasecmp(out[i],"n/a") || !strcasecmp(out[i],"/dev/null")) {
				bc[b].out[i] = NULL;
				bc[b].fout[i] = NULL;
				continue;
			}
			const char *p=strchr(out[i],'%');
			if (!p) fail("Each output file name must contain a '%%' sign, which is replaced by the barcode id or sequence\n");
			bc[b].out[i]=(char *) malloc(strlen(out[i])+nameseq_len+100);
			strncpy(bc[b].out[i], out[i], p-out[i]);
			bc[b].out[i][p-out[i]]='\0';
			if (seqnames && (b < bcnt)) {
			  strcat(bc[b].out[i], bc[b].seq.s);
                          if (bc[b].dual) {
                            strcat(bc[b].out[i], "-");
                            strcat(bc[b].out[i], bc[b].dual);
                          }
                        }
			else {
			  strcat(bc[b].out[i], bc[b].id.s);
                        }
			strcat(bc[b].out[i], p+1);
			if (!(bc[b].fout[i]=gzopen(bc[b].out[i], "w", &bc[b].gzout[i]))) {
				fprintf(stderr, "Error opening output file '%s': %s\n",bc[b].out[i], strerror(errno));
				return 1;
			}
		}
	}

	// seek back to beginning of fastq
	for (i=0;i<f_n;++i) {
		if (!gzin[i])
			fseek(fin[i],0,0);
		else {
			pclose(fin[i]);
			fin[i]=gzopen(in[i],"r",&gzin[i]);
		}
	}

    // don't trim if you're not outputting the read

	struct fq fq[8];

    meminit(fq);

	int nrec=0;
	int nerr=0;
	int nok=0;
	int ntooshort=0;
	int ntrim=0;
	int nbtrim=0;
	int read_ok;

    // ACTUAL DEMUX HAPPENS HERE
	// read in 1 record from EACH file supplied
	while (read_ok=read_fq(fin[0], nrec, &fq[0])) {
		for (i=1;i<f_n;++i) {
			int mate_ok=read_fq(fin[i], nrec, &fq[i]);
			if (read_ok != mate_ok) {
				fprintf(stderr, "# of rows in mate file '%s' doesn't match primary file, quitting!\n", in[i]);
				return 1;
			}
			if (verify) {
				// verify 1 in 100
				if (0 == (nrec % 100)) {
					char *p=strchr(fq[i].id.s,verify);
					if (!p) {
						fprintf(stderr, "File %s is missing id verification char %c at line %d", in[i], verify, nrec*4+1);
						return 1;
					}
					int l = p-fq[i].id.s;
					if (strncmp(fq[0].id.s, fq[i].id.s, l)) {
						fprintf(stderr, "File %s, id doesn't match file %s at line %d", in[0], in[i], nrec*4+1);
						return 1;
					}
				}
			}
		}
		++nrec;
		if (read_ok < 0) continue;

		int i, best=-1, bestmm=mismatch+distance+1, bestd=mismatch+distance+1, next_best=mismatch+distance*2+1;

        if (bcinheader) {
            for (i=f_n-1;i>=0;--i) {
                fq[i+(dual?2:1)]=fq[i];
            }
            meminit(fq[0]); 
            if (dual) {
                meminit(fq[1]); 
                getbcfromheader(&fq[2], &fq[0], &fq[1].seq.s, &fq[1].seq.n);
            } else {
                getbcfromheader(&fq[2], &fq[0]);
            }
        }

		if (debug) {
			if (!bcinheader) fq[0].id.s[fq[0].id.n-1] = '\0';
			fprintf(stderr, "id: %s, seq: %s %d", fq[0].id.s, fq[0].seq.s, fq[0].seq.n);
			if (dual) fprintf(stderr, ", sdual: %s %d", fq[1].seq.s, fq[1].seq.n);
			if (!bcinheader) fq[0].id.s[fq[0].id.n] = '\n';
			if (debug > 1) printf("\n");
		}

        if (quality > 0) {
            // low quality base = 'N'
            for (i=0;i<fq[0].seq.n;++i) {
                if (fq[0].qual.s[i]<quality) {
                    fq[0].seq.s[i]='N';
                }
            }
        }

        // for each barcode
        for (i =0; i < bcnt; ++i) {
            int d;
            if (end == 'e') {
                if (bc[i].shifted) {
                    if (fq[0].seq.n > bc[i].seq.n) {
                        d=hd(fq[0].seq.s+fq[0].seq.n-bc[i].seq.n-1, bc[i].seq.s, bc[i].seq.n);
                    } else {
                        d=bc[i].seq.n;
                    }
                } else {
                    if (fq[0].seq.n >= bc[i].seq.n) {
                        d=hd(fq[0].seq.s+fq[0].seq.n-bc[i].seq.n, bc[i].seq.s, bc[i].seq.n);
                    } else {
                        d=bc[i].seq.n;
                    }
                }

                if (dual) {
                    // distance is added in for duals
                    if (fq[1].seq.n >= bc[i].dual_n) {
                        d+=hd(fq[1].seq.s+fq[1].seq.n-bc[i].dual_n, bc[i].dual, bc[i].dual_n);
                    } else {
                        d+=bc[i].dual_n;
                    }
                }
            } else {
                if (bc[i].shifted) 
                    d=hd(fq[0].seq.s+1,bc[i].seq.s, bc[i].seq.n);
                else
                    d=hd(fq[0].seq.s,bc[i].seq.s, bc[i].seq.n);

                // distance is added in for duals
                if (dual) 
                    d+=hd(fq[1].seq.s,bc[i].dual, bc[i].dual_n);

                //				if (debug > 1) {
                //					fprintf(stderr, "index: %d dist: %d bc:%s n:%d", i, d, bc[i].seq.s, bc[i].seq.n);
                //					if (dual) fprintf(stderr, ", idual: %s %d", bc[i].dual, bc[i].dual_n);
                //					fprintf(stderr, "\n");
                //				}
            }
            // simple... 
            if (d < bestd) {
                next_best=bestd;
                bestd=d;
                if (debug > 1) fprintf(stderr,"next_dist: %d, best_seq: %s:%d\n", next_best, bc[i].seq.s, bestd);
            }
            // if exact match
            if (d==0) { 
                if (debug) fprintf(stderr, ", found bc: %d bc:%s n:%d, bestd: %d, next_best: %d", i, bc[i].seq.s, bc[i].seq.n, bestd, next_best);
                best=i; 
                break; 
            } else if (d <= mismatch) {
                // if ok match
                if (d == bestmm) {
                    best=-1;		// more than 1 match... bad
                } else if (d < bestmm) {
                    bestmm=d;		// best match...ok
                    best=i;
                }
            }
        }

        if ((best >= 0) && distance && (next_best-bestd) < distance) {
            if (debug) fprintf(stderr, "%d<%d, skipping", next_best-bestd, distance);
            // match is ok, but distance is poor
            ++poor_distance;
            best=-1;
        }

        bool trimmed = false;
        // only trim if you're outputting the sequence
		if (trim && best >= 0 && bc[best].fout[0]) {
			// todo: save trimmed
            trimmed = true;
			int len=bc[best].seq.n;
			if (end =='b') {
				memmove(fq[0].seq.s, fq[0].seq.s+len, fq[0].seq.n-len);
				memmove(fq[0].qual.s, fq[0].qual.s+len, fq[0].seq.n-len);
			}
			fq[0].seq.s[fq[0].seq.n-len]='\0';
			fq[0].qual.s[fq[0].qual.n-len]='\0';
		}

		if (best < 0) {
            // shuttle to unmatched file
			best=bcnt;
		}

		if (debug) fprintf(stderr, ", best: %d %s\n", best, bc[best].id.s);

		++bc[best].cnt;

        int shift_index=0;
        if (bcinheader) {
            shift_index = 1;
            if (dual) 
                shift_index = 2;
        }

		for (i=shift_index;i<f_n+shift_index;++i) {
			FILE *f=bc[best].fout[i-shift_index];
			if (!f) continue;
            if (!trimmed) {
			    // todo: capture always, not just when trim is off
                if (!debug) *strrchr(fq[i].id.s, '\n') = '\0';
                fputs(fq[i].id.s,f);
                fputc(' ', f);
                fputs(fq[0].seq.s,f);
                if (dual) {
                    fputc('-', f);
                    fputs(fq[1].seq.s,f);
                }
                fputc('\n', f);
            } else {
                // id still has chr
                fputs(fq[i].id.s,f);
            }
            fputs(fq[i].seq.s,f);
            fputc('\n',f);
            fputs(fq[i].com.s,f);
            fputs(fq[i].qual.s,f);
            fputc('\n',f);
		}
	}

    bool io_ok=1;
    for (b=0;b<=bcnt;++b) {
        for (i=0;i<f_n;++i) {
            if (bc[b].fout[i]) {
                if (bc[b].gzout[i]) {
                    io_ok = io_ok && !pclose(bc[b].fout[i]);
                } else {
                    io_ok = io_ok && !fclose(bc[b].fout[i]);
                }
            }
        }
    }


    if (poor_distance > 0)
        fprintf(stderr, "Skipped because of distance < %d : %d\n", distance, poor_distance);

    if (!io_ok)
        fprintf(stderr, "Returning error because of i/o error during file close\n");

	int j;
	printf("Id\tCount\tFile(s)\n");
	uint64_t tot=0;
	for (i=0;i<=bcnt;++i) {
		printf("%s\t%lu", bc[i].id.s, bc[i].cnt);
		tot+=bc[i].cnt;
		for (j=0;j<f_n;++j) {
			if (bc[i].out[j])
				printf("\t%s", bc[i].out[j]);
		}
		printf("\n");
	}
	printf("total\t%lu\n", tot);

    if (!io_ok)
        return 3;

	return 0;
}

struct group* getgroup(char *s) {
	int i;
	for (i=0;i<grcnt;++i) {
		if (!strcasecmp(s,grs[i].id)) {
			return &grs[i];
		}
	}
    if (grcnt >= MAX_GROUP_NUM) {
        fprintf(stderr,"Too many barcode groups, quitting\n");
        exit(1);
    }
	grs[grcnt].id=(char *)malloc(strlen(s)+1);
    strcpy(grs[grcnt].id,s);
	grs[grcnt].tcnt=0;
	grs[grcnt].i=grcnt;
	return &grs[grcnt++];
}

void pickbest(const void *nodep, const VISIT which, const int depth)
{
	if (which==endorder || which==leaf) {
		bnode *ent = *(bnode **) nodep;
//		 printf("HERE!! %s, %d, %d\n", ent->seq, ent->cnt, pickmax);
		// allow one sample to be as much as 1/10 another, possibly too conservative
		if (ent->cnt > pickmax && bcnt < MAX_BARCODE_NUM) {
			bc[bcnt].seq.s=ent->seq;
			bc[bcnt].id.s=ent->seq;
			bc[bcnt].id.n=strlen(bc[bcnt].id.s);
			bc[bcnt].seq.n=strlen(bc[bcnt].seq.s);
			++bcnt;
		} else {
			//free(ent->seq);
		}
		//free(ent);
		//tdelete((void*)ent, &picktab, scompare);
	}
}

void usage(FILE *f) {
	fprintf(f,
"Usage: fastq-multx [-g|-l|-B] <barcodes.fil> <read1.fq> -o r1.%%.fq [mate.fq -o r2.%%.fq] ...\n"
"Version: %s\n"
"\n"
"Output files must contain a '%%' sign which is replaced with the barcode id in the barcodes file.\n"
"Output file can be n/a to discard the corresponding data (use this for the barcode read)\n"
"\n"
"Barcodes file (-B) looks like this:\n"
"\n"
"<id1> <sequence1>\n"
"<id2> <sequence2> ...\n"
"\n"
"Default is to guess the -bol or -eol based on clear stats.\n"
"\n"
"If -g is used, then it's parameter is an index lane, and frequently occuring sequences are used.\n"
"\n"
"If -l is used then all barcodes in the file are tried, and the *group* with the *most* matches is chosen.\n" 
"\n"
"Grouped barcodes file (-l or -L) looks like this:\n"
"\n"
"<id1> <sequence1> <group1>\n"
"<id1> <sequence1> <group1>\n"
"<id2> <sequence2> <group2>...\n"
"\n"
"Mated reads, if supplied, are kept in-sync\n"
"\n"
"Options:\n"
"\n"
"-o FIL1     Output files (one per input, required)\n"
"-g SEQFIL   Determine barcodes from the indexed read SEQFIL\n"
"-l BCFIL    Determine barcodes from any read, using BCFIL as a master list\n"
"-L BCFIL    Determine barcodes from <read1.fq>, using BCFIL as a master list\n"
"-B BCFIL    Use barcodes from BCFIL, no determination step, codes in <read1.fq>\n"
"-H          Use barcodes from illumina's header, instead of a read\n"
"-s          Substitute barcode sequence instead of barcode label into output file names\n"
"-b          Force beginning of line (5') for barcode matching\n"
"-e          Force end of line (3') for batcode matching\n"
"-t NUM      Divide threshold for auto-determine by factor NUM (1), > 1 = more sensitive\n"
"-G NAME     Use group(s) matching NAME only\n"
"-x          Don't trim barcodes off before writing out destination\n"
"-n          Don't execute, just print likely barcode list\n"
"-v C        Verify that mated id's match up to character C (Use ' ' for illumina)\n"
"-m N        Allow up to N mismatches, as long as they are unique (1)\n"
"-d N        Require a minimum distance of N between the best and next best (2)\n"
"-q N        Require a minimum phred quality of N to accept a barcode base (0)\n"
	,VERSION);
}

void getbcfromheader(struct fq *fq, struct fq *bc, char **s2, int *ns2) {
    // reallocate bc to match fq
    // warning... result has no newline!

    // copy id to sequence
    bc->seq.s=(char *)realloc(bc->seq.s,fq->id.n+1);
    strncpy(bc->seq.s,fq->id.s,fq->id.n+1);
    bc->seq.n=fq->id.n;
    bc->seq.a=fq->id.n+1;


    // this extracts the new sequence
    getbcfromheader(bc->seq.s, &(bc->seq.n), &(bc->qual.s), s2, ns2);
    bc->qual.n=bc->seq.n;

//printf("DEBUG: seq is %s, length is %d\n", bc->seq.s, bc->seq.n);
}

// looks for barcode in s, totally replaces s with barcode only, sets ns to length
void getbcfromheader(char *s, int *ns, char **q, char **s2, int *ns2) {
    char *p=strchr(s, ' ');
    if (!p) {
        fprintf(stderr,"Barcode not in header: %s.\n", s);
        exit(1);
    }

    char *t;
    while(t=strchr(p,':')) {
        p=t+1;
    }

    // remove newline
    if (s[*ns-1] == '\n') {
        if (s[*ns-1] == '\r') {
            --*ns;
        }
        --*ns;
        s[*ns]='\0';
    }

    // result has no newline
    *ns-=(p-s);
    memmove(s,p,*ns);
    s[*ns]='\0';

    if (q) {
        *q=(char*)realloc(*q,(*ns)+1);
        memset(*q,'h',*ns);
        (*q)[*ns]='\0';
    }

    if (p=strchr(s,'+')) {
        *p='\0';
        *ns = p-s;

        if (ns2) {
            *ns2=(*ns-((int)(p-s))-1);
            *s2=p+1;
        } else {
            // ERROR: maybe die here?   Or assume the user knows what's up?
        }
    }
}

