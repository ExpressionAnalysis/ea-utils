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
#define SCANLEN 15
#define SCANMIDP ((int) SCANLEN/2)

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
	int bcnt;			// number found at beginning
	int bcntz;			// number found at beginning
	int ecnt;			// number found at end
	int ecntz;			// number found at end
	char end;			// 'b' or 'e'
	int thr;			// min-length for clip
};

int read_fa(FILE *in, int rno, struct ad *ad);		// 0=done, 1=ok, -1=err+continue
int read_fq(FILE *in, int rno, struct fq *fq);		// 0=done, 1=ok, -1=err+continue

void usage(FILE *f);
int hd(char *a, char *b, int n);
int debug=0;
int main (int argc, char **argv) {
	char c;
	bool eol;
	int nmin = 1, nkeep = 15;
	float minpct = 0.25;
	int pctdiff = 20;
	char *outfile = NULL;
	int sampcnt = 40000;
	int xmax = -1;
	float scale = 2.2;
	int noclip=0;
	char end = '\0';

	char *mate_in[5];
	char *mate_out[5];
	int mate_n=0;
	int mate_oarg=0;

	int i;
	
	char *afil = NULL, *ifil = NULL;
	while (	(c = getopt (argc, argv, "-ndbehp:o:l:s:m:t:")) != -1) {
		switch (c) {
		case '\1': 
			if (!afil) 
				afil = optarg; 
			else if (!ifil) 
				ifil = optarg; 
			else if (mate_n<3) 
				mate_in[mate_n++] = optarg; 
			else {
				usage(stderr); return 1;
			}
			break;
		case 't': minpct = atof(optarg); break;
		case 'm': nmin = atoi(optarg); break;
		case 'l': nkeep = atoi(optarg); break;
		case 'p': pctdiff = atoi(optarg); break;
		case 'h': usage(stdout); return 1; 
		case 'o': if (!outfile) 
				outfile = optarg; 
			  else if (mate_oarg < 3) 
				mate_out[mate_oarg++] = optarg;
			  break;
		case 's': scale = atof(optarg); break;
		case 'i': ifil = optarg; break;
		case 'n': noclip = 1; break;
		case 'd': debug = 1; break;
		case 'b': end = 'b'; break;
		case 'e': end = 'e'; break;
		case '?': 
		     if (strchr("lm", optopt))
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

	if (!(noclip && !outfile) && mate_n != mate_oarg) {
		fprintf(stderr, "Error: number of input files must match number of '-o' output files.\n");
		return 1;
	}

	if (argc < 3 || !afil || !ifil) {
		usage(stderr);
		return 1;
	}

	FILE *fin = fopen(ifil, "r"); 
	if (!fin) {
		fprintf(stderr, "Error opening file '%s': %s\n",ifil, strerror(errno));
		return 1;
	}

        FILE *ain = fopen(afil, "r");
        if (!ain) {
                fprintf(stderr, "Error opening file '%s': %s\n",afil, strerror(errno));
                return 1;
        }

	FILE *fout = stdout;
	FILE *fstat = stderr;
	if (outfile  && !noclip) {
		fstat = stdout;
	}
	if (noclip) {
		fstat = stdout;
	}

	FILE *mate_fin[5];
	FILE *mate_fout[5];

	for (i=0;i<mate_n;++i) {
		if (!(mate_fin[i]=fopen(mate_in[i], "r"))) {
			fprintf(stderr, "Error opening file '%s': %s\n",mate_in[i], strerror(errno));
			return 1;
		}
		if (!noclip)
                if (!(mate_fout[i]=fopen(mate_out[i], "w"))) {
                        fprintf(stderr, "Error opening output file '%s': %s\n",mate_out[i], strerror(errno));
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

	struct stat st;
	stat(ifil, &st);
	fseek(fin, st.st_size > sampcnt ? (st.st_size-sampcnt)/3 : 0, 0);

	char *s = NULL; size_t na = 0; int nr = 0, ns = 0;
        while (getline(&s, &na, fin) > 0) {
		if (*s == '@')  {
			if ((ns=getline(&s, &na, fin)) <=0) 
				break;
			--ns;				// don't count newline for read len
			++nr;
			int a;
			char buf[SCANLEN+1];
			strncpy(buf, s, SCANLEN);
			for(a=0;a<acnt;++a) {
				char *p;
				if (p = strstr(s, ad[a].escan)) { 
			//		if (debug) fprintf(stderr, "END  : A: %s (%s), P: %d, SL: %d, Z:%d\n", ad[a].id, ad[a].escan, p-s, ns, (p-s) == ns-SCANLEN);
					if ((p-s) == ns-SCANLEN) 
						++ad[a].ecntz;
					++ad[a].ecnt;
				}
					
				if ((p = strstr(ad[a].seq, buf))) { 
			//		if (debug) fprintf(stderr, "BEGIN: A: %s (%s), P: %d, SL: %d, Z:%d\n", ad[a].id, ad[a].seq, p-ad[a].seq, ns, (p-ad[a].seq )  == ad[a].nseq-SCANLEN);
					if (p-ad[a].seq == ad[a].nseq-SCANLEN) 
						++ad[a].bcntz;
					++ad[a].bcnt;
				}
			}
		}
		if (nr >= sampcnt) 
			break;
        }

	int a;
	int athr = (int) ((float)nr * minpct) / 100;
	fprintf(fstat, "Scale used: %g\n", scale);
	fprintf(fstat, "Threshold used: %d out of %d\n", athr+1, nr);
	int newc=0;
	for(a=0;a<acnt;++a) {
		if (ad[a].ecnt > athr || ad[a].bcnt > athr) {
			int cnt;
			if (debug) fprintf(stderr, "EC: %d, BC:%d, ECZ: %d, BCZ: %d\n", ad[a].ecnt, ad[a].bcnt, ad[a].ecntz, ad[a].bcntz);
			// heavily weighted toward start/end maches
			if ((ad[a].ecnt + 10*ad[a].ecntz) >= (ad[a].bcnt + 10*ad[a].bcntz)) {
				ad[a].end='e';
				cnt = ad[a].ecnt;
			} else {
				ad[a].end='b';
				cnt = ad[a].bcnt;
			}
			
			// user supplied end.... don't clip elsewhere
			if (end && ad[a].end != end)
				continue;

			ad[a].thr = max(nmin,(int) (-log(cnt / (float) nr)/log(scale)));
			fprintf(fstat, "Adapter %s (%s): counted %d at the '%s' of sequences, clip set to %d", ad[a].id, ad[a].seq, cnt, ad[a].end == 'e' ? "end" : "start", ad[a].thr);
			if (abs((ad[a].bcnt-ad[a].ecnt)) < athr/4) {
				fprintf(fstat, ", warning end was not reliable\n", ad[a].id, ad[a].seq);
			} else {
				fputc('\n', fstat);
			}

			ad[newc++]=ad[a];
		}
	}
	acnt=newc;

	if (acnt == 0) {
		fprintf(fstat, "No adapters found, no clipping needed\n");
		if (noclip) exit (1);			// for including in a test
		exit(0);				// not really an error, check size of output files
	}

	if (noclip)
		exit(0);

        if (outfile) {
                fout = fopen(outfile, "w");
                if (!fout) {
                        fprintf(stderr, "Error opening output file '%s': %s",outfile, strerror(errno));
                        return 1;
                }
        }

        for (i=0;i<mate_n;++i) {
                if (!(mate_fout[i]=fopen(mate_out[i], "w"))) {
                        fprintf(stderr, "Error opening output file '%s': %s\n",mate_out[i], strerror(errno));
                        return 1;
                }
        }

	struct fq fq;	
        memset(&fq, 0, sizeof(fq));
	struct fq mate_fq[5];	
        memset(&mate_fq, 0, 5*sizeof(fq));

	int nrec=0;
	int nerr=0;
	int nok=0;
	int ntooshort=0;
	int ntrim=0;
	int nbtrim=0;
	int read_ok;
	fseek(fin, 0, 0);
	while (read_ok=read_fq(fin, nrec, &fq)) {
		for (i=0;i<mate_n;++i) {
			read_ok=read_fq(mate_fin[i], nrec, &mate_fq[i]);
			if (!read_ok) {
				fprintf(stderr, "# of rows in mate file '%s' doesn't match primary file, quitting!\n", mate_in[i]);
				return 1;
			}
		}
		++nrec;
		if (read_ok < 0) continue;

		// chomp

		if (debug) fprintf(stderr, "seq: %s %d\n", fq.seq, fq.nseq);

		bool skip = 0;
		int bestscore_e = INT_MAX, bestoff_e = 0, bestlen_e = 0; 
		int bestscore_b = INT_MAX, bestoff_b = 0, bestlen_b = 0; 

		for (i =0; i < acnt; ++i) {
			int nmatch = ad[i].thr;
			if (!nmatch) nmatch = ad[i].nseq;			// full match required if nmin == 0
	
			// how far in to search for a match?
			int mx = ad[i].nseq;
			if (xmax) {
				 mx = fq.nseq;
				 if (xmax > 0 && (xmax+ad[i].nseq) < mx)
					mx = xmax+ad[i].nseq;			// xmax is added to adapter length
			}

			if (debug)
				fprintf(stderr, "adapter: %s, adlen: %d, nmatch: %d, mx: %d\n", ad[i].seq, ad[i].nseq, nmatch, mx);

			if (ad[i].end == 'e') {
			int off;
			for (off = nmatch; off <= mx; ++off) {			// off is distance from tail of sequence
				char *seqtail = fq.seq+fq.nseq-off; 		// search at tail
				int ncmp = off<ad[i].nseq ? off : ad[i].nseq;
				int mind = (pctdiff * ncmp) / 100;
				int d = hd(ad[i].seq,seqtail,ncmp);		// # differences
				if (debug>1)
					fprintf(stderr, "tail: %s, bestoff: %d, off: %d, ncmp: %d, mind: %d, hd %d\n", seqtail, bestoff_e, off, ncmp, mind, d);
				if (d <= mind) {
					int score = (1000*(d*d+1))/ncmp;
					if (score <= bestscore_e) {			// better score?
						bestscore_e = score;			// save max score
						bestoff_e = off;				// offset at max
						bestlen_e = ncmp;				// cmp length at max
					}
					if (d == 0 && (ncmp == ad[i].nseq)) {
						break;
					}
				}
			}
			} else {
                        int off;
                        for (off = nmatch; off <= mx; ++off) {                  // off is distance from start of sequence
                                int ncmp = off<ad[i].nseq ? off : ad[i].nseq;	// number we are comparing
                                char *matchtail = ad[i].seq+ad[i].nseq-ncmp;    // tail of adapter
                                char *seqstart = fq.seq+off-ncmp;		// offset into sequence (if any)
                                int mind = (pctdiff * ncmp) / 100;
                                int d = hd(matchtail,seqstart,ncmp);             // # differences
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

			int totclip = bestoff_e + bestoff_b;
			if (totclip > 0) {
				if ( (fq.nseq-totclip) < nkeep) {
					++ntooshort;
					skip = 1;
					break;
				}
			}
		}	

		if (!skip) {
			if (bestoff_b > 0 || bestoff_e > 0) {
				++ntrim;
				if (bestoff_e > 0) {
					if (debug) printf("trimming %d from end\n", bestoff_e);
					fq.seq [fq.nseq -bestoff_e]='\0';
					fq.qual[fq.nqual-bestoff_e]='\0';
				}
				if (bestoff_b > 0) {
					if (debug) printf("trimming %d from begin\n", bestoff_b);
					memmove(fq.seq ,fq.seq +bestoff_b,fq.nseq -=bestoff_b);
					memmove(fq.qual,fq.qual+bestoff_b,fq.nqual-=bestoff_b);
					fq.seq[fq.nseq]='\0';
					fq.qual[fq.nqual]='\0';
				}
			}
			fputs(fq.id,fout);
			fputs(fq.seq,fout);
			fputc('\n',fout);
			fputs(fq.com,fout);
			fputs(fq.qual,fout);
			fputc('\n',fout);
			for (i=0;i<mate_n;++i) {
				fputs(mate_fq[i].id,mate_fout[i]);
				fputs(mate_fq[i].seq,mate_fout[i]);
				fputc('\n',mate_fout[i]);
				fputs(mate_fq[i].com,mate_fout[i]);
				fputs(mate_fq[i].qual,mate_fout[i]);
				fputc('\n',mate_fout[i]);
			}

		}
	}
	fprintf(fstat, "Total reads: %d\n", nrec);
	fprintf(fstat, "Too short after clip: %d\n", ntooshort);
	fprintf(fstat, "Clipped reads: %d\n", ntrim);
	fprintf(fstat, "Errors: %d\n", nerr);
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

void usage(FILE *f) {
	fprintf(f, 
"usage: fastq-mcf [options] <adapters.fa> <reads.fq> [mates1.fq ...] \n"
"\n"
"Detects levels of adapter presence, computes likelihoods and\n"
"locations (start, end) of the adapters.   Removes the adapter\n"
"sequences from the fastq file, and sends it to stdout.\n\n"
"Stats go to stderr, unless -o is specified.\n\n"
"If you specify multiple 'paired-end' inputs, then a -o option is\n" 
"required for each.  IE: -o read1.clip.q -o read2.clip.fq\n"
"\n"
"Options:\n"
"	-h	This help\n"
"	-o FIL	Output file (stats to stdout)\n"
"	-s N.N	Log scale for clip pct to threshold (2.5)\n"
"	-t N	%% occurance threshold before clipping (0.25)\n"
"	-m N	Minimum clip length, overrides scaled auto (4)\n"
"	-p N	Maximum adapter difference percentage (20)\n"
"	-l N	Minimum remaining sequence length (15)\n"
"	-n	Don't clip, just output what would be done\n"
	);
}

/*
#!/usr/bin/perl

my ($f, $a) = @ARGV;

my @a = split(/,/, $a);

open (F, $f) || die;

while (my $r = <F>) {
	for my $a (@a) {
		for (my $i = 1; $i < length($a); ++$i) {
			
		}
	}
}
# http://www.perlmonks.org/?node_id=500235
sub hd{ length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }
*/
