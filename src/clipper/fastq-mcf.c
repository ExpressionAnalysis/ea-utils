#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
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
	char bscan[SCANLEN+1]; 			// scan sequence
	int bcnt;			// number found at beginning
	int ecnt;			// number found at end
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
	int nmin = 4, nkeep = 15;
	float minpct = 0.25;
	int pctdiff = 20;
	char *outfile = NULL;
	int sampcnt = 20000;
	int xmax = -1;
	float scale = 2.5;
	
	int i;
	
	char *afil = NULL, *ifil = NULL;
	while (	(c = getopt (argc, argv, "-hp:o:l:m::")) != -1) {
		switch (c) {
		case '\1': 
			if (!afil) 
				afil = optarg; 
			else if (!ifil) 
				ifil = optarg; 
			else {
				usage(stderr); return 1;
			}
			break;
		case 'm': nmin = atoi(optarg); break;
		case 'l': nkeep = atoi(optarg); break;
		case 'h': usage(stdout); return 1; 
		case 'o': outfile = optarg; break;
		case 'i': ifil = optarg; break;
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
	if (outfile) {
		fout = fopen(outfile, "w"); 
		if (!fout) {
			fprintf(stderr, "Error opening output file '%s': %s",outfile, strerror(errno));
			return 1;
		}
		fstat = stdout;
	}

	struct ad ad[MAX_ADAPTER_NUM+1];
	memset(ad, 0, sizeof(*ad)*(MAX_ADAPTER_NUM+1));

	int acnt=0, ok=0, rno=0;	// adapter count, ok flag, record number
	while (acnt < MAX_ADAPTER_NUM && (ok = read_fa(ain, rno, &ad[acnt]))) {
		++rno;
		if (ok < 0)
			break;
		// copy in truncated to max scan length
		strncpy(ad[acnt].escan, ad[acnt].seq, SCANLEN);
		strncpy(ad[acnt].bscan, ad[acnt].seq+max(0,ad[acnt].nseq-SCANLEN), SCANLEN);
		ad[acnt].escan[SCANLEN] = '\0';
		ad[acnt].bscan[SCANLEN] = '\0';
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
			for(a=0;a<acnt;++a) {
				char *p;
				if (p = strstr(s+10, ad[a].escan)) { 
					printf("P: %d, SL: %d\n", p-s, ns);
					++ad[a].ecnt;
				}
				if ((p = strstr(s, ad[a].bscan)) && ((p-s) < (ns/2))) { 
					++ad[a].bcnt;
				}
			}
		}
		if (nr >= sampcnt) 
			break;
        }

	int a;
	int athr = (int) ((float)nr * minpct) / 100;
	fprintf(stderr, "Threshold used: %d out of %d\n", athr, nr);
	int newc=0;
	for(a=0;a<acnt;++a) {
		if (ad[a].ecnt > athr || ad[a].bcnt > athr) {
			int cnt;
			printf("EC: %d, BC:%d\n", ad[a].ecnt, ad[a].bcnt);
			if (ad[a].ecnt >= ad[a].bcnt) {
				ad[a].end='e';
				cnt = ad[a].ecnt;
			} else {
				ad[a].end='b';
				cnt = ad[a].bcnt;
			}
			
			ad[a].thr = max(1,(int) (-log(cnt / (float) nr)/log(scale)));
			fprintf(stderr, "Adapter %s (%s): counted %d at the '%s' of sequences, clip set to %d", ad[a].id, ad[a].seq, cnt, ad[a].end == 'e' ? "end" : "start", ad[a].thr);
			if (abs((ad[a].bcnt-ad[a].ecnt)) < athr/4) {
				fprintf(stderr, ", warning end was not reliable\n", ad[a].id, ad[a].seq);
			} else {
				fputc('\n', stderr);
			}

			ad[newc++]=ad[a];
		}
	}
	acnt=newc;

	struct fq fq;	
        memset(&fq, 0, sizeof(fq));

	int nrec=0;
	int nerr=0;
	int nok=0;
	int ntooshort=0;
	int ntrim=0;
	int nbtrim=0;
	int read_ok;
	fseek(fin, 0, 0);
	while (read_ok=read_fq(fin, nrec, &fq)) {
		++nrec;
		if (read_ok < 0) continue;

		// chomp

		if (debug) fprintf(stderr, "seq: %s %d\n", fq.seq, fq.nseq);

		bool skip = 0;
		int bestscore = 999, bestoff = 0, bestlen = 0;

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

			int off;
			for (off = nmatch; off <= mx; ++off) {			// off is distance from tail of sequence
				char *seqtail = fq.seq+fq.nseq-off; 		// search at tail
				int ncmp = off<ad[i].nseq ? off : ad[i].nseq;
				int mind = (pctdiff * ncmp) / 100;
				int d = hd(ad[i].seq,seqtail,ncmp);		// # differences
				if (debug)
					fprintf(stderr, "tail: %s, bestoff: %d, off: %d, ncmp: %d, mind: %d, hd %d\n", seqtail, bestoff, off, ncmp, mind, d);
				if (d <= mind) {
					int score = (d*d+1)/ncmp;
					if (score <= bestscore) {			// better score?
						bestscore = score;			// save max score
						bestoff = off;				// offset at max
						bestlen = ncmp;				// cmp length at max
					}
					if (d == 0 && (ncmp == ad[i].nseq)) {
						break;
					}
				}
			}

			if (bestoff > 0) {
				if ( (fq.nseq-bestoff) < nkeep) {
					++ntooshort;
					skip = 1;
					break;
				}
			}
		}	

		if (!skip) {
			if (bestoff > 0) {
				++ntrim;
				fq.seq[fq.nseq-bestoff]='\0';
				fq.qual[fq.nseq-bestoff]='\0';
			}
			fputs(fq.id,fout);
			fputs(fq.seq,fout);
			fputc('\n',fout);
			fputs(fq.com,fout);
			fputs(fq.qual,fout);
			fputc('\n',fout);
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
	strcpy(fa->id, p);
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
"usage: fastq-mcf [options] <adapters.fa> <reads.fq>\n"
"\n"
"Detects levels of adapter presence, computes likelihoods and \n"
"locations of the adapters.   Removes the adapter sequences from\n"
"the fastq file, and sends it to stdout.\n"
"Stats go to stderr, unless -o is specified.\n"
"\n"
"Options:\n"
"	-h	This help\n"
"	-o FIL	Output file (stats to stdout)\n"
"	-p N	Maximum difference percentage (20)\n"
"	-m N	Minimum clip length (4)\n"
"	-l N	Minimum remaining sequence length (15)\n"
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
