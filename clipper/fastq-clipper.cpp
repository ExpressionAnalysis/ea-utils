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

Replaced, largely, by fastq-mcf.

See "void usage" below for usage.

*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

#include "fastq-lib.h"

#define MAX_ADAPTER_NUM 20
#define MAX_ADAPTER_LEN 160

void usage(FILE *f);
int hd(char *a, char *b, int n);
int debug=0;
int main (int argc, char **argv) {
	char c;
	bool eol;
	int nmin = 4, nkeep = 15, xmax=-1, pctdiff = 20;
	char *outfile = NULL;
	
	int i;
	
	char *a = NULL, *f = NULL;
	while (	(c = getopt (argc, argv, "-hedbp:i:o:l:m:x::")) != -1) {
		switch (c) {
		case '\1': 
			if (!f) 
				f=optarg; 
			else if (!a) 
				a=optarg; 
			else {
				usage(stderr); return 1;
			}
			break;
		case 'm': nmin = atoi(optarg); break;
		case 'p': pctdiff = atoi(optarg); break;
		case 'l': nkeep = atoi(optarg); break;
		case 'e': eol = 1; break;
		case 'h': usage(stdout); return 1; 
		case 'b': eol = 0; break;
		case 'd': debug = 1; break;
		case 'x': xmax = optarg ? atoi(optarg) : -1; break;
		case 'o': outfile = optarg; break;
		case 'i': f = optarg; break;
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

	if (argc < 3 || !a || !f) {
		usage(stderr);
		return 1;
	}

	FILE *fin = strcmp(f,"-") ? fopen(f, "r") : stdin; 
	if (!fin) {
		fprintf(stderr, "Error opening file '%s': %s\n",f, strerror(errno));
		return 1;
	}

	FILE *fout = stdout;
	FILE *fstat = stderr;
	if (outfile ) {
		fout = fopen(outfile, "w"); 
		if (!fout) {
			fprintf(stderr, "Error opening output file '%s': %s",outfile, strerror(errno));
			return 1;
		}
		fstat = stdout;
	}

	char *adapters[MAX_ADAPTER_NUM+1];
	int adapter_len[MAX_ADAPTER_NUM+1];
	char *p;
	int adapter_count=0;
	while (p=strtok(a,":")) {
		a = NULL;					// strtok requirement
                adapters[adapter_count] = p;
                adapter_len[adapter_count] = strlen(p);         // append to list
		++adapter_count;
                if (adapter_count >= MAX_ADAPTER_NUM) {
                        break;
                }
        }

	char *s[4] = {0,0,0,0}; 	// id, sequence, comment, quality
	size_t na[4] = {0,0,0,0};	// lengths of above
	int ns[4] = {0,0,0,0};	// lengths of above
	int nrec=0;
	int nerr=0;
	int nok=0;
	int ntooshort=0;
	int ntrim=0;
	int nbtrim=0;
	while (1) {
		int i;
		for (i = 0; i < 4; ++i ) {
			ns[i] = getline(&s[i], &na[i], fin);
		}

		if (ns[1] <= 0) { 
			break;
		}

		++nrec;

		// skip malformed records
		if (ns[1] != ns[3] || s[0][0] != '@' || s[2][0] != '+') {
			if (nerr < 10) {
				fprintf(stderr, "Malformed fastq record at line %d\n", nrec*4-3);
			}
			++nerr;
			continue;
		}

		// chomp
		s[1][ns[1]-1]='\0';
		--ns[1];
		s[3][ns[3]-1]='\0';
		--ns[3];

		if (debug) fprintf(stderr, "seq: %s %d\n", s[1], ns[1]);

		bool skip = 0;
		int bestscore = 999, bestoff = 0, bestlen = 0;

		for (i =0; i < adapter_count; ++i) {
			int nmatch = nmin;
			if (!nmatch) nmatch = adapter_len[i];			// full match required if nmin == 0
	
			// how far in to search for a match?
			int mx = adapter_len[i];
			if (xmax) {
				 mx = ns[1];
				 if (xmax > 0 && (xmax+adapter_len[i]) < mx)
					mx = xmax+adapter_len[i];		// xmax is added to adapter length
			}

			if (debug)
				fprintf(stderr, "adapter: %s, adlen: %d, nmatch: %d, mx: %d\n", adapters[i], adapter_len[i], nmatch, mx);

			int off;
			for (off = nmatch; off <= mx; ++off) {			// off is distance from tail of sequence
				char *seqtail = s[1]+ns[1]-off; 		// search at tail
				int ncmp = off<adapter_len[i] ? off : adapter_len[i];
				int mind = (pctdiff * ncmp) / 100;
				int d = hd(adapters[i],seqtail,ncmp);		// # differences
				if (debug)
					fprintf(stderr, "tail: %s, bestoff: %d, off: %d, ncmp: %d, mind: %d, hd %d\n", seqtail, bestoff, off, ncmp, mind, d);
				// calc squared distance over length score
				if (d <= mind) {
					int score = (d*d+1)/ncmp;
					if (score <= bestscore) {			// better score?
						bestscore = score;			// save max score
						bestoff = off;				// offset at max
						bestlen = ncmp;				// cmp length at max
					}
					if (d == 0 && (ncmp == adapter_len[i])) {
						break;
					}
				}
			}

			// assure time wasn't wasted running a comparison that couldn't matter
			assert((bestlen == 0) || (bestlen >= nmatch));

			if (bestoff > 0) {
				if ( (ns[1]-bestoff) < nkeep) {
					++ntooshort;
					skip = 1;
					break;
				}
			}
		}	

		if (!skip) {
			if (bestoff > 0) {
				++ntrim;
				s[1][ns[1]-bestoff]='\0';
				s[3][ns[1]-bestoff]='\0';
			}
			fputs(s[0],fout);
			fputs(s[1],fout);
			fputc('\n',fout);
			fputs(s[2],fout);
			fputs(s[3],fout);
			fputc('\n',fout);
		}
	}
	fprintf(fstat, "Total: %d\n", nrec);
	fprintf(fstat, "Too Short: %d\n", ntooshort);
	fprintf(fstat, "Trimmed: %d\n", ntrim);
	fprintf(fstat, "Errors: %d\n", nerr);
	return 0;
}

void usage(FILE *f) {
	fprintf(f, 
"usage: fastq-clipper [options] <fastq-file> <adapters>\n"
"\n"
"Removes one or more adapter sequences from the fastq file.\n"
"Adapter sequences are colon-delimited.\n"
"Stats go to stderr, unless -o is specified.\n"
"\n"
"Options:\n"
"	-h	This help\n"
"	-o FIL	Output file (stats to stdout)\n"
"	-p N	Maximum difference percentage (10)\n"
"	-m N	Minimum clip length (1)\n"
"	-l N	Minimum remaining sequence length (15)\n"
"	-x [N]	Extra match length past adapter length, \n"
"		 N =-1 : search all\n"
"		 N = 0 : search only up to adapter length\n"
"	-e	End-of-line (default)\n"
"	-b	Beginning-of-line (not supported yet)\n"
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
