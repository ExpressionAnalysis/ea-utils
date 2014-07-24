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

ALSO, IT WOULD BE NICE IF YOU LET ME KNOW YOU USED IT.

$Id$
*/
const char * VERSION = "2.0";

#include <ctype.h>
#include <stdio.h>

void usage(FILE *f) {
        fputs(
"Usage: bam-filter [-h] [-i IN] [-b BAD-READS] [-s STAT] -f FILTER1 [-f FITLER2...] -o OUT \n"
"Version: %s\n" 
"\n"
"	-h	help\n"
"	-i	input (stdin)\n"
"	-s	stats (stderr)\n"
"	-t	trim char (none)\n"
"	-o	output bam prefix\n"
"	-b      writes reads that were removed (bad reads) to this file, if specified\n"
"	-e	save all equal alignments\n"
"	-m	save equal only if edit distance is smaller\n"
"	-f	bam file to filter input against\n"
"\n"
"For each read in the input bam, searches all the filters.  If an analogous read in a filter has\n"
"a higher mapping quality, then the read from the input bam is discarded.\n"
"\n"
        ,f);
}

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include <sys/stat.h>

#include <string>
#include <sparsehash/sparse_hash_map> // or sparse_hash_set, dense_hash_map, ...

#include <api/BamReader.h>
#include <api/BamWriter.h>

using namespace BamTools;

class mapq {
public:
	int	mq;
	int	nm;
	mapq() {mq=0;nm=0;};
};

#define MAX_F 128

void usage(FILE *f);

int debug=0;
int main(int argc, char **argv) {
  char c, *in=NULL, *out=NULL, *err=NULL, *bad=NULL, trimchar = '\0';
	char *filter[MAX_F];
	int nfilter = 0, saveeq = 0, savenm = 0;


        while ( (c = getopt (argc, argv, "demo:i:s:f:b:t:")) != -1) {
                switch (c) {
                case 'd': ++debug; break;
                case 'b': bad=optarg; break;
                case 'o': out=optarg; break;
                case 'e': saveeq=1; break;
                case 'm': savenm=1; break;
                case 'i': in=optarg; break;
                case 't': trimchar=*optarg; break;
                case 's': err=optarg; break;
                case 'h': usage(stdout); return 0;
		case 'f': filter[nfilter++]=optarg;
			  if (nfilter >= MAX_F) {
                       		fprintf(stderr, "Too many filters: %d\n", nfilter);
                     		return 1;
			  }
			  break;
                case '?':
		     if (optopt == '?') {
                	usage(stdout); return 0;
                     } else if (optopt && strchr("oisft", optopt))
                       fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                     else if (isprint(optopt))
                       fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                     else
                       fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                     usage(stderr);
                     return 1;
		}
	}
	if (!nfilter) {
		usage(stderr);
		return 1;
	}

	if (!out) {
		usage(stderr);
		return 1;
	}

	FILE *ferr=stderr;
	if (err) {
		fprintf(stderr,"Opening %s\n",err);
		ferr=fopen(err,"w");
		if (!ferr) {
        		fprintf(stderr, "Error writing '%s': %s\n", err, strerror(errno));	
		}
	}

	// options are done
	BamReader inbam;
	if ( !inbam.Open(in) ) {
        	fprintf(stderr, "Error reading '%s': %s\n", in, strerror(errno));	
		return 1;
	}


	google::sparse_hash_map<std::string, mapq> pmap;
	int i;
	for (i=0;i<nfilter;++i) {
		BamReader fbam;
		if ( !fbam.Open(filter[i]) ) {
			fprintf(stderr, "Error reading '%s': %s\n", filter[i], strerror(errno));
			return 1;
		}
		if (debug) fprintf(stderr, "Indexing '%s'\n",filter[i]);
		google::sparse_hash_map<std::string,mapq>::iterator it; 
		BamAlignment al;
		while ( fbam.GetNextAlignment(al) ) {
			it = pmap.find(al.Name);

			mapq m;
			m.mq = al.MapQuality;
			al.GetTag("NM",m.nm);

			if (it == pmap.end()) {
				pmap[al.Name]=m;
			} else {
				if (al.MapQuality > it->second.mq) {
					pmap[al.Name]=m;
				}
			}
		}
	}

	SamHeader header = inbam.GetHeader();
	RefVector references = inbam.GetReferenceData();

	BamWriter writer;
	BamWriter badwriter;

	if ( !writer.Open(out, header, references) ) {
                fprintf(stderr, "Error writing '%s': %s", out, strerror(errno));
		return 1;
	}
	
	if ( !badwriter.Open(bad, header, references) ) {
                fprintf(stderr, "Error writing '%s': %s", bad, strerror(errno));
		return 1;
	}

	google::sparse_hash_map<std::string,mapq>::iterator it; 
	BamAlignment al;
	if (debug) fprintf(stderr, "Filtering\n");
	int na=0, gt=0, eq=0, lt=0, to=0, eq_r=0, eq_s=0;

	while ( inbam.GetNextAlignment(al) ) {
		++to;
		try {
			it = pmap.find(al.Name);
			if (it == pmap.end()) {
				// not found?
				writer.SaveAlignment(al);
				++na;
			} else if (al.MapQuality > it->second.mq) {
				// gt
				++gt;
				writer.SaveAlignment(al);
			} else if (al.MapQuality == it->second.mq) {
				// eq
				int nm;
				al.GetTag("NM",nm);
				if (nm < it->second.nm) {
					// fewer mismatches
					if (savenm) {
						writer.SaveAlignment(al);
						++eq_s;					
					} else {
					  badwriter.SaveAlignment(al);
						++eq_r;					
					}
				} else if (nm == it->second.nm) {
					if (saveeq) {
						writer.SaveAlignment(al);
						++eq_s;
					} else {
					  badwriter.SaveAlignment(al);
						++eq_r;					
					}
				} else {
				  badwriter.SaveAlignment(al);
					++eq_r;					
				}
			} else {
				// lt
			  badwriter.SaveAlignment(al);
				++lt;
			}
		} catch (...) {
		}
	}
	fprintf(ferr,"total\t%d\n",to);
	fprintf(ferr,"better\t%d\t%2.2f%%\n",gt,100.0*gt/(float)to);
	if (eq_s > 0) fprintf(ferr,"eq-saved\t%d\t%2.2f%%\n",eq_s,100.0*eq_s/(float)to);
	if (eq_r > 0) fprintf(ferr,"eq-removed\t%d\t%2.2f%%\n",eq_r,100.0*eq_r/(float)to);
	fprintf(ferr,"removed\t%d\t%2.2f%%\n",lt,100.0*lt/(float)to);
	if (na) fprintf(ferr,"missing\t%d\t%2.2f%%\n",na,100.0*na/(float)to);
	return 0;
}

