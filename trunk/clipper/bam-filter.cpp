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
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include <sys/stat.h>

#include <string>
#include <google/sparse_hash_map> // or sparse_hash_set, dense_hash_map, ...

#include <api/BamReader.h>
#include <api/BamWriter.h>

using namespace BamTools;

#define MAX_F 128

void usage(FILE *f);

int debug=0;
int main(int argc, char **argv) {
	char c, *in=NULL, *out=NULL, *err=NULL, trimchar = '\0';
	char *filter[MAX_F];
	int nfilter = 0;
        while ( (c = getopt (argc, argv, "do:i:s:f:t:")) != -1) {
                switch (c) {
                case 'd': ++debug; break;
                case 'o': out=optarg; break;
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
                     } else if (optopt && strchr("oisf", optopt))
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

	google::sparse_hash_map<std::string, int> pmap;
	int i;
	for (i=0;i<nfilter;++i) {
		BamReader fbam;
		if ( !fbam.Open(filter[i]) ) {
			fprintf(stderr, "Error reading '%s': %s\n", filter[i], strerror(errno));
			return 1;
		}
		if (debug) fprintf(stderr, "Indexing '%s'\n",filter[i]);
		google::sparse_hash_map<std::string,int>::iterator it; 
		BamAlignment al;
		while ( fbam.GetNextAlignment(al) ) {
			it = pmap.find(al.Name);
			if (it == pmap.end()) {
				pmap[al.Name]=al.MapQuality;
			} else {
				if (al.MapQuality > it->second) {
					pmap[al.Name]=al.MapQuality;
				}
			}
		}
	}

	SamHeader header = inbam.GetHeader();
	RefVector references = inbam.GetReferenceData();

	BamWriter writer;
	if ( !writer.Open(out, header, references) ) {
                fprintf(stderr, "Error writing '%s': %s", out, strerror(errno));
		return 1;
	}

	google::sparse_hash_map<std::string,int>::iterator it; 
	BamAlignment al;
	if (debug) fprintf(stderr, "Filtering\n");
	int na=0, gt=0, eq=0, lt=0, to=0;
	while ( inbam.GetNextAlignment(al) ) {
		++to;
		try {
			it = pmap.find(al.Name);
			if (it == pmap.end()) {
				// not found?
				writer.SaveAlignment(al);
				++na;
			} else if (al.MapQuality > it->second) {
				// gt
				++gt;
				writer.SaveAlignment(al);
			} else if (al.MapQuality == it->second) {
				// eq
				++eq;
				writer.SaveAlignment(al);
			} else {
				// lt
				++lt;
			}
		} catch (...) {
		}
	}
	fprintf(ferr,"total\t%d\n",to);
	fprintf(ferr,"better\t%d\t%2.2f%%\n",gt,100.0*gt/(float)to);
	fprintf(ferr,"equal\t%d\t%2.2f%%\n",eq,100.0*eq/(float)to);
	fprintf(ferr,"removed\t%d\t%2.2f%%\n",lt,100.0*lt/(float)to);
	if (na) fprintf(ferr,"Missing\t%d\t%2.2f%%\n",na,100.0*na/(float)to);
	return 0;
}

void usage(FILE *f) {
        fputs(
"Usage: bam-filter [-h] [-i IN] [-o OUT] [-s STAT] -f FILTER1 [-f FITLER2...]\n"
"\n"
"	-h	help\n"
"	-i	input (stdin)\n"
"	-o	output (stout)\n"
"	-s	stats (stderr)\n"
"	-t	trim char (none)\n"
"	-f	bam file to filter input against\n"
"\n"
"For each probe id in the input bam, searches all the filters.  If any have\n"
"a higher mapping quality, then the alignment is discared.\n"
"\n"
        ,f);
}
