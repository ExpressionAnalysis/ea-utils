#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string>
#include <vector>

#include <sys/time.h>
#include <unistd.h>

#include <sparsehash/dense_hash_map>

#include "fastq-lib.h"
#include "utils.h"
#include "tidx.h"

void usage(FILE *f);

using namespace std;
using namespace google;

int main (int argc, char **argv) {
    bool debug = false;
    bool echo = false;
    bool build = false;
    vector<const char *> vin;
    const char *ain= NULL;
	const char *sep = "\t";
    const char *msep = "^";
    const char *trim = "chr";
    char *point = NULL;
    int nchr = 1, nbeg = 2, nend = 3;
    char skip_c = '#';
    char sub_e = 0;
    int skip_i = 0;
    bool dump = 0;

    char c;
    while ( (c = getopt (argc, argv, "Dlhdt:r:c:b:T:e:p:i:s:a:nB")) != -1) {
        switch (c) {
            case 'd':
                debug=true; break;
            case 'D':
                dump=true; break;
            case 'n':
                echo=false; break;
            case 'l':
                sub_e=true; break;
            case 'B':
                build = true; break;
            case 'h':
                usage(stdout); exit(0);
            case 't':
                sep = optarg; break;
            case 'p':
                point = optarg; break;
            case 's':
                if (isdigit(*optarg))
                    skip_i = atoi(optarg);
                else
                    skip_c = *optarg; 
                break;
            case 'r':
                msep = optarg; break;
            case 'T':
                trim = optarg; break;
            case 'c':
                nchr = atoi(optarg); break;
            case 'i':
                vin.push_back(optarg); break;
            case 'a':
                ain = optarg; break;
            case 'b':
                nbeg = atoi(optarg); break;
            case 'e':
                nend = atoi(optarg); break;
            case '?':
                if (strchr("tncbe", optopt))
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

    if (!vin.size())  {
        if (argc>1) warn("Error: at least one -i index file is required\n"); usage(stderr); exit(1);
    }

    if (! build && ! ain && !point && !dump) { 
        fail("Error: one of -D -B, -p or -a is required\n");
    }

    if ((!!build + !!ain + !!point + !!dump) > 1) {
        warn("Error: only one of -B, -p or -a is allowed\n");
    }

    --nchr; --nbeg; --nend;

    if ( build ) {
        int f_i;
        for (f_i=0;f_i<vin.size();++f_i) {
            tidx x;
            if (debug) 
                x.debug=true;
            x.build(vin[f_i], sep, nchr, nbeg, nend, skip_i, skip_c, sub_e);
        }
    } else {
	    struct line l; meminit(l);
        int nl = 0;
        int f_i;
        vector<tidx *>vmap; vmap.resize(vin.size());
        for (f_i=0;f_i<vin.size();++f_i) { 
            vmap[f_i]=new tidx(vin[f_i]);
            if (debug) 
                vmap[f_i]->debug=true;
            if (dump) {
                vmap[f_i]->dump(stdout);
            }
        }
        if (dump) {
            exit(0);
        }

        if ( point ) {
            char * p = strchr(point, ':');
            if (!p) {
                fail("Error: -p requires chr:pos argument\n");
            }
            *p++ = '\0';
            long pos = atol(p);
            int found = 0;
            for (f_i=0;f_i<vin.size();++f_i) {
                string tmp = vmap[f_i]->lookup(point, pos, msep);
                if (tmp.size()) {
                    ++found;
                    fputs(tmp.c_str()+(found==1),stdout);                             // echo
                }
            }
            if (found) fputc('\n',stdout);
            return !found;
        } else {
            FILE *fin = !strcmp(ain,"-") ? stdin : fopen(ain, "r");
            if (!fin)
                fail("error '%s':%s", ain,strerror(errno));

            while (read_line(fin, l)>0) {
                ++nl;

                chomp_line(l);

                fputs(l.s,stdout);                              // echo

                vector<char *> v = split(l.s, sep);     // todo, only get the keys desired, don't destroy

                string res;
                if (v.size() > nchr && v.size() > nbeg) {
                    for (f_i=0;f_i<vin.size();++f_i) {
                        string tmp = vmap[f_i]->lookup(v[nchr], atol(v[nbeg]), msep);
                        if (tmp.size()) {
                            res = res + tmp;
                        }
                    }
                }
                fputs(res.c_str(),stdout);                             // echo
                fputc('\n',stdout);
            }
            free_line(&l);
        }

        for (f_i=0;f_i<vin.size();++f_i) 
            delete vmap[f_i];
    }
}

void usage(FILE *f) {
fputs(
"Usage: tidx [options] -i IFILE [-i IFILE2...] -a AFILE\n"
"   or: tidx [options] -B -i IFILE\n"
"\n"
"Fragments and merges overlapping regions in an file with start-stop values.\n"
"Creating a simple, fast, compressed index\n"
"\n"
"Also can load that index, and search AFILE for intersecting lines\n"
"\n"
"If the 'group' column is zero, no grouping will be used\n"
"\n"
"If just -b is present during a search, then only that column\n"
"is searched.\n"
"\n"
"If both -b and -e are present during a search, then all regions\n"
"that overlap will be returned.\n"
"\n"
"Options and (defaults):\n"
"\n"
"-i IFILE       Text file to index (can specify more than one)\n"
"-B             Build index, don't annotate\n"
"-a FILE        Read text file and annotate\n"
"-p CHR:POS     Lookup a single point (slow!)\n"
"-r STRING      Annotation response separator (^)\n"
"-t CHAR(s)     Field separator (TAB)\n"
"-c INT         Group by (chromosome) column (1)\n"
"-b INT         Begin region column (2) (or position for annot)\n"
"-e INT         End region column (3)\n"
"-s INT or CHAR Skip rows starting with CHAR (#), or skip INT rows\n"
"-l             Less than end, not less than or equal-to\n"
"-n             Don't echo input lines\n"
//"-d             Verbose debug output\n"
//"-D             Dump input table (debug)\n"
//"-p CHR:POS     Single point lookup (debug, slow)\n"
"\n"
        ,f);
}

