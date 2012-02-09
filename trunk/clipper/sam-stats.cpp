#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <sys/stat.h>

#include <string>
#include <google/sparse_hash_map> // or sparse_hash_set, dense_hash_map, ...
#include <google/dense_hash_map> // or sparse_hash_set, dense_hash_map, ...

#include <api/BamReader.h>
#include <api/BamWriter.h>

const char * VERSION = "1.2";

using namespace BamTools;
using namespace std;

void usage(FILE *f);

#define max(a,b) (a>b?a:b)
#define min(a,b) (a<b?a:b)
#define meminit(l) (memset(&l,0,sizeof(l)))
#define debugout(s,...) if (debug) fprintf(stderr,s,##__VA_ARGS__)
#define warn(s,...) ((++errs), fprintf(stderr,s,##__VA_ARGS__))
#define stdev(cnt, sum, ssq) sqrt((((double)cnt)*ssq-pow((double)sum,2)) / ((double)cnt*((double)cnt-1)))

double quantile(const std::vector<int> &vec, double p);
std::string string_format(const std::string &fmt, ...);

int debug=0;
int errs=0;
extern int optind;
int histnum=30;
/// if we use this a lot may want to make it variable size
class scoverage {
public:
	scoverage() {mapb=reflen=0; dist.resize(histnum+2);};
	long long int mapb;
	int reflen;
	vector <int> dist;
};

class sstats {
public:
	sstats() {
		memset((void*)&dat,0,sizeof(dat));
		covr.set_empty_key("-");
	}
	~sstats() {
		covr.clear();
	}
	struct {
		int n, mapn;		// # of entries, # of mapped entries, 
		int lenmin, lenmax; double lensum, lenssq;	// read length stats
		double mapsum, mapssq;	// map quality sum/ssq 
		double nmnz, nmsum;	// # of mismatched reads, sum of mismatch lengths 
		int nbase, qualmax, qualmin;	// num bases samples, min/max qual 
		double qualsum, qualssq;	// sum quals, sum-squared qual
		int nrev, nfor;		// rev reads, for reads
		double tmapb;		// number of mapped bases
		long long int basecnt[5];
		int del, ins;		// length total dels/ins found
		bool pe;		// paired-end ? 0 or 1	
		int dupmax;		// max dups found
	} dat;
	vector<int> vmapq;		// all map qualities
	vector<int> visize;		// all insert sizes
	google::dense_hash_map<std::string, scoverage> covr;	// # mapped per ref seq
	google::sparse_hash_map<std::string, int> dups;		// alignments by read-id (not necessary for some pipes)

	// file-format neutral ... called per read
	void dostats(string name, int rlen, int bits, const string &ref, int pos, int mapq, int nmate, const string &seq, const string &qual, int nm, int del, int ins);

	// read a bam/sam file and call dostats over and over
	bool parse_bam(const char *in);
	bool parse_sam(FILE *f);
};

class line {
public:
	line() {s=NULL; n=0; a=0;}
        char *s; int n; size_t a;
};
int read_line(FILE *in, line &l); 

#define T_A 0
#define T_C 1
#define T_G 2
#define T_T 3
#define T_N 4


int qualreads = 1000000;
int dupreads = 1000000;
bool trackdup=0;
int main(int argc, char **argv) {
	const char *ext = NULL;
	bool multi=0, newonly=0, inbam=0;
	char c;
	optind = 0;
    while ( (c = getopt (argc, argv, "?BDdx:MhH:")) != -1) {
                switch (c) {
                case 'd': ++debug; break;
                case 'D': ++trackdup; break;
                case 'B': inbam=1; break;
                case 'b': qualreads=atoi(optarg); break;
                case 'H': histnum=atoi(optarg); break;
                case 'x': ext=optarg; break;
                case 'M': newonly=1; break;
                case 'h': usage(stdout); return 0;
                case '?':
                     if (!optopt) {
                        usage(stdout); return 0;
                     } else if (optopt && strchr("ox", optopt))
                       fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                     else if (isprint(optopt))
                       fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                     else
                       fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                     usage(stderr);
                     return 1;
                }
    }

	// recompute argc owing to getopt
	const char *stdv[3] = {argv[0],"-",NULL}; 
	if (!argv[optind]) {
		argc=2;
		argv = (char **) stdv;
		optind=1;
	}

	multi = (argc-optind-1) > 0;

	if (multi && !ext) 
		ext = "stats";

	debugout("argc:%d, argv[1]:%s, multi:%d, ext:%s\n", argc,argv[optind],multi,ext);
	const char *p;
	for (;optind < argc;++optind) {
		sstats s;
		const char *in = argv[optind];
		FILE *f;
		FILE *o=NULL;
		bool needpclose = 0;
		if (!strcmp(in,"-")) {
			// read sam/bam from stdin
			if (ext) {
				warn("Can't use file extension with stdin\n");
				continue;
			}
			f = stdin;
			o = stdout;
		} else {
			string out;
			if ((p = strrchr(in,'.')) && !strcmp(p, ".gz")) {
				// maybe this is a gzipped sam file...
				string cmd = string_format("gunzip -c '%s'", in);
				f = popen(cmd.c_str(), "r");
				needpclose=1;
				if (f) {
					char c;
					if (!inbam) {
						// guess file format with 1 char
						c=getc(f); ungetc(c,f);
						if (c==-1) {
							warn("Can't unzip %s\n", in);
							pclose(f);
							continue;
						}
						if (c==31) {
							// bam file... reopen to reset stream... can't pass directly
							string cmd = string_format("gunzip -c '%s'", in);
							f = popen(cmd.c_str(), "r");
							inbam=1;
						}
					} else 
						c = 31;	// user forced bam, no need to check/reopen

					if (inbam) {
						// why did you gzip a bam... weird? 
						if (dup2(fileno(f),0) == -1) {
						      warn("Can't dup2 STDIN\n");
						      continue;
						}
						in = "-";
					}
				} else {
					warn("Can't unzip %s: %s\n", in, strerror(errno));
					continue;
				}
				// extension mode... output to file minus .gz
				if (ext) 
					out=string(in, p-in);
			} else {
	 			f = fopen(in, "r");
				if (!f) {
					warn("Can't open %s: %s\n", in, strerror(errno));
					continue;
				}
				// extension mode... output to file
				if (ext) 
					out=in;
			}
			if (ext) {
				( out += '.') += ext;
				o=fopen(out.c_str(), "w");
				if (!o) {
					warn("Can't write %s: %s\n", out.c_str(), strerror(errno));
					continue;
				}
			} else
				o=stdout;
		}

		debugout("file:%s, f: %lx\n", in, (long int) f);
		char c;
		if (!inbam) {
			// guess file format
			c=getc(f); ungetc(c,f);
			if (c==31 && !strcmp(in,"-")) {
				// if bamtools api allowed me to pass a stream, this wouldn't be an issue....
				warn("Specify -B to read a bam file from standard input\n");
				continue;
			}
		} else 
			c = 31;
		if (c != 31) {
			// (could be an uncompressed bam... but can't magic in 1 char)
			if (!s.parse_sam(f)) {
				if (needpclose) pclose(f); else fclose(f);
				warn("Invalid sam file %s\n", in);
				continue;
			}
		} else {
			if (!s.parse_bam(in)) {
				if (needpclose) pclose(f); else fclose(f);
				warn("Invalid bam file %s\n", in);
				continue;
			}
		}
		if (needpclose) pclose(f); else fclose(f);

		sort(s.vmapq.begin(), s.vmapq.end());
		sort(s.visize.begin(), s.visize.end());

		int phred = s.dat.qualmin < 64 ? 33 : 64;
		if (!s.dat.n) {
			warn("No reads in %s\n", in);
			continue;
		}
		fprintf(o, "reads\t%d\n", s.dat.n);
		fprintf(o, "version\t%s\n", VERSION);

		// mapped reads is the number of reads that mapped at least once (either mated or not)
		if (s.dat.mapn > 0) {
			if (trackdup && s.dat.dupmax > (s.dat.pe+1)) {
				google::sparse_hash_map<string,int>::iterator it = s.dups.begin();
				vector<int> vtmp;
				int amb = 0;
				int sing = 0;
				while(it!=s.dups.end()) {
					// *not* making the distinction between 2 singleton mappings and 1 paired here
					if (it->second > (s.dat.pe+1)) {
						++amb;
					}
					if (it->second == 1 && s.dat.pe) {
						++sing;	
					}
					++it;
				}
				fprintf(o,"mapped reads\t%d\n", (int) s.dups.size()*(s.dat.pe+1)-sing);
				if (amb > 0) {
					fprintf(o,"ambiguous\t%d\n", amb*(s.dat.pe+1));
					fprintf(o,"pct ambiguous\t%.6f\n", 100.0*((double)amb/(double)s.dups.size()));
					fprintf(o,"max dup align\t%.d\n", s.dat.dupmax-s.dat.pe);
				}
				if (sing)
					fprintf(o,"singleton mappings\t%.d\n", sing);
				// number of total mappings
				fprintf(o, "total mappings\t%d\n", s.dat.mapn);
			} else {
				// dup-id's not tracked
				fprintf(o, "mapped reads\t%d\n", s.dat.mapn);
				// todo: add support for bwa's multiple alignment tag
				// fprintf(o, "total mappings\t%d\n", s.dat.mapn);
			}
		} else {
			fprintf(o, "mapped reads\t%d\n", s.dat.mapn);
		}

		if (s.dat.pe) {
			fprintf(o, "library\tpaired-end\n");
		}

		if (s.dat.mapn > 0) {
			fprintf(o, "bsize\t%d\n", qualreads);
			fprintf(o, "phred\t%d\n", phred);
			fprintf(o, "forward\t%d\n", s.dat.nfor);
			fprintf(o, "reverse\t%d\n", s.dat.nrev);
			if (s.dat.lenmax != s.dat.lenmin) {
				fprintf(o, "len max\t%d\n", s.dat.lenmax);	
				fprintf(o, "len mean\t%.4f\n", s.dat.lensum/s.dat.mapn);	
				fprintf(o, "len stdev\t%.4f\n", stdev(s.dat.mapn, s.dat.lensum, s.dat.lenssq));	
			} else {
				fprintf(o, "len max\t%d\n", s.dat.lenmax);	
			}
			fprintf(o, "mapq mean\t%.4f\n", s.dat.mapsum/s.dat.mapn);
			fprintf(o, "mapq stdev\t%.4f\n", stdev(s.dat.mapn, s.dat.mapsum, s.dat.mapssq));
			fprintf(o, "mapq Q1\t%.2f\n", quantile(s.vmapq,.25));
			fprintf(o, "mapq median\t%.2f\n", quantile(s.vmapq,.50));
			fprintf(o, "mapq Q3\t%.2f\n", quantile(s.vmapq,.75));

			if (s.dat.lensum > 0) {
				fprintf(o, "snp rate\t%.6f\n", s.dat.nmsum/s.dat.lensum);
				if (s.dat.ins >0 ) fprintf(o, "ins rate\t%.6f\n", s.dat.ins/s.dat.lensum);
				if (s.dat.del >0 ) fprintf(o, "del rate\t%.6f\n", s.dat.del/s.dat.lensum);
				fprintf(o, "pct mismatch\t%.4f\n", 100.0*((double)s.dat.nmnz/s.dat.mapn));
			}

			if (s.visize.size() > 0) {
				double p10 = quantile(s.visize, .10);
				double p90 = quantile(s.visize, .90);
				double matsum=0, matssq=0;
				int matc = 0;
				int i;
				for(i=0;i<s.visize.size();++i) {
					int v = s.visize[i];
					if (v >= p10 && v <= p90) {
						++matc;
						matsum+=v;
						matssq+=v*v;
					}
				}
				fprintf(o, "insert mean\t%.4f\n", matsum/matc);
				if (matc > 1) {
					fprintf(o, "insert stdev\t%.4f\n", stdev(matc, matsum, matssq));
					fprintf(o, "insert Q1\t%.2f\n", quantile(s.visize, .25));
					fprintf(o, "insert median\t%.2f\n", quantile(s.visize, .50));
					fprintf(o, "insert Q3\t%.2f\n", quantile(s.visize, .75));
				}
			}

			if (s.dat.nbase >0) {
				fprintf(o,"base qual mean\t%.4f\n", (s.dat.qualsum/s.dat.nbase)-phred);
				fprintf(o,"base qual stdev\t%.4f\n", stdev(s.dat.nbase, s.dat.qualsum, s.dat.qualssq));
				fprintf(o,"%%A\t%.4f\n", 100.0*((double)s.dat.basecnt[T_A]/(double)s.dat.nbase));
				fprintf(o,"%%C\t%.4f\n", 100.0*((double)s.dat.basecnt[T_C]/(double)s.dat.nbase));
				fprintf(o,"%%G\t%.4f\n", 100.0*((double)s.dat.basecnt[T_G]/(double)s.dat.nbase));
				fprintf(o,"%%T\t%.4f\n", 100.0*((double)s.dat.basecnt[T_T]/(double)s.dat.nbase));
				if (s.dat.basecnt[T_N] > 0) {
					fprintf(o,"%%N\t%.4f\n", 100.0*((double)s.dat.basecnt[T_N]/(double)s.dat.nbase));
				}
			}
			// how many ref seqs have mapped bases?
			int mseq=0;
			google::dense_hash_map<string,scoverage>::iterator it = s.covr.begin();
			vector<string> vtmp;
			bool haverlen = 0;
			while (it != s.covr.end()) {
				if (it->second.mapb > 0) {
					++mseq;								// number of mapped refseqs
					if (mseq <= 1000) vtmp.push_back(it->first);		// don't bother if too many chrs
					if (it->second.reflen > 0) haverlen = 1;
				}
				++it;
			}
			// don't print per-seq percentages if size is huge, or is 1
			if ((haverlen || mseq > 1) && mseq <= 1000) {			// worth reporting
				// sort the id's
				sort(vtmp.begin(),vtmp.end());
				vector<string>::iterator vit=vtmp.begin();
				double logb=log(2);
				while (vit != vtmp.end()) {
					scoverage &v = s.covr[*vit];
					if (v.reflen && histnum > 0) {
						string sig;
						int d;
						for (d=0;d<histnum;++d) {
							sig += ('0' + (v.dist[d] ? (int) (log(v.dist[d])/logb) : 0));
						}
						fprintf(o,"%%%s\t%.2f\t%s\n", vit->c_str(), 100.0*((double)v.mapb/s.dat.lensum), sig.c_str());
					} else {
						fprintf(o,"%%%s\t%.2f\n", vit->c_str(), 100.0*((double)v.mapb/s.dat.lensum));
					}
					++vit;
				}
			}
			if (s.covr.size() > 1) {
				fprintf(o,"num ref seqs\t%d\n", (int) s.covr.size());
				fprintf(o,"num ref aligned\t%d\n", (int) mseq);
			}
		}
	}
	return errs ? 1 : 0;
}

#define S_ID 0
#define S_BITS 1
#define S_NMO 2
#define S_POS 3
#define S_MAPQ 4
#define S_CIG 5
#define S_MATE 8
#define S_READ 9
#define S_QUAL 10
#define S_TAG 11

void sstats::dostats(string name, int rlen, int bits, const string &ref, int pos, int mapq, int nmate, const string &seq, const string &qual, int nm, int del, int ins) {

	++dat.n;

	if (pos<=0) return;

	++dat.mapn;

	if (rlen > dat.lenmax) dat.lenmax = rlen;
	if ((rlen < dat.lenmin) || dat.lenmin==0) dat.lenmin = rlen;
	dat.lensum += rlen;
	dat.lenssq += rlen*rlen;

	if (bits & 16) 
		++dat.nrev;
	else
		++dat.nfor;

	dat.mapsum += mapq;
	dat.mapssq += mapq*mapq;

	vmapq.push_back(mapq);
	
	dat.nmsum += nm;
	if (nm > 0) dat.nmnz += 1;
	dat.del+=del;
	dat.ins+=ins;

	if (ref.length()) {
		scoverage *sc = &(covr[ref]);
		if (sc) {
			sc->mapb+=rlen;
			if (histnum > 0 && sc->reflen > 0) {
				int x = histnum * ((double)pos / sc->reflen);
				if (debug > 1) { 
					warn("chr: %s, hn: %d, pos: %d, rl: %d, x: %x\n", ref.c_str(), histnum, pos, sc->reflen, x);
				}
				if (x < histnum) {
	                                sc->dist[x]+=rlen;
				} else {
					// out of bounds.... what to do?
					sc->dist[histnum] +=rlen;
				}
			}	
		}
	}
	dat.tmapb+=rlen;
	if (nmate>0) {
		visize.push_back(nmate);
		dat.pe=1;
	}

	if (dat.mapn <= qualreads) {
		int i, j;
		for (i=0;i<qual.length();++i) {
			if (qual[i]>dat.qualmax) dat.qualmax=qual[i];
			if (qual[i]<dat.qualmin) dat.qualmin=qual[i];
			dat.qualsum+=qual[i];
			dat.qualssq+=qual[i]*qual[i];
			switch(seq[i]) {
				case 'A': case 'a':
					j=T_A; break;
				case 'C': case 'c':
					j=T_C; break;
				case 'G': case 'g':
					j=T_G; break;
				case 'T': case 't':
					j=T_T; break;
				default:
					j=T_N; break;
			}
			++dat.basecnt[j];
			++dat.nbase;
		}
	}
	if (trackdup) {
		size_t p;
		if ((p = name.find_first_of('/'))!=string::npos) 
				name.resize(p);
		int x=++dups[name];
		if (x>dat.dupmax) 
			dat.dupmax=x;
	}
}

bool sstats::parse_sam(FILE *f) {
	line l;
	while (read_line(f, l)>0)  {
		char *sp;
		if (l.s[0]=='@') {
			if (!strncmp(l.s,"@SQ\t",4)) {
				char *t=strtok_r(l.s, "\t", &sp);
				string sname; int slen=0;
				while(t) {
					if (!strncmp(t,"SN:",3)) {
						sname=&(t[3]);
						if (slen) 
							break;
					} else if (!strncmp(t,"LN:",3)) {
						slen=atoi(&t[3]);
						if (sname.length()) 
							break;
					}
					t=strtok_r(NULL, "\t", &sp);
				}
				covr[sname].reflen=slen;
			}
			continue;
		}
		char *t=strtok_r(l.s, "\t", &sp);
		char *d[100]; meminit(d);
		int n =0;
		while(t) {
			d[n++]=t;
			t=strtok_r(NULL, "\t", &sp);
		}
		int nm=0;
		int i;
		// get # mismatches
		for (i=S_TAG;i<n;++i){
			if (d[i] && !strncmp(d[i],"NM:i:",5)) {
				nm=atoi(&d[i][5]);
			}
		}

		if (!d[S_BITS] || !isdigit(d[S_BITS][0]) 
		 || !d[S_POS]  || !isdigit(d[S_POS][0])
		   ) { 
			// invalid sam
			return false;
		}

		int ins = 0, del = 0;	
		char *p=d[S_CIG];
		// sum the cig
		while (*p) {
			int n=strtod(p, &sp);
			if (sp==p) {
				break;
			}
			if (*sp == 'I') 
				ins+=n;
			else if (*sp == 'D') 
				del+=n;
			p=sp;
		}
		if (d[S_CIG][0] == '*') d[S_POS] = (char *) (char *) (char *) (char *) (char *) (char *) (char *) (char *) (char *) "-1";
		dostats(d[S_ID],strlen(d[S_READ]),atoi(d[S_BITS]),d[S_NMO],atoi(d[S_POS]),atoi(d[S_MAPQ]),atoi(d[S_MATE]),d[S_READ],d[S_QUAL],nm, ins, del);
	}
	return true;
}

bool sstats::parse_bam(const char *in) {
        BamReader inbam;
        if ( !inbam.Open(in) ) {
                warn("Error reading '%s': %s\n", in, strerror(errno));
                return false;
        }
	SamHeader theader = inbam.GetHeader();
	RefVector references = inbam.GetReferenceData();
	int i;
	for (i = 0; i < references.size(); ++i) {
		covr[references[i].RefName].reflen=references[i].RefLength;
	}
	BamAlignment al;
        while ( inbam.GetNextAlignment(al) ) {
		int nm;
		al.GetTag("NM",nm);
		int ins=0, del=0;
		int i;
		for (i=0;i<al.CigarData.size();++i) {
			if (al.CigarData[i].Type=='I') {
				ins+=al.CigarData[i].Length;
			} else if (al.CigarData[i].Type=='D') {
				del+=al.CigarData[i].Length;
			}
		}
		if (al.CigarData.size() == 0) {
			al.Position=-1;
		}
		dostats(al.Name,al.Length,al.AlignmentFlag,al.RefID>=0?references.at(al.RefID).RefName:"",al.Position+1,al.MapQuality, al.InsertSize, al.QueryBases, al.Qualities, nm, ins, del);
	}
	return true;
}

int read_line(FILE *in, line &l) {
        return (l.n = getline(&l.s, &l.a, in));
}

void usage(FILE *f) {
        fprintf(f,
"Usage: sam-stats [options] [file1] [file2...filen]\n"
"Version: %s\n"
"\n"
"Produces lots of easily digested statistics for the files listed\n"
"\n"
"Options:\n"
"\n"
"-D             Keep track of multiple alignments (slower!)\n"
"-M             Only overwrite if newer (requires -x, or multiple files)\n"
"-B             Input is bam, don't bother looking at magic\n"
"-x FIL         File extension for multiple files (stats)\n"
"-b INT         Number of reads to sample for per-base stats (1M)\n"
"-S INT         Size of ascii-signature (30)\n"
"\n"
"OUTPUT:\n"
"\n"
"If one file is specified, then the output is to standard out.  If\n"
"multiple files are specified, or if the -x option is supplied,\n"
"the output file is <filename>.<ext>.  Default extension is 'stats'.\n"
"\n"
"Complete Stats:\n"
"\n"
"  <STATS>           : mean, max, stdev, median, Q1 (25 percentile), Q3\n"
"  reads             : # of entries in the file\n"
"  phred             : phred scale used\n"
"  mapped reads      : number of aligned reads\n"
"  mapped bases      : total of the lengths of the aligned reads\n"
"  forward           : number of forward-aligned reads\n"
"  reverse           : number of reverse-aligned reads\n"
"  snp rate          : mismatched bases / total bases\n"
"  ins rate          : insert bases / total bases\n"
"  del rate          : deleted bases / total bases\n"
"  pct mismatch      : percent of reads that have mismatches\n"
"  len <STATS>       : read length stats, ignored if fixed-length\n"
"  mapq <STATS>      : stats for mapping qualities\n"
"  insert <STATS>    : stats for insert sizes\n"
"  %%<CHR>           : percentage of mapped bases per chr, followed by a signature\n"
"\n"
"Subsampled stats (1M reads max):\n"
"  base qual <STATS> : stats for base qualities\n"
"  %%A,%%T,%%C,%%G       : base percentages\n"
"\n"
"Meaning of the per-chromosome signature:\n"
"  A ascii-histogram of mapped reads by chromosome position.\n"
"  It is only output if the original SAM/BAM has a header. The values\n"
"  are the log2 of the # of mapped reads at each position + ascii '0'.\n"
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

double quantile(const std::vector<int> &vec, double p) {
        int l = vec.size();
        double t = ((double)l-1)*p;
        int it = (int) t;
        int v=vec[it];
        if (t > (double)it) {
                return (v + p * (vec[it+1] - v));
        } else {
                return v;
        }
}

