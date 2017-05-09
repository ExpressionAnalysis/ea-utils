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
#include <search.h>
#include <limits.h>
#include <sparsehash/sparse_hash_map> // or sparse_hash_set, dense_hash_map, ...
#include <string>

typedef struct line {
        char *s; int n; size_t a;
} line;

struct fq {
        line id;
        line seq;
        line com;
        line qual;
};

class to_merge {
public:
	std::string sseq;
	int cnt;
	to_merge(std::string s, int c) {sseq=s; cnt=c;};
	to_merge() {};
};

class ent {
public:
	std::string seq;
	int cnt;

	ent(const std::string &s, const int &c) { seq=s; cnt=c; };

	static bool comp_cnt (const ent &a, const ent &b) {
		return a.cnt < b.cnt;
	};
};
class idseq {
public:
        std::string id;
        std::string seq;

        idseq(const std::string &s, const std::string &i) { seq=s; id=i; };
};

// get file extension
const char *fext(const char *f);
FILE *gzopen(const char *f, const char *m, bool*isgz);

int read_line(FILE *in, struct line &l);                // 0=done, 1=ok, -1=err+continue
int read_fq(FILE *in, int &lno, struct fq *fq);          // 0=done, 1=ok, -1=err+continue
int read_fa(FILE *in, int &lno, struct fq *fq);          // 0=done, 1=ok, -1=err+continue
#ifdef _WIN32
	ssize_t getline(char **lineptr, size_t *n, FILE *stream);
#endif
void usage(FILE *f);
double quantile(std::vector<int> vec, double p);
#define meminit(l) (memset(&l,0,sizeof(l)))
FILE *openordie(const char *nam, const char * mode, FILE *def, const char *errstr, bool *isgz=NULL);

std::string string_format(const std::string &fmt, ...);
bool file_newer(const char *f1, const char *f2);
std::vector<std::string> split(char* str,const char* delim);

#define MAX_EX 4

int main (int argc, char **argv) {
	char *in = NULL;
	char *out = NULL;
	char *stat = NULL;
	char *ref = NULL;
	char *pat = NULL;
	char *vexcl[MAX_EX] = {NULL};

	int excl_n = 0;
	int targ = 1000;
	int thr = -1;
	int mergs = 1;
	int mergc = 0;
	int mergn = 0;
	bool debug = false;

	char c;
	while ( (c = getopt (argc, argv, "-hdi:o:s:p:r:n:t:x:m:")) != -1) {
		switch (c) {
			case '\1':
				if (!in)
					in=optarg;
				else {
					fprintf(stderr, "Unknown parameter '%s'.\n", optarg);
					exit(1);
				}
			case 'i':
				in = optarg; break;
			case 'o':
				out = optarg; break;
			case 's':
				stat = optarg; break;
			case 'm':
				mergs = atoi(optarg); break;
			case 'r':
				ref = optarg; break;
			case 'p':
				pat = optarg; break;
			case 'x':
				vexcl[excl_n++] = optarg; break;
			case 'n':
				targ = atoi(optarg); break;
			case 't':
				thr = atoi(optarg); break;
			case 'd':
				debug=true; break;
			case 'h':
				usage(stdout);
				exit(0);
			case '?':
				if (strchr("iosrn", optopt))
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

	FILE *fin, *fout, *fstat, *fref = NULL, *fpat=NULL, *ftmp=NULL;
	bool fingz;
	fin = openordie(in, "r", stdin, "Error opening file '%s': %s\n", &fingz);
	fout = openordie(out, "w", stdout, "Error writing to file '%s': %s\n");
	fstat = openordie(stat, "w", stderr, "Error writing to file '%s': %s\n");

	if (ref) 
		fref = openordie(ref, "r", NULL, "Error opening file '%s': %s\n");
	
	if (pat) 
		fpat = openordie(pat, "r", NULL, "Error opening file '%s': %s\n");

	google::sparse_hash_map<std::string, std::string> mmap;


	int read_ok, nref=0, nrec = 0, lno = 0, npat=0; struct fq fq; meminit(fq);

	if (fref) {
		// read fasta referernce
		lno = 0;
		while (read_ok=read_fa(fref, lno, &fq)) {
			++nref;
			char * p=strchr(fq.id.s, ' ');
			if (p) *p='\0';
			mmap[fq.seq.s] = fq.id.s+1;
		}
		fprintf(fstat,"nrefseq\t%d\n", nref);
	}

	std::vector<idseq> pvec;
	if (fpat) {
		// patterns to search and exclude
		lno = 0;
		while (read_ok=read_fa(fpat, lno, &fq)) {
			++npat;
			char * p=strchr(fq.id.s, ' ');
			if (p) *p='\0';
			idseq np(fq.seq.s,fq.id.s+1);
			pvec.push_back(np);
		}
		fprintf(fstat,"npatseq\t%d\n", npat);
	}


	// map all sequences to the hash
	google::sparse_hash_map<std::string, int> pmap;
	lno = 0;
	while (read_ok=read_fq(fin, lno, &fq)) {
		++nrec;
//		fprintf(stderr, "read %d '%s'\n", nrec, fq.seq.s);
		if (fq.seq.n > 1) {
			++pmap[fq.seq.s];
		}
	}
	if (thr<0) thr = (int)(log(1+nrec)/log(10));

	std::vector<int> vec;
	std::vector<int> lvec;
	std::vector<ent> lis;

	std::string tmp;
	if (excl_n) {
		// make a temp file
		if (in) {
			tmp=string_format("%s.tmp.fq", in);
		} else {
			tmp=string_format("/tmp/mirna-quant-%d.tmp.fq", getpid());
		}
        ftmp = openordie(tmp.c_str(), "w", NULL, "Error opening file '%s': %s\n");
	}

	if (npat > 0 || excl_n) {
		// fill mmap with pattern matches, and fill tmp fastq with sequences
		std::vector<idseq>::iterator pit;
      	google::sparse_hash_map<std::string,int>::iterator it = pmap.begin();
		while (it != pmap.end()) {
			if (excl_n) {
				fputs("@\n", ftmp);
				fputs(it->first.c_str(),ftmp);
				fputs("\n+\n", ftmp);
				int i;
				for(i=0;i<it->first.size();++i) {
					fputc('h',ftmp);
				}
				fputc('\n', ftmp);
			}
			if (npat > 0) {
				for (pit=pvec.begin(); pit != pvec.end(); ++pit) {
					if (it->first.find(pit->seq) != std::string::npos) {
						mmap[it->first]=pit->id;
						break;
					}
				}
			}
			++it;
		}
	}

	if (excl_n) {
		int i;
		for (i=0;i<excl_n;++i) {
			char *excl=vexcl[i];
			std::string ebwt = string_format("%s.1.ebwt", excl);
			std::string cmd;
			int ret;
			if (file_newer(excl,ebwt.c_str())) {
				cmd=string_format("bowtie-build %s %s > /dev/null", excl, excl).c_str();
				fprintf(stderr,"+%s\n",cmd.c_str());
				if (ret=system(cmd.c_str())) {
					exit(ret >> 8);
				}
			}
			cmd = string_format("bowtie --sam-nohead -S %s %s 2> /dev/null", excl, tmp.c_str());
			fprintf(stderr,"+%s\n",cmd.c_str());
			FILE *aln;
			if (!(aln=popen(cmd.c_str(),"r"))) {
				fprintf(stderr, "Can't run bowtie: $!\n", strerror(errno));
				exit(1);
			}
			struct line l; meminit(l);
			while (read_line(aln, l)>0) {
				std::vector<std::string> v = split(l.s,"\t");
				if (v[9].size() > 17) {
					if (atoi(v[3].c_str()) > 0) {
						// add to mmap
						mmap[v[9]]=v[2];
					}
				}
			}
			unlink(tmp.c_str());
		}
	}	

	if (mergs > 0) {	
		if (debug) fprintf(stderr, "merge\n");
		google::sparse_hash_map<std::string,int>::iterator it = pmap.begin();
		google::sparse_hash_map<std::string,std::string>::iterator mit; 
		std::string sseq;
		std::string oseq;
		std::vector<to_merge> merge_list;
		while (it != pmap.end()) {
			int i;
			for (i=1;i<=mergs;++i) {
				sseq = it->first;
				if (sseq.length() > i) {
					oseq=sseq;
					sseq.erase(sseq.length()-i);
					if (debug) fprintf(stderr, "merge %s %s\n", oseq.c_str(), sseq.c_str());

					// search for "close enough" ref seq
					mit = mmap.find(sseq);
					// merge named sequences
					if (mit != mmap.end()) {
						// not safe to search map while iterating!
						if (debug) fprintf(stderr, "found: %d+%d\n", it->second, (int)pmap[oseq]);
						++mergn;
						mergc+=it->second;
						merge_list.push_back(to_merge(sseq, it->second));
						break;
					}
				}
			}
			++it;
		}
		int i;
		for (i=0;i<merge_list.size();++i) {
			pmap[merge_list[i].sseq]+=merge_list[i].cnt;
		}
	}

	int tot = 0;
	// build lists
	google::sparse_hash_map<std::string,int>::iterator it = pmap.begin();
	while (it != pmap.end()) {
 		ent e(it->first,it->second);
		lvec.push_back(it->first.size());
		if (it->second >= thr) {
			vec.push_back(it->second);
			tot+=it->second;
		}
		lis.push_back(e);
		++it;
	}
	
	std::sort(vec.begin(), vec.end());
	std::sort(lvec.begin(), lvec.end());
	std::sort(lis.begin(), lis.end(), ent::comp_cnt);

	fprintf(fstat, "reads\t%d\n", nrec);
	fprintf(fstat, "threshold\t%d\n", thr);
	fprintf(fstat, "pass num\t%d\n", (int)vec.size());
	fprintf(fstat, "pass tot\t%d\n", tot);

	if (mergc>0) fprintf(fstat, "merge tot\t%d\n", mergc);
	if (mergn>0) fprintf(fstat, "merge num\t%d\n", mergn);

	if (vec.size() == 0) {
		exit(0);
	}

	int q1=(int)quantile(vec,.25);
	int q2=(int)quantile(vec,.50);
	double q3=quantile(vec,.75);
	int q4=(int)quantile(vec,1);

	int lq0=(int)quantile(lvec,0);
	int lq1=(int)quantile(lvec,.25);
	int lq2=(int)quantile(lvec,.50);
	int lq3=(int)quantile(lvec,.75);
	int lq4=(int)quantile(lvec,1);

	double norm = ((double)targ)/q3;
	if (targ)
		fprintf(fstat, "q3 target\t%d\n", targ);

	if (mergs)
		fprintf(fstat, "allow mismatch\t%d\n", mergs);

	fprintf(fstat, "cnt q1\t%d\n", q1);
	fprintf(fstat, "cnt med\t%d\n", q2);
	fprintf(fstat, "cnt q3\t%d\n", (int) q3);
	fprintf(fstat, "cnt max\t%d\n", q4);

	fprintf(fstat, "len min\t%d\n", lq0);
	fprintf(fstat, "len q1\t%d\n", lq1);
	fprintf(fstat, "len med\t%d\n", lq2);
	fprintf(fstat, "len q3\t%d\n", lq3);
	fprintf(fstat, "len max\t%d\n", lq4);

	int i = 0;
	std::vector<ent>::iterator vit;
	for (vit=(lis.end()-1); vit>=lis.begin();--vit) {
		++i;
		std::string anno;
                google::sparse_hash_map<std::string,std::string>::iterator mit = mmap.find(vit->seq);
                if (mit != mmap.end()) {
                        anno=mit->second;
                }
		fprintf(fstat, "top %d\t%d\t%s\t%s\n", i, vit->cnt, vit->seq.c_str(), anno.c_str());
		if (i >= 10) break;	
	}

        google::sparse_hash_map<std::string,std::string>::iterator fn;
	std::string seq;
	int ncnt, idf;
	for (vit=(lis.end()-1); vit>=lis.begin();--vit) {
		if (vit->cnt <= thr) {
			break;
		}
		seq=vit->seq;

		std::string anno;
		// merge map
		google::sparse_hash_map<std::string,std::string>::iterator mit = mmap.find(seq);
		if (mit != mmap.end()) {
			anno=mit->second;
		}
		ncnt = (int) (norm * (double) vit->cnt);
		fprintf(fout, "%s\t%d\t%d\t%s\n", seq.c_str(), vit->cnt, ncnt, anno.c_str());
	}

}

////////////// pasted library stuff ////////////////

double quantile(std::vector<int> vec, double p) {
	int l = vec.size();
	double t = ((double)l-1)*p;
	int it = (int) t;
	int v=vec[it];
	if (t > (double)it) {
		return (v + (t-it) * (vec[it+1] - v));
	} else {
		return v;
	}
}

FILE *openordie(const char *nam, const char * mode, FILE *def, const char *errstr, bool *isgz) {
	FILE *r;
	if (!nam && def) r = def;
	else {
			if (isgz) {
				r = gzopen(nam, mode, isgz);
			} else {
				r = fopen(nam, mode);
			}
			if (!r) {
					fprintf(stderr, errstr,nam, strerror(errno));
					exit(1);
			}
	}
	return r;
}


#define comp(c) ((c)=='A'?'T':(c)=='a'?'t':(c)=='C'?'G':(c)=='c'?'g':(c)=='G'?'C':(c)=='g'?'c':(c)=='T'?'A':(c)=='t'?'a':(c))
void revcomp(struct fq *d, struct fq *s) {
        if (!d->seq.s) {
                d->seq.s=(char *) malloc(d->seq.a=s->seq.n);
                d->qual.s=(char *) malloc(d->qual.a=s->qual.n);
        } else if (d->seq.a < s->seq.n) {
                d->seq.s=(char *) realloc(s->seq.s, d->seq.a=s->seq.n);
                d->qual.s=(char *) realloc(s->qual.s, d->qual.a=s->qual.n);
        }
        int i;
        for (i=0;i<s->seq.n/2;++i) {
                char b=s->seq.s[i];
                char q=s->qual.s[i];
                //printf("%d: %c, %c\n", i, comp(s->seq.s[s->seq.n-i-1]), s->qual.s[s->qual.n-i-1]);
                d->seq.s[i]=comp(s->seq.s[s->seq.n-i-1]);
                d->qual.s[i]=s->qual.s[s->qual.n-i-1];
                //printf("%d: %c, %c\n", s->seq.n-i-1, comp(b), q);
                d->seq.s[s->seq.n-i-1]=comp(b);
                d->qual.s[s->seq.n-i-1]=q;
        }
        if (s->seq.n % 2) {
                //printf("%d: %c, %c\n", 1+s->seq.n/2, comp(s->seq.s[s->seq.n/2]));
                d->seq.s[s->seq.n/2] = comp(s->seq.s[s->seq.n/2]);
                d->qual.s[s->seq.n/2] = s->qual.s[s->seq.n/2];
        }
        d->seq.n=s->seq.n;
        d->qual.n=s->qual.n;
}

int read_line(FILE *in, struct line &l) {
	l.n = getline(&l.s, &l.a, in);
	if (l.n > 1 && l.s[l.n-1] == '\n' && l.s[l.n-2] == '\r') {
		--l.n;
		l.s[l.n-1] = '\n';
		l.s[l.n] = '\0';
	}
	return l.n;
}

int read_fa(FILE *in, int &lno, struct fq *fa) {
// note: this only reads one line of sequence!
        read_line(in, fa->id);
        read_line(in, fa->seq);
	lno+=2;
	if (fa->seq.n <= 0) 
		return 0;
        if (fa->id.s[0] != '>') {
                fprintf(stderr, "Malformed fasta record at line %d\n", lno+1);
                return -1;
        }
        fa->seq.s[--fa->seq.n] = '\0';
        fa->id.s[--fa->id.n] = '\0';
}

int read_fq(FILE *in, int &lno, struct fq *fq) {
        read_line(in, fq->id);
        read_line(in, fq->seq);
        read_line(in, fq->com);
        read_line(in, fq->qual);
	lno+=4;
        if (fq->qual.n <= 0)
                return 0;
        if (fq->id.s[0] != '@' || fq->com.s[0] != '+' || fq->seq.n != fq->qual.n) {
		if (fq->seq.n) fq->seq.s[--fq->seq.n] = '\0';
		if (fq->id.n) fq->id.s[--fq->id.n] = '\0';
		if (fq->com.n) fq->com.s[--fq->com.n] = '\0';
		if (fq->qual.n) fq->qual.s[--fq->qual.n] = '\0';
                fprintf(stderr, "Malformed fastq record at line %d\n", lno+1);
                return -1;
        }
        // chomp
        fq->seq.s[--fq->seq.n] = '\0';
        fq->qual.s[--fq->qual.n] = '\0';
        return 1;
}

void usage(FILE *f) {
        fputs(
"Usage: mirna-quant [options] [reads.fq]\n"
"\n"
"Computes counts of miRNA species.\n"
"\n"
"Options and (defaults):\n"
"\n"
"-o FIL         Output file (stdout)\n"
"-i FIL         Input file (stdin)\n"
"-s FIL         Stats file (stderr)\n"
"-r FIL         miRNA reference fasta (none)\n"
"-x FIL         Annot via alignment, ie: rRNA (none)\n"
"-n INT         Upper quartile normalize with target (1000)\n"
"-t INT         Threshold (log10(record count))\n"
"-m INT         Combine known (-r) seqs with up to INT mismatches\n"
"\n"
        ,f);
}

/* getline.c -- Replacement for GNU C library function getline

Copyright (C) 1993 Free Software Foundation, Inc.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

/* Written by Jan Brittenson, bson@gnu.ai.mit.edu.  */

#include <sys/types.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

/* Read up to (and including) a TERMINATOR from STREAM into *LINEPTR
   + OFFSET (and null-terminate it). *LINEPTR is a pointer returned from
   malloc (or NULL), pointing to *N characters of space.  It is realloc'd
   as necessary.  Return the number of characters read (not including the
   null terminator), or -1 on error or EOF.  */

int getstr (char ** lineptr, size_t *n, FILE * stream, char terminator, int offset)
{
  int nchars_avail;             /* Allocated but unused chars in *LINEPTR.  */
  char *read_pos;               /* Where we're reading into *LINEPTR. */
  int ret;

  if (!lineptr || !n || !stream)
    return -1;

  if (!*lineptr)
    {
      *n = 64;
      *lineptr = (char *) malloc (*n);
      if (!*lineptr)
        return -1;
    }

  nchars_avail = *n - offset;
  read_pos = *lineptr + offset;

  for (;;)
    {
      register int c = getc (stream);

      /* We always want at least one char left in the buffer, since we
         always (unless we get an error while reading the first char)
         NUL-terminate the line buffer.  */

      assert(*n - nchars_avail == read_pos - *lineptr);
      if (nchars_avail < 1)
        {
          if (*n > 64)
            *n *= 2;
          else
            *n += 64;

          nchars_avail = *n + *lineptr - read_pos;
          *lineptr = (char *) realloc (*lineptr, *n);
          if (!*lineptr)
            return -1;
          read_pos = *n - nchars_avail + *lineptr;
          assert(*n - nchars_avail == read_pos - *lineptr);
        }

      if (c == EOF || ferror (stream))
        {
          /* Return partial line, if any.  */
          if (read_pos == *lineptr)
            return -1;
          else
            break;
        }

      *read_pos++ = c;
      nchars_avail--;

      if (c == terminator)
        /* Return the line.  */
        break;
    }

  /* Done - NUL terminate and return the number of chars read.  */
  *read_pos = '\0';

  ret = read_pos - (*lineptr + offset);
  return ret;
}

#ifdef _WIN32
ssize_t getline(char **lineptr, size_t *n, FILE *stream)
{
  return getstr(lineptr, n, stream, '\n', 0);
}
#endif

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

bool file_newer(const char *f1, const char *f2) {
	struct stat st1, st2;
	st1.st_mtime = st2.st_mtime = 0;
	stat(f1, &st1);
	stat(f2, &st2);
	return (st1.st_mtime > st2.st_mtime);
}


std::vector<std::string> split(char* str,const char* delim)
{
    char* token = strtok(str,delim);
    std::vector<std::string> result;
    while(token != NULL)
    {
        result.push_back(token);
        token = strtok(NULL,delim);
    }
    return result;
}

int gzclose(FILE *f, bool isgz) {
    return isgz ? pclose(f) : fclose(f);
}

FILE *gzopen(const char *f, const char *m, bool*isgz) {
    // maybe use zlib some day?
        FILE *h;
        if (!strcmp(fext(f),".gz")) {
                char *tmp=(char *)malloc(strlen(f)+100);
                if (strchr(m,'w')) {
                        strcpy(tmp, "gzip  > '");
                        strcat(tmp, f);
                        strcat(tmp, "'");
                } else {
                        strcpy(tmp, "gunzip -c '");
                        strcat(tmp, f);
                        strcat(tmp, "'");
                }
        h = popen(tmp, m);
        *isgz=1;
        free(tmp);
        } else {
                h = fopen(f, m);
                *isgz=0;
        }
        if (!h) {
                fprintf(stderr, "Error opening file '%s': %s\n",f, strerror(errno));
                exit(1);
        }
        return h;
}

const char *fext(const char *f) {
        const char *x=strrchr(f,'.');
        return x ? x : "";
}

