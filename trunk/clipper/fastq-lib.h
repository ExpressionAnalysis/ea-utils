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

// 32-bit o/s support
#if defined(__i386__)
	#define _FILE_OFFSET_BITS 64
#endif

// standard libs
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <assert.h>
#include <math.h>
#include <sys/stat.h>
#include <search.h>
#include <limits.h>
#include <stdint.h>
#include <stddef.h>

#if defined(__APPLE__)
	#define getopt(a,b,c) getopt_long(a,b,c,NULL,NULL)
#endif

// misc useful macros
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
#define meminit(l) (memset(&l,0,sizeof(l)))
#define fail(s,...) ((fprintf(stderr,s,##__VA_ARGS__), exit(1)))
#define warn(s,...) ((fprintf(stderr,s,##__VA_ARGS__)))
#define stdev(cnt, sum, ssq) sqrt((((double)cnt)*ssq-pow((double)sum,2)) / ((double)cnt*((double)cnt-1)))

// maximum number of files that can be tracked by poorquals lib
#define MAX_FILENO_QUALS 6

// read line, read fq
typedef struct line {
        char *s; int n; size_t a;
} line;

struct fq {
        line id;
        line seq;
        line com;
        line qual;
};


void free_line(struct line *l);
void free_fq(struct fq *fq);

// not GNU?  probably no getline & strtok_r...
#if !defined( __GNUC__) || defined(WIN32) || defined(__APPLE__)
	ssize_t getline(char **lineptr, size_t *n, FILE *stream);
	char* strtok_r(char *str, const char *delim, char **nextp);
#endif
    
// get file extension
const char *fext(const char *f);

// read fq
int read_line(FILE *in, struct line &l);                // 0=done, 1=ok, -1=err+continue
int read_fq(FILE *in, int rno, struct fq *fq, const char *name=NULL);          // 0=done, 1=ok, -1=err+continue
int read_fq_sam(FILE *in, int rno, struct fq *fq, const char *name=NULL);          // 0=done, 1=ok, -1=err+continue
void free_fq(struct fq *fq);

// open a file, possibly gzipped, exit on failure
FILE *gzopen(const char *in, const char *mode, bool *isgz);
int gzclose(FILE *f, bool isgz);

// keep track of poor quals (n == "file number", maybe should have persistent stat struct instead?)
bool poorqual(int n, int l, const char *s, const char *q);

// returns number of differences between 2 strings, where n is the "max-length to check"
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

// reverse complement an fq entry into a blank (memset 0) one
void revcomp(struct fq *dest, struct fq* src);


