#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include <sys/stat.h>
#include <limits.h>
#include <vector>
#include <string>
#include <sparsehash/sparse_hash_map> // or sparse_hash_set, dense_hash_map, ...

#define max(a,b) (a>b?a:b)
#define min(a,b) (a<b?a:b)
#define meminit(l) (memset(&l,0,sizeof(l)))
#define die(...) { fprintf(stderr, __VA_ARGS__); exit(1); }

using namespace std;

typedef struct line {
        char *s; int n; size_t a;
} line;

class counts {
public:
	counts() {am=un=0;};
	int am;
	int un;
};

int read_line(FILE *in, struct line &l);
vector<char *> split(char *s, const char *d);

int main(int argc, char **argv) {
	char *in = argv[1];
	FILE * fin=fopen(in,"r");
	if (!fin) {
		die("%s:%s", in, strerror(errno));
	}
	struct line l;
	vector<char *> v;
        google::sparse_hash_map<string, counts> cmap;
	string lastid, lastref;
	int lastloc;
	int mated = 0;
	char *p;
	while(read_line(fin, l)>0) {
		if (l.s[0] == '@') continue;

		v = split(l.s,"\t");
		if (p=strchr(v[0],' ')) {
			if (isdigit(p[1])) {
				*p='\0';
			}	
		}
		string id=v[0];
		string ref=v[2];

		if (ref[0] == '*') continue;

		int loc=atoi(v[3]);
		int mateloc=atoi(v[7]);
		int isize=atoi(v[8]);
		if (mateloc > 0) mated = 1;
		counts *pcnt;
		if (loc == 0) continue;
		string firstref;
		bool newid=0, doitonce=0;
		if (lastref.size()>0) {
			if (mated) {
				if (id == lastid) {
					// proper pair second read?
					if (loc==mateloc && ref == lastref) {
						//  new probe id
						if (newid) {
							// previous proper-pair read was a 'newid'
							cmap[ref].un++;
							firstref=ref;
							doitonce=1;
						} else {
							// current tr is ambig
							cmap[ref].am++;
							if (doitonce) {
								// previous tr was not unambig
								cmap[firstref].un--;
								doitonce=0;
							}
						}
					}
					newid=0;
				} else {
					newid=1;
				}
			} else {
				if (id == lastid) {
					// current tr is ambig
					cmap[ref].am++;
					if (doitonce) {
						// and prev was not unamb
						cmap[lastref].un--;
						doitonce=0;
					}
				} else {
					cmap[ref].un++;
					doitonce= 1;	
				}
			}
		}

                lastid = id;
                lastloc = loc;
                lastref = ref;
	}
	google::sparse_hash_map<string, counts>::iterator it;
	for (it=cmap.begin();it!=cmap.end();++it) {
		printf("%s\t%d\t%d\n", it->first.c_str(), it->second.un, it->second.am);
	}
}

int read_line(FILE *in, struct line &l) {
        return (l.n = getline(&l.s, &l.a, in));
}

vector<char *> split(char *s, const char *d) {
	char *sp;
	std::vector<char *> v;
	char *t=strtok_r(s,d,&sp);
	while(t) {
		v.push_back(t);
		t=strtok_r(NULL,d,&sp);
	}
	return v;
}

