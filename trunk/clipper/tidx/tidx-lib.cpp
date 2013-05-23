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

double xtime();

bool annot_comp (const annot &a, const annot &b) { return (a.beg < b.beg); }

template <typename L, typename R> void append(L& lhs, R const& rhs) { lhs.insert(lhs.end(), rhs.begin(), rhs.end()); }
template <typename L, typename R> void prepend(L& lhs, R const& rhs) { lhs.insert(lhs.begin(), rhs.begin(), rhs.end()); }

struct string_annot_serializer {
  bool operator()(FILE* fp, const std::pair<const string&, const vector<annot> >& value) const {

    {
    assert(value.first.length() <= UCHAR_MAX);
    const unsigned char size = value.first.length();
    if (fwrite(&size, sizeof(size), 1, fp) != 1)
      return false;
    if (fwrite(value.first.data(), size, 1, fp) != 1)
      return false;
    }

    {
    const vector<annot>&van=value.second;
    const unsigned long size = van.size();
    if (fwrite(&size, sizeof(size), 1, fp) != 1)
      return false;
    int i;
    for (i=0;i<size;++i) {
        if (fwrite(&van[i].beg, sizeof(van[i].beg), 1, fp) != 1)
          return false;
        if (fwrite(&van[i].end, sizeof(van[i].end), 1, fp) != 1)
          return false;
        assert(van[i].pos.size() <= USHRT_MAX);
        const unsigned short size = van[i].pos.size();
        if (fwrite(&size, sizeof(size), 1, fp) != 1)
          return false;
        int j;
        for (j=0;j<size;++j) {
            if (fwrite(&van[i].pos[j], sizeof(van[i].pos[j]), 1, fp) != 1)
              return false;
        }
    }
    }

    return true;
  }

  bool operator()(FILE* fp, std::pair<const string, vector<annot> >* value) const {

    {

    string buf;
    unsigned char size;    // all strings are <= 255 chars long
    if (fread(&size, sizeof(size), 1, fp) != 1)
      return false;

    if(size>buf.size()) buf.resize(size*2);

    if (fread((void *)buf.data(), size, 1, fp) != 1) {
      return false;
    }
    // necessarry to "new" the value which must be const, except during "unsearialization"
    // api shouldn't foist this on the user ... should be behind the scenes
    string * ncs = const_cast<string *>(&value->first);
    new(ncs) string(buf.data(), (size_t)size);
    
    }

    {
    
    vector<annot> &van=value->second;
    unsigned long size;
    if (fread(&size, sizeof(size), 1, fp) != 1)
      return false;
    int i;
    van.resize(size);
    for (i=0;i<size;++i) {
        if (fread(&van[i].beg, sizeof(van[i].beg), 1, fp) != 1)
          return false;
        if (fread(&van[i].end, sizeof(van[i].beg), 1, fp) != 1)
          return false;
        unsigned short size;
        if (fread(&size, sizeof(size), 1, fp) != 1)
          return false;
        int j;
        van[i].pos.resize(size);
        for(j=0;j<size;++j) {
            if (fread(&van[i].pos[j], sizeof(van[i].pos[j]), 1, fp) != 1)
              return false;
        }
    }

    }

    return true;
  }
};


void chomp_line(struct line &l) {
    if (l.s[l.n-1] == '\n') l.s[--l.n]='\0';       // chomp
    if (l.s[l.n-1] == '\r') l.s[--l.n]='\0';       // chomp
}

vector <long int> empty_vector;
const vector<long int> &tidx::lookup(const char *chr, int pos) {
    dense_hash_map<string,vector<annot> >::iterator it=map.find(chr);
    if (it == map.end()) return empty_vector;
    vector<annot> &va = it->second;
    if (debug) fprintf(stderr,"lookup: %s:%d -> %d\n", chr, pos, (int) va.size());
    int b=0, t=va.size(), c=0;
    while (t>b) {
        c=(t+b)/2;
//        printf("here1: c:%d, t:%d, b:%d, pos:%d, beg:%d, end:%d, res:%d\n", c, t, b, pos, va[c].beg, va[c].end, va[c].pos[0]);
        if (pos == va[c].beg)
            break;
        else if (pos < va[c].beg)
            t=c-1;
        else if (pos > va[c].beg) {
            if (pos <= va[c].end) {
                return va[c].pos;
            }
            b=c+1;
        }
    }
    
    if (t == b)
        c = t;
//    printf("here2: c:%d, t:%d, b:%d, pos:%d, beg:%d, end:%d, res:%d\n", c, t, b, pos, va[c].beg, va[c].end, va[c].pos[0]);
    if (pos >= va[c].beg && pos <= va[c].end) {
        return va[c].pos;
    }
    return empty_vector;
}

vector<long int> tidx::lookup_r(const char *chr, int beg, int end) {
    dense_hash_map<string,vector<annot> >::iterator it=map.find(chr);
    if (it == map.end()) return empty_vector;
    vector<annot> &va = it->second;
    if (debug) fprintf(stderr,"lookup_r: %s:%d.%d -> %d\n", chr, beg, end, (int) va.size());
    int b=0, t=va.size(), c=0;
    while (t>b) {
        c=(t+b)/2;
        if (beg == va[c].beg)
            break;
        else if (beg < va[c].beg)
            t=c-1;
        else if (beg > va[c].beg) {
            if (beg <= va[c].end) 
                break;
            b=c+1;
        }
    }
    if (t == b)
        c = t;
    vector<long int> res;
    while (c<va.size() && end >= va[c].beg && beg <= va[c].end) {
        append(res,va[c].pos);
        ++c;
    }
    return res;
}

string tidx::lookup(const char *chr, int pos, const char *msep) { 
//    printf("here2\n");
    const vector<long int> &v = lookup(chr, pos);
    string res;
    if (!fh) {
        fh=fopen(path.c_str(),"rb");
        if (!fh) 
            fail("%s:%s\n",path.c_str(),strerror(errno));
    }
    string line;
    int i;
	struct line l; meminit(l);
    for (i=0;i<v.size();++i) {
        fseek(fh,v[i],0);
    	read_line(fh, l);
        chomp_line(l);
        res += msep;
        res += string(l.s, l.n);
    }
    free_line(&l);
    return res;
}

string tidx::lookup_r(const char *chr, int beg, int end, const char *msep) {
//    printf("here2\n");
    const vector<long int> &v = lookup_r(chr, beg, end);
    string res;
    if (!fh) {
        fh=fopen(path.c_str(),"rb");
        if (!fh)
            fail("%s:%s\n",path.c_str(),strerror(errno));
    }
    string line;
    int i;
    struct line l; meminit(l);
    for (i=0;i<v.size();++i) {
        fseek(fh,v[i],0);
        read_line(fh, l);
        chomp_line(l);
        res += msep;
        res += string(l.s, l.n);
    }
    free_line(&l);
    return res;
}

string api_ret = "";
const char *tidx::lookup_c(const char *chr, int pos, const char *msep) {
    api_ret = lookup(chr, pos, msep);
    return api_ret.c_str();
}

const char *tidx::lookup_cr(const char *chr, int beg, int end, const char *msep) {
    api_ret = lookup_r(chr, beg, end, msep);
    return api_ret.c_str();
}

bool tidx::read(const char *in) {
    string uin = string_format("gunzip -c %s.tidx", in);

    if (debug) fprintf(stderr, "read %s\n", in);
    FILE *fun=popen(uin.c_str(),"r");
    if (!fun) {
        return false;
    }
    map.unserialize(string_annot_serializer(), fun);
    path=in;
    return true;
}

void tidx::init() {
    debug=false;
    fh=NULL;
    map.set_empty_key("-");
}

void tidx::dump(FILE *fh) {
    fprintf(fh,"#file\t%s\n",path.c_str());
    dense_hash_map<string,vector<annot> >::iterator it = map.begin();;
    while (it != map.end()) {
        vector<annot> &van = it->second;
        int i;
        for (i=0;i<van.size();++i) {
            fprintf(fh, "%s\t%d\t%d\t%ld\t%ld\n", it->first.c_str(), van[i].beg, van[i].end, van[i].pos.size(), van[i].pos[0]);
        }
        ++it;
    }    
}

// fun part
void tidx::build(const char *in, const char *sep, int nchr, int nbeg, int nend, int skip_i, char skip_c, bool sub_e) {
	FILE *fin=fopen(in,"r");

    if (!fin)
        fail("%s:%s\n",in,strerror(errno));

    if (nend == -1)
        nend = nbeg;

    string out = string_format("gzip -c > %s.tidx", in);
    FILE *fout=popen(out.c_str(),"w");
    if (!fout)
        fail("%s:%s\n", out.c_str(),strerror(errno));

    double xst = xtime();

	struct line l; meminit(l);
    int nlast = max(max(nbeg,nend),nchr);
    int nl = 0;

    string p_chr = "%";
    path=in;
    vector<annot> *pvan;
    long tpos = ftell(fin);
    // read in the annotation file
    if (debug) fprintf(stderr, "reading %s (%d, %d, %d)\n", in, nchr, nbeg, nend);
	while (read_line(fin, l)>0) {
        ++nl;
        if (skip_i > 0 || *l.s==skip_c) {
            --skip_i;
        } else {
            chomp_line(l);
            vector<char *> v = split(l.s, sep);
            if (nlast >= v.size()) {
                fail("error, file %s, line %d: missing info\n", in, nl);
            } 
            annot a;
            a.beg=atoi(v[nbeg]);
            a.end=atoi(v[nend]);
            if (sub_e) --a.end;
            if (a.beg > a.end) {
                fail("error, file %s, line %d: beg > end : %d > %d\n", in, nl, a.beg, a.end);
            }
            a.pos.push_back(tpos);
            if (strcmp(v[nchr], p_chr.c_str())) {       // speed up
                pvan = &map[v[nchr]];
                p_chr=v[nchr];
            }
            pvan->push_back(a);
        }
        tpos = ftell(fin);
	}
    dense_hash_map<string,vector<annot> >::iterator it;

    free_line(&l);
 
    // for each chromosome
    it = map.begin();
    while (it != map.end()) {
        vector<annot> &van = it->second;
        // sort the annotation file by beginning of region
        sort(van.begin(), van.end(), annot_comp);
        int i;
        if (debug) fprintf(stderr, "frag %s : %ld ->", it->first.c_str(), van.size());
        for (i=0;i<van.size()-1;++i) {
            if (van[i].beg >= van[i+1].beg && van[i].end == van[i+1].end) {
                // exact match
                if (debug) fprintf(stderr, " [dup %d]", van[i].beg);
                // merge annotations
                prepend(van[i+1].pos,van[i].pos);
                // skip next... (empty pos won't be serialized)
                assert(van[i].beg == van[i+1].beg);
                van[i].pos.clear();
            } else if (van[i].end >= van[i+1].beg) {
                if (debug) fprintf(stderr, " [ovr %d-%d:%ld ]", van[i].beg, van[i].end, van[i].pos[0]);
                // overlap next
                int new_st;
                int new_en;

                // forced to initialize here so we can use a reference (for efficiency)
                vector<long> new_ro = van[i].pos;
                
                if (van[i].end < van[i+1].end) {
                    // contained within next, so new frag starting after i stop
                    new_st = van[i].end + 1;
                    new_en = van[i+1].end;
                    new_ro = van[i+1].pos;                  // that only contains the other
                    van[i+1].end=van[i].end;                // shorten next to my end
                    append(van[i+1].pos,van[i].pos);        // and now the other contains me
                } else {
                    // passes next, so next contains all of me
                    new_st = van[i+1].end+1;                // new frag is after the end of next
                    new_en = van[i].end;
                    append(van[i+1].pos,van[i].pos);
                }

                van[i].end=van[i+1].beg-1;                  // shorten my end to less than the next's start
                
                if (debug) fprintf(stderr, " [i:%d:%d:%ld]", van[i].beg, van[i].end, van[i].pos[0]);
                if (debug) fprintf(stderr, " [i+1:%d:%d:%ld]", van[i+1].beg, van[i+1].end, van[i+1].pos[0]);

                if (new_en >= new_st) {                     // is this a real one?
                    if (debug) fprintf(stderr, " [n:%d:%d:%ld]", new_st, new_en, new_ro[0]);

                    int j = i+2;                            // figure out where it goes (shouldn't be far)
                    while (j < van.size() & new_st > van[j].beg) {
                        ++j;
                    }

                    annot a;
                    a.beg=new_st;
                    a.end=new_en;
                    a.pos=new_ro;
                                                            // (slow... use linked list, turn to array later for storage/bin search?)
                    van.insert(van.begin()+j, a);           // insert into the annot array
                }
            }
        }
        long j = 0;
        for (i=0;i<van.size();++i) {
            // overlap next
            if (van[i].pos.size() != 0 && van[i].beg <= van[i].end)
                if (i != j) 
                    van[j++]=van[i];
                else
                    ++j;
        }
        if (j != van.size()) {
            if (debug) fprintf(stderr, "(rm %ld) ", van.size()-j);
            van.resize(j);
        }
        if (debug) fprintf(stderr, " %ld\n", van.size());
        ++it;
    }

    double xen = xtime();
    double speed = xen-xst;

    if (debug) fprintf(stderr, "compiled in %g secs\n", speed);

    map.serialize(string_annot_serializer(), fout);
    pclose(fout);

    //
    xst = xtime();
    tidx tmap(in);
    xen = xtime();
    speed = xen-xst;
    if (debug) fprintf(stderr, "read in %g secs\n", speed);
    path=in;
}

double xtime() {
    struct timeval tm;
    gettimeofday(&tm, NULL);
    return (double) tm.tv_sec + ((double)tm.tv_usec)/1000000.0;
}

// for the api
void tidx_build(const char *file, const char *sep, int chr, int beg, int end, int skip_i, char skip_c, bool sub_e) {
    tidx n;
    n.build(file, sep, chr, beg, end, skip_i, skip_c, sub_e); 
}
