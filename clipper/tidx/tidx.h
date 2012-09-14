#include <string>
#include <vector>
#include <sparsehash/dense_hash_map>

class annot {
public:
    int beg;
    int end;
    std::vector<long> pos;
};

class tidx {
    FILE *fh;
    void init();
public:
    bool debug;
    tidx() {init();};
    tidx(const char *path)  {init(); read(path);};

    std::string path;
    google::dense_hash_map<std::string,std::vector<annot> > map;

    void dump(FILE *stream);
    bool read(const char *path);
    void build(const char *path, const char *sep, int nchr, int nbeg, int nend, int skip_i, char skip_c, bool sub_e);

    const std::vector <long int> & lookup(const char *chr, int pos);
    std::string lookup(const char *chr, int pos, const char *msep);

// range lookup
    std::vector <long int> lookup_r(const char *chr, int beg, int end);
    std::string lookup_r(const char *chr, int beg, int end, const char *msep);

// const char * return value
    const char * lookup_c(const char *chr, int pos, const char *msep);
    const char * lookup_cr(const char *chr, int beg, int end, const char *msep);
};

void chomp_line(struct line &l);

// build, with no return value, for API use
void tidx_build(const char *path, const char *sep, int nchr, int nbeg, int nend, int skip_i, char skip_c, bool sub_e);

