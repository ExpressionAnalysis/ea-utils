#ifdef SWIG
%module "Text::Tidx::SWIG"
%{
#include "tidx.h"
void tidx_build(const char *path, const char *sep, int nchr, int nbeg, int nend, int skip_i, char skip_c, bool sub_e);
%}
#endif

class annot {
public:
    int beg;
    int end;
    std::vector<long> pos;
};

class tidx {
public:
    bool debug;
    tidx() {init();};
    tidx(const char *path)  {init(); read(path);};

    std::string path;
    google::dense_hash_map<std::string,std::vector<annot> > map;

    bool read(const char *path);

    /* interface for single-position */
    const std::vector <long int> & lookup(const char *chr, int pos);
    const char * lookup_c(const char *chr, int pos, const char *msep);

    /* interface for range overlap */
    const char * lookup_cr(const char *chr, int beg, int end, const char *msep);
    const std::vector <long int> & lookup_r(const char *chr, int beg, int end);
};

void tidx_build(const char *path, const char *sep, int nchr, int nbeg, int nend, int skip_i, char skip_c, bool sub_e);
