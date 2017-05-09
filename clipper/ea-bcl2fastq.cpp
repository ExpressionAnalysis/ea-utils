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

/* 

See "void usage" below for usage.

*/

#include <stdarg.h>  // for va_start, etc
#include <memory>    // for std::unique_ptr
#include <getopt.h>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#include <time.h>
#include <vector>

using namespace std;            // bad practice

// #include "fastq-lib.h"

#define CHUNK 32768
#define MAX_ERR_FILES 10
#include "zlib.h"

#define VERSION "1.01.816"
#define warn(...) { fprintf(stderr, __VA_ARGS__); }
#define die(...) { warn(__VA_ARGS__); exit(1); }

std::string arg2cmdstr(int argc, char** argv);
std::string string_format(const std::string fmt_str, ...);
void usage(FILE *f, const char *msg=NULL);
FILE *openordie(const char *path, const char *mode);
FILE *popenordie(const char *path, const char *mode);
char* itoa(int value, char* result, int base, char **endp);

// per file/output file
typedef struct {
    bool useit;
    int cyc_offset;
    int cyc_len;
    int rnum;
    FILE *fout;
} mask;

// per-cycle info
typedef struct {
    bool useit;
    gzFile fin;
} cycle;

// tile record
typedef struct  __attribute__ ((packed)) {
    uint32_t tid;           // tile id
    uint32_t ccnt;          // cluster count
} tile_record;

int main (int argc, char **argv) {
    static struct option long_options[] = {
       {"debug", 0, 0, 0},
       {0, 0, 0, 0}
    };


    string run;                     // run path
    int lane=0;                     // lane
    string out;                     // output prefix
    vector<mask> masks;             // Y50N20Y50 format

    unsigned int cluster_start=0;                 // offset into cluster list, ZERO BASED 
    unsigned int cluster_count=0;                 // number of reads to process
    unsigned int output_cluster_count=0;                 // number of reads to process
    int tile=0;                                   // tile number
    int debug=0;                    // debug flag
    bool usegz=false;
    const char *fcid="X";

    int option_index = 0;
    int c;
    while (	(c = getopt_long(argc, argv, "zhr:l:t:o:m:s:n:f:",long_options,&option_index)) != -1) {
		switch (c) {
			case '\0':
                { 
                    const char *oname=long_options[option_index].name;
                    if(!strcmp(oname,        "debug")) {
                        debug=1;
                    }
                    break;
                }
			case 'h': usage(stdout); exit(0); break;
			case 'r': run = optarg; break;
			case 'l': lane = atoi(optarg); break;
			case 'o': out = optarg; break;
			case 't': tile = atoi(optarg); break;
			case 'f': fcid = optarg; break;
			case 'z': usegz = 1; break;
			case 'm': 
                    {
                        int typ;
                        char *p = optarg;
                        bool err=0;
                        int cur_offset=0;
                        while (*p) {
                            mask m; 
                            if (*p=='Y') {
                                m.useit=1;
                            } else if (*p=='N') {
                                m.useit=0;
                            } else {
                                err=1;
                            }
                            ++p;
                            int len;
                            if (isdigit(*p)) {
                                m.fout=NULL;
                                m.cyc_offset=cur_offset;
                                m.cyc_len=strtol(p, &p, 10);
                                cur_offset+=m.cyc_len;
                                masks.push_back(m);
                            } else {
                                err=1;
                            }
                            if (err) {
                                die("Mask should be something like: Y50Y30Y50");    
                            }
                       }
                    } 
                break;
			case 's': char *endp; cluster_start=strtoul(optarg, &endp, 10); break;
			case 'n': cluster_count=atoi(optarg); break;
			case '?': 
					  if (strchr("rltomsn", optopt))
						  fprintf(stderr, "Option -%c requires an argument.\n", optopt);
					  else if (isprint(optopt))
						  fprintf(stderr, "Unknown option `-%c'.\n", optopt);
					  else
						  fprintf(stderr,
								  "Unknown option character `\\x%x'.\n",
								  optopt);
					  usage(stderr);
					  return 1;
		}
	}

    if (!run.size() || !masks.size() || !lane || !out.length()) {
		die("Run, mask and lane are required.\n");
    }
  
    char lanestr[5];
    char lanestrnopad[5];
    sprintf(lanestr, "L%03d", lane);
    sprintf(lanestrnopad, "%d", lane);
   
    string locspath = run + "/Data/Intensities/" + lanestr + "/s_" + lanestrnopad + ".locs";
    string bcipath = run + "/Data/Intensities/BaseCalls/" + lanestr + "/s_" + lanestrnopad + ".bci";
    string filterpath = run + "/Data/Intensities/BaseCalls/" + lanestr + "/s_" + lanestrnopad + ".filter";

    // no die on flocs ... because it's ok to output no locs
    FILE *flocs = fopen(locspath.c_str(), "r");

    // die if no log, because presumably this is a write error in the output location
    FILE *flog = openordie((string(out) + ".log").c_str(), "w");

    // if no filter, then we'll just show all of them, so no die needed
    FILE *ffilter = fopen(filterpath.c_str(), "r");

	fprintf(flog, "Command Line: %s\n", arg2cmdstr(argc, argv).c_str());

    fprintf(flog,"Tile path: %s\n", bcipath.c_str());
    fprintf(flog,"Filter path: %s\n", filterpath.c_str());
    fflush(flog);

    FILE *ftnums = fopen(bcipath.c_str(), "r");

    struct {  
        uint32_t zero;
        uint32_t version;
        uint32_t numclusters;
    } filter_info;

    struct __attribute__ ((packed)) {
        // MAGIC NUMBER HEADER
        uint32_t field1;
        float field2;
        // real info
        uint32_t numclusters;
    } locs_info;

    // READ LOCS header
    if(!flocs || !fread(&locs_info,sizeof(locs_info),1,flocs)) {
        warn("Locs file is broken, no locations will be output\n");
        if (flocs) fclose(flocs);
        flocs=NULL;
    } else if (fseek(flocs,cluster_start*sizeof(float)*2,SEEK_CUR) < 0) {
        warn("Locs is there, but is no good\n");
        fclose(flocs);
        flocs=NULL;
    }

//    printf("TELL LOCS: %ld\n", ftell(flocs));

    // READ FILTER header
    bool ok=true;
    ok = ok && (fread(&filter_info.zero,4,1,ffilter)==1);
    ok = ok && (fread(&filter_info.version,4,1,ffilter)==1);
    ok = ok && (fread(&filter_info.numclusters,4,1,ffilter)==1);
    if(!ok) {
        die("Filter file is broken\n");
    }

    // ERROR/check numclusters
    if (flocs && filter_info.numclusters != locs_info.numclusters ) {
        die("Filter and locs numclusters don't match: %u\n", locs_info.numclusters);
    }

    if (!flocs) {
        fprintf(flog,"Locations invalid at: %u\n", cluster_start);
    }

    // READ TILE INFO
    vector <tile_record> tinfo;
    if (ftnums) {
        tile_record tr;
        while(fread(&tr, sizeof(tr), 1, ftnums)==1) {
//            printf("tile: %d, count: %d\n", tr.tid, tr.ccnt);
            tinfo.push_back(tr);
        }
    } else { 
        warn("Proceeding without tile info: %s\n", strerror(errno));
        fprintf(flog,"Tile numbers invalid at: %u\n", cluster_start);
    }

    // general purpose iterators
    int i,j;

    if (tile) {
        unsigned int cur=0;
        bool ok=false;
        for(i=0;i<tinfo.size();++i) {
            if (tinfo[i].tid==tile) {
                cluster_start=cur;
                cluster_count=tinfo[i].ccnt;
                ok=1;
                break;
            }
            cur+=tinfo[i].ccnt;
        }
        if (!ok) {
            die("Tile %d not found\n", tile);
        }
    }

//    warn("%d tiles read\n", (int) tinfo.size());

    // unsigned 8 bit integer
    if (fseek(ffilter,cluster_start,SEEK_CUR) < 0) {
        die("Can't seek in filter file\n");
    }
//    printf("TELL FILTER: %ld\n", ftell(ffilter));

    vector<cycle>cycles;
    int output_fnum=0;
    string outtmp;
    for(i=0;i<masks.size();++i) {
        FILE *fo=NULL;                            // null by default
        if(masks[i].useit) {
            ++output_fnum;                        // output file number is sequential
            masks[i].rnum=output_fnum;            // save file number as "read number"
            if (usegz) {
                outtmp = string_format("gzip -2  -c > %s.%d.fq.gz",out.c_str(),output_fnum); 
                fo=popenordie(outtmp.c_str(),"w");
            } else {
                outtmp = string_format("%s.%d.fq",out.c_str(),output_fnum); 
                fo=openordie(outtmp.c_str(),"w");
            }
        }
        masks[i].fout=fo;                         // pointer to output fiule for this mask level
        for (j = 0; j < masks[i].cyc_len; ++j) {
            cycle c;                            
            c.useit = masks[i].useit;
            cycles.push_back(c);
        }
    }

    string bclbase = run + "/Data/Intensities/BaseCalls/" + lanestr + "/";
    string bclpath;
    vector<gzFile> fbclv;      // vector of open files

    int err_files = 0;
    for (i=0;i<cycles.size();++i) {
        bclpath=bclbase + string_format("%04d.bcl.bgzf",i+1);
        gzFile fil=Z_NULL;
        if (cycles[i].useit) {
            fil=gzopen(bclpath.c_str(), "r");                                   // open the basecall
            if (fil == Z_NULL) {
                warn("Cycle %d file open failed, %s: %s\n", i, bclpath.c_str(), strerror(errno));
                ++err_files;
                if (err_files >= MAX_ERR_FILES) {
                    die("Too many errors, quitting\n");
                }
                cycles[i].useit=false;
            }  else {
                uint32_t numc;
                gzread(fil,&numc,4);                                                // read the header
                if (numc!=filter_info.numclusters) {
                    warn("Cycle %d num clusters mismatch/corrupt", i);
                    gzclose(fil);
                    fil=Z_NULL;
                    cycles[i].useit=false;
                } else {
    //                warn("Seek bcl %d\n",i);
                    if (gzseek(fil,cluster_start,SEEK_CUR)<0) {                            // seek to cluster_start (8 bits per record)
                        gzclose(fil);
                        fil=Z_NULL;
                        cycles[i].useit=false;
                    }
                }
            }

         
            if (!cycles[i].useit) {
                fprintf(flog,"Cycle %d invalid at %d\n", i+1, cluster_start);
            }
        }
        cycles[i].fin=fil;
    }

    struct {
        unsigned int base : 2;
        unsigned int qual : 6;
    } rec;

    char seqs[cycles.size()];
    char quals[cycles.size()];

    //ID Template:
    //@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>
    //@NS500184:5:H0K79AGXX:1:11103:20690:3982 1:N:0:ATTCAGAA+GCCTCTAT

    char read_id[1000];
    sprintf(read_id,"@NS:1:%s:",fcid);

    // pid is a pointer to the end of the id, after the lane: has been added
    char *pid_after_lane=read_id+strlen(read_id);
    itoa(lane, pid_after_lane, 10, &pid_after_lane);
    *pid_after_lane++ = ':';

    // map from aa to char
    char aa_map[4] = {'A','C','G','T'};

    // map from qual to score
    char qc_map[64];
    for(i=0;i<64;++i) 
        qc_map[i]=33+i;

    struct __attribute__ ((packed)) {
        float x;
        float y;
    } locrec;
 
    char pf;        /// purity filter (PF in illumina-speak)
   
    int tidx=(tinfo.size()>0)?0:-1; 
    int tileid=(tinfo.size()>0)?tinfo[tidx].tid:0;
    int trnum=0;

    // for each cluster requested... (we should be all seeked to the correct offsets at this point)

//    printf("TINFO: %d, %d\n", tileid, tinfo[tidx].ccnt);
    fprintf(flog,"Cluster count: %u\n", filter_info.numclusters);

    if (cluster_count == 0) {
        cluster_count = filter_info.numclusters;
    }

    fprintf(flog,"Cluster start: %u\n", cluster_start);
    fprintf(flog,"Cluster subset: %u\n", cluster_count);

    for(j=0;j<cluster_count;++j) {
        if (tidx > 0 && trnum > tinfo[tidx].ccnt) {
            ++tidx;
            trnum=0;
            if (tidx > tinfo.size()) {
                // tile numbers are invalid at this point... !
                fprintf(flog,"Tile numbers invalid at: %u\n", cluster_start+j);
                tidx=-1;
                tileid=0;
            }
            tileid=tinfo[tidx].tid;
//            printf("TINFO: %d, %d\n", tileid, tinfo[tidx].ccnt);
        }
        ++trnum;

        // read filter flag from filter file
        if(ffilter && (fread(&pf,1,1,ffilter)==1)) {
            pf = pf ? 'N' : 'Y';
        } else {
            pf = 'U';
        }

        // read x/y location from locs file
        int x, y;
        if (flocs && fread(&locrec,sizeof(locrec), 1, flocs)==1) {
            x=int(locrec.x * 10 + 1000 + 0.5);
            y=int(locrec.y * 10 + 1000 + 0.5);
        } else {
            x=0;
            y=0;
            if (flocs) {
                fprintf(flog,"Locations invalid at: %u\n", cluster_start+j);
                flocs = NULL;
            }
        }


        // read cycles
        for (i=0;i<cycles.size();++i) {
            if (cycles[i].useit) {
                // 1 byte read
                if(gzread(cycles[i].fin,&rec,1)==1) {
//                    warn("Record n:%d, cy:%d: b:%d, q:%d\n", j, i, rec.base, rec.qual);
                    if (rec.base==0 && rec.qual == 0) {
                        seqs[i]='N';
                        quals[i]='#';
                    } else {
                        char aa = aa_map[rec.base];
                        char qc = qc_map[rec.qual];
                        seqs[i]=aa;
                        quals[i]=qc;
                    }
                } else {
                    // zero out the cycle from now on... 
                    fprintf(flog,"Cycle %d invalid at %d\n", i+1, cluster_start);
                    cycles[i].useit = 0;
                    seqs[i]='N';
                    quals[i]='#';
                }
            } else {
                seqs[i]='N';
                quals[i]='#';
            }
        }
        
        // output read(s)
        if (pf != 'Y') {
            // convert tileid, x y to read header
            output_cluster_count++;
            char *tmpid = pid_after_lane;
            itoa(tileid, tmpid, 10, &tmpid);
            *tmpid++=':';
            itoa(x, tmpid, 10, &tmpid);
            *tmpid++=':';
            itoa(y, tmpid, 10, &tmpid);
            *tmpid++=' ';

            // id after the space
            char *pid_after_space = tmpid;
            for (i=0;i<masks.size();++i) {
                if (masks[i].useit) {
                    // output file number, pf flag and control flag
                    tmpid=pid_after_space;
                    itoa(masks[i].rnum, tmpid, 10, &tmpid);
                    *tmpid++=':';
                    *tmpid++=pf;
                    *tmpid++=':';
                    *tmpid++='0';
                    *tmpid='\0';

                    // output the id, sequence, and quals for the current file output
                    fputs(read_id,masks[i].fout);
                    fputc('\n',masks[i].fout);
                    fwrite(seqs+masks[i].cyc_offset,1,masks[i].cyc_len, masks[i].fout);
                    fputc('\n',masks[i].fout);
                    fputc('+',masks[i].fout);
                    fputc('\n',masks[i].fout);
                    fwrite(quals+masks[i].cyc_offset,1,masks[i].cyc_len, masks[i].fout),
                    fputc('\n',masks[i].fout);
                } 
            }
        }
    }


    for(i=0;i<masks.size();++i) {
        if(masks[i].fout && usegz) {
            if (!pclose(masks[i].fout)) {
                fprintf(flog, "Error : gzip file may be corrupt\n");
                die("Error : gzip file may be corrupt\n");
            }
        }
    }

    // all ok?
    fprintf(flog,"Cluster output: %u\n", output_cluster_count);
    exit(0);
}

void usage(FILE *f, const char *msg) {
	if(msg)
		fprintf(f, "%s\n", msg);

	fprintf(f, 
"Usage: ea-bcl2fastq -r PATH -l NUM -m MASK -o PREFIX [options] \n"
"Version: %s\n"
"\n"
"Converts Illumina bcl files to fastq\n"
"\n"
"Required:\n"
"    -r PATH     Path to run folder\n"
"    -m MASK     Y50Y6Y50\n"
"    -o PREFIX   Output file prefix\n"
"\n"
"Optional:\n"
"    -s START    Cluster offset (ZERO BASED OFFSET)\n"
"    -n COUNT    Cluster count\n"
"\n"
    ,VERSION);
}

char *arg2cmd(int argc, char** argv) {
    char *buf=NULL;
    int n = 0;
    int k, i;
    for (i=1; i <argc;++i) {
        int k=strlen(argv[i]);
        buf=( char *)realloc(buf,n+k+4);
        char *p=buf+n;
        char endq=0;
        // this is a poor mans quoting, which is good enough for anything that's not rediculous
        if (strchr(argv[i], ' ')) {
            if (!strchr(argv[i], '\'')) {
                *p++='\'';
                endq='\'';
            } else {
                *p++='\"';
                endq='\"';
            }
        }
        memcpy(p, argv[i], k);
        p+=k;
        if (i < (argc-1)) *p++=' ';
        if (endq) *p++=endq;
        *p='\0';
        n = p-buf;
    }
    return buf;
}

std::string arg2cmdstr(int argc, char **argv) {
    char *tmp=arg2cmd(argc, argv);
    std::string ret=tmp;
    free(tmp);
    return ret;
}

FILE *openordie(const char *path, const char *mode) {
    FILE *f=fopen(path, mode);
    if (!f) {
        warn("Can't open-%s %s: %s\n", mode, path, strerror(errno));
        exit(1);
    }
    return f;
}

FILE *popenordie(const char *path, const char *mode) {
    FILE *f=popen(path, mode);
    if (!f) {
        warn("Can't popen-%s %s: %s\n", mode, path, strerror(errno));
        exit(1);
    }
    return f;
}


std::string string_format(const std::string fmt, ...) {
    int size = 100;
    std::string str;
    va_list ap;
    while (1) {
        str.resize(size);
        va_start(ap, fmt);
        int n = vsnprintf((char *)str.c_str(), size, fmt.c_str(), ap);
        va_end(ap);
        if (n > -1 && n < size) {
            str.resize(n);
            return str;
        }
        if (n > -1)
            size = n + 1;
        else
            size *= 2;
    }
    return str;
}

 // fast int to string code: from http://www.jb.man.ac.uk/~slowe/cpp/itoa.html

char* itoa(int value, char* result, int base, char **endp) {
        // check that the base if valid
        if (base < 2 || base > 36) { *result = '\0'; return result; }
    
        char* ptr = result, *ptr1 = result, tmp_char;
        int tmp_value;
    
        do {
            tmp_value = value;
            value /= base;
            *ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz" [35 + (tmp_value - value * base)];
        } while ( value );
    
        // Apply negative sign
        if (tmp_value < 0) *ptr++ = '-';
        *endp=ptr;
        *ptr-- = '\0';
        while(ptr1 < ptr) {
            tmp_char = *ptr;
            *ptr--= *ptr1;
            *ptr1++ = tmp_char;
        }
        return result;
}

/* vim: set noai ts=4 sw=4: */
