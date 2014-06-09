///A tutorial about local alignments.
#include <iostream>
#include <seqan/align.h>
//#include <seqan/graph_msa.h>
//#include <seqan/consensus/consensus_calling.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <sys/stat.h>

using namespace seqan;

void revcomp(struct fq *dest, struct fq* src);

struct fq {
        char *id;   int nid;   size_t naid;
        char *seq;  int nseq;  size_t naseq;
        char *com;  int ncom;  size_t nacom;
        char *qual; int nqual; size_t naqual;
};
int read_fq(FILE *in, int rno, struct fq *fq);          // 0=done, 1=ok, -1=err+continue

const char *adapter = "TCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGA";

#define stdev(cnt, sum, ssq) sqrt((((double)cnt)*ssq-pow((double)sum,2)) / ((double)cnt*((double)cnt-1)))
int quantile(std::vector<int> vec, double p);
#define strendcmp(hay, nee) strcmp(hay+strlen(hay)-strlen(nee), nee)

int main(int argc, char **argv)
{
	if (argc < 2) {
		fprintf(stderr, "usage: pbi-clean FASTQ > CLEANED\n");
		exit(1);
	}

	int al = strlen(adapter);
	const char *in = argv[1];
	FILE *fin;

	if (!strendcmp(in,".gz")) {
		std::string gunz = "gunzip -c '";
		gunz += in;
		gunz += "'";
		fin=popen(gunz.c_str(), "r");
	} else {
		fin = !strcmp(in,"-")?stdin:fopen(in,"r");
	}

	if (!fin) {
                fprintf(stderr, "%s : %s\n",argv[1],strerror(errno));
		exit(1);
	}

        struct fq fq, rc;
        memset(&fq, 0, sizeof(fq));

	bool read_ok; int nrec=0;
	Score<int> score(3,-3,-2, -2);
	std::vector<int> cnt_vec;
	std::vector<int> len_vec;
	int cnt_sum=0, cnt_ssq=0, len_sum=0, len_ssq=0;
        while (read_ok=read_fq(fin, nrec, &fq)) {
		++nrec;
		Align< String<char> > ali;
		appendValue(rows(ali), adapter);
		appendValue(rows(ali), fq.seq);
		LocalAlignmentFinder<> finder(ali);
		if (fq.nseq > al*1.5) {
			try {
				int s;
				int cnt = 0;
				while ((s=localAlignment(ali, finder, score, (int)(al*1.8), WatermanEggert()))>0) {
					int b=clippedBeginPosition(row(ali, 1));
					int e=clippedEndPosition(row(ali, 1));
//					printf("%d<-->%d\n", b,e);
					fq.seq[b]='\1';
					fq.seq[e]='\2';
					++cnt;
				}
				if (cnt > 0) {
					cnt = 0;
					char *b = fq.seq;
					char *p = fq.seq;
					char *e = fq.seq+fq.nseq;
					char *t;
	//				Align< String<DnaQ> > sub;
					while (p <= e) {
						if (*p == '\1' || *p =='\0') {
							char *m = (char *) ((unsigned long)p/2+(unsigned long)b/2);
							for (t = b; t < p; ++t) {
								if (*t=='\1' || *t=='\2') {
									if (t < m) {
										b=t+1;
									} else {
										p=t-1;
										break;
									}
								}
							}
							*p = '\0';
							fq.qual[(p-fq.seq)]='\0';
							if ((p-b) > al && *b) {
							//	String<DnaQ> dq = b;
							//	assignQualities(dq,fq.qual+(b-fq.seq));
							//	if (cnt % 2 == 1)
							//		dq=DnaStringReverseComplement(dq);
							//	appendValue(rows(sub),dq);

								printf("%s:%d-%d\n", fq.id, (int)(b-fq.seq), (int)(p-fq.seq-1));
								printf("%s\n", b);
								printf("+\n");
								printf("%s\n", fq.qual+(b-fq.seq));
								len_vec.push_back(p-b);
								len_ssq+=(p-b)*(p-b);
								len_sum+=p-b;
								++cnt;
							}
							++p;
							while (p < e && *p != '\2') 
								++p;
							p=b=p+1;
						} else {
							++p;
						}
					}
					if (cnt > 1) {
						/*
						typedef String<DnaQ>         TSequence;  // sequence type
						typedef Align<TSequence,ArrayGaps>  TAlign;     // align type
						typedef Row<TAlign>::Type       TRow;
						typedef Iterator<TRow>::Type        TIterator;
						for (unsigned i = 1; i < length(rows(sub)); ++i) {
							TAlign pw;
							appendValue(rows(pw),row(sub,i-1));
							appendValue(rows(pw),row(sub,i));
							int sc = localAlignment(pw, score);
						printf("Score: %d\n", sc);
						::std::cout << pw << ::std::endl;
							
							TIterator it = iter(row(pw,0),0);
							int alEnd = _max(length(row(pw, 0)), length(row(pw, 1)));
							TIterator itEnd = iter(row(pw,0),alEnd);
						}*/
	//					globalMsaAlignment(sub, score);

						//computeProfiles(prof, pinfo, sub);
					}
				} else {
					cnt = 1;
					printf("%s:%d-%d\n", fq.id, 0, fq.nseq-1);
					printf("%s\n", fq.seq);
					printf("+\n");
					printf("%s\n", fq.qual);
					len_vec.push_back(fq.nseq);
					len_ssq+=fq.nseq*fq.nseq;
					len_sum+=fq.nseq;
				}
				cnt_vec.push_back(cnt);
				cnt_ssq+=cnt*cnt;
				cnt_sum+=cnt;
//				if (cnt_vec.size() > 50) {
//					goto JUMP;
//				}
			} catch (...) {
				fprintf(stderr,"Alignment error, line %d\n", nrec*4);
			}
		}
	}
//JUMP:
	fprintf(stderr,"reads\t%d\n",nrec);
	if (len_vec.size() > 0) {
		std::sort(cnt_vec.begin(), cnt_vec.end());
		std::sort(len_vec.begin(), len_vec.end());
		double sd;
		float mn = cnt_sum/(float)cnt_vec.size();
		int q1=quantile(len_vec,.25);
		int q2=quantile(len_vec,.50);
		int q3=quantile(len_vec,.75);
		fprintf(stderr,"subreads\t%d\n",cnt_sum);
		fprintf(stderr,"len max\t%d\n",len_vec[len_vec.size()-1]);
		if (cnt_vec.size() > 1) {
			fprintf(stderr,"len mean\t%2.2f\n",(double)len_sum/len_vec.size());
			sd = stdev(len_vec.size(), len_sum, len_ssq);
			fprintf(stderr,"len stdev\t%2.2f\n",sd);
			fprintf(stderr,"len q1\t%d\n",q1);
			fprintf(stderr,"len q2\t%d\n",q2);
			fprintf(stderr,"len q3\t%d\n",q3);
		}
		fprintf(stderr,"cyc max\t%d\n",cnt_vec[cnt_vec.size()-1]);
		if (cnt_vec.size() > 1) {
			int q1=quantile(cnt_vec,.25);
			int q2=quantile(cnt_vec,.50);
			int q3=quantile(cnt_vec,.75);
			fprintf(stderr,"cyc mean\t%2.2f\n",(double)cnt_sum/cnt_vec.size());
			sd = stdev(cnt_vec.size(), cnt_sum, cnt_ssq);
			fprintf(stderr,"cyc stdev\t%2.2f\n",sd);
			fprintf(stderr,"cyc q1\t%d\n",q1);
			fprintf(stderr,"cyc q2\t%d\n",q2);
			fprintf(stderr,"cyc q3\t%d\n",q3);
		}
	} else {
		fprintf(stderr,"subreads\t%d\n",0);
	}

	return 0;
}

#define MAXWARN 10
int read_fq(FILE *in, int rno, struct fq *fq) {
	static int fq_warncount=0;
        fq->nid = getline(&fq->id, &fq->naid, in);
        fq->nseq = getline(&fq->seq, &fq->naseq, in);
        fq->ncom = getline(&fq->com, &fq->nacom, in);
        fq->nqual = getline(&fq->qual, &fq->naqual, in);
        if (fq->nqual <= 0)
                return 0;
        if (fq->id[0] != '@' || fq->com[0] != '+' || fq->nseq != fq->nqual) {
                if (fq_warncount < MAXWARN) {
                        fprintf(stderr, "Malformed fastq record at line %d\n", rno*2+1);
                        ++fq_warncount;
                }
                return -1;
        }

        // chomp
        fq->id[--fq->nid] = '\0';
        fq->seq[--fq->nseq] = '\0';
        fq->qual[--fq->nqual] = '\0';
        return 1;
}

#define comp(c) ((c)=='A'?'T':(c)=='a'?'t':(c)=='C'?'G':(c)=='c'?'g':(c)=='G'?'C':(c)=='g'?'c':(c)=='T'?'A':(c)=='t'?'a':(c))
void revcomp(struct fq *d, struct fq *s) {
        if (!d->seq) {
                d->seq=(char *) malloc(d->naseq=s->nseq);
                d->qual=(char *) malloc(d->naqual=s->nqual);
        } else if (d->naseq < s->nseq) {
                d->seq=(char *) realloc(s->seq, d->naseq=s->nseq);
                d->qual=(char *) realloc(s->qual, d->naqual=s->nqual);
        }
        int i;
        for (i=0;i<s->nseq/2;++i) {
                char b=s->seq[i];
                char q=s->qual[i];
                //printf("%d: %c, %c\n", i, comp(s->seq.s[s->seq.n-i-1]), s->qual.s[s->qual.n-i-1]);
                d->seq[i]=comp(s->seq[s->nseq-i-1]);
                d->qual[i]=s->qual[s->nqual-i-1];
                //printf("%d: %c, %c\n", s->seq.n-i-1, comp(b), q);
                d->seq[s->nseq-i-1]=comp(b);
                d->qual[s->nseq-i-1]=q;
        }
        if (s->nseq % 2) {
                //printf("%d: %c, %c\n", 1+s->seq.n/2, comp(s->seq.s[s->seq.n/2]));
                d->seq[s->nseq/2] = comp(s->seq[s->nseq/2]);
                d->qual[s->nseq/2] = s->qual[s->nseq/2];
        }
        d->nseq=s->nseq;
        d->nqual=s->nqual;
}


int quantile(std::vector<int> vec, double p) {
        int l = vec.size();
        double t = ((double)l-1)*p;
	int it = (int) t;
        int v=vec[it];
        if (t > (double)it) {
                return v + (t-it) * (vec[it+1] - v);
        } else {
                return v;
        }
}

