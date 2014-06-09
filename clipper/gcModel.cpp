/*
$Id$
*/
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
#include <iostream>

#include "fastq-lib.h"
#include "gcModel.h"

// #define UNIT_TEST
 
void gcInit(int maxReadLength);
void gcProcessSequence(int l,int c);
void gcPrintDistribution(FILE *fp);
void gcClose();

using namespace std;

#define roundgt0(x) (long)(x<0.5?0:x+0.5)

typedef struct GCModelValue {
  int percentage;
  double increment;
} GC_MODEL_VALUE; 

typedef struct GCModelValues {
  GC_MODEL_VALUE * values;
  int valuesLength;
} GC_MODEL_VALUES;

typedef GC_MODEL_VALUES *GC_MODELS;


static int claimingCounts[101]; 
static double gcDistribution[101]; 
static GC_MODELS *cachedModels;
static int gMaxReadLength = -1;

GC_MODEL_VALUES *calcModels(int readLength) {

  memset(claimingCounts,0,sizeof(claimingCounts));

  GC_MODEL_VALUES *models = (GC_MODEL_VALUES *) malloc((readLength+1) * sizeof(GC_MODEL_VALUES ));
  memset(models,0,(readLength+1) * sizeof(GC_MODEL_VALUES ));

  for (int pos=0;pos<=readLength;pos++) {
    double lowCount = pos-0.5;
    double highCount = pos+0.5;
    
    if (lowCount < 0.0) lowCount = 0.0;
    if (highCount < 0.0) highCount = 0.0;
    if (highCount > readLength) highCount = readLength;
    if (lowCount > readLength) lowCount = readLength;
    
    int lowPercentage = (int)roundgt0(((lowCount*100) / readLength));
    int highPercentage = (int)roundgt0(((highCount*100) / readLength));
    
    for (int p=lowPercentage;p<=highPercentage;p++) {
      claimingCounts[p]++;
    }
  }

  // We now do a second pass to make up the model using the weightings
  // we calculated previously.
  
  for (int pos=0;pos<=readLength;pos++) {
    double lowCount = pos-0.5;
    double highCount = pos+0.5;
    
    if (lowCount < 0) lowCount = 0;
    if (highCount < 0) highCount = 0;
    if (highCount > readLength) highCount = readLength;
    if (lowCount > readLength) lowCount = readLength;
    
    int lowPercentage = (int)roundgt0((lowCount*100) / readLength);
    int highPercentage = (int)roundgt0((highCount*100) / readLength);
    
    models[pos].values = (GC_MODEL_VALUE *) malloc(((highPercentage-lowPercentage)+1) * sizeof(GC_MODEL_VALUE) );
    memset(models[pos].values,0,
	   ((highPercentage-lowPercentage)+1) * sizeof(GC_MODEL_VALUE) );
    models[pos].valuesLength = (highPercentage-lowPercentage)+1;

    for (int p=lowPercentage;p<=highPercentage;p++) {
      models[pos].values[p-lowPercentage].percentage = p;
      models[pos].values[p-lowPercentage].increment = 1.0/claimingCounts[p];
    }
  }
  
  return (models);
}


void gcProcessSequence(int l,int c) {

  if(l > gMaxReadLength) { printf("Error: read length (%d) exceeds specified maximum length(%d)\n", l, gMaxReadLength); }
  if(c > l) { printf("Error: GC-count (%d) exceeds actual read length(%d)\n", c, l) ;}

  GC_MODEL_VALUE *values = cachedModels[l][c].values;

  for(int i=0; i < cachedModels[l][c].valuesLength; i++) {
    gcDistribution[values[i].percentage] += values[i].increment;
  }

}

void printModels(int rl) {
  GC_MODEL_VALUES *m = cachedModels[rl];

  printf("## Model values for read length=%d\n",rl);

  for(int i = 0; i <= rl; i++) {
    printf("%d: ",i);
    for(int j = 0; j < m[i].valuesLength; j++) {
      printf("%d,%.2f ",m[i].values[j].percentage, m[i].values[j].increment);
    }
    printf("\n");
  }
}

void gcPrintDistribution(FILE *fp) {
  if(fp == NULL) {
    fp = stdout;
  }
  fprintf(fp, "pct_GC\tCount\n");
  for(int i=0; i<=100;i++) {
    fprintf(fp, "%d\t%.2f\n",i,gcDistribution[i]);
  }
}

void gcClose() {
  if(gMaxReadLength < 0)return; // never initialized

  for(int rl = 0; rl < gMaxReadLength; rl++) {
    GC_MODEL_VALUES * m =  cachedModels[rl];
    for(int i = 0; i <= rl; i++) {
      free(m[i].values);
    }
    free(m);
  }
  
  free(cachedModels);
}

void gcInit(int maxReadLength) {
  gMaxReadLength = maxReadLength;

  memset(gcDistribution,0,sizeof(gcDistribution));
  // Build all models for a given max readlength:
  cachedModels = (GC_MODELS*)malloc((maxReadLength+1) * sizeof(GC_MODELS));
  // original code fills this in,caching, as necessary 
  // here, we just build all models at outset:
  int pos;
  for( pos = 0; pos <= maxReadLength; pos++) {
    cachedModels[pos] = calcModels(pos);
  }
}

#ifdef UNIT_TEST
main() {

  //  int maxReadLength = 35;
  int maxReadLength = 5;

  gcInit(maxReadLength);

  // ***
  // simulate processing A sequence:
  //  int seqLength = 3; // this sequence's length
  //  int gcCount = 2; // total G's & C's -- count 'em
  /*
  for(int i = 0; i < 10000000; i++) {
    gcProcessSequence(35,15);
  }
  for(int i = 0; i < 5000000; i++) {
    gcProcessSequence(35,10);
  }
  */

    printModels(4);
  //  exit(0);

  //  for(int pos=0; pos <= maxReadLength; pos++) {
  //    gcProcessSequence(maxReadLength,pos);
    //  }

    //  gcProcessSequence(3,0);
    //  gcProcessSequence(4,2);
    //  gcProcessSequence(5,4);

    //  gcPrintDistribution(NULL);

 
  gcClose();

}
#endif
