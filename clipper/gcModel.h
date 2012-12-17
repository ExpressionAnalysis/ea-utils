/*
$Id$
*/
extern void gcInit(int maxReadLength);
extern void gcProcessSequence(int l,int c);
extern void gcPrintDistribution(FILE *fp);
void gcClose();
