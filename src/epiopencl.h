#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else 
#include <CL/cl.h>
#endif
#include "id4117.h"
#include "readplink.h"

int UPERM, perm;
int UDEVICE;
int UCHUNKSIZE;
int UTEST;
int USAFE;
char plinkbin[100];
float UTHRESHOLD1, UTHRESHOLD2;
char UPHEN[1000];


char *load_program_source(const char *filename);
void ereport(int ret, int code);
int perform_scan(int *geno, chromosome *chrstat, ped *dat, map *genmap, int nid, int nsnp, int npack, int remain, int nchr, char *filename);
void permute(int nid, ped *dat);
void resume(char *filename, int *chr1, int *chr2, int *tothits, int nid, int nsnp);
void parsefilename(char *plinkbin);
void analysismodecl(int argc, char **argv, char **binfile, char **outfile);
void readphenotype(char *filename, int nid, ped *dat);


