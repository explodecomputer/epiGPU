#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#define MINAF 0.01
#define MINCALL 0.9
#define MISS 3
#define PI 3.14159

// space or tab delimited
// -9 or 0 = missing
// 1 = male, 2 = female
// PED: family ID | Individual ID | Paternal ID | Maternal ID | Sex | Phenotype |||| Geno

typedef struct ped {
	char family[30];
	char id[30];
	char paternal[30];
	char maternal[30];
	char sex;
	float phen;
} ped;

typedef struct map {
	char chr[30];
	char snpname[30];
	double gd;
	long pd;
	char allelename[2];
} map;

typedef struct chromosome {
	int chrid;
	char chrname[30];
	int chrsize;
	int chrstart;
} chromosome;

void dims(int *rowCount, int *colCount, char *filename);
char **get_space(int rows, int cols);
void release_space(char **a);
void readped(int nid, int nsnp, char *pedfile, ped *dat, char *geno, map *genmap);
void readmap(int nsnp, map* genmap, int *nchr, chromosome *chrstat, char *mapfile);
void cleangeno(int *nid, int *nsnp, int nchr, char *geno, ped *dat, chromosome *chrstat, map *genmap);
void guessmissing(int nid, int nsnp, char *geno);
void contphen(int nid, ped *dat);
void rnorm(int n, ped *dat);
void createbinary(int nid, int nsnp, int nchr, map *genmap, ped *dat, char *geno, chromosome *chrstat, char *filename);
void readbinary(int *nid, int *nsnp, int *nchr, map **genmap, ped **dat, char **geno, chromosome **chrstat, char *filename);
void simulate(int *nid, int *nsnp, int *nchr, map **genmap, ped **dat, char **geno, chromosome **chrstat);
void simqtl(int nid, int nsnp, char *geno, ped *dat, map *genmap);
void readgenabel(int *nid1, int *nsnp1, int *nchr1, map **genmap1, ped **dat1, char **geno1, chromosome **chrstat1, char *pedfile, char *mapfile);
int datamode(int argc, char **argv, int *nid, int *nsnp, int *npack, int *remain, int *nchr, ped **dat, map **genmap, chromosome **chrstat, char **binfile, char **pedfile, char **mapfile, int **genop, int **gpack, char **geno);
void arguments(char **argv);

void readpackedbinary(int *nid, int *nsnp, int *npack, int *remain, int *nchr, map **genmap1, ped **dat1, int **gpack1, chromosome **chrstat1, char *filename);
void createpackedbinary(int nid, int nsnp, int npack, int remain, int nchr, map *genmap, ped *dat, int *gpack, chromosome *chrstat, char *filename);
void unpack(int nid, int nsnp, int npack, int remain, int *gpack, char **geno1);
void pack(int nid, int nsnp, int *npack, int *remain, int **gpack, char *geno);
void bit_print(int a);
void extractfrombinary(char *binfile, char *pedfile, char *mapfile);
char isbinary(int nid, ped *dat);
void walkersalias(int n, float *p, int nans, int *ans);

