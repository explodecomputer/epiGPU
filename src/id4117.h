#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include <sys/time.h>
#include <math.h>
#include <ctype.h>

#define MAX 15000
#define MAX2 200
#define THRESHOLD 7

int silent;

typedef struct table {
	int snp1;
	int snp2;
	int df1;
	int df2;
	float F;
	float F2;
} table;


//int timediff(double *result, struct timeval *x, struct timeval *y);
void timediff(double *result, time_t x, time_t y);
void readtime(char *timetaken, double result);

#define KERN kernelfull

