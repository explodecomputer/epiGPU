#include "epiopencl.h"



// FULLY OPTIMISED
// GEFORCE 8800GT: 1092842 tests / second
// with vectorisation: 1242173 tests / second


/*
Function description:
Called by clEnqueueNDRangeKernel, performs an 8df F-test for all SNP pairs 
in a square (or triangle) of the chromosome x chromosome search grid

Main optimisation strategies implemented:
- Warping
2D thread arrays are executed in multiples of 16. The biggest I/O constraint 
is fetching genotype data, this should occur simultaneously across threads
- Bit packing
Genotypes are stored as 0,1,2 (or 3=missing), so only 2 bits required. 
16 genotypes are packed into ints, and organised with the variables 'npack' 
and 'remain'
- Loop unrolling
Factorising the genotypes and calculating class means requires a 9-way 
switch for each genotype
- Phenotype variance and mean calculated simultaneously in O(n)
Sum of squares and mean of phenotype (SSY and mY) are precalculated at 
the beginning of the programme. When missing values are encountered the 
mY and SSY are updated by removing the corresponding phenotype values by 
converting the mean and SSW calculations to 'online algorithms' (modified 
from Welford 1962)
- Local memory
The small cache (16kb) constrains the local work size, but it is fast to 
access. It is used to store class means and class counts for all threads 
in the local work group.
- Avoiding atomic operations
When an F value exceeds the threshold it is stored in the results table 
(global). This is a pre-allocated array that is partitioned into chunks,
one chunk for each local thread. Each local thread only writes to its own
chunk to avoid race conditions. The alternative is to use atomic functions,
but this forces all global reads to follow the slow path on ATI cards, and has
some performance penalty on NVIDIA.
*/

char *kernelfull =
"typedef struct table {\n" \
"	int snp1;\n" \
"	int snp2;\n" \
"	int df1;\n" \
"	int df2;\n" \
"	float F;\n" \
"	float F2;\n" \
"} table;\n" \
"__kernel void square(\n" \
"	__global int *gchr1,  // Pointers to bitpacked genotype data\n" \
"	__global int *gchr2,  // Pointers to bitpacked genotype data\n" \
"	__global const float *phen,\n" \
"	__global table *results,\n" \
"	const float mY,\n" \
"	const float SSY,\n" \
"	const float threshold1,\n" \
"	const float threshold2,\n" \
"	const int nid,\n" \
"	const int nsnp,\n" \
"	const int npack,\n" \
"	const int remain,\n" \
"	__local float *mY_fac,\n" \
"	__local unsigned short int *factor_count,\n" \
"	__global int *nhits,\n" \
"	const int chr1,\n" \
"	const int chr2,\n" \
"	const int chr1g,\n" \
"	const int chr2g,\n" \
"	const int chr1m,\n" \
"	const int chr2m\n" \
"){\n" \
"	int\n" \
"		lsx = get_local_size(0),\n" \
"		lsy = get_local_size(1),\n" \
"		gsx = get_global_size(0),\n" \
"		gsy = get_global_size(1),\n" \
"		gx = get_global_id(0)+chr1,\n" \
"		gy = get_global_id(1)+chr2,\n" \
"		lx = get_local_id(0),\n" \
"		ly = get_local_id(1);\n" \
"	int i,j,k,l,n,npack1,gen1,gen2,locnid=nid;\n" \
"	unsigned mask;\n" \
"	int tmpi, g1, g2;\n" \
"	__global int *g1p, *g2p;\n" \
"	__global const float *phenp;\n" \
"	float tmpf;\n" \
"	float4 vtmp;\n" \
"	float locmy = mY, locssy = SSY, delta;\n" \
// chr1g and chr2g are ranks of the chromosomes' first SNP in the genome
// chr1m and chr2m are the boundaries of the chromosomes, required for diagonal regions
// Calculations only proceed if the SNPs are in the lower triangle of the search space
"	i = gx + chr1g; j = gy + chr2g;\n" \
"	if((i > j) && (gx < chr1m) && (gy < chr2m))\n" \
"	{\n" \
// lad is the local memory array start position for the thread's factor_count and mY_fac
// Each thread in a local work group has its own allocation of 9 short ints and 9 floats respectively
// 'lad' is the offset variable for the thread's access to the arrays
"	int lad = 9*(lx*lsx + ly);\n" \
"	for(i = 0; i < 9; i++)\n" \
"	{\n" \
"		tmpi = lad+i;\n" \
"		factor_count[tmpi] = 0;\n" \
"		mY_fac[tmpi] = 0;\n" \
"	}\n" \
"	k = 0;\n" \
// For each genotype fetch from global memory 16 packed genotypes are imported, packed into an int
// They must be unpacked sequentially to update the class counts / class means
"	npack1 = remain ? npack + 1 : npack;\n" \
"	g1p = gchr1 + npack1*gx;\n" \
"	g2p = gchr2 + npack1*gy;\n" \
"	for(i = 0; i < npack; i++)\n" \
"	{\n" \
// synchronise global memory fetches into warps
"		g1 = *(g1p++);\n" \
"		g2 = *(g2p++);\n" \
"		n = 30;\n" \
"		mask = 3; mask <<= n;\n" \
// begin unpacking data and accumulating factor info
// the phenotype call is vectorised, so the loop is unrolled into 4 iterations of float4 calls
"		for(j = 0; j < 4; j++)\n" \
"		{\n" \
"			phenp = phen + k; k+=4;\n" \
"			vtmp = vload4(0,phenp);\n" \
// reveal genotype data 2 bits at a time
"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.x;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
// Convert genotypes to 2x2 genotype classes (9 genotype classes from 2 SNPs)
// Count how many individuals for each class (stored locally in short int)
// Calculate genotype class means (stored locally in float)
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
// update SSY and mY in event of missing data
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \

"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.y;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \

"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.z;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \

"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.w;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \
"		}\n" \
"	}\n" \

// if the number of individuals isn't a factor of 16 then complete the factorising
"	if(remain)\n" \
"	{\n" \
"		g1 = *g1p;\n" \
"		g2 = *g2p;\n" \
"		n = (remain-1) * 2;\n" \
"		for(j = (remain-1); j >= 0; j--)\n" \
"		{\n" \
"			mask = 3;\n" \
"			mask <<= n;\n" \
"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			tmpf = phen[k++];\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
// Convert genotypes to 2x2 genotype classes (9 genotype classes from 2 SNPs)
// Count how many individuals for each class (stored locally in short int)
// Calculate genotype class means (stored locally in float)
"				tmpi = 3*gen1+gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \
"		}\n" \
"	}\n" \


"	int nfac = 0;\n" \
"	float SSB = 0;\n" \
// loop unrolling...
// calculate SSB class means
"	for(i = 0; i < 9; i++)\n" \
"	{\n" \
"		if(factor_count[lad+i] > 0)\n" \
"		{\n" \
"			nfac++;\n" \
"			mY_fac[lad+i] /= factor_count[lad+i];\n" \
"			SSB += factor_count[lad+i] * ((mY_fac[lad+i] - locmy) * (mY_fac[lad+i] - locmy));\n" \
"		}\n" \
"	}\n" \
// calculate F
"	int d1 = nfac - 1;\n" \
"	int d2 = locnid - nfac;\n" \
"	float Ftemp = (SSB/d1) / ((locssy-SSB)/d2);\n" \

// If F is above threshold then store it in global results array
// position in global results array is determined by how many have been found so far and local thread id
"	if(Ftemp > threshold1)\n" \
"	{\n" \
"		tmpi = nhits[lx*lsx + ly];\n" \
"		nhits[lx*lsx + ly] += 1;\n" \
"		tmpi += (200*(lx*lsx + ly));\n" \
"		results[tmpi].snp1 = gx+chr1g;\n" \
"		results[tmpi].snp2 = gy+chr2g;\n" \
"		results[tmpi].df1 = d1;\n" \
"		results[tmpi].df2 = d2;\n" \
"		results[tmpi].F = Ftemp;\n" \
"		results[tmpi].F2 = 0;\n" \
"	}\n" \
"}\n" \
"}";


/*
The logic of this function is identical to kernelfull
It also calculates the F test when parameterised only for interaction terms
This is done by subtracting the SSQ of the marginal genotypes from the SSQ of the full genotypes
as an extra step in the F test calculation. The full F test is also calculated and this is the
value that must surpass the threshold.

*/

char *kernelint =
"typedef struct table {\n" \
"	int snp1;\n" \
"	int snp2;\n" \
"	int df1;\n" \
"	int df2;\n" \
"	float F;\n" \
"	float F2;\n" \
"} table;\n" \
"__kernel void square(\n" \
"	__global int *gchr1,  // Pointers to bitpacked genotype data\n" \
"	__global int *gchr2,  // Pointers to bitpacked genotype data\n" \
"	__global const float *phen,\n" \
"	__global table *results,\n" \
"	const float mY,\n" \
"	const float SSY,\n" \
"	const float threshold1,\n" \
"	const float threshold2,\n" \
"	const int nid,\n" \
"	const int nsnp,\n" \
"	const int npack,\n" \
"	const int remain,\n" \
"	__local float *mY_fac,\n" \
"	__local unsigned short int *factor_count,\n" \
"	__global int *nhits,\n" \
"	const int chr1,\n" \
"	const int chr2,\n" \
"	const int chr1g,\n" \
"	const int chr2g,\n" \
"	const int chr1m,\n" \
"	const int chr2m\n" \
"){\n" \
"	int\n" \
"		lsx = get_local_size(0),\n" \
"		lsy = get_local_size(1),\n" \
"		gsx = get_global_size(0),\n" \
"		gsy = get_global_size(1),\n" \
"		gx = get_global_id(0)+chr1,\n" \
"		gy = get_global_id(1)+chr2,\n" \
"		lx = get_local_id(0),\n" \
"		ly = get_local_id(1);\n" \
"	int i,j,k,l,n,npack1,gen1,gen2,locnid=nid;\n" \
"	unsigned mask;\n" \
"	int tmpi, g1, g2;\n" \
"	__global int *g1p, *g2p;\n" \
"	__global const float *phenp;\n" \
"	float tmpf;\n" \
"	float4 vtmp;\n" \
"	float locmy = mY, locssy = SSY, delta;\n" \
// chr1g and chr2g are ranks of the chromosomes' first SNP in the genome
// chr1m and chr2m are the boundaries of the chromosomes, required for diagonal regions
// Calculations only proceed if the SNPs are in the lower triangle of the search space
"	i = gx + chr1g; j = gy + chr2g;\n" \
"	if((i > j) && (gx < chr1m) && (gy < chr2m))\n" \
"	{\n" \
// lad is the local memory array start position for the thread's factor_count and mY_fac
// Each thread in a local work group has its own allocation of 9 short ints and 9 floats respectively
// 'lad' is the offset variable for the thread's access to the arrays
"	int lad = 9*(lx*lsx + ly);\n" \
"	for(i = 0; i < 9; i++)\n" \
"	{\n" \
"		tmpi = lad+i;\n" \
"		factor_count[tmpi] = 0;\n" \
"		mY_fac[tmpi] = 0;\n" \
"	}\n" \
"	k = 0;\n" \
// For each genotype fetch from global memory 16 genotypes are imported, packed into an int
// They must be unpacked sequentially to update the class counts / class means
"	npack1 = remain ? npack + 1 : npack;\n" \
"	g1p = gchr1 + npack1*gx;\n" \
"	g2p = gchr2 + npack1*gy;\n" \
"	for(i = 0; i < npack; i++)\n" \
"	{\n" \
"		g1 = *(g1p++);\n" \
"		g2 = *(g2p++);\n" \
"		n = 30;\n" \
"		mask = 3; mask <<= n;\n" \
"		for(j = 0; j < 4; j++)\n" \
"		{\n" \
"			phenp = phen + k; k+=4;\n" \
"			vtmp = vload4(0,phenp);\n" \
"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.x;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \

"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.y;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \

"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.z;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \

"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.w;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \
"		}\n" \
"	}\n" \

"	if(remain)\n" \
"	{\n" \
"		g1 = *g1p;\n" \
"		g2 = *g2p;\n" \
"		n = (remain-1) * 2;\n" \
"		for(j = (remain-1); j >= 0; j--)\n" \
"		{\n" \
"			mask = 3;\n" \
"			mask <<= n;\n" \
"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			tmpf = phen[k++];\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
// Convert genotypes to 2x2 genotype classes (9 genotype classes from 2 SNPs)
// Count how many individuals for each class (stored locally in short int)
// Calculate genotype class means (stored locally in float)
"				tmpi = 3*gen1+gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \
"		}\n" \
"	}\n" \

"	float tmpf, mean_row[3], mean_col[3], SSI = 0;\n" \
"	int nfac = 0;\n" \
"	float SSB = 0;\n" \
"	for(i = 0; i < 9; i++)\n" \
"	{\n" \
"		if(factor_count[lad+i] > 0)\n" \
"		{\n" \
"			nfac++;\n" \
"		}\n" \
"	}\n" \
"	int d1 = nfac - 1;\n" \
"	for(i = 0; i < 3; i++)\n" \
"	{\n" \
// Calculate the independent genotype means for SNPs at each chromosome
"		mean_col[i] = (mY_fac[lad+i] + mY_fac[lad+i+3] + mY_fac[lad+i+6]) / (factor_count[lad+i] + factor_count[lad+i+3] + factor_count[lad+i+6]);\n" \
"		mean_row[i] = (mY_fac[lad+3*i] + mY_fac[lad+3*i+1] + mY_fac[lad+3*i+2]) / (factor_count[lad+3*i] + factor_count[lad+3*i+1] + factor_count[lad+3*i+2]);\n" \
"	}\n" \
"	for(i = 0; i < 9; i++)\n" \
"	{\n" \
// Calculate SSB
"		if(factor_count[lad+i]) {\n" \
"			mY_fac[lad+i] /= factor_count[lad+i];\n" \
"			SSB += factor_count[lad+i] * ((mY_fac[lad+i] - locmy) * (mY_fac[lad+i] - locmy));\n" \
"			tmpi = i % 3;\n" \
"			j = (int)floor((float)i/3); \n" \
// Calculate SSI as the residual SSB after subtracting the variance accounted for by marginal effects
"			tmpf = mY_fac[lad+i] - mean_row[j] - mean_col[tmpi] + locmy;\n" \
"			SSI += factor_count[lad+i] * tmpf * tmpf;\n" \
"		}\n" \
"	}\n" \
"	int d2 = locnid - nfac;\n" \
"	tmpf = locssy - SSB;\n" \
"	float MSW = tmpf / d2;\n" \
"	float Ftemp = (SSB/d1) / MSW;\n" \
"	float Ftemp2 = (SSI/4) / MSW;\n" \
"	if(Ftemp > threshold1 || Ftemp2 > threshold2)\n" \
"	{\n" \
"		if(nfac > 6)\n" \
"		{\n" \
"		//tmpi = atom_inc(&nhits[0]);\n" \
"		tmpi = nhits[lx*lsx + ly];\n" \
"		nhits[lx*lsx + ly] += 1;\n" \
"		tmpi += (200*(lx*lsx + ly));\n" \
"		results[tmpi].snp1 = gx;\n" \
"		results[tmpi].snp2 = gy;\n" \
"		results[tmpi].df1 = d1;\n" \
"		results[tmpi].df2 = d2;\n" \
"		results[tmpi].F = Ftemp;\n" \
"		results[tmpi].F2 = Ftemp2;\n" \
"		}\n" \
"	}\n" \
"}\n" \
"}\n";


char *kernelfullsafe =
"typedef struct table {\n" \
"	int snp1;\n" \
"	int snp2;\n" \
"	int df1;\n" \
"	int df2;\n" \
"	float F;\n" \
"	float F2;\n" \
"} table;\n" \
"__kernel void square(\n" \
"	__global int *gchr1,  // Pointers to bitpacked genotype data\n" \
"	__global int *gchr2,  // Pointers to bitpacked genotype data\n" \
"	__global const float *phen,\n" \
"	__global table *results,\n" \
"	const float mY,\n" \
"	const float SSY,\n" \
"	const float threshold1,\n" \
"	const float threshold2,\n" \
"	const int nid,\n" \
"	const int nsnp,\n" \
"	const int npack,\n" \
"	const int remain,\n" \
"	__local float *mY_fac,\n" \
"	__local unsigned int *factor_count,\n" \
"	__global int *nhits,\n" \
"	const int chr1,\n" \
"	const int chr2,\n" \
"	const int chr1g,\n" \
"	const int chr2g,\n" \
"	const int chr1m,\n" \
"	const int chr2m\n" \
"){\n" \
"	int\n" \
"		lsx = get_local_size(0),\n" \
"		lsy = get_local_size(1),\n" \
"		gsx = get_global_size(0),\n" \
"		gsy = get_global_size(1),\n" \
"		gx = get_global_id(0)+chr1,\n" \
"		gy = get_global_id(1)+chr2,\n" \
"		lx = get_local_id(0),\n" \
"		ly = get_local_id(1);\n" \
"	int i,j,k,l,n,npack1,gen1,gen2,locnid=nid;\n" \
"	unsigned mask;\n" \
"	int tmpi, g1, g2;\n" \
"	__global int *g1p, *g2p;\n" \
"	__global const float *phenp;\n" \
"	float tmpf;\n" \
"	float4 vtmp;\n" \
"	float locmy = mY, locssy = SSY, delta;\n" \
// chr1g and chr2g are ranks of the chromosomes' first SNP in the genome
// chr1m and chr2m are the boundaries of the chromosomes, required for diagonal regions
// Calculations only proceed if the SNPs are in the lower triangle of the search space
"	i = gx + chr1g; j = gy + chr2g;\n" \
"	if((i > j) && (gx < chr1m) && (gy < chr2m))\n" \
"	{\n" \
// lad is the local memory array start position for the thread's factor_count and mY_fac
// Each thread in a local work group has its own allocation of 9 short ints and 9 floats respectively
// 'lad' is the offset variable for the thread's access to the arrays
"	int lad = 9*(lx*lsx + ly);\n" \
"	for(i = 0; i < 9; i++)\n" \
"	{\n" \
"		tmpi = lad+i;\n" \
"		factor_count[tmpi] = 0;\n" \
"		mY_fac[tmpi] = 0;\n" \
"	}\n" \
"	k = 0;\n" \
// For each genotype fetch from global memory 16 packed genotypes are imported, packed into an int
// They must be unpacked sequentially to update the class counts / class means
"	npack1 = remain ? npack + 1 : npack;\n" \
"	g1p = gchr1 + npack1*gx;\n" \
"	g2p = gchr2 + npack1*gy;\n" \
"	for(i = 0; i < npack; i++)\n" \
"	{\n" \
// synchronise global memory fetches into warps
"		g1 = *(g1p++);\n" \
"		g2 = *(g2p++);\n" \
"		n = 30;\n" \
"		mask = 3; mask <<= n;\n" \
// begin unpacking data and accumulating factor info
// the phenotype call is vectorised, so the loop is unrolled into 4 iterations of float4 calls
"		for(j = 0; j < 4; j++)\n" \
"		{\n" \
"			phenp = phen + k; k+=4;\n" \
"			vtmp = vload4(0,phenp);\n" \
// reveal genotype data 2 bits at a time
"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.x;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
// Convert genotypes to 2x2 genotype classes (9 genotype classes from 2 SNPs)
// Count how many individuals for each class (stored locally in short int)
// Calculate genotype class means (stored locally in float)
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
// update SSY and mY in event of missing data
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \

"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.y;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \

"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.z;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \

"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.w;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \
"		}\n" \
"	}\n" \

// if the number of individuals isn't a factor of 16 then complete the factorising
"	if(remain)\n" \
"	{\n" \
"		g1 = *g1p;\n" \
"		g2 = *g2p;\n" \
"		n = (remain-1) * 2;\n" \
"		for(j = (remain-1); j >= 0; j--)\n" \
"		{\n" \
"			mask = 3;\n" \
"			mask <<= n;\n" \
"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			tmpf = phen[k++];\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
// Convert genotypes to 2x2 genotype classes (9 genotype classes from 2 SNPs)
// Count how many individuals for each class (stored locally in short int)
// Calculate genotype class means (stored locally in float)
"				tmpi = 3*gen1+gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \
"		}\n" \
"	}\n" \


"	int nfac = 0;\n" \
"	float SSB = 0;\n" \
// loop unrolling...
// calculate SSB class means
"	for(i = 0; i < 9; i++)\n" \
"	{\n" \
"		if(factor_count[lad+i] > 0)\n" \
"		{\n" \
"			nfac++;\n" \
"			mY_fac[lad+i] /= factor_count[lad+i];\n" \
"			SSB += factor_count[lad+i] * ((mY_fac[lad+i] - locmy) * (mY_fac[lad+i] - locmy));\n" \
"		}\n" \
"	}\n" \
// calculate F
"	int d1 = nfac - 1;\n" \
"	int d2 = locnid - nfac;\n" \
"	float Ftemp = (SSB/d1) / ((locssy-SSB)/d2);\n" \

// If F is above threshold then store it in global results array
// position in global results array is determined by how many have been found so far and local thread id
"	if(Ftemp > threshold1)\n" \
"	{\n" \
"		tmpi = nhits[lx*lsx + ly];\n" \
"		nhits[lx*lsx + ly] += 1;\n" \
"		tmpi += (200*(lx*lsx + ly));\n" \
"		results[tmpi].snp1 = gx+chr1g;\n" \
"		results[tmpi].snp2 = gy+chr2g;\n" \
"		results[tmpi].df1 = d1;\n" \
"		results[tmpi].df2 = d2;\n" \
"		results[tmpi].F = Ftemp;\n" \
"		results[tmpi].F2 = 0;\n" \
"	}\n" \
"}\n" \
"}";


/*
The logic of this function is identical to kernelfull
It also calculates the F test when parameterised only for interaction terms
This is done by subtracting the SSQ of the marginal genotypes from the SSQ of the full genotypes
as an extra step in the F test calculation. The full F test is also calculated and this is the
value that must surpass the threshold.

*/

char *kernelintsafe =
"typedef struct table {\n" \
"	int snp1;\n" \
"	int snp2;\n" \
"	int df1;\n" \
"	int df2;\n" \
"	float F;\n" \
"	float F2;\n" \
"} table;\n" \
"__kernel void square(\n" \
"	__global int *gchr1,  // Pointers to bitpacked genotype data\n" \
"	__global int *gchr2,  // Pointers to bitpacked genotype data\n" \
"	__global const float *phen,\n" \
"	__global table *results,\n" \
"	const float mY,\n" \
"	const float SSY,\n" \
"	const float threshold1,\n" \
"	const float threshold2,\n" \
"	const int nid,\n" \
"	const int nsnp,\n" \
"	const int npack,\n" \
"	const int remain,\n" \
"	__local float *mY_fac,\n" \
"	__local unsigned int *factor_count,\n" \
"	__global int *nhits,\n" \
"	const int chr1,\n" \
"	const int chr2,\n" \
"	const int chr1g,\n" \
"	const int chr2g,\n" \
"	const int chr1m,\n" \
"	const int chr2m\n" \
"){\n" \
"	int\n" \
"		lsx = get_local_size(0),\n" \
"		lsy = get_local_size(1),\n" \
"		gsx = get_global_size(0),\n" \
"		gsy = get_global_size(1),\n" \
"		gx = get_global_id(0)+chr1,\n" \
"		gy = get_global_id(1)+chr2,\n" \
"		lx = get_local_id(0),\n" \
"		ly = get_local_id(1);\n" \
"	int i,j,k,l,n,npack1,gen1,gen2,locnid=nid;\n" \
"	unsigned mask;\n" \
"	int tmpi, g1, g2;\n" \
"	__global int *g1p, *g2p;\n" \
"	__global const float *phenp;\n" \
"	float tmpf;\n" \
"	float4 vtmp;\n" \
"	float locmy = mY, locssy = SSY, delta;\n" \
// chr1g and chr2g are ranks of the chromosomes' first SNP in the genome
// chr1m and chr2m are the boundaries of the chromosomes, required for diagonal regions
// Calculations only proceed if the SNPs are in the lower triangle of the search space
"	i = gx + chr1g; j = gy + chr2g;\n" \
"	if((i > j) && (gx < chr1m) && (gy < chr2m))\n" \
"	{\n" \
// lad is the local memory array start position for the thread's factor_count and mY_fac
// Each thread in a local work group has its own allocation of 9 short ints and 9 floats respectively
// 'lad' is the offset variable for the thread's access to the arrays
"	int lad = 9*(lx*lsx + ly);\n" \
"	for(i = 0; i < 9; i++)\n" \
"	{\n" \
"		tmpi = lad+i;\n" \
"		factor_count[tmpi] = 0;\n" \
"		mY_fac[tmpi] = 0;\n" \
"	}\n" \
"	k = 0;\n" \
// For each genotype fetch from global memory 16 genotypes are imported, packed into an int
// They must be unpacked sequentially to update the class counts / class means
"	npack1 = remain ? npack + 1 : npack;\n" \
"	g1p = gchr1 + npack1*gx;\n" \
"	g2p = gchr2 + npack1*gy;\n" \
"	for(i = 0; i < npack; i++)\n" \
"	{\n" \
"		g1 = *(g1p++);\n" \
"		g2 = *(g2p++);\n" \
"		n = 30;\n" \
"		mask = 3; mask <<= n;\n" \
"		for(j = 0; j < 4; j++)\n" \
"		{\n" \
"			phenp = phen + k; k+=4;\n" \
"			vtmp = vload4(0,phenp);\n" \
"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.x;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \

"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.y;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \

"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.z;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \

"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = vtmp.w;\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \
"		}\n" \
"	}\n" \

"	if(remain)\n" \
"	{\n" \
"		g1 = *g1p;\n" \
"		g2 = *g2p;\n" \
"		n = (remain-1) * 2;\n" \
"		for(j = (remain-1); j >= 0; j--)\n" \
"		{\n" \
"			mask = 3;\n" \
"			mask <<= n;\n" \
"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			tmpf = phen[k++];\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
// Convert genotypes to 2x2 genotype classes (9 genotype classes from 2 SNPs)
// Count how many individuals for each class (stored locally in short int)
// Calculate genotype class means (stored locally in float)
"				tmpi = 3*gen1+gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \
"		}\n" \
"	}\n" \

"	float tmpf, mean_row[3], mean_col[3], SSI = 0;\n" \
"	int nfac = 0;\n" \
"	float SSB = 0;\n" \
"	for(i = 0; i < 9; i++)\n" \
"	{\n" \
"		if(factor_count[lad+i] > 0)\n" \
"		{\n" \
"			nfac++;\n" \
"		}\n" \
"	}\n" \
"	int d1 = nfac - 1;\n" \
"	for(i = 0; i < 3; i++)\n" \
"	{\n" \
// Calculate the independent genotype means for SNPs at each chromosome
"		mean_col[i] = (mY_fac[lad+i] + mY_fac[lad+i+3] + mY_fac[lad+i+6]) / (factor_count[lad+i] + factor_count[lad+i+3] + factor_count[lad+i+6]);\n" \
"		mean_row[i] = (mY_fac[lad+3*i] + mY_fac[lad+3*i+1] + mY_fac[lad+3*i+2]) / (factor_count[lad+3*i] + factor_count[lad+3*i+1] + factor_count[lad+3*i+2]);\n" \
"	}\n" \
"	for(i = 0; i < 9; i++)\n" \
"	{\n" \
// Calculate SSB
"		if(factor_count[lad+i]) {\n" \
"			mY_fac[lad+i] /= factor_count[lad+i];\n" \
"			SSB += factor_count[lad+i] * ((mY_fac[lad+i] - locmy) * (mY_fac[lad+i] - locmy));\n" \
"			tmpi = i % 3;\n" \
"			j = (int)floor((float)i/3); \n" \
// Calculate SSI as the residual SSB after subtracting the variance accounted for by marginal effects
"			tmpf = mY_fac[lad+i] - mean_row[j] - mean_col[tmpi] + locmy;\n" \
"			SSI += factor_count[lad+i] * tmpf * tmpf;\n" \
"		}\n" \
"	}\n" \
"	int d2 = locnid - nfac;\n" \
"	tmpf = locssy - SSB;\n" \
"	float MSW = tmpf / d2;\n" \
"	float Ftemp = (SSB/d1) / MSW;\n" \
"	float Ftemp2 = (SSI/4) / MSW;\n" \
"	if(Ftemp > threshold1 || Ftemp2 > threshold2)\n" \
"	{\n" \
"		//tmpi = atom_inc(&nhits[0]);\n" \
"		tmpi = nhits[lx*lsx + ly];\n" \
"		nhits[lx*lsx + ly] += 1;\n" \
"		tmpi += (200*(lx*lsx + ly));\n" \
"		results[tmpi].snp1 = gx;\n" \
"		results[tmpi].snp2 = gy;\n" \
"		results[tmpi].df1 = d1;\n" \
"		results[tmpi].df2 = d2;\n" \
"		results[tmpi].F = Ftemp;\n" \
"		results[tmpi].F2 = Ftemp2;\n" \
"	}\n" \
"}\n" \
"}\n";

char *kernelbase3 =
"typedef struct table {\n" \
"	int snp1;\n" \
"	int snp2;\n" \
"	int df1;\n" \
"	int df2;\n" \
"	float F;\n" \
"	float F2;\n" \
"} table;\n" \
"__kernel void square(\n" \
"	__global int *gchr1,  // Pointers to bitpacked genotype data\n" \
"	__global int *gchr2,  // Pointers to bitpacked genotype data\n" \
"	__global const float *phen,\n" \
"	__global table *results,\n" \
"	const float mY,\n" \
"	const float SSY,\n" \
"	const float threshold1,\n" \
"	const float threshold2,\n" \
"	const int nid,\n" \
"	const int nsnp,\n" \
"	const int npack,\n" \
"	const int remain,\n" \
"	__local float *mY_fac,\n" \
"	__local unsigned short int *factor_count,\n" \
"	__global int *nhits,\n" \
"	const int chr1,\n" \
"	const int chr2,\n" \
"	const int chr1g,\n" \
"	const int chr2g,\n" \
"	const int chr1m,\n" \
"	const int chr2m\n" \
"){\n" \
"	int\n" \
"		lsx = get_local_size(0),\n" \
"		lsy = get_local_size(1),\n" \
"		gsx = get_global_size(0),\n" \
"		gsy = get_global_size(1),\n" \
"		gx = get_global_id(0)+chr1,\n" \
"		gy = get_global_id(1)+chr2,\n" \
"		lx = get_local_id(0),\n" \
"		ly = get_local_id(1);\n" \
"	int i,j,k,l,n,npack1,gen1,gen2,locnid=nid;\n" \
"	unsigned mask;\n" \
"	int tmpi, g1, g2;\n" \
"	__global int *g1p, *g2p;\n" \
"	__global const float *phenp;\n" \
"	float tmpf;\n" \
"	//float4 vtmp;\n" \
"	float locmy = mY, locssy = SSY, delta;\n" \
// chr1g and chr2g are ranks of the chromosomes' first SNP in the genome
// chr1m and chr2m are the boundaries of the chromosomes, required for diagonal regions
// Calculations only proceed if the SNPs are in the lower triangle of the search space
"	i = gx + chr1g; j = gy + chr2g;\n" \
"	if((i > j) && (gx < chr1m) && (gy < chr2m))\n" \
"	{\n" \
// lad is the local memory array start position for the thread's factor_count and mY_fac
// Each thread in a local work group has its own allocation of 9 short ints and 9 floats respectively
// 'lad' is the offset variable for the thread's access to the arrays
"	int lad = 9*(lx*lsx + ly);\n" \
"	for(i = 0; i < 9; i++)\n" \
"	{\n" \
"		tmpi = lad+i;\n" \
"		factor_count[tmpi] = 0;\n" \
"		mY_fac[tmpi] = 0;\n" \
"	}\n" \
"	k = 0;\n" \
// For each genotype fetch from global memory 16 packed genotypes are imported, packed into an int
// They must be unpacked sequentially to update the class counts / class means
"	npack1 = remain ? npack + 1 : npack;\n" \
"	g1p = gchr1 + npack1*gx;\n" \
"	g2p = gchr2 + npack1*gy;\n" \
"	for(i = 0; i < npack; i++)\n" \
"	{\n" \
// synchronise global memory fetches into warps
"		g1 = *(g1p++);\n" \
"		g2 = *(g2p++);\n" \
"		n = 30;\n" \
"		mask = 3; mask <<= n;\n" \
// begin unpacking data and accumulating factor info
// the phenotype call is vectorised, so the loop is unrolled into 4 iterations of float4 calls
"		for(j = 0; j < 16; j++)\n" \
"		{\n" \
"			//phenp = phen + k; k+=4;\n" \
"			//vtmp = vload4(0,phenp);\n" \
// reveal genotype data 2 bits at a time
"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			mask >>= 2;\n" \
"			tmpf = phen[k++];\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
// Convert genotypes to 2x2 genotype classes (9 genotype classes from 2 SNPs)
// Count how many individuals for each class (stored locally in short int)
// Calculate genotype class means (stored locally in float)
"				tmpi = 3*gen1 + gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
// update SSY and mY in event of missing data
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \
"		}\n" \
"	}\n" \

// if the number of individuals isn't a factor of 16 then complete the factorising
"	if(remain)\n" \
"	{\n" \
"		g1 = *g1p;\n" \
"		g2 = *g2p;\n" \
"		n = (remain-1) * 2;\n" \
"		for(j = (remain-1); j >= 0; j--)\n" \
"		{\n" \
"			mask = 3;\n" \
"			mask <<= n;\n" \
"			gen1 = ((g1 & mask) >> n);\n" \
"			gen2 = ((g2 & mask) >> n);\n" \
"			n -= 2;\n" \
"			tmpf = phen[k++];\n" \
"			if(gen1 != 3 && gen2 != 3)\n" \
"			{\n" \
// Convert genotypes to 2x2 genotype classes (9 genotype classes from 2 SNPs)
// Count how many individuals for each class (stored locally in short int)
// Calculate genotype class means (stored locally in float)
"				tmpi = 3*gen1+gen2;\n" \
"				factor_count[lad+tmpi] += 1;\n" \
"				mY_fac[lad+tmpi] += tmpf;\n" \
"			} else {\n" \
"				delta = (locnid*locmy - tmpf)/(--locnid);\n" \
"				locssy -= (tmpf - locmy)*(tmpf - delta);\n" \
"				locmy = delta;\n" \
"			}\n" \
"		}\n" \
"	}\n" \


"	int nfac = 0;\n" \
"	float SSB = 0;\n" \
// loop unrolling...
// calculate SSB class means
"	for(i = 0; i < 9; i++)\n" \
"	{\n" \
"		if(factor_count[lad+i] > 0)\n" \
"		{\n" \
"			nfac++;\n" \
"			mY_fac[lad+i] /= factor_count[lad+i];\n" \
"			SSB += factor_count[lad+i] * ((mY_fac[lad+i] - locmy) * (mY_fac[lad+i] - locmy));\n" \
"		}\n" \
"	}\n" \
// calculate F
"	int d1 = nfac - 1;\n" \
"	int d2 = locnid - nfac;\n" \
"	float Ftemp = (SSB/d1) / ((locssy-SSB)/d2);\n" \

// If F is above threshold then store it in global results array
// position in global results array is determined by how many have been found so far and local thread id
"	if(Ftemp > threshold1)\n" \
"	{\n" \
"		tmpi = nhits[lx*lsx + ly];\n" \
"		nhits[lx*lsx + ly] += 1;\n" \
"		tmpi += (200*(lx*lsx + ly));\n" \
"		results[tmpi].snp1 = gx+chr1g;\n" \
"		results[tmpi].snp2 = gy+chr2g;\n" \
"		results[tmpi].df1 = d1;\n" \
"		results[tmpi].df2 = d2;\n" \
"		results[tmpi].F = Ftemp;\n" \
"		results[tmpi].F2 = 0;\n" \
"	}\n" \
"}\n" \
"}";



/////////////////////



char *load_program_source(const char *filename)
{ 
	FILE *KERNEL = fopen(filename, "r");
	char *source; 

	if (KERNEL == 0)
		return 0; 
	
	source = (char *) malloc(0x100000);
	fread(source, 1, 0x100000, KERNEL);
	fclose(KERNEL);
	return source;
} 


void ereport(int ret, int code)
{
	if(ret == CL_DEVICE_NOT_FOUND)
	{
		printf("\n\nError:\nThere is no suitable graphics device on this computer.\nOpenCL requires NVIDIA GeForce 8 series or above, or ATI Radeon 4000 series or above.\n\n");
		exit(1);
	}
	if(ret != CL_SUCCESS)
	{
		printf("\n\nThere was an error. \nError line %d\nCode %d\n",code,ret);
		exit(code);
	} else {
		if(silent == 1)
		{
			printf("%d\n",code);
			fflush(stdout);
		}
	}
}


/*
Function descripton:
If the Checks the user's input file if the analysis is being resumed. If it exists and a scan has already started then the
analysis will resume from the last chromosome x chromosome scan provided that the data in the input file hasn't changed

*/

void resume(char *filename, int *chr1, int *chr2, int *tothits, int nid, int nsnp)
{
	FILE *OUT;
	int i,j,k,c1,c2;
	float thresh1, thresh2;
	char temp1[300], x;

	if((OUT = fopen(filename,"rt")))
	{
		fgets(temp1,300,OUT);
		sscanf(temp1,"Datafile: %s\n",&*temp1);
		k = strcmp(temp1,plinkbin);
		if(k != 0)
		{
			printf("Can't resume analysis with different datafile!\n");
			exit(1);
		}
		fgets(temp1,300,OUT);
		sscanf(temp1,"Number of individuals: %d\n",&i);
		if(i != nid)
		{
			printf("Can't resume analysis with different data!\n");
			exit(1);
		}

		fgets(temp1,300,OUT);
		sscanf(temp1,"Number of SNPs: %d\n",&i);
		if(i != nsnp)
		{
			printf("Can't resume analysis with different data!\n");
			exit(1);
		}

		fgets(temp1,300,OUT);
		sscanf(temp1,"Permutation: %d\n",&i);
		if(i != UPERM)
		{
			printf("Can't resume analysis with different permutation!\n");
			exit(1);
		}

		fgets(temp1,300,OUT);
		sscanf(temp1,"Threshold: %f %f\n",&thresh1, &thresh2);
		if(thresh1 != UTHRESHOLD1)
		{
			printf("Warning: Resuming analysis with different full threshold, changed from %f to %f\n",thresh1,UTHRESHOLD1);
		}
		if(thresh2 != UTHRESHOLD2)
		{
			printf("Warning: Resuming analysis with different interaction threshold, changed from %f to %f\n",thresh2,UTHRESHOLD2);
		}

		fgets(temp1,300,OUT);
		sscanf(temp1,"Test: %c",&x);
		if(x != UTEST)
		{
			printf("Can't resume analysis performing a different test!\n");
			exit(1);
		}

		fgets(temp1,300,OUT);
		fgets(temp1,300,OUT);
		j = 0;
		while(fgets(temp1,300,OUT) != NULL)
		{
			j++;
			sscanf(temp1,"# %d x %d : %d interactions\n", &c1, &c2, &k);
			*tothits += k;
			for(i = 0; i < k; i++)
			{
				fgets(temp1,300,OUT);
			}
		}
		if(j == 0)
		{
			*chr1 = 0; *chr2 = 0;
		} else {
			if(c1 == c2)
			{
				*chr1 = c1+1; *chr2 = 0;
			} else {
				*chr1 = c1; *chr2 = c2+1;
			}
		}
	} else {
		*chr1 = 0; *chr2 = 0;
	}
}

/*

Function description:
Main analysis routine (called by main)
The full scan is represented as a lower triangular matrix. This is partitioned into chromosome segments, 
such that the triangle is diced up into rectangles that represent chromosome i x chromsome !i, and triangles
that represent chromosome i x chromosome i scans.
Each chromosome x chromosome scan is then split up into squares of size 'chunksize' (e.g. 512x512 SNPs). Each
square is a single kernel execution (must last < 5sec on GUI OS)
The local work size is set to 16x16, this is because anything larger would require more local memory to store 
genotype class means / counts (this could be improved by using an on-line algorithm to store only class means, 
removing the need for short int class counts, but the cost in speed could be prohibitive).
The global work size of chunksize x chunksize is distributed across the GPU such that streaming multiprocessors
queue to receive 16x16 blocks.
Each local work thread has its own global results array, and global result counter. If a thread calculates
an F value above some threshold then it increments its results counter and saves the result as a table struct.
These results are downloaded and reset between each chromosome x chromosome scan



1. Calculate mean / SSY of phenotype
2. Check device availability
3. Check input data and output file (resume or fresh start)
4. Create context, command queue, kernel
5. Allocate memory on device, write phenotype
6. Set kernel arguments that are constant for all kernel executions
7. Perform scan, partitioned chromosome x chromosome
 a. Write chromosome genotype data to device, set boundaries as kernel arguments
 b. Local work size is 2D {16,16}, allocate space for results arrays and initialise results counters to 0
 c. Partition chromosome x chromosome scan into chunksize x chunksize size square blocks
 d. Set block specific kernel arguments
 e. Iterate through blocks, sequentially executing kernels for these SNP ranges
 f. Append results to output file

*/

int perform_scan(int *geno, chromosome *chrstat, ped *dat, map *genmap, int nid, int nsnp, int npack, int remain, int nchr, char *filename)
{
	cl_context context;
	cl_context_properties props[3];
	cl_command_queue cmd_queue;
	cl_device_id devices;
	
	cl_device_id devs[100];
	cl_uint ndevs = 0;
	cl_bool av[10];

	cl_platform_id platform;
	cl_uint ret_num_platforms;
	cl_int err;
	cl_program program1;
	cl_kernel kernel1;
	cl_event event;
	cl_ulong start, end;
	double kerneltime;
	char timetaken[50];
	char build[2048], devname[2048];
	int tnhits = 0;
	int *nhits = (int *)malloc(sizeof(int)*256);
	int tothits = 0;
	cl_mem g_chr1, g_chr2, g_phen, g_nhits, g_results;
	float *phen = (float *)malloc(sizeof(float) * nid);
	float mY = 0, SSY = 0;
	int s1, s2, i, j, k, l, npack1, ccount, chrint, flag;
	int niterations0, niterations1, chunksize = UCHUNKSIZE, *chunks0, *chunks1;
	int chr1, chr1g, chr2, chr2g, chr1m, chr2m;
	time_t oval,nval;
	double tres;
	FILE *OUT;
	size_t gws2d[2], offset;
	size_t lws2d[2];
	size_t factor_count_size;
	size_t mY_fac_size;
	table *results = (table *)malloc(sizeof(table) * MAX);
	oval = time(NULL);
	
	if(USAFE)
	{
		lws2d[0] = lws2d[1] = 8;
		factor_count_size = lws2d[0]*lws2d[1]*9*sizeof(unsigned int);
		mY_fac_size = lws2d[0]*lws2d[1]*9*sizeof(float);
	} else {
		lws2d[0] = lws2d[1] = 16;
		factor_count_size = lws2d[0]*lws2d[1]*9*sizeof(unsigned short int);
		mY_fac_size = lws2d[0]*lws2d[1]*9*sizeof(float);
	}
	
	// Assuming no missing values the phenotype mean and ssq will remain constant
	// Therefore only needs to be calculated once.
	// If missing values exist then on-line algorithm adjusts mean and ssq at runtime in kernel
	for(i = 0; i < nid; i++)
	{
		phen[i] = dat[i].phen;
		mY += phen[i];
	}
	mY /= nid;
	for(i = 0; i < nid; i++)
	{
		SSY += (phen[i] - mY) * (phen[i] - mY);
	}
	for(i = 0; i < 256; i++)
	{
		nhits[i] = 0;
	}
	
	// Platform information
	err = clGetPlatformIDs(1,&platform,&ret_num_platforms);   ereport(err,__LINE__);
	props[0] = (cl_context_properties)CL_CONTEXT_PLATFORM;
	props[1] = (cl_context_properties)platform;
	props[2] = (cl_context_properties)0;
	err = clGetDeviceIDs(platform,CL_DEVICE_TYPE_GPU, 1, &devices, NULL);   ereport(err,__LINE__);
	
	// Device information
	clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 100, devs, &ndevs);
	printf("Devices (availability):\n");
	for(i = 0; (unsigned int)i < ndevs; i++)
	{
		clGetDeviceInfo(devs[i],CL_DEVICE_AVAILABLE,sizeof(cl_bool),&av[i],NULL);
		clGetDeviceInfo(devs[i],CL_DEVICE_NAME,sizeof(devname),devname,NULL);
		printf("%d:\t%s (%d)\n",i,devname,av[i]);
	}
	// Check if device exists
	if((ndevs-1) < (unsigned int)UDEVICE)
	{
		printf("Invalid device choice\nPlease select one of the listed devices\n\n");
		exit(1);
	}
	// Check if device is available (possible bug in OpenCL 1.0?)
	if(av[UDEVICE] != 1)
	{
		printf("The selected device is unavailable, please select a different one\n\n");
		exit(1);
	}

	printf("\n\n");

	// Check if this is a fresh analysis, or if it is resuming a previously started analysis.
	// If it is resuming then import previous results and check parameters and data are the same
	resume(filename, &s1, &s2, &tothits, nid, nsnp);
	ccount = 0;
	chrint = nchr * (nchr+1) / 2; // total number of chromosome x chromosome scans to perform
	if(s1 != 0 || s2 != 0)
	{
		ccount = s1 * (s1+1) / 2 + s2;
		if(ccount == chrint)
		{
			printf("Scan already completed!\n\n");
			return tothits;
		} else if(ccount > chrint) {
			printf("The initial analysis was performed with a different dataset!\nPlease use a different filename.\n");
			exit(1);
		}		
		printf("Resuming from %d/%d, %d hits found so far...\n\n", ccount+1, chrint, tothits);
	} else {
		// If fresh analysis then record analysis details and write table header in output file
		OUT = fopen(filename, "wt");
		fprintf(OUT,"Datafile: %s %s\n",plinkbin, UPHEN);
		fprintf(OUT,"Number of individuals: %d\nNumber of SNPs: %d\nPermutation: %d\nThreshold: %f %f\n",nid,nsnp,UPERM, UTHRESHOLD1, UTHRESHOLD2);
		if(UTEST == 'f')
		{
			fprintf(OUT,"Test: full\n\n");
			fprintf(OUT,"SNP1\tSNP1name\tChr1\tSNP2\tSNP2name\tChr2\tdf1\tdf2\tFval\n");
		} else {
			fprintf(OUT,"Test: interactions + full\n\n");
			fprintf(OUT,"SNP1\tSNP1name\tChr1\tSNP2\tSNP2name\tChr2\tdf1\tdf2\tFval\tFint\n");
		}
		fclose(OUT);
	}

	// Create the context
	context = clCreateContextFromType(props, CL_DEVICE_TYPE_GPU, NULL, NULL, &err);
	 ereport(err,__LINE__);

	clGetDeviceInfo(devs[UDEVICE],CL_DEVICE_NAME,sizeof(devname),devname,NULL);
	printf("Using device %d: %s\n\n", UDEVICE, devname);

	// Create the command queue
	cmd_queue = clCreateCommandQueue(context, devs[UDEVICE], CL_QUEUE_PROFILING_ENABLE, NULL);  ereport(err,__LINE__);

	// Compile kernel. 
	if(UTEST == 'f')
	{
		if(USAFE)
		{
			program1 = clCreateProgramWithSource(context,1, (const char **)&kernelfullsafe, NULL, &err);  ereport(err,__LINE__);			
		} else {
			program1 = clCreateProgramWithSource(context,1, (const char **)&KERN, NULL, &err);  ereport(err,__LINE__);
		}
	} else {
		if(USAFE)
		{
			program1 = clCreateProgramWithSource(context,1, (const char **)&kernelintsafe, NULL, &err);  ereport(err,__LINE__);			
		} else {
			program1 = clCreateProgramWithSource(context,1, (const char **)&kernelint, NULL, &err);  ereport(err,__LINE__);
		}
	}
	err = clBuildProgram(program1, 0, NULL, NULL, NULL, NULL);  ereport(err,__LINE__);
	clGetProgramBuildInfo(program1, devs[UDEVICE], CL_PROGRAM_BUILD_LOG, 2048, build, NULL);  ereport(err,__LINE__);
	//printf("\n%s\n\n",build);

	kernel1 = clCreateKernel(program1, "square", &err);  ereport(err,__LINE__);
	npack1 = remain ? npack + 1 : npack;

	g_phen = clCreateBuffer(context, CL_MEM_READ_ONLY, nid*sizeof(float), NULL, &err);  ereport(err,__LINE__);
	err = clEnqueueWriteBuffer(cmd_queue, g_phen, CL_TRUE, 0, nid*sizeof(float), phen, 0, NULL, NULL);  ereport(err,__LINE__);
	g_results = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(table)*MAX*lws2d[0]*lws2d[1], NULL, &err);  ereport(err,__LINE__);
	g_nhits = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int)*lws2d[0]*lws2d[1], NULL, &err);  ereport(err,__LINE__);

	err = clSetKernelArg(kernel1, 2, sizeof(cl_mem), (void *)&g_phen);  ereport(err,__LINE__);
	err = clSetKernelArg(kernel1, 3, sizeof(cl_mem), (void *)&g_results);  ereport(err,__LINE__);
	err = clSetKernelArg(kernel1, 4, sizeof(cl_float), &mY);  ereport(err,__LINE__);
	err = clSetKernelArg(kernel1, 5, sizeof(cl_float), &SSY);  ereport(err,__LINE__);
	err = clSetKernelArg(kernel1, 6, sizeof(cl_float), &UTHRESHOLD1);  ereport(err,__LINE__);
	err = clSetKernelArg(kernel1, 7, sizeof(cl_float), &UTHRESHOLD2);  ereport(err,__LINE__);
	err = clSetKernelArg(kernel1, 8, sizeof(cl_int), &nid);  ereport(err,__LINE__);
	err = clSetKernelArg(kernel1, 9, sizeof(cl_int), &nsnp);  ereport(err,__LINE__);
	err = clSetKernelArg(kernel1, 10, sizeof(cl_int), &npack);  ereport(err,__LINE__);
	err = clSetKernelArg(kernel1, 11, sizeof(cl_int), &remain);  ereport(err,__LINE__);
	err = clSetKernelArg(kernel1, 12, mY_fac_size, NULL);  ereport(err,__LINE__);
	err = clSetKernelArg(kernel1, 13, factor_count_size, NULL);  ereport(err,__LINE__);
	err = clSetKernelArg(kernel1, 14, sizeof(cl_mem), (void *)&g_nhits);  ereport(err,__LINE__);

	
	for(i = s1; i < nchr; i++)
	{
		// upload first chromosome and set argument
		g_chr1 = clCreateBuffer(context, CL_MEM_READ_ONLY, npack1*chrstat[i].chrsize*sizeof(int), NULL, &err);  ereport(err,__LINE__);
		err = clEnqueueWriteBuffer(cmd_queue, g_chr1, CL_TRUE, 0, npack1*chrstat[i].chrsize*sizeof(int), &geno[npack1*chrstat[i].chrstart], 0, NULL, NULL);  ereport(err,__LINE__);
		err = clSetKernelArg(kernel1, 0, sizeof(cl_mem), (void *)&g_chr1);  ereport(err,__LINE__);

		// define boundary for chr1
		chr1m = chrstat[i].chrsize;

		// chr1g refers chromosome start position, gws2d refers to iteration SNP position
		chr1g = chrstat[i].chrstart;

		for(j = s2; j <= i; j++)
		{
			//upload second chromosome and set argument
			g_chr2 = clCreateBuffer(context, CL_MEM_READ_ONLY, npack1*chrstat[j].chrsize*sizeof(int), NULL, &err);  ereport(err,__LINE__);
			err = clEnqueueWriteBuffer(cmd_queue, g_chr2, CL_TRUE, 0, npack1*chrstat[j].chrsize*sizeof(int), &geno[npack1*chrstat[j].chrstart], 0, NULL, NULL);  ereport(err,__LINE__);
			err = clSetKernelArg(kernel1, 1, sizeof(cl_mem), (void *)&g_chr2);  ereport(err,__LINE__);

			// define boundary for chr2
			chr2m = chrstat[j].chrsize;

			// reset hit count locally and on GPU
			for(k = 0; (unsigned int)k < lws2d[0]*lws2d[1]; k++)
			{
				nhits[k] = 0;
			}
			err = clEnqueueWriteBuffer(cmd_queue, g_nhits, CL_TRUE, 0, sizeof(int)*lws2d[0]*lws2d[1], nhits, 0, NULL, NULL);  ereport(err,__LINE__);
			printf("%d/%d\tScanning %s vs %s\n",++ccount,chrint,chrstat[i].chrname,chrstat[j].chrname);
			// begin chrxchr breakdown
			chr1 = 0;
			gws2d[0] = (size_t)ceil((double)chrstat[i].chrsize / lws2d[0]) * lws2d[0];
			gws2d[1] = (size_t)ceil((double)chrstat[j].chrsize / lws2d[1]) * lws2d[1];
			if(i != j) { printf("\tWork size: %lu x %lu square\n\n  ",gws2d[0],gws2d[1]); 
			} else { printf("\tWork size: %lu x %lu triangle\n\n  ",gws2d[0],gws2d[1]);}

			niterations0 = (int)ceil((double)gws2d[0] / chunksize);
			niterations1 = (int)ceil((double)gws2d[1] / chunksize);
			chunks0 = (int *)malloc(sizeof(int) * niterations0);
			chunks1 = (int *)malloc(sizeof(int) * niterations1);
			for(k = 0; k < niterations0; k++) chunks0[k] = chunksize;
			for(k = 0; k < niterations1; k++) chunks1[k] = chunksize;
			l = gws2d[0] % chunksize;
			if(l) chunks0[(niterations0-1)] = l;
			l = gws2d[1] % chunksize;
			if(l) chunks1[(niterations1-1)] = l;
	
			kerneltime = 0;
			for(k = 0; k < niterations0; k++)
			{
				if(i == j)
				{
					flag = k+1;
				} else {
					flag = niterations1;
				}
				chr2g = chrstat[j].chrstart;
				chr2 = 0;
				gws2d[0] = chunks0[k];
				for(l = 0; l < flag; l++)
				{
					gws2d[1] = chunks1[l];
					printf(".");fflush(stdout);
					err = clSetKernelArg(kernel1, 15, sizeof(cl_int), &chr1);  ereport(err,__LINE__);
					err = clSetKernelArg(kernel1, 16, sizeof(cl_int), &chr2);  ereport(err,__LINE__);
					err = clSetKernelArg(kernel1, 17, sizeof(cl_int), &chr1g);  ereport(err,__LINE__);
					err = clSetKernelArg(kernel1, 18, sizeof(cl_int), &chr2g);  ereport(err,__LINE__);
					err = clSetKernelArg(kernel1, 19, sizeof(cl_int), &chr1m);  ereport(err,__LINE__);
					err = clSetKernelArg(kernel1, 20, sizeof(cl_int), &chr2m);  ereport(err,__LINE__);
					err = clEnqueueNDRangeKernel(cmd_queue, kernel1, 2, NULL, gws2d, lws2d, 0, NULL, &event);  ereport(err,__LINE__);
					err = clWaitForEvents(1,&event);  ereport(err,__LINE__);
					err = clFinish(cmd_queue);  ereport(err,__LINE__);
					clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
					clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL);
					kerneltime += (double)(end - start) / 1e9;
					chr2 += chunksize;
				}
				chr1 += chunksize;
				printf("\n  ");fflush(stdout);
			}
			clReleaseMemObject(g_chr2);
			err = clEnqueueReadBuffer(cmd_queue, g_nhits, CL_TRUE, 0, sizeof(int)*lws2d[0]*lws2d[1], nhits, 0, NULL, NULL);  ereport(err,__LINE__);
			readtime(timetaken,kerneltime);
			tnhits = 0;
			for(k = 0;(unsigned int) k < lws2d[0]*lws2d[1]; k++)
			{
				tnhits += nhits[k];
			}
			OUT = fopen(filename, "a");
			printf("\n\t%d interactions found in %s\n\n",tnhits,timetaken);
			fprintf(OUT,"# %d x %d : %d interactions\n", i,j, tnhits);
			if(tnhits > 0)
			{
				for(k = 0; (unsigned int)k < lws2d[0]*lws2d[1]; k++)
				{
					if(nhits[k])
					{
						offset = MAX2 * k * sizeof(table);
						err = clEnqueueReadBuffer(cmd_queue, g_results, CL_TRUE, offset, sizeof(table)*nhits[k], results, 0, NULL, NULL);  ereport(err,__LINE__);
						err = clFinish(cmd_queue);  ereport(err,__LINE__);
						for(l = 0; l < nhits[k]; l++)
						{
							if(UTEST == 'f')
							{
	fprintf(OUT, "%d\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%f\n", \
		results[l].snp1+1 + chrstat[i].chrstart, \
		genmap[results[l].snp1 + chrstat[i].chrstart].snpname, \
		chrstat[i].chrname, \
		results[l].snp2+1 + chrstat[j].chrstart, \
		genmap[results[l].snp2 + chrstat[j].chrstart].snpname, \
		chrstat[j].chrname, \
		results[l].df1, \
		results[l].df2, \
		results[l].F);
							} else {
	fprintf(OUT, "%d\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%f\t%f\n",
		results[l].snp1+1 + chrstat[i].chrstart, \
		genmap[results[l].snp1 + chrstat[i].chrstart].snpname, \
		chrstat[i].chrname, \
		results[l].snp2+1 + chrstat[j].chrstart, \
		genmap[results[l].snp2 +  chrstat[j].chrstart].snpname, \
		chrstat[j].chrname, \
		results[l].df1, \
		results[l].df2, \
		results[l].F, \
		results[l].F2);
							}
							tothits++;
						}
					}
				}
			}
			fclose(OUT);
			free(chunks0);
			free(chunks1);
		}
		clReleaseMemObject(g_chr1);
	}

//	fclose(OUT);
	clReleaseKernel(kernel1);
	clReleaseProgram(program1);
	clReleaseCommandQueue(cmd_queue);
	clReleaseContext(context);
	
	//clReleaseMemObject(g_geno);
	clReleaseMemObject(g_phen);
	clReleaseMemObject(g_results);
	clReleaseMemObject(g_nhits);
	free(results);
	free(nhits);
	free(phen);

	nval = time(NULL);
	timediff(&tres,nval,oval);
	readtime(timetaken,tres);
	printf("\nProgramme finished in %s\n%d epistatic interactions discovered\n", timetaken, tothits);
	return tothits;
}

