#include "readplink.h"

// find text file table dimensions
// Input files should have the same number of columns for each row
// - otherwise programme exits
// Columns are counted by character spaces (ASCII space or tab)
// robust to double spaces/tabs
void dims(int *rowCount, int *colCount, char *filename)
{
	int cur, row=0, col=0, flag = 0;
	char sep;
	FILE *ifp = fopen(filename,"r");
	for(;;)
	{
		fscanf(ifp,"%c",&sep);
		if(sep == ' ' || sep == '\t' || sep == EOF) break;
	}
	fclose(ifp);
	if(sep == EOF)
	{
		printf("Problem with data, must be space or tab delimited\n\n");
		exit(1);
	}
	ifp = fopen(filename,"r");

	printf("Getting %s dimensions...",filename);fflush(stdout);
	while ((cur = fgetc(ifp)) != EOF)
	{
		if (cur == '\n')
		{
			row++;
			flag = 0;
			//col++;
			// set col to be 
			if(row == 1)
			{
				*colCount = col;
				col = 0;
			} else {
				if(col != *colCount)
				{
					printf("\nRagged matrix:\nLine %d\n%d %d\nExiting...\n",row,col,*colCount);
					exit(1);
				} else {
					col = 0;
				}
			}
		} else if (cur == sep) {
			flag = 0;
		} else if (flag == 0) {
			col++;
			flag = 1;
		}
	}
	*rowCount = row;
	*colCount = (*colCount - 6) / 2;
	printf("done!\n%d SNPs, %d Individuals\n\n",*colCount, *rowCount);fflush(stdout);

	fclose(ifp);
}

void walkersalias(int n, float *p, int nans, int *ans)

{

	float *q, rU;

	int i,j,k;

	int *HL,*H,*L;

	int *a;

	float sum = 0;



	a = (int *)malloc(n*sizeof(int));

	HL = (int *)malloc(n*sizeof(int));

	q = (float *)malloc(n*sizeof(float));

	H = HL - 1; L = HL + n;

	

	for(i = 0; i < n; i++)

	{

		sum += p[i];

	}

	for(i = 0; i < n; i++)

	{

		p[i] /= sum;

	}

	for(i = 0; i < n; i++)

	{

		q[i] = p[i] * n;

		if(q[i] < 1.) *++H = i; else *--L = i;

	}

	if(H >= HL && L < HL +n)

	{

		for(k = 0; k < n-1; k++)

		{

			i = HL[k];

			j = *L;

			a[i] = j;

			q[j] += q[i] - 1;

			if(q[j] < 1.) L++;

			if(L >= HL + n) break;

		}

	}

	for(i = 0; i < n; i++) q[i] += i;

	for(i = 0; i < nans; i++)

	{

		rU = (float) rand() / RAND_MAX * n;

		k = (int) rU;

		ans[i] = (rU < q[k]) ? k : a[k];

	}

	free(HL);

	free(q);

	free(a);

}



char isbinary(int nid, ped *dat)
{
	char flag = 0;
	int count[3], i;
	count[0] = count[1] = count[2] = 0;
	for(i = 0; i < nid; i++)
	{
		if(dat[i].phen == 1)
		{
			count[0]++;
		} else if(dat[i].phen == 2) {
			count[1]++;
		} else {
			count[2]++;
		}
	}

	if(count[2] == 0)
	{
		if((count[0] + count[1]) == nid)
		{
			flag = 1;
		} else {
			printf("Phenotype error!\n");
			exit(1);
		}
	}

	return flag;
}

// read in epiGPU binary format and write to PLINK format
void extractfrombinary(char *binfile, char *pedfile, char *mapfile)
{
	int i,j, nid, nsnp, npack, remain, nchr;
	int flag;
	map *genmap;
	ped *dat;
	int *genop;
	char *geno, g;
	chromosome *chrstat;
	FILE *map, *ped;

	// read and unpack
	readpackedbinary(&nid, &nsnp, &npack, &remain, &nchr, &genmap, &dat, &genop, &chrstat, binfile);
	printf("Unpacking...");
	unpack(nid, nsnp, npack, remain, genop, &geno);
	free(genop);
	
	// write out mapfile
	// 4 Columns: chromosome | SNP name | genetic distance | physical distance
	map = fopen(mapfile, "w");
	for(i = 0; i < nsnp; i++)
	{
		fprintf(map,"%s %s %f %ld\n", genmap[i].chr, genmap[i].snpname, genmap[i].gd, genmap[i].pd);
	}
	fclose(map);
	
	// check if data is binary or continuous
	flag = isbinary(nid,dat);

	// write out pedfile
	// first 6 columns: family | ID | father | mother | sex | phenotype
	
	if(flag)
	{
		printf("Writing to files (binary phenotype)...");
	} else {
		printf("Writing to files (continuous phenotype)...");
	}
	ped = fopen(pedfile, "w");
	for(i = 0; i < nid; i++)
	{
		if(flag)
		{
			fprintf(ped,"%s %s %s %s %d %d",dat[i].family, dat[i].id, dat[i].paternal, dat[i].maternal, (int)dat[i].sex, (int)dat[i].phen);
		} else {
			fprintf(ped,"%s %s %s %s %d %f",dat[i].family, dat[i].id, dat[i].paternal, dat[i].maternal, (int)dat[i].sex, dat[i].phen);
		}
		for(j = 0; j < nsnp; j++)
		{
			g = geno[j*nid+i];
			switch(g)
			{
				case 0:
					fprintf(ped," %c %c", genmap[j].allelename[0], genmap[j].allelename[0]);
					break;
				case 1:
					fprintf(ped," %c %c", genmap[j].allelename[0], genmap[j].allelename[1]);
					break;
				case 2:
					fprintf(ped," %c %c", genmap[j].allelename[1], genmap[j].allelename[1]);
					break;
				case 3:
					fprintf(ped," 0 0");
					break;
			}
		}
		fprintf(ped,"\n");
	}
	fclose(ped);
	printf("done\n\n");
	free(geno);
	free(genmap);
	free(dat);
	free(chrstat);
}

// Read text file in PLINK format (.ped)
// first 6 columns: family | ID | father | mother | sex | phenotype
// Remaining columns: 2 columns per SNP
// Convert 2 character allele code to single character genotype
// missing values = 0
// first homozygote allele <- allele1 <- 0
// store true allele bases in genmap
// Columns can be separated by spaces or tabs, double spaces will be ignored
void readped(int nid, int nsnp, char *pedfile, ped *dat, char *geno, map *genmap)
{
	int i,j,k,l;
	FILE *PED;
	char allele1 = 1, allele2 = 1;
	char *buffer = (char *)malloc(sizeof(char) * nsnp * 2);
	
	for(i = 0; i < nsnp; i++)
	{
		genmap[i].allelename[0] = genmap[i].allelename[1] = 0;
	}
	PED = fopen(pedfile,"r");
	printf("Reading in %s.",pedfile);fflush(stdout);
	for(i = 0; i < nid; i++)
	{
		if(i % ((int)nid/10) == 1) {printf(".");fflush(stdout);}
		(void)fscanf(PED,"%s",dat[i].family);
		(void)fscanf(PED,"%s",dat[i].id);
		(void)fscanf(PED,"%s",dat[i].paternal);
		(void)fscanf(PED,"%s",dat[i].maternal);
		(void)fscanf(PED,"%d",&j);
		(void)fscanf(PED,"%f%*[ /t]",&dat[i].phen);
		dat[i].sex = j;
		// store individual's genotype data in 'buffer'
		for(j = 0; j < (nsnp*2); j++)
		{
			(void)fscanf(PED,"%c%*[ /t]",&buffer[j]);
		}
		
		// recode 'buffer' and store in correct position in genotype array
		// in this round, recode heterozygotes as 1, and homozygotes as single instance of allele
		k = 0;
		for(j = 0; j < (nsnp*2); j+=2)
		{
			if(buffer[j] == '0')
			{
				geno[(k*nid+i)] = MISS;
			} else {
				if(buffer[j] != buffer[j+1])
				{
					genmap[(j/2)].allelename[0] = buffer[j];
					genmap[(j/2)].allelename[1] = buffer[j+1];

					geno[(k*nid+i)] = 1;
				} else {
					geno[(k*nid+i)] = buffer[j];
				}
			}
			k++;
		}
	}
	fclose(PED);
	printf("done!\nConverting..");fflush(stdout);

	// Convert homozygote codes
	for(i = 0; i < nsnp; i++)
	{
		if(i % ((int)nsnp/10) == 1) { printf(".");fflush(stdout);}
		allele1 = 1;
		allele2 = 1;
		j = 0;
		// find the first homozygote, call this allele1
		while(allele1 == 1 || allele1 == MISS)
		{
			allele1 = geno[(i*nid + j++)];
		}
		if(genmap[i].allelename[0] == 0)
		{
			genmap[i].allelename[0] = allele1;
		} else if(allele1 == genmap[i].allelename[1]) {
			genmap[i].allelename[1] = genmap[i].allelename[0];
			genmap[i].allelename[0] = allele1;
		}
		k = j;
		// find the second homozygote, call this allele2 (may not exist)
		while((allele2 == 1 || allele2 == allele1 || allele2 == MISS) && k < nid)
		{
			allele2 = geno[(i*nid + k++)];
		}
		// convert homozygotes to 0 or 2
		if(allele2 != allele1 && allele2 != 1 && allele2 != MISS)
		{
			for(l = j-1; l < nid; l++)
			{
				if(geno[(i*nid + l)] == allele1)
				{
					geno[(i*nid + l)] = 0;
				} else if (geno[(i*nid+l)] == allele2) {
					geno[(i*nid + l)] = 2;
				}
			}
		} else {
			for(l = j-1; l < nid; l++)
			{
				if(geno[(i*nid + l)] == allele1)
				{
					geno[(i*nid + l)] = 0;
				}
			}
		}
	}
	free(buffer);
	printf("done!\n");fflush(stdout);

//	for(i = 0; i < 10; i++)
//	{
//		for(j = 0; j < 10; j++)
//		{
//			printf("%d ", geno[(i*nid+j)]);
//		}
//		printf("\n");
//	}
	
}

// Read in PLINK map file (.map)
// 4 Columns: chromosome | SNP name | genetic distance | physical distance
// still to do: sort map by physical distance
void readmap(int nsnp, map* genmap, int *nchr, chromosome *chrstat, char *mapfile)
{
	int i, j, k;
	FILE *MAP;
	char sep;
	MAP = fopen(mapfile,"r");
	printf("Reading in %s..",mapfile);fflush(stdout);
	*nchr = 0;
	for(;;)
	{
		fscanf(MAP,"%c",&sep);
		if(sep == ' ' || sep == '\t' || sep == EOF) break;
	}
	fclose(MAP);
	MAP = fopen(mapfile,"r");
	
	for(i = 0; i < nsnp; i++)
	{
		if(i % ((int)nsnp/10) == 1){ printf(".");fflush(stdout);}

		(void)fscanf(MAP,"%s",genmap[i].chr);			
		k = 0;
		for(j = 0; j < *nchr; j++)
		{
			if(strcmp(genmap[i].chr, chrstat[j].chrname) == 0)
			{
				(chrstat[j].chrsize)++;
				k = 1;
				break;
			}
		}
		if(k == 0)
		{
			chrstat[(*nchr)].chrsize = 1;
			chrstat[(*nchr)].chrid = (*nchr);
			while(genmap[i].chr[k] != '\0')
			{
				chrstat[(*nchr)].chrname[k] = genmap[i].chr[k];
				k++;
			}
			chrstat[(*nchr)].chrname[k] = '\0';
			(*nchr)++;
		}

		(void)fscanf(MAP,"%s",genmap[i].snpname);
		(void)fscanf(MAP,"%lf ",&genmap[i].gd);
		(void)fscanf(MAP,"%li ",&genmap[i].pd);
	}
	fclose(MAP);

	chrstat[0].chrstart = 0;
	for(i = 1; i < *nchr; i++)
	{
		chrstat[i].chrstart = chrstat[i-1].chrstart + chrstat[i-1].chrsize;
	}
	
	printf("done!\n");fflush(stdout);
	printf("%d chromosomes\n",*nchr); fflush(stdout);
	for(i = 0; i < *nchr; i++)
	{
		printf("%s\t%d\t%d\n",chrstat[i].chrname,chrstat[i].chrsize,chrstat[i].chrstart);fflush(stdout);
	}
}

// Remove SNPs/IDs with low callrates (defined in readplink.h as MINAF and MINCALL)
void cleangeno(int *nid, int *nsnp, int nchr, char *geno, ped *dat, chromosome *chrstat, map *genmap)
{
	float tmpf1;
	int i,j,k,nsnp1, nid1, sum1;
	float *maf = (float *)malloc(*nsnp * sizeof(float));
	int *SNPs = (int *)malloc(*nsnp * sizeof(int));
	int *IDs = (int *)malloc(*nid * sizeof(int));
	printf("Cleaning genotypes:\n%0.3f maf cutoff\n%0.3f per cent minimum callrate\nCleaning SNPs..",MINAF,MINCALL*100); fflush(stdout);
	
	nsnp1 = 0;
	for(i = 0; i < *nsnp; i++)
	{
		if(i % ((int)*nsnp/10) == 1){ printf(".");fflush(stdout); }
		sum1 = 0;
		k = 0;
		for(j = 0; j < *nid; j++)
		{
			if(geno[(*nid * i + j)] != MISS)
			{
				sum1 += (geno[(*nid * i + j)]);
				k++;
			}
		}
		tmpf1 = (float)k / (float)*nid;
		maf[i] = (float)sum1 / (float)(2*k);
		if((maf[i] > MINAF) && (tmpf1 > MINCALL) && (genmap[i].pd > 0)) SNPs[nsnp1++] = i;
	}
	printf("done!\n%d SNPs remain\nChecking individuals..",nsnp1);fflush(stdout);

	nid1 = 0;
	for(i = 0; i < *nid; i++)
	{
		if(i % ((int)*nid/10) == 1) { printf(".");fflush(stdout);}
		k = 0;
		for(j = 0; j < nsnp1; j++)
		{
			if(geno[(*nid * SNPs[j] + i)] != MISS) k++;
		}
		tmpf1 = (float)k / (float)nsnp1;
		if(tmpf1 > MINCALL && dat[i].phen != -9) IDs[nid1++] = i;
	}
	printf("done!\n%d individuals remain\n",nid1);fflush(stdout);

	k = 0;
	for(i = 0; i < nsnp1; i++)
	{
		genmap[i] = genmap[SNPs[i]];
		for(j = 0; j < nid1; j++)
		{
			geno[k++] = geno[(*nid * SNPs[i] + IDs[j])];
		}
	}
	for(i = 0; i < nid1; i++)
		dat[i] = dat[IDs[i]];

	*nid = nid1;
	*nsnp = nsnp1;

	for(i = 0; i < nchr; i++)
	{
		chrstat[i].chrsize = 0;
	}

	for(i = 0; i < nsnp1; i++)
	{
		for(j = 0; j < nchr; j++)
		{
			if(strcmp(genmap[i].chr, chrstat[j].chrname) == 0)
			{
				chrstat[j].chrsize++;
				break;
			}
		}
	}
	chrstat[0].chrstart = 0;
	for(i = 1; i < nchr; i++)
	{
		chrstat[i].chrstart = chrstat[(i-1)].chrstart + chrstat[(i-1)].chrsize;
	}
}

// replace all 'MISS's with genotype, sample proportional to allele frequency
void guessmissing(int nid, int nsnp, char *geno)
{
	int i,j,temp,count=0;
	printf("Guessing missing values..");fflush(stdout);
	for(i = 0; i < nsnp; i++)
	{
		if(i % ((int)nsnp/10) == 1) { printf(".");fflush(stdout);}
		for(j = 0; j < nid; j++)
		{
			if(geno[(nid * i + j)] == MISS)
			{
				count++;
				while(geno[(nid * i + j)] == MISS)
				{
					temp = (int)rand() % nid;
					geno[(nid * i + j)] = geno[(nid * i + temp)];
				}
			}
		}
	}
	printf("done!\n%d missing values removed\n",count);
}

void contphen(int nid, ped *dat)
{
	int i;
	for(i = 0; i < nid; i++)
	{
		dat[i].phen = (float) rand() / RAND_MAX;
	}
}

// Box-Muller method approximation. mean = 0; var = 1
// Marsaglia polar method may be faster
void rnorm(int n, ped *dat)
{
	int
		i;
	float
		tmp1f, tmp2f,
		x1,x2;

	for(i = 0; i < n; i++)
	{
		x1 = (float)rand() / RAND_MAX;
		x2 = (float)rand() / RAND_MAX;
		tmp1f = sqrtf(-2 * logf(x1));
		tmp2f = cosf(2*PI*x2);
		dat[i].phen = tmp1f * tmp2f;
	}
	printf("Generating random phenotype...\n");
}


//int nid | int nsnp | int nchr | ped[nid] | map[nsnp] | char[nid*nsnp]
void createbinary(int nid, int nsnp, int nchr, map *genmap, ped *dat, char *geno, chromosome *chrstat, char *filename)
{
	FILE *BIN = fopen(filename,"wb");
	
	printf("Storing\n\t%d individuals\n\t%d SNPs\n\t%d chromosomes\nas binary in %s..",nid,nsnp,nchr,filename);fflush(stdout);

	fwrite(&nid,sizeof(int),1,BIN);
	fwrite(&nsnp,sizeof(int),1,BIN);
	fwrite(&nchr,sizeof(int),1,BIN);
	fwrite(dat,sizeof(ped),nid,BIN);
	fwrite(genmap,sizeof(map),nsnp,BIN);
	fwrite(chrstat,sizeof(chromosome),nchr,BIN);
	fwrite(geno,sizeof(char),nid*nsnp,BIN);

	fclose(BIN);
	printf("done!\n");
}

// to be used in conjunction with convgenabel.R
void readgenabel(int *nid1, int *nsnp1, int *nchr1, map **genmap1, ped **dat1, char **geno1, chromosome **chrstat1, char *pedfile, char *mapfile)
{
	char *geno;
	int nid,nsnp,nchr;
	int i,j;
	float *phen;
	map *genmap;
	ped *dat;
	chromosome *chrstat;
	FILE *WHOLE;

	WHOLE = fopen(pedfile,"rb");
	fread(&nchr,sizeof(int),1,WHOLE);
	fread(&nid,sizeof(int),1,WHOLE);
	fread(&nsnp,sizeof(int),1,WHOLE);
	printf("%d IDs\n%dSNPs\n%d chromsomes\n",nid,nsnp,nchr);
	geno = (char *)malloc(sizeof(char) * nid * nsnp);
	phen = (float *)malloc(sizeof(float) * nid);
	dat = (ped *)malloc(sizeof(ped) * nid);
	fread(geno,sizeof(char),nid*nsnp,WHOLE);
	fread(phen,sizeof(float),nid,WHOLE);
	fclose(WHOLE);
	for(i = 0; i < nid; i++)
	{
		sprintf(dat[i].family,"%d",i);
		sprintf(dat[i].id,"1");
		sprintf(dat[i].paternal,"0");
		sprintf(dat[i].maternal,"0");
		//dat[i].id = i;
		dat[i].phen = phen[i];
	}
	free(phen);

	// map / chromosomes
	chrstat = (chromosome *) malloc(sizeof(chromosome) * nchr);
	genmap = (map *) malloc(sizeof(map) * nsnp);
	readmap(nsnp, genmap, &nchr, chrstat, mapfile);
	for(i = 0; i < nsnp; i++)
	{
		genmap[i].allelename[0] = 49;
		genmap[i].allelename[1] = 50;
	}

	printf("\nMap info:\n");
	for(i = 0; i < 10; i++)
	{
		printf("%s\t%s\t%f\t%lu\n",genmap[i].chr,genmap[i].snpname,genmap[i].gd,genmap[i].pd);
	}
	
	printf("\nGenotype sample:\n");
	for(i = 0; i < 20; i++)
	{
		for(j = 0; j < 20; j++)
		{
			printf("%d ",geno[(i + j*nid)]);
		}
		printf("\n");
	}

	printf("\nPhenotypes:\n");
	for(i = 0; i < 20; i++)
	{
		printf("%s:\t%f\n",dat[i].id, dat[i].phen);
	}
	printf("\n\n");
	
	*nid1 = nid;
	*nsnp1 = nsnp;
	*nchr1 = nchr;
	*chrstat1 = chrstat;
	*geno1 = geno;
	*genmap1 = genmap;
	*dat1 = dat;
}

// import binary file
void readbinary(int *nid, int *nsnp, int *nchr, map **genmap1, ped **dat1, char **geno1, chromosome **chrstat1, char *filename)
{
	int i,j;
	ped *dat;
	char *geno;
	map *genmap;
	chromosome *chrstat;

	FILE *BIN = fopen(filename,"rb");

	printf("\nReading in data...\n");
	fread(&*nid,sizeof(int),1,BIN);
	fread(&*nsnp,sizeof(int),1,BIN);
	fread(&*nchr,sizeof(int),1,BIN);
	printf("%d SNPs, %d individuals, %d chromosomes\n\n",*nsnp,*nid,*nchr); fflush(stdout);
	
	genmap = (map *)malloc(sizeof(map) * (*nsnp));
	dat = (ped *)malloc(sizeof(ped) * (*nid));
	chrstat = (chromosome *)malloc(sizeof(chromosome) * (*nchr));
	
	fread(&*dat,sizeof(ped),(*nid),BIN);
	fread(&*genmap,sizeof(map),(*nsnp),BIN);
	fread(&*chrstat,sizeof(chromosome),(*nchr),BIN);

	geno = malloc(sizeof(char) * (*nid) * (*nsnp));
	fread(geno,sizeof(char), (*nid) * (*nsnp), BIN);
	fclose(BIN);
	
	for(i = 0; i < 10; i++)
	{
		for(j = 0; j < 10; j++)
		{
			printf("%d ", geno[((*nid) * i + j)]);
		}
		printf("\n");
	}
	printf("\n");

	if(*nchr > 1){
		printf("%d chromosomes:\n",*nchr); fflush(stdout);
	} else {
		printf("%d chromosome:\n",*nchr); fflush(stdout);
	}
	for(i = 0; i < *nchr; i++)
	{
		printf("%s\t%d\t%d\n",chrstat[i].chrname,chrstat[i].chrsize,chrstat[i].chrstart);fflush(stdout);
	}
	printf("\nInfo:\n");
	for(i = 0; i < 10; i++)
	{
		printf("%s\t%s\t%s\t%s\t%d\t%f\n",dat[i].family,dat[i].id,dat[i].paternal,dat[i].maternal,(int)dat[i].sex, dat[i].phen);
	}
	printf("\nMap info:\n");
	
	for(i = 0; i < 10; i++)
	{
		printf("%s\t%s\t%f\t%lu\n",genmap[i].chr,genmap[i].snpname,genmap[i].gd,genmap[i].pd);
	}
	printf("\n");
	
	*chrstat1 = chrstat;
	*geno1 = geno;
	*genmap1 = genmap;
	*dat1 = dat;
}


// simulate data in correct format
// specify number if individuals, number of chromosomes, number of SNPs for each chromosome
void simulate(int *nid, int *nsnp, int *nchr, map **genmap1, ped **dat1, char **geno1, chromosome **chrstat1)
{
	int i, j, k, l, nmiss, ntot, binary, *buffer;
	float pmiss,p,gp[3];
	map *genmap;
	ped *dat;
	char *geno, userread[20], nucs[4] = {'A','T','C','G'}, s1, s2;
	chromosome *chrstat;
	*nid = 0;
	while(!(*nid))
	{
		printf("\nHow many individuals? --> ");
		scanf("%s",&*userread);
		*nid = atoi(userread);
		if(!(*nid))
		{
			printf("Invalid value!\n");
		}
	}
	*nchr = 0;
	while(!(*nchr) || *nchr > 200)
	{
		printf("\nHow many chromosomes? --> ");
		scanf("%s",&*userread);
		*nchr = atoi(userread);
		if(!(*nchr) || (*nchr) > 200)
		{
			printf("Invalid value!\n");
		}
	}
	chrstat = (chromosome *)malloc(sizeof(chromosome) * *nchr);
	*nsnp = 0;
	for(i = 0; i < *nchr; i++)
	{
		j = 0;
		while(!j)
		{
			printf("How many SNPs for chromosome %d --> ",i+1);
			scanf("%s",&*userread);
			j = atoi(userread);
			if(!j)
			{
				printf("Invalid value!\n");
			}
		}
		chrstat[i].chrid = i;
		sprintf(chrstat[i].chrname, "chr%d",i+1);
		chrstat[i].chrsize = j;
		(*nsnp)+=j;
	}
	chrstat[0].chrstart = 0;
	for(i = 1; i < *nchr; i++)
	{
		chrstat[i].chrstart = chrstat[i-1].chrstart + chrstat[i-1].chrsize;
	}

	printf("\n%d IDs, %d SNPs, %d chromosomes\n",*nid,*nsnp,*nchr);

	pmiss = -1;
	while(pmiss < 0 || pmiss > 0.5)
	{
		printf("\nProportion missing (0-0.5)? --> ");
		scanf("%s",&*userread);
		pmiss = (float)atof(userread);
		if(pmiss < 0 || pmiss > 0.5)
		{
			printf("Invalid value!\n");
		}
	}

	geno = (char *)malloc(sizeof(char) * *nid * *nsnp);
	genmap = (map *)malloc(sizeof(map) * *nsnp);
	dat = (ped *)malloc(sizeof(ped) * *nid);	
	printf("Memory allocated...\n");

	k = 0;
	buffer = (int *) malloc(sizeof(int) * *nid);
	for(i = 0; i < *nsnp; i++)
	{
		p = ((rand() * 0.9) / RAND_MAX + 0.05);
		gp[0] = p*p;
		gp[1] = 2*p*(1-p);
		gp[2] = (1-p)*(1-p);
		walkersalias(3, gp, *nid, buffer);
		for(j = 0; j < *nid; j++)
		{
			geno[j+k] = (char)buffer[j];
		}
		k += *nid;
//		printf("%f : ",p);
//		for(j = 0; j < 10; j++)
//		{
//			printf("%d ",buffer[j]);
//		}
//		printf("\n");
	}
	free(buffer);

	for(i = 0; i < 20; i++)
	{
		for(j = 0; j < 20; j++)
		{
			printf("%d ",geno[i* *nid + j]);
		}
		printf("\n");
	}

//	for(i = 0; i < (*nid * *nsnp); i++)
//	{
//		geno[i] = rand() % 3;
//	}

	ntot = *nid * *nsnp;
	nmiss = ntot * pmiss;
	if(nmiss)
	{
		printf("Generating random missing values...\n");
		k = 0;
		for(i = 0; i < nmiss; i++)
		{
			j = rand() % *nid;
			geno[k * *nid + j] = MISS;
			k++;
			if(k == *nsnp) k = 0;
		}
	}

	for(i = 0; i < *nid; i++)
	{
		sprintf(dat[i].family,"%d",i+1);
		sprintf(dat[i].id,"1");
		sprintf(dat[i].paternal,"0");
		sprintf(dat[i].maternal,"0");
		dat[i].sex = rand() % 2 + 1;
	}

	
	binary = -1;
	while(binary < 0 || binary > 1)
	{
		printf("\nContinuous (Recommended: 1) or binary (0) phenotype? --> ");
		scanf("%s",&*userread);
		binary = atoi(userread);
		if(binary < 0 || binary > 1)
		{
			printf("Invalid value!\n");
		}
	}

	if(binary)
	{
		rnorm(*nid,dat); // create random normal phenotype
	} else {
		for(i = 0; i < *nid; i++)
		{
			dat[i].phen = rand() % 2 + 1;
		}
	}

//	printf("\nInfo:\n");
//	for(i = 0; i < 10; i++)
//	{
//		printf("%d\t%d\t%d\t%d\t%d\t%f\n", dat[i].family, dat[i].id, dat[i].paternal, dat[i].maternal, dat[i].sex, dat[i].phen);
//	}
	
	k = 0;
	for(i = 0; i < *nchr; i++)
	{
		for(j = 0; j < chrstat[i].chrsize; j++)
		{
			for(l = 0; l < 30; l++)
			{
				genmap[k].chr[l] = chrstat[i].chrname[l];
				if(!genmap[k].chr[l]) break;
			}
			sprintf(genmap[k].snpname,"SNP%d",k+1);
			genmap[k].gd = 0;
			genmap[k].pd = 0;
			s1 = rand() % 4;
			s2 = s1;
			while(s2==s1)
			{
				s2 = rand() % 4;
			}
			genmap[k].allelename[0] = nucs[(int)s1];
			genmap[k].allelename[1] = nucs[(int)s2];
			k++;
		}
	}

//	for(i = 0; i < 20; i++)
//	{
//		for(j = 0; j < 20; j++)
//		{
//			printf("%d ", geno[((*nid) * i + j)]);
//		}
//		printf("\n");
//	}
//	printf("\n");

//	if(*nchr > 1){
//		printf("%d chromosomes:\n",*nchr); fflush(stdout);
//	} else {
//		printf("%d chromosome:\n",*nchr); fflush(stdout);
//	}
//	for(i = 0; i < *nchr; i++)
//	{
//		printf("%s\t%d\t%d\n",chrstat[i].chrname,chrstat[i].chrsize,chrstat[i].chrstart);fflush(stdout);
//	}
//	printf("\nMap info:\n");
	
//	for(i = 0; i < 10; i++)
//	{
//		printf("%s\t%s\t%lu\t%lu\n",genmap[i].chr,genmap[i].snpname,genmap[i].gd,genmap[i].pd);
//	}
//	printf("\n");

	*geno1 = geno;
	*dat1 = dat;
	*genmap1 = genmap;
	*chrstat1 = chrstat;
}

// simulate epistatic QTL pairs in phenotype
void simqtl(int nid, int nsnp, char *geno, ped *dat, map *genmap)
{
	int i, j, k, *X, nqtl, *qtl1, *qtl2, gcount[9], *buffer;
	float gmean[9], gp[3] = {0.25,0.5,0.25};
	char userread[20];
	int pat[9] = {-1,0,1,0,0,0,1,0,-1};
	X = (int *)malloc(sizeof(int) * nid);
	for(i = 0; i < 9; i++)
	{
		gcount[i] = 0;
		gmean[i] = 0;
	}

	printf("How many QTLs --> ");
	scanf("%s",&*userread);
	nqtl = atoi(userread);
	if(!nqtl)
	{
		printf("No QTLs simulated\n");
		return;
	}
	
	qtl1 = (int *)malloc(sizeof(int) * nqtl);
	qtl2 = (int *)malloc(sizeof(int) * nqtl);
	for(i = 0; i < nqtl; i++)
	{
		printf("QTL %d locus 1 [1:%d] --> ",i,nsnp);
		qtl1[i] = nsnp;
		while(qtl1[i] >= nsnp || qtl1[i] < 0)
		{
			scanf("%s",&*userread);
			qtl1[i] = atoi(userread) - 1;
			if(qtl1[i] >= nsnp || qtl1[i] < 0)
			{
				printf("Outside range!\nQTL %d locus 1 [1:%d] --> ",i,nsnp);
			}
		}
		printf("QTL %d locus 2 [1:%d] --> ",i,nsnp);
		qtl2[i] = nsnp;
		while(qtl2[i] >= nsnp || qtl2[i] < 0)
		{
			scanf("%s",&*userread);
			qtl2[i] = atoi(userread) - 1;
			if(qtl2[i] >= nsnp || qtl2[i] < 0)
			{
				printf("Outside range!\nQTL %d locus 2 [1:%d] --> ",i,nsnp);
			}
		}
		printf("\nSimulating AxA between %s and %s\n",genmap[qtl1[i]].snpname,genmap[qtl2[i]].snpname);
		
		buffer = (int *) malloc(sizeof(int) * nid);
		walkersalias(3, gp, nid, buffer);
		for(j = 0; j < nid; j++)
		{
			geno[j+nid*qtl1[i]] = (char)buffer[j];
		}
		walkersalias(3, gp, nid, buffer);
		for(j = 0; j < nid; j++)
		{
			geno[j+nid*qtl2[i]] = (char)buffer[j];
		}
		free(buffer);
		
		for(j = 0; j < nid; j++)
		{
			X[j] = 3*geno[(nid * qtl1[i] + j)] + geno[(nid * qtl2[i] + j)];
			for(k = 0; k < 9; k++)
			{
				if(X[j] == k)
				{
					(gcount[k])++;
					dat[j].phen += pat[k];
					(gmean[k]) += dat[j].phen;
					break;
				}
			}
		}
		printf("\nPairwise genotype frequencies:");
		for(j = 0; j < 9; j++)
		{
			if(j % 3 == 0) printf("\n");
			printf("%d\t",gcount[j]);
		}
		printf("\n\n");
		printf("G-P map:");
		for(j = 0; j < 9; j++)
		{
			if(j % 3 == 0) printf("\n");
			printf("%0.2f\t",gmean[j]);
		}
		printf("\n\n");
	}
	printf("Simulated %d QTL pair(s)\n\n",nqtl);
}

void bit_print(int a)
{
	int i, n = sizeof(int) * CHAR_BIT, mask = 1 << (n - 1);
	for(i = 1; i <= n; i++)
	{
		putchar(((a & mask) == 0) ? '0' : '1');
		a <<= 1;
		if(i % CHAR_BIT == 0 && i < n) putchar (' ');
	}
	putchar('\n');
}

void pack(int nid, int nsnp, int *npack, int *remain, int **gpack, char *geno)
{
	int i, j, k, l, n;
	int p, *gp;
	int size = 16;
	*npack = (int)floor((double)nid / size);
	*remain = nid % size;
//	printf("%d %d\n",*npack, *remain);
	n = *remain ? (*npack + 1) : *npack;
	gp = (int *)malloc(sizeof(int) * n * nsnp);
	
	k = 0;
	for(i = 0; i < nsnp; i++)
	{
		for(j = 0; j < *npack; j++)
		{
			p = geno[k++];
			for(l = 0; l < (size-1); l++)
			{
				p = (p << (32/size)) | geno[k++];
			}
			gp[(i*n + j)] = p;
			//bit_print(p);
		}
		if(*remain)
		{
			p = geno[k++];
			for(l = 0; l < (*remain-1); l++)
			{
				p = (p << (32/size)) | geno[k++];
			}
			gp[(i*n + j)] = p;
			//bit_print(p);
		}
	}
	*gpack = gp;
}

void unpack(int nid, int nsnp, int npack, int remain, int *gpack, char **geno1)
{
	int i, j, k, l, m, n;
	unsigned mask;
	char *geno;
	int size = 16;
	geno = (char *)malloc(sizeof(char) * nid * nsnp);
	l = 0;
	m = 0;
	for(i = 0; i < nsnp; i++)
	{
		for(j = 0; j < npack; j++)
		{
			n = (size-1) * 32/size;
			for(k = (size-1); k >= 0; k--)
			{
				//n = k * (32/size);
				mask = 3;
				mask <<= n;
				geno[l] = ((gpack[m] & mask) >> n);
				l++;
				n -= (32/size);
			}
			m++;
		}
		if(remain)
		{
			for(k = (remain-1); k >= 0; k--)
			{
				n = k * (32/size);
				mask = 3;
				mask <<= n;
				geno[l] = ((gpack[m] & mask) >> n);
				l++;
			}
			m++;
		}
	}
	*geno1 = geno;
}

// Create the epiGPU binary file
// Stores the following information in order:
// char "epiGPU"      - binary signature - to identify as epiGPU binary file
// int <nid>
// int <nsnp>
// int <npack>
// int <remain>
// int <nchr>
// ped array[nid]     - individual information (family, id, phenotype etc)
// map array[nsnp]    - map data
// chrstat array[nchr] - chromosome summary information
// int array[npack*nsnp] - packed genotype data

void createpackedbinary(int nid, int nsnp, int npack, int remain, int nchr, map *genmap, ped *dat, int *gpack, chromosome *chrstat, char *filename)
{
	int i;
	FILE *BIN = fopen(filename,"wb");
	char sig[10] = "epiGPU";

	printf("\nStoring\n- %d individuals\n- %d SNPs\n- %d chromosomes\nin epigpu binary format as file \"%s\"..",nid,nsnp,nchr,filename);fflush(stdout);

	fwrite(sig,sizeof(char),6,BIN);
	fwrite(&nid,sizeof(int),1,BIN);
	fwrite(&nsnp,sizeof(int),1,BIN);
	fwrite(&npack,sizeof(int),1,BIN);
	fwrite(&remain,sizeof(int),1,BIN);
	fwrite(&nchr,sizeof(int),1,BIN);
	fwrite(dat,sizeof(ped),nid,BIN);
	fwrite(genmap,sizeof(map),nsnp,BIN);
	fwrite(chrstat,sizeof(chromosome),nchr,BIN);
	i = remain ? (npack+1) : npack;
	
	fwrite(gpack,sizeof(int),i*nsnp,BIN);

	fclose(BIN);
	printf("done!\n");
}

void readpackedbinary(int *nid, int *nsnp, int *npack, int *remain, int *nchr, map **genmap1, ped **dat1, int **gpack1, chromosome **chrstat1, char *filename)
{
	int i,ii;
	ped *dat;
	int *gpack;
	map *genmap;
	char sigread[7], sig[7] = "epiGPU";
//	int n,mask,j;

	chromosome *chrstat;

	FILE *BIN = fopen(filename,"rb");
	if(BIN == NULL)
	{
		printf("No input file found\n");
		exit(1);
	}

	printf("\nReading in data...\n");

	// check that first 6 chars are "epiGPU"
	fread(sigread,sizeof(char),6,BIN);
	sigread[6] = '\0';
	i = strcmp(sigread,sig);
	if(i != 0)
	{
		printf("%d Invalid input file!\nPlease convert plink .ped and .map data files using the plinkbin programme.\n",i);
		exit(1);
	}
	fread(&*nid,sizeof(int),1,BIN);
	fread(&*nsnp,sizeof(int),1,BIN);
	fread(&*npack,sizeof(int),1,BIN);
	fread(&*remain,sizeof(int),1,BIN);
	fread(&*nchr,sizeof(int),1,BIN);
	printf("%d SNPs, %d individuals, %d chromosomes\n\n",*nsnp,*nid,*nchr); fflush(stdout);
	
	genmap = (map *)malloc(sizeof(map) * (*nsnp));
	dat = (ped *)malloc(sizeof(ped) * (*nid));
	chrstat = (chromosome *)malloc(sizeof(chromosome) * (*nchr));
	
	fread(&*dat,sizeof(ped),(*nid),BIN);
	fread(&*genmap,sizeof(map),(*nsnp),BIN);
	fread(&*chrstat,sizeof(chromosome),(*nchr),BIN);

	ii = *remain ? (*npack + 1) : *npack;

	gpack = (int *)malloc(sizeof(int) * ii * (*nsnp));
	fread(gpack,sizeof(int), ii * (*nsnp), BIN);
	fclose(BIN);

//	for(i = 0; i < 16; i++)
//	{
//		n = 30;
//		for(j = 15; j >= 0; j--)
//		{
//			mask = 3;
//			mask <<= n;
//			printf("%d ",(char)((gpack[(i*ii)] & mask) >> n));
//			n -= 2;
//		}
//		printf("\n");
//	}
//	printf("\n");

//	if(*nchr > 1){
//		printf("%d chromosomes:\n",*nchr); fflush(stdout);
//	} else {
//		printf("%d chromosome:\n",*nchr); fflush(stdout);
//	}
//	for(i = 0; i < *nchr; i++)
//	{
//		printf("%s\t%d\t%d\n",chrstat[i].chrname,chrstat[i].chrsize,chrstat[i].chrstart);fflush(stdout);
//	}
//	printf("\nInfo:\n");
//	for(i = 0; i < 10; i++)
//	{
//		printf("%d\t%d\t%d\t%d\t%d\t%f\n",dat[i].family,dat[i].id,dat[i].paternal,dat[i].maternal,dat[i].sex, dat[i].phen);
//	}
//	printf("\nMap info:\n");
	
//	for(i = 0; i < 10; i++)
//	{
//		printf("%s\t%s\t%f\t%lu\n",genmap[i].chr,genmap[i].snpname,genmap[i].gd,genmap[i].pd);
//	}
//	printf("\n");
	
	*chrstat1 = chrstat;
	*gpack1 = gpack;
	*genmap1 = genmap;
	*dat1 = dat;
}

void arguments(char **argv)
{
	printf("DATA MANAGEMENT MODE:\n\n" \
	"%s -D[ arguments ] [ filenames ... ]\n\n\n" \
	"Arguments:\n\n" \
	"r\tRead PLINK data\n" \
	"c\tClean data, remove low call rate SNPs and individuals\n" \
	"m\tImpute missing genotypes based on allele frequency\n" \
	"n\tReplace phenotype with random normally distributed sample\n" \
	"q\tSimulate AxA interactions at specific SNP pairs\n" \
	"s\tSimulate entire dataset\n" \
	"e\tExtract binary data to PLINK format\n\n" \
	"Example:\n" \
	"%s -Drm [ .ped file ] [ .map file ] [ epigpu file ]\n" \
	"%s -De [ .ped file ] [ .map file ] [ epigpu file ]\n" \
	"%s -Dsq [ epigpu file ]\n\n\n\n", argv[0], argv[0], argv[0], argv[0]);
}



// interpret arguments and perform data manipulation operations
int datamode(int argc, char **argv, int *nid, int *nsnp, int *npack, int *remain, int *nchr, ped **dat, map **genmap, chromosome **chrstat, char **binfile, char **pedfile, char **mapfile, int **genop, int **gpack, char **geno)
{
	int arg, flag, i, j;
	
	if(argc < 2) { printf("Wrong arguments\n"); arguments(argv); exit(1); }
	arg = (int)strlen(argv[1]); if(arg > 6) { printf("Wrong flags\n"); arguments(argv); exit(1); }

	flag = 0;
	for(j = 2; j < arg; j++)
	{
		if(argv[1][j] == 's') // simulate
		{
			if(argc != 3) { printf("Wrong arguments\n"); arguments(argv); exit(1); }
			*binfile = argv[2];
			printf("Simulating data and storing in binary file \"%s\"\n",*binfile);
			simulate(&*nid,&*nsnp,&*nchr,&*genmap,&*dat,&*geno,&*chrstat);
			flag = 1;
		}
	}

	for(j = 2; j < arg; j++)
	{
		if(argv[1][j] == 'e') // extract from binary
		{
			if(argc != 5) { printf("Wrong arguments\n"); arguments(argv); exit(1); }
			*pedfile = argv[2];
			*mapfile = argv[3];
			*binfile = argv[4];
			printf("Extracting data stored in \"%s\", and writing to \"%s\" and \"%s\"\n",*binfile,*pedfile,*mapfile);
			extractfrombinary(*binfile,*pedfile,*mapfile);
			return 0;
		}
	}

	for(j = 2; j < arg; j++)
	{
		if(argv[1][j] == 'r') // read plink data
		{
			flag = 1;
			if(argc != 5) { arguments(argv); exit(1); }
			*pedfile = argv[2];
			*mapfile = argv[3];
			*binfile = argv[4];
			dims(&*nid, &*nsnp, *pedfile);
			*dat = (ped *)malloc(sizeof(ped) * *nid);
			*genmap = (map *)malloc(sizeof(map) * *nsnp);
			*geno = (char *)malloc(sizeof(char) * *nid * *nsnp);
			*chrstat = (chromosome *)malloc(sizeof(chromosome) * 50);
			readped(*nid, *nsnp, *pedfile, *dat, *geno, *genmap);
			readmap(*nsnp, *genmap, &*nchr, *chrstat, *mapfile);
			for(j = 2; j < arg; j++)
			{
				if(argv[1][j] == 'c') // clean data
				{
					printf("Cleaning data\n");
					cleangeno(&*nid, &*nsnp, *nchr, *geno, *dat, *chrstat, *genmap);
				}
			}
			for(j = 1; j < arg; j++)
			{
				if(argv[1][j] == 'm') // impute
				{
					printf("Guessing missing genotypes from frequencies\n");
					guessmissing(*nid, *nsnp, *geno);
				}
			}
			for(j = 1; j < arg; j++)
			{
				if(argv[1][j] == 'n') // random normal phenotype
				{
					printf("Generating random continuous phenotype\n");
					rnorm(*nid, *dat);
				}
			}
		}
	}

	for(i = 2; i < arg; i++)
	{
		if(argv[1][i] == 'g') // GenABEL format - first run R --no-save < comvgenabel.R
		{
			if(argc != 5) { printf("Wrong arguments - pedfile, mapfile, binfile\n"); arguments(argv); exit(1); }
			*pedfile = argv[2];
			*mapfile = argv[3];
			*binfile = argv[4];
			printf("Reading GenABEL output and storing in binary file\n");
			readgenabel(&*nid, &*nsnp, &*nchr, &*genmap, &*dat, &*geno, &*chrstat, *pedfile, *mapfile);
			flag = 1;
			for(j = 2; j < arg; j++)
			{
				if(argv[1][j] == 'm') // impute
				{
					printf("Guessing missing genotypes from frequencies\n");
					guessmissing(*nid, *nsnp, *geno);
				}
			}
			for(j = 2; j < arg; j++)
			{
				if(argv[1][j] == 'n') // random normal phenotype
				{
					printf("Generating random continuous phenotype\n");
					rnorm(*nid, *dat);
				}
			}
		}
	}

	if(flag == 0)
	{
		printf("Neither simulating nor reading\n"); arguments(argv); exit(1);
	}

	for(j = 2; j < arg; j++)
	{
		if(argv[1][j] == 'q') // simulate QTLs
		{
			simqtl(*nid, *nsnp, *geno, *dat, *genmap);
		}
	}

	flag = 0;
	for(j = 2; j < arg; j++)
		if(argv[1][j] == 'u') // don't pack genotypes (this should never be needed)
			flag = 1;

	if(!flag)
	{
		pack(*nid, *nsnp, &*npack, &*remain, &*gpack, *geno);
		createpackedbinary(*nid, *nsnp, *npack, *remain, *nchr, *genmap, *dat, *gpack, *chrstat, *binfile);
		free(*gpack);
	} else {
		createbinary(*nid, *nsnp, *nchr, *genmap, *dat, *geno, *chrstat, *binfile);
	}

	free(*geno);
	free(*chrstat);
	free(*dat);
	free(*genmap);
	return flag;
}
