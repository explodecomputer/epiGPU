#include "readplink.h"


char **get_space(int m, int n)
{
        char i, *p, **a;
        p = malloc(m * n * sizeof(int));
        a = malloc(m * sizeof(int *));
        for(i = 0; i < m; i++) a[i] = p + (i * n);
        return a;
}

void release_space(char** X)
{
        int * p;
        p = (int *) X[0];
        free(p);
        free(X);
}

void general_dims(int *rowCount, int *colCount, char *filename)
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
	printf("done!\n");fflush(stdout);

	fclose(ifp);
}

void readplinkfam(int nid, ped *dat, char *famfile)
{
	int i,j;
	FILE *FAM;

	FAM = fopen(famfile,"r");
	printf("Reading in %s",famfile);fflush(stdout);
	for(i = 0; i < nid; i++)
	{
		if(i % ((int)nid/10) == 1) {printf(".");fflush(stdout);}
		(void)fscanf(FAM,"%s",dat[i].family);
		(void)fscanf(FAM,"%s",dat[i].id);
		(void)fscanf(FAM,"%s",dat[i].paternal);
		(void)fscanf(FAM,"%s",dat[i].maternal);
		(void)fscanf(FAM,"%d",&j);
		(void)fscanf(FAM,"%f",&dat[i].phen);
		dat[i].sex = j;
	}
	fclose(FAM);
	printf("done!\n");
}

void readplinkbim(int nsnp, map *genmap, char *bimfile)
{
	int i;
	FILE *MAP;
	char sep;
	MAP = fopen(bimfile,"r");
	printf("Reading in %s",bimfile);fflush(stdout);
	for(;;)
	{
		fscanf(MAP,"%c",&sep);
		if(sep == ' ' || sep == '\t' || sep == EOF) break;
	}
	fclose(MAP);

	if(sep == EOF)
	{
		printf("Problem with data, must be space or tab delimited\n\n");
		exit(1);
	}

	MAP = fopen(bimfile,"r");
	
	for(i = 0; i < nsnp; i++)
	{
		if((i % (nsnp / 10)) == 0) { printf("."); fflush(stdout); }
		(void)fscanf(MAP,"%s",genmap[i].chr);
		(void)fscanf(MAP,"%s",genmap[i].snpname);
		(void)fscanf(MAP,"%lf",&genmap[i].gd);
		(void)fscanf(MAP,"%li",&genmap[i].pd);
		(void)fscanf(MAP,"%s",&genmap[i].allelename[0]);
		(void)fscanf(MAP,"%s",&genmap[i].allelename[1]);
	}
	fclose(MAP);
	printf("done!\n");fflush(stdout);
}

void chromosomedetails(int nsnp, map *genmap, int *nchr, chromosome *chrstat)
{
	int i, j, k;

	*nchr = 0;
	for(i = 0; i < nsnp; i++)
	{
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

		if(k == 0) // found new chromosome
		{
			chrstat[(*nchr)].chrsize = 1;
			chrstat[(*nchr)].chrid = (*nchr);
			strcpy(chrstat[(*nchr)].chrname, genmap[i].chr);
			(*nchr)++;
		}
	}
	
	chrstat[0].chrstart = 0;
	for(i = 1; i < *nchr; i++)
	{
		chrstat[i].chrstart = chrstat[i-1].chrstart + chrstat[i-1].chrsize;
	}

	if(*nchr == 1) {
		printf("%d chromosome:\n",*nchr); fflush(stdout);
	} else {
		printf("%d chromosomes:\n",*nchr); fflush(stdout);
	}

	printf("Name\tSNPs\tStart\n");
	for(i = 0; i < *nchr; i++)
	{
		printf("%s\t%d\t%d\n",chrstat[i].chrname,chrstat[i].chrsize,chrstat[i].chrstart);fflush(stdout);
	}
	printf("\n");
}


// read in SNP major mode binary plink file (.bed)
// store as char
// convert to int
// genotypes are reverse order

void readplinkbedchar(int nid, int nsnp, int npack, int remain, char *bped, char *bedfile)
{
	char tmp;
	FILE *BED = fopen(bedfile, "rb");
	
	fread(&tmp, sizeof(char), 1, BED);
	if(tmp != 108) {
		printf("%s is not a binary plink file", bedfile);
		exit(1);
	}
	
	fread(&tmp, sizeof(char), 1, BED);
	if(tmp != 27) {
		printf("%s is not a binary plink file", bedfile);
		exit(1);
	}
		
	fread(&tmp, sizeof(char), 1, BED);
	if(tmp != 1) {
		printf("%s must be in SNP-major order", bedfile);
		exit(1);
	}
	
	fread(bped, sizeof(char), npack * nsnp, BED);
	
	fclose(BED);
}
	

void unpackint(int genop, int *geno)
{
	int k, l;
	unsigned int mask;
	
	l = 0;
	mask = 3;
	for(k = 0; k < 16; k++)
	{
		geno[l] = ((genop & mask) >> (k*2));
		l++;
		mask <<= 2;
	}
}

void readplinkbedint(int nid, int nsnp, int *genop, char *bedfile)
{
	int i, j, l;
	int npack_int, remain_int, npack_char, remain_char, npack_char2;
	unsigned char tmp, tmp4[4];
	FILE *BED = fopen(bedfile, "rb");

	fread(&tmp, sizeof(char), 1, BED);
	if(tmp != 108) {
		printf("%s is not a binary plink file", bedfile);
		exit(1);
	}

	fread(&tmp, sizeof(char), 1, BED);
	if(tmp != 27) {
		printf("%s is not a binary plink file", bedfile);
		exit(1);
	}

	fread(&tmp, sizeof(char), 1, BED);
	if(tmp != 1) {
		printf("%s must be in SNP-major order", bedfile);
		exit(1);
	}


	npack_int = (int) floor((float)nid / 16);
	remain_int = nid % 16;
	npack_char = (int) floor((float)(nid - npack_int * 16) / 4);
	remain_char = (nid - npack_int * 16) % 4;
	npack_char2 = remain_char ? npack_char + 1 : npack_char;

	printf("npack_int = %d\nremain_int = %d\nnpack_char = %d\nremain_char=%d\nnpack_char2=%d\n", npack_int, remain_int, npack_char, remain_char, npack_char2);
	
//  	printf("%d\tindividuals\n", nid);
//  	printf("%d\tnpack_int\n", npack_int);
//  	printf("%d\tremain_int\n", remain_int);
//  	printf("%d\tnpack_char\n", npack_char);
//  	printf("%d\tremain_char\n", remain_char);
//  	printf("%d\tnpack_char2\n", npack_char2);
//  	printf("%d\tindividuals\n", npack_int * 16 + npack_char * 4 + remain_char);

	printf("Reading in %s", bedfile);
	l = 0;
	for(i = 0; i < nsnp; i++)
	{
		if((i % (nsnp / 10)) == 0) { printf("."); fflush(stdout); }

		genop[l] = 0;
		for(j = 0; j < npack_int; j++)
		{
			fread(tmp4, sizeof(char), 4, BED);
			genop[l] |= ( tmp4[3] << 24 );
			genop[l] |= ( tmp4[2] << 16 );
			genop[l] |= ( tmp4[1] <<  8 );
			genop[l] |= ( tmp4[0]       );
			l++;
		}

		for(j = npack_char2; j > 0; j--)
		{
			fread(&tmp, sizeof(char), 1, BED);
			genop[l] |= ( tmp  << (j - 1)*8 );
		}
		if(npack_char2 > 0)
		{
			l++;
		}
	}
	printf("done!\n"); fflush(stdout);

}

// This function takes the readplinkbedint output and checks that it looks about right
void bplink2geno(int nid, int nsnp, int *genop)
{
	int i, j, k, l, m, n, npack, remain;
	int size=16;
	char *geno;

	geno = (char *)malloc(sizeof(char) * nid * nsnp);

	npack = (int) floor((float)nid / 16);	// print first few SNPs and individuals
	remain = nid % 16;
	l = 0; // number of individuals
	m = 0; // number of packs

	printf("npack = %d\nremain = %d\n\n", npack, remain);

	// SNP major order
	for(i = 0; i < nsnp; i++)
	{
		// there are npack packets of individuals, each with 16 individuals
		for(j = 0; j < npack; j++)
		{
			n = 3; // start at the end
			for(k = 0; k < size; k++)
			{
				geno[l] = ((genop[m] & n) >> (k*2));
				l++;
				n <<= 2;
			}
			m++;
		}
		if(remain)
		{
			n = 3;
			for(k = 0; k < remain; k++)
			{
				geno[l] = ((genop[m] & n) >> (k*2));
				l++;
				n <<= 2;
			}
			m++;
		}
	}

	for(i = 0; i < 20; i++) // individual
	{
		for(j = 0; j < 20; j++) // snp
		{
			printf("%d ", geno[j*nid+i]);
		}
		printf("\n");
	}

	printf("\n");

	// print first few SNPs and individuals
	for(i = nid-20; i < nid; i++) // individual
	{
		for(j = nsnp-20; j < nsnp; j++) // snp
		{
			printf("%d ", geno[j*nid+i]);
		}
		printf("\n");
	}

	free(geno);
}


void unpack2(int nid, int nsnp, int *gpack)
{
	int i, j, k, l, m, n, npack, remain;
	unsigned mask;
	char *geno;
	int size = 16;
	geno = (char *)malloc(sizeof(char) * nid * nsnp);

	npack = (int) floor((float)nid / 16);	// print first few SNPs and individuals
	remain = nid % 16;

	l = 0;
	m = 0;


	printf("recoding... again\n");
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
	for(i = 0; i < 20; i++) // individual
	{
		for(j = 0; j < 20; j++) // snp
		{
			printf("%d ", geno[j*nid+i]);
		}
		printf("\n");
	}

	printf("\n");

	// print first few SNPs and individuals
	for(i = nid-20; i < nid; i++) // individual
	{
		for(j = nsnp-20; j < nsnp; j++) // snp
		{
			printf("%d ", geno[j*nid+i]);
		}
		printf("\n");
	}
	printf("%d %d\n", geno[1000*nid+500],geno[10000*nid+200]);

	free(geno);
}


void recodeplinkbedint(int nid, int nsnp, int *genop)
{
	int i, j, npack_int, pack, mask, recode, snp1, snp2;
	
	npack_int = (int) ceil((float)nid / 16);

	printf("\nRecoding");
	printf("npack_int = %d\n", npack_int);

	for(i = 0; i < npack_int*nsnp; i++)
	{
		pack = genop[i];
		mask = 3;
		recode = 0;
		for(j = 0; j < 16; j++)
		{
			snp1 = ((pack & mask) >> (j*2));
//			if(i < 100) printf("%d ", snp1);
			mask <<= 2;
			if(snp1 == 1)
			{
				snp2 = 3;
			} else if(snp1 == 2) {
				snp2 = 1;
			} else if(snp1 == 3) {
				snp2 = 0;
			} else {
				snp2 = 2;
			}
			recode <<= 2;
			recode |= snp2;
//			if(i < 100) printf("%d ", recode);
		}
		if(i < 100)
		{
//			printf("%d %d %d; ", genop[i], recode, pack);
		}
//		if(i < 100) printf("\n");
		genop[i] = recode;
	}
	printf("done!\n\n");
}


void dec2bin(unsigned int decimal, char *binary)
{
  int  k = 0, n = 0;
  int  neg_flag = 0;
  int  remain;
  int  old_decimal;  // for test
  char temp[80];
  // take care of negative input
  if (decimal < 0)
  {      
    decimal = -decimal;
    neg_flag = 1;
  }
  do 
  {
    old_decimal = decimal;   // for test
    remain    = decimal % 2;
    // whittle down the decimal number
    decimal   = decimal / 2;
    // this is a test to show the action
    // printf("%d/2 = %d  remainder = %d\n", old_decimal, decimal, remain);
    // converts digit 0 or 1 to character '0' or '1'
    temp[k++] = remain + '0';
  } while (decimal > 0);
  if (neg_flag)
    temp[k++] = '-';       // add - sign
  else
    temp[k++] = ' ';       // space
  // reverse the spelling
  while (k >= 0)
    binary[n++] = temp[--k];
  binary[n-1] = 0;         // end with NULL
}

char file_exists(const char * filename)
{
	FILE *file = fopen(filename, "r");
	if (file)
	{
		fclose(file);
		return 1;
	}
	printf("%s does not exist!\n\n", filename);
	return 0;
}


void memory_calc(long bytes, char *size)
{
	if(bytes / 1024 < 1024)
	{
		sprintf(size, "%0.1f Kb", (double) bytes / 1024);
	} else if(bytes / 1024 / 1024 < 1024) {
		sprintf(size, "%0.1f Mb", (double) bytes / 1024 / 1024);
	} else {
		sprintf(size, "%0.1f Gb", (double) bytes / 1024 / 1024 / 1024);
	}
}


void readbinaryplink(int *nid, int *nsnp, int *npack, int *remain, int *nchr, map **genmap1, ped **dat1, int **gpack1, chromosome **chrstat1, char **argv)
{
	int tmp;
	char bimfile[1000], famfile[1000], bedfile[1000];
	int *genop;
	chromosome *chrstat;
	map *genmap;
	ped *dat;
	char size[1000];

//	strcpy(bedfile, filename);
//	strcat(bedfile, ".bed");
	strcpy(bedfile, argv[2]);
	if(!file_exists(bedfile)) exit(1);
	
	strcpy(bimfile, argv[3]);
//	strcat(bimfile, ".bim");
	if(!file_exists(bimfile)) exit(1);

	strcpy(famfile, argv[4]);
	if(!file_exists(famfile)) exit(1);

	general_dims(&*nid, &tmp, famfile);
	general_dims(&*nsnp, &tmp, bimfile);

	printf("\n%d individuals\n%d SNPs\n\n", *nid, *nsnp);

	*npack = (int) ceil((float) *nid / 16);
	*remain = *nid % 16;
	memory_calc((long)(*npack * *nsnp * 4), size);
	printf("Allocating CPU memory (%s)..........", size);
	genop = malloc(sizeof(int) * *npack * *nsnp);
	printf("Done!\n\n");
	
	dat = malloc(sizeof(ped) * *nid);
	readplinkfam(*nid, dat, famfile);
	
	genmap = malloc(sizeof(map) * *nsnp);
	readplinkbim(*nsnp, genmap, bimfile);

	
	// read binary format genotypes

	readplinkbedint(*nid, *nsnp, genop, bedfile);
	bplink2geno(*nid, *nsnp, genop);
	recodeplinkbedint(*nid, *nsnp, genop);
	unpack2(*nid, *nsnp, genop);

	chrstat = malloc(sizeof(chromosome) * 1000);
	chromosomedetails(*nsnp, genmap, &*nchr, chrstat);

	int i;
	for(i = 0; i < 20; i++)
	{
		printf("%s\t%s\t%s\t%s\t%d\t%f\n", dat[i].family, dat[i].id, dat[i].paternal, dat[i].maternal, dat[i].sex, dat[i].phen);
	}
	for(i = *nid-1; i > *nid-20; i--)
	{
		printf("%s\t%s\t%s\t%s\t%d\t%f\n", dat[i].family, dat[i].id, dat[i].paternal, dat[i].maternal, dat[i].sex, dat[i].phen);
	}


	for(i = 0; i < 20; i++)
	{
		printf("%s\t%s\t%lf\t%li\t%d\t%d\n", genmap[i].chr, genmap[i].snpname, genmap[i].gd, genmap[i].pd, genmap[i].allelename[0], genmap[i].allelename[1]);
	}
	for(i = *nsnp-1; i > *nsnp-20; i--)
	{
		printf("%s\t%s\t%lf\t%li\t%d\t%d\n", genmap[i].chr, genmap[i].snpname, genmap[i].gd, genmap[i].pd, genmap[i].allelename[0], genmap[i].allelename[1]);
	}


// 	for(i = 0; i < 10; i++)
// 	{
// 		printf("%s\t%s\t%s\t%s\t%d\t%f\n",
// 			dat[i].family, 
// 			dat[i].id, 
// 			dat[i].paternal, 
// 			dat[i].maternal, 
// 			dat[i].sex, 
// 			dat[i].phen);
// 	}
// 	printf("\n\n\n");
// 	for(i = 0; i < 10; i++)
// 	{
// 		printf("%s\t%s\t%f\t%lu\t%c\t%c\n",
// 			genmap[i].chr, 
// 			genmap[i].snpname, 
// 			genmap[i].gd, 
// 			genmap[i].pd, 
// 			genmap[i].allelename[0], 
// 			genmap[i].allelename[1]);
// 	}
	
	*chrstat1 = chrstat;
	*gpack1 = genop;
	*genmap1 = genmap;
	*dat1 = dat;
}

/*
int main(int argc, char **argv)
{
	int nid, nsnp, nchr, npack, remain;
	char bimfile[1000], famfile[1000], bedfile[1000];
	int *genop;
	chromosome *chrstat;
	map *genmap;
	ped *dat;


	readbinaryplink(&nid, &nsnp, &npack, &remain, &nchr, &genmap, &dat, &genop, &chrstat, argv[1]);


	free(genop);
	free(dat);
	free(genmap);
	free(chrstat);

	return 0;
}
*/

