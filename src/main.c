#include "epiopencl.h"
#define ARGS printf("%s",instructions)
// usage: 
// epigpu -D [data stuff...]
// epigpu -A [analysis stuff...]

// data
// -s
// -q
// -p - change to r
// -g
// -c
// -i - change to m
// -r - change to n
// -u

// analysis
// -p [n]	permutation
// -i [n]	chunksize
// -d [n]	device to use
// -t [f/i]	full test or interactions test
// -F [n]	threshold1
// -I [n]	threshold2

char *instructions =
	"\n\nepiGPU functions in two modes, data management (D) and analysis (A)\n\n" \
	"DATA MANAGEMENT MODE:\n\n" \
	"<epiGPU> -D[ arguments ] [ filenames ... ]\n\n\n" \
	"Arguments:\n\n" \
	"r\tRead PLINK data\n" \
	"c\tClean data, remove low call rate SNPs and individuals\n" \
	"m\tImpute missing genotypes based on allele frequency\n" \
	"n\tReplace phenotype with random normally distributed sample\n" \
	"q\tSimulate AxA interactions at specific SNP pairs\n" \
	"s\tSimulate entire dataset\n" \
	"e\tExtract binary data to PLINK format\n\n" \
	"Example:\n" \
	"<epigpu> -Drm [ .ped file ] [ .map file ] [ epigpu file ]\n" \
	"<epigpu> -De [ .ped file ] [ .map file ] [ epigpu file ]\n" \
	"<epigpu> -Dsq [ epigpu file ]\n\n\n\n" \
	"ANALYSIS MODE:\n\n" \
	"<epigpu> -A [ epigpu file] [ output file] [ Options ... ]\n\n" \
	"Options:\n\n" \
	"p [n]\tPermutation, default [n] = 0 (no permutation)\n" \
	"i [n]\tIteration size, default [n] = 128\n" \
	"d [n]\tDevice, default [n] = 0\n" \
	"t [f/i]\tFull test or full test with interaction test, default = f\n" \
	"F [n]\tF value threshold for full test, default [n] = 6.5\n" \
	"I [n]\tF value threshold for interaction test, default [n] = 10.5\n" \
	"s [0/1]\tNormal / Safe mode. Default = 0 (Normal mode)\n\n" \
	"Example:\n" \
	"<epigpu> -A [ epigpu file ] [ output file ]\n" \
	"<epigpu> -A [ epigpu file ] [ output file ] -i 512 -p 100 -F 7\n\n\n\n" \
	"For additional help please see README.pdf\n\n\n";


// Parse analysis mode argument flags
void analysismodecl(int argc, char **argv, char **binfile, char **outfile)
{
	int i, j, k;
	FILE *file;

	if(argc < 4 || (((int)argc & 1) != 0))
	{
		printf("Incorrect arguments\n");
		ARGS;
		exit(1);
	}
	parsefilename(argv[2]);
	*binfile = argv[2];
	*outfile = argv[3];
	silent = 0;
	if(argc >= 4)
	{
		for(i = 4; i < argc; i+=2)
		{
			if(argv[i][0] != '-')
			{
				printf("Incorrect optional arguments %c\n",argv[i][0]);
				ARGS;
				exit(1);
			}
			switch(argv[i][1])
			{
				case 's' :
					j = 0; k = 0;
					while(argv[i+1][j])
					{
						k += isdigit((int)argv[i+1][j++]) ? 0 : 1;
					}
					if(k != 0)
					{
						printf("Safe mode option: Some graphics cards may not support all optimisation techniques. To disable them set safe mode to 1.\nDefault = 0 (normal mode)\n\n"); exit(1);
					}
					USAFE = atoi(argv[i+1]);
					if(USAFE != 0 && USAFE != 1)
					{
						printf("Safe mode option: Some graphics cards may not support all optimisation techniques. To disable them set safe mode to 1.\nDefault = 0 (normal mode)\n\n"); exit(1);
					}
					break;
				case 'p' :
					j = 0; k = 0;
					while(argv[i+1][j])
					{
						k += isdigit((int)argv[i+1][j++]) ? 0 : 1;
					}
					if(k != 0)
					{
						printf("Please enter a numeric value for permutations\nDefault = 0 (no permutation)\n\n"); exit(1);
					}
					UPERM = atoi(argv[i+1]);
					break;
				case 'i' :
					j = 0; k = 0;
					while(argv[i+1][j])
					{
						k += isdigit((int)argv[i+1][j++]) ? 0 : 1;
					}
					if(k != 0)
					{
						printf("Please enter a numeric value for iteration size\nDefault = 128\n\n"); exit(1);
					}
					k = atoi(argv[i+1]);
					j = k % 16;
					if(k < 16 || j != 0)
					{
						printf("Iteration size must be a multiple of 16\nDefault = 128\n\n");
						exit(1);
					}
					UCHUNKSIZE = atoi(argv[i+1]);
					break;
				case 'd' :
					j = 0; k = 0;
					while(argv[i+1][j])
					{
						k += isdigit((int)argv[i+1][j++]) ? 0 : 1;
					}
					if(k != 0)
					{
						printf("Please enter a numeric value for device\nDefault = 0\nOnly change this if you have multiple graphics cards\n\n"); exit(1);
					}
					UDEVICE = atoi(argv[i+1]);
					break;
				case 't' :
					if(argv[i+1][0] != 'f' && argv[i+1][0] != 'i')
					{
						printf("Test selection:\nf\tFull 8df test\ni\tFull 8df test + 4df interaction test\n\n");
						exit(1);
					}
					UTEST = argv[i+1][0];
					break;
				case 'F' :
					UTHRESHOLD1 = (float)atof(argv[i+1]);
					if(UTHRESHOLD1 == 0)
					{
						printf("Please enter a numeric F-test threshold value.\nRecommended values 5.5 - 10\nDefault = 7\n\n");
						exit(1);
					}
					break;
				case 'I' :
					UTHRESHOLD2 = (float)atof(argv[i+1]);
					if(UTHRESHOLD2 == 0)
					{
						printf("Please enter a numeric F-test threshold value.\nRecommended values 5.5 - 10\nDefault = 7\n\n");
						exit(1);
					}
					if(UTEST == 'f')
					{
						printf("\n\nWarning: Interaction test threshold has been set, but interaction test isn't being performed\n\n\n");
					}
					break;
				case 'P' :
					if((file = fopen(argv[i+1], "rt")))
					{
						fclose(file);
						strcpy(UPHEN, argv[i+1]);
					} else {
						printf("Separate phenotype file does not exist\n");
						exit(1);
					}
					break;
				default :
					printf("Incorrect optional arguments %c\n",argv[i][0]);
					exit(1);
			}
		}
	}
}


void readphenotype(char *filename, int nid, ped *dat)
{
	FILE *file;
	int i;
	float temp;
	file = fopen(filename, "rt");
	printf("Reading in separate phenotype file\n");
	for(i = 0; i < nid; i++)
	{
		if(i == 0)
		{
			printf("First 3 values: "); fflush(stdout);
		}
		fscanf(file, "%f\n", &temp);
		dat[i].phen = temp;
		//printf("%f ", dat[i].phen);
		if(i < 3)
		{
			printf("%f ", dat[i].phen); fflush(stdout);
		}
		if(i == (nid - 3))
		{
			printf("\nLast 3 values: "); fflush(stdout);
		}
		if(i >= (nid - 3))
		{
			printf("%f ", dat[i].phen); fflush(stdout);
		}
	}
	fclose(file);
	printf("\nDone!\n\n");
}


int main(int argc, char **argv)
{
	char *binfile, *outfile;
	int nsnp, npack, remain, nid, nchr, tothits;
	ped *dat;
	map *genmap;
	int *genop;
	char *geno;
	chromosome *chrstat;
	int i;

	FILE *file;

	//plinkbin
	int flag;
	char *pedfile;
	char *mapfile;
	int *gpack;

	UPERM = 0;
	UDEVICE = 0;
	UCHUNKSIZE = 128;
	UTEST = 'f';
	UTHRESHOLD1 = 6.5;
	UTHRESHOLD2 = 10.5;
	if(argc < 3)
	{
		ARGS;
		exit(1);
	}
	if(argv[1][0] != '-')
	{
		printf("Incorrect arguments\n");
		ARGS;
		exit(1);
	}
	if(argv[1][1] == 'A') // analysis mode
	{
		flag = 2;
	} else if(argv[1][1] == 'D') // data mode
	{
		flag = 3;
	} else {
		ARGS;
		exit(1);
	}

	// DATA MODE
	if(flag == 3)
	{
		flag = datamode(argc, argv, &nid, &nsnp, &npack, &remain, &nchr, &dat, &genmap, &chrstat, &binfile, &pedfile, &mapfile, &genop, &gpack, &geno);
		return 0;
	}

	// ANALYSIS MODE
	// must have at least input and output files
	// must use 2 arguments per option e.g. number of options is always odd

	analysismodecl(argc, argv, &binfile, &outfile);
	
	printf("\n\nepiGPU\n\n" \
		"This programme will perform an exhaustive 2D scan for\n" \
		"epistasis, using the graphics card to improve speed\n\n" \
		"WARNING: Users may experience difficulty in performing\n" \
		"other tasks while this programme runs.\n\n" \
		"To interupt the scan press ctrl + c, and you will be able\n" \
		"to resume later by using the same output file.\n\n\n");

	readpackedbinary(&nid, &nsnp, &npack, &remain, &nchr, &genmap, &dat, &genop, &chrstat, binfile);

	if((file = fopen(UPHEN, "r")))
	{
		readphenotype(UPHEN, nid, dat);
	}

	flag = isbinary(nid,dat);
	if(flag)
	{
		printf("Error: Binary phenotype!\n" \
			"epiGPU can only use quantitative phenotypes for analysis.\n\n");
		exit(1);
	}

	// if permutation is set then set seed to permutation number and print rearranged phenotype to phen<n>.txt
	if(UPERM > 0)
	{
		printf("Permutation %d\n",UPERM);
		srand((size_t)UPERM);
		permute(nid, dat);
		for(i = 0; i < 10; i++)
		{
			printf("%f\n",dat[i].phen);
		}
	}

	if(USAFE)
	{
		printf("Safe mode\n");
	}
	if(UTEST == 'i')
	{
		printf("Performing 8 d.f. + 4 d.f. tests\n");
	} else {
		printf("Performing 8 d.f. test\n");
	}

	tothits = perform_scan(genop, chrstat, dat, genmap, nid, nsnp, npack, remain, nchr, outfile);

	return 0;
}

void permute(int nid, ped *dat)
{
	int randnum, i;
	ped temp;
	char filename[50];
	FILE *out;
	for(i = nid; i > 1; i--)
	{
		randnum = rand() % i;
		temp = dat[randnum];
		dat[randnum] = dat[i - 1];
		dat[i - 1] = temp;
	}
	sprintf(filename,"perm%d.txt",UPERM);
	out = fopen(filename,"w");
	fprintf(out,"Family ID Phenotype\n");
	for(i = 0; i < nid; i++)
	{
		fprintf(out,"%s %s %f\n",dat[i].family, dat[i].id, dat[i].phen);
	}
	fclose(out);
}

void parsefilename(char *argv1)
{
	int i, j, pos = 0;
	char temp = 1;
	i = 0;
	while(temp)
	{
		i++;
		temp = argv1[i];
		if(temp == '\\' || temp == '/') pos = i+1;
	}
	for(j = 0; j < (i-pos); j++)
	{
		plinkbin[j] = argv1[j+pos];
	}
	plinkbin[i] = '\0';
}

