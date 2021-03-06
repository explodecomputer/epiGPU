epiGPU v2.0
===========

The original release is located here: [http://sourceforge.net/projects/epigpu/](http://sourceforge.net/projects/epigpu/). All the documentation and binaries for Windows, Mac and Linux can be found there.

This is work in progress for version 2. The following changes have been made:
- Allow missing phenotypic values
- Allow phenotypes to be specified in separate file
- Allow multiple platforms

Still to be done:
- Make major improvements on data handling side
- Deal exclusively with binary plink files (no conversion to epiGPU format required)
- Restructure command-line arguments to be similar to `plink` and `gcta`

Master branch is stable. To build, just run:

    git clone git://github.com/explodecomputer/epiGPU.git
    cd epiGPU
    make


## NOTE

Unfortunately the data handling side of this software needs to be substantially improved. 
- If you are testing the code, please try using at least 2000 SNPs per chromosome, and at least 500 individuals. Using a very small number of SNPs or individuals may lead to problems.
- There are problems converting larger datasets, i.e. above 4000 individuals and 600,000 SNPs. This will be dealt with as soon as I have time..!
