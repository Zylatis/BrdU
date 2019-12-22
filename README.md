Package for curve fitting and analysis of 
BrdU + Ki67 experimental data

G. H. Gossel 2016
graeme.gossel@gmail.com

## Repo contents

### R files used as imports

- analysisHeaders.R
- formatOutput.R
- getBootSettings.R
- getData.R
- getFunctions.R
- getModelConfigs.R
- getODEs.R
- getODEsolver.R
- getSettings.R
- singleSetCI.R

### R files run directly by user
- runBestFits.R
- runAnalysis.R
- runCombineBootData.R
- runBootstrap.R
- runExamine.R

### Mathematica files used to prepare analytic expressions

- makeCfns.nb
- MakeRFunctions.nb

### Raw data files
- 29AprilRegate.xlsx

### Formatted data files

4Tem and 4Tcm versions of:
- 4TxxBRDU.xls
- 4TxxAllKi67Pos.xls
- 4TxxMeanKi67Pos.xls

## Requirements
The following R packages are required

- iterators
- parallel
- doParallel
- foreach
- Rcpp
- xlsx*
- odeintr**
- gdata
- stats
- nleqslv
- deSolve
- fBasics
- GenSA
- ggplot2
- reshape
- gridExtra


* This can be difficult to get working because of rJava issues. See "rJavaNotes.txt" for more information.

** This is an in-development package. To resolve several issues, the most up to date (development)
version from Github was used and can be found here: https://github.com/thk686/odeintr. In order to
obtain this version directly from git via commandline, additional packages will be required (httr, git2r, curl, openssl)

## Instructions

1. Modify workDir.R to reflect correct local package directory.

2. Modify getModelConfigs.R to define kBox/bBox min/max and desired source terms to loop over.
    For now one must use only configurations for which the ODE expressions already exist in Functions/
    and CFunctions/. Further information on how to use the Mathematica files to generate further expressions
    will be forthcoming.

3. Run runBestFits.R using Rscript which takes command line arguments according to:

	Rscript runBestFits.R <cell> <source switch> <model> <ncores>

    where the arguments are chosen from sets according to:
    
    cell - {4Tem, 4Tcm}
    source.switch -{immediate, delayStep}
    model - {kinHet, kinHetExtended, kinHetExtended2,tempHet}
    ncores - integer > 0

    Example:	
	Rscript runBestFits.R 4Tem delayStep kinHet 64
 
    Note: There is no check to ensure the number of cores is reasonable,
    it is currently up to the user to ensure this for their machine.

4. Modify runAnalysis.R according to which source terms the user wants included in analysis by
    re-defining 'target.source.list'
    
5. Run Rscript runAnalysis.R <model> <cell> <source switch>. This will output the summary files for the best fits.
    
    Example:
	Rscript runAnalysis.R kinHet 4Tem delayStep 
    
6. Run 'runBootstrap.R' with command line arguments for cell type, source switch function, sigma value, model, and number of cores
    
    Example:
	Rscript runBootStrap.R 4Tem delayStep 0.4 kinHet 64

    The above command will run the bootstrap code for the case where the source term is 40% of
    the memory pool per week. This presupposes the bestfits have been run for this source. It will
    automatically do the bootstraps on the 'best best-fit' box configuration available

7. Run 'runCombinedBootData.R' to collect all the bootstrap outputs into a useful form (CIs). This takes some time. 
    Only takes three command line arguments, <model>, <cell>, and <ncores>
