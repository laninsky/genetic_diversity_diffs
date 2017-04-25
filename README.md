# genetic_diversity_diffs v1.0.3
Tests for significant differences in genetic diversity between populations based on DNA sequence data

# How it works
This R script carries out a permutation procedure, reasmpling from the combined haplotype frequencies over all your populations, in order to test whether the observed differences in haplotype/nucleotide diversity between specific populations are greater than expected by chance.

# Things to note
1) This program "eats" arlequin files that have a haplotype block at the top of the file. Check out the example data if you are unsure on this format.

2) You need to make sure the library stringr is loaded in to your R environment

3) The program will look for and consolidate any identical haplotypes

4) It counts indels as a state for the calculation of nucleotide diversity (so make sure sequences of different lengths are 'padded' by missing data, not dashes).

5) If there are haplotypes that do not reciprocally match (i.e. a haplotype with missing data can match to more than one haplotype, and the haplotypes it matches to do not match to each other), the program will identify the haplotype with the largest amount of missing data. In most cases this will be the culprit, but a better idea is to make sure all of your haplotypes are unique before running this program.

6) Nucleotide diversity will differ slightly from arlequin if there is missing data, as genetic_diversity_diffs calculates the proportional difference between haplotypes only based on non-missing data, whereas arlequin averages the differences over the entire length of alignment to calculate proportional differences. 

7) Nucleotide diversity will also differ slightly between the programs if there are indels, because Arlequin only counts indels as part of the haplotype alignment length if they are variable within the population in question, and genetic_diversity_diffs counts them as part of the length if they are variable across all haplotype definitions (over all of your populations).

8) The program can handle spaces and tabs (or a mixture of these) as delimiters in your arlequin file. However, make sure that closing braces i.e. "}" do not occur on the same line as your data. This might give you a weird non-specific error.

If arlequin can accept the file, and genetic_diversity_diffs can't, then email alana dot alexander at ku dot edu. If you are having trouble with arlequin, make sure your hyphens haven't been made into em-dashes by word.

# How to run it
The easiest way for people less familiar with R, is to paste the entire function into R. You can then call the function by:

genetic_diversity_diffs(working_dir,file_name,n_iter,test_for_hap,test_for_nucl)

where:

working_dir == pathway to the folder with your arlequin file e.g. "C:/blahblahblah" 

file_name == the name of your arlequin file (in haplotype list format) e.g. "data.txt"

n_iter == the number of iterations you wish to perform e.g. 1000

test_for_hap == "Y"/"N" - do you wish to test for differences in haplotype diversity?

test_for_nucl == "Y"/"N" - do you wish to test for difference in nucleotide diversity?

# Example of input
genetic_diversity_diffs("C:/Users/Folder/","ATL_by_region_394.arp",1000,"Y","Y")

Demo arlequin file "ATL_by_region_394.arp" located in example folder

# Suggested citation
genetic_diversity_diffs was first published in:

Alexander, A., Steel, D., Hoekzema, K., Mesnick, S.L., Engelhaupt, D., Kerr, I., Payne, R. and Baker, C.S., 2016. What influences the worldwide genetic structure of sperm whales (Physeter macrocephalus)?. Molecular ecology. 25: 2754–2772

If you want to cite the specific version of the script you used, I suggest the following:

Alexander, A. 2017. genetic_diversity_diffs v1.0.3. Available from https://github.com/laninsky/genetic_diversity_diffs

This pipeline also wouldn't be possible without:

R: R Core Team. 2015. R: A language and environment for statistical computing. URL http://www.R-project.org/. R Foundation for Statistical Computing, Vienna, Austria. https://www.r-project.org/

Stringr:  Hadley Wickham (2012). stringr: Make it easier to work with strings..
  R package version 0.6.2. http://CRAN.R-project.org/package=stringr (for up-to-date citation information run citation("stringr" in R)

# Version history
1.0.4: Fixed a bug introduced in 1.0.3 that stopped the code from working for anything other than four populations

1.0.3: Vectorized some of the code to make it quicker. It now carries out nucleotide diversity calculations only counting differences at sites that do not involve ambiguous base calls, which might lead to slight differences in comparison with previous code.

1.0.2: gave some more commented-out notes on the arlequin file in the examples folder

1.0.1: made the output tables a little prettier for opening in excel

