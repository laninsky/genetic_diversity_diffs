genetic_diversity_diffs <- function(working_dir,file_name,n_iter,test_for_hap,test_for_nucl) {

# The stringr library is required
library(stringr)

# Throwing out error messages if any of the inputs are missing from the command line
x <- 0
error_one <- 0
error_two <- 0
error_three <- 0
error_four <- 0
error_five <- 0

killswitch <- "no"

if(missing(working_dir)) {
x <- 1
error_one <- 1
killswitch <- "yes"
}

if(missing(file_name)) {
x <- 1
error_two <- 2
killswitch <- "yes"
}

if(missing(n_iter)) {
x <- 1
error_three <- 3
killswitch <- "yes"
}

if(missing(test_for_hap)) {
x <- 1
error_four <- 4
killswitch <- "yes"
}

if(missing(test_for_nucl)) {
x <- 1
error_five <- 5
killswitch <- "yes"
}

if(x==1) {
cat("Call the program by genetic_diversity_diffs(working_dir,file_name,n_iter,test_for_hap,test_for_nucl), where:\nworking_dir == pathway to the folder with your arlequin file e.g. \"C:/blahblahblah\" \nfile_name == the name of your arlequin file (in haplotype list format) e.g. \"data.txt\"\nn_iter == the number of iterations you wish to perform e.g. 1000\ntest_for_hap == \"Y\"/\"N\" - do you wish to test for differences in haplotype diversity?\ntest_for_nucl == \"Y\"/\"N\" - do you wish to test for difference in nucleotide diversity?\n\nThe program can handle spaces and tabs (or a mixture of these) as delimiters in your arlequin file\nHowever, make sure that closing braces i.e. \"}\" do not occur on the same line as your data. This might give you a weird non-specific error.\n\nExample of input:\ngenetic_diversity_diffs(\"C:/Users/Folder/\",\"ATL_by_region_394.arp\",1000,\"Y\",\"Y\")\n\nSpecific errors/missing inputs:\n")
}
if(error_one==1) {
cat("Sorry, I am missing a working directory pathway\nworking_dir == pathway to the folder with your arlequin file e.g. \"C:/blahblahblah\" \n\n")
}
if(error_two==2) {
cat("Sorry, I am missing a filename for your input arlequin file with haplotype list\nfile_name == the name of your arlequin file (in haplotype list format) e.g. \"data.txt\"\n\n")
}
if(error_three==3) {
cat("How many iterations would you like to complete?\nn_iter == the number of iterations you wish to perform e.g. 1000\n\n")
}
if(error_four==4) {
cat("Would you like me to test for differences in haplotype diversity?\ntest_for_hap == \"Y\"/\"N\"\n\n")
}
if(error_five==5) {
cat("Would you like me to test for differences in nucleotide diversity?\ntest_for_nucl == \"Y\"/\"N\"\n\n")
}

# Error message tested and functioning
if (killswitch=="yes") {
stop("Please fill in the missing info in your function call")
}

#Checking status of working directory
print(noquote("STEP ONE: Loading in all the variables"))
print(noquote(""))
print(noquote("An error message after this indicates your working directory is not valid"))
flush.console()
setwd(working_dir)
print(noquote("Not to worry, your working directory IS valid! I've successfully set the working directory"))
print(noquote(""))
flush.console()

#Checking status of arlequin file
print(noquote("An error message after this indicates your file is not located in the directory you listed"))
flush.console()
input <- readLines(file_name)
print(noquote("Not to worry, your file IS located in the directory! I'm pulling it into my memory to extract the parts we are interested in"))
print(noquote(""))
flush.console()

inputlength <- length(input)
matrixinput <- matrix(NA, nrow=inputlength, ncol=2)
matrixinput[,1] <- input
rm(input)

i <- 0
j <- 1

missingdata <- NULL

for (j in 1:inputlength) {
if (grepl("MissingData[[:blank:]]*=", matrixinput[j,1], fixed=FALSE)==TRUE) {
if (grepl("'", matrixinput[j,1])==TRUE) {
missingdata <- (unlist(strsplit(matrixinput[j,1],"'"))[2])
} else {
missingdata <- (unlist(strsplit(matrixinput[j,1],'"'))[2])
}
break
}
}

j <- 1

# Carrying out some manipulations on the input arlequin file to extract the bits we are interested in
while (j <= inputlength)  {
if (grepl("HaplList[[:blank:]]*=[[:blank:]]*\\{", matrixinput[j,1], fixed=FALSE)==TRUE) {
i <- i + 1
}
if  (grepl("SampleName[[:blank:]]*=", matrixinput[j,1], fixed=FALSE)==TRUE) {
i <- i + 1
}
if  (grepl("\\}", matrixinput[j,1])==TRUE) {
i <- i + 1
}
matrixinput[j,2] <- i
j <- j+1
}

i <- NULL
j <- 1
protohaplist <- NULL
matchbad <- c("",NA)

while (j <= inputlength)  {
if (grepl("HaplList[[:blank:]]*=[[:blank:]]*\\{", matrixinput[j,1],ignore.case=TRUE)==TRUE) {
   k <- j + 1
   while (matrixinput[k,2]==matrixinput[j,2]) {
   if(!(matrixinput[k,1] %in% matchbad)){
   protohaplist <- rbind(protohaplist, matrixinput[k,])
   }
   k <- k + 1
   }
   }
j <- j+1
}

# Error message tested and functioning
if(is.null(protohaplist)) {
stop("\n\nYour arlequin file does not seem to have a HaplList embedded in the file under [[HaplotypeDefinitions]].\nPlease reformat your arlequin file, so the haplotypes are defined in a separate [[HaplotypesDefinitions]] block\n")
} else {
print(noquote("I have found the haplotype definition block in your arlequin file"))
flush.console()
}

# Checking how many populations there are in the file
numpop <- length(grep("SampleName[[:blank:]]*=", matrixinput))

protohaplistlength <- dim(protohaplist)[1]
haplist <- matrix(NA, nrow=protohaplistlength+1, ncol=2+numpop)
haplist[1,1] <- c("hap_name")
haplist[1,2] <- c("hap_sequence")

for (j in 1:protohaplistlength) {
haplist[j+1,1] <- unlist(strsplit((str_trim(protohaplist[j,1],side="left")),"[[:blank:]]+",fixed=FALSE))[1]
haplist[j+1,2] <- unlist(strsplit((str_trim(protohaplist[j,1],side="left")),"[[:blank:]]+",fixed=FALSE))[2]
}

haplistlength <- dim(haplist)[1]

checkingforduphapnames <- as.data.frame(table(haplist[2:haplistlength,1]))
checkingforduphapnames <- checkingforduphapnames[order(checkingforduphapnames[,2]),]
checkingforduphapnameslen <- dim(checkingforduphapnames)[1]
if(checkingforduphapnames[checkingforduphapnameslen,2]>1) {
print(noquote(""))
print(noquote("Multiple haplotype definitions with the same haplotype name are present"))
print(noquote("Please check the following haplotype name, and make sure all haplotype definitions"))
print(noquote("have unique haplotype names:"))
print(droplevels(checkingforduphapnames[checkingforduphapnameslen,1]))
flush.console()
stop("Please rename the multiple haplotype definitions that have the above name so that they have unique names")
}

minhaplength <- min(nchar(haplist[2:haplistlength,2]))
maxhaplength <- max(nchar(haplist[2:haplistlength,2]))

#Error messages tested and funcitons
if(minhaplength==maxhaplength) {
   if(all(is.na(haplist[2:haplistlength,2]))) {
stop("\n\nNo sequence given for haplotypes. Please make sure you have sequence defined for each of your haplotypes in your haplist\n")
}
} else  {
stop("\n\nYour arlequin file does not have equal length haplotypes. Please check your alignment and try again\n")
}

print(noquote("Your haplotypes appear to be of equal lengths"))
print(noquote(""))
flush.console()

# Filling in the number of haplotypes by population in our file 'haplist'.  Error message tested and functioning
print(noquote("An error message following this suggests there is something wrong with your haplist block"))
print(noquote("If arlequin opens your file correctly, but this program doesn't, please contact alana.alexander@ku.edu"))
flush.console()

popno <- 2
j <- 1
matcharray <- NULL

while (j <= inputlength)  {
if (grepl("SampleName[[:blank:]]*=", matrixinput[j,1])==TRUE) {
   popno <- popno + 1
   haplist[1,popno] <-  (unlist(strsplit(matrixinput[j,1],'"'))[2])
   k <- j + 1
   while (matrixinput[k,2]==matrixinput[j,2]) {
   isitamatch <-  (unlist(strsplit((str_trim(matrixinput[k,1],side="left")),"[[:blank:]]+",fixed=FALSE))[1])
   matchbad <- c("",NA)
   if(!(isitamatch %in% matchbad)){
   if (grepl("sample", matrixinput[k,1],ignore.case=TRUE)==FALSE) {
   matcharray <- append(matcharray,isitamatch)
   for (m in 2:haplistlength) {
   if (haplist[m,1]==isitamatch) {
   haplist[m,popno] <-  (unlist(strsplit((str_trim(matrixinput[k,1],side="left")),"[[:blank:]]+",fixed=FALSE))[2])
    }
    }
    }
    }
   k <- k + 1                                  }
                                                  }
j <- j + 1
                          }

# Error message tested and functioning
lenarray <- length(matcharray)
k <- 0
for (j in 1:lenarray) {
if (!(matcharray[j] %in% haplist[2:haplistlength,1])) {
print(matcharray[j])
k <- k + 1
}
}

if(k>0) {
stop("\nThe haplotype(s) printed above are defined under your population data but aren't present in your haplotype list.\nPlease fix the file and try again.")
}

print(noquote("Not to worry: I have successfully parsed the haplotypes from your haplotype definition block"))
print(noquote(""))
flush.console()

haplist[is.na(haplist)] <- 0

# Error message tested and functioning
if (dim(haplist)[2] >= 4) {
print(noquote("Your arlequin file has more than one population, so permutations can be performed"))
print(noquote(""))
flush.console()
} else {
stop("\n\nYour arlequin file does not appear to contain more than one population, so permutations cannot be performed.\n Please try a file with more than one pop.\n\n")
}

# Getting some parameters that will be used whether haplotype and/or nucleotide diversity is being tested
i <- NULL
j <- NULL
k <- NULL
pop_size <- NULL # defining variables used in the for loop below

no_haps <- dim(haplist)[1] - 1  #Getting the number of haplotypes from the input file
no_pops <- dim(haplist)[2] - 2  #Getting the number of populations from the input file

for (i in 1:no_pops) {
j <- sum(as.numeric(haplist[1:no_haps+1, i+2])) # summing the total number of samples in population i
pop_size[i] <- j # adding the pop size for i to 'pop_size', an array recording all population sizes
}

# Error message tested and functioning
if(0 %in% pop_size) {
stop("\n\nYour arlequin file contains a population with no data. Please remove/correct this population and try again\n\n")
}

print(noquote("I am now going to check for duplicate haplotypes in your haplotype definitions"))
flush.console()
 
# We need to check if any of the haplotypes are identical and pool these if so
# Defining our variables for the loop below
pattern <- NULL
propdiffs <- matrix(NA,nrow = (no_haps + 1),ncol = (no_haps + 1))
propdiffs[1,] <- t(haplist[,1])
propdiffs[,1] <- haplist[,1]
if(!(is.null(missingdata))) {
ambigs <- c("R","Y","S","W","K","M","B","D","H","V","N",missingdata,"r","y","s","w","k","m","b","d","h","v","n")
} else {
ambigs <- c("R","Y","S","W","K","M","B","D","H","V","N","r","y","s","w","k","m","b","d","h","v","n")
}

OKbases <- c("A","C","G","T","-","a","c","g","t")

# Calculating the proportion of differences between haplotypes for all base positions that do not have an ambiguous nucleotide
for (j in 2:(no_haps)) {
   first <- unlist(strsplit(haplist[j,2],pattern))
   for (k in (1+j):(no_haps+1)) {
      second <- unlist(strsplit(haplist[k,2],pattern))
      mismatch <- which(toupper(first)!=toupper(second)) 
      for (m in mismatch) {
         if (!(toupper(first[m]) %in% toupper(ambigs))) {
            if (!(toupper(first[m]) %in% toupper(OKbases))) {
               stop("\n\nThere are non-IUPAC codes in your DNA sequence data. Please check this and try again.\n\n")
            }
         }
         if (!(toupper(second[m]) %in% toupper(ambigs))) {
            if (!(toupper(second[m]) %in% toupper(OKbases))) {
               stop("\n\nThere are non-IUPAC codes in your DNA sequence data. Please check this and try again.\n\n")
            }
         }
      }
     ambig_sites <- unique(c(which((toupper(first[mismatch]) %in% toupper(ambigs))),which((toupper(second[mismatch]) %in% toupper(ambigs)))))
     mismatch <- mismatch[-ambig_sites]
     propdiffs[k,j] <- length(mismatch)/length(first)
   }
}

# Calculating the amount of missing data between haplotypes   
j <- 2
namearray <- NULL
missingamount <- NULL

while (j <= no_haps) {
namearray[j] <- haplist[j,1]
if (!(is.null(missingdata))) {
missingamount[j] <- nchar(haplist[j,2])-nchar(gsub(missingdata,"",haplist[j,2],fixed=TRUE))
} else {
missingamount[j] <- 0
}
k <- j + 1
m <- 0
while (k <= (no_haps+1)) {
if (propdiffs[k,j]==0) {
m <- m + 1
namearray[j] <- paste(namearray[j],haplist[k,1])
if (!(is.null(missingdata))) {
missingamount[j] <- paste(missingamount[j],(nchar(haplist[k,2])-nchar(gsub(missingdata,"",haplist[k,2],fixed=TRUE))))
} else {
missingamount[j] <- paste(missingamount[j],0)
}
}
k <- k + 1
}
j <- j + 1
}
   
newnamearray <- NULL
newmissingamount <- NULL

namearraylen <- length(namearray)
for (j in 1:namearraylen) {
if ((length(unlist(strsplit(namearray[j]," "))))>1) {
newnamearray <- rbind(newnamearray,namearray[j])
newmissingamount <- rbind(newmissingamount,missingamount[j])
}
}

#This for loop is figuring out which haplotype might be identical, given ambiguous sites
if(length(newnamearray)>0) {
newmissingamount <- cbind(newmissingamount,NA)
newnamearray <- cbind(newnamearray,NA)
newnamearraylen <- dim(newnamearray)[1]
j <- 1
i <- 1

finnamearray <- NULL
finmissingamount <- NULL

if(newnamearraylen==1) {
finnamearray[i] <- newnamearray[j]
finmissingamount[i] <- newmissingamount[j]
} else {
for (j in 1:(newnamearraylen)) {
if(is.na(newnamearray[j,2])) {
newnamearray[j,2] <- i
newmissingamount[j,2] <- i
rowlen <- length(unlist(strsplit(newnamearray[j,1], " ")))
for (m in 1:rowlen) {
k <- j + 1
while (k <= newnamearraylen) {
if(is.na(newnamearray[k,2])) {
if(((unlist(strsplit(newnamearray[j,1]," "))) [m]) %in% (unlist(strsplit(newnamearray[k,1]," ")))) {
newnamearray[k,2] <- i
newmissingamount[k,2] <- i
}
}
k <- k + 1
}
}
i <- i + 1
}
j <- j + 1
}

if(is.na(newnamearray[newnamearraylen,2])) {
newnamearray[newnamearraylen,2] <- i + 1
newmissingamount[newnamearraylen,2] <- i + 1
}

newnamearray <- cbind(newnamearray,newmissingamount)
newnamearray <- newnamearray[order(newnamearray[,2]),]
i <- 1
j <- 1

while (j < newnamearraylen) {
k <- j + 1
finnamearray[i] <- newnamearray[j,1]
finmissingamount[i] <- newnamearray[j,3]
while ((k <= newnamearraylen) && (newnamearray[j,2]==newnamearray[k,2])) {
finnamearray[i] <- paste(finnamearray[i],newnamearray[k,1])
finmissingamount[i] <- paste(finmissingamount[i],newnamearray[k,3])
k <- k + 1
}
i <- i + 1
j <- k
}

if(!(newnamearray[newnamearraylen,2]==newnamearray[(newnamearraylen-1),2])) {
finmissingamount[length(finnamearray)+1] <- newnamearray[newnamearraylen,3] 
finnamearray[length(finnamearray)+1] <- newnamearray[newnamearraylen,1]
}
}


temphaplist <- haplist[1,]
uniques_for_adding <- NULL
j <- 1

while (j <= length(finnamearray)) {
addsies <- NULL
secondcols <- rep.int(0,no_pops)
vectorizing <- unlist(strsplit(finnamearray[j], " "))
checkformults <- as.data.frame(table(vectorizing))
checkformults <- checkformults[order(checkformults[,2]),]
checkformults <- droplevels(checkformults)
checkformultslen <- dim(checkformults)[1]

if(checkformultslen>2) {
for (k in 1:(checkformultslen-2)) {
if (!((checkformults[k,2]+1)==checkformults[(k+1),2])) {
print(noquote(""))
print(noquote("You have haplotype(s) with so much missing data they could be a match"))
print(noquote("to multiple different haplotypes. Please remove these haplotypes and try again"))
toprint <- vectorizing[which.max(unlist(strsplit(finmissingamount[j], " ")))]
print(noquote(toprint))
flush.console()
stop("The above haplotype has the most missing data of the matching haplotypes. Remove it and try again")
}
}

if(!(checkformults[checkformultslen,2]<=(checkformults[(checkformultslen-1),2]+1))) {
print(noquote(""))
print(noquote("You have a haplotype with so much missing data it could be a match"))
print(noquote("to multiple different haplotypes. Please remove this haplotype and try again"))
toprint <- vectorizing[which.max(unlist(strsplit(finmissingamount[j], " ")))]
print(noquote(toprint))
flush.console()
stop("Please remove the above haplotype which has too much missing data and try again")
}
}

onlyunique <- vectorizing[!duplicated(vectorizing)]
print(noquote("The following set of haplotypes is identical"))
print(noquote(onlyunique))
flush.console()
uniques_for_adding <- append(uniques_for_adding,onlyunique)
for (k in 1:length(onlyunique)) {
for (m in 2:haplistlength) {
if (haplist[m,1]==onlyunique[k]) {
firstcols <- haplist[m,1:2]
secondcols <- secondcols + as.numeric(haplist[m,3:(no_pops+2)])
}
}
}
addsies <- cbind((t(as.matrix(firstcols))),(t(as.matrix(secondcols))))
temphaplist <- rbind(temphaplist,addsies)
j <- j + 1
}

for (m in 2:haplistlength) {
if (!(haplist[m,1] %in% uniques_for_adding)) {
temphaplist <- rbind(temphaplist,haplist[m,])
}
}

haplist <- temphaplist
print(noquote("The frequency for sets of identical haplotypes has been merged"))
print(noquote("The modified input table is in the process of being written to:"))
print("haplist.txt")
print(noquote(""))
flush.console()
write.table(haplist, "haplist.txt", sep="\t",quote=FALSE, row.names=FALSE,col.names=FALSE)

no_haps <- dim(haplist)[1] - 1 
pattern <- NULL
propdiffs <- matrix(NA,nrow = (no_haps + 1),ncol = (no_haps + 1))
propdiffs[1,] <- t(haplist[,1])
propdiffs[,1] <- haplist[,1]

if(!(is.null(missingdata))) {
ambigs <- c("R","Y","S","W","K","M","B","D","H","V","N",missingdata,"r","y","s","w","k","m","b","d","h","v","n")
} else {
ambigs <- c("R","Y","S","W","K","M","B","D","H","V","N","r","y","s","w","k","m","b","d","h","v","n")
}

OKbases <- c("A","C","G","T","-","a","c","g","t")

for (j in 2:(no_haps)) {
   first <- unlist(strsplit(haplist[j,2],pattern))
   for (k in (1+j):(no_haps+1)) {
      second <- unlist(strsplit(haplist[k,2],pattern))
      mismatch <- which(toupper(first)!=toupper(second)) 
      for (m in mismatch) {
         if (!(toupper(first[m]) %in% toupper(ambigs))) {
            if (!(toupper(first[m]) %in% toupper(OKbases))) {
               stop("\n\nThere are non-IUPAC codes in your DNA sequence data. Please check this and try again.\n\n")
            }
         }
         if (!(toupper(second[m]) %in% toupper(ambigs))) {
            if (!(toupper(second[m]) %in% toupper(OKbases))) {
               stop("\n\nThere are non-IUPAC codes in your DNA sequence data. Please check this and try again.\n\n")
            }
         }
      }
     ambig_sites <- unique(c(which((toupper(first[mismatch]) %in% toupper(ambigs))),which((toupper(second[mismatch]) %in% toupper(ambigs)))))
     mismatch <- mismatch[-ambig_sites]
     propdiffs[k,j] <- length(mismatch)/length(first)
   }
}

   
   
   
} else {
print(noquote("No identical haplotypes were found in your haplotype definition"))
print(noquote("The input table has been written to:"))
print("haplist.txt")
print(noquote(""))
flush.console()
write.table(haplist, "haplist.txt", sep="\t",quote=FALSE, row.names=FALSE,col.names=FALSE)
}

i <- NULL
j <- NULL
k <- NULL
haplotypes <- NULL  # defining variables used in the for loop below

for (i in 1:no_haps) {
j <- sum(as.numeric(haplist[i+1,3:(no_pops+2)])) # summing the total number of haplotype i over all populations
k <- rep.int((i+1),j) # coding haplotype i as a numeric value
haplotypes <- c(haplotypes, k) # adding the values for haplotype i to the previous array of haplotypes
}

total_no_samples <- length(haplotypes)

#Testing for haplotype diversity if that option was selected
if (test_for_hap=="Y") {
print(noquote("STEP TWO: You've said you would like to test for haplotype diversity, so doing that now"))
print(noquote(""))
flush.console()

i <- NULL # defining variables used in the for loop below
obshd <- dim(no_pops)
obshi <- dim(no_pops)

for (i in 1:no_pops) {# This loop is calculating the observed haplotype diversity for our populations
obshi[i] <- 0
y <- 1
while (y < no_haps+1) {
# While y is less than or equal to the total number of haplotypes in the population
obshi[i] <- obshi[i] + (as.numeric(haplist[y+1,i+2])/pop_size[i])^2
# Count the number of haplotypes in each population that equals
# haplotype code 'y+1'. Divide them by the total sample size (to get the
# frequency of that haplotype in the sample) and square this
# difference. Add this to the value calculated for other haplotypes.
y <- y+1
}
obshd[i] <- (pop_size[i]/(pop_size[i] - 1))*(1-obshi[i])
# Calculating haplotype diversity from the haplotype identity (hi_) calculated above
}

i <- NULL
j <- NULL
k <- NULL
temp <- NULL
x <- NULL
y <- NULL
hd_ <- NULL
hi_ <- matrix(0,nrow = n_iter,ncol = no_pops)
hd_ <- matrix(0,nrow = n_iter,ncol = no_pops)  # defining variables for below

print(noquote("I am going to give you progress updates at 1%, 10%, 25%, 50%, 75% and 90%"))
flush.console()

for (i in 1:n_iter) {
if((i/n_iter)==0.01) {
print(noquote("About 1% done with haplotype diversity permutations"))
flush.console()
}
if((i/n_iter)==0.1) {
print(noquote("About 10% done with haplotype diversity permutations"))
flush.console()
}
if((i/n_iter)==0.25) {
print(noquote("About 25% done with haplotype diversity permutations"))
flush.console()
}
if((i/n_iter)==0.5) {
print(noquote("About 50% done with haplotype diversity permutations"))
flush.console()
}
if((i/n_iter)==0.75) {
print(noquote("About 75% done with haplotype diversity permutations"))
flush.console()
}
if((i/n_iter)==0.9) {
print(noquote("About 90% done with haplotype diversity permutations"))
flush.console()
}

# This sampling process will be completed n_iter times
x <- 1
temp <- sample(haplotypes, total_no_samples, replace=FALSE, prob=NULL) # rearranging order of haplotype array
for (j in 1:no_pops) {
if (!(j == 1)) {x <- x + pop_size[j-1]}
pop_haps <- temp[x:(x+pop_size[j]-1)] # getting the simulated haplotypes that will correspond to each population

y <- 1

while (y < no_haps+1) {
# While y is less than or equal to the total number of haplotypes in the population
hi_[i,j] <- hi_[i,j] + ((sum(pop_haps==y))/pop_size[j])^2
# Count the number of haplotypes in each population that equals
# haplotype code 'j'. Divide them by the total sample size (to get the
# frequency of that haplotype in the sample) and square this
# difference. Add this to the value calculated for other haplotypes.


y <- y+1
# The counter for the while loop goes up one (we move on and do these 
# calculations for the next haplotype)
}
hd_[i,j] <- pop_size[j]/(pop_size[j] - 1)*(1-hi_[i,j])
# Calculating haplotype diversity from the haplotype identity (hi_) calculated above
}
}

# Defining variables for the loop where we will record our results
x <- (no_pops*(no_pops - 1))/2
results <- matrix(,nrow = (n_iter+1),ncol = x)

i <- 1
j <- 1
k <- NULL
x <- 1
obs_diff <- NULL 

summary_results <-  matrix("",nrow = (no_pops+2),ncol = (no_pops+1))
summary_results[1,1] <-  c("\"Haplotype diversity for each population given to right\"")
summary_results[2,1] <-  c("\"Haplotype diversity difference between populations below diagonal, p-values above diagonal\"")
summary_results[2,2:(no_pops+1)] <- haplist[1,3:(no_pops+2)]
summary_results[3:(no_pops+2),1] <- t(haplist[1,3:(no_pops+2)])

print(noquote(""))
print(noquote("####RESULTS#### Haplotype diversity permutation test"))
print(noquote(""))
flush.console()

while (j <= no_pops) {
summary_results[j+2,j+1] <- "--"
summary_results[1,j+1] <- obshd[j]
k <- j + 1
while (k <= no_pops){
results[1,x] <- paste(haplist[1,(j+2)],haplist[1,(k+2)])
for (i in 1:n_iter) {
results[i+1,x] <- abs(hd_[i,j] - hd_[i,k])  # calculAting the simulated haplotype diversity differences between each population pair
}
obs_diff <- abs(obshd[j]-obshd[k]) # Comparing the simulated values to our observed (printing into STDOUT)
p_val <- sum(as.numeric(results[2:(n_iter+1),x])>=obs_diff)
print1 <- paste("Comparing",paste(haplist[1,(j+2)]),"(h =",paste(obshd[j]),") &",paste(haplist[1,(k+2)]),"(h =",paste(obshd[k]),"):")
print(noquote(print1))
print2 <- paste(p_val,"out of ",n_iter," simulated results had a haplotype diversity difference")
print(noquote(print2))
print3 <- paste("equal or greater than the observed difference of", obs_diff)
print(noquote(print3))
print(noquote(""))
flush.console()
summary_results[k+2,j+1] <- obs_diff
summary_results[j+2,k+1] <- as.numeric(p_val)/n_iter
x <- x + 1
k <- k + 1
}
j <- j + 1
}

write.table(summary_results, "hapl_diffs_summary.txt", sep="\t",quote=FALSE, row.names=FALSE,col.names=FALSE)
write.table(results, "haplotype_diversity_differences.txt", sep="\t",quote=FALSE, row.names=FALSE,col.names=FALSE)

print(noquote("I have also written out the simulated haplotype diversity for each population, for each iteration, to:"))
print("haplotype_diversity_differences.txt")
print(noquote("I have written a summarized table of differences and p-values to:"))
print("hapl_diffs_summary.txt")
print(noquote(""))
flush.console()

} else {
print(noquote("STEP TWO: You want to skip a test for haplotype diversity... moving on"))
print(noquote(""))
flush.console()
}


#Testing for nucleotide diversity if that option was selected
if (test_for_nucl=="Y") {
print(noquote("STEP THREE: You've said you would like to test for differences in nucleotide diversity, so doing that now"))
print(noquote(""))
flush.console()

obspi <- NULL

for (j in 1:no_pops) {# This loop is calculating the observed nucleotide diversity for our populations
obspi[j] <- 0
for (k in 1:no_haps) {# For every haplotype...
if (haplist[k+1,j+2]!=0) {
m <- k + 1
while (m <= no_haps) {
if (haplist[m+1,j+2]!=0) {
obspi[j] <- obspi[j] + 2*as.numeric(propdiffs[m+1,k+1])*(as.numeric(haplist[(k+1),(j+2)]))*(as.numeric(haplist[(m+1),(j+2)]))/(pop_size[j]*(pop_size[j]-1))
# For every other haplotype in population j, multiplying the proportion of differences between the haplotypes
# by the frequency of those haplotypes in population j. 
# Adding this to the value calculated for other haplotypes.


}
m <- m + 1
}
}
}
}

simpi <- matrix(0,nrow = n_iter,ncol = no_pops) 

print(noquote("I am going to give you progress updates at 1%, 10%, 25%, 50%, 75% and 90%"))
flush.console()

for (i in 1:n_iter) {
if((i/n_iter)==0.01) {
print(noquote("About 1% done with nucleotide diversity permutations"))
flush.console()
}
if((i/n_iter)==0.1) {
print(noquote("About 10% done with nucleotide diversity permutations"))
flush.console()
}
if((i/n_iter)==0.25) {
print(noquote("About 25% done with nucleotide diversity permutations"))
flush.console()
}
if((i/n_iter)==0.5) {
print(noquote("About 50% done with nucleotide diversity permutations"))
flush.console()
}
if((i/n_iter)==0.75) {
print(noquote("About 75% done with nucleotide diversity permutations"))
flush.console()
}
if((i/n_iter)==0.9) {
print(noquote("About 90% done with nucleotide diversity permutations"))
flush.console()
}

# This sampling process will be completed n_iter times
x <- 1
temp <- sample(haplotypes, total_no_samples, replace=FALSE, prob=NULL) # rearranging order of haplotype array
for (j in 1:no_pops) {
if (!(j == 1)) {x <- x + pop_size[j-1]}
pop_haps <- temp[x:(x+pop_size[j]-1)] # getting the simulated haplotypes that will correspond to each population
summedhaps <- as.matrix(as.data.frame(table(pop_haps)))
lensummedhaps <- dim(summedhaps)[1]

if(lensummedhaps>1){
for (k in 1:(lensummedhaps-1)) {
hapid1 <- as.numeric(summedhaps[k,1])
m <- k + 1
while (m <= lensummedhaps) {
hapid2 <- as.numeric(summedhaps[m,1])
simpi[i,j] <- simpi[i,j] + (2*as.numeric(propdiffs[hapid2,hapid1])*as.numeric(summedhaps[k,2])*as.numeric(summedhaps[m,2]))/(pop_size[j]*(pop_size[j]-1))
m <- m + 1
# For eeach haplotype in simulated population j, multiplying the proportion of differences between the haplotypes
# by the frequency of those haplotypes in population j. 
# Adding this to the value calculated for other haplotypes.
}
}
} else {
simpi[i,j] <- 0
}
}
}

# Defining variables for the loop where we will record our results
x <- (no_pops*(no_pops - 1))/2
results <- matrix(,nrow = (n_iter+1),ncol = x)
i <- 1
j <- 1
k <- NULL
x <- 1
obs_diff <- NULL

summary_results <-  matrix("",nrow = (no_pops+2),ncol = (no_pops+1))
summary_results[1,1] <-  c("\"Nucleotide diversity for each population given to right\"")
summary_results[2,1] <-  c("\"Nucleotide diversity difference between populations below diagonal, p-values above diagonal\"")
summary_results[2,2:(no_pops+1)] <- haplist[1,3:(no_pops+2)]
summary_results[3:(no_pops+2),1] <- t(haplist[1,3:(no_pops+2)])


print(noquote(""))
print(noquote("####RESULTS#### Nucleotide diversity permutation test"))
print(noquote(""))
flush.console()


while (j <= no_pops) {
k <- j + 1
summary_results[j+2,j+1] <- "--"
summary_results[1,j+1] <- obspi[j]
while (k <= no_pops){
results[1,x] <- paste(haplist[1,(j+2)],haplist[1,(k+2)])
for (i in 1:n_iter) {
results[i+1,x] <- abs(simpi[i,j] - simpi[i,k])  # calcualting the simulated haplotype diversity differences between each population pair
}
obs_diff <- abs(obspi[j]-obspi[k]) # Comparing the simulated values to our observed (printing into STDOUT)
p_val <- sum(as.numeric(results[2:(n_iter+1),x])>=obs_diff)
print1 <- paste("Comparing",paste(haplist[1,(j+2)]),"(pi =",paste(obspi[j]),") &",paste(haplist[1,(k+2)]),"(pi =",paste(obspi[k]),"):")
print(noquote(print1))
print2 <- paste(p_val,"out of ",n_iter," simulated results had a nucleotide diversity difference")
print(noquote(print2))
print3 <- paste("equal or greater than the observed difference of", obs_diff)
print(noquote(print3))
print(noquote(""))
flush.console()
summary_results[k+2,j+1] <- obs_diff
summary_results[j+2,k+1] <- as.numeric(p_val)/n_iter
x <- x + 1
k <- k + 1
}
j <- j + 1
}

write.table(summary_results, "nucl_diffs_summary.txt", sep="\t",quote=FALSE, row.names=FALSE,col.names=FALSE)
write.table(results, "nucleotide_diversity_differences.txt", sep="\t",quote=FALSE, row.names=FALSE,col.names=FALSE)

print(noquote("I have also written out the simulated nucleotide diversity for each population, for each iteration, to:"))
print("nucleotide_diversity_differences.txt")
print(noquote("I have written a summarized table of differences and p-values to:"))
print("nucl_diffs_summary.txt")
print(noquote(""))
flush.console()

} else {
print(noquote("STEP THREE: You want to skip a test for nucleotide diversity... program is all done!"))
flush.console()
}
}
