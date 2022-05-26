#setwd("1_QC_GWAS/")
library(ggplot2)

### Step 1 ### 

# Investigate missingness per individual and per SNP and make histograms.
system("./plink --bfile HapMap_3_r3_1 --missing")
# output: ./plink.imiss and ./plink.lmiss, these files show respectively the proportion of missing SNPs per individual and the proportion of missing individuals per SNP.


# Generate plots to visualize the missingness results.
#system("Rscript --no-save hist_miss.R")

indmiss<-read.table(file="plink.imiss", header=TRUE)
snpmiss<-read.table(file="plink.lmiss", header=TRUE)

# read data into R 

# pdf("histimiss.pdf") #indicates pdf format and gives title to file
# hist(indmiss[,6],main="Histogram individual missingness") #selects column 6, names header of file
# dev.off()

head(indmiss)
ggplot()+
  geom_histogram(data=indmiss,aes(x=F_MISS))+
  theme_light()

head(snpmiss)
ggplot()+
  geom_histogram(data=snpmiss,aes(x=F_MISS))+
  theme_light()


# Delete SNPs and individuals with high levels of missingness, explanation of this and all following steps can be found in box 1 and table 1 of the article mentioned in the comments of this script.
# The following two QC commands will not remove any SNPs or individuals. However, it is good practice to start the QC with these non-stringent thresholds.  
# Delete SNPs with missingness >0.2.
system("./plink --bfile HapMap_3_r3_1 --geno 0.2 --make-bed --out HapMap_3_r3_2")

# Delete individuals with missingness >0.2.
system("./plink --bfile HapMap_3_r3_2 --mind 0.2 --make-bed --out HapMap_3_r3_3")

# Delete SNPs with missingness >0.02.
system("./plink --bfile HapMap_3_r3_3 --geno 0.02 --make-bed --out HapMap_3_r3_4")

# Delete individuals with missingness >0.02.
system("./plink --bfile HapMap_3_r3_4 --mind 0.02 --make-bed --out HapMap_3_r3_5")

#
###################################################
### Step 3 ### 

# Generate a bfile with autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF).

# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).

system("./plink --bfile HapMap_3_r3_6 --chr 1-22 --make-bed --out HapMap_3_r3_7")

# Generate a plot of the MAF distribution.

system("./plink --bfile HapMap_3_r3_7 --freq --out MAF_check")

#system("Rscript --no-save MAF_check.R")

maf_freq <- read.table("MAF_check.frq", header =TRUE, as.is=T)

head (maf_freq)

ggplot()+
  geom_histogram(data=maf_freq,aes(x=MAF))

# Remove SNPs with a low MAF frequency.
system("./plink --bfile HapMap_3_r3_7 --maf 0.05 --make-bed --out HapMap_3_r3_8")
# 1073226 SNPs are left
# A conventional MAF threshold for a regular GWAS is between 0.01 or 0.05, depending on sample size.

####################################################
### Step 4 ###

# Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).
# Check the distribution of HWE p-values of all SNPs.

# By default the --hwe option in ./plink only filters for controls.
# Therefore, we use two steps, first we use a stringent HWE threshold for controls, followed by a more stringent threshold for the case data.
system("./plink --bfile HapMap_3_r3_8 --hwe 1e-6 --make-bed --out HapMap_hwe_filter_step1")

# The HWE threshold for the cases filters out only SNPs which deviate extremely from HWE. 
# This second HWE step only focusses on cases because in the controls all SNPs with a HWE p-value < hwe 1e-6 were already removed
system("./plink --bfile HapMap_hwe_filter_step1 --hwe 1e-10 --hwe-all --make-bed --out HapMap_3_r3_9")

# Theoretical background for this step is given in our accompanying article: https://www.ncbi.nlm.nih.gov/pubmed/29484742 .

############################################################
### step 5 ###

# Generate a plot of the distribution of the heterozygosity rate of your subjects.
# And remove individuals with a heterozygosity rate deviating more than 3 sd from the mean.

# Checks for heterozygosity are performed on a set of SNPs which are not highly correlated.
# Therefore, to generate a list of non-(highly)correlated SNPs, we exclude high inversion regions (inversion.txt [High LD regions]) and prune the SNPs using the command --indep-pairwise?.
# The parameters ?50 5 0.2? stand respectively for: the window size, the number of SNPs to shift the window at each step, and the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously.

system("./plink --bfile HapMap_3_r3_9 --exclude inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP")
# Note, don't delete the file indepSNP.prune.in, we will use this file in later steps of the tutorial.

system("./plink --bfile HapMap_3_r3_9 --extract indepSNP.prune.in --het --out R_check")
# This file contains your pruned data set.

# Plot of the heterozygosity rate distribution

het <- read.table("R_check.het", head=TRUE)
pdf("heterozygosity.pdf")
head(het)
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."

ggplot()+
  geom_histogram(data=het,aes(x=HET_RATE),)



# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
# For data manipulation we recommend using UNIX. However, when performing statistical calculations R might be more convenient, hence the use of the Rscript for this step:

outliers <- het[abs(scale(het$HET_RATE))>3,]

write.table(outliers[,1:2],file = "het_fail_ind.txt",quote = F, row.names = F)

# Output of the command above: fail-het-qc.txt .
# When using our example data/the HapMap data this list contains 2 individuals (i.e., two individuals have a heterozygosity rate deviating more than 3 SD's from the mean).
# Adapt this file to make it compatible for ./plink, by removing all quotation marks from the file and selecting only the first two columns.

# Remove heterozygosity rate outliers.
system("./plink --bfile HapMap_3_r3_9 --remove het_fail_ind.txt --make-bed --out HapMap_3_r3_10")


############################################################
### step 6 ###

# It is essential to check datasets you analyse for cryptic relatedness.
# Assuming a random population sample we are going to exclude all individuals above the pihat threshold of 0.2 in this tutorial.

# # Check for relationships between individuals with a pihat > 0.2.
system("./plink --bfile HapMap_3_r3_10 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2")
# 
# # The HapMap dataset is known to contain parent-offspring relations. 
# # The following commands will visualize specifically these parent-offspring relations, using the z values. 

# The generated plots show a considerable amount of related individuals (explentation plot; PO = parent-offspring, UN = unrelated individuals) in the Hapmap data, this is expected since the dataset was constructed as such.
# Normally, family based data should be analyzed using specific family based methods. In this tutorial, for demonstrative purposes, we treat the relatedness as cryptic relatedness in a random population sample.
# In this tutorial, we aim to remove all 'relatedness' from our dataset.
# To demonstrate that the majority of the relatedness was due to parent-offspring we only include founders (individuals without parents in the dataset).

system("./plink --bfile HapMap_3_r3_10 --filter-founders --make-bed --out HapMap_3_r3_11")

# Now we will look again for individuals with a pihat >0.2.
system("./plink --bfile HapMap_3_r3_11 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders")
# The file 'pihat_min0.2_in_founders.genome' shows that, after exclusion of all non-founders, only 1 individual pair with a pihat greater than 0.2 remains in the HapMap data.
# This is likely to be a full sib or DZ twin pair based on the Z values. Noteworthy, they were not given the same family identity (FID) in the HapMap data.

# For each pair of 'related' individuals with a pihat > 0.2, we recommend to remove the individual with the lowest call rate. 
system("./plink --bfile HapMap_3_r3_11 --missing")
# Use an UNIX text editor (e.g., vi(m) ) to check which individual has the highest call rate in the 'related pair'. 

# Generate a list of FID and IID of the individual(s) with a Pihat above 0.2, to check who had the lower call rate of the pair.
# In our dataset the individual 13291  NA07045 had the lower call rate.
system("echo '13291  NA07045' > 0.2_low_call_rate_pihat.txt")

# Delete the individuals with the lowest call rate in 'related' pairs with a pihat > 0.2 
system("./plink --bfile HapMap_3_r3_11 --remove 0.2_low_call_rate_pihat.txt --make-bed --out HapMap_3_r3_12")

################################################################################################################################

# CONGRATULATIONS!! You've just succesfully completed the first tutorial! You are now able to conduct a proper genetic QC. 

# For the next tutorial, using the script: 2_Main_script_MDS.txt, you need the following files:
# - The bfile HapMap_3_r3_12 (i.e., HapMap_3_r3_12.fam,HapMap_3_r3_12.bed, and HapMap_3_r3_12.bim
# - indepSNP.prune.in

