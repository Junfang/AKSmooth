AKSmooth
========

Ajusted local kernel smoothing for low-coverage bisulfite sequencing data.


1. INSTALL

# Download and unzip this file, and run the following codes in R:

install.packages("devtools") 
library(devtools)

# Find the directory where AKSmooth-master file is located and set the following new working directory.
setwd("./AKSmooth-master") 

# Load the package into memory
load_all()

# Install the package 
install()

2. USAGE

# After installation one can reload the AKSmooth package.
library(AKSmooth)

# Example
load("./data/bn1chr21.rda")
data <- bn1chr21
fitChr21gau <- AKSmooth(data, 30, "Gaussian")

# For more details.
?AKSmooth


Reference

Chen, J., Lutsik, P., Akulenko, R., Walter, J., & Helms, V. (2014). AKSmooth: Enhancing low-coverage bisulfite sequencing data via kernel-based smoothing. Journal of bioinformatics and computational biology, 12(06), 1442005.


