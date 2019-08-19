#!/usr/bin/env Rscript

# ---------- Preparations ----------
# Load libraries
library(parallel)               # Detect number of cpu cores
library(foreach)                # For multicore parallel
library(doMC)                   # For multicore parallel
library(RColorBrewer)           # For colors
library(multtest)               # For diffreport
library(MSnbase)                # MS features
library(xcms)                   # Swiss army knife for metabolomics
library(CAMERA)                 # Metabolite Profile Annotation
library(mixOmics)               # Statistics for various Omics
# ERROR: vegan conflicts with mda and klaR, unload packages before using any of the analyses !!!
if ("package:mda" %in% search()) detach(package:mda, unload=TRUE)
if ("package:klaR" %in% search()) detach(package:klaR, unload=TRUE)
library(vegan)
library(multcomp)               # For Tukey test
library(Hmisc)                  # For correlation test
library(gplots)                 # For fancy heatmaps
library(lme4)                   # Linear Models with Random Effects
library(lmerTest)               # Create p-values for LMER
library(ape)                    # Phylogeny
library(pvclust)                # Phylogeny
library(dendextend)             # Phylogeny
library(cba)                    # Phylogeny
library(phangorn)               # Phylogeny
library(ontologyIndex)          # Reading obo ontology files
library(webchem)                # Converting InChI to InChIKey

# Set locales and encoding
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
loc <- Sys.setlocale(category="LC_ALL", locale="C")
options(encoding="UTF-8")

# Set options
options(stringAsfactors=FALSE, useFancyQuotes=FALSE)

# Multicore parallel
nSlaves <- detectCores(all.tests=FALSE, logical=TRUE)
registerDoMC(nSlaves)



# ---------- User variables ----------
# Working directory
setwd("~/mtbls709/statistics")

# Data directory
mzml_dir <- "~/mtbls709/data/raw"

# MS1 variables
polarity <- "positive"
pol <- substr(x=polarity, start=1, stop=3)
ppm <- 35
peakwidth <- c(4,21)
snthresh <- 10
prefilter <- c(5,50)
fitgauss <- FALSE
verbose.columns <- FALSE
mzwidth <- 0.01
minfrac <- 0.5
bwindow <- 4

# MS2 variables
mzabs <- 0.01 #0.01                        # Absolute mass error (in seconds) used for merging MS/MS spectra
mzppm <- 5 #5                         # ppm error used for merging MS/MS spectra
rtabs <- 5 #20                           # Retention time error (in seconds) used for merging MS/MS spectra
max.rt.range <- 20 #20                    # Permitted retention time window (in seconds) of grouped MS1 precursors
max.mz.range <- 0.01 #0.01                 # Permitted m/z window of grouped MS1 precursors
min.rt <- 10                          # Minimum retention time for selected precursors
max.rt <- 1020                        # Maximum retention time for selected precursors
min.mz <- 50                          # Minimum m/z value for selected precursors
max.mz <- 1500                        # Maximum m/z value for selected precursors
msms.intensity.threshold <- 100       # Minimum intensity value for MS/MS peaks

# Preparations for plotting
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)



# ---------- MS1 Preparations ----------
source("ms1_extraction.r")

# Load files and define classes
f.load_mzml()
f.sample_classes()

# MS1 Peak picking
f.ms1_pick_features()
f.ms1_preprocess_features(sclass=species)

# Export MAF
f.export_maf(xcam=peak_xcam, maf_filename="ms1.maf.tsv")

# Annotate MAF
f.annotate_maf(csv_compounds="annotation.csv", maf_input="ms1.maf.tsv", maf_output="ms2.maf.tsv")



# ---------- MS1: Ecological analyses ----------
source("ms1_diversity.r")

# Diversity
f.ms1_div_model()

# Import traits
f.ms1_import_traits(filename="../data/traits.csv")

# Species
f.species_div_unique()
f.species_div_features()
f.species_div_shannon()
f.species_div_pielou()
f.species_div_concentration()

# Seasons
f.seasons_div_unique()
f.seasons_div_features()
f.seasons_div_shannon()
f.seasons_div_pielou()
f.seasons_div_concentration()



# ---------- MS2 Spectra annotation & classification ----------
source("ms2_merge_filter.r")
source("ms2_classification.r")
source("ms2_annotation.r")

# MS2 Spectra Filtering and merging
f.ms2_find_spectra()
f.ms2_plot_spectra()

# Export MS2 spectra to MSP
f.ms2_create_msp()

# Calculate Relationships to phylogeny
f.ms2_classify_phylo()

# Generalized MS2 Spectra classification
f.ms2_classify_features()

# Variation partitioning
f.ms2_classify_varpart()

# Which ecological traits have influence on compound classes?
f.ms2_classify_eco()


