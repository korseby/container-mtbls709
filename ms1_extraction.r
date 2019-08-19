# ---------- Preparations ----------
# General variables
mzml_files <- NULL
mzml_names <- NULL
mzml_times <- NULL
qc_files <- NULL
qc_times <- NULL

# Specific variables
species <- NULL
seasons <- NULL
seasonal_species <- NULL
spesearep <- NULL
species_names <- NULL
species_labels <- c("Brachythecium rutabulum", "Calliergonella cuspidata",
                    "Fissidens taxifolius", "Grimmia pulvinata",
                    "Hypnum cupressiforme", "Marchantia polymorpha",
                    "Plagiomnium undulatum", "Polytrichum strictum",
                    "Rhytidiadelphus squarrosus")
species_colors <- NULL
species_symbols <- NULL
seasons_names <- NULL
seasons_colors <- NULL
seasons_symbols <- NULL
species_samples_colors <- NULL
seasons_samples_colors <- NULL
species_samples_symbols <- NULL
seasons_samples_symbols <- NULL
samples <- NULL
presence <- 12 * (10/12)

# MS1 Matrices
peak_xset <- NULL
peak_xcam <- NULL
peak_qcset <- NULL
peak_qccam <- NULL
diff_list <- NULL
peak_list <- NULL
feat_list <- NULL
bina_list <- NULL
uniq_list <- NULL

# MS2 Matrices
msms_spectra <- NULL
msms_merged <- NULL



# ---------- Load mzML files ----------
f.load_mzml <- function() {
    # Load files
    mzml_files <- list.files(mzml_dir, pattern="*.mzML", recursive=T, full.names=T)
    
    # Exclude blanks
    mzml_files <- mzml_files[grep("MM8", mzml_files, invert=T)]
    mzml_files <- mzml_files[grep("ACN", mzml_files, invert=T)]

    # Basenames of files without path and without extension
    mzml_names <- gsub('(.*)\\..*', '\\1', gsub('( |-|,)', '.', basename(mzml_files)))
    
    # Load files
    qc_files <- list.files(mzml_dir, pattern="*.mzML", recursive=T, full.names=T)
    qc_files <- qc_files[grep("MM8", qc_files, invert=FALSE)]
    
    # Save timestamps of samples
    for (i in 1:length(mzml_files)) {
        fl <- mzR::openMSfile(mzml_files[i])
        run_info <- mzR::runInfo(fl)
        mzR::close(fl)
        mzml_times <- c(mzml_times, run_info$startTimeStamp)
    }
    for (i in 1:length(qc_files)) {
        fl <- mzR::openMSfile(qc_files[i])
        run_info <- mzR::runInfo(fl)
        mzR::close(fl)
        qc_times <- c(qc_times, run_info$startTimeStamp)
    }
    
    # Return global variables
    mzml_files <<- mzml_files
    mzml_names <<- mzml_names
    mzml_times <<- mzml_times
    qc_files <<- qc_files
    qc_times <<- qc_times
}



# ---------- Define sample classes ----------
f.sample_classes <- function() {
    # Sample classes: species
    species <<- as.factor(sapply(strsplit(as.character(mzml_names), "_"), function(x) {
        nam <- x[4];
        nam;
    }))
    
    # Sample classes: seasons
    seasons <<- as.factor(sapply(strsplit(as.character(mzml_names), "_"), function(x) {
        nam <- x[3];
        nam;
    }))
    
    # Sample classes: seasonal species
    seasonal_species <<- as.factor(sapply(strsplit(as.character(mzml_names), "_"), function(x) {
        se <- x[3];
        sp <- x[4];
        nam <- paste(se, '_', sp, sep='')
        nam;
    }))
    
    # Sample classes: unique species-seasons-replicate
    spesearep <<- as.factor(sapply(strsplit(as.character(mzml_files), "_"), function(x) {
        se <- x[3];
        sp <- x[4];
        nam <- paste(sp, '_', se, sep='')
        nam;
    }))
    spesearep <<- make.names(as.character(spesearep), unique=TRUE)
    
    # Define species names, colors, symbols
    species_names <<- levels(species)
    species_colors <<- c("yellowgreen", "mediumseagreen", "darkorange1", "firebrick3", "darkolivegreen4", "dodgerblue4", "chocolate", "darkviolet", "darkkhaki")
    species_symbols <<- c(15, 16, 0, 1, 17, 8, 2, 5, 18)
    
    # Define seasons names, colors, symbols
    seasons_names <<- c("summer", "autumn", "winter", "spring")
    seasons_colors <<- c("darkgoldenrod3", "firebrick3", "deepskyblue3", "chartreuse3")
    seasons_symbols <<- c(1, 15, 16, 0)
    
    # Define samples colors
    species_samples_colors <<- sapply(species, function(x) { x <- species_colors[which(x==species_names)] } )
    seasons_samples_colors <<- sapply(seasons, function(x) { x <- seasons_colors[which(x==seasons_names)] } )
    
    # Define samples symbols
    species_samples_symbols <<- sapply(species, function(x) { x <- species_symbols[which(x==species_names)] } )
    seasons_samples_symbols <<- sapply(seasons, function(x) { x <- seasons_symbols[which(x==seasons_names)] } )
    
    # Define number of samples
    samples <<- length(species)/length(species_names)
}



# ---------- Pick MS1 features ----------
f.ms1_pick_features <- function() {
    # Find Peaks in grouped phenoData according to directory structure
    xset <- xcmsSet(files=mzml_files, method="centWave", BPPARAM=MulticoreParam(nSlaves),
                    ppm=ppm, peakwidth=peakwidth, snthresh=snthresh, prefilter=prefilter,
                    fitgauss=fitgauss, verbose.columns=verbose.columns)
    
    # Normalize filenames
    phenodata_old <- xset@phenoData
    phenodata_new <- xset@phenoData
    phenodata_new[,1] <- as.factor(sapply(strsplit(rownames(phenodata_new), "_"), function(x) { nam <- x[4]; nam; }))
    phenodata_new[,2] <- as.factor(sapply(strsplit(rownames(phenodata_new), "_"), function(x) { nam <- x[3]; nam; }))
    colnames(phenodata_new) <- c("V1", "V2")
    xset@phenoData <- phenodata_new
    
    # Group peaks from different samples together
    xset2 <- group(xset, mzwid=mzwidth, minfrac=minfrac, bw=bwindow)
    
    # Filling in missing peaks
    #xset3 <- fillPeaks(xset2, method="chrom", BPPARAM=MulticoreParam(nSlaves))
    xset3 <- xset2
    
    # Retention time correction
    xset4 <- retcor(xset3, method="loess", family="gaussian", plottype="mdevden", missing=10, extra=1, span=2)
    
    # Peak re-grouping
    xset5 <- group(xset4, mzwid=mzwidth, minfrac=minfrac, bw=bwindow)
    
    # Peak picking with CAMERA
    xcam <- xsAnnotate(xset5, polarity=polarity)
    xcam <- groupFWHM(xcam, perfwhm=0.6)
    xcam <- findIsotopes(xcam, ppm=5, mzabs=0.005)
    xcam <- groupCorr(xcam, calcIso=TRUE, calcCiS=TRUE, calcCaS=TRUE, graphMethod="lpc", pval=0.05, cor_eic_th=0.75)
    xcam <- findAdducts(xcam, polarity=polarity)
    
    # Peak picking including QC
    qcset <- xcmsSet(files=c(mzml_files,qc_files), method="centWave", BPPARAM=MulticoreParam(nSlaves),
                     ppm=ppm, peakwidth=peakwidth, snthresh=snthresh, prefilter=prefilter,
                     fitgauss=fitgauss, verbose.columns=verbose.columns)
    qcset@phenoData$class <- as.factor(c(rep("sample",length(1:length(mzml_files))), rep("qc",length(1:length(qc_files)))))
    qcset2 <- group(qcset, mzwid=mzwidth, minfrac=minfrac, bw=bwindow)
    qcset4 <- retcor(qcset2, method="loess", family="gaussian", plottype="mdevden", missing=10, extra=1, span=2)
    qcset5 <- group(qcset4, mzwid=mzwidth, minfrac=minfrac, bw=bwindow)
    qccam <- xsAnnotate(qcset5, polarity=polarity)
    qccam <- groupFWHM(qccam, perfwhm=0.6)
    qccam <- findIsotopes(qccam, ppm=5, mzabs=0.005)
    qccam <- groupCorr(qccam, calcIso=TRUE, calcCiS=TRUE, calcCaS=TRUE, graphMethod="lpc", pval=0.05, cor_eic_th=0.75)
    qccam <- findAdducts(qccam, polarity=polarity)
    
    # Return xcmsSet and CAMERA objects
    peak_xset <<- xset5
    peak_xcam <<- xcam
    peak_qcset <<- qcset5
    peak_qccam <<- qccam
}



# ---------- Create CSV for each baseline corrected chromatogram ----------
f.ms1_chromatograms_save <- function(mzml_files, mzml_names) {
    # List containing the chromatograms
    xchroms <- list()
    
    for (i in 1:length(mzml_files)) {
        # Prepare chromatogram
        chroma <- xcmsRaw(mzml_files[i])
        x <- chroma@scantime
        y <- scale(chroma@tic, center=FALSE)
        xchrom <- data.frame(x, y)
        colnames(xchrom) <- c("rt", "tic")
        
        # Create list
        xchroms[[i]] <- xchrom
        
        # Write CSV
        write.csv(xchrom, file=paste0("../data/chromatograms/",mzml_names[i],".csv"), row.names=FALSE)
    }
}



# ---------- Group chromatograms together ----------
f.ms1_chromatograms_group_save <- function(mzml_files, mzml_names) {
    # List containing the chromatograms
    xchroms <- list()
    for (i in 1:length(mzml_files)) {
        xchrom <- read.table(paste0("../data/chromatograms/",mzml_names[i],".csv"), header=TRUE, sep=",", quote="\"", fill=TRUE, dec=".", stringsAsFactors=FALSE)
        xchroms[[i]] <- xchrom
    }
    
    # Determine max. retention time of all chromatograms
    rt_max <- 0
    for (i in 1:length(mzml_files)) if (max(xchroms[[i]][,1]) > rt_max) rt_max <- max(xchroms[[i]][,1])
    
    # Determine max. number of measurements of all chromatograms
    m_max <- 0
    for (i in 1:length(mzml_files)) if (max(nrow(xchroms[[i]])) > m_max) m_max <- max(nrow(xchroms[[i]]))
    
    # Determine max. TIC of all chromatograms
    tic_max <- 0
    for (i in 1:length(mzml_files)) if (max(xchroms[[i]]) > tic_max) tic_max <- max(xchroms[[i]])
    
    # Create grouped alignment data frame for chromatograms
    groups <- seq(0, ceiling(rt_max), by=round(rt_max/m_max, 3))
    gchroms <- data.frame(rt=groups)
    
    for (i in 1:length(mzml_files)) {
        gchroms <- cbind(gchroms, 0)
        colnames(gchroms)[ncol(gchroms)] <- strsplit(basename(mzml_files[i]),'.',fixed=T)[[1]][1]
        fint <- findInterval(xchroms[[i]][,1], groups)
        gchroms[fint,i+1] <- xchroms[[i]][,2]
    }
    
    # Create CSV for all grouped chromatograms
    write.csv(gchroms, file=paste0("../data/chromatograms/","all_grouped.csv"), row.names=FALSE)
}



# ---------- Plot grouped species chromatograms ----------
f.ms1_chromatograms_group_plot <- function(mzml_files, mzml_names) {
    # Define plot limits
    sp_xlim <- c(100,1020)
    sp_ylim <- c(0,10)
    
    # Species full names
    species_full_names <- c(
        "Brachythecium rutabulum",
        "Calliergonella cuspidata",
        "Fissidens taxifolius",
        "Grimmia pulvinata",
        "Hypnum cupressiforme",
        "Marchantia polymorpha",
        "Plagiomnium undulatum",
        "Polytrichum strictum",
        "Rhytidiadelphus squarrosus"
    )
    
    # Abbreviated species names
    species_abbr_names <- c(
        "B. rutabulum",
        "C. cuspidata",
        "F. taxifolius",
        "G. pulvinata",
        "H. cupressiforme",
        "M. polymorpha",
        "P. undulatum",
        "P. strictum",
        "R. squarrosus"
    )
    
    # Plot
    for (sp in levels(species)) {
        pdf(paste0("../data/chromatograms/_",sp,".pdf"), encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
        plot(0, 0, type="n", xlim=sp_xlim, ylim=sp_ylim, main=paste("Grouped chromatograms of",species_full_names[which(species_names==sp)]),
             xlab="rt", ylab="TIC")
        for (i in which(species==sp)) {
            chroma <- xcmsRaw(mzml_files[i])
            x <- chroma@scantime
            y <- scale(chroma@tic, center=FALSE)
            lines(x, y, col=seasons_colors[seasons[i]], xlim=sp_xlim, ylim=sp_ylim)
        }
        legend("topleft", bty="n", pch=16, col=seasons_colors, pt.cex=0.8, cex=0.8, y.intersp=0.7, text.width=0.5, legend=seasons_names)
        dev.off()
    }
}



# ---------- Preprocess binary features list ----------
f.ms1_preprocess_features <- function(sclass) {
    # Get Reduced Peaklist
    xcam_report <- getReducedPeaklist(peak_xcam, method="median", default.adduct.info="maxint", cleanup=FALSE)
    
    # Diff report
    diff_list <- xcam_report
    diff_list <- diff_list[order(diff_list$pcgroup, decreasing=FALSE),]
    diff_list$pcgroup <- paste("pos_", diff_list$pcgroup, sep="")
    
    # Create peak list
    peak_list <- xcam_report
    peak_list <- peak_list[order(peak_list$pcgroup, decreasing=FALSE),]
    pcgroup <- peak_list$pcgroup
    peak_list <- peak_list[, which(colnames(peak_list) == mzml_names[1]) : which(colnames(peak_list) == mzml_names[length(mzml_names)])]
    
    # Create feature list
    feat_list <- peak_list
    rownames(feat_list) <- paste("pos_", unique(pcgroup), sep="")
    
    # Cleanup values in peak list: Remove NAs, negative abundances, constant features (rows)
    feat_list[is.na(feat_list)] <- 0
    feat_list[feat_list < 0] <- 0
    #feat_list <- feat_list[!apply(feat_list, MARGIN=1, function(x) max(x,na.rm=TRUE) == min(x,na.rm=TRUE)),]
    
    # Create single 0/1 matrix
    bina_list <- peak_list
    bina_list[is.na(bina_list)] <- 0
    bina_list[bina_list != 0] <- 1
    rownames(bina_list) <- paste("pos_", unique(pcgroup), sep="")
    
    # Only unique compounds in one group and not the others
    if (all(as.character(sclass) == as.character(species))) { #species
        uniq_list <- apply(X=bina_list, MARGIN=1,
                           FUN=function(x) { if (length(unique(species[grepl("1", x)])) == 1) x else rep(0, length(x)) } )
    } else { #seasons
        uniq_list <- apply(X=bina_list, MARGIN=1,
                           FUN=function(x) { if (length(unique(seasons[grepl("1", x)])) == 1) x else rep(0, length(x)) } )
    }
    uniq_list <- t(uniq_list)
    colnames(uniq_list) <- colnames(bina_list)
    
    # Return global variables
    diff_list <<- t(diff_list)
    peak_list <<- t(peak_list)
    feat_list <<- t(feat_list)
    bina_list <<- t(bina_list)
    uniq_list <<- t(uniq_list)
}



# ---------- Export ReducedPeaklist used for statistics as MAF ----------
f.export_maf <- function(xcam, maf_filename) {
    # Export ReducedPeaklist
    xcam_report <- getPeaklist(xcam)
    l <- nrow(xcam_report)
    
    # These columns are defined by MetaboLights mzTab
    maf <- apply(data.frame(database_identifier = character(l),
                            chemical_formula = character(l),
                            smiles = character(l),
                            inchi = character(l),
                            metabolite_identification = character(l),
                            mass_to_charge = xcam_report$mz,
                            fragmentation = character(l),
                            modifications = character(l),
                            charge = character(l),
                            retention_time = xcam_report$rt,
                            taxid = character(l),
                            species = character(l),
                            database = character(l),
                            database_version = character(l),
                            reliability = character(l),
                            uri = character(l),
                            search_engine = character(l),
                            search_engine_score = character(l),
                            smallmolecule_abundance_sub = character(l),
                            smallmolecule_abundance_stdev_sub = character(l),
                            smallmolecule_abundance_std_error_sub = character(l),
                            xcam_report,
                            stringsAsFactors=FALSE),
                 2, as.character)
    
    # Export MAF
    write.table(maf, file=maf_filename, row.names=FALSE, col.names=colnames(maf), quote=TRUE, sep="\t", na="\"\"")
    
    # Return variables (cols: features, rows: samples)
    #peak_maf <<- maf
}



# ---------- Update existing MAF with annotated compounds ----------
f.annotate_maf <- function(csv_compounds, maf_input, maf_output) {
    # Annotate identified compounds
    max.mz.range <- 0.5
    max.rt.range <- 20
    
    # Import new MAF
    maf_in <- read.table(file=maf_input, quote="\"", sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
    maf_in[is.na(maf_in)] <- as.character("")
    maf_out <- maf_in
    
    # Read list of annotated compounds from in-house library
    csv_in <- read.table(file=paste0(mzml_dir, "/../", csv_compounds), quote="\"", sep=";", dec=",", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
    csv_in[is.na(csv_in)] <- as.character("")
    
    # Extract InChI keys 
    inchis <- data.frame(InChI=csv_in$inchi[which(csv_in$inchi != "")])
    inchis_filename <- "ms1.inchis.tsv"
    write.table(inchis, file=inchis_filename, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", na="\"\"")
    
    # Convert InChI to InChIKey
    inchikeys <- apply(X=inchis, MARGIN=1, FUN=function(x) { x <- cs_inchi_inchikey(inchi=x, verbose=FALSE) })
    inchikeys <- data.frame(InChIKeys=inchikeys)
    inchikeys_filename <- "ms1.inchikeys.tsv"
    write.table(inchikeys, file=inchikeys_filename, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", na="\"\"")
    
    # Convert InChIKeys to ChEBI with: http://cts.fiehnlab.ucdavis.edu/batch
    chebi_filename <- "ms1.chebiids.csv"
    chebis <- read.table(file=chebi_filename, quote="\"", sep=",", na.strings="", header=TRUE, stringsAsFactors=TRUE)
    chebis[which(chebis$ChEBI == "No result"),"ChEBI"] <- NA
    
    # Integrate ChEBIs to database_identifier field
    csv_in[which(csv_in$inchi != ""), "database_identifier"] <- as.character(chebis$ChEBI)
    csv_in[which(is.na(csv_in$database_identifier)), "database_identifier"] <- ""
    
    # Calculate fixed-digit character representation for comparison with annotated compound list
    val_in <- cbind(mz=as.numeric(formatC(maf_in[,"mass_to_charge"], format="f", digits=4, flag="#")),
                    rt=as.numeric(formatC(maf_in[,"retention_time"], format="f", digits=1, flag="#")))
    com_in <- cbind(mz=as.numeric(formatC(csv_in[,"mass_to_charge"], format="f", digits=4, flag="#")),
                    rt=as.numeric(formatC(csv_in[,"RT.sec"], format="f", digits=1, flag="#")))
    
    com_list <- apply(X=com_in, MARGIN=1, FUN=function(x) { which ( ((x["mz"] >= val_in[,"mz"]-max.mz.range) & (x["mz"] <= val_in[,"mz"]+max.mz.range)) & ((x["rt"] >= val_in[,"rt"]-max.rt.range) & (x["rt"] <= val_in[,"rt"]+max.rt.range)) ) } )
    val_index <- NULL; for (i in 1:length(com_list)) { if (length(com_list[[i]])>0) val_index <- c(val_index, com_list[[i]][1]) }
    com_index <- NULL; for (i in 1:length(com_list)) { if (length(com_list[[i]])>0) com_index <- c(com_index, i) }
    
    # Annotation
    maf_out[val_index, "chemical_formula"] <- as.character(csv_in[com_index, "chemical_formula"])
    maf_out[val_index, "metabolite_identification"] <- as.character(csv_in[com_index, "metabolite_identification"])
    maf_out[val_index, "inchi"] <- as.character(csv_in[com_index, "inchi"])
    maf_out[val_index, "database_identifier"] <- as.character(csv_in[com_index, "database_identifier"])
    
    # Export MAF
    write.table(maf_out, file=maf_output, row.names=FALSE, col.names=colnames(maf_out), quote=TRUE, sep="\t", na="\"\"")
    
    
    # Annotate compound classes
    # Re-Import MAF
    maf_in <- maf_out
    
    # Read classified samples
    classifiers <- list()
    for (i in 1:length(mzml_names)) {
        classifiers[[i]] <- read.table(file=paste(mzml_dir,"/../classifier/",mzml_names[i],".tsv", sep=""), header=TRUE, sep="\t", quote="\"", fill=FALSE, dec=".", stringsAsFactors=FALSE)
    }
    
    # Read CHEMONT ontology
    obo <- get_ontology(file=paste0(mzml_dir, "/../ChemOnt_2_1.obo"), extract_tags="minimal")
    
    # Query CHEMONTIDs for determined classes
    classes_id <- NULL
    for (i in 1:length(classes)) {
        classes_id <- c(classes_id, as.character(obo$id[which(as.character(obo$name) == gsub(x=classes[i], pattern=".*; ", replacement=""))]))
    }
    
    # Process annotated samples
    for (i in 1:length(mzml_names)) {
        for (j in 1:nrow(classifiers[[i]])) {
            # Compound info
            mz <- classifiers[[i]][j,"m.z"]
            rt <- classifiers[[i]][j,"RT"]
            cl <- gsub(x=classifiers[[i]][j,"Annotation..putative."], pattern=".*; ", replacement="")
            id <- as.character(obo$id[which(as.character(obo$name) == as.character(cl))])
            
            # Add info in MAF
            li <- which( (maf_out$mass_to_charge > mz - max.mz.range) & (maf_out$mass_to_charge < mz + max.mz.range) &
                         (maf_out$retention_time > rt - max.rt.range) & (maf_out$retention_time < rt + max.rt.range) )
            if (length(li) == 0) {
                print(paste("No candidates in peak list found for",classifiers[[i]][j,"Label"],"in",classifiers[[i]][j,"Metabolite.name"]))
            } else {
                # Take the closest mz
                if (length(li) > 1) {
                    li <- li[which.min(abs(maf_out[li,"mass_to_charge"]-mz))]
                }
                # Write database_identifier
                if (nchar(maf_in[li, "database_identifier"]) > 0)
                    maf_out[li, "database_identifier"] <- paste(id, maf_in[li, "database_identifier"], sep=" | ")
                else
                    maf_out[li, "database_identifier"] <- id
                # Write metabolite_identification
                if (nchar(maf_in[li, "metabolite_identification"]) > 0)
                    maf_out[li, "metabolite_identification"] <- paste(cl, maf_in[li, "metabolite_identification"], sep=" | ")
                else
                    maf_out[li, "metabolite_identification"] <- cl
            }
        }
    }
    
    # Export MAF
    write.table(maf_out, file=maf_output, row.names=FALSE, col.names=colnames(maf_out), quote=TRUE, sep="\t", na="\"\"")
}


