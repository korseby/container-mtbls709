# ---------- Preparations ----------
# MS2 Matrices
msms_spectra <- NULL
msms_merged <- NULL



# ---------- Find MS/MS spectra ----------
f.ms2_find_spectra <- function() {
    # Method to use
    MergeSpectra <- FALSE
    
    # Find all MS/MS spectra in each sample
    msms_spectra <- list()
    msms_spectra <- foreach(i=1:length(mzml_files)) %dopar% {
        readMSData(files=mzml_files[i], msLevel=2, verbose=TRUE)
    }
    
    # Merge MS/MS spectra in each sample
    msms_merged <- list()
    if (MergeSpectra == TRUE) {
        msms_merged <- foreach(i=1:length(mzml_files)) %dopar% {
            merge.spectra(msms_spectra[[i]], mzabs, mzppm, rtabs, max.rt.range, max.mz.range, min.rt, max.rt, min.mz, max.mz, msms.intensity.threshold)
        }
    } else {
        msms_merged <- foreach(i=1:length(mzml_files)) %dopar% {
            collect.spectra.lists(msms_spectra[[i]], mzabs, mzppm, rtabs)
        }
    }
    
    # Return variables
    msms_spectra <<- msms_spectra
    msms_merged <<- msms_merged
}



# ---------- Plot MS/MS spectra ----------
f.ms2_plot_spectra <- function() {
    # Plot number of MS/MS spectra in each sample
    msms_num_spectra <- data.frame(species=species, sample=mzml_names, nspectra=0)
    for (i in 1:length(msms_merged))
        msms_num_spectra[i,"nspectra"] <- length(msms_merged[[i]])
    
    # For species
    boxplot(nspectra ~ species, data=msms_num_spectra, col=species_colors, names=NA,
            main="Number of acquired Auto-MS spectra per species", xlab="Species", ylab="Number of spectra")
    text(1:length(species_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=species_names, xpd=TRUE, cex=0.9)
    
    # For seasons
    # R-bug: prevent sorting when using formula
    model_boxplot <- data.frame(summer=msms_num_spectra$nspectra[seasons=="summer"],
                                autumn=msms_num_spectra$nspectra[seasons=="autumn"],
                                winter=msms_num_spectra$nspectra[seasons=="winter"],
                                spring=msms_num_spectra$nspectra[seasons=="spring"])
    boxplot(model_boxplot, data=msms_num_spectra, col=seasons_colors,
            main="Number of acquired Auto-MS spectra in the seasons", xlab="Seasons", ylab="Number of spectra")

    # For each profile
    msms_plot <- list()
    msms_plot <- foreach(i=1:length(mzml_files)) %dopar% {
        temp <- as.data.frame(matrix(unlist(lapply(msms_merged[[i]], FUN=function(x) { x <- data.frame(mz=x@precursorMz, rt=x@rt) } )), ncol=2))
        colnames(temp) <- c("mz", "rt")
        temp <- temp[order(temp$rt), ]
        return(temp)
    }
    
    pdf("ms_merged_spectra.pdf", encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    for (i in c(1:length(mzml_files))) {
        x_min <- floor(min(msms_plot[[i]][,"rt"]))
        x_max <- ceiling(max(msms_plot[[i]][,"rt"]))
        y_min <- floor(min(msms_plot[[i]][,"mz"]))
        y_max <- ceiling(max(msms_plot[[i]][,"mz"]))
        plot(x=msms_plot[[i]][,"rt"], y=msms_plot[[i]][,"mz"],
             xlim=c(x_min, x_max), ylim=c(y_min, y_max),
             main=mzml_names[i], xlab="retention time [rt]", ylab="merged m/z [mz]",
             pch=16, col=rainbow(n=nrow(msms_plot[[i]])))
    }
    dev.off()
}



# ---------- Write out MS/MS spectra in MSP text format ----------
f.ms2_create_msp <- function() {
    for (i in 1:length(msms_merged)) {
        msp_text <- NULL
        for (j in 1:length(msms_merged[[i]])) {
            NAME <- paste(mzml_names[i], ":", j, sep='')
            AlignmentID <- j
            RETENTIONTIME <- msms_merged[[i]][[j]]@rt
            PRECURSORMZ <- msms_merged[[i]][[j]]@precursorMz
            METABOLITENAME <- "Unknown"
            ADDUCTIONNAME <- "[M+H]+"
            NumPeaks <- msms_merged[[i]][[j]]@peaksCount
            Peaks <- NULL
            for (k in 1:length(msms_merged[[i]][[j]]@mz)) {
                Peaks <- c(Peaks, paste(msms_merged[[i]][[j]]@mz[k], msms_merged[[i]][[j]]@intensity[k], sep="\t") )
            }
            msp_text <- c(msp_text, paste("NAME:",NAME))
            msp_text <- c(msp_text, paste("AlignmentID:",AlignmentID))
            msp_text <- c(msp_text, paste("RETENTIONTIME:",RETENTIONTIME))
            msp_text <- c(msp_text, paste("PRECURSORMZ:",PRECURSORMZ))
            msp_text <- c(msp_text, paste("METABOLITENAME:",METABOLITENAME))
            msp_text <- c(msp_text, paste("ADDUCTIONNAME:",ADDUCTIONNAME))
            msp_text <- c(msp_text, paste("NumPeaks:",NumPeaks))
            msp_text <- c(msp_text, Peaks)
            msp_text <- c(msp_text, "")
        }
        cat(msp_text, file=paste(mzml_dir,"/../msp/",mzml_names[i],".msp",sep=""), sep="\n")
    }
}



# ---------- Tukey-Test ----------
tukey.test <- function(response, term) {
    if (sum(response) <= 0)
        return(data.frame("tukey_groups"=rep("z",length(term))))
    model_anova <- aov(formula(response ~ term))
    model_mc <- multcomp::glht(model_anova, multcomp::mcp(term="Tukey"))
    model_cld <- multcomp::cld(summary(model_mc), decreasing=TRUE)
    model_tukey <- data.frame("tukey_groups"=model_cld$mcletters$Letters)
    return(model_tukey)
}



# ---------- P-Values ----------
print_p.values <- function(p.values) {
    p.values[1] <- 1
    p.values[p.values < 0.001] <- '***'
    p.values[(p.values >= 0.001 & p.values < 0.01)] <- '**'
    p.values[(p.values >= 0.01 & p.values < 0.05)] <- '*'
    p.values[(p.values >= 0.05 & p.values < 0.1)] <- '.'
    p.values[p.values >= 0.1] <- ' '
    return(p.values)
}



# ---------- Manually annotate some MS/MS features ----------
f.ms2_classify_features <- function() {
    #classifier_name <- "2018-05-30_MSMS_pos_28558_MoNA_IPB_GNPS_Classifier"
    classifier_name <- "2019-04-17_MSMS_pos_21908_MoNA_Classifier"
    #classifier_name <- "2019-04-27_MSMS_pos_23938_410767_MoNA_WeizMASS_NIST2017_IPB_GNPS_PlaSMA_LignSesq_Classifier"
    #classifier_name <- "2019-06-16_MSMS_HR_pos_free__c_13314__s_71354__mFam_MoNA_GNPS_PlaSMA_Classifier"

    # Apply classifier on sample spectra
    classifiers <- list()
    classifiers <- foreach(i=1 : length(mzml_names)) %dopar% {
        #print(paste("Processing file",mzml_names[i],"..."))
        applyClassifierMs2(classifierFile = paste(mzml_dir, "/../classifier/", classifier_name, ".RData", sep=""),
                           propertiesFile = paste(mzml_dir, "/../classifier/", classifier_name, ".txt", sep=""),
                           fileMs1Path = NULL,
                           fileMs2Path = paste(mzml_dir, "/../msp/", mzml_names[i], ".msp", sep=""),
                           fileClasses = paste(mzml_dir, "/../classifier/classes.txt", sep=""),
                           minimumIntensityOfMaximalMS2peak = msms.intensity.threshold, #100
                           minimumProportionOfMS2peaks = 0.05,
                           mzDeviationAbsolute_grouping = mzabs, #0.01
                           mzDeviationInPPM_grouping = mzppm ) #10
    }
    
    # Save classifiers
    for (i in 1:length(mzml_names)) {
        write.table(x=classifiers[[i]], file=paste(mzml_dir,"/../classifier/",mzml_names[i],".tsv", sep=""), sep="\t", quote=TRUE, row.names=FALSE, dec=".")
    }
    
    # Load classifier
    load(paste(mzml_dir, "/../classifier/", classifier_name, ".RData", sep=""))
    
    # Read classified samples
    classifiers <- list()
    for (i in 1:length(mzml_names)) {
        classifiers[[i]] <- read.table(file=paste(mzml_dir,"/../classifier/",mzml_names[i],".tsv", sep=""), header=TRUE, sep="\t", quote="\"", fill=FALSE, dec=".", stringsAsFactors=FALSE)
    }
    
    # Diversity of classes in all samples
    div_classes_samples <- data.frame()
    for (i in 1:length(mzml_names)) {
        obj <- data.frame(sample=mzml_names[[i]], classes=unique(classifiers[[i]][,"Annotation..putative."]), frequency=0)
        for (j in 1:nrow(obj)) obj[j,"frequency"] <- length(which(obj$classes[j] == classifiers[[i]][,"Annotation..putative."]))
        div_classes_samples <- rbind(div_classes_samples, obj)
    }
    
    # Diversity of classes per species
    div_classes_species <- data.frame()
    for (i in 1:length(species_names)) {
        k <- c(which(species==species_names[i]))
        obj <- data.frame(species=species_names[i], classes=unique(div_classes_samples[which(div_classes_samples$sample %in% mzml_names[k]), "classes"]), frequency=0)
        for (j in 1:nrow(obj)) obj[j,"frequency"] <- length(which(obj$classes[j] == div_classes_samples[which(div_classes_samples$sample %in% mzml_names[k]), "classes"]))
        obj <- obj[order(obj$frequency, decreasing=TRUE),]
        div_classes_species <- rbind(div_classes_species, obj)
    }
    
    # Diversity of classes per seasons
    div_classes_seasons <- data.frame()
    for (i in 1:length(seasons_names)) {
        k <- c(which(seasons==seasons_names[i]))
        obj <- data.frame(seasons=seasons_names[i], classes=unique(div_classes_samples[which(div_classes_samples$sample %in% mzml_names[k]), "classes"]), frequency=0)
        for (j in 1:nrow(obj)) obj[j,"frequency"] <- length(which(obj$classes[j] == div_classes_samples[which(div_classes_samples$sample %in% mzml_names[k]), "classes"]))
        obj <- obj[order(obj$frequency, decreasing=TRUE),]
        div_classes_seasons <- rbind(div_classes_seasons, obj)
    }
    
    # Diversity of classes
    obj <- data.frame()
    for (i in 1:length(mzml_names)) obj <- rbind(obj, classifiers[[i]])
    obj <- data.frame(classes=obj[,"Annotation..putative."], frequency=0)
    div_classes <- data.frame(classes=unique(obj[,"classes"]), frequency=0)
    for (j in 1:nrow(div_classes)) div_classes[j,"frequency"] <- length(which(div_classes$classes[j] == obj[,"classes"]))
    div_classes <- div_classes[order(div_classes$frequency, decreasing=TRUE),]
    
    # Plot most abundant classes
    pdf(file="classes_counts.pdf", encoding="ISOLatin1", pointsize=10, width=10, height=10, family="Helvetica")
    par(mfrow=c(1,1), mar=c(38,4,4,1), oma=c(0,0,0,0), cex.axis=0.8, cex=0.9)
    barplot(div_classes[1:nrow(div_classes),"frequency"], names.arg=div_classes[1:nrow(div_classes),"classes"], las=3, ylab="frequency", main="Most abundant compound classes")
    dev.off()
    
    # Determine how many spectra were classified
    spectra_number <- sum(unlist(lapply(X=as.numeric(1:length(mzml_names)), FUN=function(x) { x <- length(msms_merged[[x]]) } )))
    spectra_classified <- sum(div_classes_samples$frequency)
    print(paste("Number of merged spectra:", spectra_number))
    print(paste("Number of spectra classified:", spectra_classified))
    print(paste("Number of unclassified spectra:", spectra_number - spectra_classified))
    print(paste("Number of classes:", length(classes_order)))
    print(paste("Number of classes with entities:", length(classes)))
    print(paste("Number of classes without entities:", length(classes_order) - length(classes)))
    print("Classes with entities:")
    print(gsub(x=classes_order[which( (classes_order %in% classes))], pattern=".*; ", replacement=""))
    print("Classes without entities:")
    print(gsub(x=classes_order[which( ! (classes_order %in% classes))], pattern=".*; ", replacement=""))

    # Classes
    classes_order <- readLines(con=paste(mzml_dir, "/../classifier/classes.txt", sep=""))
    classes_order <- classes_order[which(grepl(pattern="^Organic compounds", x=classes_order))]
    classes_order <- gsub(x=classes_order, pattern=":", replacement="")
    classes_order <- gsub(x=classes_order, pattern="/", replacement="; ")
    
    classes <- div_classes$classes
    classes <- classes[which(grepl(pattern="^Organic compounds", x=classes))]
    classes <- gsub(x=classes, pattern=":", replacement="")
    classes <- gsub(x=classes, pattern="/", replacement="; ")
    
    classes <- classes_order[which(classes_order %in% classes)]
    
    # Classifiers
    classifiers_class <- get(load(paste(mzml_dir, "/../classifier/", classifier_name, ".RData", sep="")))
    
    # Area under Precision Recall Curve
    classifier_auc <- as.numeric(unlist(lapply(X=classes, FUN = function(x) { x <- classifiers_class[[x]]$AUC_PR } )))
    
    # True Positive Rate for False Positive Rate of 5 Percent
    classifier_fpr <- as.numeric(unlist(lapply(X=classes, FUN = function(x) { x <- classifiers_class[[x]]$TPR_for_FPR_of_5Percent } )))
    
    # Save table with AUC-PR and TPR-FPR rates
    #write.table(x=data.frame(compound_class=classes, AUC=classifier_auc, FPR=classifier_fpr),
    #            file=paste(mzml_dir,"/../classifier/","correctness",".csv", sep=""),
    #            sep=";", quote=TRUE, row.names=FALSE, dec=".")
    write.table(x=data.frame(compound_class=classes, "AUC-PR"=classifier_auc, "TPR-FPR"=classifier_fpr),
                file="table_1.csv",
                sep=";", quote=TRUE, row.names=FALSE, dec=".")
    
    ######################################
    
    # Count numbers of matched classes in each sample
    class_list <- data.frame(species=species)
    for (i in classes) {
        vec <- NULL
        for (j in 1:length(mzml_names)) {
            vec <- as.numeric(c(vec, sum(grepl(pattern=i, x=classifiers[[j]][,"Annotation..putative."]), na.rm=TRUE)))
        }
        class_list <- cbind(class_list, vec)
    }
    colnames(class_list) <- c("species", gsub(x=classes, pattern=".*\\; ", replacement=""))
    rownames(class_list) <- mzml_names
    
    # Sunburst plot of classes
    source("sunburst.r")
    numberOfSpectra <- as.numeric(unlist(apply(X=class_list[,2:ncol(class_list)], MARGIN=2, FUN = function(x) { sum(x) })))
    numberOfSpectra <- numberOfSpectra[which(numberOfSpectra > 0)]
    classifierClasses <- classes[which(numberOfSpectra > 0)]
    #pdf(file="classes_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
    pdf(file="fig_1.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
    par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=1, cex=1)
    sunBurstPlotFromSubstanceClasses(classifierClasses, numberOfSpectra, colorStart=0.5, colorAlpha=0.5)
    dev.off()
    write.csv(data.frame("Classifier classes"=classifierClasses, "Number of spectra"=numberOfSpectra), file="fig_1.csv", row.names=FALSE)
    
    heatmap.2(x=as.matrix(class_list[,2:ncol(class_list)]), cexRow=0.3, cexCol=0.4,
              Rowv=as.dendrogram(hclust(dist(class_list[,2:ncol(class_list)]))), offsetRow=0,
              Colv=as.dendrogram(hclust(dist(t(class_list[,2:ncol(class_list)])))), offsetCol=0,
              col=colorRampPalette(c('lightgrey','white','darkblue'), alpha=0.1, bias=8)(256),
              trace="none", margins=c(7.2,6),
              key=TRUE, key.title="Color key")
    
    # Sum of matched classes in species
    model_class <- data.frame(matrix(0, ncol=length(classes), nrow=length(species_names)))
    colnames(model_class) <- gsub(x=classes, pattern=".*\\; ", replacement="")
    rownames(model_class) <- species_names
    for (i in 1:length(classes)) {
        model_class[,i] <- as.numeric(lapply(X=species_names, FUN=function(x) { sum(as.numeric(class_list[which(class_list[,1]==x), i+1])) } ))
    }
    
    # Boxplots
    #pdf(file="species_classes_boxplots.pdf", encoding="ISOLatin1", pointsize=15, width=5*4, height=4*6, family="Helvetica")
    pdf(file="fig_s4.pdf", encoding="ISOLatin1", pointsize=15, width=5*4, height=4*6, family="Helvetica")
    par(mfrow=c(6,4), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=1.0, cex=0.9)
    for (i in 1:length(classes)) {
        if (sum(class_list[,i+1]) > 0) {
            boxplot(class_list[,i+1] ~ species, col=species_colors, names=NA, main=colnames(class_list)[i+1], cex.main=0.9, xlab="Species", ylab="Chemical richness")
            text(1:length(species_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/10, srt=-22.5, adj=0.5, labels=species_names, xpd=TRUE, cex=0.9)
            div_tukey <- tukey.test(response=class_list[,i+1], term=species)
            text(1:length(species_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/15, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
            mtext(paste0("(",letters[i],")"), outer=TRUE, adj=0.01 + 0.2535*((i-1) %% 4), padj=2 + ((ceil(i/4) - 1) * 22.25), font=2, cex=1.2)
        }
    }
    dev.off()
    write.csv(class_list, file="fig_s4.csv", row.names=FALSE)
    
    #pdf(file="species_classes_boxplots_chosen.pdf", encoding="ISOLatin1", pointsize=12, width=5, height=4, family="Helvetica")
    pdf(file="fig_3.pdf", encoding="ISOLatin1", pointsize=12, width=5, height=4, family="Helvetica")
    for (i in 1:length(classes)) {
        if (sum(class_list[,i+1]) > 0) {
            if ( (colnames(class_list)[i+1] == "Sesquiterpenoids") | 
                 (colnames(class_list)[i+1] == "Flavonoids") |
                 (colnames(class_list)[i+1] == "Glycosyl compounds") |
                 (colnames(class_list)[i+1] == "Anthocyanins") |
                 (colnames(class_list)[i+1] == "Lignans, neolignans and related compounds") | 
                 (colnames(class_list)[i+1] == "Stilbenes") ) {
                boxplot(class_list[,i+1] ~ species, col=species_colors, names=NA, main=colnames(class_list)[i+1], cex.main=0.9, xlab="Species", ylab="Chemical richness")
                text(1:length(species_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/10, srt=-22.5, adj=0.5, labels=species_names, xpd=TRUE, cex=0.9)
                div_tukey <- tukey.test(response=class_list[,i+1], term=species)
                text(1:length(species_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/15, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
                if (colnames(class_list)[i+1] == "Sesquiterpenoids") mtext("(a)", outer=TRUE, adj=0.01, line=-1.4, font=2, cex=1.4)
                if (colnames(class_list)[i+1] == "Flavonoids") mtext("(b)", outer=TRUE, adj=0.01, line=-1.4, font=2, cex=1.4)
                if (colnames(class_list)[i+1] == "Glycosyl compounds") mtext("(c)", outer=TRUE, adj=0.01, line=-1.4, font=2, cex=1.4)
                if (colnames(class_list)[i+1] == "Anthocyanins") mtext("(d)", outer=TRUE, adj=0.01, line=-1.4, font=2, cex=1.4)
                if (colnames(class_list)[i+1] == "Lignans, neolignans and related compounds") mtext("(e)", outer=TRUE, adj=0.01, line=-1.4, font=2, cex=1.4)
                if (colnames(class_list)[i+1] == "Stilbenes") mtext("(f)", outer=TRUE, adj=0.01, line=-1.4, font=2, cex=1.4)
            }
        }
    }
    dev.off()
    write.csv(class_list[,c("species","Sesquiterpenoids","Flavonoids","Glycosyl compounds","Anthocyanins","Lignans, neolignans and related compounds","Stilbenes")], file="fig_3.csv", row.names=FALSE)
    
    # Heatmap of classes per species
    model_heat <- as.matrix(log(model_class))
    model_heat[which(is.infinite(model_heat))] <- 0
    
    class_rowclust <- as.dendrogram(hclust(vegdist(model_class, method="bray"), method="single"))
    class_colclust <- as.dendrogram(hclust(dist(t(model_class), method="maximum"), method="complete"))
    class_colclust <- reorder(as.dendrogram(class_colclust), seq(1:ncol(model_class)))
    #gsub(x=classes, pattern=".*\\; ", replacement="")
    #which(gsub(x=classes, pattern=".*\\; ", replacement="") %in% as.character(class_colclust$labels[class_colclust$order]))
    
    #pdf(file="species_classes_heatmap.pdf", encoding="ISOLatin1", pointsize=15, width=14, height=10, family="Helvetica")
    pdf(file="fig_2.pdf", encoding="ISOLatin1", pointsize=15, width=14, height=10, family="Helvetica")
    par(oma=c(0,0,0,15))
    heatmap.2(x=as.matrix(model_heat), cexRow=1.2, cexCol=1.2,
                         Rowv=class_rowclust, offsetRow=0, colRow=species_colors,
                         Colv=class_colclust, offsetCol=0,
                         #cellnote=as.matrix(round(ceiling(model_class/(length(mzml_names)/length(species_names))),0)), notecol="black", notecex=0.6,
                         #col=colorRampPalette(c('lightgrey','white','darkblue'), alpha=0.1, bias=8)(256),
                         col=colorRampPalette(c('white','lightgrey','darkgreen'), alpha=0.01, bias=1.6)(256),
                         trace="none", margins=c(21,5),
                         key=TRUE, key.title="Color key", density.info='density', denscol="black")
    mtext(text="(a)", side=2, adj=0, padj=-25, las=1, line=3, font=2, cex=1.2)
    
    par(new=TRUE, oma=c(12.5,40,6.9,0))
    plot(phylo_tree, type="phylogram", direction="leftwards", x.lim=c(0,12),
                    label.offset=0.5, use.edge.length=TRUE, show.tip.label=TRUE,
                    tip.color=species_colors[phylo_index], font=1, main="")
    mtext(text="(b)", side=2, adj=0, padj=-21, las=1, line=-1, font=2, cex=1.2)
    dev.off()
    write.csv(model_heat, file="fig_2.csv", row.names=TRUE)
    
    # Linear model for each class in the species, seasons: random factor
    model_class_lmer <- data.frame(matrix(0, ncol=length(classes), nrow=length(species_names)))
    colnames(model_class_lmer) <- gsub(x=classes, pattern=".*\\; ", replacement="")
    rownames(model_class_lmer) <- species_names

    model_class_lmer_p <- data.frame(matrix(0, ncol=length(classes), nrow=length(species_names)))
    colnames(model_class_lmer_p) <- gsub(x=classes, pattern=".*\\; ", replacement="")
    rownames(model_class_lmer_p) <- species_names
    
    # LMER for each class
    for (i in 1:length(classes)) {
        if (sum(class_list[,i+1]) > 0) {
            model_lmer <- lmer(class_list[,i+1] ~ species + (1|seasons))
            model_class_lmer[,i] <- coef(summary(model_lmer))[,5]
            model_class_lmer_p[,i] <- print_p.values(coef(summary(model_lmer))[,5])
        } else {
            model_class_lmer[,i] <- rep(1, length(species_names))
            model_class_lmer_p[,i] <- rep("0", length(species_names))
        }
    }

    # Save models of statistics
    write.table(x=cbind(species=species_names,model_class), file="species_class.tsv", sep="\t", quote=TRUE, row.names=FALSE, dec=".")
    write.table(x=cbind(species=species_names,model_class_lmer), file="species_lmer.tsv", sep="\t", quote=TRUE, row.names=FALSE, dec=".")
    write.table(x=cbind(species=species_names,model_class_lmer_p), file="species_lmer_p.tsv", sep="\t", quote=TRUE, row.names=FALSE, dec=".")
    
    ############################
    
    # Count numbers of matched classes in each sample
    class_list <- data.frame(seasons=seasons)
    for (i in classes) {
        vec <- NULL
        for (j in 1:length(mzml_names)) {
            vec <- as.numeric(c(vec, sum(grepl(pattern=i, x=classifiers[[j]][,"Annotation..putative."]), na.rm=TRUE)))
        }
        class_list <- cbind(class_list, vec)
    }
    colnames(class_list) <- c("seasons", gsub(x=classes, pattern=".*\\; ", replacement=""))
    rownames(class_list) <- mzml_names
    
    # Sum of matched classes in seasons
    model_class <- data.frame(matrix(0, ncol=length(classes), nrow=length(seasons_names)))
    colnames(model_class) <- gsub(x=classes, pattern=".*\\; ", replacement="")
    rownames(model_class) <- seasons_names
    for (i in 1:length(classes)) {
        model_class[,i] <- as.numeric(lapply(X=seasons_names, FUN=function(x) { sum(as.numeric(class_list[which(class_list[,1]==x), i+1])) } ))
    }
    
    # Boxplots
    #pdf(file="seasons_classes_boxplots.pdf", encoding="ISOLatin1", pointsize=15, width=4*4, height=4*6, family="Helvetica")
    pdf(file="fig_s5.pdf", encoding="ISOLatin1", pointsize=15, width=4*4, height=4*6, family="Helvetica")
    par(mfrow=c(6,4), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=1.0, cex=0.9)
    for (i in 1:length(classes)) {
        if (sum(class_list[,i+1]) > 0) {
            model_boxplot <- data.frame(summer=class_list[,i+1][seasons=="summer"],
                                        autumn=class_list[,i+1][seasons=="autumn"],
                                        winter=class_list[,i+1][seasons=="winter"],
                                        spring=class_list[,i+1][seasons=="spring"])
            boxplot(x=model_boxplot, col=seasons_colors, names=NA, main=colnames(class_list)[i+1], cex.main=0.9, xlab="Seasons", ylab="Chemical richness")
            text(1:length(seasons_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/10, srt=-22.5, adj=0.5, labels=seasons_names, xpd=TRUE, cex=0.9)
            div_tukey <- tukey.test(response=as.numeric(apply(X=model_boxplot, MARGIN=1, FUN=function(x) { x })), term=seasons)
            text(1:length(seasons_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/20, adj=0.5, labels=div_tukey[match(seasons_names, rownames(div_tukey)),1], xpd=TRUE, cex=0.8)
            mtext(paste0("(",letters[i],")"), outer=TRUE, adj=0.01 + 0.256*((i-1) %% 4), padj=2 + ((ceil(i/4) - 1) * 22.25), font=2, cex=1.2)
        }
    }
    dev.off()
    write.csv(class_list, file="fig_s5.csv", row.names=FALSE)
    
    #pdf(file="seasons_classes_boxplots_chosen.pdf", encoding="ISOLatin1", pointsize=12, width=4, height=4, family="Helvetica")
    pdf(file="fig_5.pdf", encoding="ISOLatin1", pointsize=12, width=4, height=4, family="Helvetica")
    for (i in 1:length(classes)) {
        if (sum(class_list[,i+1]) > 0) {
            if ( (colnames(class_list)[i+1] == "Fatty acids and conjugates") | 
                 (colnames(class_list)[i+1] == "Carbohydrates and carbohydrate conjugates") |
                 (colnames(class_list)[i+1] == "Phenylpropanoids and polyketides") |
                 (colnames(class_list)[i+1] == "Flavonoids") |
                 (colnames(class_list)[i+1] == "Anthocyanins") ) {
                model_boxplot <- data.frame(summer=class_list[,i+1][seasons=="summer"],
                                            autumn=class_list[,i+1][seasons=="autumn"],
                                            winter=class_list[,i+1][seasons=="winter"],
                                            spring=class_list[,i+1][seasons=="spring"])
                boxplot(x=model_boxplot, col=seasons_colors, names=NA, main=colnames(class_list)[i+1], cex.main=0.9, xlab="Seasons", ylab="Chemical richness")
                text(1:length(seasons_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/10, srt=-22.5, adj=0.5, labels=seasons_names, xpd=TRUE, cex=0.9)
                div_tukey <- tukey.test(response=as.numeric(apply(X=model_boxplot, MARGIN=1, FUN=function(x) { x })), term=seasons)
                text(1:length(seasons_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/20, adj=0.5, labels=div_tukey[match(seasons_names, rownames(div_tukey)),1], xpd=TRUE, cex=0.8)
                if (colnames(class_list)[i+1] == "Fatty acids and conjugates") mtext("(a)", outer=TRUE, adj=0.01, line=-1.4, font=2, cex=1.4)
                if (colnames(class_list)[i+1] == "Carbohydrates and carbohydrate conjugates") mtext("(b)", outer=TRUE, adj=0.01, line=-1.4, font=2, cex=1.4)
                if (colnames(class_list)[i+1] == "Phenylpropanoids and polyketides") mtext("(c)", outer=TRUE, adj=0.01, line=-1.4, font=2, cex=1.4)
                if (colnames(class_list)[i+1] == "Flavonoids") mtext("(d)", outer=TRUE, adj=0.01, line=-1.4, font=2, cex=1.4)
                if (colnames(class_list)[i+1] == "Anthocyanins") mtext("(e)", outer=TRUE, adj=0.01, line=-1.4, font=2, cex=1.4)
            }
        }
    }
    dev.off()
    write.csv(class_list[,c("seasons","Fatty acids and conjugates","Carbohydrates and carbohydrate conjugates","Phenylpropanoids and polyketides","Flavonoids","Anthocyanins")], file="fig_5.csv", row.names=FALSE)

    # Heatmap of classes per seasons
    model_heat <- as.matrix(log(model_class))
    model_heat[which(is.infinite(model_heat))] <- 0
    
    #class_rowclust <- hclust(dist(model_class, method="maximum"), method="ward.D2")
    class_rowclust <- as.dendrogram(hclust(vegdist(model_class, method="cao"), method="ward.D2"))
    class_rowclust <- reorder(as.dendrogram(class_rowclust), seq(1:nrow(model_class)))
    class_colclust <- as.dendrogram(hclust(dist(t(model_class), method="maximum"), method="complete"))
    class_colclust <- reorder(as.dendrogram(class_colclust), seq(1:ncol(model_class)))
    
    #pdf(file="seasons_classes_heatmap.pdf", encoding="ISOLatin1", pointsize=15, width=10, height=8, family="Helvetica")
    pdf(file="fig_4.pdf", encoding="ISOLatin1", pointsize=15, width=10, height=8, family="Helvetica")
    heatmap.2(x=as.matrix(model_heat), cexRow=1.2, cexCol=1.2,
              Rowv=class_rowclust, offsetRow=0,
              Colv=class_colclust, offsetCol=0,
              #cellnote=as.matrix(model_class), notecol="black", notecex=0.6,
              col=colorRampPalette(c('white','lightgrey','darkgreen'), alpha=0.01, bias=1.6)(256),
              trace="none", margins=c(21,5),
              key=TRUE, key.title="Color key", density.info='density', denscol="black")
    dev.off()
    write.csv(model_heat, file="fig_4.csv", row.names=TRUE)
    
    # Linear model for each class in the seasons, seasons: random factor
    model_class_lmer <- data.frame(matrix(0, ncol=length(classes), nrow=length(seasons_names)))
    colnames(model_class_lmer) <- gsub(x=classes, pattern=".*\\; ", replacement="")
    rownames(model_class_lmer) <- seasons_names
    
    model_class_lmer_p <- data.frame(matrix(0, ncol=length(classes), nrow=length(seasons_names)))
    colnames(model_class_lmer_p) <- gsub(x=classes, pattern=".*\\; ", replacement="")
    rownames(model_class_lmer_p) <- seasons_names
    
    # LMER for each class
    for (i in 1:length(classes)) {
        if (sum(class_list[,i+1]) > 0) {
            model_lmer <- lmer(class_list[,i+1] ~ seasons + (1|species))
            model_class_lmer[,i] <- coef(summary(model_lmer))[,5]
            model_class_lmer_p[,i] <- print_p.values(coef(summary(model_lmer))[,5])
        } else {
            model_class_lmer[,i] <- rep(1, length(seasons_names))
            model_class_lmer_p[,i] <- rep("0", length(seasons_names))
        }
    }
    
    # Save models of statistics
    write.table(x=cbind(seasons=seasons_names,model_class), file="seasons_class.tsv", sep="\t", quote=TRUE, row.names=FALSE, dec=".")
    write.table(x=cbind(seasons=seasons_names,model_class_lmer), file="seasons_lmer.tsv", sep="\t", quote=TRUE, row.names=FALSE, dec=".")
    write.table(x=cbind(seasons=seasons_names,model_class_lmer_p), file="seasons_lmer_p.tsv", sep="\t", quote=TRUE, row.names=FALSE, dec=".")
    
    ######################################
    
    # Count numbers of matched classes in each sample
    class_list <- data.frame(species=species, seasons=seasons)
    for (i in classes) {
        vec <- NULL
        for (j in 1:length(mzml_names)) {
            vec <- as.numeric(c(vec, sum(grepl(pattern=i, x=classifiers[[j]][,"Annotation..putative."]), na.rm=TRUE)))
        }
        class_list <- cbind(class_list, vec)
    }
    colnames(class_list) <- c("species", "seasons", gsub(x=classes, pattern=".*\\; ", replacement=""))
    rownames(class_list) <- mzml_names
    class_list <- class_list[,c(3:ncol(class_list))]
    for (i in 1:ncol(class_list)) { class_list[,i] <- as.numeric(class_list[,i]) }
    
    # Boxplots to show seasonal effect in each species
    #pdf(file="seasons_classes_boxplots_species.pdf", encoding="ISOLatin1", pointsize=15, width=4*4, height=4*6, family="Helvetica")
    pdf(file="fig_s6.pdf", encoding="ISOLatin1", pointsize=15, width=4*4, height=4*6, family="Helvetica")
    par(mfrow=c(6,4), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=1.0, cex=0.9)
    for (sp in species_names) {
        for (i in 1:length(classes)) {
            if (sum(class_list[,i]) > 0) {
                model_boxplot <- data.frame(summer=class_list[,i][seasons=="summer" & species==sp],
                                            autumn=class_list[,i][seasons=="autumn" & species==sp],
                                            winter=class_list[,i][seasons=="winter" & species==sp],
                                            spring=class_list[,i][seasons=="spring" & species==sp])
                boxplot(x=model_boxplot, col=seasons_colors, names=NA, main=colnames(class_list)[i], cex.main=0.9, xlab="Seasons", ylab="Chemical richness")
                text(1:length(seasons_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/10, srt=-22.5, adj=0.5, labels=seasons_names, xpd=TRUE, cex=0.9)
                div_tukey <- tukey.test(response=as.numeric(apply(X=model_boxplot, MARGIN=1, FUN=function(x) { x })), term=as.factor(rep(colnames(model_boxplot),nrow(model_boxplot))))
                text(1:length(seasons_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/20, adj=0.5, labels=div_tukey[match(seasons_names, rownames(div_tukey)),1], xpd=TRUE, cex=0.8)
            }
        }
        mtext(paste0("(",letters[which(species_names == sp)],") ",species_full_names[which(species_names == sp)]), outer=TRUE, adj=0.01, padj=1.6, font=4, cex=1.2)
    }
    
    sp <- c("Brarut", "Calcus", "Hypcup", "Rhysqu")
    for (i in 1:length(classes)) {
        if (sum(class_list[,i]) > 0) {
            model_boxplot <- data.frame(summer=class_list[,i][seasons=="summer" & (species==sp[1] | species==sp[2] | species==sp[3] | species==sp[4])],
                                        autumn=class_list[,i][seasons=="autumn" & (species==sp[1] | species==sp[2] | species==sp[3] | species==sp[4])],
                                        winter=class_list[,i][seasons=="winter" & (species==sp[1] | species==sp[2] | species==sp[3] | species==sp[4])],
                                        spring=class_list[,i][seasons=="spring" & (species==sp[1] | species==sp[2] | species==sp[3] | species==sp[4])])
            boxplot(x=model_boxplot, col=seasons_colors, names=NA, main=colnames(class_list)[i], cex.main=0.9, xlab="Seasons", ylab="Chemical richness")
            text(1:length(seasons_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/10, srt=-22.5, adj=0.5, labels=seasons_names, xpd=TRUE, cex=0.9)
            div_tukey <- tukey.test(response=as.numeric(apply(X=model_boxplot, MARGIN=1, FUN=function(x) { x })), term=as.factor(rep(colnames(model_boxplot),nrow(model_boxplot))))
            text(1:length(seasons_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/20, adj=0.5, labels=div_tukey[match(seasons_names, rownames(div_tukey)),1], xpd=TRUE, cex=0.8)
        }
    }
    mtext("(j) Pleurocarpous species", outer=TRUE, adj=0.01, padj=1.6, font=2, cex=1.2)
    
    sp <- c("Fistax", "Gripul", "Plaund", "Polstr")
    for (i in 1:length(classes)) {
        if (sum(class_list[,i]) > 0) {
            model_boxplot <- data.frame(summer=class_list[,i][seasons=="summer" & (species==sp[1] | species==sp[2] | species==sp[3] | species==sp[4])],
                                        autumn=class_list[,i][seasons=="autumn" & (species==sp[1] | species==sp[2] | species==sp[3] | species==sp[4])],
                                        winter=class_list[,i][seasons=="winter" & (species==sp[1] | species==sp[2] | species==sp[3] | species==sp[4])],
                                        spring=class_list[,i][seasons=="spring" & (species==sp[1] | species==sp[2] | species==sp[3] | species==sp[4])])
            boxplot(x=model_boxplot, col=seasons_colors, names=NA, main=colnames(class_list)[i], cex.main=0.9, xlab="Seasons", ylab="Chemical richness")
            text(1:length(seasons_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/10, srt=-22.5, adj=0.5, labels=seasons_names, xpd=TRUE, cex=0.9)
            div_tukey <- tukey.test(response=as.numeric(apply(X=model_boxplot, MARGIN=1, FUN=function(x) { x })), term=as.factor(rep(colnames(model_boxplot),nrow(model_boxplot))))
            text(1:length(seasons_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/20, adj=0.5, labels=div_tukey[match(seasons_names, rownames(div_tukey)),1], xpd=TRUE, cex=0.8)
        }
    }
    mtext("(k) Acrocarpous species", outer=TRUE, adj=0.01, padj=1.6, font=2, cex=1.2)
    dev.off()
    write.csv(cbind(species=species,season=seasons,class_list), file="fig_s6.csv", row.names=FALSE)
    
    # Boxplots to show seasonal effect in each species
    #pdf(file=paste0("seasons_classes_boxplots_species_chosen.pdf"), encoding="ISOLatin1", pointsize=12, width=4, height=4, family="Helvetica")
    pdf(file="fig_6.pdf", encoding="ISOLatin1", pointsize=12, width=4, height=4, family="Helvetica")
    model_csv <- data.frame(species=NULL, summer=NULL, autumn=NULL, winter=NULL, spring=NULL, compound_class=NULL)
    for (i in 1:10) {
        if (i == 1) { sp <- "Brarut"; cl <- "Lactones" }
        if (i == 2) { sp <- "Fistax"; cl <- "Anthocyanins" }
        if (i == 3) { sp <- "Gripul"; cl <- "Sesquiterpenoids" }
        if (i == 4) { sp <- "Marpol"; cl <- "Lignans, neolignans and related compounds" }
        if (i == 5) { sp <- "Marpol"; cl <- "Stilbenes" }
        if (i == 6) { sp <- "Marpol"; cl <- "Methoxyphenols" }
        if (i == 7) { sp <- "Marpol"; cl <- "Anthocyanins" }
        if (i == 8) { sp <- "Marpol"; cl <- "Lactones" }
        if (i == 9) { sp <- "Plaund"; cl <- "Lactones" }
        if (i == 10){ sp <- "Polstr"; cl <- "Carbohydrates and carbohydrate conjugates" }
        model_boxplot <- data.frame(summer=class_list[,cl][seasons=="summer" & species==sp],
                                    autumn=class_list[,cl][seasons=="autumn" & species==sp],
                                    winter=class_list[,cl][seasons=="winter" & species==sp],
                                    spring=class_list[,cl][seasons=="spring" & species==sp])
        model_csv <- rbind(model_csv, cbind(data.frame(species=c(sp,sp,sp)), model_boxplot, data.frame(compound_class=c(cl,cl,cl))))
        boxplot(x=model_boxplot, col=seasons_colors, names=NA, main=cl, cex.main=0.9, xlab="seasons", ylab="")
        mtext(paste0(species_full_names[which(species_names == sp)]), outer=FALSE, side=2, line=3, font=3, cex=1)
        text(1:length(seasons_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/10, srt=-22.5, adj=0.5, labels=seasons_names, xpd=TRUE, cex=0.9)
        div_tukey <- tukey.test(response=as.numeric(apply(X=model_boxplot, MARGIN=1, FUN=function(x) { x })), term=as.factor(rep(colnames(model_boxplot),nrow(model_boxplot))))
        text(1:length(seasons_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/20, adj=0.5, labels=div_tukey[match(seasons_names, rownames(div_tukey)),1], xpd=TRUE, cex=0.8)
        mtext(paste0("(",letters[i],")"), outer=TRUE, adj=0.01, line=-2, font=2, cex=1.4)
    }
    dev.off()
    write.csv(model_csv, file="fig_6.csv", row.names=FALSE)
}



# ---------- Variation partitioning ----------
f.ms2_classify_varpart <- function() {
    # Count numbers of matched classes in each sample
    class_list <- data.frame(species=species, seasons=seasons)
    for (i in classes) {
        vec <- NULL
        for (j in 1:length(mzml_names)) {
            vec <- as.numeric(c(vec, sum(grepl(pattern=i, x=classifiers[[j]][,"Annotation..putative."]), na.rm=TRUE)))
        }
        class_list <- cbind(class_list, vec)
    }
    colnames(class_list) <- c("species", "seasons", gsub(x=classes, pattern=".*\\; ", replacement=""))
    rownames(class_list) <- mzml_names
    class_list <- class_list[,c(3:ncol(class_list))]
    for (i in 1:ncol(class_list)) { class_list[,i] <- as.numeric(class_list[,i]) }
    
    # Variance partitioning (result depending on Curtis distance)
    vp_list <- class_list
    #vp_list[is.na(vp_list)] <- 0
    #vp_list[is.infinite(vp_list)] <- 0
    model_vp <- varpart(vp_list, ~ seasons, ~ species)
    
    # Plot results
    #pdf(file="ms2_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
    pdf(file="fig_s3b.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
    plot(model_vp, Xnames=c("seasons","species"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("darkgreen","darkblue"))
    dev.off()
    write.csv(vp_list, file="fig_s3b.csv", row.names=TRUE)
    
    
    
    # Trait of interest
    constrain <- as.data.frame(model.matrix(~ 0 + seasons))
    trait <- as.data.frame(model.matrix(~ 0 + species))
    colnames(trait) <- species_names
    
    # Model with chosen explanatory variables of interest
    model_dbrda <- dbrda(formula=class_list ~ ., data=constrain, distance="bray", metaMDSdist=FALSE, add=TRUE, sqrt.dist=FALSE)
    model_dbrda_constraints <- summary(model_dbrda)$constr.chi / summary(model_dbrda)$tot.chi
    model_dbrda_scores <- scores(model_dbrda, choices=c(1:20))
    
    model_dbrda_ef <- envfit(formula=model_dbrda ~ ., data=trait, choices=c(1:20), perm=1000)
    model_fit <- data.frame(r2=c(model_dbrda_ef$vectors$r,model_dbrda_ef$factors$r),
                            pvals=c(model_dbrda_ef$vectors$pvals,model_dbrda_ef$factors$pvals) )
    rownames(model_fit) <- c(names(model_dbrda_ef$vectors$r),names(model_dbrda_ef$factors$r))
    #write.csv(model_fit, file="seasons_classes_varpart.csv", row.names=TRUE)
    write.csv(model_fit, file="table_s1.csv", row.names=TRUE)
}



# ---------- Test the influence of traits on the compound classes ----------
f.ms2_classify_eco <- function() {
    # Count numbers of matched classes in each sample
    class_list <- data.frame(species=species, seasons=seasons)
    for (i in classes) {
        vec <- NULL
        for (j in 1:length(mzml_names)) {
            vec <- as.numeric(c(vec, sum(grepl(pattern=i, x=classifiers[[j]][,"Annotation..putative."]), na.rm=TRUE)))
        }
        class_list <- cbind(class_list, vec)
    }
    colnames(class_list) <- c("species", "seasons", gsub(x=classes, pattern=".*\\; ", replacement=""))
    rownames(class_list) <- mzml_names
    class_list <- class_list[,c(3:ncol(class_list))]
    class_list <- log(class_list)
    for (i in 1:ncol(class_list)) {
        class_list[,i] <- as.numeric(class_list[,i])
        class_list[which(is.infinite(class_list[,i])), i] <- 0
    }
    
    # Choose a model by permutation tests in constrained ordination
    # Model with intercept only
    model_0 <- dbrda(formula=class_list ~ 1, comm=class_list, data=traits, distance="bray", metaMDSdist=FALSE, add=TRUE, sqrt.dist=FALSE)
    
    # Model with all explanatory variables
    model_1 <- dbrda(formula=class_list ~ ., comm=class_list, data=traits, distance="bray", metaMDSdist=FALSE, add=TRUE, sqrt.dist=FALSE)
    
    # Goodness of fit statistic: Squared correlation coefficient
    model_dbrda_ef <- envfit(formula=model_1 ~ ., data=traits, perm=1000)
    model_fit <- data.frame(r2=c(model_dbrda_ef$vectors$r,model_dbrda_ef$factors$r),
                            pvals=c(model_dbrda_ef$vectors$pvals,model_dbrda_ef$factors$pvals) )
    rownames(model_fit) <- c(names(model_dbrda_ef$vectors$r),names(model_dbrda_ef$factors$r))
    #write.csv(model_fit, file="classes_traits_stat.csv", row.names=TRUE)
    write.csv(model_fit, file="table_s2.csv", row.names=TRUE)
    
    # Find best model using stepwise model selection
    model_step <- ordistep(object=model_0, scope=formula(model_1), R2scope=TRUE, Pin=0.2, Pout=5, direction="both", perm.max=1000, trace=TRUE)
    model_step_scores <- scores(model_step)
    model_step_constraints <- summary(model_step)$constr.chi / summary(model_step)$tot.chi
    write.csv(as.data.frame(model_step$anova), file="classes_traits_model.csv", row.names=TRUE)
    
    # Remove species and seasons as only ecological characteristics are interesting to us
    #model_formula <- as.formula(model_step$terms)
    model_formula <- as.character(as.formula(model_step$terms))
    model_formula <- as.character(paste0(model_formula[2], model_formula[1], model_formula[3]))
    model_formula <- as.formula(gsub(x=model_formula, pattern="(species \\+ |season \\+ )", replacement=""))
    
    # Ordistepped dbRDA
    model_dbrda <- dbrda(formula=model_formula, data=traits, distance="bray", metaMDSdist=TRUE, sqrt.dist=FALSE, choices=c(1:20))
    ef_formula <- stats::update(model_formula, model_dbrda ~ .)
    ef_factors <- as.factor(sapply(strsplit(as.character(ef_formula)[[3]], "\\+"), function(x) { x <- gsub("(\\`|^ | $)","",x) }))
    model_dbrda_ef <- envfit(formula=ef_formula, data=traits, choices=c(1:20), perm=1000)
    model_dbrda_constraints <- summary(model_dbrda)$constr.chi / summary(model_dbrda)$tot.chi
    model_dbrda_scores <- scores(model_dbrda, choices=c(1:20))
    
    # Plot results
    #pdf(file=args[3], encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    plot(0, 0, xlim=c(min(model_dbrda_scores$sites[,1])-1, max(model_dbrda_scores$sites[,1])+1),
         ylim=c(min(model_dbrda_scores$sites[,2]), max(model_dbrda_scores$sites[,2])),
         xlab="dbRDA1", ylab="dbRDA2",
         main=paste("dbRDA: Traits"," (explained variance = ",round(model_cap_constraints,3),")",sep=''), type="n")
    points(model_dbrda_scores$sites[,c(1,2)], pch=16, col=species_samples_colors)
    plot(model_dbrda_ef, cex=0.6, p.max=1, col="black")
    legend("topleft", bty="n", pch=16, col=species_colors, pt.cex=0.8, cex=0.8, y.intersp=0.7, text.width=0.5, legend=sort(unique(species)))
    #dev.off()
    
    
    
    # Traits of interest
    trait <- traits[, c("onsite_moisture.wet", "onsite_moisture.dry", "onsite_moisture.damp", "onsite_exposition.N", "onsite_exposition.O")]
    
    # Model with chosen explanatory variables of interest
    model_dbrda <- dbrda(formula=feat_list ~ ., comm=class_list, data=traits, distance="bray", metaMDSdist=FALSE, add=TRUE, sqrt.dist=FALSE)
    model_dbrda_ef <- envfit(formula=model_dbrda ~ ., data=traits, choices=c(1:20), perm=1000)
    model_fit <- data.frame(r2=c(model_dbrda_ef$vectors$r,model_dbrda_ef$factors$r),
                            pvals=c(model_dbrda_ef$vectors$pvals,model_dbrda_ef$factors$pvals) )
    rownames(model_fit) <- c(names(model_dbrda_ef$vectors$r),names(model_dbrda_ef$factors$r))
    #write.csv(model_fit, file="classes_traits_stat.csv", row.names=TRUE)
    trait <- traits[, c(rownames(model_fit[which(model_fit[,"pvals"] <= 0.06),]))]
    
    model_dbrda <- dbrda(formula=feat_list ~ ., comm=class_list, data=trait, distance="bray", metaMDSdist=FALSE, add=TRUE, sqrt.dist=FALSE)
    model_dbrda_ef <- envfit(formula=model_dbrda ~ ., data=trait, choices=c(1:20), perm=1000)
    model_dbrda_constraints <- summary(model_dbrda)$constr.chi / summary(model_dbrda)$tot.chi
    model_dbrda_scores <- scores(model_dbrda, choices=c(1:20))
    
    # Plot results
    #pdf(file=args[3], encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    plot(0, 0, xlim=c(min(model_dbrda_scores$sites[,3])-1, max(model_dbrda_scores$sites[,3])+1),
         ylim=c(min(model_dbrda_scores$sites[,4]), max(model_dbrda_scores$sites[,4])),
         xlab="dbRDA1", ylab="dbRDA2",
         main=paste("dbRDA: Traits"," (explained variance = ",round(model_cap_constraints,3),")",sep=''), type="n")
    points(model_dbrda_scores$sites[,c(3,4)], pch=16, col=seasons_samples_colors)
    plot(model_dbrda_ef, choices=c(3,4), cex=0.6, p.max=1, col="black")
    legend("topleft", bty="n", pch=16, col=species_colors, pt.cex=0.8, cex=0.8, y.intersp=0.7, text.width=0.5, legend=sort(unique(species)))
    #dev.off()
    
    
    
    # Trait of interest
    trait <- as.factor(as.character(traits$onsite_moisture))

    # Define colors and symbols for trait
    #classes_names <- gsub(x=classes, pattern=".*\\; ", replacement="")
    trait_names <- as.character(unique(trait))
    trait_colors <- rainbow(n=length(unique(trait)))
    trait_symbols <- seq(1, length(unique(trait)))
    trait_samples_colors <- sapply(trait, function(x) { x <- trait_colors[which(x==unique(trait))] } )
    
    # Count numbers of matched classes in each sample
    class_list <- data.frame()
    for (i in classes) {
        vec <- NULL
        for (j in 1:length(mzml_names)) {
            vec <- c(vec, sum(grepl(pattern=i, x=classifiers[[j]][,"Annotation..putative."]), na.rm=TRUE))
        }
        class_list <- rbind(class_list, vec)
    }
    class_list <- t(class_list)
    colnames(class_list) <- c(gsub(x=classes, pattern=".*\\; ", replacement=""))
    rownames(class_list) <- mzml_names
    
    # Sum of matched classes in a chosen trait
    model_class <- data.frame(matrix(0, ncol=length(classes), nrow=length(unique(trait))))
    colnames(model_class) <- gsub(x=classes, pattern=".*\\; ", replacement="")
    rownames(model_class) <- unique(trait)
    for (i in 1:length(classes)) {
        model_class[,i] <- as.numeric(lapply(X=unique(trait), FUN=function(x) { sum(as.numeric(class_list[which(trait==x), i])) } ))
    }
    
    # Boxplots
    pdf(file="classes_traits_boxplots.pdf", encoding="ISOLatin1", pointsize=10, width=16, height=12, family="Helvetica")
    par(mfrow=c(4,4), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=1.0, cex=0.9)
    for (i in 1:ncol(class_list)) {
        if (sum(class_list[,i]) > 0) {
            boxplot(class_list[,i] ~ trait, col=trait_colors, names=NA, main=colnames(class_list)[i], xlab="trait (onsite moisture)", ylab=colnames(class_list)[i])
            text(1:length(trait_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/10, srt=-22.5, adj=0.5, labels=trait_names, xpd=TRUE, cex=0.9)
            div_tukey <- tukey.test(response=class_list[,i], term=trait)
            text(1:length(trait_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/20, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
        }
    }
    dev.off()
    
    #model_anova <- aov(formula(response ~ term))
    #model_mc <- multcomp::glht(model_anova, multcomp::mcp(term="Tukey"))
    #model_cld <- multcomp::cld(summary(model_mc), decreasing=TRUE)
    #model_tukey <- data.frame("tukey_groups"=model_cld$mcletters$Letters)
    
    # Linear model for each class in the seasons, seasons: random factor
    model_class_lmer <- data.frame(matrix(0, ncol=length(classes), nrow=length(trait_names)))
    colnames(model_class_lmer) <- gsub(x=classes, pattern=".*\\; ", replacement="")
    rownames(model_class_lmer) <- trait_names
    
    model_class_lmer_p <- data.frame(matrix(0, ncol=length(classes), nrow=length(trait_names)))
    colnames(model_class_lmer_p) <- gsub(x=classes, pattern=".*\\; ", replacement="")
    rownames(model_class_lmer_p) <- trait_names
    
    # LMER for each class
    for (i in 1:length(classes)) {
        if (sum(class_list[,i]) > 0) {
            model_lmer <- lm(class_list[,i] ~ trait)
            model_class_lmer[,i] <- coef(summary(model_lmer))[,4]
            model_class_lmer_p[,i] <- print_p.values(coef(summary(model_lmer))[,4])
        } else {
            model_class_lmer[,i] <- rep(1, length(seasons_names))
            model_class_lmer_p[,i] <- rep("0", length(seasons_names))
        }
    }
    
    # Save models of statistics
    write.table(x=cbind(trait=trait_names,model_class), file="classes_traits.tsv", sep="\t", quote=TRUE, row.names=FALSE, dec=".")
    write.table(x=cbind(trait=trait_names,model_class_lmer), file="classes_traits_lm.tsv", sep="\t", quote=TRUE, row.names=FALSE, dec=".")
    write.table(x=cbind(trait=trait_names,model_class_lmer_p), file="classes_traits_lm_p.tsv", sep="\t", quote=TRUE, row.names=FALSE, dec=".")

    # Plot parallel coordinates
    pdf(file="classes_traits_parcoord.pdf", encoding="ISOLatin1", pointsize=8, width=16, height=8, family="Helvetica")
    parcoord(x=class_list, col=trait_colors, lty=1, lwd=2, las=3, xaxt="n", var.label=TRUE)
    dev.off()
}



# ---------- Test the influence of phylogeny on the compound classes ----------
f.ms2_classify_phylo <- function() {
    # Read phylogenetic tree
    phylo_tree <- read.tree("../data/moss_phylo.tre")
    
    # Replace names with species codes
    phylo_index <- match(substr(phylo_tree$tip.label,1,3), as.character(lapply(X=species_names, FUN=function(x) { x <- substr(x,1,3) })))
    phylo_tree$tip.label <- as.character(species_names[phylo_index])
    
    # Cophenetic distance matrix, needed for mpd and mntd calculations below
    phylo_dist <- cophenetic.phylo(phylo_tree)
    
    # Distance matrix of phylogenetic tree using Bray-Curtis
    phylo_dist <- vegdist(phylo_dist, method="bray")
    
    # Hierarchical clustering
    phylo_hclust <- hclust(phylo_dist, method="complete")
    
    
    # Sum of matched classes in species
    model_class <- data.frame(matrix(0, ncol=length(classes), nrow=length(species_names)))
    colnames(model_class) <- gsub(x=classes, pattern=".*\\; ", replacement="")
    rownames(model_class) <- species_names
    for (i in 1:length(classes)) {
        model_class[,i] <- as.numeric(lapply(X=species_names, FUN=function(x) { sum(as.numeric(class_list[which(class_list[,1]==x), i+1])) } ))
    }
    #class_list_species <- NULL
    #for (i in species_names) class_list_species <- rbind(class_list_species, apply(X=class_list[species==i,], MARGIN=2, FUN=function(x) { median(x) } ))
    #rownames(class_list_species) <- species_names
    
    # Reorder rows according to phylogenetic tree order
    model_class <- model_class[phylo_index,]
    
    # Distance matrix of feat_list using Bray-Curtis
    class_dist <- vegdist(model_class, method="bray")
    
    # Hierarchical clustering
    class_hclust <- hclust(class_dist, method="ward.D2")
    
    # Optimal order
    class_oclust <- reorder(class_hclust, seq(1:ncol(model_class)))

    # Procrustes analysis
    model_procrust <- protest(X=phylo_dist, Y=class_dist, permutations=10000)
    plot(model_procrust)
    
    # Mantel test
    model_mantel <- mantel(xdis=phylo_dist, ydis=class_dist, method="pearson", permutations=10000)
    
    # Correlation tests
    model_cor <- cor(phylo_dist, class_dist, method="pearson")
    model_cop <- cor_cophenetic(hclust(phylo_dist), hclust(class_dist), method="pearson")
    
    # Robinson-Foulds metric
    model_rfm <- RF.dist(phylo_tree, as.phylo(class_oclust), normalize=TRUE, check.labels=TRUE, rooted=TRUE)
    
    
    # Plot phylogenetic tree
    pdf("classes_phylo_corr.pdf", encoding="ISOLatin1", pointsize=12, width=8, height=5, family="Helvetica")
    par(mfrow=c(1,2), mar=c(1,1,2,1), cex=1.0)
    plot(phylo_tree, type="phylogram", direction="rightwards", x.lim=c(0,11), label.offset=0.4, use.edge.length=TRUE, show.tip.label=TRUE, tip.color=species_colors[phylo_index], font=2, main="")
    mtext(text="(a)", adj=0, line=0.5, font=2, cex=1.2)
    plot(as.phylo(class_oclust), type="phylogram", direction="leftwards", x.lim=c(0,0.35), label.offset=0.01, use.edge.length=TRUE, show.tip.label=TRUE, tip.color=species_colors[phylo_index], font=2, main="")
    mtext(text="(b)", adj=0, line=0.5, font=2, cex=1.2)
    dev.off()
    
    
    # r = Mantel
    print(paste("Mantel statistic:", round(model_mantel$statistic,3)))
    
    # c = cor_cophenetic
    print(paste("Correlation:", round(model_cop,3)))
    
    # rf = Robinson-Foulds
    print(paste("Robinson-Foulds metric:", round(model_rfm, 3)))
}


