# ---------- Preparations ----------
model_div <- NULL
traits <- NULL



# ---------- Shannon Diversity ----------
shannon.diversity <- function(p) {
    # Based on Li et al. (2016)
    # Function is obsolete, as it returns same values than vegan::diversity(x, index="shannon")
    # Since we do not expect features with negative intensity,
    # we exclude negative values and take the natural logarithm
    if (min(p) < 0 || sum(p) <= 0) 
        return(NA)
    pij <- p[p>0] / sum(p) 
    -sum(pij * log(pij)) 
}



# ---------- Menhinick Diversity ----------
menhinick.diversity <- function(p) {
    # Based on: http://www.coastalwiki.org/wiki/Measurements_of_biodiversity#Species_richness_indices
    D_Mn <- length(p) / sqrt(vegan::specnumber(p))
}



# ---------- Tukey-Test ----------
tukey.test <- function(response, term) {
    model_anova <- aov(formula(response ~ term))
    model_mc <- multcomp::glht(model_anova, multcomp::mcp(term="Tukey"))
    model_cld <- multcomp::cld(summary(model_mc), decreasing=TRUE)
    model_tukey <- data.frame("tukey_groups"=model_cld$mcletters$Letters)
    return(model_tukey)
}



# ---------- Create diversity data frame ----------
f.ms1_div_model <- function() {
    # Create data frame
    model_div               <- data.frame(features=apply(X=bina_list, MARGIN=1, FUN=function(x) { sum(x) } ))
    model_div$unique        <- apply(X=uniq_list, MARGIN=1, FUN=function(x) { sum(x) } )
    model_div$richness      <- apply(X=feat_list, MARGIN=1, FUN=function(x) { menhinick.diversity(x) })
    model_div$diversity     <- apply(X=uniq_list, MARGIN=1, FUN=function(x) { shannon.diversity(x) })
    model_div$shannon       <- apply(X=bina_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
    model_div$pielou        <- apply(X=scale(feat_list, center=F), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
    model_div$chao          <- vegan::specpool2vect(X=vegan::specpool(feat_list, species), index="chao")
    model_div$simpson       <- apply(X=bina_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
    model_div$inverse       <- apply(X=bina_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
    model_div$fisher        <- apply(X=bina_list, MARGIN=1, FUN=function(x) { fisher.alpha(x) })
    model_div$concentration <- apply(X=as.data.frame(mzml_files), MARGIN=1, FUN=function(x) { xr <- xcmsRaw(x); x <- sum(xr@env$intensity) } )
    
    # Remove NAs if present
    model_div[is.na(model_div)] <- 0
    
    # Write csv with results
    write.csv(model_div, file=paste("diversity.csv",sep=""), row.names=TRUE)
    
    # Return results
    model_div <<- model_div
}



# ---------- Convert trait of a group to presence/abscence matrix ----------
make.presence.abscence.matrix <- function(trait, name) {
    # Construct presence/abscence matrix
    tmp <- model.matrix(~ 0 + sapply(samp, function(x) { x <- traits_single[which(x==code), trait] } ))
    
    # Handle colnames correctly
    colnames(tmp) <- unique(traits_single[,trait])
    
    # Remove blank column
    if (length(which(colnames(tmp) == "blank value")) > 0)
        tmp <- tmp[, - which(colnames(tmp) == "blank value")]
    
    # Attach trait name prior to levels
    colnames(tmp) <- paste(name, colnames(tmp), sep=".")
    
    return(tmp)
}



# ---------- Import traits ----------
f.import_traits <- function(filename) {
    # Load traits from Peters et al. (2018) 10.1002/ece3.4361 paper
    traits_single <- read.table(file=filename, sep=";", header=TRUE, stringsAsFactors=TRUE)
    
    # Apply traits to sample matrix
    traits <- data.frame(sample=mzml_names, species=species, season=seasons)
    traits$code <- sapply(species, function(x) { x <- traits_single$Code[which(x==traits_single$Code)] } )
    traits$type <- sapply(species, function(x) { x <- traits_single$Type[which(x==traits_single$Code)] } )
    traits$growth_form <- sapply(species, function(x) { x <- traits_single$Growth.form[which(x==traits_single$Code)] } )
    traits$habitat_type <- sapply(species, function(x) { x <- traits_single$Habitat.type[which(x==traits_single$Code)] } )
    traits$substrate <- sapply(species, function(x) { x <- traits_single$Substrate[which(x==traits_single$Code)] } )
    traits$life_strategy <- sapply(species, function(x) { x <- traits_single$Life.strategy[which(x==traits_single$Code)] } )
    traits$gametangia <- sapply(species, function(x) { x <- traits_single$Gametangia.distribution[which(x==traits_single$Code)] } )
    traits$mean_spore_size <- as.numeric(sapply(species, function(x) { x <- traits_single$Mean.spore.size[which(x==traits_single$Code)] } ))
    traits$sex_frequency <- sapply(species, function(x) { x <- traits_single$Sexual.reproduction.frequency[which(x==traits_single$Code)] } )
    traits$light_index <- as.numeric(sapply(species, function(x) { x <- traits_single$Light.index[which(x==traits_single$Code)] } ))
    traits$temperature_index <- as.numeric(sapply(species, function(x) { x <- traits_single$Temperature.index[which(x==traits_single$Code)] } ))
    traits$continentality_index <- as.numeric(sapply(species, function(x) { x <- traits_single$Continentality.index[which(x==traits_single$Code)] } ))
    traits$moisture_index <- as.numeric(sapply(species, function(x) { x <- traits_single$Moisture.index[which(x==traits_single$Code)] } ))
    traits$reaction_index <- as.numeric(sapply(species, function(x) { x <- traits_single$Reaction.index[which(x==traits_single$Code)] } ))
    traits$nitrogen_index <- as.numeric(sapply(species, function(x) { x <- traits_single$Nitrogen.index[which(x==traits_single$Code)] } ))
    traits$life_form_index <- sapply(species, function(x) { x <- traits_single$Life.form.index[which(x==traits_single$Code)] } )
    
    # Load traits from ISA-Tab
    traits_single <- read.table(file=filename, quote="\"", sep="\t", header=TRUE, stringsAsFactors=TRUE)
    
    # Apply traits to sample matrix
    traits <- data.frame(sample=mzml_names, species=species, season=seasons)
    samp <- mzml_names
    for (i in species_names) samp <- gsub(x=samp, pattern=paste0(i,"_.*"), replacement=paste0(i))
    code <- traits_single[c(1:length(samp)), "Source.Name"]
    
    # Apply traits to sample matrix
    #as.data.frame(model.matrix(~ 0 + seasons))
    #traits$code <- sapply(samp, function(x) { x <- traits_single[which(x==code), "Characteristics.Species.Code."] } )
    traits <- cbind(traits, make.presence.abscence.matrix(name="type", trait="Characteristics.Type."))
    traits <- cbind(traits, make.presence.abscence.matrix(name="growth_form", trait="Characteristics.Growth.form."))
    traits <- cbind(traits, make.presence.abscence.matrix(name="habitat_type", trait="Characteristics.Habitat.type."))
    traits <- cbind(traits, make.presence.abscence.matrix(name="substrate", trait="Characteristics.Substrate."))
    traits <- cbind(traits, make.presence.abscence.matrix(name="life_strategy", trait="Characteristics.Life.strategy."))
    traits <- cbind(traits, make.presence.abscence.matrix(name="gametangia", trait="Characteristics.Gametangia.distribution."))
    traits$mean_spore_size <- as.numeric(sapply(samp, function(x) { x <- traits_single[which(x==code), "Characteristics.Mean.spore.size."] } ))
    traits <- cbind(traits, make.presence.abscence.matrix(name="sex_frequency", trait="Characteristics.Sexual.reproduction.frequency."))
    traits$light_index <- as.numeric(sapply(samp, function(x) { x <- traits_single[which(x==code), "Characteristics.Ellenberg.Light.index."] } ))
    traits$temperature_index <- as.numeric(sapply(samp, function(x) { x <- traits_single[which(x==code), "Characteristics.Ellenberg.Temperature.index."] } ))
    traits$continentality_index <- as.numeric(sapply(samp, function(x) { x <- traits_single[which(x==code), "Characteristics.Ellenberg.Continentality.index."] } ))
    traits$moisture_index <- as.numeric(sapply(samp, function(x) { x <- traits_single[which(x==code), "Characteristics.Ellenberg.Moisture.index."] } ))
    traits$reaction_index <- as.numeric(sapply(samp, function(x) { x <- traits_single[which(x==code), "Characteristics.Ellenberg.Reaction.index."] } ))
    traits$nitrogen_index <- as.numeric(sapply(samp, function(x) { x <- traits_single[which(x==code), "Characteristics.Ellenberg.Nitrogen.index."] } ))
    traits <- cbind(traits, make.presence.abscence.matrix(name="life_form_index", trait="Characteristics.Ellenberg.Life.form.index."))
    traits <- cbind(traits, make.presence.abscence.matrix(name="onsite_substrate", trait="Characteristics.Onsite.Substrate."))
    traits <- cbind(traits, make.presence.abscence.matrix(name="onsite_light", trait="Characteristics.Onsite.Light."))
    traits <- cbind(traits, make.presence.abscence.matrix(name="onsite_moisture", trait="Characteristics.Onsite.Moisture."))
    traits <- cbind(traits, make.presence.abscence.matrix(name="onsite_exposition", trait="Characteristics.Onsite.Exposition."))

    # Return results
    traits <<- traits
}



# ---------- Variation partitioning ----------
f.species_varpart <- function() {
    # Variance partitioning (result depending on Curtis distance)
    vp_list <- log(feat_list)
    vp_list[is.na(vp_list)] <- 0
    vp_list[is.infinite(vp_list)] <- 0
    model_vp <- varpart(vp_list, ~ seasons, ~ species)
    
    # Plot results
    #pdf(file="ms1_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
    pdf(file="fig_s3a.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
    plot(model_vp, Xnames=c("seasons","species"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("darkgreen","darkblue"))
    dev.off()
    write.csv(vp_list, file="fig_s3a.csv", row.names=TRUE)
}



# ---------- Plot unique features ----------
f.species_div_unique <- function() {
    #pdf(paste("species_features_unique.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    pdf(paste("fig_s1a.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    boxplot(model_div$unique ~ species, col=species_colors, names=NA, main="Number of unique features", xlab="Species", ylab="number of unique features")
    text(1:length(species_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=species_names, xpd=TRUE, cex=0.9)
    div_tukey <- tukey.test(response=model_div$unique, term=species)
    text(1:length(species_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
    dev.off()
    write.csv(cbind(data.frame(species=species), model_div[,c("unique","shannon","pielou","concentration")]), file="fig_s1.csv", row.names=FALSE)
}



# ---------- Plot total features ----------
f.species_div_features <- function() {
    pdf(paste("species_features_total.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    boxplot(model_div$features ~ species, col=species_colors, names=NA, main="Number of features", xlab="Species", ylab="number of features")
    text(1:length(species_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=species_names, xpd=TRUE, cex=0.9)
    div_tukey <- tukey.test(response=model_div$features, term=species)
    text(1:length(species_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
    dev.off()
}



# ---------- Plot Shannon index ----------
f.species_div_shannon <- function() {
    #pdf(paste("species_features_div_shannon.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    pdf(paste("fig_s1b.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    boxplot(model_div$shannon ~ species, col=species_colors, names=NA, main="Shannon diversity (H\')", xlab="Species", ylab="Shannon diversity index (H\')")
    text(1:length(species_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=species_names, xpd=TRUE, cex=0.9)
    div_tukey <- tukey.test(response=model_div$shannon, term=species)
    text(1:length(species_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
    dev.off()
}



# ---------- Plot Pielou index ----------
f.species_div_pielou <- function() {
    #pdf(paste("species_features_div_pilou.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    pdf(paste("fig_s1c.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    boxplot(model_div$pielou ~ species, col=species_colors, names=NA, main="Pielou\'s evenness", xlab="Species", ylab="Pielou diversity index (J)")
    text(1:length(species_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=species_names, xpd=TRUE, cex=0.9)
    div_tukey <- tukey.test(response=model_div$pielou, term=species)
    text(1:length(species_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
    dev.off()
}



# ---------- Plot concentration ----------
f.species_div_concentration <- function() {
    #pdf(paste("species_features_div_concentration.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    pdf(paste("fig_s1d.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    boxplot(model_div$concentration ~ species, col=species_colors, names=NA, main="Concentration", xlab="Species", ylab="concentration")
    text(1:length(species_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=species_names, xpd=TRUE, cex=0.9)
    div_tukey <- tukey.test(response=model_div$concentration, term=species)
    text(1:length(species_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
    dev.off()
}



# ---------- Plot unique features ----------
f.seasons_div_unique <- function() {
    # R-bug: prevent sorting when using formula
    model_boxplot <- data.frame(summer=model_div$unique[seasons=="summer"],
                                autumn=model_div$unique[seasons=="autumn"],
                                winter=model_div$unique[seasons=="winter"],
                                spring=model_div$unique[seasons=="spring"])
    # plot
    #pdf(paste("seasons_features_unique.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    pdf(paste("fig_s2a.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    boxplot(x=model_boxplot, col=seasons_colors, main="Number of unique features", xlab="Seasons", ylab="number of unique features")
    #text(1:length(seasons_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=seasons_names, xpd=TRUE, cex=0.9)
    div_tukey <- tukey.test(response=as.numeric(apply(X=model_boxplot, MARGIN=1, FUN=function(x) { x })), term=seasons)
    text(1:length(seasons_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
    dev.off()
    write.csv(cbind(data.frame(seasons=seasons), model_div[,c("unique","shannon","pielou","concentration")]), file="fig_s2.csv", row.names=FALSE)
}



# ---------- Plot total features ----------
f.seasons_div_features <- function() {
    # R-bug: prevent sorting when using formula
    model_boxplot <- data.frame(summer=model_div$features[seasons=="summer"],
                                autumn=model_div$features[seasons=="autumn"],
                                winter=model_div$features[seasons=="winter"],
                                spring=model_div$features[seasons=="spring"])
    # plot
    pdf(paste("seasons_features_total.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    boxplot(x=model_boxplot, col=seasons_colors, main="Number of features", xlab="Seasons", ylab="number of features")
    #text(1:length(seasons_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=seasons_names, xpd=TRUE, cex=0.9)
    div_tukey <- tukey.test(response=as.numeric(apply(X=model_boxplot, MARGIN=1, FUN=function(x) { x })), term=seasons)
    text(1:length(seasons_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
    dev.off()
}



# ---------- Plot Shannon index ----------
f.seasons_div_shannon <- function() {
    # R-bug: prevent sorting when using formula
    model_boxplot <- data.frame(summer=model_div$shannon[seasons=="summer"],
                                autumn=model_div$shannon[seasons=="autumn"],
                                winter=model_div$shannon[seasons=="winter"],
                                spring=model_div$shannon[seasons=="spring"])
    # plot
    #pdf(paste("seasons_features_shannon.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    pdf(paste("fig_s2b.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    boxplot(x=model_boxplot, col=seasons_colors,  main="Shannon diversity (H\')", xlab="Seasons", ylab="Shannon diversity index (H\')")
    #text(1:length(seasons_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=seasons_names, xpd=TRUE, cex=0.9)
    div_tukey <- tukey.test(response=as.numeric(apply(X=model_boxplot, MARGIN=1, FUN=function(x) { x })), term=seasons)
    text(1:length(seasons_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
    dev.off()
}



# ---------- Plot Pielou index ----------
f.seasons_div_pielou <- function() {
    # R-bug: prevent sorting when using formula
    model_boxplot <- data.frame(summer=model_div$pielou[seasons=="summer"],
                                autumn=model_div$pielou[seasons=="autumn"],
                                winter=model_div$pielou[seasons=="winter"],
                                spring=model_div$pielou[seasons=="spring"])
    # plot
    #pdf(paste("seasons_features_pielou.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    pdf(paste("fig_s2c.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    boxplot(x=model_boxplot, col=seasons_colors, main="Pielou\'s evenness", xlab="Seasons", ylab="Pielou diversity index (J)")
    #text(1:length(seasons_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=seasons_names, xpd=TRUE, cex=0.9)
    div_tukey <- tukey.test(response=as.numeric(apply(X=model_boxplot, MARGIN=1, FUN=function(x) { x })), term=seasons)
    text(1:length(seasons_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
    dev.off()
}



# ---------- Plot concentration ----------
f.seasons_div_concentration <- function() {
    # R-bug: prevent sorting when using formula
    model_boxplot <- data.frame(summer=model_div$concentration[seasons=="summer"],
                                autumn=model_div$concentration[seasons=="autumn"],
                                winter=model_div$concentration[seasons=="winter"],
                                spring=model_div$concentration[seasons=="spring"])
    # plot
    #pdf(paste("seasons_features_concentration.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    pdf(paste("fig_s2d.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    boxplot(x=model_boxplot, col=seasons_colors, main="Concentration", xlab="Seasons", ylab="concentration diversity index (J)")
    #text(1:length(seasons_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=seasons_names, xpd=TRUE, cex=0.9)
    div_tukey <- tukey.test(response=as.numeric(apply(X=model_boxplot, MARGIN=1, FUN=function(x) { x })), term=seasons)
    text(1:length(seasons_names), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
    dev.off()
}


