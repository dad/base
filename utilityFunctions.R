
no.legend <- theme(legend.position="none")

geom_diagline <- function(linetype='solid',size=0.1,colour="grey20",...) {
    geom_abline(slope=1,intercept=0,linetype=linetype,colour=colour)
}

scientific_10 <- function(x) {
    xout <- gsub("1e", "10^{", format(x),fixed=TRUE)
    xout <- gsub("{-0", "{-", xout,fixed=TRUE)
    xout <- gsub("{+", "{", xout,fixed=TRUE)
    xout <- gsub("{0", "{", xout,fixed=TRUE)
    xout <- paste(xout,"}",sep="")
    return(parse(text=xout))
}

scale_x_log10nice <- function(name=waiver(),omag=seq(-10,20),...) {
    breaks10 <- 10^omag
    scale_x_log10(name,breaks=breaks10,labels=scientific_10(breaks10),...)
}

scale_y_log10nice <- function(name=waiver(),omag=seq(-10,20),...) {
    breaks10 <- 10^omag
    scale_y_log10(name,breaks=breaks10,labels=scientific_10(breaks10),...)
}

scale_loglog <- function(xname=waiver(), yname=waiver(), ...) {
    list(scale_x_log10nice(name=xname, ...),scale_y_log10nice(name=yname,...))
}

scale_x_log2nice <- function(name=waiver(),omag=seq(-6,6),...) {
    breaks2 <- 2^omag
    scale_x_log10(name,breaks=breaks2,
                  labels=parse(text=paste("2^{",omag,"}")),...)
}

scale_y_log2nice <- function(name=waiver(),omag=seq(-6,6),...) {
    breaks2 <- 2^omag
    scale_y_log10(name,breaks=breaks2,
                  labels=parse(text=paste("2^{",omag,"}")),...)
}

scale_loglog2 <- function(...) {
    list(scale_x_log2nice(...),scale_y_log2nice(...))
}

theme_density <- function(...) {
    list(scale_y_continuous(expand=c(0.01,0.01)),
         theme(panel.border=element_blank(),
               panel.background=element_blank(),
               axis.line=element_blank(),
               axis.text.y=element_blank(),
               axis.title.y=element_blank(),
               axis.ticks.y=element_blank()),...)
}

# No legend for plots
# Usage: ggplot(d, aes(...)) + no.legend
no.legend <- theme(legend.position="none")


logfinite <- function(x) {
  xl <- log(x)
  xl[!is.finite(xl)] <- NA
  xl
}

corfinite <- function(x,y=NULL,use='pairwise.complete.obs',method="pearson") {
    # wrapper to calculate correlation coefficients for only finite values
    # useful for correlations of log-transformed values
  if (!is.null(y)) {
    niceinds <- which(is.finite(x) & is.finite(y))
    res <- cor(x[niceinds],y[niceinds],method=method)
  } else {
    x[!is.finite(x)] <- NA
    res <- cor(x,use=use,method=method)
    
  }
  res
}

logcor <- function(x,y=NULL,use='pairwise.complete.obs',method="pearson") {
    # wrapper for correlations of log-transformed values
  x <- as.matrix(x)
  ly <- y
  if (!is.null(y)) {ly <- log(as.matrix(y))}
  corfinite(log(x),ly,use=use,method=method)
}

odds <- function(p) {
  p/(1-p)
}

odds2p <- function(o) {
  o/(1+o)
}

p2odds <- function(p) {
  res <- odds(p)
  res[!is.finite(res)] <- NA
  res
}

logodds <- function(x) {
  # log odds, a shortcut 
  res <- log(odds(x))
  res[!is.finite(res)] <- NA
  res
}

invlogodds <- function(x) {
  # inverse log odds: a = invlogodds(logodds(a))
  y <- exp(x)
  y/(1+y)
}


##########
# Compute intracellular concentration of proteins and protein complexes
# for haploid budding yeast grown in rich medium
# 
# Assumes input data x has columns "gene" and "prot": common gene name, and protein level in molecules/cell
# 
# Typical call: 
# x <- read_tsv("scer-mrna-protein-absolute-estimate.txt",comment="#") # From Csardi et al. 2015
# genes <- c('SSA1','SSA2','SSA3','SSA4','HSP26','GBP2','NUG1','LSG1','PAB1','HSP104','SSE1','SSE2','YDJ1','SIS1')
# stoichiometry <- c('HSP104'=6,'SIS1'=2, 'YDJ1'=2)
# localization <- c("GBP2"="nucleus", "NUG1"="nucleolus", "LSG1"="cytosol","PAB1"="cytosol")
# conc <- intracellular.concentration(bg, genes, stoichiometry, localization)
# 
intracellular.concentration <- function(x, genes, stoichiometry=NULL, localization=NULL) {
  # Haploid cell volume ~40 fL in G1, S phase
  # https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=15&id=100427
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3404734/
  cell.volume <- 40e-15 # L (42 fL) cell volume
  # Cytosol volume is ~50% of cell volume
  cytosol.volume <- 0.5 * cell.volume
  # Nuclear volume is ~10% of cell volume
  nuclear.volume <- 0.1 * cell.volume
  # Nucleolar volume is 20% of nucleus volume
  nucleolar.volume <- 0.2 * nuclear.volume
  subcellular.volumes <- list("cell"=cell.volume, "cytosol"=cytosol.volume, "nucleus"=nuclear.volume, "nucleolus"=nucleolar.volume)
  # prot = Protein level in molecules/cell
  N_avo <- 6.022e23 # Avogadro's number, molecules/mole
  # Stoichiometry = number of monomers in complex whose concentration is desired
  # E.g. Hsp104 = 6, provided as list("HSP104"=6, "SIS1"=2) etc.
  stoichs <- tibble(gene=x$gene, stoich=1) # Assume monomers
  if (!is.null(stoichiometry)) {
    # Update proteins with known stoichiometries
    stoichs[match(names(stoichiometry), stoichs$gene),'stoich'] <- as.numeric(stoichiometry)
  }

  vols <- tibble(gene=x$gene, subcellular.localization="cell", subcellular.volume=cell.volume) # Assume localization throughout cell
  if (!is.null(localization)) {
    # Update proteins with known subcellular localizations
    vols[match(names(localization), vols$gene),'subcellular.localization'] <- as.character(localization)
    vols[match(names(localization), vols$gene),'subcellular.volume'] <- as.numeric(subcellular.volumes[as.character(localization)])
  }
  
  # Concentration in micromolar given number of molecules and volume
  # concentration (uM) = molecules/cell * 1 mol/N_avo molecules * (1/volume(L)) * 1e6 uM/M
  conc.uM <- function(prot, vol) {
    prot*(1/N_avo)*(1/vol)*1e6
  }

  # Cellular concentrations:
  res <- x %>% filter(gene %in% genes) %>% left_join(stoichs, by='gene') %>% left_join(vols, by='gene') %>%
    mutate(
    # Monomer
    cell.monomer.conc.uM=conc.uM(prot,cell.volume),
    # Complex
    cell.complex.conc.uM=conc.uM(prot,cell.volume)/stoich,
    # Complex within subcellular location
    subcellular.complex.conc.uM=conc.uM(prot,subcellular.volume)/stoich,
    )
  res <- res %>% select(gene, prot, stoich, subcellular.localization, cell.monomer.conc.uM, cell.complex.conc.uM, subcellular.complex.conc.uM)
  if (is.null(localization)) {
    res <- res %>% select(-subcellular.localization, -subcellular.complex.conc.uM)
  }
  if (is.null(stoichiometry)) {
    res <- res %>% select(-stoich, -cell.complex.conc.uM)
  }
  # Arrange in alphabetical order by gene name
  res %>% arrange(-desc(gene))
}
