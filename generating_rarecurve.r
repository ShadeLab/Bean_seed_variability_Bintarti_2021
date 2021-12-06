################################################################################################################################
###########  Generating Rarecurve using ggrare Function from gauravsk/ranacapa Package   #######################################


#devtools::install_github("gauravsk/ranacapa")
#install.packages("ranacapa")
library(ranacapa)


#run the ggrare function

ggrare <- function(physeq_object, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {

  x <- as(otu_table(physeq_object), "matrix")
  if (taxa_are_rows(physeq_object)) { x <- t(x) }

  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)

  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)

  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }

  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }

  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }

  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))

  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")

  #if (!is.null(label)) {
    #p <- p + ggplot2::geom_text(data = labels,
                                #ggplot2::aes_string(x = "x",
                                                    #y = "y",
                                                    #label = label,
                                                    #color = color),
                       #size = 4, hjust = 0)
  #}

  p <- p + ggplot2::geom_line(ggplot2::aes_string(), size=1)
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}


#set seed
set.seed(42)

####################################  Bacteria/Archaea  ###################################################

#rarefy the data
# data = phyloseq object of decontaminated non normalized otu table
p.rare <- ggrare(phyl.obj1, step = 1, color = "Plant", label = "Sample", se = FALSE)

#set up your own color palette
Palette <- c("#440154FF","#1F968BFF","#FDE725FF")
names(Palette) <- levels(sample_data(phyl.obj1)$Plant)
Palette

#plot the rarecurve

p <- p.rare + 
 #facet_wrap(~Plant, labeller = label_both)+
 theme_bw()+
 scale_color_manual(values = Palette)+
 scale_size_manual(values = 60)+
 labs(title = "(a)")+
 theme( strip.text.x = element_text(size=14, face='bold'),
        axis.text.x=element_text(size = 13),
        axis.text.y = element_text(size = 13),
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size =20 ,face='bold'),
        axis.title = element_text(size=15,face="bold"),
        legend.position = "none",
        #legend.title = element_text(size=15),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
        xlab("Sequencing depth (Reads)") +  ylab("Number of Bacterial/archaeal OTUs")
 

plot(p)

####################################  Fungi  ###################################################

#set seed
set.seed(42)

#rarefy the data
# data = phyloseq object of decontaminated non normalized otu table
p.rare.its <- ggrare(fil.obj, step = 1, color = "Plant", label = "Sample", se = FALSE)

#set up your own color palette
Palette <- c("#440154FF","#1F968BFF","#FDE725FF")
names(Palette) <- levels(sample_data(fil.obj)$Plant)
Palette

#plot the rarecurve
#p <- ggrare(psdata, step = 1000, color = "SampleType", label = "Sample", se = FALSE)
p.its <- p.rare.its + 
 #facet_wrap(~Plant)+
 theme_bw()+
 scale_color_manual(values = Palette)+
 labs(title = "(d)")+
 theme( strip.text.x = element_text(size=14, face='bold'),
        axis.text.x=element_text(size = 13),
        axis.text.y = element_text(size = 13),
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size =20 ,face='bold'),
        axis.title = element_text(size=15,face="bold"),
        legend.position = "none",
        #legend.title = element_text(size=15),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
        xlab("Sequencing depth (Reads)") +  ylab("Number of Fungal OTUs")
 

plot(p.its)


df_1 = data.frame(lat = rnorm(20), 
                  lon = rnorm(20), 
                  cor = c(rep('positive', 7), rep('negative', 13)), 
                  sign = c(rep(99, 5), rep(95, 6), rep(90,9)))





















