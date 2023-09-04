hack_Ideogram <- function(obj, subchr = NULL, which = NULL, xlabel = FALSE, cytobands = TRUE,
                          color = "red", fill = "red", alpha = 0.7,
                          zoom.region = NULL,
                          zoom.offset = 0.2, size = 1,
                          aspect.ratio = 1/20, ..., genome){
  if(missing(obj)){
    data(ideoCyto, package = "biovizBase")
    if(genome %in% names(ideoCyto)){
      obj <- ideoCyto[[genome]]
    }else{
      obj <- biovizBase::getIdeogram(genome = genome, subchr = subchr, cytobands = cytobands)
    }
  }
  ## do we need subchr here
  obj.ori <- obj
  if(!length(subchr)){
    subchr <- sort(unique(as.character(seqnames(obj))))[1]
    message("use ", subchr, " automatically")
    obj <- obj[seqnames(obj) == subchr]
    obj <- keepSeqlevels(obj, subchr)
  } else {
    obj <- hack_selectChromosome(obj, subchr)
  }
  if(!biovizBase::isIdeogram(obj))
    cytobands <- FALSE

  p <- ggplot2::ggplot() +
    ggbio::layout_karyogram(obj, cytobands = cytobands, geom = NULL)
  p <- hack_adjustZoom(obj, p, zoom.region, zoom.offset, color, fill, size, alpha)
  p <- hack_applyTheme(p, xlabel, subchr, aspect.ratio)

  p
}

hack_selectChromosome <- function(obj, subchr) {
  if(length(subchr)) {
    obj <- obj[seqnames(obj) == subchr]
    obj <- keepSeqlevels(obj, subchr)
  } else {
    subchr <- sort(unique(as.character(seqnames(obj))))[1]
    message("use ", subchr, " automatically")
    obj <- obj[seqnames(obj) == subchr]
    obj <- keepSeqlevels(obj, subchr)
  }
  if(length(unique(as.character(seqnames(obj)))) > 1)
    stop("Mulptiple chromosome information found")
  obj
}

hack_adjustZoom <- function(obj, plot, zoom.region, zoom.offset, color, fill, size, alpha) {
  if(length(zoom.region)) {
    if(length(zoom.region) != 2)
      stop("zoom.region must be a numeric vector of length 2")
    zoom.df <- data.frame(x1 = zoom.region[1],
                          x2 = zoom.region[2],
                          y1 = 0 - zoom.offset,
                          y2 = 10 + zoom.offset,
                          seqnames = unique(as.character(seqnames(obj))))
    plot <- plot + ggplot2::geom_rect(data = zoom.df,
                                      do.call(aes, list(xmin = substitute(x1),
                                                        xmax = substitute(x2),
                                                        ymin = substitute(y1),
                                                        ymax = substitute(y2))),
                                      color = color, fill = fill, size = size,
                                      alpha = alpha)
  } else {
    plot
  }
}

hack_applyTheme <- function(plot, xlabel, subchr, aspect.ratio) {
  plot <- plot + ggbio::theme_alignment(grid = FALSE, ylabel = TRUE, border = FALSE) +
    ggplot2::scale_y_continuous(breaks = 5, labels = subchr) +
    ggplot2::theme(strip.background = ggplot2::element_rect(colour = 'NA', fill = 'NA')) +
    ggplot2::theme(strip.text.y = ggplot2::element_text(colour = 'white')) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab("")
  plot <- plot + ggplot2::theme(aspect.ratio = aspect.ratio, axis.ticks.y = ggplot2::element_blank())
  if(!xlabel)
    plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  plot
}
