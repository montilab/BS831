#AlignPlots: solution to vertically align ggplot2 plots
#source: http://stackoverflow.com/questions/26159495/align-multiple-ggplot-graphs-with-and-without-legends

#' @import ggplot2
#' @import data.table
AlignPlots <- function(...) {
  LegendWidth <- function(x) x$grobs[[8]]$grobs[[1]]$widths[[4]]

  plots.grobs <- lapply(list(...), ggplotGrob)

  max.widths <- do.call(unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })

  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]],
                                      unit(abs(diff(c(LegendWidth(x),
                                                      max.legends.width))),
                                           "mm"))
    }
    x
  })

  plots.grobs.eq.widths.aligned
}

#' @import Biobase
#' @export
scale_row<-function(eset){
	rowz<-t(apply(Biobase::exprs(eset), 1, function(z)
			scale(z)))
	Biobase::exprs(eset)<-rowz
	return(eset)
}

#merge a list of color labels,breaks,values into single list
#' @export
merge_labels<-function(x){
	l1 <- do.call('c', x)
	l2<-lapply(split(l1,sub('.*\\.', '', names(l1))),
                      function(i) do.call(c, i))
	inds<-which(!duplicated(l2$col_breaks))
	l3<-lapply(l2, function(i) i[inds])
}

#helper function for extracting ggplot legend
#' @import ggplot2
#' @export
g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

#make a single color legend
#' @import ggplot2
#' @export
make_legend<-function(hmcolors,
	... #other parameters in theme()
	){
	names(col_values)<-col_breaks
	df<-data.frame(num=1:length(col_breaks), breaks=col_breaks, values=col_values, labels=col_labels)
	p<-ggplot(df, aes(x = breaks, y = breaks, fill = factor(breaks)))+geom_tile(size = 1)+
	scale_fill_manual(values = col_values, breaks = col_breaks, labels = col_labels,
		guide = guide_legend(title = ""))+
	guides(fill = guide_legend(title = "",
 			override.aes = list(colour = "black"))) +
	theme(legend.text.align = 0,
			legend.justification = c(0,0), ...)
	p.legend<-g_legend(p)
	return(p.legend)
}

#make a list of color legends
#' @import ggplot2
#' @export
make_legend_list<-function(x, #x in the format of named list(x1, x2...)
	#e.g. x1 = list(col_breaks, col_values, col_labels)
	... #other parameters in theme()
	){

	g_legend<-function(a.gplot){
	    tmp <- ggplot_gtable(ggplot_build(a.gplot))
	    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	    legend <- tmp$grobs[[leg]]
	    legend
	}

	p.legend<-list()

	for(i in 1:length(x)){
		legend.title<-names(x)[i]
		col_breaks<-x[[i]]$col_breaks
		col_values<-x[[i]]$col_values
		col_labels<-x[[i]]$col_labels

		names(col_values)<-col_breaks
		df<-data.frame(num=1:length(col_breaks), breaks=col_breaks,
			values=col_values, labels=col_labels)
		p<-ggplot(df, aes(x = breaks, y = breaks, fill = factor(breaks))) +
		geom_tile(size = 1) +
		scale_fill_manual(values = col_values, breaks = col_breaks, labels = col_labels,
			guide = guide_legend(title = legend.title)) +
		guides(fill = guide_legend(title = legend.title,
	 			override.aes = list(colour = "black"))) +
		theme(legend.text.align = 0,
			legend.justification = c(0,0), ...)
		p.legend[[i]]<-g_legend(p)
	}

	p.legend$ncol <-1
	p<-do.call(arrangeGrob, p.legend)
	return(p)
}

#make a single color legend
#' @import ggplot2
#' @export
make_legend_continuous<-function(hmcolors,
	clow = -3, chigh = 3, n = 25, label = "expression"
	){

	col_breaks<-clow:chigh

	hmc<-hmcolors(name = label)

	df<-lapply(1:n, function(x){
		res<-lapply(1:n, function(y){
			data.frame(y = y, breaks = runif(n, clow, chigh))
			})
		res<-do.call(rbind, res)
		data.frame(x = x, res)
		})

	df<-do.call(rbind, df)
	#force dummy to be span range of clow and chigh
	df$breaks[1]<-clow
	df$breaks[2]<-chigh
	df$breaks<-as.numeric(df$breaks)

	p<-ggplot(df, aes(x = x, y = y, fill = breaks)) +
		geom_tile(size=1) +
		hmc +
		theme(legend.text.align = 0, legend.justification = c(0,0))
	p.legend<-g_legend(p)
	return(p.legend)
}

#' @import Biobase
#' @export
clust_eset<-function(eset){
	mat<-Biobase::exprs(eset)

	#column clustering
	#using euclidean distance, ward.D agglomeration
	distC <- function(x) dist(t(x), method="euclidean")
	dist_c<-distC(mat)
	hc<-hcopt(dist_c, method = "ward.D")

	#row clustering
	#using 1-cor as distance, ward.D agglomeration
	distR <- function(x) stats::as.dist(1- cor(t(x)))
	dist_r<-distR(mat)
	hr<-hcopt(dist_r, method = "ward.D")

	return(list(hc = hc, hr = hr))
}

#' \code{ggheat.make.groups} splits eset into list of esets based on phenotype label, do clustering within groups
#' @param eset expression set
#' @param labelcol column name for grouping in pData(eset)
#' @param labelvals values to group on (e.g. factor levels of pData(eset)[, labelcol])
#' clustFUN clustering function for eset, e.g. clust_eset
#' @import ggplot2
#' @export
ggheat.make.groups<-function(eset,
	labelcol,  #column name for grouping in pData(eset)
	labelvals, #values to group on (e.g. factor levels of pData(eset)[, labelcol])
	clustFUN, #clustering output for eset
	fixRowOrd = TRUE
	){

	#recommended fix to keep horizontal alignment
	nmax<-max(sapply(colnames(eset), function(i) nchar(i)))
	colnames(eset)<-sapply(colnames(eset), function(i)
	paste( paste(rep("  ", nmax - nchar(i)+1), collapse = ""), i, sep = ""))

	labelvec<-as.character(pData(eset)[, labelcol])
	numvals<-sapply(labelvals, function(i) sum(which(labelvec %in% i)))
	if(any(numvals < 3))
		stop("one or more groups is missing or has less than 3 members")

	reslist<-lapply(labelvals, function(i){
		eseti<-eset[, labelvec %in% i]
		clusti<-clustFUN(eseti)
		return(list(eseti = eseti, clusti = clusti))
		})
	esetlist<-lapply(reslist, function(i)
		i$eseti)
	hclist<-lapply(reslist, function(i)
		i$clusti$hc)

	hrlist<-lapply(reslist, function(i)
		i$clusti$hr)

	if(fixRowOrd){
		##fixed row clustering for all samples
		hr.all<-clustFUN(eset)$hr
		hrlist<-lapply(1:length(labelvals), function(i){hr.all})
	}

	return(list(esetlist = esetlist, hclist = hclist, hrlist = hrlist))
}

#helper function for plotting discretized heatmap
#' \code{ggheat.continuous} helper function for plotting ggheatmap, see \code{ggheat.continuous.single} and
#' \code{ggheat.continuous.group} for examples of usage
#' @import ggplot2
#' @export
ggheat.continuous<-function(eset,
	hc = NA, #hcopt for column leave NA for no ordering
	hr = NA, #hcopt for row leave NA for no ordering
	hmcolors = NA,
	col_lab,
	col_values,
	col_breaks,
	col_labels,
	ylabstr = "",
	type = c("left", "right", "middle", "regular"),
	fout = NA,
	p.heights = c(1.5, 0.5, 5),
	xsize = 4,
	ysize = 4,
	ysizelab = 7,
	override.hc = NA
	){

	theme_none <- theme(
	  panel.grid=element_blank(),
	  panel.grid.major=element_blank(),
	  panel.grid.minor=element_blank(),
	  panel.background=element_blank(),
	  axis.title.x=element_blank(),
	  axis.title.y=element_blank(),
	  axis.text.x=element_blank(),
	  axis.text.y=element_blank(),
	  axis.line=element_blank(),
	  axis.ticks.x=element_blank(),
	  axis.ticks.y=element_blank(),
	  plot.margin = unit(c(0,0.1,0,0), "lines"),
	  legend.margin = margin(6,6,6,6),
	  legend.key = element_rect(colour = "black"),
	  strip.background=element_blank(),
	  panel.spacing=unit(0, "cm"),
	  panel.border=element_blank(),
	  plot.background=element_blank()
	)

	mat<-Biobase::exprs(eset)

	#column dendrogram
	#default no clustering
	col_ord<-1:ncol(mat)
	row_ord<-1:nrow(mat)

	if(length(hc) > 1){
		dd_col<-as.dendrogram(hc)
		col_ord<-order.dendrogram(dd_col)
		data_col <- dendro_data(dd_col, draw=FALSE)
		HC <- ggplot(segment(data_col)) +
		geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
		scale_x_continuous( expand=c(0,0),
			limits = c(min(data_col$labels$x)-0.5, max(data_col$labels$x)+0.5)) +
		scale_y_continuous(expand=c(0.0,0.0))+  theme_none+
		theme(plot.margin = unit(c(0.4,0.1,0,0), "lines")) #extra padding on the top margin
	} else
		HC<-NA

	#add override hc if specified
	if(suppressWarnings(!is.na(override.hc))){
		dd_col.override<-as.dendrogram(override.hc)
		data_col <- dendro_data(dd_col.override, draw=FALSE)
		HC <- ggplot(segment(data_col)) +
		geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
		scale_x_continuous( expand=c(0,0),
			limits = c(min(data_col$labels$x)-0.5, max(data_col$labels$x)+0.5)) +
		scale_y_continuous(expand=c(0.0,0.0))+  theme_none+
		theme(plot.margin = unit(c(0.4,0.1,0,0), "lines")) #extra padding on the top margin
	}

	#order rows
	#by number of non-zero elements
	if(length(hr) > 1 ){
		dd_row<-as.dendrogram(hr)
		row_ord<-order.dendrogram(dd_row)
	}

	matcopy<-mat

	#quick fix to avoid bug in ordering
	rn<-rownames(mat)
	cn<-colnames(mat)
	if(is.null(rn)) rn <-1:nrow(mat)
	if(is.null(cn)) cn <-1:ncol(mat)
	rownames(matcopy)<-paste("R", rn, sep = "")
	colnames(matcopy)<-paste("C", cn, sep = "")
	matcopy<-matcopy[row_ord, col_ord]
	mat<-matcopy
	rownames(mat)<-rn[row_ord]
	colnames(mat)<-cn[col_ord]

	dt <- data.table(melt(mat))
	#factor, otherwise order might screw up if they can be coerced into numbers
	#row levels
	Var1levels<-rownames(mat)
	#col levels
	Var2levels<-colnames(mat)

	dt$Var1<-factor(dt$Var1, levels= Var1levels)
	dt$Var2<-factor(dt$Var2, levels= Var2levels)

	#main heatmap
	text.y<-element_text(size = ysize)
	if(type %in% c("middle", "right")){
		text.y<-element_blank()
	}

	#default heatmap fill gradient
	if(suppressWarnings(is.na(hmcolors)[1])){
		warning("heatmap color gradient not specified, setting to default hmcolors")
		hmcolors<-function(... ) scale_fill_gradient2(low = "blue", mid = "white",
       high = "red", midpoint = 0, limits=c(-3,3), oob=squish, ...)
	}

	if(type %in% c("left", "middle", "regular"))
		scfill<-hmcolors(guide =FALSE)
	else
		scfill<-hmcolors(guide = guide_legend(title = ""))

	p<-ggplot(dt, aes(Var2,y=Var1, fill = value ))+
		geom_tile( size=1) +
		scfill +
		theme(axis.text.x = element_text(angle = 90, size = xsize, hjust = 1,
			margin=margin(0,0,0,0)),
			axis.text.y = text.y,
			plot.margin = unit(c(0,0,1,0), "lines"),
			axis.title.x = element_blank(),
			panel.grid.minor.x = element_blank(),
			panel.grid.minor.y = element_blank(),
			panel.grid.major.x = element_blank(),
			panel.grid.major.y = element_blank())+
		scale_x_discrete(expand=c(0,0)) +
	    scale_y_discrete(expand=c(0,0))

	if(type %in% c("left", "regular")){
		p<-p + ylab(ylabstr)
		if(ysize ==0)
		p<-p+theme(axis.ticks.y = element_blank())
	} else if (type %in% c("middle")){
		p<-p + theme(axis.title.y = element_blank(),
			axis.ticks.y = element_blank(),
			axis.text.y = element_blank(),
		)
	} else { #right, account for legend
 		p<-p + guides(fill = guide_legend(title = "",
 			override.aes = list(colour = "black")))+
 		theme(axis.title.y = element_blank(),
			axis.ticks.y = element_blank(),
			axis.text.y = element_blank(),
			legend.key = element_rect(colour="black", size=0.5),
			legend.position = "right")
	}

	#column labels
	columnlab<-pData(eset)[col_ord, col_lab]
	# refactor column labels
	for(i in colnames(columnlab)){
		j<-as.character(columnlab[,i])
		columnlab[,i]<-factor(j, levels = unique(j))
	}

	if(length(col_lab)>1)
		rownames(columnlab)<-1:nrow(columnlab)

	dtcol<-data.table(melt(as.matrix(columnlab)))
	names(col_values)<-col_breaks

	text.lab.y<-element_text(size = ysizelab)
	if(type %in% c("middle", "right")){
		text.lab.y<-element_blank()
	}

	if(type %in% c("left", "middle", "regular"))
		scfilllab<-scale_fill_manual(values = col_values, breaks = col_breaks, labels = col_labels,
			guide = FALSE)
	else
		scfilllab<-scale_fill_manual(values = col_values, breaks = col_breaks, labels = col_labels,
			guide = guide_legend(title = ""))

	if(length(col_lab) == 1){
		#edge case
		dtcol$Var2 <-rep(col_lab, nrow(dtcol))
		lims<-c(min(dtcol$Var1)-0.5, max(dtcol$Var1)+0.5)
		pcol<-ggplot(dtcol, aes(Var1,y=Var2, fill = value )) +
		geom_tile( size=1) +
		scfilllab  +
		theme_none +
		scale_x_discrete(expand=c(0,0)) +
    	scale_y_discrete(expand=c(0,0))

	} else {
		pcol<-ggplot(dtcol, aes(Var1,y=Var2, fill = value )) +
			geom_tile( size=1) +
			scfilllab  +
			theme_none	+
			scale_x_discrete(expand=c(0,0)) +
    		scale_y_discrete(expand=c(0,0))
	}

	if(type %in% c("left", "regular")){
		pcol<-pcol+theme(axis.text.y = element_text(size = ysizelab, hjust = 0))
	} else if(type %in% c("right")){
		pcol<-pcol+theme(legend.key = element_rect(colour="black", size=0.5),
			legend.position = "right") +
			guides(fill = guide_legend(title = "",
			override.aes = list(colour = "black")))
	}
	if( any(is.na(HC)) ) {
		p.heights <-p.heights[-1]
		plist<-suppressWarnings(AlignPlots(pcol, p))
	}
	else plist<-suppressWarnings(AlignPlots(HC, pcol, p))

	plist$ncol <-1
	plist$heights <- p.heights
	p.combined<-do.call(arrangeGrob, plist)
	if(!is.na(fout))
		ggsave(p.combined, file = fout)

	return(p.combined)
}

#' \code{ggheat.continuous.single} ggheatmap function for plotting a single eset
#' @param eset expression set
#' @param	hc column clustering from hclust or hcopt
#' @param	hr row clustering from hclust or hcopt
#' @param	hmcolors function for heatmap color gradient, default:
#'  hmcolors<-function(... ) scale_fill_gradient2(low = "blue",
#'	mid = "white", high = "red", midpoint = 0, limits=c(-3,3), oob=squish, ...)
#' @param	hmtitle title of heatmap color label, default: "expression"
#' @param	col_lab column color label to display, one or more of columns of pData(eset)
#' @param	col_legend named list of column color legends, names correspond to col_lab, see format:
#' col_legend<-list(COL1 = list(col_breaks = COL1levels,
#'								col_values = brewer.pal(length(COL1levels),"Set1"),
#'								col_labels = COL1levels),
#'		COL2 = list(col_breaks = COL2levels,
#'								col_values = brewer.pal(length(COL2levels),"Set2"),
#'								col_labels = COL2levels))
#' @param	ylabstr y axis title string,
#' @param	fout filename to save plot, default NA, save as object, to display to device, use grid.arrange
#' @param	p.heights panel heights default c(1.5, 0.5, 5),
#' @param	xsize default 4 x axis labels (uses colnames), set to 0 to turn off display
#' @param	ysize default 4, y axis labels (uses rownames), set to 0 to turn off display
#' @param	ysizelab default 7, text size for y axis title (specified in ylabstr)
#' @param	xright default 0.24, displayed width for right margin as fraction of whole plot width
#' @examples
#' #load expression data
#' data(tcga.subset.400g.200s)
#' dat
#'
#' #hclust for rows and columns
#' hc<-clust_eset(dat)
#'
#' pData(dat)$hclust.groups<-as.factor(cutree(hc$hc, k = 4))
#'
#' #scale expression by row
#' dat.scaled<-scale_row(dat)
#'
#' subtypelevels<-levels(dat$subtype)
#'
#' #color legends for column labels
#' col_legend<-list(subtype = list(col_breaks = subtypelevels,
#' 		col_values = brewer.pal(length(subtypelevels),"Set1"),
#' 		col_labels = subtypelevels),
#' 		hclust.groups = list(col_breaks = levels(pData(dat)$hclust.groups),
#' 		col_values = sapply(c("pink", "orange", "yellow", "cyan"), to.hex),
#' 		col_labels = levels(pData(dat)$hclust.groups)))
#'
#' #heatmap fill gradient
#' hmcolors<-function(... ) scale_fill_gradient2(low = "blue", mid = "white",
#'        high = "red", midpoint = 0, limits=c(-3,3), oob=squish, ...)
#'
#' p<-ggheat.continuous.single(eset = dat.scaled,
#' 	hc = hc$hc,
#' 	hr = hc$hr,
#' 	hmcolors = hmcolors,
#' 	hmtitle = "row-zscore GE",
#' 	col_lab = c("subtype", "hclust.groups"),
#' 	col_legend = col_legend,
#' 	ylabstr = "",
#' 	fout = NA,
#' 	p.heights = c(1.5, 0.5, 5),
#' 	xsize = 0,
#' 	ysize = 0,
#' 	ysizelab = 7,
#' 	xright = 0.18
#' 
#' @import ggplot2
#' @export
ggheat.continuous.single<-function(eset,
	hc,
	hr,
	hmcolors = NA,
	hmtitle = "expression",
	col_lab,
	col_legend,
	ylabstr = "",
	fout = NA,
	p.heights = c(1.5, 0.5, 5),
	xsize = 4,
	ysize = 4,
	ysizelab = 7,
	xright = 0.24,
	override.hc = NA){

	#default heatmap fill gradient
	if(suppressWarnings(is.na(hmcolors)[1])){
		warning("heatmap color gradient not specified, setting to default hmcolors")
		hmcolors<-function(... ) scale_fill_gradient2(low = "blue", mid = "white",
       high = "red", midpoint = 0, limits=c(-3,3), oob=squish, ...)
	}

	col_legend_vec<-merge_labels(col_legend)
	col_values<-col_legend_vec$col_values
	col_breaks<-col_legend_vec$col_breaks
	col_labels<-col_legend_vec$col_labels

	p1<-ggheat.continuous(eset,
		hc, #hcopt for column leave NA for no ordering
		hr, #hcopt for row leave NA for no ordering
		hmcolors,
		col_lab,
		col_values,
		col_breaks,
		col_labels,
		ylabstr,
		type="regular",
		fout =NA,
		p.heights,
		xsize,
		ysize,
		ysizelab,
		override.hc
		)

	pcol.legend<-make_legend_list(col_legend,
		legend.key.size =  unit(0.2, "in"),
		legend.text = element_text(size=10),
		legend.title = element_text(colour = 'black', face = "bold", size = 10))

	clow<-min(Biobase::exprs(eset))
	chigh<-max(Biobase::exprs(eset))

	hm.legend<-make_legend_continuous(hmcolors, clow = clow, chigh = chigh, n = 25, label = hmtitle)

	plist<-list()
	plist[[1]]<-p1
	plist[[2]]<-arrangeGrob(pcol.legend, hm.legend, nrow = 2, ncol = 1)
	plist$nrow <-1
	plist$widths<-c(1-xright, xright)

	plistdev<-do.call(grid.arrange, plist)
	if(!is.na(fout))
		ggsave(plistdev, file = fout)
	return(plist)
}

#' \code{ggheat.continuous.group} ggheatmap function for plotting groups of esets, within-group column clustering
#' @param esetlist list of expression sets, see ?heatmap.make.groups for preprocessing
#' @param	hclist list of column clustering from hclust or hcopt
#' @param	hrlist list of row clustering from hclust or hcopt
#' @param	hmcolors function for heatmap color gradient, default:
#'  hmcolors<-function(... ) scale_fill_gradient2(low = "blue",
#'	mid = "white", high = "red", midpoint = 0, limits=c(-3,3), oob=squish, ...)
#' @param	hmtitle title of heatmap color label, default: "expression"
#' @param	col_lab column color label to display, one or more of columns of pData(eset)
#' @param	col_legend named list of column color legends, names correspond to col_lab, see format:
#' col_legend<-list(COL1 = list(col_breaks = COL1levels,
#'								col_values = brewer.pal(length(COL1levels),"Set1"),
#'								col_labels = COL1levels),
#'		COL2 = list(col_breaks = COL2levels,
#'								col_values = brewer.pal(length(COL2levels),"Set2"),
#'								col_labels = COL2levels))
#' @param	ylabstr y axis title string,
#' @param	fout filename to save plot, default NA, save as object, to display to device, use grid.arrange
#' @param	p.heights panel heights default c(1.5, 0.5, 5),
#' @param	xsize default 4 x axis labels (uses colnames), set to 0 to turn off display
#' @param	ysize default 4, y axis labels (uses rownames), set to 0 to turn off display
#' @param	ysizelab default 7, text size for y axis title (specified in ylabstr)
#' @param	xleft default 0.15, displayed width for left panel + margin as fraction of whole plot width
#' @param	xright default 0.24, displayed width for right margin as fraction of whole plot width
#' @examples
#' data(tcga.subset.400g.200s)
#' dat
#'
#' subtypelevels<-levels(dat$subtype)
#'
#' grps<-ggheat.make.groups(eset = dat,
#' 	labelcol = "subtype",  #column name for grouping in pData(eset)
#' 	labelvals = subtypelevels, #values to group on (e.g. factor levels of pData(eset)[, labelcol])
#' 	clustFUN = clust_eset #clustering function for eset
#' 	)
#'
#' col_legend<-list(subtype = list(col_breaks = subtypelevels,
#' 		col_values = brewer.pal(length(subtypelevels),"Set1"),
#' 		col_labels = subtypelevels))
#'
#' esetlist<-lapply(grps$esetlist, function(i) scale_row(i))
#'
#' hmcolors<-function(... ) scale_fill_gradient2(low = "blue", mid = "white",
#'        high = "red", midpoint = 0, limits=c(-3,3), oob=squish, ...)
#'
#' p<-ggheat.continuous.group(esetlist,
#' 	grps$hclist,
#' 	grps$hrlist,
#' 	hmcolors,
#' 	hmtitle = "row-zscore GE",
#' 	col_lab = "subtype",
#' 	col_legend = col_legend,
#' 	ylabstr = "",
#' 	fout  = NA,
#' 	p.heights = c(1.5, 0.5, 5),
#' 	xsize = 0,
#' 	ysize = 0,
#' 	ysizelab = 7,
#' 	xleft = 0.10,
#' 	xright = 0.24)
#' 
#' @import ggplot2
#' @export
ggheat.continuous.group<-function(esetlist,
	hclist,
	hrlist,
	hmcolors,
	hmtitle = "expression",
	col_lab,
	col_legend,
	ylabstr = "",
	fout,
	p.heights = c(1.5, 0.5, 5),
	xsize = 4,
	ysize = 4,
	ysizelab = 7,
	xleft = 0.15,
	xright = 0.24){

	n<-length(esetlist)

	if(length(hclist) != length(hrlist) | n != length(hclist))
		stop("esetlist, hclist, hrlist lengths must be equal")
	if(n<2)
		stop("esetlist length must be greater than 1")

	#default heatmap fill gradient
	if(suppressWarnings(is.na(hmcolors)[1])){
		warning("heatmap color gradient not specified, setting to default hmcolors")
		hmcolors<-function(... ) scale_fill_gradient2(low = "blue", mid = "white",
       high = "red", midpoint = 0, limits=c(-3,3), oob=squish, ...)
	}

	col_legend_vec<-merge_labels(col_legend)
	col_values<-col_legend_vec$col_values
	col_breaks<-col_legend_vec$col_breaks
	col_labels<-col_legend_vec$col_labels

	p1<-ggheat.continuous(esetlist[[1]],
		hclist[[1]], #hcopt for column leave NA for no ordering
		hrlist[[1]], #hcopt for row leave NA for no ordering
		hmcolors,
		col_lab,
		col_values,
		col_breaks,
		col_labels,
		ylabstr,
		type="left",
		fout =NA,
		p.heights,
		xsize,
		ysize,
		ysizelab
		)

	nmid<-0

	pmid<-lapply(2:n, function(i){
		ggheat.continuous(esetlist[[i]],
		hclist[[i]], #hcopt for column leave NA for no ordering
		hrlist[[i]], #hcopt for row leave NA for no ordering
		hmcolors,
		col_lab,
		col_values,
		col_breaks,
		col_labels,
		ylabstr,
		type ="middle",
		fout =NA,
		p.heights,
		xsize,
		ysize,
		ysizelab)
		})

	nmid<-sapply(2:n, function(i){
		ncol(esetlist[[i]])
		})

	n1<-ncol(esetlist[[1]])
	ntot<-n1+sum(nmid)

	xrem<-1-xright
	f1<-(n1/ntot)*xrem
	fmid<-sapply(nmid, function(i) (i/ntot)*xrem)

	pcol.legend<-make_legend_list(col_legend,
		legend.key.size =  unit(0.2, "in"),
		legend.text = element_text(size=10),
		legend.title = element_text(colour = 'black', face = "bold", size = 10))
	clow<-min(unlist(lapply(esetlist, function(i){
		min(Biobase::exprs(i))
		})))
	chigh<-max(unlist(lapply(esetlist, function(i){
		max(Biobase::exprs(i))
		})))
	hm.legend<-make_legend_continuous(hmcolors,
	clow = clow, chigh = chigh, label = hmtitle)

	if(n == 2) plist<-list(p1, pmid[[1]])
	else {
		plist<-list()
		plist[[1]]<-p1
		for(i in 1:length(pmid)){
			plist[[i+1]]<-pmid[[i]]
		}
	}
	plist[[n+1]]<-arrangeGrob(pcol.legend, hm.legend, nrow = 2, ncol = 1)
	plist$nrow <-1
	plist$widths<-c(xleft+f1, fmid, xright)

	plistdev<-do.call(grid.arrange, plist)
	if(!is.na(fout))
		ggsave(plistdev, file = fout)
	return(plist)
}
