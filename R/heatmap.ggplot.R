#####
# This function plots a heatmap from expression set object using ggplot2
# Customization for heatmap color range (options: row.scaling and z.norm)
# Customization for RowSideColor and ColSideColor, may include multiple features (options: col.lab, row.lab)
# required packaged not installed on scc4: ggdendro, url: http://cran.r-project.org/web/packages/ggdendro/index.html
# to install:
# 1) download ggdendro_0.1-15.tar.gz from url
# 2) R CMD INSTALL -l LOCAL_R_DIR ggdendro_0.1-15.tar.gz
# 3) make sure LOCAL_R_DIR is included in R package path: .libPaths( c( .libPaths(), "/path/to/LOCAL_R_DIR") )
#####
#' heatmap.ggplot2
#' 
#' \code{heatmap.ggplot2} is the main function to draw and/or save the heatmap from an 
#'expression set object
#'
#' @import Biobase ggplot2 reshape2 ggdendro grid gridExtra gtable RColorBrewer scales stats
#' @param eSet expression set object to plot
#' @param brewer.pal.name name of colorbrewer, see options at RColorBrewer::display.brewer.all()
#' @param brewer.pal.rev reverse colorbrewer (TRUE or FALSE)
#' @param brewer.numColors number of colors for heatmap palette default 11
#' @param col.legend.brewer string vector of hexdecimal values for column labels, see example p3 and p4 for how to determine annotation labels and assigning colors
#' @param row.legend.brewer string vector of hexdecimal values for row labels
#' @param col.clust perform column-wise hierarchical clustering (TRUE or FALSE)
#' @param row.clust perform row-wise hierarchical clustering (TRUE or FALSE)
#' @param col.clust.hc hc object to be passed to as.dendrogram, available when col.clust = TRUE, default NA: hc.col = hcopt(stats::dist(t(x)), method="ward.D")
#' @param row.clust.hc hc object to be passed to as.dendrogram, available when row.clust = TRUE, default NA: hc.row = hcopt(stats::as.dist(1-cor(t(x))),method="ward.D")
#' @param col.lab column labels to include: subset of pData colnames character vector
#' @param row.lab row labels to include: subset of fData colnames character vector
#' @param heatmap.y.text include y axis labels for heatmap, uses rownames(eSet) (TRUE or FALSE)
#' @param heatmap.x.text include x axis labels for heatmap, uses colnames(eSet)) (TRUE or FALSE)
#' @param heatmap.y.text.size text size for y axis labels
#' @param heatmap.x.text.size text size for x axis labels
#' @param heatmap.colorlegend.name name for heatmap color legend
#' @param title.text main title for the plot
#' @param col.legend.name character vector for subset of col.lab to include in color legend
#' @param row.legend.name character vector for subset of row.lab to include in color legend
#' @param legend.lab.max.char number of characters limit for legend labels, default = 15
#' @param row.scaling how should rows be scaled ("none", "quantile", "z-score.all", or "z-score.capped")
#' @param z.norm heatmap colors reflect z-scores rather than original values, can be TRUE if row.scaling is "none"
#' @param cuttree.col number of clusters for columns, default 0: do not show cluster assignment
#' @param cuttree.row number of clusters for rows, default 0: do not show cluster assignment
#' @param verbose return additional clustred annotation information as well as heatmap
#' @param grid.heights numeric vector of length 7 representing heights of subpanels, leave NA for default
#' @param grid.widths numeric vector of length 4 represeting weights of subpanels, leave NA for default
#' @param show prints the heatmap within execution of the function (TRUE or FALSE)
#' 
#' @examples
#' 
#' #Use example data #1, for data set information: ?eSet1
#' data(eSet.brca.100)
#' eSet1<-eSet.brca.100
#' eSet1<-eSet1[1:10,1:25]
#'
#' p1<-heatmap.ggplot2(eSet=eSet1, col.clust = TRUE, row.clust = TRUE, 
#' col.clust.hc = NA, row.clust.hc = NA,
#' col.lab = c("HER2_status", "ER_status", "PR_status", "TN_status"), row.lab = "", 
#' heatmap.y.text = TRUE, heatmap.x.text = TRUE,
#' heatmap.colorlegend.name = "RNASeq_expression",
#' title.text = "TCGA BRCA log2 RNA-seq expression, z-score row normalized",
#' col.legend.name = c("HER2_status", "ER_status", "PR_status", "TN_status"), 
#' row.legend.name = "", 
#' row.scaling = "z-score.capped", 
#' z.norm = FALSE, 
#' cuttree.col = 4, cuttree.row = 3,
#' verbose = FALSE, show = FALSE)
#' grid.arrange(p1)
#' 
#'
#' x<-exprs(eSet1)
#' hc.row<-hcopt(stats::as.dist(1-cor(t(x))),method="ward.D")
#' hc.col <- hcopt(stats::dist(t(x), method = "euclidean"), method="ward.D") 
#'
#' #Adding custom hclust object in col.clust.hc and row.clust.hc
#' p2 <- heatmap.ggplot2(eSet=eSet1, col.clust = TRUE, row.clust = TRUE, 
#' col.clust.hc = hc.col, row.clust.hc = hc.row,
#' col.lab = c("HER2_status", "ER_status", "PR_status", "TN_status"), row.lab = "", 
#' heatmap.y.text = TRUE, heatmap.x.text = TRUE,
#' heatmap.colorlegend.name = "RNASeq_expression",
#' title.text = "TCGA BRCA log2 RNA-seq expression, z-score row normalized",
#' col.legend.name = c("HER2_status", "ER_status", "PR_status", "TN_status"), 
#' row.legend.name = "", 
#' row.scaling = "z-score.capped", 
#' z.norm = FALSE, 
#' cuttree.col = 4, cuttree.row = 3,
#' verbose = FALSE, show = FALSE)
#' grid.newpage()
#' grid.arrange(p2)
#' 
#' 
#' #Saving plot in verbose format
#' p3 <- heatmap.ggplot2(eSet=eSet1, col.clust = TRUE, row.clust = TRUE, 
#' col.clust.hc = hc.col, row.clust.hc = hc.row,
#' col.lab = c("HER2_status", "ER_status", "PR_status", "TN_status"), row.lab = "", 
#' heatmap.y.text = TRUE, heatmap.x.text = TRUE,
#' heatmap.colorlegend.name = "RNASeq_expression",
#' title.text = "TCGA BRCA log2 RNA-seq expression, z-score row normalized",
#' col.legend.name = c("HER2_status", "ER_status", "PR_status", "TN_status"), 
#' row.legend.name = "", 
#' row.scaling = "z-score.capped", 
#' z.norm = FALSE, 
#' cuttree.col = 4, cuttree.row = 3,
#' verbose = TRUE, show = FALSE)
#' grid.newpage()
#' grid.arrange(p3$heatmap)
#' 
#'
#' #Adding custom colors to column and row annotation labels
#' print(p3$meta.c$id)
#' meta.c.color.string<-c("yellow", "khaki3", "gold", "chocolate", "darkred", "cyan")
#' meta.c.color<-as.character(sapply(meta.c.color.string, to.hex))
#' names(meta.c.color)<-c("Negative", "Positive", "1", "2", "3", "4")
#' print(p3$meta.r$id)
#' meta.r.color.string<-c("pink", "azure", "green")
#' meta.r.color<-as.character(sapply(meta.r.color.string, to.hex))
#' names(meta.r.color)<-c("1", "2", "3")
#'
#' p4<-heatmap.ggplot2(eSet=eSet1, 
#' col.legend.brewer = meta.c.color,
#' row.legend.brewer = meta.r.color,
#' col.clust = TRUE, row.clust = TRUE, 
#' col.clust.hc = hc.col, row.clust.hc = hc.row,
#' col.lab = c("HER2_status", "ER_status", "PR_status", "TN_status"), row.lab = "cluster.row", 
#' heatmap.y.text = TRUE, heatmap.x.text = TRUE,
#' heatmap.colorlegend.name = "RNASeq_expression",
#' title.text = "TCGA BRCA log2 RNA-seq expression, z-score row normalized",
#' col.legend.name = c("HER2_status", "ER_status", "PR_status", "TN_status", "cluster.col"), 
#' row.legend.name = "cluster.row", 
#' row.scaling = "z-score.capped", 
#' z.norm = FALSE, 
#' cuttree.col = 4, cuttree.row = 3,
#' verbose = FALSE, show = FALSE)
#' grid.newpage()
#' grid.arrange(p4)
#'
#' 
#' 
#' 
#' @import Biobase
#' @import ggplot2
#' @import reshape2
#' @import ggdendro
#' @import grid
#' @import gridExtra
#' @import gtable
#' @import RColorBrewer
#' @import scales
#' @import stats
#' @export 
heatmap.ggplot2<-function(eSet, 
	brewer.pal.name = "RdBu",#color gradient, see options at RColorBrewer::display.brewer.all()
	brewer.pal.rev = TRUE, #reverse color gradient
	brewer.numColors = 11,
	col.legend.brewer = "",
	row.legend.brewer = "",
	col.clust = TRUE, #column clustering
	row.clust = TRUE, #row clustering
	col.clust.hc = NA, #hc object to be passed to as.dendrogram, available when col.clust = TRUE, default NA: hc = hcopt(stats::dist(t(x)), method="ward.D")
	row.clust.hc = NA,#hc object to be passed to as.dendrogram, available when row.clust = TRUE, default NA: hc01.row <- hcopt(stats::as.dist(1-cor(t(x))),method="ward.D")
	col.lab = "", #column side labels, character vector of pData colnames subset to include in label
	row.lab = "", #row side labels, character vector of fData colnames subset to include in label
	heatmap.y.text = FALSE, #text labels for heatmap rows: rownames of heatmap matrix
	heatmap.x.text = FALSE, #text labels for heatmap columns: colnames of heatmap matrix
	heatmap.y.text.size = 7, #text labels size for heatmap rows
	heatmap.x.text.size = 6, #text labels size for heatmap columns
	heatmap.colorlegend.name = "val", #name for heatmap color legend
	title.text = "",
	col.legend.name = "", #character vector for colname to plot as col.lab legend
	row.legend.name = "", #character vector for colname to plot as row.lab legend
	legend.lab.max.char = 15, #number of characters limit for legend labels
	row.scaling = "none", #one of c("none", "quantile", "z-score.all", "z-score.capped") 
	z.norm = FALSE, #heatmap colors reflect z-scores rather than original values
					#only can be TRUE if row.scaling == "none" 
	cuttree.col = 0, #TO BE IMPLEMENTED: shows cluster assignment color label for columns
					#number of clusters for columns, default 0: do not show cluster assignment
	cuttree.row = 0, #TO BE IMPLEMENTED: shows cluster assignment color label for rows
					#number of clusters for rows, default 0: do not show cluster assignment
	verbose = FALSE,  #return addition information as well as heatmap
	grid.heights =NA, #numeric vector of length 7 representing heights of subpanels, leave NA for default
	grid.widths = NA, #numeric vector of length 4 represeting weights of subpanels, leave NA for default
	show = FALSE #prints the heatmap within execution of the function
	){

  #make blank ggplot2 theme
	theme_none <- theme(
	  panel.grid = element_blank(),
	  panel.grid.major = element_blank(),
	  panel.grid.minor = element_blank(),
	  panel.background = element_blank(),
	  axis.title.x = element_text(colour=NA),
	  axis.title.y = element_blank(),
	  axis.text.x = element_blank(),
	  axis.text.y = element_blank(),
	  axis.line = element_blank(), 
	  axis.ticks.x = element_blank(),
	  axis.ticks.y = element_blank(),
	  axis.ticks.margin = unit(0, "cm"),
	  plot.margin = unit(c(0.1,0.1,-1.2,-0.6), "lines"), # top,right, bottom, left
	  legend.margin =unit(0,"cm"),
	 # strip.background = element_blank(), 
	  panel.margin = unit(0, "cm"),
	  panel.border = element_blank(),
	  plot.background = element_blank(),
	  strip.background = element_blank()
	)
	
	meta.c<-NA
	meta.r<-NA
	col.meta<-NA
	row.meta<-NA

	#--blank panel for paddings in grid layout--
	p<-grid.rect(gp=gpar(col="white"), draw = F)

	#--main title--
	main.title<-textGrob(label=title.text,just=c("center","center"))

	#set default ordering (no clustering)
	x<-as.matrix(Biobase::exprs(eSet))

	row.ord<-1:dim(x)[1]
	col.ord<-1:dim(x)[2]

	clusMember.col<-NA
	clusMember.row<-NA

  	#--column dendrogram--
  	HC<-p
  	if (col.clust == TRUE){ 
  		#column distance: euclidean, clustering method: ward
  		if (is.na(col.clust.hc[1])) {
  			hc01.col <- hcopt(stats::dist(t(x)), method="ward.D") 
  			dd.col<-as.dendrogram(hc01.col)
  		}
  		else {
  			hc01.col <- col.clust.hc
  			dd.col<-as.dendrogram(hc01.col)
  		}
        
        col.ord<-order.dendrogram(dd.col)
	  	data_col <- dendro_data(dd.col, draw = F)
	  	HC <- ggplot(segment(data_col)) + 
	  	geom_segment(aes_string(x = "x", y = "y", xend = "xend", yend = "yend"))+
	  	scale_x_continuous( expand = c(0,0)) + 
	  	scale_y_continuous(expand = c(0.01,0.01))+
	  	theme_none 

	  	if (cuttree.col > 1){
	  		clusMember.col<-cutree(hc01.col, cuttree.col)
	  	}
  	}

  	#--row dendrogram--
  	HR<-p
  	if (row.clust ==TRUE){
  		if (is.na(row.clust.hc[1])) {
  			hc01.row <- hcopt(stats::as.dist(1-cor(t(x))),method="ward.D")
  			dd.row<-as.dendrogram(hc01.row)
  		}
  		else {
  			hc01.row <-row.clust.hc
  			dd.row<-as.dendrogram(hc01.row)
  		}
  		
  		row.ord<-order.dendrogram(dd.row)

  		##reverse row ordering

	  	data_row <- dendro_data(dd.row)
	  	HR <- ggplot(segment(data_row)) + 
	  	geom_segment(aes_string(x = "x", y = "y", xend = "xend", yend = "yend"))+
	    scale_x_continuous( expand = c(0,0)) + 
	    scale_y_continuous( expand = c(0.01, 0.01)) +
	    theme_none+
	    coord_flip()

	    if (cuttree.row > 1){
	  		clusMember.row<-cutree(hc01.row, cuttree.row)
	  	}
	}

	#--main heatmap--
  	
  	x.ordered <- x[row.ord, col.ord]

  	if (row.scaling == "quantile"){
		  x.ordered.scaled<-t(apply(x.ordered, 1, function(z) 
			cut(z, 
          breaks = quantile(z, (0:11)/11 ), 
          dig.lab =10, 
          include.lowest=T, 
          labels =F))) - 6
		  rownames(x.ordered.scaled)<-rownames(x.ordered)
		  colnames(x.ordered.scaled)<-colnames(x.ordered)
		  x.ordered<-x.ordered.scaled
	  } else if (row.scaling == "z-score.all" || row.scaling == "z-score.capped"){
		  x.ordered.scaled<-t(apply(x.ordered, 1, function(z) 
			scale(z))) 
		  rownames(x.ordered.scaled)<-rownames(x.ordered)
		  colnames(x.ordered.scaled)<-colnames(x.ordered)
		  x.ordered<-x.ordered.scaled
	  } 

	df<-as.data.frame(x.ordered)
	df$gene<-rownames(x.ordered)
	df$gene<-factor(df$gene, levels=unique(df$gene))
  
	dfm <- melt(df, id.vars = "gene")
	colnames(dfm) <- c("gene", "sample", heatmap.colorlegend.name)

	if (row.scaling =="none" & z.norm == TRUE){
    dfm[,heatmap.colorlegend.name]<-scale(dfm[,heatmap.colorlegend.name])
  }

 	if (brewer.pal.rev == TRUE)
	heatmap_color_scale<-scale_fill_gradientn(colours=rev(brewer.pal(brewer.numColors,brewer.pal.name)))
	else 
	heatmap_color_scale<-scale_fill_gradientn(colours=brewer.pal(brewer.numColors,brewer.pal.name))



	if (row.scaling == "z-score.capped"){
		if (brewer.pal.rev == TRUE)
		heatmap_color_scale<-scale_fill_gradientn(colours=rev(brewer.pal(brewer.numColors,brewer.pal.name)), limits=c(-3,3), oob=squish)
		else 
		heatmap_color_scale<-scale_fill_gradientn(colours=brewer.pal(brewer.numColors,brewer.pal.name), limits=c(-3,3), oob=squish)
	}

	M <- ggplot(dfm, aes_string(x="gene", y="sample")) + 
  	geom_tile(aes_string(fill=heatmap.colorlegend.name)) + 
  	heatmap_color_scale +
  	theme_none +
  	theme(legend.position="bottom") +
  	coord_flip() + 
  	theme(axis.text.x = element_text(angle = 90, hjust = 1, size = heatmap.x.text.size), 
	  			axis.text.y=element_text(hjust = 1, size =heatmap.y.text.size) ) 

  	#--color legend for heatmap--
  	ML <- gtable_filter(ggplot_gtable(ggplot_build(M)), "guide-box")

  	#--heatmap y axis label text
  	if(heatmap.y.text == TRUE){
  		MCT <- gtable_filter(ggplot_gtable(ggplot_build(M)), "axis.l") 
  	}
  	else{
  		MCT <- p
  	}

  	 if(heatmap.x.text == TRUE){
  		MRT <- gtable_filter(ggplot_gtable(ggplot_build(M)), "axis.b") 
  	}
  	else{
  		MRT <- p
  	}

  	M <- M + theme_none + theme(legend.position = "none" )

  	#--column side label--
  	if (cuttree.col>1){
  		pData(eSet)$cluster.col<- clusMember.col
  		if(col.lab[1] == ""){
  			col.lab <-"cluster.col"
  		} else {
  			col.lab <-c("cluster.col", col.lab)
  		}
  	}

  	if (col.lab[1] == ""){
  		#no specificied column side label
  		LC <- p
  		LCT <- p
  		SLC <- list(p)
  		names(SLC)<-"blank"

  	} else{

	  	col.meta<-pData(eSet)[col.ord, ]
	  	for (i in col.lab){
	  		col.meta[,i]<-strtrim(col.meta[,i], legend.lab.max.char)
	  	}
	 	meta.c<-data.frame(type = vector(), id = vector(), num = vector() )
		for (x in col.lab){
			meta.c<-rbind(meta.c, 
				data.frame(cbind(type = x, 
							id = as.character(col.meta[,x]), 
							num =  1:length(col.meta[,x]))) )
		}
		meta.c$num <- factor(meta.c$num, levels=unique(meta.c$num))
		meta.c$id <-factor(meta.c$id, levels = unique(meta.c$id))

		palette.old<- brewer.pal(11, "Spectral")
		getPalette <- colorRampPalette(palette.old)
		numColors <- length(unique(meta.c$id))
		#generate number of colored seed for all unique factor levels
		palette.all<-getPalette(numColors) 
		set.seed(57)
		palette.all.permute<-sample(palette.all, replace = FALSE, size = length(palette.all))

		if (col.legend.brewer[1] != ""){
			metacolunq<-as.character(unique(meta.c$id))

			if(is.null(names(col.legend.brewer))){
				stop("col.legend.brewer must have names representing ids of color labels")
			}
			if (any(duplicated(names(col.legend.brewer)))){
				stop("duplicated column legend colors in col.lgend.brewer")
			} 
			missingcols<-setdiff(names(col.legend.brewer), metacolunq)
			if (length(missingcols)!= 0){
				stop(paste("missing column legend colors", 
					paste(missingcols, collapse = ","), sep = ":"))
			}
			col.legend.brewer.ordered<- col.legend.brewer[ match( metacolunq,
				names(col.legend.brewer))]
			col.legend.brewer.ordered[which(is.na(col.legend.brewer.ordered))]<- to.hex("white")
			palette.all.permute<-col.legend.brewer.ordered
			#palette.all.permute <-col.legend.brewer
		}

		LC<-ggplot(meta.c, aes_string(x = "num", y = "type", fill = "id")) + 
			geom_tile() + 
			scale_y_discrete(expand =c(0,0)) + 
			scale_x_discrete(expand=c(0,0)) + 
	  		theme(legend.position="bottom", axis.ticks.y = element_blank(),
	  			legend.key.size = unit(0.15,"cm"),
	  			legend.text.align = 0, legend.title.align = 0, 
	  			axis.text.y=element_text(hjust = 1, size =7, vjust = 0)) 

	  	if (col.legend.name[1] != ""){
	  		SLC <-list()
	  		for (ind in rev(col.legend.name)){
	  			LC.temp<-LC + scale_fill_manual(name = ind, 
	  					breaks = levels(factor(col.meta[, ind])),
	                    values = palette.all.permute, 
	                    guide = guide_legend(direction = "vertical", 
                                           title.position = "top", 
                                           label.position="right")
	                    )
	            SLC[[ind]]<- gtable_filter(ggplot_gtable(ggplot_build(LC.temp)), "guide-box")

	  		}
	  		LC<-LC.temp
	  		
	  	} else {
	  		LC<-LC + scale_fill_manual(values = palette.all.permute)
	  		SLC <-list(p)
	  		names(SLC)<-"blank"
	  	}

	  	#--column side label names--
	    LCT <- gtable_filter(ggplot_gtable(ggplot_build(LC)), "axis.l") 

	    LC <-LC + theme_none + theme(legend.position = "none" )
	}

	if (cuttree.row>1){
  		fData(eSet)$cluster.row<- clusMember.row
  		if(row.lab[1] == ""){
  			row.lab <-"cluster.row"
  		} else{
  			row.lab <-c("cluster.row", row.lab)
  		}
  }

	if (row.lab[1] == ""){
  		#no specificied column side label
  		LR <- p
  		LRT <- p
  		SLR<-list(p)
	  	names(SLR)<-"blank"

  	} else {

	  	row.meta<-fData(eSet)[row.ord,]
	  	for (i in row.lab){
	  		row.meta[,i]<-strtrim(row.meta[,i], legend.lab.max.char)
	  	}
	 	meta.r<-data.frame(type = vector(), id = vector(), num = vector() )
		for (x in row.lab){
			meta.r<-rbind(meta.r, 
				data.frame(cbind(type = x, 
							id = row.meta[,x], 
							num =  1:length(row.meta[,x]))) )
		}

		meta.r$num <- factor(meta.r$num, levels=unique(meta.r$num))
		meta.r$id<- factor(meta.r$id, levels = unique(meta.r$id))
		palette.old<- brewer.pal(9, "Set1")
		getPalette <- colorRampPalette(palette.old)
		numColors <- length(unique(meta.r$id))
		#generate number of colored need for all unique factor levels
		palette.all<-getPalette(numColors) 

		set.seed(57)
		palette.all.permute<-sample(palette.all, replace = FALSE, size = length(palette.all))

		if (row.legend.brewer[1] != ""){
			metarowunq<-as.character(unique(meta.r$id))

			if(is.null(names(row.legend.brewer))){
				stop("row.legend.brewer must have names representing ids of color labels")
			}
			if (any(duplicated(names(row.legend.brewer)))){
				stop("duplicated column legend colors in col.lgend.brewer")
			} 
			missingrows<-setdiff(names(row.legend.brewer), metarowunq)
			if (length(missingrows)!= 0){
				stop(paste("missing column legend colors", 
					paste(missingrows, collapse = ","), sep = ":"))
			}
			row.legend.brewer.ordered<- row.legend.brewer[ match( metarowunq,
				names(row.legend.brewer))]
			row.legend.brewer.ordered[which(is.na(row.legend.brewer.ordered))]<- to.hex("white")
			palette.all.permute<-row.legend.brewer.ordered

			#palette.all.permute <-row.legend.brewer
		}

		LR<-ggplot(meta.r, aes_string(x = "num", y = "type", fill = "id")) + 
			geom_tile() + 
			scale_y_discrete(expand =c(0,0)) + 
			scale_x_discrete(expand=c(0,0)) + 
	  	theme(legend.position="bottom", 
        axis.ticks.x = element_blank(),
	  		legend.key.size = unit(0.15,"cm"),
	 			legend.text.align = 0, 
        legend.title.align = 0,
        axis.text.x = element_text(angle = -90, hjust = 0, size =7, vjust =0)) + 
        coord_flip()

	  	if (row.legend.name[1] != ""){
	  		SLR <-list()
	  		for (ind in rev(row.legend.name)){
	  			LR.temp<-LR + scale_fill_manual(name = ind, 
	  					breaks = levels(factor(row.meta[, ind])),
	                    values = palette.all.permute, 
	                    guide = guide_legend(direction = "vertical", 
                                           title.position = "top", 
                                           label.position="right")
	                    )
	            SLR[[ind]]<- gtable_filter(ggplot_gtable(ggplot_build(LR.temp)), "guide-box")
	  		}
	  		LR<-LR.temp

	  	} else {
	  		LR<-LR + scale_fill_manual(values = palette.all.permute)
	  		SLR<-list(p)
	  		names(SLR)<-"blank"
	  	}
	  	#--row side label names--
	    LRT <- gtable_filter(ggplot_gtable(ggplot_build(LR)), "axis.b") 

	    LR <-LR + theme_none + theme(legend.position = "none" )
	}
  	#make heatmap with layout
	dx1<-1/12
	if (is.na(title.text[1])){
		dx1<-0
	}

	dx2<-1.5/12
	if (col.clust ==FALSE){
		dx2<-0
	}
	dy4<-1.5/12
	if (row.clust == FALSE){
		dy4<-0
	}

	dx3<-(0.3 * length(col.lab))/12
	if (col.lab[1] == ""){
		dx3<-0
	}

	dy3 <-(0.3 * length(row.lab))/12 
	if (row.lab[1] == ""){
		dy3<-0
	}

	if (length(grid.heights) <7){
		dx <-c(dx1, dx2, (1-dx1 -dx2)* dx3 , 8/12, 3.2/12, 1.5/12,1/12)
	} else{
		dx<-grid.heights 
	}

	if (length(grid.widths)<4){
		dy<-c( (1- dy4) * c(2/12, 8/12, dy3), dy4)
	} else{
		dy<-grid.widths
	}

	p1<-arrangeGrob(p, HC, p, p, widths = dy, ncol = 4)
	p2<-arrangeGrob(LCT, LC, p, p, widths = dy, ncol = 4)
	p3<-arrangeGrob(MCT, M, LR, HR, widths = dy, ncol = 4)

	if (heatmap.x.text == FALSE){
		p4<-arrangeGrob(p, ML, LRT, p, widths = dy, ncol = 4)
		p6<-p
	} else{
		p4<-arrangeGrob(p, MRT, LRT, p, widths = dy, ncol = 4)
		p6<-ML
	}



	p5.args<-c(SLC, SLR , heights=1/2, nrow=1)
	p5<-do.call(arrangeGrob, p5.args)

	g.main<-arrangeGrob(p1,p2,p3,p4, p5, p6,heights = dx[2:7], nrow = 6)
 	g<-arrangeGrob(main.title, g.main, heights= c(dx[1], 1-dx[1]), nrow=2)

 	
 	if (show == TRUE){
 		print(g)
 	}
 	if (verbose == FALSE){
  		return(g)
  	}
  	else {
  		res<-list(heatmap = g, meta.c = meta.c, meta.r = meta.r, col.meta = col.meta, row.meta = row.meta)
  		return(res)
  	}
  
}

