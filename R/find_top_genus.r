#' Determine the top Genus (by abundance, prevalence or centrality)
#'
# This function takes in the list output from the avt.data_organisation
# and returns the top genus for abundance, prevalence and centrality for each sample grouping in the given input
#' @param avt_dat Numeric vector.
#' @param out_dir Numeric vector.
#' @param top_number Numeric vector.
#' @param centrality_p Numeric vector.
#' @param centrality_r_upper Numeric vector.
#' @param centrality_r_lower Numeric vector.
#'
#' @export

topGenus <- function(avt_dat, out_dir, top_number=50, centrality_p=0.05, centrality_r_upper=0.7, centrality_r_lower=0){
	# define output directory
	dir.create(out_dir, showWarnings = FALSE)

	# go through each of the grouping
	sample_grouping <- avt_dat[["grouping"]]
	count <- avt_dat[["count"]]
	taxa <- avt_dat[["taxa"]]
	count$OUT_ID <- row.names(count)
	taxa$OUT_ID <- row.names(count)
	print(paste("Total number of features:", nrow(count), sep=" "))
	print(paste("Total number of taxonomy:", nrow(taxa), sep=" "))
	count_annotated <- merge(count, taxa, by="OUT_ID")
	count_annotated <- count_annotated[,!(colnames(count_annotated) %in% c("Kingdom","Phylum","Class","Order","Family"))]
	top_genus <- list()
	# Get the top X genus for each grouping using each method of assessment (abundance, prevalence and centrality) and save the results in a list
	for(grouping_idx in 1:length(sample_grouping)){
		current_grouping_name <- names(sample_grouping)[grouping_idx]
		print(paste("Processing grouping:", current_grouping_name, sep=" ")) 
		current_count <- count_annotated[,colnames(count_annotated) %in% c(sample_grouping[[grouping_idx]], "Genus")] 
		print(paste("current_count nrow " , nrow(current_count), sep=""))
		if(nrow(current_count) > 0){
			current_count_mean_genus <- stats::aggregate( .~Genus, data=current_count, mean)
			current_count_mean_genus <- current_count_mean_genus[!is.na(current_count_mean_genus$Genus),]
			rownames(current_count_mean_genus) <- current_count_mean_genus$Genus
			print("Top abundance")
			top_genus[[paste("ABUN", current_grouping_name, sep="_")]] = topByPrevalence(mean_genus_count = current_count_mean_genus, top_number = top_number)
			print("Top prevalence")
			top_genus[[paste("PREV", current_grouping_name, sep="_")]] = topByPrevalence(mean_genus_count = current_count_mean_genus, top_number = top_number)
			print("Top centrality")
			if(TRUE){
				top_genus[[paste("CENT", current_grouping_name, sep="_")]] = topByCentrality(mean_genus_count = current_count_mean_genus, 
																								taxa = taxa, 
																								out_dir = out_dir, 
																								current_grouping_name = current_grouping_name,
																								top_number = top_number ,
																								centrality_p = centrality_p, 
																								centrality_r_upper = centrality_r_upper, 
																								centrality_r_lower = centrality_r_lower)

			}
		}
	}
	# calculate the variance
	return(top_genus)
	
}

#' Find the top Genus by abundance
#'
#' Find the top Genus by abundance
#' @param mean_genus_count Numeric vector.
#' @param top_number Numeric vector.
#'
#' @export

topByAbundance <- function(mean_genus_count, top_number){
	mean_genus_count_long <- reshape2::melt(mean_genus_count)
	colnames(mean_genus_count_long) <- c("Genus", "sample_id", "value")
	out <- plyr::ddply(mean_genus_count_long, .(Genus), summarize, base_mean_genus=(mean(value, na.rm = TRUE)))
	out <- out[order(out$base_mean_genus,decreasing = T),]
	out <- utils::head(out, top_number)
	names(out)<- c("Genus","mean_genus")
	out <- out[out$mean_genus > 0,]
	return(out)
}


#' Find the top Genus by prevalence
#'
#' Find the top Genus by prevalence
#' @param mean_genus_count Numeric vector.
#' @param top_number Numeric vector.
#'
#' @export

topByPrevalence <- function(mean_genus_count, top_number){
	mean_genus_count_long<- reshape2::melt(mean_genus_count)
	colnames(mean_genus_count_long) <- c("Genus", "sample_id", "value")
	total_Base <- length(unique(mean_genus_count_long$sample_id))
	mean_genus_count_long %>%
	  filter(value != 0) %>%
	  group_by(Genus) %>%
	  tally() %>%
	  mutate(prev = n/total_Base) %>%
	  arrange(-prev) -> out
	out <- out[order(out$prev,decreasing = T),]
	out <- utils::head(out, top_number)
	names(out)<- c("Genus","n","mean_genus")
	out <- out[out$mean_genus > 0,]
	return(out)
}


#' Find the top Genus by centrality
#'
#' Find the top Genus by centrality
#' @param mean_genus_count Numeric vector.
#' @param taxa Numeric vector.
#' @param out_dir Numeric vector.
#' @param current_grouping_name Numeric vector.
#' @param top_number Numeric vector.
#' @param centrality_p Numeric vector.
#' @param centrality_r_upper Numeric vector.
#' @param centrality_r_lower Numeric vector.
#'
#' @export

topByCentrality <- function(mean_genus_count, taxa, out_dir, current_grouping_name, top_number, centrality_p, centrality_r_upper, centrality_r_lower){
	# organise output graphs file paths
	out_dendrogram_fp <- paste(out_dir, "/", "dendrogram_", current_grouping_name, ".pdf", sep="")
	out_genera_co_occurance_fp <- paste(out_dir, "/", "genera_co_occurance_", current_grouping_name, ".pdf", sep="")
	# Pearson correlation
	mean_genus_count <- as.matrix(mean_genus_count[,!colnames(mean_genus_count) %in% c("Genus")])
	correlation_table <- microbiome::associate(t(mean_genus_count), t(mean_genus_count), method = "spearman", mode = "table", p.adj.threshold = centrality_p, n.signif = 1)
	correlation_table<-as.data.frame(correlation_table)
	print("pearson-corr done")
	names(correlation_table) <- c("from","to","Correlation","p.adj")
	correlation_table <- correlation_table[correlation_table$Correlation<1,]
	correlation_table <- correlation_table[!duplicated(correlation_table[c(1:2)]),]
	# Filter with p-value and r value thresholds
	correlation_table$SIG[correlation_table$p.adj<centrality_p] <- "YES"
	correlation_table$SIG[is.na(correlation_table$SIG)] <- "NO"
	correlation_table <- correlation_table[correlation_table$SIG == "YES",]
	correlation_table <-correlation_table[((correlation_table$Correlation>centrality_r_upper) | (correlation_table$Correlation<centrality_r_lower)),] ## Filter
	correlation_table$SIG<-as.factor(correlation_table$SIG)
	correlation_table$SigColors[correlation_table$SIG=="YES"]<-"#FF0000"
	correlation_table$SigColors[correlation_table$SIG=="NO"]<-"#808000"
	# Label phylum
	correlation_table$Phylum <-taxa[match(correlation_table$from, taxa$Genus),]$Phylum
	correlation_table$Phylum<-as.factor(correlation_table$Phylum)
	# Draw graphs

	correlation_graph =  igraph::graph_from_data_frame(correlation_table, directed = FALSE, vertices = NULL)

	if(TRUE){
		comm_base <- igraph::cluster_edge_betweenness(correlation_graph)
		comm_base <- igraph::edge.betweenness.community(correlation_graph)
		hclust_error <- FALSE
		tryCatch(dendPlot(comm_base, mode="hclust"), error=function(e){ hclust_error <<- TRUE; print("denPlot failed")})
		if(! hclust_error){ 
			dendPlot(comm_base, mode="hclust")
			grDevices::pdf(out_dendrogram_fp,width = 10,height = 21)
			plot_dendrogram(comm_base,mode = "phylo", main = paste("Communities identified for", current_grouping_name, sep=" "))
			grDevices::dev.off()
		}
		grDevices::pdf(out_genera_co_occurance_fp,width = 8,height = 8)
		plot(comm_base,correlation_graph, 
			 
			 # === vertex
			 layout=igraph::layout_with_fr,
			 vertex.color = "blue",          				# Node color
			 vertex.frame.color = "white",                 # Node border color
			 vertex.shape="circle",                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
			 vertex.size=5,                               # Size of the node (default is 15)
			 # vertex.size2=NA,                              # The second size of the node (e.g. for a rectangle)
			 
			 # === vertex label
			 # Character vector used to label the nodes
			 vertex.label.color="black",
			 vertex.label.family="Times",                  # Font family of the label (e.g.“Times”, “Helvetica”)
			 vertex.label.font=0.9,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
			 vertex.label.cex=0.3,                           # Font size (multiplication factor, device-dependent)
			 vertex.label.dist=0.2,                          # Distance between the label and the vertex
			 vertex.label.degree=0.2 ,                       # The position of the label in relation to the vertex (use pi)
			 
			 # === Edge
			 edge.color=correlation_table$SigColors,                           # Edge color
			 edge.width=(correlation_table$Correlation),                                 # Edge width, defaults to 1
			 edge.arrow.size=0.1,                            # Arrow size, defaults to 1
			 edge.arrow.width=1.5,                           # Arrow width, defaults to 1
			 edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
			 edge.curved=0.6                         # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
			 
		)


		 graphics::legend(x=-0.9, y=-1, unique(correlation_table$SIG), pch=21,
					 col="#777777", pt.bg=unique(correlation_table$SigColors), pt.cex=2, cex=.8, bty="n", ncol=1)
		 
		graphics::title(paste("Microbial communities for", current_grouping_name, sep=" "),cex.main=1,col.main="black")
		grDevices::dev.off()
	}
	print("before betweenness")
	centrality <- as.data.frame(igraph::betweenness(correlation_graph))
	centrality$Genus <- row.names(centrality)
	row.names(centrality) <- NULL
	centrality <- dplyr::rename(centrality, "betweeness" = `igraph::betweenness(correlation_graph)`)
	centrality <- dplyr::select(centrality, Genus, betweeness)
	out <- centrality[order(centrality$betweeness,decreasing = T),]
	out <- utils::head(out , top_number)
	out <-as.data.frame(out )
	names(out)<- c("Genus","mean_genus")
	return(out)
}


