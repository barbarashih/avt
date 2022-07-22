#' Calculate change in variance 
#'
#' this function takes the output from organised data and top genus list
#' and calculate the variance in the top genus
#' and returns the sum of the top eigen values

#' @param avt_dat List. Output from avt_data().
#' @param top_genera List. Output from topGenera().
#' @param out_dir String. File path to output directory.
#' @param top_genera_num2 Numeric. Number of genus to consider. Default set to 200.
#' @param top_eigen Numeric. Number of top eigen vectors to consider. Default set to 10.
#' @param seed Numeric. 
#'
#'
#' @export

calculateVariance <- function(avt_dat, top_genera, out_dir = "avt_out", top_genera_num2=200, top_eigen=10, seed=123){
	# create output directory
	dir.create(out_dir, showWarnings = FALSE)
	# split the top genus by groups
	grouping <- sapply(strsplit(names(top_genera), "_"), FUN=function(x)paste(x[2:length(x)], collapse="_"))
	out_variance_list <- list()
	# Go through the top geneus in each comparison
	all_comparisons = names(top_genera)
	all_comparisons = all_comparisons[order(all_comparisons)]
	for(top_geneus_comoparison_idx in 1:length(all_comparisons)){
		current_comparison_name <- all_comparisons[top_geneus_comoparison_idx]
		current_comparison_type <- strsplit(current_comparison_name, "_")[[1]][1]
		current_comparison_type_longname <- list( "ABUN"= "abundance",
													"PREV" = "prevalence",
													"CENT" = "centrality")
		out_fp <- paste(out_dir, "/", current_comparison_name, "_taxaConfidence_", avt_dat[["taxa_confidence"]], ".csv", sep="")
		current_grouping = grouping[top_geneus_comoparison_idx]
		#print(paste("Current grouop:", current_comparison_name, sep=" "))
		# filter the data with relavent samples and genus
		samples <- avt_dat[["grouping"]][[current_grouping]]
		counts <- avt_dat[["count"]][,samples]
		taxa <- avt_dat[["taxa"]]
		meta <- avt_dat[["meta"]]

		OTU_tree <- avt_dat[["otu_tree"]]
		OTU_tree <- ape::multi2di(OTU_tree)
		phyloseq_obj <- phyloseq::phyloseq(phyloseq::otu_table(counts,taxa_are_rows = T),phyloseq::tax_table(as.matrix(taxa)),meta,OTU_tree)
		taxa_sums = tapply(phyloseq::taxa_sums(phyloseq_obj), phyloseq::tax_table(phyloseq_obj)[, "Genus"], sum, na.rm=T)
		#taxa_top = names(sort(taxa_sums, TRUE))[1:top_genera_num2]
		#phyloseq_obj_base <- prune_taxa((tax_table(phyloseq_obj)[, "Genus"] %in% taxa_top), phyloseq_obj)
		phyloseq_obj_base <- phyloseq_obj
		# try to catch ordinate error
		ordinate_error <- FALSE
		set.seed(seed)
		tryCatch(phyloseq::ordinate(physeq=phyloseq_obj_base,  "PCoA", "unifrac", weighted=TRUE), error=function(e){ ordinate_error <<- TRUE; print(paste(current_comparison_name, ": ordinate step failed when no genius was removed.", sep=" "))})
		
		if (ordinate_error){
			#print("There may be too few mapped reads for the ordinate step to be carried out.")
			#print(colSums(phyloseq::otu_table(phyloseq_obj_base)))
			next
		} else {
			set.seed(seed)
			ordinate_all = phyloseq::ordinate(physeq=phyloseq_obj_base,  "PCoA", "unifrac", weighted=TRUE)
			ordinate_top_eig_nofilt <- sum(ordinate_all$values$Relative_eig[1:top_eigen])
			base_ordinate_df <- data.frame(Genus="All", Percentage = ordinate_top_eig_nofilt, Removed= 0, stringsAsFactors=FALSE)
			# remove the top genus step-wise
			top_genera_list <- lapply(1:nrow(top_genera[[top_geneus_comoparison_idx]]), FUN=function(x)top_genera[[top_geneus_comoparison_idx]]$Genus[1:x])
			# Loop through the inclusion criteria from least number of excluded genes to most
			skip_loop <- FALSE
			all_ordinate_list <- list()
			for (genus_list_idx in 1:length(top_genera_list)){
				current_excluding_genus <- top_genera_list[[genus_list_idx]]
				taxa_filtered <- as.character(taxa$Genus[!(taxa$Genus %in% current_excluding_genus)])
				phyloseq_obj_filtered = phyloseq::prune_taxa((phyloseq::tax_table(phyloseq_obj_base)[, "Genus"] %in% taxa_filtered), phyloseq_obj_base)
				taxa_sums_filtered = tapply(phyloseq::taxa_sums(phyloseq_obj_filtered), phyloseq::tax_table(phyloseq_obj_filtered)[, "Genus"], sum, na.rm=T)
				# keep only the top top_genera2 genus with the highest count sum
				#top_genera2 = names(sort(taxa_sums_filtered, TRUE))[1:top_genera_num2]
				#taxa_filatered2 <- prune_taxa((tax_table(phyloseq_obj_base)[, "Genus"] %in% top_genera2), phyloseq_obj_base) # use taxa_filatered2 instead of phyloseq_obj_filtered if filtering by top_genera_num2 
				ordinate_error <- FALSE
				tryCatch(phyloseq::ordinate(phyloseq_obj_filtered, "PCoA", "unifrac", weighted=TRUE), error=function(e){ ordinate_error <<- TRUE; print(paste(current_comparison_name, ": ordinate step failed when the top",  length(top_genera_list[[genus_list_idx]]), "genera are removed", sep=" "))})
				if(ordinate_error){ 
					break
				} else {
					ordinate_filatered <- phyloseq::ordinate(phyloseq_obj_filtered, "PCoA", "unifrac", weighted=TRUE)
					all_ordinate_list[[genus_list_idx]] <- data.frame(Genus= paste(top_genera_list[[genus_list_idx]], collapse=","),  Percentage=sum(ordinate_filatered$values$Relative_eig[1:top_eigen]), Removed=length(top_genera_list[[genus_list_idx]]))
				}
			}
			if(length(top_genera_list[[genus_list_idx]])  == length(top_genera_list)){
				print(paste(current_comparison_name, ": variance calculation up to",  length(top_genera_list[[genus_list_idx]]), "genera removed", sep=" "))
			}
			all_ordinate_df <- do.call(rbind,all_ordinate_list)
			current_top_genera_eigen_sum <- rbind(base_ordinate_df, all_ordinate_df)
			out_variance_list[[current_comparison_name]] <- current_top_genera_eigen_sum
			utils::write.csv(current_top_genera_eigen_sum, out_fp, row.names=FALSE)		
		}
	}
	out_df <- lapply(out_variance_list, FUN=function(x){out_x <- x; out_x$VC = out_x$Percentage[1] - out_x$Percentage; return(out_x)})
	out_df <- lapply(names(out_df), FUN=function(x){out_x <- out_df[[x]]; out_x$Type = x; return(out_x)})
	names(out_df) <- names(out_variance_list)
	return(out_df)
}
