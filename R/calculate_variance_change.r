#' Calculate change in variance 
#'
#' this function takes the output from organised data and top genus list
#' and calculate the variance in the top genus
#' and returns the sum of the top eigen values

#' @param avt_dat Numeric vector.
#' @param top_genus Numeric vector.
#' @param out_dir Numeric vector.
#' @param top_genus_num2 Numeric vector.
#' @param top_eigen Numeric vector.
#'
#'
#' @export

calculateVariance <- function(avt_dat, top_genus, out_dir, top_genus_num2=200, top_eigen=5){
	# split the top genus by groups
	grouping <- sapply(strsplit(names(top_genus), "_"), FUN=function(x)paste(x[2:length(x)], collapse="_"))
	out_variance_list <- list()
	# Go through the top geneus in each comparison
	for(top_geneus_comoparison_idx in 1:length(top_genus)){
		current_comparison_name <- names(top_genus)[top_geneus_comoparison_idx]
		current_comparison_type <- strsplit(current_comparison_name, "_")[[1]][1]
		current_comparison_type_longname <- list( "ABUN"= "abundance",
													"PREV" = "prevalence",
													"CENT" = "centrality")
		out_fp <- paste(out_dir, "/", current_comparison_name, "_taxaConfidence_", avt_dat[["taxa_confidence"]], ".csv", sep="")
		current_grouping = grouping[top_geneus_comoparison_idx]
		print(paste("Calculating variance when top",  current_comparison_type_longname[[current_comparison_type]], "genus are removed. Group:", current_comparison_name, sep=" "))
		# filter the data with relavent samples and genus
		samples <- avt_dat[["grouping"]][[current_grouping]]
		counts <- avt_dat[["count"]][,samples]
		taxa <- avt_dat[["taxa"]]
		meta <- avt_dat[["meta"]]

		OTU_tree <- avt_dat[["otu_tree"]]
		OTU_tree <- ape::multi2di(OTU_tree)
		phyloseq_obj <- phyloseq::phyloseq(phyloseq::otu_table(counts,taxa_are_rows = T),phyloseq::tax_table(as.matrix(taxa)),meta,OTU_tree)
		taxa_sums = tapply(phyloseq::taxa_sums(phyloseq_obj), phyloseq::tax_table(phyloseq_obj)[, "Genus"], sum, na.rm=T)
		#taxa_top = names(sort(taxa_sums, TRUE))[1:top_genus_num2]
		#phyloseq_obj_base <- prune_taxa((tax_table(phyloseq_obj)[, "Genus"] %in% taxa_top), phyloseq_obj)
		phyloseq_obj_base <- phyloseq_obj
		# try to catch ordinate error
		ordinate_error <- FALSE
		tryCatch(phyloseq::ordinate(physeq=phyloseq_obj_base,  "PCoA", "unifrac", weighted=TRUE), error=function(e){ ordinate_error <<- TRUE; print("Ordinate step failed when no genius was removed")})
		
		if (ordinate_error){
			print("There may be too few mapped reads for the ordinate step to be carried out.")
			#print(colSums(phyloseq::otu_table(phyloseq_obj_base)))
			next
		} else {
			ordinate_all = phyloseq::ordinate(physeq=phyloseq_obj_base,  "PCoA", "unifrac", weighted=TRUE)
			ordinate_top_eig_nofilt <- sum(ordinate_all$values$Relative_eig[1:top_eigen])
			base_ordinate_df <- data.frame(Genus="All", Percentage = ordinate_top_eig_nofilt, Removed= 0, stringsAsFactors=FALSE)
			# remove the top genus step-wise
			top_genus_list <- lapply(1:nrow(top_genus[[top_geneus_comoparison_idx]]), FUN=function(x)top_genus[[top_geneus_comoparison_idx]]$Genus[1:x])
			# Loop through the inclusion criteria from least number of excluded genes to most
			skip_loop <- FALSE
			all_ordinate_list <- list()
			for (genus_list_idx in 1:length(top_genus_list)){
				current_excluding_genus <- top_genus_list[[genus_list_idx]]
				taxa_filtered <- as.character(taxa$Genus[!(taxa$Genus %in% current_excluding_genus)])
				phyloseq_obj_filtered = phyloseq::prune_taxa((phyloseq::tax_table(phyloseq_obj_base)[, "Genus"] %in% taxa_filtered), phyloseq_obj_base)
				taxa_sums_filtered = tapply(phyloseq::taxa_sums(phyloseq_obj_filtered), phyloseq::tax_table(phyloseq_obj_filtered)[, "Genus"], sum, na.rm=T)
				# keep only the top top_genus2 genus with the highest count sum
				#top_genus2 = names(sort(taxa_sums_filtered, TRUE))[1:top_genus_num2]
				#taxa_filatered2 <- prune_taxa((tax_table(phyloseq_obj_base)[, "Genus"] %in% top_genus2), phyloseq_obj_base) # use taxa_filatered2 instead of phyloseq_obj_filtered if filtering by top_genus_num2 
				ordinate_error <- FALSE
				tryCatch(phyloseq::ordinate(phyloseq_obj_filtered, "PCoA", "unifrac", weighted=TRUE), error=function(e){ ordinate_error <<- TRUE; print(paste("ordinate step failed when the top",  length(top_genus_list[[genus_list_idx]]), "genus are removed", sep="  "))})
				if(ordinate_error){ 
					break
				} else {
					ordinate_filatered <- phyloseq::ordinate(phyloseq_obj_filtered, "PCoA", "unifrac", weighted=TRUE)
					all_ordinate_list[[genus_list_idx]] <- data.frame(Genus= paste(top_genus_list[[genus_list_idx]], collapse=","),  Percentage=sum(ordinate_filatered$values$Relative_eig[1:top_eigen]), Removed=length(top_genus_list[[genus_list_idx]]))
				}
			}
			all_ordinate_df <- do.call(rbind,all_ordinate_list)
			current_top_genus_eigen_sum <- rbind(base_ordinate_df, all_ordinate_df)
			out_variance_list[[current_comparison_name]] <- current_top_genus_eigen_sum
			utils::write.csv(current_top_genus_eigen_sum, out_fp, row.names=FALSE)		
		}
	}
	out_df <- lapply(out_variance_list, FUN=function(x){out_x <- x; out_x$VC = out_x$Percentage[1] - out_x$Percentage; return(out_x)})
	out_df <- lapply(names(out_df), FUN=function(x){out_x <- out_df[[x]]; out_x$Type = x; return(out_x)})
	names(out_df) <- names(out_variance_list)
	return(out_df)
}
