#' Organise input data into AVT format
#'
#' This tool tracks microbiome variance and identifies bacterial genera responsible for significant variance shifts using 16S rRNA gene metabarcoding for longitudinal data. 
#' @import dplyr
#' @param data_dir Numeric vector.
#' @param grouping_column Numeric vector.
#' @param out_dir Numeric vector.
#' @param count_filename Numeric vector.
#' @param taxa_filename Numeric vector.
#' @param meta_filename Numeric vector.
#' @param tree_filename Numeric vector.
#' @param f_taxaconfidence Numeric vector.
#' @param f_grouping_inclusion Numeric vector.
#' @param min_req_sample_number Numeric vector.
#'
#' @export

avt_data <- function(data_dir, grouping_column, out_dir ,count_filename="OTU.txt", taxa_filename="taxonomy.tsv", meta_filename="MetaDB.txt", tree_filename="tree.nwk", f_taxaconfidence=0.9, f_grouping_inclusion=vector(), min_req_sample_number =8){
	 
	# read in data
	f_study_name=basename(data_dir)
	dir.create(out_dir, showWarnings = FALSE)

	f_taxa <- utils::read.delim(paste(data_dir,taxa_filename,sep="/"),  header = T, sep="\t")
	# if taxa table is not already formatted, format it so that it contains the genus column
	if(!("Kingdom" %in% colnames(f_taxa))){
		row.names(f_taxa) <- f_taxa$Feature.ID
		f_taxa$Taxon <- as.character(f_taxa$Taxon)
		f_taxa$Kingdom <- gsub(";", "", gsub("D_0__", "", stringr::str_extract(paste0(f_taxa$Taxon, ";"), "D_0__(.*?);.*?")))
		f_taxa$Phylum <- gsub(";", "", gsub("D_1__", "", stringr::str_extract(paste0(f_taxa$Taxon, ";"), "D_1__(.*?);.*?")))
		f_taxa$Class <- gsub(";", "", gsub("D_2__", "", stringr::str_extract(paste0(f_taxa$Taxon, ";"), "D_2__(.*?);.*?")))
		f_taxa$Order <- gsub(";", "", gsub("D_3__", "", stringr::str_extract(paste0(f_taxa$Taxon, ";"), "D_3__(.*?);.*?")))
		f_taxa$Family <- gsub(";", "", gsub("D_4__", "", stringr::str_extract(paste0(f_taxa$Taxon, ";"), "D_4__(.*?);.*?")))
		f_taxa$Genus <- gsub(";", "", gsub("D_5__", "", stringr::str_extract(paste0(f_taxa$Taxon, ";"), "D_5__(.*?);|D_5__(.*?)$]")))
		f_taxa$Genus  <- ifelse(f_taxa$Genus %in% c("", " "), NA, f_taxa$Genus)
		# Given Genus are sometimes messy - this bit cleans it up
		# find the Genus labels that are 1 word + space + number
		f_taxa$Genus_word_count <- sapply(strsplit(f_taxa$Genus, " "), function(x)length(x))
		f_taxa$Genus_word_and_num <- stringr::str_detect(f_taxa$Genus,"[:alpha:] [:digit:]")
		f_taxa$Genus_word_first <- sapply(strsplit(f_taxa$Genus, " "), function(x)x[1])
		f_taxa$Genus <- ifelse((f_taxa$Genus_word_count ==2) & f_taxa$Genus_word_and_num, f_taxa$Genus_word_first, f_taxa$Genus)
		f_taxa <- cleanGenus(f_taxa)
		f_taxa$Genus <- ifelse(f_taxa$Confidence < f_taxaconfidence, "Low_confidence", f_taxa$Genus)
		f_taxa <- f_taxa[!(f_taxa$Genus %in% "Low_confidence"),]
		f_taxa <- f_taxa[!is.na(f_taxa$Genus),]
		row.names(f_taxa) <- f_taxa$Feature.ID
		# print out the number of read in each taxa 
		genus_count <- stats::aggregate(data=f_taxa, Feature.ID~Genus, length)
		genus_count_out_fp <- paste(out_dir, "/", "taxa_used_with_taxaConfidence", f_taxaconfidence, ".csv", sep="")
		utils::write.csv(genus_count, genus_count_out_fp, row.names=FALSE)
		
		f_taxa <- f_taxa[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")]
	}
	
	# read in each of the input file
	f_counts <- utils::read.delim(paste(data_dir,count_filename,sep="/"), sep="\t", row.names=1)
	f_meta<- utils::read.delim(paste(data_dir,meta_filename,sep="/"), header = T, comment.char="#", row.names=1)
	f_meta <- f_meta[row.names(f_meta) %in% colnames(f_counts),colnames(f_meta), drop=FALSE]
	f_tree <- phyloseq::read_tree(paste(data_dir,tree_filename, sep="/"),errorIfNULL=FALSE)
	# remove any line that starts with # in the meta data
	
	
	# Convert objects for phyloseq object
	f_otu <- phyloseq::otu_table(f_counts,taxa_are_rows = T)
	f_meta<-phyloseq::sample_data(f_meta)
	f_OTU_tree <- ape::compute.brlen(f_tree, method = "Grafen")

	# Create Phyloseq experiment object
	f_EXPname <- phyloseq::phyloseq(f_otu,f_taxa,f_meta,f_OTU_tree)
	f_OTU <- phyloseq::otu_table(as.matrix(f_otu), taxa_are_rows = FALSE)
	# Sample groups
	if(length(f_grouping_inclusion)>0){
		f_meta <- f_meta[f_meta[[grouping_column]] %in% f_grouping_inclusion,]
	}
	f_sample_groups <- split(row.names(f_meta), f_meta[[grouping_column]])
	for(group_idx in 1:length(f_sample_groups)){
		if(length(f_sample_groups[[group_idx]]) < min_req_sample_number){
			print(paste("Not enough observations in group - ", names(f_sample_groups)[group_idx], " - at least ", min_req_sample_number, " samples required per group.", sep=""))
		}
	}
	f_sample_groups <- f_sample_groups[lapply(f_sample_groups,length)>=min_req_sample_number]
	if(length(f_sample_groups) < 2){
		print("Less than 2 groups with enough samples. The pipeline cannot be performed.")
	} else {
		out_list <- list(
				study_name=f_study_name,
				count =f_counts[row.names(f_counts) %in% row.names(f_taxa),],
				phyloseq_count = f_EXPname,
				taxa = f_taxa,
				meta = phyloseq::sample_data(f_meta),
				otu = f_OTU,
				otu_tree = f_OTU_tree,
				grouping=f_sample_groups,
				taxa_confidence=f_taxaconfidence
			)
		return(out_list)
	}
}
