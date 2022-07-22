#' Summarise the change in variance 
#'
# This function prints out a table with ordered variance change within each grouping
#' @param avt_var List. Output from calculateVariance().
#' @param out_dir String. File path to output directory.
#'
#' @export

orderedVarianceChange <- function(avt_var, out_dir="avt_out"){
	# create output directory
	dir.create(out_dir, showWarnings = FALSE)
	# function for calculating the difference in vairance
	fun_var_diff <- function(current_df){
		# calculate the variancce change between pre- and post- removal of each genus
		# calculate the variance change if there are two or more entries for genus removal
		if(nrow(current_df)>1){
			out_diff <- NA
			# make sure it's in order
			current_df <- current_df[order(current_df$Removed),]
			# calculate the change in variance between before/after genus removal
			for(i in 2:nrow(current_df)){
				current_diff <- current_df$VC[i] - current_df$VC[(i-1)]
				out_diff <- c(out_diff, abs(current_diff))
			}
			out_df <- current_df
			out_df$variance_change <- out_diff
			out_df$genus_removed <- sapply(strsplit(out_df$Genus, ","), FUN=function(x)x[length(x)])
			out_df <- out_df[order(out_df$variance_change, decreasing=TRUE),]
			
			return(out_df[,c("genus_removed", "Type", "variance_change", "Removed")])
		# If there are only one or less rows for the variance analysis, return NA
		} else if (nrow(current_df) == 1) {
			out_df <- current_df
			out_df$variance_change <- NA
			out_df$genus_removed <- sapply(strsplit(out_df$Genus, ","), FUN=function(x)x[length(x)])
			return(out_df[,c("genus_removed", "Type", "variance_change", "Removed")])
		}
	}
	# organise the output, remove NAs 
	out_ordered <- lapply(avt_var, FUN=function(x)fun_var_diff(x))
	out_ordered <- do.call(rbind, out_ordered)
	out_ordered <- out_ordered[!is.na(out_ordered$variance_change),]
	# organise output columns
	out_ordered$ranking_method <- sapply(strsplit(out_ordered$Type,"_"),FUN=function(x)x[1])
	out_ordered$grouping <- sapply(strsplit(out_ordered$Type,"_"),FUN=function(x)paste(x[2:length(x)], collapse="_"))
	out_ordered$Type = NULL
	row.names(out_ordered) <- 1:nrow(out_ordered)
	# write output
	utils::write.csv(out_ordered, paste0(out_dir, "/variance_change_summary.csv"), row.names=FALSE)		
	return(out_ordered)
}
