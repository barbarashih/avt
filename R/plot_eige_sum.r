#' Plot the change in variance 
#'
# This function plots the output from calculateVariance 
#' @import ggplot2
#' @param avt_var List. Output from calculateVariance().
#' @param out_dir String. File path to output directory.
#' @param height Numeric. Height of the output PDF plot.
#' @param width Numeric. Width of the output PDF plot.
#'
#' @export

plotEigeSum <- function(avt_var, out_dir="avt_out", height = 6, width = 10){
	# create output directory
	dir.create(out_dir, showWarnings = FALSE)
	f_out_fp <- paste(out_dir, "/",  "Variance_analysis.pdf", sep="")
	f_plot_df <- do.call(rbind, avt_var)
	f_plot_df$Type <- as.character(f_plot_df$Type)
	f_plot_df$grouping <- sapply(strsplit(f_plot_df$Type,"_"),FUN=function(x)paste(x[2:length(x)], collapse="_"))
	f_plot_df$ranking <- sapply(strsplit(f_plot_df$Type,"_"),FUN=function(x)x[1])
	p <- ggplot(f_plot_df,aes_string(x="Removed", y = "VC")) + geom_area(alpha=0.5) + theme_bw() +
		geom_bar(stat = "identity") + xlab("Genera removed") + ylab("Impact on variance") + facet_grid(vars(f_plot_df$ranking), vars(f_plot_df$grouping))
	ggsave( f_out_fp, width = width, height = height)
}
