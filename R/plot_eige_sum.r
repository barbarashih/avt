#' Plot the change in variance 
#'
# This function plots the output from calculateVariance 
#' @import ggplot2
#' @param avt_var Numeric vector.
#' @param out_dir Numeric vector.
#' @param height Numeric vector.
#' @param width Numeric vector.
#'
#' @export

plotEigeSum <- function(avt_var, out_dir, height = 6, width = 10){
	f_out_fp <- paste(out_dir, "/",  "Variance_analysis.pdf", sep="")
	f_plot_df <- do.call(rbind, avt_var)
	f_plot_df$Type <- as.character(f_plot_df$Type)
	f_plot_df$grouping <- sapply(strsplit(f_plot_df$Type,"_"),FUN=function(x)paste(x[2:length(x)], collapse="_"))
	ggplot(f_plot_df,aes(Removed,VC)) + geom_area(alpha=0.5) + theme_bw() +
		geom_bar(stat = "identity") + xlab("Genera removed") + ylab("Impact on variance") + facet_wrap(~Type, ncol=length(unique(f_plot_df$grouping)))
	ggsave( f_out_fp, width = width, height = height)
}
