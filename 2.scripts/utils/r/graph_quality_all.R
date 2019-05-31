library("optparse")
library( "RColorBrewer")

colors <- c( rep( brewer.pal( 5, 'RdYlGn')[ 1], 2), 
               brewer.pal( 5, 'Greens')[ -1])

#==========================================================================
# Parsing argument to get data
#==========================================================================


## Collect arguments
args <- commandArgs(TRUE)

option_list = list(
  make_option( c( "-f", "--file"), type="character", default=NULL,
               help="path to tab file create by snakemake (experiment_name	NSC	RSC	FRiP	nb_peaks	score_NSC	score_RSC	score_FRiP	score_total)", metavar="character"),
  make_option( c( "-o", "--outfile"), type="character", default=NULL,
                help="name to outfile", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#==========================================================================
# Function
#==========================================================================

## Function to create a graph with all data
large.graph.score <- function(df.score, graph.title) {
  
  plot(df.score[,'NSC'], 
       df.score[,'RSC'], 
       xlim = c(1,8),
       ylim = c(0,10), 
       pch = 20, 
       col = colors[ df.score[,"score_total"]+1], 
       main = graph.title, 
       xlab = 'NSC', 
       ylab = 'RSC')
  
  abline(h = 0.8,
         col='grey')
  
  abline(h = 1.0,
         col='grey')
  
  abline(v = 1.10,
         col='grey')
  
  abline(v = 1.05,
         col='grey')
  
  legend('topright',
         legend = 0:5,
         fill = colors,
         title = 'Score', 
         bty='n')
}

## Function to create a graph with data near origins of graph
zoom.graph.score <- function(df.score, graph.title) {
  
  plot(df.score[,'NSC'],
       df.score[,'RSC'], 
       xlim = c(1,2), 
       ylim = c(0,2), 
       pch = 20, 
       col = colors[df.score[,"score_total"]+1], 
       main = graph.title, 
       xlab = 'NSC', 
       ylab = 'RSC')
  
  abline(h = 0.8,
         col = 'grey')
  
  abline(h = 1.0,
         col = 'grey')
  
  abline(v = 1.10,
         col = 'grey')
  
  abline(v = 1.05,
         col = 'grey')
  
  legend('bottomright',
         legend = 0:5,
         fill = colors,
         title = 'Score', 
         bty = 'n')
  
}


#==========================================================================
# Begining of script
#==========================================================================

# importing data
df.quality.all = read.table( opt$file,
                              header=TRUE,
                              sep = "\t",
                              fill = TRUE
                              )

# df.quality.all <- read.table( "/home/cheneby/sacapus/tagc-remap3/src/jeanne/8.quality/results/macs2.quality_all",
#                               header=TRUE,
#                               sep = "\t",
#                               fill = TRUE
#                             )
                              
                              


pdf( opt$outfile)
  large.graph.score( df.quality.all,
                     "MACS2 all datasets quality")
  
  zoom.graph.score( df.quality.all,
                   "MACS2 all datasets quality")


dev.off()