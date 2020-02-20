# Plot a profile of the estimated free energy difference along domain sequence
# 
# Aleix Lafita - April 2018

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))

# Parse the file names from the CLI
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to the input dG_cut.tsv file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Path for the output plot in PDF", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output)){
  print_help(opt_parser)
  stop("Required -i input and -o output file paths.\n", call.=FALSE)
}

data = read.csv(opt$input, sep="\t", header=T)

p = ggplot(data, aes(y=dGcut, x=position)) + 
  geom_step(stat = "identity") + 
  geom_hline(yintercept = 0, alpha = 0.4)+
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(
    #limits = c(1, max(data$position) + min(data$position) - 1),
    breaks = c(1, seq(10, max(data$position) + min(data$position) - 1, 10))
  ) +
  xlab("Position") + ylab(expression(paste(Delta,"Gc (kcal/mol)", sep="")))

pdf(opt$output, 5, 3)
print(p)
log = dev.off()
