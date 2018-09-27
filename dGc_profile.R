
library(optparse)
library(ggplot2)

# Parse the IO files from the CLI
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to the input dG_cut.tsv file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Path for the output profile plot in PDF", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output)){
  print_help(opt_parser)
  stop("Required -i input and -o output file paths.\n", call.=FALSE)
}

data = read.csv(opt$input, sep="\t", header=T)

p = ggplot(data, aes(y=dGcut, x=position)) + 
  geom_line(stat = "identity") + 
  geom_hline(yintercept = 0) +
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(min(data$position), max(data$position), 10)) +
  xlab("Position") + ylab(expression(paste(Delta,"Gc (kcal/mol)", sep="")))

pdf(opt$output, 4, 3)
print(p)
log = dev.off()
