
library(optparse)
library(ggplot2)

# Parse the IO files from the CLI
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to the input GO contact list file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Path for the output contact matrix plot in PDF", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output)){
  print_help(opt_parser)
  stop("Required -i input and -o output file paths.\n", call.=FALSE)
}

golist = read.csv(opt$input, header = F, sep = "")
names(golist) = c("i", "j", "distance", "energy")

golist_rev = golist
golist_rev$x = golist_rev$i
golist_rev$i = golist_rev$j
golist_rev$j = golist_rev$x
golist_rev$x = NULL

golist = rbind(golist, golist_rev)

# Contact matrix plot
p = ggplot(golist, aes(x = i, y = j, fill=-energy)) + 
  geom_tile() + 
  theme_bw() + theme(panel.grid = element_blank(), 
                     legend.position = "top") + 
  xlab("") + ylab("") +
  scale_fill_continuous(name = "Go energy (kcal/mol)")

pdf(opt$output, 5, 5.5)
print(p)
dev.off()