# Plot the Go contact energy matrix of a protein
# 
# Aleix Lafita - April 2018

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))

contmat_theme = theme_bw() +theme(
  panel.grid.minor = element_blank(),
  #axis.text.x = element_text(angle = 90, hjust = 1),
  legend.position = "right",
  axis.title = element_blank()
)
theme_set(contmat_theme)

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

golist = rbind(golist, golist_rev) %>%
  mutate(energy = -energy)

# Contact matrix plot
p = ggplot(golist, aes(x = i, y = j, fill=energy)) + 
  geom_tile() +
  geom_abline(slope = -1, alpha = 0.5) +
  scale_fill_continuous(name = "kcal/mol") +
  scale_y_reverse(
    expand = c(0, 0)
  ) + scale_x_continuous(
    expand = c(0, 0)
  ) + coord_fixed()

pdf(opt$output, 6, 5)
print(p)
log = dev.off()
