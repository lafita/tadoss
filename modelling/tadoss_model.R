# Generate tandem domain swap models from a TADOSS run
# Runs external scripts based on edmc distance matrix transforms
# 
# Aleix Lafita - July 2020

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(bio3d))
suppressPackageStartupMessages(library(edmcr))

######################## Parameters #########################

domain = "sh3"
path = "~/tadoss"

# Parse the file names from the CLI
option_list = list(
  make_option(c("-p", "--path"), type="character", default=path,
              help="Path to the TADOSS source code directory", metavar="character"),
  make_option(c("-d", "--domain"), type="character", default=domain,
              help="Directory and domain name for input dG_cut, dG_join and output PDB files", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

########################### Parsing ############################

# Cut ddG estimates
cut = read.csv(
  sprintf("%s-dG_cut.tsv", opt$domain), 
  sep="\t", 
  header=T
)

# Join ddG estimate
join = read.csv(
  sprintf("%s-dG_join.tsv", opt$domain), 
  sep="\t", 
  header=T
)

########################### DS indices ############################

# Index position of the swap
max_cut = cut %>% top_n(1, dGcut)

cmd = sprintf(
  "Rscript %s/modelling/model_tswap.R -s tmp/%s_hadded.pdb -b %s/modelling/calpha_bounds.tsv -o %s_tswap%i.pdb -i %i -g %i -f %i",
  opt$path,
  opt$domain,
  opt$path,
  opt$domain,
  max_cut$position,
  max_cut$position,
  ceiling(max_cut$Lhinge / 2),
  ceiling((join$unfold_N + join$unfold_C) / 2)
)

message(cmd)

system(cmd)


