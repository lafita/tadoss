# C-alpha model of a tandem domain swap from TADOSS
# Using Euclidean Distance Matrices
# More info: https://github.com/lafita/protein-edm-demo
# 
# Aleix Lafita - July 2020

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))

suppressPackageStartupMessages(library(bio3d))
suppressPackageStartupMessages(library(edmcr))

######################## Parameters #############################

seed = 0
#seed = as.integer(runif(1) * 1000000000)
tol = 0.1

######################## Input and Options #############################

# File names
pdb = "example/e1shgA1.pdb"
output = "example/e1shgA1_tswap33.pdb"
bounds = "modelling/calpha_bounds.tsv"
index = 33 # residue index to for the hinge loop
off = 1 # offset at each side of the old termini to form new loop
hinge = 1 # residues to unfold in the hinge loop

# create parser object
parser = ArgumentParser(description = 'C-alpha model of a tandem domain swap for TADOSS')

parser$add_argument("-s", "--pdb", default=pdb,
                    help="Structure file in PDB format [default \"%(default)s\"]")
parser$add_argument("-b", "--bounds", default=bounds,
                    help="Path to the file of distance bounds to use [default \"%(default)s\"]")
parser$add_argument("-i", "--index", default=index,
                    help="Residue index for the domain swap (hinge loop) [default \"%(default)s\"]")
parser$add_argument("-f", "--off", default=off,
                    help="Offset residues to unfold to join the termini [default \"%(default)s\"]")
parser$add_argument("-g", "--hinge", default=hinge,
                    help="Offset residues to unfold to form the hinge loop [default \"%(default)s\"]")
parser$add_argument("-d", "--seed", default=seed,
                    help="Random seed to use [default \"%(default)s\"]")
parser$add_argument("-t", "--toler", default=tol,
                    help="Precision tolerance for EDMC convergence [default \"%(default)s\"]")
parser$add_argument("-o", "--output", default=output,
                    help="File name for the output structural model [default \"%(default)s\"]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args = parser$parse_args()

# do some operations based on user input
pdb = args$pdb
bounds = args$bounds
index = as.integer(args$index)
off = as.integer(args$off)
hinge = as.integer(args$hinge)
output = args$output

seed = as.integer(args$seed)
set.seed(seed)
tol = as.numeric(args$toler)

######################## Parsing #############################
message("## Parsing PDB structure...")

# Use bio3d to parse the PDB structure
data.pdb.raw = read.pdb2(pdb, multi = F)

# Clean the PDB
data.pdb = clean.pdb(
  data.pdb.raw,
  consecutive = F,
  force.renumber = T,
  fix.chain = T,
  fix.aa = T,
  rm.wat = T,
  rm.lig = T,
  rm.h = T
)

# Extract atoms
pdb.atoms = data.pdb$atom %>% 
  filter(chain == "A") %>%
  filter(elety == "CA") %>%
  mutate(eleno = 1:length(eleno))

atomnum = nrow(pdb.atoms)
resnum = atomnum
revindex = resnum - index

# The atom details of the second domain - just change the IDs
pdb.atoms.2 = pdb.atoms %>%
  mutate(
    eleno = eleno + atomnum,
    resno = resno + resnum
  )

pdb.atoms.tandem = rbind(pdb.atoms, pdb.atoms.2)

# Parse the protein distance bounds constraints
calpha.bounds = read.csv(
  bounds,
  sep = "\t"
) %>% mutate(
  # Distance bounds were limited to 20A
  d.u = ifelse(d.u < 20, d.u, Inf)
)

######################## Distance matrix #############################
message("## Calculating distance matrix...")

# Compute interatomic distances
distmat.dist = dist(pdb.atoms.tandem %>% select(x, y, z), diag = T, upper = T)

distmat.df = data.frame(
  d = matrix(as.matrix(distmat.dist), ncol = 1),
  eleno.x = 1:(atomnum*2)
) %>% mutate(
  eleno.y = as.integer(0:(length(d)-1) / (atomnum*2)) + 1
)

distmat.df = merge(distmat.df, pdb.atoms.tandem, by.x = "eleno.x", by.y = "eleno")
distmat.df = merge(distmat.df, pdb.atoms.tandem, by.x = "eleno.y", by.y = "eleno")

distmat.df = distmat.df[order(distmat.df$eleno.x, distmat.df$eleno.y),]

######################## Transformation #############################
message("## Creating tandem domain distance matrix...")

# The distance matrix of the native domain
D.nat = matrix(distmat.df$d, nrow = atomnum*2)

# Tandem domain transform
distmat.df.swap = distmat.df %>%
  mutate(
    d.na = d,
    # Delete inter-domain distances from non-swap - top, bottom, right and left
    d.na = ifelse(is.element(resno.x, (index+1):(resnum+index-1)) & is.element(resno.y, 1:(index-1)), NA, d.na),
    d.na = ifelse(is.element(resno.y, (index+1):(resnum+index-1)) & is.element(resno.x, 1:(index-1)), NA, d.na),
    d.na = ifelse(is.element(resno.x, (index+1):(resnum+index-1)) & is.element(resno.y, (resnum+index):(2*resnum)), NA, d.na),
    d.na = ifelse(is.element(resno.y, (index+1):(resnum+index-1)) & is.element(resno.x, (resnum+index):(2*resnum)), NA, d.na),
    # Unfold the old terminal residues to form the loop
    d.na = ifelse(off > 0 & is.element(resno.x, (resnum-off+1):(resnum+off)) & resno.x != resno.y, NA, d.na),
    d.na = ifelse(off > 0 & is.element(resno.y, (resnum-off+1):(resnum+off)) & resno.x != resno.y, NA, d.na),
    # Unfold the hinge loop residues
    d.na = ifelse(is.element(resno.x, c((index-hinge):(index+hinge-1), (resnum+index-hinge):(resnum+index+hinge-1))) & resno.x != resno.y, NA, d.na),
    d.na = ifelse(is.element(resno.y, c((index-hinge):(index+hinge-1), (resnum+index-hinge):(resnum+index+hinge-1))) & resno.x != resno.y, NA, d.na)
  )

# Apply protein bounds to unknown distances
distmat.df = distmat.df.swap %>% mutate(
  resno.diff = resno.y - resno.x,
  # Take care of the last residue bin being representative of any difference higher
  resno.diff.bkbn = ifelse(
    abs(resno.diff) > max(calpha.bounds$resno.diff), 
    sign(resno.diff)*max(calpha.bounds$resno.diff), 
    resno.diff)
)

# Apply the upper and lower bounds
distmat.df.bounded = merge(
  distmat.df, 
  calpha.bounds, 
  by.x = c("elety.x", "elety.y", "resno.diff.bkbn"), 
  by.y = c("elety.x", "elety.y", "resno.diff"),
  all.x = T
) %>% mutate(
  d.lower = ifelse(is.na(d.na), d.l, d.na),
  d.upper = ifelse(is.na(d.na), d.u, d.na),
  # Set adjacent CA atoms to fixed distances
  d.na = ifelse(is.na(d.na) & d.avg < 4, d.avg, d.na)
)

# Make sure it is properly sorted after the merge
distmat.df.bounded = distmat.df.bounded[order(distmat.df.bounded$eleno.x, distmat.df.bounded$eleno.y),]

# Matrices of upper and lower bounds
D = matrix(distmat.df.bounded$d.na, nrow = atomnum*2)
L = matrix(distmat.df.bounded$d.lower, nrow = atomnum*2)
U = matrix(distmat.df.bounded$d.upper, nrow = atomnum*2)

# Complete upper and lower bound matrices
U = mstUB(U)
U = ceiling(U*10000000) / 10000000

message("## Completing distance matrix...")

# Print some information before running the completion
message(sprintf("##  Distance matrix: %i / %i missing entries", length(D[is.na(D)]), length(D)))

# Complete the distance matrix using EDM completion
edmc.start = Sys.time()
D.edmc = edmc(D, method = "dpf", d=3, lower=L, upper=U, toler = tol)
edmc.end = Sys.time()

D.c = D.edmc$D

message(sprintf("##  EDMC convergence value: %.1f", D.edmc$optval))
z = edmc.end - edmc.start
message(sprintf("##  EDMC running time: %.0fs", as.numeric(z, units="secs")))


######################## 3D coords #############################
message("## Multidimensional scaling to convert to 3D structure...")

# Use multidimensional scaling
fit = getConfig(D.c, d=3)
C = fit$X

message(sprintf("##  MDS accuracy: %.2f", fit$Accuracy))

######################## Save PDB #############################

# Modify the atom table and coordinates from original
data.pdb.mds = data.pdb
data.pdb.mds$atom = pdb.atoms.tandem
data.pdb.mds$xyz = matrix(t(C), nrow = 1)

# Save the model
write.pdb(
  data.pdb.mds,
  output
)

message(sprintf("## Model saved to '%s'", output))
message("## Done!")


