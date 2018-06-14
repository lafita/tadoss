'''
==============================================================

                     PUBLIC DOMAIN NOTICE
                National Institutes of Health

This software is a "United States Government Work" under the
terms of the United States Copyright Act.  It was written as
part of the authors' official duties as United States
Government employees and thus cannot be copyrighted.  This
software is freely available to the public for use. The
National Institutes of Health and the U.S. Government have not
placed any restriction on its use or reproduction.

Although all reasonable efforts have been taken to ensure
the accuracy and reliability of the software and data, the
NIH and the U.S. Government do not and cannot warrant the
performance or results that may be obtained by using this
software or data. The NIH and the U.S. Government disclaim
all warranties, express or implied, including warranties of
performance, merchantability or fitness for any particular
purpose.

Please cite the authors in any work or product based on this
material.

==============================================================
'''

import numpy as np
import argparse as ap
import warnings

from math import acos, sin, sqrt

import Bio.PDB as bp
from Bio import BiopythonWarning

################### Custom functions

def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))

def length(v):
    return sqrt(dotproduct(v, v))

def angle(v1, v2):
    return acos( dotproduct(v1, v2) / (length(v1) * length(v2)) )

################### Inputs

parser = ap.ArgumentParser()

parser.add_argument("domain", help="code of the domain to analyze", type=str)
parser.add_argument("file", help="path to the structure file", type=str)

parser.add_argument("-m", "--ext", help="residue extension needed between N anc C terminal",
                    type=int, default=-1)
parser.add_argument("-l", "--linker", help="linker length between the N and C terminal (default: 0)",
                    type=int, default=0)
parser.add_argument("-p", "--hinge", help="minimum hinge length, the loop between domains (default: 3)",
                    type=int, default=3)
parser.add_argument("-w", "--warn", help="print warning messages",
                    action='store_true')

args = parser.parse_args()


# Code of the domain
domain = args.domain

# Path to the structure file in PDB format and
path = args.file

# Note that if N- and C- termini point in opposite directions, the two termini will form the
# turn of the joint loop which does not contribute to the the effective length.
M = args.ext

# Linker length in residues (equivalent to 21A of effective length)
linker = args.linker

# Size in resdiues of the loop extension (even number) - minimum number of residues to peel to form the connecting loop.
hinge = args.hinge

# Offset in residues at the N and C terminal to start trying the circular permutant cuts
offset_cp = 10

# Silent BioPython warnings
if not args.warn:
    warnings.simplefilter('ignore', BiopythonWarning)


################### Parameters

# Temperature of reference
Tref = 350.0

# kcal/mol*K scaling factor for the go model
dS = 0.0054

# Bringing the N- and C-terminus closer together will require peeling off a small part of the native structure, starting from the termini.
# break_L is the maximum number of residues to peel from each termini.
break_L = 10

# Cutting a loop will require peeling off a small part of the native structure, starting from the break point.
# break_L is the maximum number of residues to peel from each end of the loop cut.
break_L_cut = 4

# Average length contribution of each residue.
r0 = 3.5


################### Calculation

print "1  Parsing PDB file..."

# Use BioPython to parse the PDB structure
parser = bp.PDBParser()
structure = parser.get_structure(domain, path)

# Extract all the residues of the structure in a list
res_list = bp.Selection.unfold_entities(structure, 'R')

# Total number of residues in the structure
L_p = len(res_list)

print "2  Parsing contact energies between residues..."

# Parse the GO model contact energy file
contacts = np.genfromtxt("go_{}/go_{}_gomodel_golist.dat".format(domain,domain))

# Number of free residues at N- and C-termini, defined as the number of residues which do not have native contacts.
free_N = L_p
free_C = 0

for c in contacts:
    if c[0] < free_N:
        free_N = c[0]
    if c[1] > free_C:
        free_C = c[1]

free_N_real = int(free_N) - 1
free_C_real = L_p - int(free_C)

free_N_list = []
free_C_list = []
free_NC_indx = 0
dist_NC_list = []
angle_NC_list = []
dG_joinNC_list = []
M_list = []

for free_N in range(0,max(10, free_N_real)):
    for free_C in range(0,max(10, free_C_real)):

        print "   Free resdiues at the N and C terminal: {}, {}".format(free_N, free_C)

        dist_NC = res_list[L_p - free_C - 1]['CA'] - res_list[free_N]['CA']
        dist_NC4 = res_list[L_p - free_C - 5]['CA'] - res_list[free_N+4]['CA']

        print "   Distance N to C teminal: {} A".format(round(dist_NC,1))
        print "   Distance N+4 to C-5 teminal: {} A".format(round(dist_NC4,1))

        # Calculate the angle between the N and the C termini
        c_ter_res1 = res_list[L_p - free_C - 1]['CA'].get_vector()  # get the coordinates of residue C_i
        c_ter_res2 = res_list[L_p - free_C - 5]['CA'].get_vector()  # get the coordinates of the residue C_i - 4

        n_ter_res1 = res_list[free_N]['CA'].get_vector()  # get the coordinates of residue N_j
        n_ter_res2 = res_list[free_N + 4]['CA'].get_vector()  # get the coordinates of the residue N_j + 4

        v1 = np.array(c_ter_res1) - np.array(c_ter_res2)  # get the direction of the C-terminus
        v2 = np.array(n_ter_res1) - np.array(n_ter_res2)  # get the direction of the N-terminus

        angle_NC = angle(v1, v2)  # calculate the angle between N and C - terminus

        print "   Angle between N and C temini: {} rad".format(round(angle_NC, 1))

        if M < 0:
            # Case 1, they are pointing away to each other (loop needed)
            if dist_NC4 < dist_NC:
                # instead of M = 0 or 6 as the origional paper, use the effective residue distance M' = 6*sin(ang/2)
                M = 6 * sin(angle_NC / 2)

            # Case 2, they are pointing to each other (no loop needed)
            else:
                M = 0

            print "    Parameter M set to {}".format(round(M,1))


        print "3  Calculating the dG of the tandem domain swap..."

        # Free energy change (dG_J) of joining the N- and C- termini

        energy = [0] * break_L * break_L  # List of enthalpy contribution from different i,j combination for dG_J
        entropy = [0] * break_L * break_L  # List of entropy contribution

        zero = False

        for i in range(break_L):
            for j in range(break_L):

                # Approximation of the entropy per residue on N-termini (lost in WT, gained in SWAP)
                entropy[i * break_L + j] += dS * Tref * (i)

                # Approximation of the entropy per residue on C-termini.
                entropy[i * break_L + j] += dS * Tref * (j)

                # Adjust the residue indices to the first and last residues in the domain
                n_i = i + free_N + 1
                c_i = L_p - j - free_C
                uniq_ij = []

                # The topological requirement of joining the two termini by peeling off residues
                if (i + j - M + linker) * r0 > (res_list[n_i - 1]['CA'] - res_list[c_i - 1]['CA']):
                    if i == 0 and j == 0:
                        zero = True # Indicate that the (0,0) solution is good

                    for t in range(len(contacts)):
                        for ii in range(1, n_i) + range(c_i + 1, L_p + 1):
                            if ii == contacts[t, 0] or ii == contacts[t, 1]:
                                if (contacts[t, 0], contacts[t, 1]) not in uniq_ij:

                                    # Energy of native contacts (formed in WT, lost in SWAP)
                                    energy[i * break_L + j] -= contacts[t, 3]
                                    uniq_ij.append((contacts[t, 0], contacts[t, 1]))

        # This automatically removes the (0,0) soution even if it is valid, which gives a 0 kcal/mol value
        non_zero_ind = np.nonzero(np.array(energy))
        energy_array = np.array(energy)[non_zero_ind]
        entropy_array = np.array(entropy)[non_zero_ind]

        # The minimum change of stabilities
        dG_J = entropy_array + energy_array

        dG_joinNC = max(dG_J)

        if zero:
            dG_joinNC = max(max(dG_J),0)


        print "   The dG contribution of joining termini is {} kcal/mol".format(round(dG_joinNC,1))

        if free_N == free_N_real & free_C == free_C_real:
            free_NC_indx = len(free_N_list+1)

        free_N_list.append(free_N)
        free_C_list.append(free_C)
        dist_NC_list.append(dist_NC)
        angle_NC_list.append(angle_NC)
        dG_joinNC_list.append(dG_joinNC)
        M_list.append(M)


# Free energy change (dG_C) of cutting at certain positions of the structure

dG_cuts = {}
hinge_cuts ={}

# Number of residues at each side of the cut point to extend (unfold)
ext_N = min(hinge / 2, break_L_cut)
ext_C = min((hinge + 1) / 2, break_L_cut)
break_N_cut = break_L_cut - ext_N
break_C_cut = break_L_cut - ext_C


for cp in range(free_N + offset_cp,L_p-offset_cp-free_C+1):

    energy_tot = [0] * break_N_cut * break_C_cut
    entropy_tot = [0] * break_N_cut * break_C_cut
    hinge_len = [0] * break_N_cut * break_C_cut

    for i in range(break_N_cut):
        for j in range(break_C_cut):
            uniq_ij = []
            entropy_tot[i * (break_C_cut) + j] = dS * Tref * (i + j)  # Loss entropy per residue
            hinge_len[i * (break_C_cut) + j] = i+ext_N + j+ext_C
            n_i = cp - i - ext_N
            c_i = cp + j + ext_C

            for t in range(len(contacts)):
                for ii in range(n_i, c_i):
                    if ii == contacts[t, 0] or ii == contacts[t, 1]:
                        if (contacts[t, 0], contacts[t, 1]) not in uniq_ij:

                            # Energy of native contacts (formed in WT, lost in SWAP)
                            energy_tot[i * (break_C_cut) + j] -= contacts[t, 3]
                            uniq_ij.append((contacts[t, 0], contacts[t, 1]))

    energy_tot = np.array(energy_tot)
    entropy_tot = np.array(entropy_tot)
    hinge_len = np.array(hinge_len)
    dgtot = energy_tot + entropy_tot

    dG_cuts[cp] = max(dgtot)
    hinge_cuts[cp] = hinge_len[np.argmax(dgtot)]


max_dGc = max(dG_cuts.values())
print "   The max dG cut position in the structure is {} kcal/mol".format(round(max_dGc),1)


print "4  Saving results for {}".format(domain)


print "   Saving dGjoin and domain information to '{}-dG_join.tsv'".format(domain)

with open('{}-dG_join.tsv'.format(domain), 'w') as fout:

    fout.write('domain\tfree_N\tfree_C\tdist_NC\tLlinker\tangle_NC\tM\tdGjoin\n')
    for i in range(0, len(free_N_list)):
        fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            domain, free_N_list[i], free_C_list[i], dist_NC_list[i],
            linker, angle_NC_list[i], M_list[i], dG_joinNC_list[i]))


print "   Saving dGcut profile to '{}-dG_cut.tsv'".format(domain)

with open('{}-dG_cut.tsv'.format(domain), 'w') as fout:

    fout.write('domain\tposition\tLhinge\tdGcut\n')

    for cp in dG_cuts.keys():
        fout.write('{}\t{}\t{}\t{}\n'.format(domain, cp+1, hinge_cuts[cp], dG_cuts[cp]))


print "   Saving PDB structure with dGcut in B-factors to '{}-dG_cut.pdb'".format(domain)

for p in range(0,len(res_list)):
    if p in dG_cuts.keys():
        res_list[p]['CA'].set_bfactor(dG_cuts[p])
    else:
        res_list[p]['CA'].set_bfactor(min(dG_cuts.values()))

io = bp.PDBIO()
io.set_structure(structure)
io.save('{}-dG_cut.pdb'.format(domain))


print "   Saving total alchemical ddG = dGjoin + max(dGcut) to '{}-ddG_tot.tsv'".format(domain)

with open('{}-ddG_tot.tsv'.format(domain), 'w') as fout:

    fout.write('domain\tlength\tdGjoin\tmax_dGc\tddGtot\n')
    fout.write('{}\t{}\t{}\t{}\t{}\n'.format(domain, L_p, dG_joinNC_list[free_NC_indx],
                                             max_dGc, dG_joinNC_list[free_NC_indx]+max_dGc))

