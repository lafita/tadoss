# Entry point to construct GO models

import copy,getopt,os,sys
import go_lib

# ======================================================================

Usage="""

Usage:
	go_builder.py [-h] [--help] [-t Tf] [--tf=Tf] [-k key] [--key=key]
	[-g gamma] [-n nrep] [-d nongo_dist] [--nrepeat=nrep] [-l list] [--gmx] 
	[-S scscale] [-B bbscale] [-G] [-A] [-D dslist] [--disulfides=dslist] 
	[-O] [--fixrep] [-p] [-s seq] [-H] [-c A:B:C] [-q] [--kh]
	--fixalpha --zerodihe --nativedihe --aacut=CUT file.pdb --ignore-inter --bind-n-go

	in which:
	        -q add charges
		-t Tf is the desired folding temperature (default = 350 K)
		-k key is a handy key for naming the output files (default = junk)
		-n nrep is the number of NC-linked repeats to make (default = 1)
		-g gamma is the fraction of the native contact energy assigned to
			non-native contacts (default = 0.0)
		-d nongo_dist is the sigma to use for non-go contacts when gamma > 0.0
		-l list is a list of go contacts and relative weights to be 
			used instead of the list that would usually be
			constructed.
		--gmx will tell the script to build a force-field for gromacs
			in addition to that for CHARMM [still in testing!].
		-3 threebodyexcl    will compute a 3-body map
		--water=N | -w N will put the protein in the centre of a
		box of "water" made of NxN cubes of %8.3f A each containing
		%i W4 particles, representing four waters (an equilibrated
		box I happen to have!). If you want to use a different size
		box that is relatively simple to do by hand.
		[--seq=SEQ ] will build a model based only on sequence, with
		all possible hydrogen bond partners defined, no tertiary
		contacts (?)
		-p : the file is in standard RCSB PDB format (default is 
		     to assume CHARMM pdb format
		-S scscale : scale all side-chain contacts by scscale
		-B bbscale : scale all backbone contacts by bbscale
		-G : use geometric rather than electrostatic hbond definition
		-A : use statistical angle potentials rather than Go angles
		-D dslist ... make disulfides in dslist; dslist has format
		A20-35:A50-82:B42-73 which means bonds between residues 20 & 35, in chain A etc.
		-O : make Go models compatible with Orland analytical solution
		-H: construct and save second-derivative matrix (Hessian) - two
			body terms only! 
		-c A:B:C  - use chains A,B,C only
		--fixrep: include 1-4 interactions when making nb repulsion list
		--fixalpha: correct alpha/beta balance in torsion angles
		--zerodihe: negligible dihedrals
		--nativedihe=1: native dihedral potential with native bias
		--nativedihe=2: native dihedral potential without native bias
		--kh: make parameter file with Kim-Hummer potential
		--ignore-inter : ignore intermolecular contacts
		--aacut=X : construct native contacts using a simple all-atom cutoff
			at X Angstroms

	Only the pdb file is a required argument. It must have all backbone heavy
	atoms and amide hydrogens defined, though side-chain heavy atoms are
	used too.

	All lengths are in Angstroms (even if output to gromacs will be in nm).

""" % (go_lib.std_boxlen, go_lib.std_nwat)

# ======================================================================
# PARSE INPUT
# ======================================================================

if (len(sys.argv)==1):
	print Usage
	sys.exit(0)

E0 = -2.27 # kT
LAMBDA = 0.159
kh =0
Tf = [ 350.0] 	# Kelvin rescale 350
#Tf = [ 334.9462689] 	# Kelvin rescale 350
n_repeat = 1
key = "junk"
use_list = 0
use_siglist = 0
gamma = 0.0
average = 0
gmx = 0
charmm = 1
symm = 0
water = 0
defbon = 0
scscale = 1.0
bbscale = 1.0
ngdist = 5.5	# default distance (in A) for non-native CA-CA contacts
geomh = 0
stata = 0
orland = 0
disulf_str = ""
sequence = ""
fixrep = 0
charmmpdb = 1
hessian = 0
addcharges = 0
threebodyexcl = 0
rcrowd = 6.0
rion = 3.0
ncrowd = 0
crowdbox = 80.0
crowdprot = ""
fixalpha = 0
zerodihe = 0
nativedihe = 0
ignore_inter = 0
aacut = -1.0

optlist, files = getopt.getopt(sys.argv[1:], 'hqpc:t:k:n:g:d:l:r:w:sS:B:GAD:Os:H3:', \
		['help','tf=','key=','nrepeat=','gamma=','list=',\
		'symmetrize','nongodist=','gmx','water=','scscale=',\
		'scscale=','defbon','replist=','geomh','stat-angle','disulfides=',\
		'orland','seq=','fixrep','pdb','kh', 'nocharmm',
		'ncrowd=','rcrowd=','crowdbox=','crowdprot=','fixalpha',
		'ignore-inter','bind-n-go','nativedihe=','zerodihe','aacut='])

usechains = []
for opt, val in optlist:
	if opt in [ "-t", "--tf" ]:
		Tf = map(float, val.split(':'))
	elif opt in [ "-h", "--help" ]:
		print Usage
		sys.exit(0)
	elif opt in [ "-k", "--key" ]:
		key = val
	elif opt in [ "-q", "--charge" ]:
		addcharges = 1
	elif opt in [ "-p", "--pdb" ]:
		charmmpdb = 0
	elif opt in [ "-n", "--nrepeat" ]:
		n_repeat = int(val)
	elif opt in [ "-s", "--symmetrize" ]:
		symm = 1
	elif opt in [ "-g", "--gamma" ]:
		gamma = float(val)
	elif opt in [ "-d", "--nongodist" ]:
		ngdist = float(val)
	elif opt in [ "-l", "--list" ]:	
		go_list_file = val
		use_list = 1
	elif opt in [ "-r", "--replist" ]:	
		sig_list_file = val
		use_siglist = 1
	elif opt in [ "--ignore-inter" ]:	
		ignore_inter = 1
	elif opt in [ "-w", "--water" ]:	
		water = int(val)
#	elif opt in [ "-a", "--average" ]:	
#		average = 1
	elif opt in [ "-S", "--scscale" ]:	
		scscale = float(val)
	elif opt in [ "-B", "--bbscale" ]:	
		bbscale = float(val)
	elif opt in [ "-G", "--geomh" ]:	
		geomh = 1
	elif opt in [ "-A", "--stat-angle" ]:	
		stata = 1
	elif opt in [ "--fixalpha" ]:	
		fixalpha = 1
	elif opt in [ "--zerodihe" ]:	
		zerodihe = 1
	elif opt in [ "--nativedihe" ]:	
		nativedihe = int(val)
	elif opt in [ "-D", "--disulfides" ]:	
		disulf_str = val
	elif opt in [ "--ncrowd" ]:	
		ncrowd = int(val)
	elif opt in [ "--rcrowd" ]:	
		rcrowd = float(val)
	elif opt in [ "--crowdbox" ]:	
		crowdbox = float(val)
	elif opt in [ "--crowdprot" ]:	
		crowdprot = val
	elif opt in ["--gmx" ]:	
		gmx = 1
	elif opt in ["--nocharmm" ]:	
		charmm = 0
	elif opt in ["--defbon" ]:	
		defbon = 1
	elif opt in ["-O", "--orland" ]:	
		orland = 1
	elif opt in ["--fixrep" ]:	
		fixrep = 1
	elif opt in ["--seq" ]:	
		sequence = val.upper()
	elif opt in ["--kh" ]:	
		kh = 1
	elif opt in ["-H" ]:	
		hessian = 1;
	elif opt in ["-3" ]:	
		threebodyexcl = int(val);
	elif opt in ["-c" ]:	
		usechains = val.split(':')
	elif opt in ["--aacut" ]:	
		aacut = float(val)
	elif opt in ["--bind-n-go"]:
		print "Bind-n-Go, setting:"
		print "		Kim-Hummer, charges, ignore-inter"
		kh = 1
		addcharges = 1
		ignore_inter = 1
	else:
		print "\nUNKNOWN OPTION: %s\n\n" % opt
		print Usage
		sys.exit(1)
if len(files) < 1 and len(sequence) == 0:
	print Usage
	sys.exit(1)
elif len(files) > 1:
	average = 1

# ======================================================================
# CREATE OUTPUT DIRECTORIES AND BASE FILE NAMES
# ======================================================================

if gamma < 0.00001:
	if n_repeat == 1:
		outp_base = "go_%s" % (key)
	else:
		outp_base = "go_%s_nr%i" % (key,n_repeat)
else:
	if n_repeat == 1:
		outp_base = "go_%s_ng%s" % (key,str(gamma))
	else:
		outp_base = "go_%s_ng%s_nr%i" % (key,str(gamma),n_repeat)

if symm:
	outp_base += "_sym"
if water>0:
	outp_base += "_w"

if kh == 1 and addcharges == 0:
	sys.stdout.write("Using Kim-Hummer implies adding charges\n");
	sys.stdout.write("... setting addcharges to true\n");
	addcharges = 1

outp_dir = outp_base
if not os.path.exists(outp_dir):
	os.mkdir(outp_dir)

outp_base = outp_dir + "/" + outp_base

# ======================================================================
# PRINT HANDY INFORMATION TO *_INFO.DAT FILE AND STDOUT
# ======================================================================

info = []
info.append("=========================================================\n")
command = reduce(lambda x, y: x+" "+y,sys.argv)
command = command[command.find("go_builder.py"):]
info.append("This go model generated with the command:\n\"%s\"\n" % (command))

if len(sequence) > 0:
	info.append("Building go model from sequence:\n")
	info.append("\t%s\n" % (sequence))
else:
	info.append("Building go model from PDB file %s\n" % (files[0]))
info.append("Target folding temperatures = ")
for t in Tf: 
	info.append("%8.3f " % (t))
info.append(" K\n")

info.append("Key for naming output files = %s\n" % (key))
info.append("Number of NC-linked repeats = %i\n" % (n_repeat))
info.append("Non-go weighting ('gamma') = %8.3f\n" % (gamma))
info.append("Default non-go contact distance = %8.3f\n" % (ngdist))
info.append("Scaling side-chain contacts by %8.3f\n" % (scscale))
info.append("Scaling backbone contacts by %8.3f\n" % (bbscale))
if (gmx):
	info.append("Writing GROMACS files\n")
if symm:
	info.append("Will attempt to symmetrize interactions in oligomer\n")
if average:
	info.append("Will average over multiple structures\n")
else:
	info.append("Will build model from single structure\n")
if water>0:
	info.append("Will place protein in water box of %ix%i cubes of\n"%(water,water))
	info.append("dimension %8.3f A containing %i pseudo waters\n" \
                % (go_lib.std_boxlen, go_lib.std_nwat))
	info.append("i.e. total box side = %8.3f A\n" % (float(water) * go_lib.std_boxlen))

if use_list:
	info.append("Taking go contacts from the following file = %s\n" \
			% (go_list_file))

info.append("Base for output file names = %s\n" % (outp_base))
info.append("=========================================================\n")

info_out = open(outp_base+"_info.dat","w")
for i in info:
	sys.stdout.write("%s" %(i))
	info_out.write("%s" %(i))
info_out.close()

# ======================================================================
# PARSE PDB AND WRITE OUT CA-ONLY COORDS -- IF THERE IS A PDB FILE
# ======================================================================

if len(sequence)==0:	# we are building from pdb
	if average:
		meta_list = []
		xyz_list = []
		for file in files:
			meta_crd, xyz = go_lib.make_coords(open(file).readlines(), charmmpdb)
			meta_list.append(meta_crd)
			xyz_list.append(xyz)
	else:
		meta_crd, xyz = go_lib.make_coords(open(files[0]).readlines(), charmmpdb)
else:
	# if building from sequence ... make 'linear' coords
	# this allows us to continue as we normally would for PDB data
	meta_crd, xyz = go_lib.make_sequ_coords(sequence)

chains = meta_crd.keys()
chains.sort()
nchain = len(chains)
protein_chains = copy.deepcopy(chains)

sys.stdout.write("Chains in structure file:\n")
for c in chains:
	sys.stdout.write("\t%s\n" % (c))

if len(crowdprot) != 0:	# crowder protein
	meta_crd_crowd, xyz_crowd = go_lib.make_coords(open(crowdprot).readlines(), charmmpdb)
	crowd_chain = meta_crd_crowd.keys()[0] # only allow one for now
	sys.stdout.write("Crowder chain: %s\n"%(crowd_chain))

if len(usechains)>0:
	meta_crd_nu = {}
	xyz_nu = {}
	sys.stdout.write("Using chains:\n")
	for c in usechains:
		if c in chains:
			sys.stdout.write("\t%s\n" % (c))
			meta_crd_nu[c] = meta_crd[c]
			xyz_nu[c] = xyz[c]
		else:
			sys.stdout.write("\tOops, couldn't find chain %s\n" % (c))
			sys.exit(1)
	chains = usechains
	xyz = xyz_nu
	meta_crd = meta_crd_nu
	chains.sort()
	nchain = len(chains)
	
if len(chains) >1 and n_repeat>1:
	sys.stdout.write("\n******************************************************\n")
	sys.stdout.write("WARNING!\n")
	sys.stdout.write("There is more than one chain in your system and you are\n")
	sys.stdout.write("building an N-C linked repeat protein;\n")
	sys.stdout.write("only the first chain will be used!\n")
	sys.stdout.write("******************************************************\n")

if (n_repeat==1):
	oxyz = xyz
	ochains = chains
else:
	oxyz = {chains[0]: go_lib.make_CA_nclink(xyz, n_repeat)}
	ochains = [ chains[0] ]

go_lib.write_CA_pdb(outp_base + ".pdb", oxyz, ochains)

if (gmx):
	go_lib.write_CA_gro(outp_base + ".gro", oxyz, ochains, crowdbox)

if symm:
	# all chains must AT LEAST be the same length!
	l0 = len(xyz[chains[0]])
	for chaini in chains:
		li = len(xyz[chaini])
		for chainj in chains:
			lj = len(xyz[chainj])
			if li != lj:
				sys.stderr.write("  !!! ERROR !!!  \n")
				sys.stderr.write("Cannot symmetrize contacts!\n")
				sys.stderr.write("Length of chain %s = %i\n"\
						% ( chaini, li ))
				sys.stderr.write("Length of chain %s = %i\n"\
						% ( chainj, lj ))
				sys.exit(1)

# ======================================================================
# SOLVATE PROTEIN AND WRITE COORDS IF NECC.
# ======================================================================

if water > 0:
	# (i) centre protein coordinates
	xyz,meta_crd = go_lib.center_coords(xyz, meta_crd)
	#go_lib.write_CA_pdb(outp_base + "_center.pdb", xyz, chains)
	# (ii) create BIG water box!
	wbox = go_lib.make_wbox(water)
	# (iii) attempt to solvate by deleting waters
	go_lib.solvate(xyz, wbox)
	go_lib.write_CA_pdb(outp_base + "_solv.pdb", xyz, xyz.keys())

if ncrowd > 0:
	# (i) centre protein coordinates
	xyz,meta_crd = go_lib.center_coords(xyz, meta_crd)
	#go_lib.write_CA_pdb(outp_base + "_center.pdb", xyz, chains)
	# (ii) add crowders
	if len(crowdprot) == 0:
		go_lib.crowd(xyz, ncrowd, crowdbox, rcrowd)	# use spherical crowders
	else:
		xyz_crowd,meta_crd_crowd = go_lib.center_coords(xyz_crowd, meta_crd_crowd)
		go_lib.protein_crowd(xyz, xyz_crowd, ncrowd, crowdbox, rcrowd)	# use protein crowders
	go_lib.write_CA_pdb(outp_base + "_crowd.pdb", xyz, xyz.keys())


# ======================================================================
# WRITE OUT SEQUENCES
# ======================================================================
seqs = {}
if symm:
	for c in chains:
		seqs[c] = map(lambda x: "G%s"%(x), range(1,len(xyz[c])+1))
	nres_tot = len(xyz[chains[0]])
else:
	off = 0
	nres_tot = 0
	for c in chains:
		seqs[c] = map(lambda x: "G%s"%(x), range(1+off,len(xyz[c])+1+off))
		off += len(xyz[c])
		nres_tot += len(xyz[c])


if water == 0:
	nsolv = 0
else:
	nsolv = len(xyz["SOLV"])

if charmm:
	go_lib.write_sequences(outp_base + "_seq.inp", seqs, disulf_str, key, \
                           nsolv, ncrowd, n_repeat)
	if ncrowd > 0:
		if len(crowdprot) == 0:
			go_lib.write_sequences(outp_base + "_crowd_seq.inp", \
                                   seqs, disulf_str, key, nsolv, ncrowd)
		else:
			for i in range(ncrowd):
				seqs["C%i"%(i)]= map(lambda x: "G%s"%(x), \
					range(1+off,len(xyz_crowd[crowd_chain])+1+off))
				off += len(xyz_crowd[crowd_chain])
				nres_tot += len(xyz_crowd[crowd_chain])
			go_lib.write_sequences(outp_base + "_crowd_seq.inp", \
                                   seqs, disulf_str, key, nsolv)


# ======================================================================
# WRITE OUT TOPOLOGY
# ======================================================================
if symm:
	chain = chains[0]
	nres = len(xyz[chain])
	residues = []
	for i in range(1,nres+1):
		name = meta_crd[chain][i]["name"]
		residues.append(("G%i"%(i),name))
else:
	residues = []
	idx = 0
	for c in chains:
		nres = len(xyz[c])
		for i in range(1,nres+1):
			idx += 1
			#print "chain = '%s'; i='%i'"%(c, i)
			name = meta_crd[c][i]["name"]
			residues.append(("G%i"%(idx),name))

if charmm:
	go_lib.write_topology(outp_base + "_top.inp", residues, addcharges, key)

if ncrowd > 0:
	if len(crowdprot) == 0:
		for i in range(1,ncrowd+1):
			name = "CRO"
			residues.append(("C%i"%(i),name))
	else:
		nres = len(xyz_crowd[crowd_chain])
		for x in range(1,ncrowd+1):
			for i in range(1,nres+1):
				#print x, i
				idx += 1
				name = meta_crd_crowd[crowd_chain][i]["name"]
				residues.append(("G%i"%(idx),name))

	go_lib.write_topology(outp_base + "_crowd_top.inp", residues, addcharges, key)

# ======================================================================
# CALCULATE PARAMETERS
# ======================================================================

if n_repeat > 1:
	link = 1
else:
	link = 0

ntf = 1+nchain*(nchain+1)/2
tfmap = {}
tfbb = 350 # rescal kelvin
#tfbb = 334.9462689 # rescal kelvin

#print nchain
#print Tf
#print len(Tf),ntf
if len(Tf) == ntf:
	idx = 0
	tfbb = Tf[idx]
	idx += 1 
	for i in range(nchain):
		chaini = chains[i]
		tfmap[chaini] = {}
		tfmap[chaini][chaini] = Tf[idx]
		idx+=1
		for j in range(i+1,nchain):
			chainj = chains[j]
			tfmap[chaini][chainj] = Tf[idx]
			idx+=1
else:
	if len(Tf) == 1:
		tfbb = Tf[0]
	else:
		print "wrong number of Tf's specified!"
		print len(Tf),ntf
		sys.exit(1)
	for i in range(nchain):
		chaini = chains[i]
		tfmap[chaini] = {}
		tfmap[chaini][chaini] = tfbb
		for j in range(i+1,nchain):
			chainj = chains[j]
			tfmap[chaini][chainj] = tfbb

# new ----------------------------------------------
bonds = {}
angles = {}
dihedrals = {}
nonbondedx = {}
threebody ={}
if kh==1:
	khmap = {}

#print meta_crd[chaini]
import time
time_a = time.time()
NSFAC = 1.0
disulfides = disulf_str.split(":")
#nresprot = 0
for i in range(nchain):
	chaini = chains[i]
	print chaini
	#print meta_crd[chaini][142]
	bb_bonds = go_lib.calc_bonds(meta_crd[chaini], link)
	if (defbon==1 or len(sequence)>0):
		bb_bonds = map(lambda x: (x[0], x[1], go_lib.std_bond), bb_bonds)

	bb_angles = go_lib.calc_angles(meta_crd[chaini], link, stata)
	bb_dihedrals = go_lib.setup_dihe(meta_crd[chaini], link, fixalpha, zerodihe, nativedihe)
	nonbondedx[chaini] = {}
	if kh==1:
		khmap[chaini] = {}
	#nonbonded[chaini][chaini],tmp_nsfac = go_lib.calc_ncon_intra(meta_crd[chaini],tfmap[chaini][chaini],ngdist,gamma,geomh,fixrep)
	if (aacut < 0.0): # do Karanicolas&Brooks
		nonbondedx[chaini][chaini] = go_lib.calc_ncon_intra(meta_crd[chaini], ngdist, gamma, geomh, fixrep)
	else:
		nonbondedx[chaini][chaini] = go_lib.calc_ncon_intra_simplecut(meta_crd[chaini], ngdist, \
                                                                      gamma, fixrep, aacut)
	#
	if kh==1:
		khmap[chaini][chaini] = go_lib.calc_kh_intra(meta_crd[chaini], E0, LAMBDA)

	# include disulfides if necc...
	clen = len(chaini)
	#nresprot += clen
	my_disulfides = map(lambda x: x[clen:], filter(lambda x: x.find(chaini)==0,disulfides))
	#
	if len(my_disulfides)>0:
		ds_bonds, ds_angles, ds_dihedrals = go_lib.calc_disulf_bonded(meta_crd[chaini], tfbb, my_disulfides)
		bonds[chaini] = bb_bonds+ds_bonds
		angles[chaini] = bb_angles+ds_angles
		dihedrals[chaini] = bb_dihedrals+ds_dihedrals
	else:
		bonds[chaini] = bb_bonds
		angles[chaini] = bb_angles
		dihedrals[chaini] = bb_dihedrals
	#
	# do 3 body if asked
	if threebodyexcl > 0:
		nresi = len(meta_crd[chaini])
		#nnatpair = len(nonbonded[chaini][chaini]['nc'].keys())+\
		#		len(nonbonded[chaini][chaini]['hb'].keys())
		threebody[chaini] = {}
		inte3 = tfmap[chaini][chaini] * go_lib.big_fat_fudge_factor * float(nresi)
		threebody[chaini][chaini] = go_lib.calc_threebody_intra(meta_crd[chaini], \
                                                                inte3, threebodyexcl)
		print len(threebody[chaini][chaini].keys())
	#
	# inter-subunit contacts
	for j in range(i+1,nchain):
		chainj = chains[j]
		if ignore_inter:
			# null intermolecular contacts
			nonbondedx[chaini][chainj] = { 'nc': {}, 'nnc': {}, 'hb': {}, 'siga': {}, \
			'sigb': {} } 
		else:
			if aacut < 0.0: # do K&B
				nonbondedx[chaini][chainj] = go_lib.calc_ncon_inter(meta_crd[chaini], \
                                                                    meta_crd[chainj], \
                                                                    ngdist, gamma, geomh, fixrep)
			else:
				nonbondedx[chaini][chainj] = go_lib.calc_ncon_inter_simplecut(meta_crd[chaini], \
                                                                              meta_crd[chainj], \
                                                                              ngdist, gamma, fixrep, aacut)
		if kh==1:
			khmap[chaini][chainj] = go_lib.calc_kh_inter(meta_crd[chaini], meta_crd[chainj], E0, LAMBDA)

if len(crowdprot) > 0:
	bb_bonds_crow = go_lib.calc_bonds(meta_crd_crowd[crowd_chain], link)
	if (defbon==1 or len(sequence)>0):
		bb_bonds_crow = map(lambda x: (x[0], x[1], go_lib.std_bond), bb_bonds)
	bb_angles_crow = go_lib.calc_angles(meta_crd_crowd[crowd_chain], link, stata)
	bb_dihedrals_crow = go_lib.setup_dihe(meta_crd_crowd[crowd_chain], link, fixalpha)
	tmpnb = go_lib.calc_ncon_intra(meta_crd_crowd[crowd_chain], ngdist, gamma, geomh, fixrep)
	#if kh==1:
	#	khmap[chaini][chaini] = go_lib.calc_kh_intra(meta_crd[chaini],E0,LAMBDA)
	# include disulfides if necc...
	clen = len(crowd_chain)
	for c in range(ncrowd):
		chain = "C%i"%(c)
		bonds[chain] = bb_bonds_crow
		angles[chain] = bb_angles_crow
		dihedrals[chain] = bb_dihedrals_crow
		nonbondedx[chain] = {}
		nonbondedx[chain][chain] = tmpnb
		for d in range(c+1,ncrowd):
			chaind = "C%i"%(d)
			nonbondedx[chain][chaind] = {}
			nonbondedx[chain][chaind]['nc'] = {}
			nonbondedx[chain][chaind]['hb'] = {}
			nonbondedx[chain][chaind]['nnc'] = {}
			nonbondedx[chain][chaind]['siga'] = {}
			nonbondedx[chain][chaind]['sigb'] = {}
	for chainp in chains:
		for c in range(ncrowd):
			chain = "C%i"%(c)
			nonbondedx[chainp][chain] = {}
			nonbondedx[chainp][chain]['nc'] = {}
			nonbondedx[chainp][chain]['hb'] = {}
			nonbondedx[chainp][chain]['nnc'] = {}
			nonbondedx[chainp][chain]['siga'] = {}
			nonbondedx[chainp][chain]['sigb'] = {}
	## inter-subunit contacts
	#for j in range(nchain):
	#	chainj = chains[j]
	#	tmpnb[chainj] = go_lib.calc_ncon_inter(meta_crd[chaini],meta_crd[chainj],\
	#			ngdist,gamma,geomh,fixrep)

time_b = time.time()
print time_b - time_a
#print khmap['C'].keys()
#print khmap['B'].keys()

# ======================================================================
# combine non-bonded interactions into a single list (rather than in chain-chain matrix)
# ======================================================================
if nchain > 1:
	if kh==1:
		globalkh_inter = go_lib.merge_kh_inter(protein_chains, xyz, khmap)

if kh==1:
	globalkh = go_lib.merge_kh(protein_chains, xyz, khmap)

time_c = time.time()
print time_c - time_a
#print nonbonded['PROT']
orig_angles = angles
orig_bonds = bonds
if symm:
	bonds = go_lib.average_bonded_list(bonds)
	angles = go_lib.average_bonded_list(angles)
	nonbonded = go_lib.average_nonbonded_matrix(nonbondedx)
	dihedrals = dihedrals[chains[0]]
else:
	unmerged_bonds = bonds
	unmerged_angles = angles
	unmerged_dihedrals = dihedrals
	bonds = go_lib.merge_bonded(bonds)
	angles = go_lib.merge_bonded(angles)
	dihedrals = go_lib.merge_dihe(dihedrals)
	if len(Tf) == 1:
		nonbonded,NSFAC = go_lib.merge_nonbonded_matrix(nonbondedx, xyz, Tf)
	else:
		nonbonded,NSFAC = go_lib.merge_nonbonded_matrix(nonbondedx, xyz, Tf[1:])

time_d = time.time()
print time_d - time_a
#print len(nonbonded['nnc'])
if threebodyexcl > 0:
	threebody = go_lib.merge_threebody_matrix(threebody, xyz)

if len(sequence) > 0:
	nonbonded = go_lib.calc_ncon_intra_helical(meta_crd[chains[0]], gamma)

if use_list:
	nonbonded['nc'] = go_lib.make_nc_from_list(meta_crd[chains[0]], go_list_file)
	nonbonded['hb'] = {}
if use_siglist:
	nonbonded['sig'] = go_lib.make_sig_from_list(meta_crd[chains[0]], sig_list_file)

time_e = time.time()
print time_e - time_a
#go_map,nongo_map,qlist = go_lib.make_go_nongo_map(nonbonded, \
#		nres_tot,scscale,bbscale)
go_map,qlist = go_lib.make_go_map(nonbonded['nc'], nonbonded['hb'], \
                                  nres_tot, scscale, bbscale)

if ncrowd>0 and charmm:
	if len(crowdprot) == 0:
		crowd_map,crowd_sig = go_lib.make_crowd_map(nonbonded['sig'], \
                                                    ncrowd, rcrowd)
		nres_crowd = nres_tot + ncrowd
	else:
		lcrowd = len(xyz_crowd[crowd_chain])
		ncrowd_res = lcrowd * ncrowd
		crowd_map = go_lib.make_crowd_map_protein(nonbonded['sig'], \
                                                  ncrowd_res, ncrowd, rcrowd)
		nres_crowd = nres_tot + ncrowd_res

if nchain > 1:
	go_map_intra,qlist_intra = go_lib.make_go_map(nonbonded['ncintra'], nonbonded['hbintra'], \
                                                  nres_tot, scscale, bbscale)

nongo_map = go_lib.make_nongo_map(nonbonded['nnc'], \
                                  nres_tot)

go_and_nongo = go_lib.mergemap(go_map, nongo_map)
GOMODEL_nbfix = go_lib.gotozero(go_map, nongo_map, nres_tot)
if ncrowd > 0 and charmm:
	crowd_nbfix = go_lib.gotozero(crowd_map, go_map, nres_crowd)

time_f = time.time()
print time_f - time_a

print "got here 00"
if ( kh == 1 ) and charmm:
	print "got here 01"
	GOMODEL_nbfix_kh = go_lib.gotozero(go_map, globalkh, nres_tot)
	print "got here 1"
	if ncrowd > 0:
		crowd_nbfix_kh = go_lib.gotozero(crowd_map, globalkh, nres_crowd)
	if nchain > 1:
		print "got here 2"
		GOMODEL_nbfix_kh_inter = go_lib.gotozero(go_map_intra, globalkh_inter, nres_tot)
		print "got here 3"
		GOMODEL_nbfix_kh_inter_nongo = go_lib.gotozero(go_map, globalkh_inter, nres_tot)
		if ncrowd > 0:
			GOMODEL_nbfix_kh_crowd = go_lib.gotozero(crowd_map, globalkh_inter, nres_crowd)

if water > 0:
	winte = go_lib.calc_nbond_w(residues)
else:
	winte = {}

time_g = time.time()
print time_g - time_a

# ======================================================================
# WRITE LINEAR PDB FOR EACH CHAIN.
# THESE WILL PROBABLY OVERLAP, SO WILL NEED TO BE TRANSLATED/ROTATED
# BEFORE ENERGY/FORCE EVALUATION
# ======================================================================

if len(sequence) == 0:
	for chain in chains:
		go_angles = go_lib.calc_angles(meta_crd[chain], link, 0, 0)
		bb_bonds = go_lib.calc_bonds(meta_crd[chain], link)
		tmpxyz = {chain: go_lib.make_CA_linear(bb_bonds, go_angles)}
		go_lib.write_CA_pdb(outp_base + "_linear_%s.pdb" % (chain.lower()), tmpxyz, [chain])

if n_repeat > 1:
	for chain in chains:
		#angles = go_lib.calc_angles(meta_crd[chain],link,0,0)
		#bonds = go_lib.calc_bonds(meta_crd[chain],link)
		tmpxyz = {chain: go_lib.make_CA_linear(bonds, angles, n_repeat)}
		go_lib.write_CA_pdb(outp_base + "_linear_%s.pdb" % (chain.lower()), tmpxyz, [chain])

time_h = time.time()
print time_h - time_a

# ======================================================================
# WRITE OUT PARAMETERS
# ======================================================================
go_lib.write_parms(outp_base + "_parm.inp", key, bonds, angles, dihedrals,
                   nonbonded['sig'], go_and_nongo, NSFAC, tfbb, winte, 0)

morse_sig = {}
for k in nonbonded['sig'].keys():
	morse_sig[k] = nonbonded['sig'][k]*1.025/1.1224

if charmm:
	go_lib.write_parms(outp_base + "_morse_parm.inp", key, bonds, angles, dihedrals,
                       morse_sig, go_and_nongo, NSFAC, tfbb, winte, 1)
	go_lib.write_parms(outp_base + "_gomodel_parm.inp", key, bonds, angles, dihedrals,
                       nonbonded['sig'], GOMODEL_nbfix, NSFAC, tfbb, winte, 0)
	go_lib.write_golist(outp_base + "_gomodel_golist.dat", go_map,
                        nres_tot, n_repeat)

time_i = time.time()
print time_i - time_a

if ncrowd > 0:
	if len(crowdprot) == 0:
		go_lib.write_parms(outp_base + "_crowd_parm.inp", key, bonds,
                           angles, dihedrals, crowd_sig, crowd_nbfix,
                           NSFAC, tfbb, winte, 0, ncrowd)
	else:
		go_lib.write_parms(outp_base + "_crowd_parm.inp", key, bonds,
                           angles, dihedrals, nonbonded['sig'], crowd_nbfix,
                           NSFAC, tfbb, winte, 0, 0)
		time_i = time.time(); 	print time_i - time_a
		go_lib.write_golist(outp_base + "_crowd_golist.dat", crowd_map,
                            nres_crowd, 1)
		time_i = time.time(); 	print time_i - time_a


if kh==1 and charmm:
	go_lib.write_parms(outp_base + "_kh_parm.inp", key, bonds, angles, dihedrals,
                       nonbonded['sig'], GOMODEL_nbfix_kh, NSFAC, tfbb, winte, 0)
	time_i = time.time(); 	print time_i - time_a
	go_lib.write_parms(outp_base + "_onlykh_parm.inp", key, bonds, angles, dihedrals,
                       nonbonded['sig'], globalkh, NSFAC, tfbb, winte, 0)

if nchain > 1 and charmm:
	if kh==1:
		go_lib.write_parms(outp_base + "_kh_parm_inter.inp", key, bonds, angles, dihedrals,
                           nonbonded['sig'], GOMODEL_nbfix_kh_inter, NSFAC, tfbb, winte, 0)
		time_i = time.time(); 	print time_i - time_a
		go_lib.write_parms(outp_base + "_kh_parm_inter_nongo.inp", key, bonds, angles, dihedrals,
                           nonbonded['sig'], GOMODEL_nbfix_kh_inter, NSFAC, tfbb, winte, 0)
		time_i = time.time(); 	print time_i - time_a
		if ncrowd > 0:
			go_lib.write_parms(outp_base + "_kh_parm_inter_crowd.inp", key, bonds,
                               angles, dihedrals,
                               nonbonded['sig'], GOMODEL_nbfix_kh_crowd, NSFAC, tfbb, winte, 0,
                               0)

	go_lib.write_golist(outp_base + "_gomodel_golist_intra.dat", go_map_intra,
                        nres_tot, n_repeat)
	time_i = time.time(); 	print time_i - time_a

if threebodyexcl > 0:
	go_lib.write_threebody("%s_3body_excl%i.dat" % (outp_base, threebodyexcl), threebody)

# ======================================================================
# CALCULATE AND SAVE HESSIAN IF REQUIRED
# ======================================================================

if hessian:
	sys.stdout.write("Dumping dVdr matrix\n")
	go_lib.dump_dVdr(outp_base + "_dVdr.dat", bonds, go_and_nongo, tfbb)


time_j = time.time()
print time_j - time_a


# ======================================================================
# WRITE Q-LISTS
# ======================================================================
# global qlist
go_lib.write_qlist(outp_base + "_qlist.dat", qlist, nres_tot, n_repeat)
# write GROMACS topology
# ======================================================================
if (gmx):
	# KB for all pairs
	go_lib.write_gro_top(outp_base + "_gromacs_go.top", key, residues, unmerged_bonds,
                         unmerged_angles, unmerged_dihedrals, nonbonded['sig'],
                         go_and_nongo, NSFAC, tfbb, seqs, symm, winte,"12_10_6", {}, 0, 1, n_repeat)
	if nchain > 1:
		# Go for only intramolecular contacts
		go_lib.write_gro_top(outp_base + "_gromacs_intragointerrep.top", key, residues, unmerged_bonds,
                             unmerged_angles, unmerged_dihedrals, nonbonded['sigintra'],
                             go_map_intra, NSFAC, tfbb, seqs, symm, winte,"12_10_6", {}, 0, 1, n_repeat)
	if kh==1:
		# KB for all go and KH for all other pairs
		go_lib.write_gro_top(outp_base + "_gromacs_gokh.top", key, residues, unmerged_bonds,
                             unmerged_angles, unmerged_dihedrals, nonbonded['sig'],
                             go_map, NSFAC, tfbb, seqs, symm, winte,"12_10_6", globalkh, addcharges, 1, n_repeat)
		# KH for all pairs
		go_lib.write_gro_top(outp_base + "_gromacs_onlykh.top", key, residues, unmerged_bonds,
                             unmerged_angles, unmerged_dihedrals, nonbonded['sig'],
                             {}, NSFAC, tfbb, seqs, symm, winte,"12_10_6", globalkh, addcharges, 1, n_repeat)
		if nchain > 1:
			# Go for all intramolecular go pairs and KH for all other pairs
			go_lib.write_gro_top(outp_base + "_gromacs_intragokh.top", key, residues,
                                 unmerged_bonds,
                                 unmerged_angles, unmerged_dihedrals, nonbonded['sig'],
                                 go_map_intra, NSFAC, tfbb, seqs, symm, winte,"12_10_6", globalkh,
                                 addcharges, 1, n_repeat)
			# Go for all intramolecular go pairs, repulsive intramolecular
			# nongo pairs and KH for all other pairs
			go_lib.write_gro_top(outp_base + "_gromacs_intragorepkh.top", key, residues,
                                 unmerged_bonds,
                                 unmerged_angles, unmerged_dihedrals, nonbonded['sig'],
                                 go_map_intra, NSFAC, tfbb, seqs, symm, winte,"12_10_6", globalkh,
                                 addcharges, 0, n_repeat)

		
time_k = time.time()
print "Time of q-list {}".format(time_k - time_a)

sys.exit(0)
# ======================================================================
# write 'qdetails' file
# ======================================================================
go_lib.write_qdetails(outp_base+"_qdetails.dat",meta_crd[chains[0]],nonbonded['details'],
		go_map)
# for AFM concatemers
if n_repeat>1:
	# write inter-module qlists
	go_lib.write_qlist(outp_base+"_qlist.dat",qlist,nres,n_repeat)
	go_lib.write_qlist_intra(outp_base+"_intra_qlist.dat",qlist,nres,n_repeat)
	go_lib.write_qlist_inter(outp_base+"_inter_qlist.dat",qlist,nres,n_repeat)
	for r in range(n_repeat):
		go_lib.write_module_qlist("%s_r%i_qlist.dat"%(outp_base,r),qlist,nres,r)

# for multiple chains
if nchain>1:
	for i in range(nchain):
		chaini = chains[i]
		for j in range(i,nchain):
			chainj = chains[j]
			go_lib.dump_qlist("%s_%s-%s_qlist.dat"%(outp_base,chaini,chainj), \
				qlists[chaini][chainj])

