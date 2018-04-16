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

# Library of useful functions for constructing GO models

import sys,math,os,random

# globals ==========================================================================

DATA_DIR = os.path.dirname(__file__) # absolute dir the script is in

std_bond = 3.81		# Angstroms
std_hbond = 6.3		# Angstroms
std_angle = 120.0	# degrees
std_helix_angle = 91.7	# degrees
std_urey = 6.0		# just as much a fudge as std_angle
ecrowd_std = -1.0

# little water box:
std_boxlen = 29.6	# Angstroms
std_nwat = 216
std_wbox = "%s/littlebox.pdb" % DATA_DIR
std_wcut = 7.0
nnc_dist = 5.5		# minimum E distance for non-native contacts
big_fat_fudge_factor = 0.0054 	# (kcal/(mol.K.residue)) dS per residue of unfolded state
std_sigrep = 3.0
std_bb_sigrep = 3.0
std_sc_sigrep = 3.0
hcut_theta = 135.0
hcut_dist = 3.4
qtol = 1.2	# lambda

gcrit = 0.00001

res2khsig = { 'ALA': 5.0, 'GLY': 4.5, 'THR': 5.6, 'TYR': 6.5,
	  'VAL': 5.9, 'LEU': 6.2, 'ILE': 6.2, 'TRP': 6.8,
	  'GLU': 5.9, 'ASP':5.6, 'SER': 5.2, 'ASN': 5.7,
	  'GLN': 6.0, 'PRO': 5.6, 'PHE': 6.4, 'ARG': 6.6,
	  'CYS': 5.5, 'HIS': 6.1, 'LYS': 6.4, 'MET': 6.2 }

# nb kh use +0.5 for his
res2charge = { 'ALA': 0, 'GLY': 0, 'THR': 0, 'TYR':0,
	  'VAL':0, 'LEU':0, 'ILE':0, 'TRP':0,
	  'GLU':-1, 'ASP':-1, 'SER':0, 'ASN':0,
	  'GLN':0, 'PRO':0, 'PHE':0, 'ARG':+1,
	  'CYS':0, 'HIS':+0.5, 'LYS':+1, 'MET':0, 'CRO': 0 }

aamap = { 'A':'ALA', 'G':'GLY', 'T':'THR', 'Y':'TYR',
	  'V':'VAL', 'L':'LEU', 'I':'ILE', 'W':'TRP',
	  'E':'GLU', 'D':'ASP', 'S':'SER', 'N':'ASN',
	  'Q':'GLN', 'P':'PRO', 'F':'PHE', 'R':'ARG',
	  'C':'CYS', 'H':'HIS', 'K':'LYS', 'M':'MET' }

massmap = { 'ALA': 71.0, 'GLY': 57.0,'THR': 101.0,'TYR': 163.0,
	'VAL': 99.0, 'LEU': 113.0,'ILE': 113.0,'TRP': 186.0,
	'GLU': 128.0, 'ASP': 114.0,'SER': 87.0,'ASN': 114.0,
	'GLN': 128.0, 'PRO': 97.0,'PHE': 147.0,'ARG': 157.0,
	'CYS': 103.0, 'HIS': 138.0,'LYS': 128.0,'MET': 131.0,
	'CRO': 500.0, 'DA': 328., 'DT': 319., 'DC': 304., 'DG':344.,
	'A': 328., 'T': 319., 'C': 304., 'G':344., 'U':305. }


# * if residue in scmassmap it has side-chain with mass give;
#   backbone CA has mass 56.0
# * if residue not in massmap it has a backbone CA with same mass
#   as in CA-only go model and no side-chain
# * for nucleic acids, backbone (ribose+phosphate) has mass 195, sc is rest
scmassmap = { 'TYR': 107.0, 'VAL': 43.0, 'LEU': 57.0,'ILE': 57.0,'TRP': 130.0,
	'GLU': 72.0, 'ASP': 58.0,'ASN': 58.0,'GLN': 72.0,'PHE': 91.0,'ARG': 101.0,
	'HIS': 82.0,'LYS': 72.0,'MET': 75.0,
	'A': 134, 'T': 125, 'C': 110, 'G': 150, 'U': 111,
	'DA': 134, 'DT': 125, 'DC': 110, 'DG': 150}

# water - aa interactions (eij in kcal/mol)
# charged = 1.2
# hphobic = 0.43
# hphilic = 0.81

wmap = { 'ALA': 0.43, 'GLY': 0.43,'THR': 0.81,'TYR': 0.43,
	'VAL': 0.43, 'LEU': 0.43,'ILE': 0.43,'TRP': 0.43,
	'GLU': 1.2, 'ASP': 1.2,'SER': 0.81,'ASN': 0.81,
	'GLN': 0.81, 'PRO': 0.81,'PHE': 0.43,'ARG': 1.2,
	'CYS': 0.43, 'HIS': 1.2,'LYS': 1.2,'MET': 0.43 }

# "sc_water6"
wmap_sc = { 'CA': 0.4, 'TYR': 0.3,
	'VAL': 0.3, 'LEU': 0.3,'ILE': 0.3,'TRP': 0.2,
	'GLU': 1.2, 'ASP': 1.2, 'ASN': 0.81,
	'GLN': 0.81, 'PHE': 0.3,'ARG': 1.2,
	'HIS': 1.2,'LYS': 1.2,'MET': 0.3 }

scbondmap = { 'TYR': 3.7, 'LEU': 2.6, 'ILE':2.3, 'TRP': 3.9,
		'GLU': 2.8, 'ASP':2.5, 'ASN': 2.5, 'GLN': 2.8,
		'PHE': 3.4, 'ARG': 4.1, 'HIS': 3.1, 'LYS': 3.5,
		'MET': 2.7, 'VAL': 2.0, 'DA': 4.5, 'DT': 4.0,
		'DC': 4.0, 'DG': 4.5,'A': 4.5, 'T': 4.0,
		'C': 4.0, 'G': 4.5, 'U': 4.0 }

def keysort(k1,k2):
	i1,j1 = k1
	i2,j2 = k2
	if i2!=i1:
		return i1-i2
	else:
		return j1-j2


#              angle      k1     th1   k2    th2   beta  epsilon1
statangle = { 'BACKBONE': [106.4, 91.7, 26.3, 130.0, 0.1,  4.3],
	      'TYR1':  [20.0, 95.1,  20.0, 162.5, 0.1, -0.405],
	      'TYR2': [ 20.0, 88.1,  20.0, 152.8, 0.1, 0.118],
	      'VAL1': [20.947, 117.4],
	      'VAL2': [20.947, 117.4],
	      'LEU1': [20.0, 110.5, 20.0, 144.7, 0.1, -0.667],
	      'LEU2': [20.0, 103.9, 24.0, 137.1, 0.1, 0.510],
	      'ILE1': [20.0, 120.1],
	      'ILE2': [20.0, 120.1],
	      'TRP1': [ 15.0, 101.5,  15.0, 156.5, 0.1, -0.363],
	      'TRP2': [ 15.0,  85.1,  15.0, 141.7, 0.1, -0.002],
	      'GLU1': [ 20.0, 121.5],
	      'GLU2': [ 20.0, 121.5],
	      'ASP1': [ 20.0, 106.1,  20.0, 146.5, 0.1, -0.562],
	      'ASP2': [ 20.0, 100.1,  20.0, 142.6, 0.1, 0.144],
	      'SER1': [  20.0, 113.3],
	      'SER2': [  20.0, 113.3],
	      'ASN1': [ 20.0, 104.5,  20.0, 145.7, 0.1, -0.591],
	      'ASN2': [ 20.0, 100.3,  20.0, 143.2, 0.1, 0.190],
	      'GLN1': [  20.0, 122.2],
	      'GLN2': [  20.0, 122.2],
	      'PHE1': [ 20.0, 97.1, 20.0, 158.7, 0.1, -0.564],
	      'PHE2': [ 20.0, 91.6, 20.0, 151.0, 0.1, 0.198],
	      'ARG1': [ 15.0, 121.0],
	      'ARG2': [ 15.0, 121.0],
	      'HIS1': [ 20.0, 96.4, 20.0, 153.3, 0.1, -0.657],
	      'HIS2': [ 20.0, 93.6, 20.0, 151.1, 0.1, 0.341],
	      'LYS1': [ 20.0, 124.4],
	      'LYS2': [ 20.0, 124.4],
	      'MET1': [ 20.0, 121.1],
	      'MET2': [ 20.0, 121.1] }

# functions ========================================================================

# HELPER FUNCTIONS:
#----------------------------------------------------------------------
def read_karanicolas_dihe():
	"""read in and parse Karanicolas parameters for CA-CA-CA-CA pseudo-dihedral
	angles
	REF: Karanicolas&Brooks, Prot. Sci., 11, 2351-2361"""
	inp = map(lambda x: x.split(), open("%s/karanicolas_dihe_parm.dat"%(DATA_DIR)).readlines())
	kdihe = {}
	for line in inp:
		if line[0] not in kdihe.keys():
			kdihe[line[0]] = {}
		if line[1] not in kdihe[line[0]].keys():
			kdihe[line[0]][line[1]] = []
		kdihe[line[0]][line[1]].append((float(line[2]),int(line[3]),float(line[4])))
	return kdihe

def tweak_dihe(kdihe,k1,d1):
	newdihe = {}
	d1_rad = d1*math.pi/180.

	for k in kdihe.keys():
		if k not in newdihe.keys():
			newdihe[k] = {}
			for j in kdihe[k].keys():
				newdihe[k][j] = []
				for d in kdihe[k][j]:
					if d[1]!=1:
						newdihe[k][j].append(d)
					else:
						kold = d[0]
						dold = d[2]
						dold_rad = dold*math.pi/180.
						tan_dnew = (kold*math.sin(dold_rad)+k1*math.sin(d1_rad))/ \
								(kold*math.cos(dold_rad)+k1*math.cos(d1_rad))
						dnew_rad = math.atan(tan_dnew)
						knew = (kold*math.cos(dold_rad)+k1*math.cos(d1_rad))/math.cos(dnew_rad)
						if knew < 0.0:
							dnew_rad += math.pi
							knew = -knew
						dnew = dnew_rad*180./math.pi
						#print k,j, knew, dnew
						newdihe[k][j].append((knew,1,dnew))
	return newdihe

def read_miyazawa_jernigan():
	"""read in and parse attractive Miyazawa-Jernigan parameters
	REF: Miyazawa&Jernigan, JMB, 256, 623-644 (1996)"""
	mjdat = map(lambda x: x.split(), open("%s/miyazawa_jernigan.dat"%(DATA_DIR)).readlines())
	mjmap = {}
	ave_mj = 0.0
	for d in mjdat:
		i,j,mij = d[0],d[1],float(d[2])
		if i not in mjmap.keys():
			mjmap[i] = {}
		if j not in mjmap.keys():
			mjmap[j] = {}
		mjmap[i][j] = mij
		mjmap[j][i] = mij
		ave_mj += mij
	ave_mj = abs(ave_mj/float(len(mjdat)))
	return mjmap,ave_mj

def parsecrd(x):
	return [float(x[30:38]),float(x[38:46]),float(x[46:54])]

def isNA(crds):
	if crds[1]['name'].strip() in [ 'DA','DC', 'DT','DG','A','C', 'T','G','U' ]:
		return 1
	else:
		return 0

def dist(a,b):
	dx = a[0]-b[0]
	dy = a[1]-b[1]
	dz = a[2]-b[2]
	return math.sqrt(dx**2+dy**2+dz**2)


def pbcdist(a,b,L):
	dx = a[0]-b[0]
	dy = a[1]-b[1]
	dz = a[2]-b[2]
	dx -= L*round(dx/L)
	dy -= L*round(dy/L)
	dz -= L*round(dz/L)
	return math.sqrt(dx**2+dy**2+dz**2)

def vecd(a,b):
	dx = a[0]-b[0]
	dy = a[1]-b[1]
	dz = a[2]-b[2]
	R = math.sqrt(dx**2+dy**2+dz**2)
	return dx/R,dy/R,dz/R,R

def angle(a,b,c):
	dx_ba = a[0]-b[0]
	dy_ba = a[1]-b[1]
	dz_ba = a[2]-b[2]
	r_ba = math.sqrt(dx_ba**2 + dy_ba**2 + dz_ba**2)
	dx_ba /= r_ba
	dy_ba /= r_ba
	dz_ba /= r_ba
	dx_bc = c[0]-b[0]
	dy_bc = c[1]-b[1]
	dz_bc = c[2]-b[2]
	r_bc = math.sqrt(dx_bc**2 + dy_bc**2 + dz_bc**2)
	dx_bc /= r_bc
	dy_bc /= r_bc
	dz_bc /= r_bc
	dotp = (dx_ba*dx_bc + dy_ba*dy_bc + dz_ba*dz_bc )
	return math.acos(dotp)*180.0/math.pi

def torsion(i,j,k,l):
	"return torsion angle between specified atoms"
	xij = i[0]-j[0] #self.X[i]-self.X[j]
	yij = i[1]-j[1] #self.Y[i]-self.Y[j]
	zij = i[2]-j[2] #self.Z[i]-self.Z[j]
	xjk = j[0]-k[0] #self.X[j]-self.X[k]
	yjk = j[1]-k[1] #self.Y[j]-self.Y[k]
	zjk = j[2]-k[2] #self.Z[j]-self.Z[k]
	xlk = l[0]-k[0] #self.X[l]-self.X[k]
	ylk = l[1]-k[1] #self.Y[l]-self.Y[k]
	zlk = l[2]-k[2] #self.Z[l]-self.Z[k]
	ax = yij*zjk-zij*yjk
	ay = zij*xjk-xij*zjk
	az = xij*yjk-yij*xjk
	bx = ylk*zjk-zlk*yjk
	by = zlk*xjk-xlk*zjk
	bz = xlk*yjk-ylk*xjk
	ra2 = ax**2 + ay**2 + az**2
	rb2 = bx**2 + by**2 + bz**2
	rjk2 = xjk**2 + yjk**2 + zjk**2
	rjk = math.sqrt(rjk2)
	rjkr = 1.0/rjk
	ra2r = 1.0/ra2
	rb2r = 1.0/rb2
	rabr = math.sqrt(ra2r*rb2r)
	cosphi = ( ax*bx + ay*by + az*bz ) * rabr
	#print 'cosphi = ', cosphi
	sinphi = rjk*rabr*(ax*xlk + ay*ylk + az*zlk)
	#print 'sinphi = ', sinphi
	acos = math.acos(cosphi)
	asin = math.asin(sinphi)
	#print cosphi,sinphi,acos, asin
	if acos < math.pi/2.0:
		return asin/math.pi*180.0
	elif asin >= 0:
		return acos/math.pi*180.0	#(math.pi-asin)/math.pi*180.0
	else:
		return -acos/math.pi*180.0	#(asin-math.pi)/math.pi*180.0

# PARSE COORDINATES
#----------------------------------------------------------------------
def make_coords(pdblines,charmm=1):
	"""read in pdb file and create two data structures:
	(i) for each residue, coordinates of atoms important
	for hydrogen bonding and side-chain heavy atom coordinates,
	plus residue name
	(ii) matrix of CA coordinates only"""
	crd_dat = {}
	xyz = {}
	pchain = "abcde"
	#atomrec = [ "ATOM" ]
	atomrec = [ "ATOM", "HETATM" ]
	pdblines = filter(lambda x: x[0:4] in atomrec, pdblines)
	#for line in  pdblines:
	#	print line
	for line in pdblines:
		res_name = line[17:20]

		if res_name in [ "HEM", "O2 " ]:
			continue
		res = int(line[22:26])
		if charmm:
			chain = line[72:76]	# CHARMM segid
		else:
			chain = line[21]
		if len(chain.strip()) == 0:
			chain = 'A' 	# set a default if no chain there!
		if chain != pchain:
			isNA = 0
			if res_name.strip() in [ 'DA','DT','DC','DG','A', 'T', 'C', 'G', 'U' ]:
				isNA = 1
			pchain = chain
			pres = 0
			residue = 0
		if res != pres:
			pres=res
			residue+=1
		atom = line[12:16]
		if chain not in crd_dat.keys():
			crd_dat[chain] = {}
		if chain not in xyz.keys():
			xyz[chain] = []
		if residue not in crd_dat[chain].keys():
			crd_dat[chain][residue] = {}
		crd_dat[chain][residue]['name'] = res_name
		ast = atom.strip()

		if isNA:
			if ast == "C1'":
				xyz[chain].append((residue,float(line[30:38]),float(line[38:46]),\
						float(line[46:54])))
			if ast in [ "C1'", "C2'", "C3'", "O4'", "C5'",
					"O3'","O5'","P","OP1","OP2","O2'" ]: # special atoms
				crd_dat[chain][residue][ast] = \
						(float(line[30:38]),float(line[38:46]),float(line[46:54]))
			elif ast[0] != "H":
				if (len(ast)>1):
					if ( ast[1]=="H" and ast[0].isdigit() ):
						continue
				if "SC" not in crd_dat[chain][residue].keys():
					crd_dat[chain][residue]["SC"] = []
				crd_dat[chain][residue]["SC"].append(\
						(float(line[30:38]),float(line[38:46]),\
						float(line[46:54])))
		else:
			if ast == "CA":
				xyz[chain].append((residue,float(line[30:38]),float(line[38:46]),\
						float(line[46:54])))
			if ast in [ "CA", "N", "H", "O", "C" ]: # special atoms
				crd_dat[chain][residue][ast] = \
						(float(line[30:38]),float(line[38:46]),float(line[46:54]))
			elif ast[0] != "H":
				if (len(ast)>1):
					if ( ast[1]=="H" and ast[0].isdigit() ):
						continue
				if "SC" not in crd_dat[chain][residue].keys():
					crd_dat[chain][residue]["SC"] = []
				crd_dat[chain][residue]["SC"].append(\
						(float(line[30:38]),float(line[38:46]),\
						float(line[46:54])))
	for chain in crd_dat.keys():
		for residue in crd_dat[chain].keys():
			if "SC" not in crd_dat[chain][residue].keys():
				continue
			resname = crd_dat[chain][residue]["name"].strip()
			if resname not in scbondmap.keys():
				continue
			if "CA" in crd_dat[chain][residue].keys():
				bbx,bby,bbz = crd_dat[chain][residue]["CA"]
			elif "C1'" in crd_dat[chain][residue].keys():
				bbx,bby,bbz = crd_dat[chain][residue]["C1'"]
			scx,scy,scz = scavg(crd_dat[chain][residue]["SC"])
			dx = scx-bbx
			dy = scy-bby
			dz = scz-bbz
			dr = math.sqrt(dx*dx+dy*dy+dz*dz)
			std_blen = scbondmap[resname]
			scx = bbx+dx/dr*std_blen
			scy = bby+dy/dr*std_blen
			scz = bbz+dz/dr*std_blen
			crd_dat[chain][residue]["SCCENT"] = (scx,scy,scz)
	return crd_dat,xyz

def make_sequ_coords(sequence):
	""" create linear coordinate data as for pdb file reading, but from sequence"""
	chain = 'PROT'
	crd_dat = { chain: {} }
	xyz = { chain : [] }
#
	sign = 1.0
	for i in range(len(sequence)):
		residue = i+1
		if sequence[i] not in aamap.keys():
			print "Don't know about amino acid type \"%s\"" % sequence[i]
			sys.exit(1)
		name = aamap[sequence[i]]
		crd_dat[chain][residue] = {}
		crd_dat[chain][residue]['name'] = name
		if i==0:
			x,y,z = 0.0,0.0,0.0
		elif i==1:
			x,y,z = std_bond,0.0,0.0
			#director = calc_dir(xyz[chain][0][1:],xyz[chain][1][1:])
		else:
			bondlen = std_bond
			tmpv = sc_mult(bondlen,director)
			tmpv = rotate(sign*(math.pi-std_helix_angle*math.pi/180.0),tmpv)
			x,y,z = addv(xyz[chain][-1][1:],tmpv)
			sign *= -1.0
		xyz[chain].append((residue,x,y,z))
		crd_dat[chain][residue]["CA"] = (x,y,z)
		if i>0:
			director = calc_dir(xyz[chain][-2][1:],xyz[chain][-1][1:])
	return crd_dat,xyz

def compute_indices(crd_dat):
	residues = crd_dat.keys()
	residues.sort()
	indices = {}
	na = isNA(crd_dat)
	idx = 1
	for r in residues:
		if r not in indices.keys():
			indices[r] = {}
		if na:
			indices[r]['CA'] = idx
		indices[r]['CA'] = idx
		idx+=1
		if 'SCCENT' in crd_dat[r].keys():
			indices[r]['SCCENT'] = idx
			idx+=1
	return indices

def make_wbox(N):
	"""read in water box and translate NxNxN to make big box"""
	xyz1 = []
	file = std_wbox
	ipdb = filter(lambda x: x.find("ATOM")==0, open(file).readlines())
	for line in ipdb:
		xyz1.append((float(line[30:38]),float(line[38:46]),\
					float(line[46:54])))
	OFFSET = -0.5*float(N-1)*std_boxlen
	XYZ = []
	for i in range(N):
		xoff = OFFSET+float(i)*std_boxlen
		for j in range(N):
			yoff = OFFSET+float(j)*std_boxlen
			for k in range(N):
				zoff = OFFSET+float(k)*std_boxlen
				for c in xyz1:
					XYZ.append((c[0]+xoff,c[1]+yoff,c[2]+zoff))
	return XYZ

def center_coords(xyz,meta_crd):
	comx, comy, comz = 0.0, 0.0, 0.0
	natom = 0
	chains = xyz.keys()
	# only using CA coordinates to compute center of geometry
	for chain in chains:
		for res in xyz[chain]:
			comx += res[1]
			comy += res[2]
			comz += res[3]
			natom += 1
	comx /= natom
	comy /= natom
	comz /= natom
	# center CA's in xyz
	for chain in chains:
		nres = len(xyz[chain])
		for r in range(nres):
			res, x, y, z = xyz[chain][r]
			xyz[chain][r] = (res,x-comx,y-comy,z-comz)

	#comx, comy, comz = 0.0, 0.0, 0.0
	#natom = 0
	#for chain in chains:
	#	for res in xyz[chain]:
	#		comx += res[1]
	#		comy += res[2]
	#		comz += res[3]
	#		natom += 1
	# center all atoms in meta_crd
	for chain in chains:
		nres = len(xyz[chain])
		prevca = 0
		for r in range(1,nres+1):
			#resname = pdetails[chain][r]["name"]
			cax,cay,caz = meta_crd[chain][r]["CA"]
			meta_crd[chain][r]["CA"] = \
					( cax-comx,cay-comy,caz-comz)
			if "SCCENT" in meta_crd[chain][r].keys():
				scx,scy,scz = meta_crd[chain][r]["SCCENT"]
				meta_crd[chain][r]["SCCENT"] = \
						( scx-comx,scy-comy,scz-comz)
	return xyz,meta_crd

def protein_crowd(xyz,xyz_crowd,ncrowd,crowdbox,rcrowd):
	chains = xyz.keys()
	crowd_chain = xyz_crowd.keys()[0]
	crowd_i = 1
	sys.stdout.write("Crowding protein....\n")
	sys.stdout.write("...with %i crowders\n"%(ncrowd))
	for i in range(ncrowd):
		xyz["C%i"%(i)] = []	# IMPLICIT MAX 999 CROWDERS
		good = 0
		while not good:
			xtrial = crowdbox*(random.random()-0.5)
			ytrial = crowdbox*(random.random()-0.5)
			ztrial = crowdbox*(random.random()-0.5)
			xyz_trial = [xtrial,ytrial,ztrial]
			good  = 1
			for chain in chains:
				for res in xyz[chain]:
					xyz_p = res[1:]
					for res in xyz_crowd[crowd_chain]:
						tmp = res[1:]
						xyz_c = [xtrial+tmp[0],ytrial+tmp[1],ztrial+tmp[2]]
						if pbcdist(xyz_c,xyz_p,crowdbox) \
								< 5.0 + rcrowd:
							good = 0
							print "PC OVERLAP"
							break
					if not good:
						break
				if not good:
					break
			for c in range(crowd_i-1):
				for res in xyz["C%i"%(c)]:
					xyz_c_a = res[1:]
					for res in xyz_crowd[crowd_chain]:
						tmp = res[1:]
						xyz_c_b = [xtrial+tmp[0],ytrial+tmp[1],ztrial+tmp[2]]
						if pbcdist(xyz_c_a,xyz_c_b,crowdbox) \
								< 5.0 + rcrowd:
							good = 0
							print "CC OVERLAP"
							break
					if not good:
						break
				if not good:
					break

		#nres = len(xyz_crowd[crowd_chain])
		for res in xyz_crowd[crowd_chain]:
			#tmp = xyz_crowd[crowd_chain][crowd_cres[1:]
			xyz["C%i"%(i)].append((res[0],res[1]+xtrial,res[2]+ytrial,res[3]+ztrial))
		crowd_i+=1

def crowd(xyz,ncrowd,crowdbox,rcrowd):
	chains = xyz.keys()
	xyz["CROW"] = []
	crowd_i = 1
	sys.stdout.write("Crowding protein....\n")
	sys.stdout.write("...with %i crowders\n"%(ncrowd))
	for i in range(ncrowd):
		good = 0
		while not good:
			xtrial = crowdbox*(random.random()-0.5)
			ytrial = crowdbox*(random.random()-0.5)
			ztrial = crowdbox*(random.random()-0.5)
			xyz_trial = [xtrial,ytrial,ztrial]
			good  = 1
			for chain in chains:
				for res in xyz[chain]:
					xyz_p = res[1:]
					if pbcdist(xyz_trial,xyz_p,crowdbox) \
							< 5.0 + rcrowd:
						good = 0
						print "PC OVERLAP"
						break
				if not good:
					break
			for c in range(crowd_i-1):
				xyz_c = xyz["CROW"][c][1:]
				if pbcdist(xyz_trial,xyz_c,crowdbox)\
						< 2.*rcrowd:
					good = 0
					print "CC OVERLAP"
					break

		xyz["CROW"].append([crowd_i,xtrial,ytrial,ztrial])
		crowd_i+=1

def solvate(xyz,wbox,wcut=std_wcut):
	nwat_ini = len(wbox)
	nwat = nwat_ini
	chains = xyz.keys()
	xyz["SOLV"] = []
	solv_i = 1
	sys.stdout.write("Solvating protein....\n")
	sys.stdout.write("...initial waters = %i\n"%(nwat_ini))
	for i in range(nwat_ini):
		xyz_w = wbox[i]
		delw = 0
		for chain in chains:
			for res in xyz[chain]:
				xyz_p = res[1:]
				if dist(xyz_w,xyz_p) < wcut:
					delw = 1
					break
			if delw == 1:
				break
		if delw == 0:
			xyz["SOLV"].append((solv_i,xyz_w[0],xyz_w[1],xyz_w[2]))
			solv_i+=1
	sys.stdout.write("...final waters = %i\n"%(solv_i-1))


# FUNCTIONS FOR CALCULATING PARAMETERS
# ----------------------------------------------------------------------
# intramolecular contacts - CA version

def calc_ncon_intra(coords,ngdist,gamma,geom=0,fixrep=0,hb_ecut=-0.5):
	"""The most important function: create go/non-go lists for inclusion
	in parameter file or golist file
	This version calculates contacts within a single chain"""
	ncmap = {} # native contacts
	nncmap = {} # non-native contacts
	hbmap = {} # hydrogen bond contacts
	sig_rep = {} # repulsive radii (for non-native)
	ave_mj = 0.0
	#eps_res = big_fat_fudge_factor*Tf
	nres = len(coords.keys())
	details = {}
	mjmap, ave_mj = read_miyazawa_jernigan()
	#esum = 0.0
	ca_only = 1
	na = isNA(coords)

	for i in range(1,nres+1):
		if "SC" in coords[i].keys():
			ca_only = 0
			break

	if na:
		return calc_ncon_intra_na(coords,ngdist,gamma,geom=0,fixrep=0,hb_ecut=-0.5)

	#if not ca_only:
	for i in range(1,nres-2):
		if "SC" not in coords[i].keys():
			continue
		for j in range(i+3,nres+1):
			if "SC" not in coords[j].keys():
				continue
			nc = 0
			for ai in coords[i]["SC"]:
				for aj in coords[j]["SC"]:
					if dist(ai,aj) <= 4.5:
						nc += 1
			if nc>0:
				rij = dist(coords[i]["CA"],coords[j]["CA"])
				mij = mjmap[coords[i]["name"]][coords[j]["name"]]/(ave_mj*2.0)
				ncmap[(i,j)] = [rij,mij]
				#esum += mij
	 			if (i,j) not in details.keys():
	 				details[(i,j)] = ""
				details[(i,j)] += "S"

	# count hydrogen bonds
	q1 = 0.42
	q2 = 0.20
	f = 332.0
	F = q1*q2*f
#	if not ca_only:
	for i in range(1,nres-2):
		for j in range(i+3,nres+1):
			hbon = 0
			if coords[j]["name"] != "PRO":
				if geom==1:
					xoh,yoh,zoh,roh=vecd(coords[i]["O"],coords[j]["H"])
					xnh,ynh,znh,rnh=vecd(coords[j]["N"],coords[j]["H"])
					xno,yno,zno,rno=vecd(coords[j]["N"],coords[i]["O"])
					theta=math.acos(xoh*xnh+yoh*ynh+zoh*znh)*180.0/math.pi
					if theta > hcut_theta and rno < hcut_dist:
						hbon += 1
				else:
					r_ON = dist(coords[i]["O"],coords[j]["N"])
					r_CH = dist(coords[i]["C"],coords[j]["H"])
					r_OH = dist(coords[i]["O"],coords[j]["H"])
					r_CN = dist(coords[i]["C"],coords[j]["N"])
					e_iacc = F*(1.0/r_ON+1.0/r_CH-1.0/r_OH-1.0/r_CN)
					if e_iacc < hb_ecut:
						hbon += 1
			if coords[i]["name"] != "PRO" and i!=1 and j!=nres:
				if geom==1:
					xoh,yoh,zoh,roh=vecd(coords[j]["O"],coords[i]["H"])
					xnh,ynh,znh,rnh=vecd(coords[i]["N"],coords[i]["H"])
					xno,yno,zno,rno=vecd(coords[i]["N"],coords[j]["O"])
					theta=math.acos(xoh*xnh+yoh*ynh+zoh*znh)*180.0/math.pi
					if theta > hcut_theta and rno < hcut_dist:
						hbon += 1
				else:
					r_ON = dist(coords[j]["O"],coords[i]["N"])
					r_CH = dist(coords[j]["C"],coords[i]["H"])
					r_OH = dist(coords[j]["O"],coords[i]["H"])
					r_CN = dist(coords[j]["C"],coords[i]["N"])
					e_iacc = F*(1.0/r_ON+1.0/r_CH-1.0/r_OH-1.0/r_CN)
					if e_iacc < hb_ecut:
						hbon += 1
			if hbon >0:
				if (i,j) not in details.keys():
					details[(i,j)] = ""
				if hbon == 1:
					details[(i,j)] += "H"
				elif hbon == 2:
					details[(i,j)] += "HH"
				## only actually count the hbond if no sc contact
				if (i,j) not in ncmap.keys():
					if (i,j) in hbmap.keys():
						hbmap[(i,j)][1] += -1.0
					else:
						rij = dist(coords[i]["CA"],coords[j]["CA"])
						hbmap[(i,j)] = [rij,-1.0]
					#esum += -1.0
				if hbon == 2 or (i,j) in ncmap.keys():
					if (i,j) in ncmap.keys():
						nadd = hbon
						if nadd == 1:
							dstr = "O"
						elif nadd == 2:
							dstr = "OO"
					else:
						nadd = 1
						dstr = "O"
					pqlist = []
					if i-1 >= 1:
						pqlist.append((i-1,j))
					if j-1-i >2:
						pqlist.append((i,j-1))
						pqlist.append((i+1,j))
					if j+1 <= nres:
						pqlist.append((i,j+1))
					fpq = 1.0/float(len(pqlist))*float(nadd)
					for pq in pqlist:
						p, q = pq
						if p>q:
							tmp = p
							p=q
							p=tmp
						if (p,q) not in details.keys():
							details[(p,q)] = ""
						details[(p,q)] += dstr
						if (p,q) in hbmap.keys():
							hbmap[(p,q)][1] += -fpq
						else:
							rpq = dist(coords[p]["CA"],coords[q]["CA"])
							hbmap[(p,q)] = [rpq, -fpq ]
					#esum += -fpq

	# this scaling now done in merge_nonbonded_matrix
	#nat_scale_fac = abs(float(nres)*eps_res/esum)
	#for k in ncmap.keys():
	#	ncmap[k][1] *= nat_scale_fac
	#for k in hbmap.keys():
	#	hbmap[k][1] *= nat_scale_fac

	# count non-native contacts
	if gamma > gcrit:
		for i in range(1,nres-2):
			for j in range(i+3,nres+1):
				if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys():
					mij = mjmap[coords[i]["name"]][coords[j]["name"]]/(ave_mj*2.0)*gamma
					nncmap[(i,j)] = [ngdist,mij]
	if fixrep == 0:
		excl = 4
	else:
		excl = 3	# because that's what's considered in NB interactions!
	for i in range(1,nres-excl+1):
		for j in range(i+excl,nres+1):
			if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys(): #\
				rij = dist(coords[i]["CA"],coords[j]["CA"])
				if i not in sig_rep.keys():
					sig_rep[i] = rij
				elif rij < sig_rep[i]:
					sig_rep[i] = rij
				if j not in sig_rep.keys():
					sig_rep[j] = rij
				elif rij < sig_rep[j]:
					sig_rep[j] = rij

	for i in range(1,nres+1):
		if i not in sig_rep.keys():
			sig_rep[i] = std_sigrep
		else:
			sig_rep[i] *= 1.1224

	print "num of nc(before) :",len(ncmap.keys())
	print "num of hb(before) :",len(hbmap.keys())

	nonbonded = { 'nc': ncmap, 'nnc': nncmap, 'hb': hbmap, 'sig': sig_rep, 'details': details}
	return nonbonded

def calc_ncon_intra_simplecut(coords,ngdist,gamma,fixrep,cut=4.5):
	"""The most important function: create go/non-go lists for inclusion
	in parameter file or golist file
	This version calculates contacts within a single chain"""
	ncmap = {} # native contacts
	nncmap = {} # non-native contacts
	hbmap = {} # hydrogen bond contacts
	sig_rep = {} # repulsive radii (for non-native)
	ave_mj = 0.0
	#eps_res = big_fat_fudge_factor*Tf
	nres = len(coords.keys())
	details = {}
	#esum = 0.0

	for i in range(1,nres-2):
		if "SC" in coords[i].keys():
			ilist = coords[i]["SC"]+[coords[i]['CA']]+[coords[i]['C']] \
					+[coords[i]['N']]+[coords[i]['O']]
		else:
			ilist = [coords[i]['CA']]+[coords[i]['C']] \
					+[coords[i]['N']]+[coords[i]['O']]
		for j in range(i+3,nres+1):
			if "SC" in coords[j].keys():
				jlist = coords[j]["SC"]+[coords[j]['CA']]+[coords[j]['C']] \
						+[coords[j]['N']]+[coords[j]['O']]
			else:
				jlist = [coords[j]['CA']]+[coords[j]['C']] \
						+[coords[j]['N']]+[coords[j]['O']]
			nc = 0
			for ai in ilist:
				for aj in jlist:
					if dist(ai,aj) <= cut:
						nc += 1
			if nc>0:
				rij = dist(coords[i]["CA"],coords[j]["CA"])
				mij = -1.0 #mjmap[coords[i]["name"]][coords[j]["name"]]/(ave_mj*2.0)
				ncmap[(i,j)] = [rij,mij]
	 			if (i,j) not in details.keys():
	 				details[(i,j)] = ""
	 			details[(i,j)] += "S"


	# count non-native contacts
	if gamma > gcrit:
		for i in range(1,nres-2):
			for j in range(i+3,nres+1):
				if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys():
					mij = mjmap[coords[i]["name"]][coords[j]["name"]]/(ave_mj*2.0)*gamma
					nncmap[(i,j)] = [ngdist,mij]
	if fixrep == 0:
		excl = 4
	else:
		excl = 3	# because that's what's considered in NB interactions!
	for i in range(1,nres-excl+1):
		for j in range(i+excl,nres+1):
			if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys(): #\
				rij = dist(coords[i]["CA"],coords[j]["CA"])
				if i not in sig_rep.keys():
					sig_rep[i] = rij
				elif rij < sig_rep[i]:
					sig_rep[i] = rij
				if j not in sig_rep.keys():
					sig_rep[j] = rij
				elif rij < sig_rep[j]:
					sig_rep[j] = rij

	for i in range(1,nres+1):
		if i not in sig_rep.keys():
			sig_rep[i] = std_sigrep
		else:
			sig_rep[i] *= 1.1224

	nonbonded = { 'nc': ncmap, 'nnc': nncmap, 'hb': hbmap, 'sig': sig_rep, 'details': details}
	return nonbonded

def calc_ncon_intra_na(coords,ngdist,gamma,geom=0,fixrep=0,hb_ecut=-0.5): #,sc=1):
	"""The most important function: create go/non-go lists for inclusion
	in parameter file or golist file
	This version calculates contacts within a single chain
	FOR NUCLEIC ACIDS ONLY"""
	ncmap = {} # native contacts
	nncmap = {} # non-native contacts
	hbmap = {} # hydrogen bond contacts
	sig_rep = {} # repulsive radii (for non-native)
	ave_mj = 0.0
	nres = len(coords.keys())
	mjmap = {}
	details = {}
	ca_only = 1
	na_cut = 4.5
	#dt_cut = 14.0
	dt_cut = 16.0

	# method 1: use side-chain contacts less than na_cut A
	# (if side-chains are present)
	for i in range(1,nres-2):
		if "SC" not in coords[i].keys():
			continue
			for j in range(i+3,nres+1):
				if "SC" not in coords[j].keys():
					continue
				nc = 0
				for ai in coords[i]["SC"]:
					for aj in coords[j]["SC"]:
						if dist(ai,aj) <= na_cut:
							nc += 1
				if nc>0:
					rij = dist(coords[i]["C1'"],coords[j]["C1'"])
					#mij = mjmap[coords[i]["name"]][coords[j]["name"]]/(ave_mj*2.0)
					ncmap[(i,j)] = [rij,-1.0]
					#esum += mij
		 			if (i,j) not in details.keys():
		 				details[(i,j)] = ""
		 			details[(i,j)] += "S"

	## method 2: thirumalai: bb contacts < dt_cut A
	#for i in range(1,nres-2):
	#	for j in range(i+3,nres+1):
	#		rij = dist(coords[i]["C1'"],coords[j]["C1'"])
	#		if rij < dt_cut:
	#			#mij = mjmap[coords[i]["name"]][coords[j]["name"]]/(ave_mj*2.0)
	#			ncmap[(i,j)] = [rij,-1.0]
	#			#esum += mij
	# 			if (i,j) not in details.keys():
	# 				details[(i,j)] = ""
	# 			details[(i,j)] += "S"

	# this scaling now done in merge_nonbonded_matrix
	#nat_scale_fac = abs(float(nres)*eps_res/esum)
	#for k in ncmap.keys():
	#	ncmap[k][1] *= nat_scale_fac
	#for k in hbmap.keys():
	#	hbmap[k][1] *= nat_scale_fac

	# count non-native contacts
	if gamma > gcrit:
		for i in range(1,nres-2):
			for j in range(i+3,nres+1):
				if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys():
					mij = mjmap[coords[i]["name"]][coords[j]["name"]]/(ave_mj*2.0)*gamma
					nncmap[(i,j)] = [ngdist,mij]
	if fixrep == 0:
		excl = 4
	else:
		excl = 3	# because that's what's considered in NB interactions!

	#
	# WE MAY WANT TO RECONSIDER THIS - IT COULD MAKE THE RNA VERY STIFF
	#
	for i in range(1,nres-excl+1):
		for j in range(i+excl,nres+1):
			if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys(): #\
				rij = dist(coords[i]["C1'"],coords[j]["C1'"])
				if i not in sig_rep.keys():
					sig_rep[i] = rij
				elif rij < sig_rep[i]:
					sig_rep[i] = rij
				if j not in sig_rep.keys():
					sig_rep[j] = rij
				elif rij < sig_rep[j]:
					sig_rep[j] = rij

	for i in range(1,nres+1):
		if i not in sig_rep.keys():
			sig_rep[i] = std_sigrep
		else:
			sig_rep[i] *= 1.1224

	nonbonded = { 'nc': ncmap, 'nnc': nncmap, 'hb': hbmap, 'sig': sig_rep, 'details': details}
	return nonbonded

def calc_threebody_intra(coords,inte,excl=3):
	"""Calculate a three-body interaction list from native structure -
	uniform interaction energy inte"""
	threebody_map = {}
	nres = len(coords.keys())
	#esum = 0.0
	for i in range(1,nres-2*excl+1):
		if "SC" not in coords[i].keys():
			continue
		for j in range(i+excl,nres-excl+1):
			if "SC" not in coords[j].keys():
				continue
			ncij = 0
			for ai in coords[i]["SC"]:
				for aj in coords[j]["SC"]:
					if dist(ai,aj) <= 4.5:
						ncij += 1
			# not even ij interaction ...
			if ncij == 0:
				continue
			for k in range(j+excl,nres+1):
				if "SC" not in coords[k].keys():
					continue
				ncjk = 0
				for aj in coords[j]["SC"]:
					for ak in coords[k]["SC"]:
						if dist(aj,ak) <= 4.5:
							ncjk += 1
				# no jk interaction ...
				if ncjk == 0:
					continue
				ncik = 0
				for ai in coords[i]["SC"]:
					for ak in coords[k]["SC"]:
						if dist(ai,ak) <= 4.5:
							ncik += 1
				# no ik interaction ...
				if ncik == 0:
					continue
				# NOW: we have 3-body term
				rij = dist(coords[i]["CA"],coords[j]["CA"])*1.2
				rjk = dist(coords[j]["CA"],coords[k]["CA"])*1.2
				rik = dist(coords[i]["CA"],coords[k]["CA"])*1.2
				threebody_map[(i,j,k)] = [rij,rjk,rik,inte]
	keys = threebody_map.keys()
	nkey = len(keys)
	for k in threebody_map.keys():
		threebody_map[k][3] /= float(nkey)
	return threebody_map

def calc_ncon_intra_helical(coords,gamma):
	"""Create helical hbond contacts"""
	ncmap = {}
	nncmap = {}
	hbmap = {}
	sig_rep = {}
	ave_mj = 0.0
	mjmap = {}
	details = {}
	mjmap, ave_mj = read_miyazawa_jernigan()
	nres = len(coords.keys())
	#for i in range(1,nres-3):
	#	j = i+4
	#	hbmap[(i,j)] = [std_hbond,-1.0]
	# count non-native contacts
	if gamma > gcrit:
		for i in range(1,nres-2):
			for j in range(i+3,nres+1):
				if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys():
					mij = mjmap[coords[i]["name"]][coords[j]["name"]]/(ave_mj*2.0)*gamma
					nncmap[(i,j)] = [nnc_dist,mij]
	for i in range(1,nres+1):
		sig_rep[i] = 6.0*1.1224
	nonbonded = { 'nc': ncmap, 'nnc': nncmap, 'hb': hbmap, 'sig': sig_rep, 'details': details}
	return nonbonded

def ene_sub(ncmap, energyfile):
	einp = open(energyfile)
	emap = {}
	for line in einp.readlines():
		ls = line.split()
		i, j, eij = int(ls[0]), int(ls[1]), float(ls[2])
		if j > i:
			tmp = i
			i = j
			j = tmp
		if i not in emap.keys():
			emap[i] = {}
		emap[i][j] = eij
	for k in ncmap.keys():
		i,j = k
		if j>i:
			tmp = i
			i = j
			j = tmp
		eij = emap[i][j]
		ncmap[k] = [ ncmap[k][0], eij ]
	return ncmap

def ncon_sub(ncmap, nconmap):
	avenc = 0.0
	avencon = 0.0
	nnc = 0
	for k in ncmap.keys():
		avenc += ncmap[k][1]
		avencon += nconmap[k][1]
		nnc += 1
	avenc /= float(nnc)
	avencon /= float(nnc)
	for k in ncmap.keys():
		ncmap[k] = [ ncmap[k][0], nconmap[k][1]*avenc/avencon ]
	return ncmap

# intramolecular contacts - CA+sidechain version
#def calc_sc_ncon_intra(coords,Tf,ngdist,gamma,geom=0,hb_ecut=-0.5):
def calc_sc_ncon_intra(coords,ngdist,gamma,geom=0,fixrep=0,hb_ecut=-0.5):
	"""The most important function: create go/non-go lists for inclusion
	in parameter file or golist file
	This version calculates contacts within a single chain
	Designed for go models WITH EXPLICIT SIDE-CHAINS
	"""
	if isNA(coords):
		return calc_sc_ncon_intra_na(coords,ngdist,gamma,geom,fixrep,hb_ecut)
	ncmap = {}
	nconmap = {}
	nncmap = {}
	hbmap = {}
	sig_rep = {}
	ave_mj = 0.0
	nres = len(coords.keys())
	mjmap = {}
	#eps_res = big_fat_fudge_factor*Tf
	mjmap, ave_mj = read_miyazawa_jernigan()
	#esum = 0.0
	for i in range(1,nres-2):
		if "SC" not in coords[i].keys():
			continue
		for j in range(i+3,nres+1):
			if "SC" not in coords[j].keys():
				continue
			nc = 0
			for ai in coords[i]["SC"]:
				for aj in coords[j]["SC"]:
					if dist(ai,aj) <= 4.5:
						nc += 1
			if nc>0:
				typei = coords[i]["name"]
				typej = coords[j]["name"]
				mij = mjmap[typei][typej]/(ave_mj*2.0)
				if typei in scmassmap.keys():
					sc_i = "SCCENT"
				else:
					sc_i = "CA"
				if typej in scmassmap.keys():
					sc_j = "SCCENT"
				else:
					sc_j = "CA"
				rij = dist(coords[i][sc_i],coords[j][sc_j])
				ncmap[(i,j)] = [rij,mij]
				nconmap[(i,j)] = [rij,nc]
				#esum += mij

	# count hydrogen bonds
	q1 = 0.42
	q2 = 0.20
	f = 332.0
	F = q1*q2*f
	for i in range(1,nres-2):
		for j in range(i+3,nres+1):
			hbon = 0
			if coords[j]["name"] != "PRO":
				if geom==1:
					xoh,yoh,zoh,roh=vecd(coords[i]["O"],coords[j]["H"])
					xnh,ynh,znh,rnh=vecd(coords[j]["N"],coords[j]["H"])
					xno,yno,zno,rno=vecd(coords[j]["N"],coords[i]["O"])
					theta=math.acos(xoh*xnh+yoh*ynh+zoh*znh)*180.0/math.pi
					if theta > hcut_theta and rno < hcut_dist:
						hbon += 1
				else:
					r_ON = dist(coords[i]["O"],coords[j]["N"])
					r_CH = dist(coords[i]["C"],coords[j]["H"])
					r_OH = dist(coords[i]["O"],coords[j]["H"])
					r_CN = dist(coords[i]["C"],coords[j]["N"])
					e_iacc = F*(1.0/r_ON+1.0/r_CH-1.0/r_OH-1.0/r_CN)
					if e_iacc < hb_ecut:
						hbon += 1
			if coords[i]["name"] != "PRO" and i!=1 and j!=nres:
				if geom==1:
					xoh,yoh,zoh,roh=vecd(coords[j]["O"],coords[i]["H"])
					xnh,ynh,znh,rnh=vecd(coords[i]["N"],coords[i]["H"])
					xno,yno,zno,rno=vecd(coords[i]["N"],coords[j]["O"])
					theta=math.acos(xoh*xnh+yoh*ynh+zoh*znh)*180.0/math.pi
					if theta > hcut_theta and rno < hcut_dist:
						hbon += 1
				else:
					r_ON = dist(coords[j]["O"],coords[i]["N"])
					r_CH = dist(coords[j]["C"],coords[i]["H"])
					r_OH = dist(coords[j]["O"],coords[i]["H"])
					r_CN = dist(coords[j]["C"],coords[i]["N"])
					e_iacc = F*(1.0/r_ON+1.0/r_CH-1.0/r_OH-1.0/r_CN)
					if e_iacc < hb_ecut:
						hbon += 1
			if hbon >0:
			#	## only actually count the hbond if no sc contact
				if (i,j) not in ncmap.keys():
					if (i,j) in hbmap.keys():
						hbmap[(i,j)][1] += -1.0
					else:
						rij = dist(coords[i]["CA"],coords[j]["CA"])
						hbmap[(i,j)] = [rij,-1.0]
					#esum += -1.0
				if hbon == 2 or (i,j) in ncmap.keys():
					if (i,j) in ncmap.keys():
						nadd = hbon
						if nadd == 1:
							dstr = "O"
						elif nadd == 2:
							dstr = "OO"
					else:
						nadd = 1
						dstr = "O"
					pqlist = []
					if i-1 >= 1:
						pqlist.append((i-1,j))
					if j-1-i >2:
						pqlist.append((i,j-1))
						pqlist.append((i+1,j))
					if j+1 <= nres:
						pqlist.append((i,j+1))
					fpq = 1.0/float(len(pqlist))*float(nadd)
					for pq in pqlist:
						p, q = pq
						if p>q:
							tmp = p
							p=q
							p=tmp
						if (p,q) in hbmap.keys():
							hbmap[(p,q)][1] += -fpq
						else:
							rpq = dist(coords[p]["CA"],coords[q]["CA"])
							hbmap[(p,q)] = [rpq, -fpq ]
					#esum += -fpq

	#nat_scale_fac = abs(float(nres)*eps_res/esum)
	#for k in ncmap.keys():
	#	ncmap[k][1] *= nat_scale_fac
	#for k in hbmap.keys():
	#	hbmap[k][1] *= nat_scale_fac

	# count non-native contacts
	if gamma > gcrit:
		for i in range(1,nres-2):
			for j in range(i+3,nres+1):
				if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys():
					mij = mjmap[coords[i]["name"]][coords[j]["name"]]/(ave_mj*2.0)*gamma
					nncmap[(i,j)] = [ngdist,mij*nat_scale_fac]
	if fixrep == 0:
		excl = 4
	else:
		excl = 3	# because that's what's considered in NB interactions!

	for i in range(1,nres-excl+1):
		for j in range(i+excl,nres+1):
			if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys(): #\
				rij = dist(coords[i]["CA"],coords[j]["CA"])
				if i not in sig_rep.keys():
					sig_rep[i] = rij
				elif rij < sig_rep[i]:
					sig_rep[i] = rij
				if j not in sig_rep.keys():
					sig_rep[j] = rij
				elif rij < sig_rep[j]:
					sig_rep[j] = rij

	for i in range(1,nres+1):
		sig_rep[i] *= 1.1224

	nonbonded = { 'nc': ncmap, 'nnc': nncmap, 'hb': hbmap, 'sig': sig_rep, 'ncon': nconmap }
	return nonbonded

def calc_sc_ncon_intra_na(coords,ngdist,gamma,geom=0,fixrep=0,hb_ecut=-0.5):
	"""The most important function: create go/non-go lists for inclusion
	in parameter file or golist file
	This version calculates contacts within a single chain
	Designed for go models WITH EXPLICIT SIDE-CHAINS
	"""
	ncmap = {}
	nconmap = {}
	nncmap = {}
	hbmap = {}
	sig_rep = {}
	ave_mj = 0.0
	nres = len(coords.keys())
	mjmap = {}
	#eps_res = big_fat_fudge_factor*Tf
	mjmap, ave_mj = read_miyazawa_jernigan()
	#esum = 0.0
	#for i in range(1,nres-2):
	#	for j in range(i+3,nres+1):
	for i in range(1,nres):
		for j in range(i+1,nres+1):
			nc = 0
			for ai in coords[i]["SC"]:
				for aj in coords[j]["SC"]:
					if dist(ai,aj) <= 5.0:
						nc += 1
			if nc>0:
				typei = coords[i]["name"]
				typej = coords[j]["name"]
				#mij = mjmap[typei][typej]/(ave_mj*2.0)
				#if typei in scmassmap.keys():
				sc_i = "SCCENT"
				#else:
				#	sc_i = "C1'"
				#if typej in scmassmap.keys():
				sc_j = "SCCENT"
				#else:
				#	sc_j = "C1'"
				rij = dist(coords[i][sc_i],coords[j][sc_j])
				ncmap[(i,j)] = [rij,-1]
				nconmap[(i,j)] = [rij,nc]
				#esum += mij

	# count non-native contacts
	if gamma > gcrit:
		for i in range(1,nres-2):
			for j in range(i+3,nres+1):
				if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys():
					mij = mjmap[coords[i]["name"]][coords[j]["name"]]/(ave_mj*2.0)*gamma
					nncmap[(i,j)] = [ngdist,mij*nat_scale_fac]
	if fixrep == 0:
		excl = 4
	else:
		excl = 3	# because that's what's considered in NB interactions!

	for i in range(1,nres-excl+1):
		for j in range(i+excl,nres+1):
			if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys(): #\
				rij = dist(coords[i]["C1'"],coords[j]["C1'"])
				if i not in sig_rep.keys():
					sig_rep[i] = rij
				elif rij < sig_rep[i]:
					sig_rep[i] = rij
				if j not in sig_rep.keys():
					sig_rep[j] = rij
				elif rij < sig_rep[j]:
					sig_rep[j] = rij

	for i in range(1,nres+1):
		sig_rep[i] *= 1.1224

	nonbonded = { 'nc': ncmap, 'nnc': nncmap, 'hb': hbmap, 'sig': sig_rep, 'ncon': nconmap }
	return nonbonded

def scale_tertiary(emap,w,dij):
	"""scale all contacts with |i-j|>dij by w"""
	numap = {}
	for k in emap.keys():
		i = k[0]; j = k[1]
		if abs(i-j) > dij:
			numap[k] = [emap[k][0],emap[k][1]*w]
		else:
			numap[k] = emap[k]
	return numap

# calculate intermolecular contacts - CA version
def calc_ncon_inter(crda,crdb,ngdist,gamma,geom=0,fixrep=0,hb_ecut=-0.5):
	"""The most important function: create go/non-go lists for inclusion
	in parameter file or golist file
	This version calculates contacts between two chains"""
	ncmap = {}
	nncmap = {}
	hbmap = {}
	sig_rep_a = {}
	sig_rep_b = {}
	ave_mj = 0.0
	#eps_res = big_fat_fudge_factor*Tf

	na_a = isNA(crda)
	na_b = isNA(crdb)

	if na_a and na_b:
		return calc_ncon_inter_na(crda,crdb,ngdist,gamma,geom=0,fixrep=0,hb_ecut=-0.5)

	if ( na_a and not na_b ) or ( na_b and not na_a ):
		return calc_ncon_inter_protna(crda,crdb,ngdist,gamma,geom=0,fixrep=0,hb_ecut=-0.5)

	nresa = len(crda.keys())
	nresb = len(crdb.keys())
	mjmap = {}
	mjmap, ave_mj = read_miyazawa_jernigan()
	#esum = 0.0
	iresa, iresb = [], []
	for i in range(1,nresa+1):
		if "SC" not in crda[i].keys():
			continue
		for j in range(1,nresb+1):
			if "SC" not in crdb[j].keys():
				continue
			nc = 0
			for ai in crda[i]["SC"]:
				for aj in crdb[j]["SC"]:
					if dist(ai,aj) <= 4.5:
						nc += 1
			if nc>0:
				rij = dist(crda[i]["CA"],crdb[j]["CA"])
				mij = mjmap[crda[i]["name"]][crdb[j]["name"]]/(ave_mj*2.0)
				ncmap[(i,j)] = [rij,mij]
				#esum += mij
				if i not in iresa:
					iresa.append(i)
				if j not in iresb:
					iresb.append(j)

	# count hydrogen bonds
	q1 = 0.42
	q2 = 0.20
	f = 332.0
	F = q1*q2*f
	for i in range(1,nresa+1):
		for j in range(1,nresb+1):
			hbon = 0
			if crdb[j]["name"] != "PRO" and j!=1 and i!=nresa:
				if geom==1:
					xoh,yoh,zoh,roh=vecd(crda[i]["O"],crdb[j]["H"])
					xnh,ynh,znh,rnh=vecd(crdb[j]["N"],crdb[j]["H"])
					xno,yno,zno,rno=vecd(crdb[j]["N"],crda[i]["O"])
					theta=math.acos(xoh*xnh+yoh*ynh+zoh*znh)*180.0/math.pi
					if theta > hcut_theta and rno < hcut_dist:
						hbon += 1
				else:
					r_ON = dist(crda[i]["O"],crdb[j]["N"])
					r_CH = dist(crda[i]["C"],crdb[j]["H"])
					r_OH = dist(crda[i]["O"],crdb[j]["H"])
					r_CN = dist(crda[i]["C"],crdb[j]["N"])
					e_iacc = F*(1.0/r_ON+1.0/r_CH-1.0/r_OH-1.0/r_CN)
					if e_iacc < hb_ecut:
						hbon += 1
			if crda[i]["name"] != "PRO" and i!=1 and j!=nresb:
				if geom==1:
					xoh,yoh,zoh,roh=vecd(crda[j]["O"],crdb[i]["H"])
					xnh,ynh,znh,rnh=vecd(crdb[i]["N"],crdb[i]["H"])
					xno,yno,zno,rno=vecd(crdb[i]["N"],crda[j]["O"])
					theta=math.acos(xoh*xnh+yoh*ynh+zoh*znh)*180.0/math.pi
					if theta > hcut_theta and rno < hcut_dist:
						hbon += 1
				else:
					r_ON = dist(crda[i]["N"],crdb[j]["O"])
					r_CH = dist(crda[i]["H"],crdb[j]["C"])
					r_OH = dist(crda[i]["H"],crdb[j]["O"])
					r_CN = dist(crda[i]["N"],crdb[j]["C"])
					e_iacc = F*(1.0/r_ON+1.0/r_CH-1.0/r_OH-1.0/r_CN)
					if e_iacc < hb_ecut:
						hbon += 1
			if hbon >0:
				## only actually count the hbond if no sc contact
				if (i,j) not in ncmap.keys():
					if (i,j) in hbmap.keys():
						hbmap[(i,j)][1] += -1.0
					else:
						rij = dist(crda[i]["CA"],crdb[j]["CA"])
						hbmap[(i,j)] = [rij,-1.0]
					#esum += -1.0
					if i not in iresa:
						iresa.append(i)
					if j not in iresb:
						iresb.append(j)

				if hbon == 2 or (i,j) in ncmap.keys():
					if (i,j) in ncmap.keys():
						nadd = hbon
						if nadd == 1:
							dstr = "O"
						elif nadd == 2:
							dstr = "OO"
					else:
						nadd = 1
						dstr = "O"
					pqlist = []
					if i-1 >= 1:
						pqlist.append((i-1,j))
					if j-1>=1:
						pqlist.append((i,j-1))
					if i+1<nresa:
						pqlist.append((i+1,j))
					if j+1 <= nresb:
						pqlist.append((i,j+1))
					fpq = 1.0/float(len(pqlist))*float(nadd)
					for pq in pqlist:
						p, q = pq
						if p>q:
							tmp = p
							p=q
							p=tmp
						if (p,q) in hbmap.keys():
							hbmap[(p,q)][1] += -fpq
						else:
							rpq = dist(crda[p]["CA"],crdb[q]["CA"])
							hbmap[(p,q)] = [rpq, -fpq ]
						#esum += -fpq
						if p not in iresa:
							iresa.append(p)
						if q not in iresb:
							iresb.append(q)

	nres = len(iresa) + len(iresb)
	print "number interface residues = ", nres
	#nat_scale_fac = abs(float(nres)*eps_res/esum)
	#for k in ncmap.keys():
	#	ncmap[k][1] *= nat_scale_fac
	#for k in hbmap.keys():
	#	hbmap[k][1] *= nat_scale_fac

	# count non-native contacts
	if gamma > gcrit:
		for i in range(1,nresa+1):
			for j in range(1,nresb+1):
				if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys():
					mij = mjmap[crda[i]["name"]][crdb[j]["name"]]/(ave_mj*2.0)*gamma
					nncmap[(i,j)] = [ngdist,mij]
	for i in range(1,nresa+1):
		for j in range(1,nresb+1):
			if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys(): #\
				rij = dist(crda[i]["CA"],crdb[j]["CA"])
				if i not in sig_rep_a.keys():
					sig_rep_a[i] = rij
				elif rij < sig_rep_a[i]:
					sig_rep_a[i] = rij
				if j not in sig_rep_b.keys():
					sig_rep_b[j] = rij
				elif rij < sig_rep_b[j]:
					sig_rep_b[j] = rij

	for i in range(1,nresa+1):
		if i in sig_rep_a.keys():
			sig_rep_a[i] *= 1.1224
		else:
			sig_rep_a[i] = std_sigrep
	for i in range(1,nresb+1):
		if i in sig_rep_b.keys():
			sig_rep_b[i] *= 1.1224
		else:
			sig_rep_b[i] = std_sigrep

	nonbonded = { 'nc': ncmap, 'nnc': nncmap, 'hb': hbmap, 'siga': sig_rep_a, \
			'sigb': sig_rep_b }
	return nonbonded




def calc_ncon_inter_na(crda,crdb,ngdist,gamma,geom=0,fixrep=0,hb_ecut=-0.5):
	"""The most important function: create go/non-go lists for inclusion
	in parameter file or golist file
	This version calculates contacts between two chains"""
	ncmap = {}
	nncmap = {}
	hbmap = {}
	sig_rep_a = {}
	sig_rep_b = {}
	ave_mj = 0.0
	#eps_res = big_fat_fudge_factor*Tf
	ca_only = 1
	na_cut = 4.5
	#dt_cut = 14.0
	dt_cut = 16.0

	nresa = len(crda.keys())
	nresb = len(crdb.keys())
	mjmap = {}
	mjmap, ave_mj = read_miyazawa_jernigan()
	#esum = 0.0
	iresa, iresb = [], []
	for i in range(1,nresa+1):
		if "SC" not in crda[i].keys():
			continue
		for j in range(1,nresb+1):
			if "SC" not in crdb[j].keys():
				continue
			nc = 0
			for ai in crda[i]["SC"]:
				for aj in crdb[j]["SC"]:
					if dist(ai,aj) <= na_cut:
						nc += 1
			if nc>0:
				print crda[i]
				print crdb[i]
				rij = dist(crda[i]["C1'"],crdb[j]["C1'"])
				#mij = mjmap[crda[i]["name"]][crdb[j]["name"]]/(ave_mj*2.0)
				ncmap[(i,j)] = [rij,-1.0]
				#esum += mij
				if i not in iresa:
					iresa.append(i)
				if j not in iresb:
					iresb.append(j)


	nres = len(iresa) + len(iresb)
	print "number interface residues = ", nres
	#nat_scale_fac = abs(float(nres)*eps_res/esum)
	#for k in ncmap.keys():
	#	ncmap[k][1] *= nat_scale_fac
	#for k in hbmap.keys():
	#	hbmap[k][1] *= nat_scale_fac

	# count non-native contacts (not working!!)
	if gamma > gcrit:
		for i in range(1,nresa+1):
			for j in range(1,nresb+1):
				if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys():
					mij = mjmap[crda[i]["name"]][crdb[j]["name"]]/(ave_mj*2.0)*gamma
					nncmap[(i,j)] = [ngdist,mij]
	for i in range(1,nresa+1):
		for j in range(1,nresb+1):
			if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys(): #\
				rij = dist(crda[i]["C1'"],crdb[j]["C1'"])
				if i not in sig_rep_a.keys():
					sig_rep_a[i] = rij
				elif rij < sig_rep_a[i]:
					sig_rep_a[i] = rij
				if j not in sig_rep_b.keys():
					sig_rep_b[j] = rij
				elif rij < sig_rep_b[j]:
					sig_rep_b[j] = rij

	for i in range(1,nresa+1):
		if i in sig_rep_a.keys():
			sig_rep_a[i] *= 1.1224
		else:
			sig_rep_a[i] = std_sigrep
	for i in range(1,nresb+1):
		if i in sig_rep_b.keys():
			sig_rep_b[i] *= 1.1224
		else:
			sig_rep_b[i] = std_sigrep

	nonbonded = { 'nc': ncmap, 'nnc': nncmap, 'hb': hbmap, 'siga': sig_rep_a, \
			'sigb': sig_rep_b }
	return nonbonded

def calc_ncon_inter_protna(crda,crdb,ngdist,gamma,geom=0,fixrep=0,hb_ecut=-0.5):
	"""The most important function: create go/non-go lists for inclusion
	in parameter file or golist file
	This version calculates contacts between two chains"""
	ncmap = {}
	nncmap = {}
	hbmap = {}
	sig_rep_a = {}
	sig_rep_b = {}
	ave_mj = 0.0
	#eps_res = big_fat_fudge_factor*Tf
	ca_only = 1
	na_prot_cut = 6.0

	if isNA(crda):
		cakey_a = "C1'"
	else:
		cakey_a = "CA"

	if isNA(crdb):
		cakey_b = "C1'"
	else:
		cakey_b = "CA"

	nresa = len(crda.keys())
	nresb = len(crdb.keys())
	mjmap = {}
	mjmap, ave_mj = read_miyazawa_jernigan()
	#esum = 0.0
	iresa, iresb = [], []
	for i in range(1,nresa+1):
		if "SC" not in crda[i].keys():
			continue
		for j in range(1,nresb+1):
			if "SC" not in crdb[j].keys():
				continue
			nc = 0
			for ai in crda[i]["SC"]:
				for aj in crdb[j]["SC"]:
					if dist(ai,aj) <= na_prot_cut:
						nc += 1
			if nc>0:
				rij = dist(crda[i][cakey_a],crdb[j][cakey_b])
				#mij = mjmap[crda[i]["name"]][crdb[j]["name"]]/(ave_mj*2.0)
				ncmap[(i,j)] = [rij,-1.0]
				#esum += mij
				if i not in iresa:
					iresa.append(i)
				if j not in iresb:
					iresb.append(j)


	nres = len(iresa) + len(iresb)
	print "number interface residues = ", nres
	#nat_scale_fac = abs(float(nres)*eps_res/esum)
	#for k in ncmap.keys():
	#	ncmap[k][1] *= nat_scale_fac
	#for k in hbmap.keys():
	#	hbmap[k][1] *= nat_scale_fac

	# count non-native contacts (not working!!)
	if gamma > gcrit:
		for i in range(1,nresa+1):
			for j in range(1,nresb+1):
				if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys():
					mij = mjmap[crda[i]["name"]][crdb[j]["name"]]/(ave_mj*2.0)*gamma
					nncmap[(i,j)] = [ngdist,mij]
	for i in range(1,nresa+1):
		for j in range(1,nresb+1):
			if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys(): #\
				rij = dist(crda[i][cakey_a],crdb[j][cakey_b])
				if i not in sig_rep_a.keys():
					sig_rep_a[i] = rij
				elif rij < sig_rep_a[i]:
					sig_rep_a[i] = rij
				if j not in sig_rep_b.keys():
					sig_rep_b[j] = rij
				elif rij < sig_rep_b[j]:
					sig_rep_b[j] = rij

	for i in range(1,nresa+1):
		if i in sig_rep_a.keys():
			sig_rep_a[i] *= 1.1224
		else:
			sig_rep_a[i] = std_sigrep
	for i in range(1,nresb+1):
		if i in sig_rep_b.keys():
			sig_rep_b[i] *= 1.1224
		else:
			sig_rep_b[i] = std_sigrep

	nonbonded = { 'nc': ncmap, 'nnc': nncmap, 'hb': hbmap, 'siga': sig_rep_a, \
			'sigb': sig_rep_b }
	return nonbonded



def calc_ncon_inter_simplecut(crda,crdb,ngdist,gamma,fixrep,cut=4.5):
	"""The most important function: create go/non-go lists for inclusion
	in parameter file or golist file
	This version calculates contacts between two chains"""
	ncmap = {}
	nncmap = {}
	hbmap = {}
	sig_rep_a = {}
	sig_rep_b = {}
	ave_mj = 0.0
	nresa = len(crda.keys())
	nresb = len(crdb.keys())
	iresa, iresb = [], []
	for i in range(1,nresa+1):
		ilist = crda[i]["SC"]+[crda[i]['CA']]+[crda[i]['C']] \
			+[crda[i]['N']]+[crda[i]['O']]
		for j in range(1,nresb+1):
			nc = 0
			jlist = crdb[j]["SC"]+[crdb[j]['CA']]+[crdb[j]['C']] \
				+[crdb[j]['N']]+[crdb[j]['O']]
			for ai in ilist:
				for aj in jlist:
					if dist(ai,aj) <= 4.5:
						nc += 1
			if nc>0:
				rij = dist(crda[i]["CA"],crdb[j]["CA"])
				mij = -1.0 #mij = mjmap[crda[i]["name"]][crdb[j]["name"]]/(ave_mj*2.0)
				ncmap[(i,j)] = [rij,mij]
				#esum += mij
				if i not in iresa:
					iresa.append(i)
				if j not in iresb:
					iresb.append(j)


	nres = len(iresa) + len(iresb)
	print "number interface residues = ", nres
	if gamma > gcrit:
		for i in range(1,nresa+1):
			for j in range(1,nresb+1):
				if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys():
					mij = mjmap[crda[i]["name"]][crdb[j]["name"]]/(ave_mj*2.0)*gamma
					nncmap[(i,j)] = [ngdist,mij]
	for i in range(1,nresa+1):
		for j in range(1,nresb+1):
			if (i,j) not in ncmap.keys() and (i,j) not in hbmap.keys(): #\
				rij = dist(crda[i]["CA"],crdb[j]["CA"])
				if i not in sig_rep_a.keys():
					sig_rep_a[i] = rij
				elif rij < sig_rep_a[i]:
					sig_rep_a[i] = rij
				if j not in sig_rep_b.keys():
					sig_rep_b[j] = rij
				elif rij < sig_rep_b[j]:
					sig_rep_b[j] = rij

	for i in range(1,nresa+1):
		if i in sig_rep_a.keys():
			sig_rep_a[i] *= 1.1224
		else:
			sig_rep_a[i] = std_sigrep
	for i in range(1,nresb+1):
		if i in sig_rep_b.keys():
			sig_rep_b[i] *= 1.1224
		else:
			sig_rep_b[i] = std_sigrep

	nonbonded = { 'nc': ncmap, 'nnc': nncmap, 'hb': hbmap, 'siga': sig_rep_a, \
			'sigb': sig_rep_b }
	return nonbonded

def calc_kh_intra(crda,E0,LAMBDA):
	"""Calculate kh matrix of interactions within a chain"""
	khmap = {}
	nresa = len(crda.keys())
	mjmap = {}
	mjmap, ave_mj = read_miyazawa_jernigan()
	for i in range(1,nresa+1):
		for j in range(i,nresa+1):
			itype = crda[i]["name"]
			jtype = crda[j]["name"]
			sigi = res2khsig[itype]
			sigj = res2khsig[jtype]
			sigij = (sigi+sigj)/2.
			epsij = LAMBDA*(mjmap[itype][jtype]-E0)*0.6 # to get kcal/mol from kT
			if abs(epsij) < 1e-20:	# i.e. zero! - would be no pair interaction calculated
				epsij = 1.e-4	# introduce small hard-core repulsion
			khmap[(i,j)] = [sigij,epsij]

	return khmap

def calc_kh_inter(crda,crdb,E0,LAMBDA):
	"""Calculate kh matrix of interactions between two chains"""
	khmap = {}
	nresa = len(crda.keys())
	nresb = len(crdb.keys())
	mjmap = {}
	mjmap, ave_mj = read_miyazawa_jernigan()
	for i in range(1,nresa+1):
		for j in range(1,nresb+1):
			itype = crda[i]["name"]
			jtype = crdb[j]["name"]
			sigi = res2khsig[itype]
			sigj = res2khsig[jtype]
			sigij = (sigi+sigj)/2.
			epsij = LAMBDA*(mjmap[itype][jtype]-E0)*0.6 # to get kcal/mol from kT
			if abs(epsij) < 1e-20:	# i.e. zero! - would be no pair interaction calculated
				epsij = 1.e-4	# introduce small hard-core repulsion
			khmap[(i,j)] = [sigij,epsij]

	return khmap

def make_nc_from_list(pdb, list_file):
	"""make a go contact list read in from a file with lines:
	i   j   e_ij
	Note that no scaling of energies is done as there is a second
	rescaling later!
	Also attractive energies in the list should be negative"""
	listdat = map(lambda x: x.split(), open(list_file).readlines())
	ncmap = {}
	for d in listdat:
		i, j, e_ij = int(d[0]), int(d[1]), float(d[2])
		if len(d) == 4:
			rij = float(d[3])
		else:
			rij = dist(pdb[i]["CA"],pdb[j]["CA"])
		ncmap[(i,j)] = [rij,e_ij]
	return ncmap

def make_sig_from_list(pdb, list_file):
	"""make a sig_rep list read in from a file with lines:
	i   ri"""
	listlist = map(lambda x: x.split(), open(list_file).readlines())
	listdat = map(lambda x: (int(x[0]),float(x[1])), listlist)
	sig_rep = {}
	for d in listdat:
		i, ri = d[0], d[1]
		sig_rep[i] = ri
	return sig_rep

def calc_nbond_w(residues):
	"""calculate residue-water nbfix"""
	winte = {}
	for r in residues:
		res, type = r
		winte[(res,"W4")] = ( 5.28, -wmap[type] )

	return winte

def calc_nbond_w_sc(residues):
	"""calculate residue-water nbfix"""
	winte = {}
	for r in residues:
		res, type = r
		if type in wmap_sc.keys(): # there is an explicit side-chain
			winte[(res,"W4")] = ( 5.28, -wmap_sc['CA'] )
			winte[("%sS"%res,"W4")] = ( 5.28, -wmap_sc[type] )
		else:
			winte[(res,"W4")] = ( 5.28, -wmap[type] )

	return winte

def symm_emap(emap,nres):
	symmap = {}
	for i in range(1,nres):
		for j in range(i+1,nres+1):
			k = (i,j)
			kk = (j,i)
			if k in emap.keys():
				if kk in emap.keys():
					symmap[k] = [(emap[k][0]+emap[kk][0])/2.0,\
							(emap[k][1]+emap[kk][1])/2.0]
				else:
					symmap[k] = emap[k]
			elif kk in emap.keys():
				symmap[k] = emap[kk]
	return symmap

def symmetrize_nb(nb,nres):
	"""Make symmetric non-bonded list for symmetric dimers"""
	ncmap = nb['nc']
	nncmap = nb['nnc']
	hbmap = nb['hb']
	sig_rep_a = nb['siga']
	sig_rep_b = nb['sigb']
	#details = nb['details']
	nres = len(sig_rep_a)
	ncmap_s = symm_emap(ncmap,nres)
	nncmap_s = symm_emap(nncmap,nres)
	hbmap_s = symm_emap(hbmap,nres)
	sig = {}
	for k in sig_rep_a.keys():
		sig[k] = min(sig_rep_a[k],sig_rep_b[k])
	nb_s = { 'nc': ncmap_s, 'nnc': nncmap_s, 'hb': hbmap_s, 'sig': sig } #, 'details': details }
	return nb_s

def average_bonded(list1,list2):
	if len(list1) != len(list2):	# we definitely have a problem
		sys.stderr.write("# bonds/angles in chain a != # bonds/angles in chain b!\n")
		sys.exit(1)
	average = []
	nres1 = len(list1)+len(list1[0])-2
	for i in range(len(list1)):
		average.append(list1[i][:-1]+((list1[i][-1]+list2[i][-1])/2.0,))
	for i in range(len(list1)):
		average.append(tuple(map(lambda x: x+nres1,list1[i][:-1]))+((list1[i][-1]+list2[i][-1])/2.0,))
	return average

def average_bonded_list(lists):
	chains = lists.keys()
	chains.sort()
	nchain = len(chains)

	#if len(list1) != len(list2):	# we definitely have a problem
	#	sys.stderr.write("# bonds/angles in chain a != # bonds/angles in chain b!\n")
	#	sys.exit(1)
	average = []
	list1 = lists[chains[0]]
	nres1 = len(list1)+len(list1[0])-2
	for i in range(len(list1)):
		ave = 0.0
		for c in chains:
			ave += lists[c][i][-1]
		ave /= float(nchain)
		average.append(list1[i][:-1]+(ave,))
	return average

def merge_bonded(lists):
	merged = []
	keys = lists.keys()
	keys.sort()

	nresp = 0
	for k in keys:
		for i in lists[k]:
			merged.append(tuple(map(lambda x: x+nresp,i[:-1]))+(i[-1],))
		nresp += len(lists[k])+len(lists[k][0])-2
	return merged

def merge_dihe(lists):
	merged = []
	keys = lists.keys()
	keys.sort()

	nresp = 0
	for k in keys:
		for i in lists[k]:
			merged.append(tuple(map(lambda x: "G"+str(int(x[1:])+nresp),i[:-1]))+(i[-1],))
		nresp += len(lists[k])+len(lists[k][0])-2
	return merged

def nmer_merge_bonded(lists):
	merged = []
	nres1 = len(lists[0])+len(lists[0][0])-2
	for l in range(len(lists)):
		for i in l:
			merged.append(tuple(map(lambda x: x+l*nres1,i[:-1]))+(i[-1],))
	return merged

def min_sig(siga,sigb):
	sigm = {}
	for k in siga.keys():
		if sigb[k] < siga[k]:
			sigm[k] = sigb[k]
		else:
			sigm[k] = siga[k]
	return sigm

def cat_sig(siga,sigb):
	sigm = {}
	len_a = len(siga)
	for k in siga.keys():
		sigm[k] = siga[k]
	for k in sigb.keys():
		kk = k+len_a
		sigm[kk] = sigb[k]
	return sigm

def merge_sig(siga,sigb,sig_ab_a,sig_ab_b):
	sigm = {}
	len_a = len(siga)
	for k in siga.keys():
		if siga[k] < sig_ab_a[k]:
			sigm[k] = siga[k]
		else:
			sigm[k] = sig_ab_a[k]
	for k in sigb.keys():
		kk = k+len_a
		if sigb[k] < sig_ab_b[k]:
			sigm[kk] = sigb[k]
		else:
			sigm[kk] = sig_ab_b[k]
	return sigm

def minsig(siga,sigb):
	msig = {}
	for k in siga.keys():
		msig[k] = min(siga[k],sigb[k])
	return msig

def merge_sig_matrix(nbmat):
	chains = nbmat.keys()
	chains.sort()
	nchain = len(chains)
	msig = {}
	for i in range(nchain):
		for j in range(i,nchain):
			ci = chains[i]
			cj = chains[j]
			if ci == cj:
				matsig = nbmat[ci][cj]['sig']
			else:
				matsig = minsig(nbmat[ci][cj]['siga'],\
						nbmat[ci][cj]['sigb'])
			for k in matsig.keys():
				if k not in msig.keys():
					msig[k] = matsig[k]
				elif matsig[k] < msig[k]:
					msig[k] = matsig[k]
	return msig

def merge_sig_hetero(siga,sigb):
	sigm = {}
	len_a = len(siga)
	for k in siga.keys():
		sigm[k] = siga[k]
	for k in sigb.keys():
		kk = k+len_a
		sigm[kk] = sigb[k]
	return sigm

def common_sig(sig):
	sigc = {}
	nres = len(sig.keys())/2.0
	for i in range(1,nres+1):
		if sig[i] < sig[i+nres]:
			d = sig[i]
		else:
			d = sig[i+nres]
		sigc[i] = d
		sigc[i+nres] = d
	return sigc

def merge_emaps_matrix(maps,xyz):
	chains = xyz.keys()
	chains.sort()
	nchain = len(chains)
	offsets = {}
	off = 0
	for chain in chains:
		l = len(xyz[chain])
		offsets[chain] = off
		off += l
	mapm = {}
	for i in range(nchain):
		chaini = chains[i]
		offi = offsets[chaini]
		#for k in maps[chaini][chaini].keys():
		#	mapm[k+offi,k+offi] = maps[chaini][chaini][k]
		for j in range(i,nchain):
			chainj = chains[j]
			offj = offsets[chainj]
			for k in maps[chaini][chainj].keys():
				mapm[k+offi,k+offj] = maps[chaini][chainj][k]
	return mapm

def merge_emaps(mapa,mapb,mapab,nres_a):
	mapm = {}
	for k in mapa.keys():
		mapm[k] = mapa[k]
	for k in mapb.keys():
		kk = (k[0]+nres_a,k[1]+nres_a)
		mapm[kk] = mapb[k]
	for k in mapab:
		kk = (k[0],k[1]+nres_a)
		mapm[kk] = mapab[k]
	return mapm

def merge_emap_hetero(mapa,mapb,nres_a):
	mapm = {}
	for k in mapa.keys():
		mapm[k] = mapa[k]
	for k in mapb.keys():
		kk = (k[0]+nres_a,k[1]+nres_a)
		mapm[kk] = mapb[k]
	return mapm

def average_emaps(mapa,mapb,mapab,nres_a):
	mapm = {}
	# average intramolecular contacts
	jointkeys = mapa.keys()
	for k in mapb.keys():
		if k not in jointkeys:
			jointkeys.append(k)
	for k in jointkeys:
		kk = (k[0]+nres_a,k[1]+nres_a)
		if k in mapa.keys():
			if k in mapb.keys():
				val = [(mapa[k][0]+mapb[k][0])/2.0,(mapa[k][1]+mapb[k][1])/2.0]
			else:
				val = mapa[k]
		else:
			val = mapb[k]
		mapm[k] = val
		mapm[kk] = val
	# average intermolecular contacts
	for i in range(1,nres_a+1):
		for j in range(1,i+1):
			k = (i,j)
			kk = (j,i)
			p = (i,j+nres_a)
			q = (j,i+nres_a)
			if k not in mapab.keys() and kk not in mapab.keys():
				continue
			if k in mapab.keys():
				if kk in mapab.keys():
					val = [(mapab[k][0]+mapab[kk][0])/2.0,(mapab[k][1]+mapab[kk][1])/2.0]
				else:
					val = mapab[k]
			else:
				val = mapab[kk]
			mapm[p] =  val
			mapm[q] =  val
	return mapm

def average_emaps_matrix(maps,type):
	mapm = {}
	chains = maps.keys()
	chains.sort()
	nchain = len(chains)
	for i in range(nchain):
		ci = chains[i]
		for j in range(i,nchain):
			cj = chains[j]
			emap = maps[ci][cj][type]
			for k in emap.keys():
				if k not in mapm.keys():
					mapm[k] = [ emap[k][0],emap[k][1], 1]
				else:
					mapm[k][0] += emap[k][0]
					mapm[k][1] += emap[k][1]
					mapm[k][2] += 1
	for k in mapm.keys():
		nn = float(mapm[k][2])
		mapm[k] = [ mapm[k][0]/nn, mapm[k][1]/nn ]
	return mapm

def merge_emaps_matrix(maps,xyz):
	chains = xyz.keys()
	chains.sort()
	nchain = len(chains)
	offsets = {}
	off = 0
	for chain in chains:
		l = len(xyz[chain])
		offsets[chain] = off
		off += l
	mapm = {}
	for i in range(nchain):
		chaini = chains[i]
		offi = offsets[chaini]
		#for k in maps[chaini][chaini].keys():
		#	mapm[k+offi,k+offi] = maps[chaini][chaini][k]
		for j in range(i,nchain):
			chainj = chains[j]
			offj = offsets[chainj]
			for k in maps[chaini][chainj].keys():
				mapm[k+offi,k+offj] = maps[chaini][chainj][k]
	return mapm

def merge_kh(chains,xyz,khmap):
	globalkh = {}
	#chains = xyz.keys()
	chains.sort()
	offsets = {}
	off = 0
	nchain = len(chains)
	for chain in chains:
		l = len(xyz[chain])
		offsets[chain] = off
		off += l
	idx = -1
	for i in range(nchain):
		chaini = chains[i]
		if chaini in [ 'CROW', 'SOLV' ]:
			continue
		offi = offsets[chaini]
		#for k in khmap[chaini][chainj].keys():
		#	globalkh[(k[0]+offi,k[1]+offj)] = khmap[chaini][chainj][k]
		for j in range(i,nchain):
			idx += 1
			chainj = chains[j]
			if chainj in [ 'CROW', 'SOLV' ]:
				continue
			idx += 1
			offj = offsets[chainj]
			#
			for k in khmap[chaini][chainj].keys():
				globalkh[(k[0]+offi,k[1]+offj)] = khmap[chaini][chainj][k]
	return globalkh

def merge_kh_inter(chains,xyz,khmap):
	globalkh = {}
	#chains = xyz.keys()
	chains.sort()
	offsets = {}
	off = 0
	nchain = len(chains)
	for chain in chains:
		l = len(xyz[chain])
		offsets[chain] = off
		off += l
	idx = -1
	for i in range(nchain):
		chaini = chains[i]
		# should not be necessary...
		if chaini in [ 'CROW', 'SOLV' ]:
			continue
		offi = offsets[chaini]
		for j in range(i+1,nchain):
			idx += 1
			chainj = chains[j]
			if chainj in [ 'CROW', 'SOLV' ]:
				continue
			offj = offsets[chainj]
			#print chaini, chainj
			#
			for k in khmap[chaini][chainj].keys():
				globalkh[(k[0]+offi,k[1]+offj)] = khmap[chaini][chainj][k]
	return globalkh

def merge_nonbonded_matrix(nb,xyz,Tfs):
	#
	# Tfs is a list of folding temperatures, one for each chain
	# i.e. excluding backbone. If len(Tfs)==1, a global scaling
	# is done. Otherwise, len(Tfs) must = nchain*(nchain+1)/2, and a per
	# chain scaling is done
	chains = filter(lambda x: x not in [ 'CROW', 'SOLV' ],xyz.keys())
	#print chains
	chains.sort()
	nchain = len(chains)
	offsets = {}
	off = 0
	#
	NSFAC = 0.0
	# determine whether scaling to a global or local Tf
	if len(Tfs) == 1:
		globalscale = 1
	elif len(Tfs) != nchain*(nchain+1)/2:
		sys.stderr.write("merge_nonbonded_matrix: len(Tfs) != nchain*(nchain+1)/2\n")
	else:
		globalscale = 0
	#
	for chain in chains:
		l = len(xyz[chain])
		offsets[chain] = off
		off += l
	if globalscale:
		esum = 0.0
		totres = 0
		# first determine sum of native pair interactions
		for i in range(nchain):
			chaini = chains[i]
			li = len(xyz[chaini])
			totres += li
			for k in nb[chaini][chaini]['nc'].keys():
				esum += nb[chaini][chaini]['nc'][k][1]
			for k in nb[chaini][chaini]['hb'].keys():
				esum += nb[chaini][chaini]['hb'][k][1]
			for j in range(i+1,nchain):
				chainj = chains[j]
				#print chaini, chainj
				for k in nb[chaini][chainj]['nc'].keys():
					esum += nb[chaini][chainj]['nc'][k][1]
				for k in nb[chaini][chainj]['hb'].keys():
					esum += nb[chaini][chainj]['hb'][k][1]
		eps_res = big_fat_fudge_factor*Tfs[0]
		if abs(esum) < 1.e-15:
			nat_scale_fac = 1.0
		else:
			nat_scale_fac = abs(float(totres)*eps_res/esum)
		NSFAC = nat_scale_fac
		#print NSFAC
	#
	#
	gkhmap = {}
	ncmap = {}
	ncintra = {}
	hbintra = {}
	nncmap = {}
	hbmap = {}
	sig = {}
	sigintra = {}
	idx = -1
	for i in range(nchain):
		idx += 1
		chaini = chains[i]
		nres = len(xyz[chaini])
		offi = offsets[chaini]
		for k in nb[chaini][chaini]['sig'].keys():
			kk = k + offi
			sigk = nb[chaini][chaini]['sig'][k]
			if kk not in sig:
				sig[kk] = sigk
			elif sigk < sig[kk]:
				sig[kk] = sigk
			sigintra[kk] = sig[kk]
		#
		if not globalscale:
			# calc per-chain scale factor
			esum = 0.0
			for k in nb[chaini][chaini]['nc'].keys():
				esum += nb[chaini][chaini]['nc'][k][1]
			for k in nb[chaini][chaini]['hb'].keys():
				esum += nb[chaini][chaini]['hb'][k][1]
			eps_res = big_fat_fudge_factor*Tfs[idx]
			nat_scale_fac = abs(float(nres)*eps_res/esum)
			if i == 0:
				NSFAC= nat_scale_fac
		#
		for k in nb[chaini][chaini]['nc'].keys():
			rij,eij = nb[chaini][chaini]['nc'][k]
			eij *= nat_scale_fac
			ncmap[k[0]+offi,k[1]+offi] = [ rij, eij ]
			#ncmap[k[0]+offi,k[1]+offi][1] *= nat_scale_fac
			#ncintra[k[0]+offi,k[1]+offi] = nb[chaini][chaini]['nc'][k]
			ncintra[k[0]+offi,k[1]+offi] = [ rij, eij ]
			#ncintra[k[0]+offi,k[1]+offi][1] *= nat_scale_fac
		for k in nb[chaini][chaini]['nnc'].keys():
			nncmap[k[0]+offi,k[1]+offi] = nb[chaini][chaini]['nnc'][k]
			nncmap[k[0]+offi,k[1]+offi][1] *= nat_scale_fac
		for k in nb[chaini][chaini]['hb'].keys():
			rij,eij = nb[chaini][chaini]['hb'][k]
			eij *= nat_scale_fac
			hbmap[k[0]+offi,k[1]+offi] = [rij,eij]
			#hbmap[k[0]+offi,k[1]+offi][1] *= nat_scale_fac
			hbintra[k[0]+offi,k[1]+offi] = [rij,eij]
			#hbintra[k[0]+offi,k[1]+offi] = nb[chaini][chaini]['hb'][k]
			#hbintra[k[0]+offi,k[1]+offi][1] *= nat_scale_fac
		for j in range(i+1,nchain):
			idx += 1
			chainj = chains[j]
			offj = offsets[chainj]
			#print chaini, chainj
			#
			if not globalscale:
				# calc per-interface scale factor
				esum = 0.0
				pkeys = []
				qkeys = []
				for k in nb[chaini][chaini]['nc'].keys():
					p,q = k
					if p not in pkeys:
						pkeys.append(p)
					if q not in qkeys:
						qkeys.append(q)
					esum += nb[chaini][chaini]['nc'][k][1]
				for k in nb[chaini][chaini]['hb'].keys():
					p,q = k
					if p not in pkeys:
						pkeys.append(p)
					if q not in qkeys:
						qkeys.append(q)
					esum += nb[chaini][chaini]['hb'][k][1]
				eps_res = big_fat_fudge_factor*Tfs[idx]
				nres_interface = len(pkeys)+len(qkeys)
				nat_scale_fac = abs(float(nres_interface)*eps_res/esum)
			#
			for k in nb[chaini][chainj]['nc'].keys():
				ncmap[k[0]+offi,k[1]+offj] = nb[chaini][chainj]['nc'][k]
				ncmap[k[0]+offi,k[1]+offj][1] *= nat_scale_fac
			for k in nb[chaini][chainj]['nnc'].keys():
				nncmap[k[0]+offi,k[1]+offj] = nb[chaini][chainj]['nnc'][k]
				nncmap[k[0]+offi,k[1]+offj][1] *= nat_scale_fac
			for k in nb[chaini][chainj]['hb'].keys():
				hbmap[k[0]+offi,k[1]+offj] = nb[chaini][chainj]['hb'][k]
				hbmap[k[0]+offi,k[1]+offj][1] *= nat_scale_fac
			for k in nb[chaini][chainj]['siga'].keys():
				kk = k + offi
				sigk = nb[chaini][chainj]['siga'][k]
				if kk not in sig:
					sig[kk] = sigk
				elif sigk < sig[kk]:
					sig[kk] = sigk
			for k in nb[chaini][chainj]['sigb'].keys():
				kk = k + offj
				sigk = nb[chaini][chainj]['sigb'][k]
				if kk not in sig:
					sig[kk] = sigk
				elif sigk < sig[kk]:
					sig[kk] = sigk

	nonbonded = { 'nc': ncmap, 'ncintra': ncintra, 'nnc': nncmap, 'hb': hbmap,
			'hbintra': hbintra, 'sig': sig, 'sigintra': sigintra }
	return nonbonded,NSFAC

def merge_threebody_matrix(threebody_in,xyz):
	#
	chains = filter(lambda x: x!='SOLV',xyz.keys())
	chains.sort()
	nchain = len(chains)
	offsets = {}
	off = 0
	for chain in chains:
		l = len(xyz[chain])
		offsets[chain] = off
		off += l
	threebody = {}
	for i in range(nchain):
		chaini = chains[i]
		offi = offsets[chaini]
		for k in threebody_in[chaini][chaini].keys():
			threebody[k[0]+offi,k[1]+offi,k[2]+offi] = threebody_in[chaini][chaini][k]
	return threebody

def write_threebody(file,threebody):
	outp = open(file,"w")
	for key in threebody.keys():
		i,j,k = key
		outp.write("%5i %5i %5i %8.3f %8.3f %8.3f %8.3f\n"%(i,j,k,threebody[key][0],
			threebody[key][1],threebody[key][2],threebody[key][3]))
	outp.close()

def hetero_mergenb(nba,nbb,crda,crdb,ngdist,gamma):
	#
	ncmap_a = nba['nc']
	nncmap_a = nba['nnc']
	hbmap_a = nba['hb']
	sig_a = nba['sig']
	#details_a = nba['details']
	#
	ncmap_b = nbb['nc']
	nncmap_b = nbb['nnc']
	hbmap_b = nbb['hb']
	sig_b = nbb['sig']
	#details_b = nbb['details']
	#
	len_a = len(sig_a)
	len_b = len(sig_b)

	ncmap = merge_emap_hetero(ncmap_a,ncmap_b,len_a)
	nncmap = merge_emap_hetero(nncmap_a,nncmap_b,len_a)
	hbmap = merge_emap_hetero(hbmap_a,hbmap_b,len_a)
	sig = merge_sig_hetero(sig_a,sig_b)
	#details = details_a
	mjmap, ave_mj = read_miyazawa_jernigan()
	if gamma > 0.001:
		for i in range(1,len_a+1):
			for j in range(1,len_b+1):
				mij = mjmap[crda[i]["name"]][crdb[j]["name"]]/(ave_mj*2.0)*gamma
				nncmap[(i,j+len_a)] = [ngdist,mij]
	nonbonded = { 'nc': ncmap, 'nnc': nncmap, 'hb': hbmap, 'sig': sig }#, 'details': details }
	return nonbonded

def average_nonbonded_matrix(nbmat):
	#
	ncmap = average_emaps_matrix(nbmat,'nc')
	nncmap = average_emaps_matrix(nbmat,'nnc')
	hbmap = average_emaps_matrix(nbmat,'hb')
	sigm = merge_sig_matrix(nbmat)
	nonbonded = { 'nc': ncmap, 'nnc': nncmap, 'hb': hbmap, 'sig': sigm }#, 'details': details }
	return nonbonded

def average_nonbonded(nba,nbb,nbab):
	#
	ncmap_a = nba['nc']
	nncmap_a = nba['nnc']
	hbmap_a = nba['hb']
	sig_a = nba['sig']
	#details_a = nba['details']
	#
	ncmap_b = nbb['nc']
	nncmap_b = nbb['nnc']
	hbmap_b = nbb['hb']
	sig_b = nbb['sig']
	#details_b = nbb['details']
	#
	ncmap_ab = nbab['nc']
	nncmap_ab = nbab['nnc']
	hbmap_ab = nbab['hb']
	sig_ab_a = nbab['siga']
	sig_ab_b = nbab['sigb']
	#
	len_a = len(sig_a)

	ncmap = average_emaps(ncmap_a,ncmap_b,ncmap_ab,len_a)
	nncmap = average_emaps(nncmap_a,nncmap_b,nncmap_ab,len_a)
	hbmap = average_emaps(hbmap_a,hbmap_b,hbmap_ab,len_a)
	sigm = merge_sig(sig_a,sig_b,sig_ab_a,sig_ab_b)
	sig = common_sig(sigm)
	#details = details_a
	nonbonded = { 'nc': ncmap, 'nnc': nncmap, 'hb': hbmap, 'sig': sig }#, 'details': details }
	return nonbonded

def calc_bonds(crds,link):
	bondlist = []
	nres = len(crds.keys())
	na = isNA(crds)
	if na:
		for b in range(nres-1):
			i,j = b+1,b+2
			rij = dist(crds[i]["C1'"],crds[j]["C1'"])
			bondlist.append((i,j,rij))
	else:
		for b in range(nres-1):
			i,j = b+1,b+2
			rij = dist(crds[i]["CA"],crds[j]["CA"])
			bondlist.append((i,j,rij))
	if link:
		bondlist.append((nres,1,std_bond))
	return bondlist

def calc_sc_bonds(crds,link):
	bondlist = []
	nres = len(crds.keys())
	residues = crds.keys()
	residues.sort()
	na = isNA(crds)
	if na:
		# backbone bonds
		for r in residues[:-1]:
			rij = dist(crds[r]["C1'"],crds[r+1]["C1'"])
			bondlist.append(("G%i"%(r),"G%i"%(r+1),rij))
		if link:
			bondlist.append(("G%i"%(nres),"G1",std_bond))
		# sidechain bonds
		for r in residues:
			type = crds[r]['name'].strip()
			if type in scmassmap.keys():
				bondlist.append(("G%i"%(r),"G%iS"%(r),scbondmap[type]))
	else:
		# backbone bonds
		for r in residues[:-1]:
			rij = dist(crds[r]["CA"],crds[r+1]["CA"])
			bondlist.append(("G%i"%(r),"G%i"%(r+1),rij))
		if link:
			bondlist.append(("G%i"%(nres),"G1",std_bond))
		# sidechain bonds
		for r in residues:
			type = crds[r]['name']
			if type in scmassmap.keys():
				bondlist.append(("G%i"%(r),"G%iS"%(r),scbondmap[type]))
	return bondlist

def calc_angles(crds,link,stata=0,urey=0):
	anglelist = []
	nres = len(crds.keys())
	na = isNA(crds)
	if na:
		for b in range(nres-2):
			i,j,k = b+1,b+2,b+3
			if stata:
				anglelist.append((i,j,k,statangle['BACKBONE']))
			elif urey:
				theta_ijk = angle(crds[i]["C1'"],crds[j]["C1'"],crds[k]["C1'"])
				r_ik = dist(crds[i]["C1'"],crds[k]["C1'"])
				anglelist.append((i,j,k,[0.0,theta_ijk,50.0,r_ik]))
			else:
				theta_ijk = angle(crds[i]["C1'"],crds[j]["C1'"],crds[k]["C1'"])
				anglelist.append((i,j,k,[theta_ijk]))
	else:
		for b in range(nres-2):
			i,j,k = b+1,b+2,b+3
			if stata:
				anglelist.append((i,j,k,statangle['BACKBONE']))
			elif urey:
				theta_ijk = angle(crds[i]["CA"],crds[j]["CA"],crds[k]["CA"])
				r_ik = dist(crds[i]["CA"],crds[k]["CA"])
				anglelist.append((i,j,k,[0.0,theta_ijk,50.0,r_ik]))
			else:
				theta_ijk = angle(crds[i]["CA"],crds[j]["CA"],crds[k]["CA"])
				anglelist.append((i,j,k,[theta_ijk]))
	if link:
		if stata:
			anglelist.append((nres-1,nres,1,statangle['BACKBONE']))
			anglelist.append((nres,1,2,statangle['BACKBONE']))
		elif urey:
			anglelist.append((i,j,k,[0.0,std_angle,50.0,std_urey]))
			anglelist.append((nres-1,nres,1,[0.0,std_angle,50.0,std_urey]))
			anglelist.append((nres,1,2,[0.0,std_angle,50.0,std_urey]))
		else:
			anglelist.append((nres-1,nres,1,[std_angle]))
			anglelist.append((nres,1,2,[std_angle]))
	return anglelist

def calc_sc_angles(crds,link,stata=0):
	anglelist = []
	nres = len(crds.keys())
	if isNA(crds):
		return calc_sc_angles_na(crds,link,stata)
	# backbone angles
	for b in range(nres-2):
		i,j,k = b+1,b+2,b+3
		if stata:
			anglelist.append(("G%i"%(i),"G%i"%(j),"G%i"%(k),statangle['BACKBONE']))
		else:
			theta_ijk = angle(crds[i]["CA"],crds[j]["CA"],crds[k]["CA"])
			anglelist.append(("G%i"%(i),"G%i"%(j),"G%i"%(k),[theta_ijk]))
	if link:
		if stata:
			anglelist.append(("G%i"%(nres-1),"G%i"%(nres),"G%i"%(1),statangle['BACKBONE']))
			anglelist.append(("G%i"%(nres),"G%i"%(1),"G%i"%(2),statangle['BACKBONE']))
		else:
			anglelist.append(("G%i"%(nres-1),"G%i"%(nres),"G%i"%(1),[std_angle]))
			anglelist.append(("G%i"%(nres),"G%i"%(1),"G%i"%(2),[std_angle]))
	# sidechain angles
	for b in range(nres-1):		# forward angles
		i = b+1
		j = b+2
		if stata:
			if "SCCENT" in crds[j].keys():
				name = crds[j]['name']
				anglelist.append(("G%i"%(i),"G%i"%(j),"G%iS"%(j),statangle[name+str(1)]))
			if "SCCENT" in crds[i].keys():
				name = crds[i]['name']
				anglelist.append(("G%iS"%(i),"G%i"%(i),"G%i"%(j),statangle[name+str(2)]))
		else:
			if "SCCENT" in crds[j].keys():
				theta = angle(crds[i]["CA"],crds[j]["CA"],crds[j]["SCCENT"])
				anglelist.append(("G%i"%(i),"G%i"%(j),"G%iS"%(j),[theta]))
			if "SCCENT" in crds[i].keys():
				theta = angle(crds[i]["SCCENT"],crds[i]["CA"],crds[j]["CA"])
				anglelist.append(("G%iS"%(i),"G%i"%(i),"G%i"%(j),[theta]))
	if link:
		if stata:
			if "SCCENT" in crds[1].keys():
				name = crds[1]['name']
				anglelist.append(("G%i"%(nres),"G%i"%(1),"G%iS"%(1),statangle[name+str(1)]))
			if "SCCENT" in crds[nres].keys():
				name = crds[nres]['name']
				anglelist.append(("G%iS"%(nres),"G%i"%(nres),"G%i"%(1),statangle[name+str(2)]))
		else:
			anglelist.append(("G%iS"%(nres),"G%i"%(nres),"G1",[std_angle]))
			anglelist.append(("G%i"%(nres),"G1","G1S",[std_angle]))
	return anglelist

def calc_sc_angles_na(crds,link,stata=0):
	anglelist = []
	nres = len(crds.keys())
	# backbone angles
	for b in range(nres-2):
		i,j,k = b+1,b+2,b+3
		if stata:
			anglelist.append(("G%i"%(i),"G%i"%(j),"G%i"%(k),statangle['BACKBONE']))
		else:
			theta_ijk = angle(crds[i]["C1'"],crds[j]["C1'"],crds[k]["C1'"])
			anglelist.append(("G%i"%(i),"G%i"%(j),"G%i"%(k),[theta_ijk]))
	if link:
		if stata:
			anglelist.append(("G%i"%(nres-1),"G%i"%(nres),"G%i"%(1),statangle['BACKBONE']))
			anglelist.append(("G%i"%(nres),"G%i"%(1),"G%i"%(2),statangle['BACKBONE']))
		else:
			anglelist.append(("G%i"%(nres-1),"G%i"%(nres),"G%i"%(1),[std_angle]))
			anglelist.append(("G%i"%(nres),"G%i"%(1),"G%i"%(2),[std_angle]))
	# sidechain angles
	for b in range(nres-1):		# forward angles
		i = b+1
		j = b+2
		if stata:
			if "SCCENT" in crds[j].keys():
				name = crds[j]['name']
				anglelist.append(("G%i"%(i),"G%i"%(j),"G%iS"%(j),statangle[name+str(1)]))
			if "SCCENT" in crds[i].keys():
				name = crds[i]['name']
				anglelist.append(("G%iS"%(i),"G%i"%(i),"G%i"%(j),statangle[name+str(2)]))
		else:
			if "SCCENT" in crds[j].keys():
				theta = angle(crds[i]["C1'"],crds[j]["C1'"],crds[j]["SCCENT"])
				anglelist.append(("G%i"%(i),"G%i"%(j),"G%iS"%(j),[theta]))
			if "SCCENT" in crds[i].keys():
				theta = angle(crds[i]["SCCENT"],crds[i]["C1'"],crds[j]["C1'"])
				anglelist.append(("G%iS"%(i),"G%i"%(i),"G%i"%(j),[theta]))
	if link:
		if stata:
			if "SCCENT" in crds[1].keys():
				name = crds[1]['name']
				anglelist.append(("G%i"%(nres),"G%i"%(1),"G%iS"%(1),statangle[name+str(1)]))
			if "SCCENT" in crds[nres].keys():
				name = crds[nres]['name']
				anglelist.append(("G%iS"%(nres),"G%i"%(nres),"G%i"%(1),statangle[name+str(2)]))
		else:
			anglelist.append(("G%iS"%(nres),"G%i"%(nres),"G1",[std_angle]))
			anglelist.append(("G%i"%(nres),"G1","G1S",[std_angle]))
	return anglelist

def strip_seg(ds):
	ds

#def calc_disulf_bonded(crds,disulf_str):
def calc_disulf_bonded(crds,disulfides):
	ss_bonds = []
	ss_angles = []
	ss_dihe = []
	#disulfides = disulf_str.split(':')

	nres = max(crds.keys())
	for ss in disulfides:
		ii,jj = ss.split('-')
		i, j = int(ii), int(jj)
		rij = dist(crds[i]["CA"],crds[j]["CA"])
		ss_bonds.append((i,j,rij))
		im1 = i-1
		ip1 = i+1
		jm1 = j-1
		jp1 = j+1
		# im1 - i - j
		if im1 >= 1:
			theta = angle(crds[im1]["CA"],crds[i]["CA"],crds[j]["CA"])
			ss_angles.append((im1,i,j,[theta]))
		# ip1 - i - j
		if ip1 <= nres:
			theta = angle(crds[ip1]["CA"],crds[i]["CA"],crds[j]["CA"])
			ss_angles.append((ip1,i,j,[theta]))
		# jm1 - j - i
		if jm1 >= 1:
			theta = angle(crds[jm1]["CA"],crds[j]["CA"],crds[i]["CA"])
			ss_angles.append((jm1,j,i,[theta]))
		# jp1 - j - i
		if jp1 <= nres:
			theta = angle(crds[jp1]["CA"],crds[j]["CA"],crds[i]["CA"])
			ss_angles.append((jp1,j,i,[theta]))
		ss_dihe.append(("X","G%i"%(i),"G%i"%(j),"X",\
				[(0.001,3,0.0)]))
	return ss_bonds, ss_angles, ss_dihe

def setup_dihe(crds,link,fixalpha=0,zero=0,native=0):
	dihelist = []
	nres = len(crds.keys())
	kdihe = read_karanicolas_dihe()
	if fixalpha:
		kdihe = tweak_dihe(kdihe,-1.158288,287.35483)
	if isNA(crds):
		zero=1
	for b in range(nres-3):
		if zero:
			i,j,k,l = b+1,b+2,b+3,b+4
			dihelist.append(("G%i"%i,"G%i"%j,"G%i"%k,"G%i"%l,\
					[[1.0e-15,3,0.0]]))
		#elif isNA(crds):
		#	i,j,k,l = b+1,b+2,b+3,b+4
		#	ri = crds[i]["C1'"]
		#	rj = crds[j]["C1'"]
		#	rk = crds[k]["C1'"]
		#	rl = crds[l]["C1'"]
		#	phi = torsion(ri,rj,rk,rl)
		#	dihelist.append(("G%i"%i,"G%i"%j,"G%i"%k,"G%i"%l,\
		#			[[2.0,2,phi]]))
		elif native==1:
			# note won't work for nucleic acids!!!
			i,j,k,l = b+1,b+2,b+3,b+4
			ri = crds[i]["CA"]
			rj = crds[j]["CA"]
			rk = crds[k]["CA"]
			rl = crds[l]["CA"]
			phi = torsion(ri,rj,rk,rl)
			phi2 = 2.*phi
			if phi>0.:
				phi-=180.
			else:
				phi+=180.
			if phi2>0.:
				phi2-=180.
			else:
				phi2+=180.
			dihelist.append(("G%i"%i,"G%i"%j,"G%i"%k,"G%i"%l,\
					[[1.0,1,phi],[2.0,2,phi2]]))
		elif native==2:
			# note won't work for nucleic acids!!!
			i,j,k,l = b+1,b+2,b+3,b+4
			ri = crds[i]["CA"]
			rj = crds[j]["CA"]
			rk = crds[k]["CA"]
			rl = crds[l]["CA"]
			phi = torsion(ri,rj,rk,rl)
			phi2 = 2.*phi
			if phi>0.:
				phi-=180.
			else:
				phi+=180.
			if phi2>0.:
				phi2-=180.
			else:
				phi2+=180.
			dihelist.append(("G%i"%i,"G%i"%j,"G%i"%k,"G%i"%l,\
					[[2.0,2,phi2]]))
		else:
			i,j,k,l = b+1,b+2,b+3,b+4
			dihelist.append(("G%i"%i,"G%i"%j,"G%i"%k,"G%i"%l,\
					kdihe[crds[j]['name']][crds[k]['name']]))
	if link:
		for b in range(3):
			i,j,k,l = b-2, b-1, b, b+1
			if i<1:
				i = nres + i
			if j<1:
				j = nres + j
			if k<1:
				k = nres + k
			if l<1:
				l = nres + l
			if zero:
				dihelist.append(("G%i"%i,"G%i"%j,"G%i"%k,"G%i"%l,
					[[1.0e-15,3,0.0]]))
			else:
				dihelist.append(("G%i"%i,"G%i"%j,"G%i"%k,"G%i"%l,
					kdihe[crds[j]['name']][crds[k]['name']]))
	return dihelist

def setup_sc_dihe(crds,link,fixalpha=0):
	dihelist = []
	na = isNA(crds)
	nres = len(crds.keys())
	kdihe = read_karanicolas_dihe()
	if fixalpha:
		kdihe = tweak_dihe(kdihe,-1.158288,287.35483)
	for b in range(nres-3):
		i,j,k,l = b+1,b+2,b+3,b+4
		if na:
			ri = crds[i]["C1'"]
			rj = crds[j]["C1'"]
			rk = crds[k]["C1'"]
			rl = crds[l]["C1'"]
			phi = torsion(ri,rj,rk,rl)
			if phi>0.:
				phase = phi-180.
				phase2 = 2.*phi-180
				if phase2 < -180.:
					phase2+=360.
				elif phase2 > 180.:
					phase2-=360.
			else:
				phase = phi+180.
				phase2 = 2.*phi+180
				if phase2 < -180.:
					phase2+=360.
				elif phase2 > 180.:
					phase2-=360.
			dihelist.append(("G%i"%i,"G%i"%j,"G%i"%k,"G%i"%l,\
					[[0.35,1,phase]]))
			dihelist.append(("G%i"%i,"G%i"%j,"G%i"%k,"G%i"%l,\
					[[0.5,2,phase2]]))
			#dihelist.append(("G%i"%i,"G%i"%j,"G%i"%k,"G%i"%l,\
			#		[[1.0e-15,3,0.0]]))
			#print "phi",i,j,k,l,phi
		else:
			dihelist.append(("G%i"%(i),"G%i"%(j),"G%i"%(k),"G%i"%(l),\
					kdihe[crds[j]['name']][crds[k]['name']]))
	if link:
		for b in range(3):
			i,j,k,l = b-2, b-1, b, b+1
			if i<1:
				i = nres + i
			if j<1:
				j = nres + j
			if k<1:
				k = nres + k
			if l<1:
				l = nres + l
			dihelist.append(("G%i"%(i),"G%i"%(j),"G%i"%(k),"G%i"%(l),\
					kdihe[crds[j]['name']][crds[k]['name']]))
	for i in range(nres-1):
		j = i+1
		k = i+2
		dihelist.append(("X","G%i"%(j),"G%i"%(k),"X",\
				[(0.001,3,0.0)]))
		if na:
			ri = crds[j]["SCCENT"]
			rj = crds[j]["C1'"]
			rk = crds[k]["C1'"]
			rl = crds[k]["SCCENT"]
			#phi = torsion(ri,rj,rk,rl)
			#dihelist.append(("G%iS"%j,"G%i"%j,"G%i"%k,"G%iS"%k,\
			#		[[2.0,2,phi]]))

	return dihelist

def setup_dihe_hlink(crds_a,crds_b,fixalpha=0):
	dihelist = []
	nres_a = len(crds_a.keys())
	nres_b = len(crds_b.keys())
	kdihe = read_karanicolas_dihe()
	if fixalpha:
		kdihe = tweak_dihe(kdihe,-1.158288,287.35483)
	types_a = map(lambda x: crds_a[x]['name'], range(1,nres_a+1))
	types_b = map(lambda x: crds_b[x]['name'], range(1,nres_b+1))
	types = types_a+types_b
	for b in range(3):
		i,j,k,l = nres_a+b-2, nres_a+b-1, nres_a+b, nres_a+b+1
		dihelist.append((i,j,k,l,kdihe[types[j-1]][types[k-1]]))
		#print i,j,k,l,types[j-1],types[k-1]
	for b in range(3):
		i,j,k,l = b-2, b-1, b, b+1
		if i<1:
			i = nres_a+nres_b + i
		if j<1:
			j = nres_a+nres_b + j
		if k<1:
			k = nres_a+nres_b + k
		if l<1:
			l = nres_a+nres_b + l
		#print i,j,k,l,types[j-1],types[k-1]
		dihelist.append((i,j,k,l,kdihe[types[j-1]][types[k-1]]))
	return dihelist

def nb2qlist(nonbonded):
	ncmap = nonbonded['nc']
	qlist = []
	for k in ncmap.keys():
		i,j=k
		qlist.append((i,j,ncmap[k][0]*qtol))
	return qlist

def nb2qlists(nonbonded,chains):
	"""create a Qlist map from a nc map"""
	nchain = len(chains)
	qlists = {}
	for i in range(nchain):
		chaini = chains[i]
		qlists[chaini] = {}
		qlists[chaini][chaini] = nb2qlist(nonbonded[chaini][chaini])
		for j in range(i+1,nchain):
			chainj = chains[j]
			qlists[chaini][chainj] = nb2qlist(nonbonded[chaini][chainj])
	return qlists

def make_go_nongo_map(nonbonded,nres):
	"""pool native contacts and hydrogen bonds into one list and
	non-native contacts into another; calculate native scaling
	factor"""
	ncmap = nonbonded['nc']
	nncmap = nonbonded['nnc']
	hbmap = nonbonded['hb']
	#details = nonbonded['details']
	go_map = {}
	nongo_map  = {}
	qlist = []
	tot_enat = 0.0
	for i in range(1,nres):
		for j in range(i,nres+1):
			if (i,j) in nncmap.keys():
				if (i,j) not in nongo_map.keys():
					nongo_map[(i,j)] = nncmap[(i,j)]
				else:
					nongo_map[(i,j)][1] += nncmap[(i,j)][1]

	#nat_scale_fac = abs(float(nres)*eps_res/tot_enat)
	#
	return go_map, nongo_map, qlist

def make_nongo_map(nncmap,nres,scscale=1.0,bbscale=1.0):
	"""pool native contacts and hydrogen bonds into one list and
	non-native contacts into another; calculate native scaling
	factor"""
	#details = nonbonded['details']
	go_map = {}
	nongo_map  = {}
	qlist = []
	tot_enat = 0.0
	for i in range(1,nres):
		for j in range(i,nres+1):
			if (i,j) in nncmap.keys():
				if (i,j) not in nongo_map.keys():
					nongo_map[(i,j)] = nncmap[(i,j)]
				else:
					nongo_map[(i,j)][1] += nncmap[(i,j)][1]

	#nat_scale_fac = abs(float(nres)*eps_res/tot_enat)
	#
	return nongo_map

def make_go_map(ncmap,hbmap,nres,scscale=1.0,bbscale=1.0):
	"""pool native contacts and hydrogen bonds into one list; calculate native scaling
	factor"""
	#details = nonbonded['details']
	go_map = {}
	nongo_map  = {}
	qlist = []
	tot_enat = 0.0
	for i in range(1,nres):
		for j in range(i,nres+1):
			if (i,j) in ncmap.keys():
				qlist.append((i,j,ncmap[(i,j)][0]*qtol))
				if (i,j) not in go_map.keys():
					go_map[(i,j)] = [ncmap[(i,j)][0],ncmap[(i,j)][1]*scscale]
				else:
					go_map[(i,j)][1] += ncmap[(i,j)][1]*scscale
				tot_enat += ncmap[(i,j)][1]
			if (i,j) in hbmap.keys():
				if (i,j) not in go_map.keys():
					go_map[(i,j)] = [hbmap[(i,j)][0],hbmap[(i,j)][1]*bbscale]
				else:
					go_map[(i,j)][1] += hbmap[(i,j)][1]*bbscale
				tot_enat += hbmap[(i,j)][1]

	#nat_scale_fac = abs(float(nres)*eps_res/tot_enat)
	#
	return go_map, qlist

def sc_mergemaps(nonbonded,eps_res,nres,atomind,scscale=1.0,bbscale=1.0):
	"""create nbfix for side-chain go-models"""
	ncmap = nonbonded['nc']
	nncmap = nonbonded['nnc']
	hbmap = nonbonded['hb']
	go_map = {}
	nongo_map  = {}
	qtol = 1.2	# lambda
	qlist = []
	allqlist = []
	tot_enat = 0.0
	nbfix = {}
	gomod_nbfix = {}	# for my GOMODEL module in CHARMM
	gomod_golist = {}
	#
	eps_res_sc = 1.0 * eps_res   # ad hoc correction for larger sc entropy
	#
	# first compute nat_scale_fac
	for i in range(1,nres):
		for j in range(i,nres+1):
			if (i,j) in ncmap.keys():
				tot_enat += ncmap[(i,j)][1]
			if (i,j) in hbmap.keys():
				tot_enat += hbmap[(i,j)][1]
	nat_scale_fac = abs(float(nres)*eps_res_sc/tot_enat)
	# then construct nbfix lists
	nsc, nhb = 0,0
	for i in range(1,nres):
		for j in range(i,nres+1):
			if (i,j) in ncmap.keys():
				nsc +=1
				if 'SCCENT' in atomind[i].keys():
					sc_i = atomind[i]['SCCENT']
					sca_i = "G%iS" % (i)
				else:
					sc_i = atomind[i]['CA']
					sca_i = "G%i" % (i)
				if 'SCCENT' in atomind[j].keys():
					sc_j = atomind[j]['SCCENT']
					sca_j = "G%iS" % (j)
				else:
					sc_j = atomind[j]['CA']
					sca_j = "G%i" % (j)
				qlist.append((sc_i,sc_j,ncmap[(i,j)][0]*qtol))
				allqlist.append((sc_i,sc_j,ncmap[(i,j)][0]*qtol))
				nbfix[(sca_i,sca_j)] = [ncmap[(i,j)][0],ncmap[(i,j)][1]*scscale*nat_scale_fac]
				gomod_golist[(sc_i,sc_j)] = [ncmap[(i,j)][0],ncmap[(i,j)][1]*scscale*nat_scale_fac]
				gomod_nbfix[(sca_i,sca_j)] = [ncmap[(i,j)][0],0.0]
			if (i,j) in hbmap.keys():
				nhb += 1
				c_i = atomind[i]['CA']
				c_j = atomind[j]['CA']
				ca_i = "G%i" %i
				ca_j = "G%i" %j
				if (ca_i,ca_j) not in nbfix.keys():
					nbfix[(ca_i,ca_j)] = [hbmap[(i,j)][0],hbmap[(i,j)][1]*bbscale*nat_scale_fac]
					gomod_nbfix[(ca_i,ca_j)] = [hbmap[(i,j)][0],0.0]
					gomod_golist[(c_i,c_j)] = [hbmap[(i,j)][0],hbmap[(i,j)][1]*bbscale*nat_scale_fac]
					allqlist.append((c_i,c_j,hbmap[(i,j)][0]*qtol))
				else:
					nbfix[(ca_i,ca_j)][1] += hbmap[(i,j)][1]*bbscale
					gomod_golist[(c_i,c_j)][1] += [hbmap[(i,j)][1]*bbscale*nat_scale_fac]
			if (i,j) in nncmap.keys():
				if 'SCCENT' in atomind[i].keys():
					sc_i = atomind[i]['SCCENT']
					sca_i = "G%iS" % (i)
				else:
					sc_i = atomind[i]['CA']
					sca_i = "G%i" % (i)
				if 'SCCENT' in atomind[j].keys():
					sc_j = atomind[j]['SCCENT']
					sca_j = "G%iS" % (j)
				else:
					sc_j = atomind[j]['CA']
					sca_j = "G%i" % (j)
				nbfix[(sca_i,sca_j)] = [nncmap[(i,j)][0],nncmap[(i,j)][0]*nat_scale_face]
				gomod_nbfix[(sca_i,sca_j)] = [nncmap[(i,j)][0],nncmap[(i,j)][0]*nat_scale_face]
	print "# side-chain contacts = ", nsc
	print "# hydrogen bonds = ", nhb
	#
	return nbfix,gomod_nbfix,gomod_golist,qlist,nat_scale_fac,allqlist

def mergemap(mapa,mapb):
	merged = {}
	for k in mapa.keys():
		#print k
		merged[k] = mapa[k]
	#print "done with a"
	for k in mapb.keys():
		#print k
		merged[k] = mapb[k]
	#print "done with b"
	return merged

def gotozero(go_map,nongo_map,nres):
	zerogo = {}
	go_keys = go_map.keys()
	nongo_keys = nongo_map.keys()
	for k in nongo_keys:
		zerogo[k] = nongo_map[k]
	for k in go_keys:
		zerogo[k] = (go_map[k][0],0.0)
	return zerogo


def make_sc_sigrep(crd_dat,bbrepdist=std_bb_sigrep,screpdist=std_sc_sigrep):
	nblist = []
	residues = crd_dat.keys()
	residues.sort()
	for r in residues:
		nblist.append(("G%i"%r,bbrepdist))
		if 'SCCENT' in crd_dat[r].keys():
			nblist.append(("G%iS"%r,screpdist))
	return nblist

# OUTPUT FUNCTIONS:
#----------------------------------------------------------------------

def write_CA_pdb2(filename, xyz, crd_dat, chains):
	"""write CA-only coordinates, possibly with multiple chains"""
	outp = open(filename,"w")
	at = 0
	chains.sort()
	if "SOLV" in chains:	# water at end
		chains = filter(lambda x: x!="SOLV",chains) + [ "SOLV" ]
	for chain in chains:
		pdbc = xyz[chain]
		for p in pdbc:
			res = p[0]
			type = crd_dat[chain][res]['name']
			at+=1
			if chain == "SOLV":
				# this distinction for the benefit of CHARMM
				if res < 1000:
					outp.write("ATOM  %5i  W4  SOL  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
							% (at,res,p[1],p[2],p[3],1.00,0.00,chain))
				else:
					outp.write("ATOM  %5i  W4  SOL   %4i   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
						% (at,res,p[1],p[2],p[3],1.00,0.00,chain))
			else:
				if res < 1000:
					outp.write("ATOM  %5i  CA  %3s  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
						% (at,type,res,p[1],p[2],p[3],1.00,0.00,chain))
				else:
					outp.write("ATOM  %5i  CA  %3s   %4i   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
						% (at,type,res,p[1],p[2],p[3],1.00,0.00,chain))
	outp.write("END\n")
	outp.close()


def make_crowd_map(siglist,ncrowd,rcrowd,ecrowd=ecrowd_std):
	crowd_map= {}
	crowd_siglist = {}
	skeys = siglist.keys()
	skeys.sort()
	nprot = len(skeys)
	for s in skeys:
		crowd_siglist[s] = siglist[s]
		for c in range(1,ncrowd+1):
			crowd_map[(s,nprot+c)] = [rcrowd+siglist[s]/2.,ecrowd]
	for c in range(1,ncrowd+1):
		for d in range(c+1,ncrowd+1):
			crowd_map[(nprot+c,nprot+d)] =[2*rcrowd,ecrowd]
	for c in range(1,ncrowd+1):
		crowd_siglist[nprot+c] = rcrowd
	return crowd_map, crowd_siglist

def make_ion_map(siglist,resi,resj,ncrowd,rcrowd,ecrowd=ecrowd_std):
	crowd_map= {}
	crowd_siglist = {}
	skeys = siglist.keys()
	skeys.sort()
	nprot = len(skeys)
	for s in skeys:
		crowd_siglist[s] = siglist[s]
		for c in range(1,ncrowd+1):
			crowd_map[(s,nprot+c)] = [rcrowd+siglist[s]/2.,ecrowd]
	for c in range(1,ncrowd+1):
		for d in range(c+1,ncrowd+1):
			crowd_map[(nprot+c,nprot+d)] =[2*rcrowd,ecrowd]
	for c in range(1,ncrowd+1):
		crowd_siglist[nprot+c] = rcrowd
	return crowd_map, crowd_siglist

def make_crowd_map_protein(siglist,ncrowdres,ncrowd,rcrowd,ecrowd=ecrowd_std):
	crowd_map= {}
	skeys = siglist.keys()
	skeys.sort()
	nres = len(skeys)
	nprotres = nres - ncrowdres
	lcrowd = ncrowdres/ncrowd
	# protein-crowder
	for p in range(1,nprotres+1):
		for q in range(nprotres+1,nres+1):
			crowd_map[(p,q)] = [2*rcrowd,ecrowd]
	# crowder-crowder
	for i in range(ncrowd):
		ii = nprotres + i*lcrowd + 1
		for j in range(i+1,ncrowd):
			jj = nprotres + j*lcrowd + 1
			crowd_map[(ii,jj)] = [2*rcrowd,ecrowd]
	return crowd_map


def write_CA_pdb(filename, pdbdat, chains):
	"""write CA-only coordinates, possibly with multiple chains
	or C1' coords for nucleic acids"""
	outp = open(filename,"w")
	at = 0
	chains.sort()
	if "SOLV" in chains:	# water at end
		chains = filter(lambda x: x!="SOLV",chains) + [ "SOLV" ]
	if "CROW" in chains:	# crowders at end
		chains = filter(lambda x: x!="CROW",chains) + [ "CROW" ]
	for chain in chains:
		pdbc = pdbdat[chain]
		for p in pdbc:
			res = p[0]
			at+=1
			if chain == "SOLV":
				# this distinction for the benefit of CHARMM
				if res < 1000:
					outp.write("ATOM  %5i  W4  SOL  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
							% (at,res,p[1],p[2],p[3],1.00,0.00,chain))
				else:
					outp.write("ATOM  %5i  W4  SOL   %4i   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
						% (at,res,p[1],p[2],p[3],1.00,0.00,chain))
			elif chain == "CROW":
				# this distinction for the benefit of CHARMM
				if res < 1000:
					outp.write("ATOM  %5i  CA  CRO  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
							% (at,res,p[1],p[2],p[3],1.00,0.00,chain))
				else:
					outp.write("ATOM  %5i  CA  CRO   %4i   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
						% (at,res,p[1],p[2],p[3],1.00,0.00,chain))
			else:
				if res < 1000:
					outp.write("ATOM  %5i  CA  ALA  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
						% (at,res,p[1],p[2],p[3],1.00,0.00,chain))
				else:
					outp.write("ATOM  %5i  CA  ALA   %4i   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
						% (at,res,p[1],p[2],p[3],1.00,0.00,chain))
	outp.write("END\n")
	outp.close()

def write_CA_gro(filename, pdbdat, chains, boxlen=10.0):
	"""write CA-only coordinates, possibly with multiple chains"""
	outp = open(filename,"w")
	at = 0
	# title
	outp.write("Go model coordinates\n")
	nres = 0
	for chain in chains:
		nres += len(pdbdat[chain])
	outp.write("%5i\n" % (nres))
	cumres = 0
	for chain in chains:
		pdbc = pdbdat[chain]
		for p in pdbc:
			cumres += 1
			res = p[0]
			resname = "G%i" % cumres
			atname = "G%i" % cumres
			at+=1
			outp.write("%5i%5s%5s%5i%8.3f%8.3f%8.3f\n" \
					% (at,resname,atname,res,p[1]/10.0,p[2]/10.0,p[3]/10.0))
			# /10 since .gro uses nm not A!!
	outp.write("%10.5f%10.5f%10.5f\n" % (boxlen, boxlen, boxlen))
	outp.close()

def scavg(list):
	xav,yav,zav=0.0,0.0,0.0
	for l in list:
		xav += l[0]
		yav += l[1]
		zav += l[2]
	flen = float(len(list))
	return xav/flen, yav/flen, zav/flen

def write_CA_SC_pdb(filename, pdbdat, pdetails, chains, regularize):
	"""write CA and SC centroid coordinates, possibly with multiple chains"""
	if filename == "stdout":
		outp = sys.stdout
	else:
		outp = open(filename,"w")
	at = 0
	conect_list = []
	chains.sort()
	if "SOLV" in chains:	# water at end
		chains = filter(lambda x: x!="SOLV",chains) + [ "SOLV" ]
	for chain in chains:
		pdbc = pdbdat[chain]
		nres = len(pdbc)
		prevca = 0
		for r in range(1,nres+1):
			p = pdbc[r-1]
			#resname = pdetails[chain][r]["name"]
			writesc = 0
			bbx, bby, bbz = p[1],p[2],p[3]
			if chain == "SOLV":
				# this distinction for the benefit of CHARMM
				at+=1
				if r < 1000:
					outp.write("ATOM  %5i  W4  SOL  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
							% (at,r,bbx,bby,bbz,1.00,0.00,chain))
				else:
					outp.write("ATOM  %5i  W4  SOL   %4i   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
						% (at,r,bbx,bby,bbz,1.00,0.00,chain))
			else:
				if r < 1000:
					at+=1
					outp.write("ATOM  %5i  CA  %3s  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
						% (at,'ALA',r,bbx,bby,bbz,1.00,0.00,chain))
						#% (at,resname,r,bbx,bby,bbz,1.00,0.00,chain))
					if prevca > 0:
						conect_list.append((prevca,at))
					prevca = at
					if "SCCENT" in pdetails[chain][r].keys():
						at+=1
						scx,scy,scz = pdetails[chain][r]["SCCENT"]
						outp.write("ATOM  %5i  SC  %3s  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
								% (at,'ALA',r,scx,scy,scz,1.00,0.00,chain))
								#% (at,resname,r,scx,scy,scz,1.00,0.00,chain))
						conect_list.append((at,at-1))
				else:
					at+=1
					outp.write("ATOM  %5i  CA  %3s   %4i   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
						% (at,'ALA',r,bbx,bby,bbz,1.00,0.00,chain))
						#% (at,resname,r,bbx,bby,bbz,1.00,0.00,chain))
					if prevca > 0:
						conect_list.append((prevca,at))
					prevca = at
					if "SCCENT" in pdetails[chain][r].keys():
						at+=1
						scx,scy,scz = pdetails[chain][r]["SCCENT"]
						outp.write("ATOM  %5i  SC  %3s  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
								% (at,'ALA',r,scx,scy,scz,1.00,0.00,chain))
								#% (at,resname,r,scx,scy,scz,1.00,0.00,chain))
						conect_list.append((at,at-1))
	for x in conect_list:
		p,q=x
		if p>q:
			t=q
			q=p
			p=t
		outp.write("%6s%5i%5i\n"%("CONECT",p,q))
	outp.write("END\n")
	outp.close()

def calc_dir(a,b):
	dx = b[0]-a[0]
	dy = b[1]-a[1]
	dz = b[2]-a[2]
	dr = math.sqrt(dx**2+dy**2+dz**2)
	return (dx/dr,dy/dr,dz/dr)

def sc_mult(bondlen,director):
	return (bondlen*director[0],bondlen*director[1],bondlen*director[2])

def rotate(angle,vector):
	x = vector[0]
	y = vector[1]
	cosa = math.cos(angle)
	sina = math.sin(angle)
	nu_x = cosa*x+sina*y
	nu_y = -sina*x+cosa*y
	return (nu_x,nu_y,0.0)

def addv(a,b):
	return (a[0]+b[0],a[1]+b[1],a[2]+b[2])

def make_CA_linear(bonds, angles, n_repeat=1):
	"""write CA-only coordinates, possibly with multiple chains"""
	# first atom goes at origin
	# second atom on x-axis
	coords = [(1,0.0,0.0,0.0),(2,bonds[0][-1],0.0,0.0)]
	# remaining atoms all in x,y plane
	director = calc_dir(coords[0][1:],coords[1][1:])
	sign = 1.0
	idx = 3
	for bond in range(1,len(bonds)):
		bondlen = bonds[bond][-1]
		tmpv = sc_mult(bondlen,director)
		tmpv = rotate(sign*(math.pi-angles[bond-1][-1][0]*math.pi/180.0),tmpv)
		coords.append(tuple([idx]) + addv(coords[-1][1:],tmpv))
		idx+=1
		director = calc_dir(coords[-2][1:],coords[-1][1:])
		sign *= -1.0

	for rep in range(1,n_repeat):
		for bond in range(0,len(bonds)-1):
			bondlen = bonds[bond][-1]
			tmpv = sc_mult(bondlen,director)
			tmpv = rotate(sign*(math.pi-angles[bond-1][-1][0]*math.pi/180.0),tmpv)
			coords.append(tuple([idx]) + addv(coords[-1][1:],tmpv))
			idx+=1
			director = calc_dir(coords[-2][1:],coords[-1][1:])
			sign *= -1.0
		if rep < n_repeat-1: # not last repeat, then link to next
			bond = len(bonds)-1
			bondlen = bonds[bond][-1]
			tmpv = sc_mult(bondlen,director)
			tmpv = rotate(sign*(math.pi-angles[bond-1][-1][0]*math.pi/180.0),tmpv)
			coords.append(tuple([idx]) + addv(coords[-1][1:],tmpv))
			idx+=1
			director = calc_dir(coords[-2][1:],coords[-1][1:])
			sign *= -1.0
	return coords

# deprecated, but retained for compatibility with older scripts
def write_CA_pdb_linear(filename, bonds, angles, chain):
	"""write CA-only coordinates, possibly with multiple chains"""
	xyz = { chain: make_CA_linear(bonds, angles) }
	write_CA_pdb(filename, xyz, [ chain ])


def make_CA_nclink(pdbdat, nrep, chain="NULL"):
	"""create pseudo-native coordinates for n-c linked polyproteins 
	(assumed only one chain)"""
	at = 0
	if chain == "NULL":
		chain = pdbdat.keys()[0]
	pdbc = pdbdat[chain]
	nres = len(pdbc)
	x_nc = pdbc[-1][1] - pdbc[0][1]
	y_nc = pdbc[-1][2] - pdbc[0][2]
	z_nc = pdbc[-1][3] - pdbc[0][3]
	r_nc = math.sqrt(x_nc**2 + y_nc**2 + z_nc**2)
	dx = x_nc/r_nc
	dy = y_nc/r_nc
	dz = z_nc/r_nc
	DX = x_nc + dx * std_bond
	DY = y_nc + dy * std_bond
	DZ = z_nc + dz * std_bond
	coords = []
	for r in range(nrep):
		xoff = r*DX
		yoff = r*DY
		zoff = r*DZ
		for p in pdbc:
			res = p[0] + nres*r
			coords.append((res, p[1]+xoff,p[2]+yoff,p[3]+zoff))
	return coords

# deprecated, but retained for compatibility
def write_CA_pdb_nclink(filename, pdbdat, nrep):
	"""create pseudo-native coordinates for n-c linked polyproteins 
	(assumed only one chain)"""
	chain = pdbdat.keys()[0]
	xyz = { chain: make_CA_nclink(pdbdat,nrep,chain) }
	write_CA_pdb(filename, xyz, [ chain ])

# should be deprecated, but not yet!
def write_CA_pdb_hetero(filename, pdbdat_a, pdbdat_b,na):
	"""create pseudo-native coordinates for n-c linked polyproteins 
	(assumed only one chain)"""
	outp = open(filename,"w")
	at = 0
	pdb_a = pdbdat_a[pdbdat_a.keys()[0]]
	pdb_b = pdbdat_b[pdbdat_b.keys()[0]]
	nres_a = len(pdb_a)
	nres_b = len(pdb_b)
	x_nc_a = pdb_a[-1][1] - pdb_a[0][1]
	y_nc_a = pdb_a[-1][2] - pdb_a[0][2]
	z_nc_a = pdb_a[-1][3] - pdb_a[0][3]
	x_nc_b = pdb_b[-1][1] - pdb_b[0][1]
	y_nc_b = pdb_b[-1][2] - pdb_b[0][2]
	z_nc_b = pdb_b[-1][3] - pdb_b[0][3]
	r_nc_a = math.sqrt(x_nc_a**2 + y_nc_a**2 + z_nc_a**2)
	r_nc_b = math.sqrt(x_nc_b**2 + y_nc_b**2 + z_nc_b**2)
	dx_a = x_nc_a/r_nc_a
	dy_a = y_nc_a/r_nc_a
	dz_a = z_nc_a/r_nc_a
	dx_b = x_nc_b/r_nc_b
	dy_b = y_nc_b/r_nc_b
	dz_b = z_nc_b/r_nc_b
	DX_a = x_nc_a + dx_a * std_bond
	DY_a = y_nc_a + dy_a * std_bond
	DZ_a = z_nc_a + dz_a * std_bond
	DX_b = x_nc_b + dx_b * std_bond
	DY_b = y_nc_b + dy_b * std_bond
	DZ_b = z_nc_b + dz_b * std_bond
	oa_x = pdb_a[0][1]
	oa_y = pdb_a[0][2]
	oa_z = pdb_a[0][3]
	ob_x = pdb_b[0][1]
	ob_y = pdb_b[0][2]
	ob_z = pdb_b[0][3]
	for r in range(na):
		xoff = r*DX_a
		yoff = r*DY_a
		zoff = r*DZ_a
		for p in pdb_a:
			res = p[0] + nres_a*r
			outp.write("ATOM  %5i  CA  ALA  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n" \
					% (res,res,p[1]+xoff-oa_x,p[2]+yoff-oa_y,\
						p[3]+zoff-oa_z,1.00,0.00,""))
	# insert
	xoff = na*DX_a
	yoff = na*DY_a
	zoff = na*DZ_a
	for p in pdb_b:
		res = p[0] + nres_a*na
		outp.write("ATOM  %5i  CA  ALA  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n" \
				% (res,res,p[1]+xoff-ob_x,p[2]+yoff-ob_y,p[3]+zoff-ob_z,\
				1.00,0.00,""))
	outp.write("END\n")
	for r in range(na):
		xoff = (na+r)*DX_a+DX_b
		yoff = (na+r)*DY_a+DY_b
		zoff = (na+r)*DZ_a+DZ_b
		for p in pdb_a:
			res = p[0] + nres_a*(na+r) + nres_b
			outp.write("ATOM  %5i  CA  ALA  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n" \
					% (res,res,p[1]+xoff-oa_x,p[2]+yoff-oa_y,\
					p[3]+zoff-oa_z,1.00,0.00,""))
	outp.close()

def write_sequences(filename,chains,disulf_str,key="a protein",nsolv=0,ncrowd=0,n_repeat=1):
	"""write sequences for streaming into CHARMM
	The argument 'chains' is a map from chain name -> list containing
	the actual residue names"""
	oseq = open(filename,"w")
	oseq.write("* This CHARMM .seq file describes a Go model of %s\n*\n\n"%(key))
	chain_keys = chains.keys()
	chain_keys.sort()
	for chain in chain_keys:
		seq = chains[chain]
		nres = len(seq)*n_repeat
		oseq.write("read sequence card\n")
		oseq.write("* Sequence for Go model of %s\n*\n%i\n" % (key,nres))
		cres = 0
		for repeat in range(n_repeat):
			for r in seq:
				oseq.write( "%-6s " % (r))
				cres += 1
				if cres % 8 == 0:
					oseq.write("\n")
		oseq.write("\n\n")

		oseq.write("generate %s setup\nrename resname ala sele all end\n\n"%(chain))
		if len(disulf_str)>0:	# need to do per-chain disulfides at some point...
			ssbonds = disulf_str.split(':')
			for ss in ssbonds:
				ii,jj = ss.split('-')
				i,j=int(ii),int(jj)
		 		oseq.write("patch disu %s %i %s %i setup\n"%(chain,i,chain,j))
	if ncrowd > 0:
		oseq.write("read sequence card\n")
		oseq.write("* Sequence for %i crowders\n*\n%i\n" % (ncrowd,ncrowd))
		cres = 0
		for r in range(ncrowd):
			oseq.write( "C%i " % (r+1))
			cres += 1
			if cres % 8 == 0:
				oseq.write("\n")
		oseq.write("\n\n")
		oseq.write("generate CROW setup noangle nodihe\n\n")
		oseq.write("rename resname CRO sele segi CROW end\n\n")


	if nsolv > 0:
		#oseq.write("read sequence SOL %i\n"%(nsolv))
		cres
		oseq.write("read sequence card\n")
		oseq.write("* Sequence for %i waters\n*\n%i\n" % (nsolv,nsolv))
		cres = 0
		for r in range(nsolv):
			oseq.write( "%-6s " % ("SOL"))
			cres += 1
			if cres % 8 == 0:
				oseq.write("\n")
		oseq.write("\n\n")
		oseq.write("generate SOLV setup noangle nodihe\n\n")

	oseq.close()

def write_topology(filename,residues,addcharges,key="a protein"):
	"""write topology for streaming into CHARMM
	residues is a list of residues and their masses to write"""
	otop = open(filename,"w")
	otop.write("* This CHARMM .top file describes a Go model of %s\n"%(key))
	otop.write("*\n")
	otop.write("\n")
	otop.write("read rtf card\n")
	otop.write("* Topology for Go model of %s\n"%(key))
	otop.write("*\n")
	otop.write("   20   1\n")

	nres = len(residues)
	for i in range(nres):
		r = residues[i][0]
		mass = massmap[residues[i][1].strip()]
		otop.write("MASS %-5i %4s  %12.6f\n" % (i+1,r,mass))
	# ... for water
	otop.write("MASS %-5i %4s  %12.6f\n" % (nres+1,"W4",72.0))

	otop.write("\nDECL +CA\n\nAUTOGENERATE ANGLES DIHEDRAL\n\n")

	for res in residues:
		if addcharges:
			charge = res2charge[res[1]]
		else:
			charge = 0.0
		otop.write("RESI %4s    %5.3f\n" % (res[0],charge))
		otop.write("GROU\n")
		otop.write("Atom CA %4s   %5.3f\n" % (res[0],charge))
		#otop.write("Bond CA +CA\n\n")
		if res[0][0] != 'C':
			otop.write("Bond CA +CA\n")
		otop.write("\n")
	# ... for water
	otop.write("RESI SOL     0.0\n")
	otop.write("GROU\n")
	otop.write("Atom W4 W4   0.0\n\n")
	otop.write("PRES DISU    0.0\n")
	otop.write("BOND 1CA 2CA\n\n")

	otop.write("END\n\n")
	otop.close()

def write_sc_topology(filename,residues,key="a protein"):
	"""write topology for streaming into CHARMM - WITH SIDE-CHAINS
	residues is a list of residues and their masses to write"""
	otop = open(filename,"w")
	otop.write("* This CHARMM .top file describes a Go model of %s\n"%(key))
	otop.write("* ...with side-chains\n")
	otop.write("*\n")
	otop.write("\n")
	otop.write("read rtf card\n")
	otop.write("* Topology for Go model of %s\n"%(key))
	otop.write("*\n")
	otop.write("   20   1\n")

	nres = len(residues)
	idx = 1
	for r in residues:
		res = r[0]
		type = r[1].strip()
		if type in ['A','T','C','G','U']:
				bbmass = 195.0
				scmass = scmassmap[type]
				otop.write("MASS %-5i %4s  %12.6f\n" % (idx,res,bbmass))
				idx+=1
				otop.write("MASS %-5i %4s  %12.6f\n" % (idx,res+'S',scmass))
				idx+=1
		else:
			if type in scmassmap.keys():	# has side-chain
				bbmass = 56.0
				scmass = scmassmap[type]
				otop.write("MASS %-5i %4s  %12.6f\n" % (idx,res,bbmass))
				idx+=1
				otop.write("MASS %-5i %4s  %12.6f\n" % (idx,res+'S',scmass))
				idx+=1
			else:
				mass = massmap[type.strip()]
				otop.write("MASS %-5i %4s  %12.6f\n" % (idx,res,mass))
				idx+=1
	# ... for water
	otop.write("MASS %-5i %4s  %12.6f\n" % (idx,"W4",72.0))
	idx+=1

	otop.write("\nDECL +CA\n\nAUTOGENERATE ANGLES DIHEDRAL\n\n")

	for r in residues:
		res, type = r[0], r[1].strip()
		otop.write("RESI %4s     0.0\n" % (res))
		otop.write("GROU\n")
		otop.write("Atom CA %4s   0.0\n" % (res))
		if type in scmassmap.keys():
			otop.write("Atom SC %4s   0.0\n" % (res+'S'))
			otop.write("Bond CA +CA   CA SC\n\n")
		else:
			otop.write("Bond CA +CA\n\n")
	# ... for water
	otop.write("RESI SOL     0.0\n")
	otop.write("GROU\n")
	otop.write("Atom W4 W4   0.0\n\n")

	otop.write("END\n\n")
	otop.close()

def write_gro_top(filename,key,residues, bonds, angles, dihedrals,
		sig_rep,nbfix,nsfac,tfbb,seqs,symm,winte,potential,khmap,addcharges,
		intrakh,n_repeat=1):
	"""write GROMACS topology file including
	everything !!! (parameters,etc.)"""
	if len(khmap.keys()) > 0:
		do_kh = 1
	else:
		do_kh = 0
	nres = len(residues)
	chains = seqs.keys()
	eps_res = big_fat_fudge_factor*tfbb
	#Note epsilon must be positive for GROMACS!
	eij_rep = +1.5e-3 * eps_res / 21.5 * 4.184	# fudge
	kbond = eps_res * 200.0
	kangle = eps_res * 40.0
	dihe_fac = eps_res * 0.4
	gmx_kbond = 2.0*kbond*4.184*100.0	# convert from kcal/(mol.A^2) to kJ/(mol.nm^2)
	                                # factor of two because gmx E(dx)=0.5kbond*dx**2
	                                #                    charmm E(dx)=kbond*dx**2
	gmx_kangle = 2.0*kangle*4.184	# convert from kcal/(mol.rad^2) to kJ/(mol.rad^2)
	                                # factor of two because gmx E(dtheta)=0.5kangle*dtheta**2
	                                #                    charmm E(dtheta)=kangle*dtheta**2
	gmx_dihe_fac = dihe_fac*4.184	# convert from kcal/mol to kJ/mol
	#
	otop = open(filename,"w")
	#============================================================
	# HEADER
	#============================================================
	otop.write("; This is a GROMACS topology file describing a Go model of %s\n" % (key))
	otop.write("; Parameters scaled to obtain a Tf of: %12.6f\n" % (tfbb))
	otop.write("; Scaling all native contacts by a factor of %12.6f\n\n" % (nsfac))
	otop.write("; NB: all lengths in nm !!!\n\n")
	#============================================================
	# PARAMETER LEVEL
	#============================================================
	# NB DEFAULTS
	otop.write("[ defaults ]\n")
	otop.write("; nbfunc = 2 ... implies Buckingham, which is hacked to 12-10-6\n")
	otop.write("; nbfunc	comb-rule	gen-pairs	fudgeLJ	fudgeQQ\n")
	otop.write("2		3		no		1.0	1.0\n\n")
	# ATOM TYPES + DEFAULT (REPULSIVE) NON-BONDED INTERACTIONS
	otop.write("[ atomtypes ]\n")
	otop.write("; name  bond_type    mass    charge   ptype          sigma      epsilon\n")
	otop.write("; !! note well: sigma in nm !!!\n")
	nres = len(sig_rep.keys())
	massm = {}
	for i in range(nres):
		res = residues[i][0]
		if addcharges:
			charge = res2charge[residues[i][1]]
		else:
			charge = 0.0
		mass = massmap[residues[i][1].strip()]
		massm[res] = mass	# save for later
		#sigma = sig_rep[i+1]/(10.0)
                sigma_fei=0.9
		otop.write("%4s   %4s   %8.3f   %8.3f    A    %12.5e   %12.5e %12.5e\n" \
                           #% (res,res,mass,0.0,0.4,eij_rep,0.0)) first test simulation 27/Apr
                % (res,res, mass, 0.0, 0.001*sigma_fei**12, 0.001*sigma_fei**10, 0.001*sigma_fei**6)) #re write combination rule paras, 0.001 is epslion
		#otop.write("%4s   %4s   %8.3f   %8.3f    A    0.0 0.0 0.0\n" \
		#		% (res,res,mass,charge))
	if len(winte.keys()) > 0:	# water
		#otop.write("%4s   %4s   %8.3f   %8.3f    A    %12.5e   %12.5e\n" \
		#		% ("W4","W4",72.0,0.0,0.528,eij_rep))
		otop.write("%4s   %4s   %8.3f   %8.3f    A    0.0 0.0 0.0\n" \
				% ("W4","W4",72.0,0.0))
	# EQUIV OF NBFIX
	otop.write("\n\n[ nonbond_params ]\n")
	otop.write(";%-4s %-5s func %-12s %-12s %-12s\n" \
			% ('i','j','13*eps*R^12', '18*eps*R^10', '4*eps*R^6'))


	for i in range(1,nres+1):
		for c in chains:
			if residues[i-1][0] in seqs[c]:
				chaini = c
		for j in range(i,nres+1):
			for c in chains:
				if residues[j-1][0] in seqs[c]:
					chainj = c
			if (i,j) in nbfix.keys():
				eij = abs(nbfix[(i,j)][1])	# positive eij!!
				rij = nbfix[(i,j)][0] / (10.0) 	# A to nm
				#eij *= nsfac*4.184
				eij *= 4.184 # nsfac should be incorporated already
				A = 13.0*eij*rij**12
				B = 18.0*eij*rij**10
				C = 4.0*eij*rij**6
                                nerf_res58_fei = 1.0 #only for titin, design exp from Jane
                                if i == 58 or j == 58:
                                        A *= nerf_res58_fei
                                        B *= nerf_res58_fei
                                        C *= nerf_res58_fei
                                print i,j,rij,A,B,C, eij
				otop.write("G%-5i G%-5i 2  %14.6e %14.6e %14.6e\n" \
					% (i,j,A,B,C))
			elif do_kh == 1 and not \
					( not intrakh and (chaini==chainj) ):
				eij = khmap[(i,j)][1]*4.184	# kcal to kJ
				rij = khmap[(i,j)][0]/10.	# A to nm
				A = -1./rij**6
				B = 4.0*eij*rij**12
				C = 4.0*eij*rij**6
				otop.write("G%-5i G%-5i 2  %14.6e %14.6e %14.6e\n" \
					% (i,j,A,B,C))
			else:
				rij = (sig_rep[i]+sig_rep[j])/(2.0*10.0)
				A = 13.0*eij_rep*rij**12
				B = 18.0*eij_rep*rij**10
				C = 4.0*eij_rep*rij**6
                                #print "##",i,j,rij,A,B,C, eij_rep
				otop.write("G%-5i G%-5i 2  %14.6e %14.6e %14.6e\n" \
					% (i,j,A,B,C))
	# this is for 1-4 interations
	# seems to cause problems for GMX 3.3
	#otop.write("\n\n[ pairtypes ] 	; 1-4 interactions\n")
	#if potential == "12_10_6":
	#	otop.write(";%-4s %-5s func %-12s %-12s %-12s\n" \
	#			% ('i','j','13*eps*R^12', '18*eps*R^10', '4*eps*R^6'))
	#elif potential == "12_10":
	#	otop.write(";%-4s %-5s func %-12s %-12s\n" \
	#			% ('i','j','6*eps*R^10', '5*eps*R^12'))
	#else:	# 12-6
	#	otop.write(";%-4s %-5s func %-12s %-12s\n" \
	#			% ('i','j','2*eps*R^6', 'eps*R^12'))
	#for i in range(1,nres+1):
	#	for j in range(i,nres+1):
	#		if (i,j) in nbfix.keys():
	#			eij = abs(nbfix[(i,j)][1])	# positive eij!!
	#			rij = nbfix[(i,j)][0] / (10.0) 	# A to nm
	#			eij *= nsfac*4.184
	#			A = 13.0*eij*rij**12
	#			B = 18.0*eij*rij**10
	#			C = 4.0*eij*rij**6
	#			otop.write("G%-5i G%-5i 2  %12.6e %12.6e %12.6e\n" \
	#				% (i,j,A,B,C))
	#		else:
	#			rij = (sig_rep[i]+sig_rep[j])/(2.0*10.0)
	#			A = 13.0*eij_rep*rij**12
	#			B = 18.0*eij_rep*rij**10
	#			C = 4.0*eij_rep*rij**6
	#			otop.write("G%-5i G%-5i 2  %12.6e %12.6e %12.6e\n" \
	#				% (i,j,A,B,C))

	## BONDS
	#otop.write("\n\n[ bondtypes ]\n")
	#otop.write("; i    j  func       b0          kb\n")
	#for bond in bonds:
	#	r1 = "G%i" % (bond[0])
	#	r2 = "G%i" % (bond[1])
	#	R0 = bond[2]/10.0	# A to nm
	#	otop.write("%-6s %-6s   1   %12.6f %12.6f\n" % (r1,r2,R0,gmx_kbond))

	## ANGLES
	#otop.write("\n\n[ angletypes ]\n")
	#otop.write(";  i    j    k  func       th0       cth\n")
	#for angle in angles:
	#	r1 = "G%i" % (angle[0])
	#	r2 = "G%i" % (angle[1])
	#	r3 = "G%i" % (angle[2])
	#	theta0 = angle[3][0]
	#	otop.write("%-6s %-6s %-6s   1   %12.6f %12.6f\n" % (r1,r2,r3,theta0,gmx_kangle))

	# DIHEDRALS
	#otop.write("\n\n[ dihedraltypes ]\n")
	#otop.write(";%-5s %-6s %-6s %-6s %2s %-12s %-12s %-2s\n" % \
	#		('i','j','k','l','fn','phi0','kphi','n'))
	#for dihe in dihedrals:
	#	r1 = "G%i" % (dihe[0])
	#	r2 = "G%i" % (dihe[1])
	#	r3 = "G%i" % (dihe[2])
	#	r4 = "G%i" % (dihe[3])
	#	for p in dihe[4]:
	#		k = p[0]*gmx_dihe_fac
	#		mult = p[1]
	#		phi0 = p[2]
	#		otop.write("%-6s %-6s %-6s %-6s %2s %12.6f %12.6f %2i\n" % \
	#				(r1,r2,r3,r4,1,phi0,k,mult))

	#============================================================
	# MOLECULE LEVEL
	#============================================================
	# MOLECULES
	#if symm:
	#	otop.write("%4s   %i\n"%(chains[0],len(chains)))
	#else:
	for chain in chains:
		otop.write("\n\n[ moleculetype ]\n")
		otop.write("; molname 	nrexcl\n")
		otop.write("; nrexcl = 2 means all atoms separated by 2 bonds or less excluded\n")
		if chain == " ":
			otop.write("%4s   2\n"%("A"))
		else:
			otop.write("%4s   2\n"%(chain))
		# ATOMS
		otop.write("\n\n[ atoms ]\n")
		otop.write("; id	at type	res nr 	residu name	at name		cg nr	charge	mass\n")
		idx = 0
		seq = seqs[chain]
		seq_len= len(seq)
		for repeat in range(n_repeat):
			seq_idx = 0
			for s in seq:
				idx+=1
				seq_idx+=1
				mass = massm[s]
				if addcharges:
					for r in residues:
						if r[0] == s:
							charge = res2charge[r[1]]
				else:
					charge = 0.0
				otop.write("%5i %4s    %5i       %4s         %4s        %5i    %8.3f   %8.3f\n" \
					%(idx,s,idx,"ALA",s,idx,charge,mass))
		# BONDS
		otop.write("\n\n[ bonds ]\n")
		otop.write("; i    j  func       b0          kb\n")
		#print bonds[chain]
		for repeat in range(n_repeat):
			offset = repeat*seq_len
			for bond in bonds[chain]:
				r1, r2 = bond[0], bond[1]
				r1+=offset
				r2+=offset
				if r1>r2:
					if repeat != n_repeat-1:
						r2+=seq_len
					else:
						continue
				R0 = bond[2]/10.0	# A to nm
				otop.write("%6i %6i   1   %12.6f %12.6f\n" % (r1,r2,R0,gmx_kbond))

		otop.write("\n\n[ angles ]\n")
		otop.write(";  i    j    k  func       th0       cth\n")
		for repeat in range(n_repeat):
			offset = repeat*seq_len
			for angle in angles[chain]:
				r1, r2, r3 = angle[0],angle[1],angle[2]
				r1+=offset
				r2+=offset
				r3+=offset
				if r1>r2:
					if repeat != n_repeat-1:
						r2+=seq_len
					else:
						continue
				if r1>r3:
					if repeat != n_repeat-1:
						r3+=seq_len
					else:
						continue
				theta0 = angle[3][0]
				otop.write("%6i %6i %6i   1   %12.6f %12.6f\n" % (r1,r2,r3,theta0,gmx_kangle))

		otop.write("\n\n[ dihedrals ]\n")
		otop.write(";%-5s %-6s %-6s %-6s %2s %-12s %-12s %-2s\n" % \
				('i','j','k','l','fn','phi0','kphi','n'))
		for repeat in range(n_repeat):
			offset = repeat*seq_len
			for dihe in dihedrals[chain]:
				rr1,rr2,rr3,rr4 = dihe[0],dihe[1],dihe[2],dihe[3]
				r1 = int(rr1[1:])
				r2 = int(rr2[1:])
				r3 = int(rr3[1:])
				r4 = int(rr4[1:])
				r1+=offset
				r2+=offset
				r3+=offset
				r4+=offset
				if r1>r2:
					if repeat != n_repeat-1:
						r2+=seq_len
					else:
						continue
				if r1>r3:
					if repeat != n_repeat-1:
						r3+=seq_len
					else:
						continue
				if r1>r4:
					if repeat != n_repeat-1:
						r4+=seq_len
					else:
						continue
				for p in dihe[4]:
					k = p[0]*gmx_dihe_fac
					mult = p[1]
					phi0 = p[2]
					otop.write("%6i %6i %6i %6i %2i %12.6f %12.6f %2i\n" % \
							(r1,r2,r3,r4,1,phi0,k,mult))

	# PAIRS - GO CONTACTS COME IN HERE!
	#============================================================
	# SYSTEM LEVEL
	#============================================================
	otop.write("\n\n[ system ]\n")
	otop.write("Another crappy Go model\n")
	otop.write("\n\n[ molecules ]\n")
	if symm:
		otop.write("%s       %i\n"%(chains[0],len(chains)))
	else:
		for chain in chains:
			if chain == " ":
				otop.write("%4s      1\n"%("A"))
			else:
				otop.write("%4s      1\n"%(chain))

def write_parm_header(filep,key,Tf,nsfac):
	"""write header of parameter file"""
	filep.write("* This CHARMM .param file describes a Go model of %s\n" % (key))
	filep.write("* Parameters scaled to obtain a Tf of: %12.6f\n" % (Tf))
	filep.write("* Scaling all native contacts by a factor of %12.6f\n\n" % (nsfac))
	filep.write("read param card\n")
	filep.write("* Parameters for Go model of %s\n\n" % (key))

def write_bonded(outp,bonds,angles,dihedrals,eps_res):
	"""write bonded terms only to parameter file"""
	kbond = eps_res * 200.0
	kangle = eps_res * 40.0
	dihe_fac = eps_res * 0.4
	outp.write("BOND\n")
	for bond in bonds:
		r1 = "G%i" % (bond[0])
		r2 = "G%i" % (bond[1])
		R0 = bond[2]
		outp.write("%-6s %-6s %12.6f %12.6f\n" % (r1,r2,kbond,R0))
	outp.write("\n")
	outp.write("ANGLE\n")
	for angle in angles:
		r1 = "G%i" % (angle[0])
		r2 = "G%i" % (angle[1])
		r3 = "G%i" % (angle[2])
		theta0 = angle[3]
		if len(theta0) == 1:	# pure Go case
			outp.write("%-6s %-6s %-6s %12.6f %12.6f\n" % (r1,r2,r3,kangle,theta0[0]))
		elif len(theta0) == 2:	# harmonic statistical potential
			outp.write("%-6s %-6s %-6s %12.6f %12.6f\n" % (r1,r2,r3,theta0[0],theta0[1]))
		elif len(theta0) == 4:	# urey-bradley
			outp.write("%-6s %-6s %-6s %12.6f %12.6f %12.6f %12.6f\n" % \
					(r1,r2,r3,theta0[0],theta0[1],theta0[2],theta0[3]))
		else:			# double harmonic statistical potential
			outp.write("%-6s %-6s %-6s %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n" \
					% (r1,r2,r3,theta0[0],theta0[1],theta0[2],theta0[3],
						theta0[4],theta0[5]))
	outp.write("\n")
	outp.write("DIHEDRAL\n")
	for dihe in dihedrals:
		#r1 = "G%i" % (dihe[0])
		#r2 = "G%i" % (dihe[1])
		#r3 = "G%i" % (dihe[2])
		#r4 = "G%i" % (dihe[3])
		r1 = dihe[0]
		r2 = dihe[1]
		r3 = dihe[2]
		r4 = dihe[3]
		for p in dihe[4]:
			k = p[0]*dihe_fac
			mult = p[1]
			phi0 = p[2]
			outp.write("%-6s %-6s %-6s %-6s %12.6f %2i %12.6f\n" % \
					(r1,r2,r3,r4,k,mult,phi0))

def write_sc_bonded(outp,bonds,angles,dihedrals,eps_res):
	"""write bonded terms only to parameter file"""
	kbond = eps_res * 200.0
	kangle = eps_res * 40.0
	dihe_fac = eps_res * 0.4
	outp.write("BOND\n")
	for bond in bonds:
		r1 = bond[0]
		r2 = bond[1]
		R0 = bond[2]
		outp.write("%-6s %-6s %12.6f %12.6f\n" % (r1,r2,kbond,R0))
	outp.write("\n")
	outp.write("ANGLE\n")
	for angle in angles:
		r1 = angle[0]
		r2 = angle[1]
		r3 = angle[2]
		theta0 = angle[3]
		if len(theta0) == 1:	# pure Go case
			outp.write("%-6s %-6s %-6s %12.6f %12.6f\n" % (r1,r2,r3,kangle,theta0[0]))
		elif len(theta0) == 2:	# harmonic statistical potential
			outp.write("%-6s %-6s %-6s %12.6f %12.6f\n" % (r1,r2,r3,theta0[0],theta0[1]))
		else:			# double harmonic statistical potential
			outp.write("%-6s %-6s %-6s %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n" \
					% (r1,r2,r3,theta0[0],theta0[1],theta0[2],theta0[3],
						theta0[4],theta0[5]))
	outp.write("\n")
	outp.write("DIHEDRAL\n")
	for dihe in dihedrals:
		r1 = dihe[0]
		r2 = dihe[1]
		r3 = dihe[2]
		r4 = dihe[3]
		for p in dihe[4]:
			k = p[0]*dihe_fac
			mult = p[1]
			phi0 = p[2]
			outp.write("%-6s %-6s %-6s %-6s %12.6f %2i %12.6f\n" % \
					(r1,r2,r3,r4,k,mult,phi0))

def write_defnb_list(outp,sig_rep,eps_res,ncrowd=0):
	outp.write("\n")
	outp.write("""NONBONDED NBXMOD 3 ATOM CDIEL SHIFT VATOM VDISTANCE VSWITCH -
  CUTNB 399.0 CTOFNB 398.5 CTONNB 395.5 EPS 1.0 WMIN 1.5\n\n""")
	eij_rep = -1.5e-3 * eps_res / 21.5 	# fudge
	nres = len(sig_rep.keys())
	for i in range(1,nres-ncrowd+1):
		outp.write("G%-5i  0.0  %12.6f %12.6f\n" % (i,eij_rep,sig_rep[i]/2.0))
	for i in range(nres-ncrowd+1,nres+1):
		outp.write("C%-5i  0.0  %12.6f %12.6f\n" % (i-(nres-ncrowd),eij_rep,sig_rep[i]/2.0))
	outp.write("W4  0.0  %12.6f %12.6f\n" % (-1.2,2.64))
	outp.write('\n')

def write_sc_defnb_list(outp,sig_rep,eps_res):
	outp.write("\n")
	outp.write("""NONBONDED NBXMOD 3 ATOM CDIEL SHIFT VATOM VDISTANCE VSWITCH -
  CUTNB 399.0 CTOFNB 398.5 CTONNB 395.5 EPS 1.0 WMIN 1.5\n\n""")
	eij_rep = -1.5e-3 * eps_res / 21.5 	# fudge
	for s in sig_rep:
		outp.write("%4s  0.0  %12.6f %12.6f\n" % (s[0],eij_rep,s[1]))
	outp.write("W4  0.0  %12.6f %12.6f\n" % (-1.2,2.64))
	outp.write('\n')

def write_defnb_list_fixed(outp,nres,sig_rep,eps_res):
	outp.write("\n")
	outp.write("""NONBONDED NBXMOD 3 ATOM CDIEL SHIFT VATOM VDISTANCE VSWITCH -
  CUTNB 399.0 CTOFNB 398.5 CTONNB 395.5 EPS 1.0 WMIN 1.5\n\n""")
	eij_rep = -1.5e-3 * eps_res / 21.5 	# fudge
	for i in range(1,nres+1):
		outp.write("G%-5i  0.0  %12.6f %12.6f\n" % (i,eij_rep,sig_rep/2.0))
	outp.write('\n')

def write_qdetails(file,crds,details,go_map):
	outp = open(file,"w")
	nres = len(crds.keys())
	for i in range(1,nres):
		for j in range(i,nres+1):
			if (i,j) in details.keys():
				outp.write("%3s %5i %3s %5i %12.6f %s\n" % \
						(crds[i]['name'],i,crds[j]['name'],\
						j,go_map[(i,j)][1], \
						details[(i,j)]))

def write_nbfix(outp,nbfix,nsfac,nres,ncrowd=0):
	if len(nbfix.keys())==0:
		return
	outp.write("\nNBFIX\n")
	nbkeys = nbfix.keys()
	nbkeys.sort(keysort)
	for k in nbkeys:
		i,j = k
		if i < nres-ncrowd+1:
			iname = "G%i" % i
		else:
			iname = "C%i" % (i-nres+ncrowd)
		if j < nres-ncrowd+1:
			jname = "G%i" % j
		else:
			jname = "C%i" % (j-nres+ncrowd)
		eij = nbfix[(i,j)][1]
		if abs(eij) < 1e-5:
			outp.write("%5s %5s %12.6e %12.6f\n" \
				% (iname,jname,eij,nbfix[(i,j)][0]))
		else:
			eij *= nsfac
			outp.write("%5s %5s %12.6f %12.6f\n" \
				% (iname,jname,eij,nbfix[(i,j)][0]))

def write_sc_nbfix(outp,nbfix,nsfac):
	outp.write("\nNBFIX\n")
	for k in nbfix.keys():
		i,j = k[0],k[1]
		outp.write("%4s %4s %12.6e %12.6f\n" \
						% (i,j,nbfix[k][1],nbfix[k][0]))

def write_winte(outp,winte):
	for k in winte.keys():
		sig, eps = winte[k]
		outp.write("%4s  %4s  %12.6f %12.6f\n"%(k[0],k[1],eps,sig))

def write_parms(file,key,bonds,angles,dihedrals,sig_rep,nbfix,
		nsfac,tfbb,winte,morse,ncrowd=0):
	"""top level parameter writing function"""
	nres = len(sig_rep.keys())
	outp = open(file,"w")
	write_parm_header(outp,key,tfbb,nsfac)
	eps_res = big_fat_fudge_factor*tfbb
	write_bonded(outp,bonds,angles,dihedrals,eps_res)	# common bit
	#if repdist > 1:
	#	write_defnb_list_fixed(outp,nres,repdist,eps_res)
	#else:
	if morse == 1:
		write_defnb_list(outp,sig_rep,eps_res*3.)
	else:
		write_defnb_list(outp,sig_rep,eps_res,ncrowd)
	write_nbfix(outp,nbfix,1.0,nres,ncrowd)
	if len(winte.keys()) > 0:
		write_winte(outp,winte)
	outp.write("\nEND\n\n")
	outp.close()

def write_sc_parms(file,key,bonds,angles,dihedrals,sig_rep,nbfix,
		nsfac,Tf,winte):
	"""top level parameter writing function - with side-chains"""
	outp = open(file,"w")
	write_parm_header(outp,key,Tf,nsfac)
	eps_res = 0.0054*Tf
	write_sc_bonded(outp,bonds,angles,dihedrals,eps_res)	# common bit
	write_sc_defnb_list(outp,sig_rep,eps_res)
	write_sc_nbfix(outp,nbfix,nsfac)
	if len(winte.keys()) > 0:
		write_winte(outp,winte)
	outp.write("\nEND\n\n")
	outp.close()

def write_qlist(file,qlist,nres,n_repeat):
	outp = open(file,"w")
	for r in range(n_repeat):
		i_off = r*nres
		for s in range(n_repeat):
			j_off = s*nres
			for item in qlist:
				outp.write("%5i %5i %8.3f\n" % (item[0]+i_off,
					item[1]+j_off,item[2]))

def sortq(q1,q2):
	q1i, q1j, q1r = q1
	q2i, q2j, q2r = q2
	if k1i>k2i:
		return 1
	elif k1i<k2i:
		return -1
	else:
		if k1j>k2j:
			return 1
		elif k1j<k2j:
			return -1
		else:
			return 0

def dump_qlist(file,qlist):
	outp = open(file,"w")
	qlist.sort(sortq)
	for q in qlist:
		outp.write("%5i %5i %8.3f\n" % q)
	outp.close()

def write_qlist_intra(file,qlist,nres,n_repeat):
	outp = open(file,"w")
	for r in range(n_repeat):
		off = r*nres
		for item in qlist:
			outp.write("%5i %5i %8.3f\n" % (item[0]+off,
				item[1]+off,item[2]))

def write_qlist_inter(file,qlist,nres,n_repeat):
	outp = open(file,"w")
	for r in range(n_repeat):
		i_off = r*nres
		for s in range(n_repeat):
			if s == r:
				continue
			j_off = s*nres
		for item in qlist:
			outp.write("%5i %5i %8.3f\n" % (item[0]+i_off,
				item[1]+j_off,item[2]))

def write_module_qlist(file,qlist,nres,r):
	outp = open(file,"w")
	off = r*nres
	for item in qlist:
		outp.write("%5i %5i %8.3f\n" % (item[0]+off,
			item[1]+off,item[2]))

def write_golist(file,go_map,nres,n_repeat):
	ogolist = open(file,"w")
	go_keys = go_map.keys()
	go_keys.sort(keysort)
	if n_repeat == 1:	# usual case
		for k in go_keys:
			i,j = k
			rij,eij = go_map[k]
			ogolist.write("%5i %5i %12.6f %12.6f\n" \
					% (i,j,rij,-eij))
	else:
		for k in go_keys:
			i,j = k
			rij,eij = go_map[k]
			for r in range(n_repeat):
				for s in range(n_repeat):
					ogolist.write("%5i %5i %12.6f %12.6f\n" \
							% (r*nres+i,s*nres+j,\
							rij,-eij))

def write_sc_golist(file,go_map,nat_scale_fac,nres,n_repeat):
	ogolist = open(file,"w")
	for k in go_map.keys():
		i,j = k
		eij = go_map[(i,j)][1]
		rij = go_map[(i,j)][0]
		ogolist.write("%5s %5s %12.6f %12.6f\n" % (i,j,	rij,-eij))

# HESSIAN CALCULATION
#======================================================================

def dump_dVdr(file,bonds, go_and_nongo, tfbb):
	"""Write second-derivative matrix for two-body terms only"""
	nres = 0
	eps_res = big_fat_fudge_factor*tfbb
	k_bonded = eps_res*200.0
	for b in bonds:
		nres = max(b[0],b[1],nres)
	H = {}
	for b in bonds:
		i,j,rij = b
		if j<i:
			tmp = i
			i = j
			j = tmp
		for ii in (-2,0,2):
			for jj in (-2,0,2):
				iii = i+ii
				jjj = j+jj
				if iii<1 or jjj<1:
					continue
				if iii>nres or jjj >= nres:
					continue
				if iii == jjj:
					continue
				if abs(iii-jjj) == 1:
					kij = k_bonded
				else:
					kij = 0.3*k_bonded
				k = (iii,jjj)
				if k not in H.keys():
					H[k] = 0.0
				H[k] += kij
	for k in go_and_nongo.keys():
		i,j=k
		sig_ij,eij =go_and_nongo[k]
		if j>i:
			kk = k
		else:
			kk = (j,i)
		if kk not in H.keys():
			H[kk] = 0.0
		#H[kk] += 10. #-216.0*eij/(sig_ij)**2
		H[kk] += -216.0*eij/(sig_ij)**2

	keys = H.keys()
	keys.sort(keysort)
	dvdr_out = open(file,"w")
	for k in keys:
		dvdr_out.write("%5i %5i %12.6e\n" % (k[0],k[1],H[k]))
	dvdr_out.close()



