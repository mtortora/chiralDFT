import os
import sys

import numpy as np


if len(sys.argv) != 6:
	print("Usage is %s pqrFile nRes nRep shift twist" % sys.argv[0])
	sys.exit()


pqrPath = os.path.realpath(sys.argv[1])

nRes = int(sys.argv[2])
nRep = int(sys.argv[3])

shift = float(sys.argv[4])
twist = float(sys.argv[5])

fName = os.path.splitext(pqrPath)[0]

xyzPath = "%s.xyz" % fName

pqrFile = open(pqrPath, mode='r')
xyzFile = open(xyzPath, mode='w')

nAtoms = 0
xyzData = []

offsets = np.arange(nRep) * shift
angles = np.arange(nRep) * twist * np.pi/180.

pqrData = pqrFile.readlines()

for i in range(nRep):
	angle = angles[i]
	offset = offsets[i]

	rot = np.asarray([[np.cos(angle), -np.sin(angle), 0],
					 [np.sin(angle), np.cos(angle), 0],
					 [0,0,1]])
					 
	for pqrLine in pqrData:
		invalidLine = pqrLine.startswith('REMARK') | pqrLine.startswith('TER') | pqrLine.startswith('END')

		if not invalidLine:
			lineData = pqrLine.split()[2:]
			
			atype = lineData[0]
			rtype = lineData[1]
			idRes = int(lineData[2])
			
			if atype in ["C","CA"]:
				atype = "CH1"
					
			elif atype == "CB":
				if rtype == "ILE":
					atype = "CH1"
				else:
					atype = "CH2"
					
			elif atype[:2] == "CD":
				if (rtype == "TRP") & (atype == "CD2"):
					atype = "CR1"
				elif rtype in ["GLU","GLN"]:
					atype = "CH1"
				elif rtype in ["ILE","LEU"]:
					atype = "CH3"
				elif rtype in ["PHE","TRP","TYR"]:
					atype = "CH2r"
				else:
					atype = "CH2"
					
			elif atype[:2] == "CE":
				if (rtype == "TRP") & (atype == "CE2"):
					atype = "CR1"
				elif rtype == "MET":
					atype = "CH3"
				elif rtype in ["PHE","TRP","TYR"]:
					atype = "CH2r"
				else:
					atype = "CH2"
					
			elif atype[:2] == "CG":
				if (rtype in ["ILE","THR"]) & (atype == "CG2"):
					atype = "CH3"
				elif rtype == "VAL":
					atype = "CH3"
				elif rtype in ["ASP","LEU"]:
					atype = "CH1"
				elif rtype in ["PHE","TRP","TYR"]:
					atype = "CR1"
				else:
					atype = "CH2"
			
			elif atype[:2] == "CH":
				atype = "CH2r"
					
			elif atype[:2] == "CZ":
				if rtype == "TYR":
					atype = "CR1"
				elif rtype in ["PHE","TRP"]:
					atype = "CH2r"
				else:
					atype = "CH2"
			
			elif atype == "NE1":
				atype = "NR"
			
			elif atype in ["NZ","ND2","NE2"]:
				atype = "NT"
				
			elif (atype == "N") & (idRes % nRes == 1):
				atype = "NT"
			
			elif atype[:2] in ["OD","OE"]:
				atype = "OM"
				
			elif atype == "OXT":
				atype = "OM"

			elif atype == "OH":
				atype = "OA"
				
			else:
				atype = atype[0]
			
			pos = np.asarray([float(x) for x in lineData[3:6]])
			pos = np.dot(rot, pos) + np.asarray([0,0,offset])
			
			xyzLine = "%s " % atype
			xyzLine +=  " ".join("%.3f" % x for x in pos)
			xyzLine += " %s" % lineData[6]
		
			xyzData.append(xyzLine)
			nAtoms += 1

xyzFile.write("%d\n" % nAtoms)
xyzFile.write("%s\n" % os.path.basename(fName))

for xyzLine in xyzData:
	xyzFile.write("%s\n" % xyzLine)

pqrFile.close()
xyzFile.close()
