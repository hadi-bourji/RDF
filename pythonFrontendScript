import numpy as np
import fort_analysis as fort
import helper as hlp
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import ctypes
from ctypes import CDLL, c_int, c_double, c_float, byref
import time
startTime = time.time()
#
from numpy.ctypeslib import ndpointer
###############
# Specify cpp file
###############
pair_corr_lib = CDLL("./libTest.so")

###############
# Specify trajectory and simulation data
###############
fname = "npt.dcd"
groname = "npt.gro"
ndxname="npt.ndx"
funit = 500
CV = []

###############
# Get molecular data
###############
atomid = hlp.getnameandid_fromndx(ndxname)

###############
# Read dcd file, extract data
###############
tatom, tframe, box = fort.readheaderdcd(fname)

###############
# Specify variables to control calculation
#
# rc2: [Angstrom] Largest shell for rdf
# dr: [Angstrom] Distancce between rdf shells
# cushion: extra shells to avoid if-else condition during binning
# nB: Total number of bins
# H: Array of histogram of type integer
# vol: Volume of simulation box
# ncount: Number of frames analyzed
# N: Number of atom centers of interest
# box2: Box vectors
# N_eq: Number of equilibration frames
# N_prod: Number of production frames
###############
rc2 = 12


dr=0.02
dr_c = ctypes.c_double(dr)

cushion=2000

nB = int((rc2/dr))+cushion 
nB_c = ctypes.c_int(nB)

rc2 = rc2*rc2
rc2_c = ctypes.c_double(rc2)

H = np.zeros(nB, dtype=int) 
H_c = (ctypes.c_int * len(H))(*H)
 
vol = 0
ncount = 0
N = 0
box2 = np.zeros(3,dtype=np.float64)
N_eq = 2750
N_prod = 100

###############
# Open dcd for reading, analyzing
###############
msg = fort.opendcd(funit,fname)
msg = fort.skipheaderdcd(funit)

###############
# Skip equilibration frames
###############
for step in range(N_eq):
    fort.skipframedcd(funit)

###############
# Loop over production frames
###############
#for step in range(min(N_prod,tframe-N_eq)):
for step in range(2750,2800):
    ### Update iteration counter
    ncount += 1 
    
    ### Read frame
    rcoord,box2 = fort.readframedcd(funit,tatom)
    box2_c = np.ctypeslib.as_ctypes(box2)

    ### Store volume cumulatively
    vol = vol + (box2[0]*box2[1]*box2[2])

    ### Extract atom coordinates of interest
    r = []
    for i in atomid['OW']:
        r.append(rcoord[i-1])
    r = np.array(r) 

    N = len(atomid['OW'])
    N_c = ctypes.c_int(N)

    rx = np.array(r[:,0])
    ry = np.array(r[:,1])
    rz = np.array(r[:,2]) 
    
    ### Convert all inputs to ctypes
    rx_c = (ctypes.c_double * len(rx))(*rx)
    ry_c = (ctypes.c_double * len(ry))(*ry)
    rz_c = (ctypes.c_double * len(rz))(*rz)
    pair_corr_lib.My_Function(rx_c, ry_c, rz_c, N_c, box2_c, rc2_c, dr_c, nB_c, H_c)

# Convert back to numpy array for analysis in python   
H = np.ctypeslib.as_array(H_c)
vol = vol/ncount
rho = (N*(N-1))/(2*vol)
###############
# Close dcd file
###############
fort.closedcd(funit)

###############
# Normalize rdf
###############
x = []
y = []
for i in range(1,nB-cushion):
    r3=((i+.5)*(i+.5)*(i+.5)-(i-.5)*(i-.5)*(i-.5))*dr*dr*dr
    volr=(4./3.)*np.pi*r3
    x.append((i*dr/10))
    y.append((H[i])/(ncount*volr*rho))
   
###############
# Plot and save data
###############
plt.plot(x,y)
plt.savefig('rdfc.png')
np.savetxt('xc.txt',x)
np.savetxt('yc.txt',y)
np.savetxt('Hc.txt',H)
print ('The script took {0} second !'.format(time.time() - startTime))
