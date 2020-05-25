#================================================================#
#                                                                #
#                       maise-net  TOOLBOX                       #
#                                                                #
#                    Python interface to maise                   #
#                  for automated data generation                 #
#                                                                #
#  Script for generation of the K-mesh for existing POSCAR file  #
#                                                                #
#                             -----                              #
#                                                                #
#         With any use of this script, please cite:              #
#                                                                #
#         'Stratified construction of neural network based       #
#          interatomic models for multicomponent materials'      #
#                                                                #
#          DOI: 10.1103/PhysRevB.95.014114                       #
#                                                                #
#                             -----                              #
#                                                                #
#                   Questions and bug reports:                   #
#                                                                #
#        Samad  Hajinazar         hajinazar@binghamton.edu       #
#        Alexey Kolmogorov        kolmogorov@binghamton.edu      #
#                                                                #
#================================================================#

import sys
import errno
import math

# ======================================================================================
# Parameters and structures
# ======================================================================================

# Constants
MAX_ATOMS = 500
DIM       = 3

#----------------------------------------------------------------------------------

# POSCAR data class
class poscar:
    def __init__(self):
        self.name       = "maise"                # name of structure
        self.scale      = 1.0                    # scaling factor
        self.ENE        = 0.0                    # structure energy
        self.Natoms     = 0                      # number of atoms
        self.Ntypes     = 0                      # number of species
        self.Ntypeatoms = []                     # number of atoms of each species
        self.SPC        = []                     # symbol of species
        self.LAT        = array2(DIM,DIM)        # lattice parameters
        self.POS        = array2(MAX_ATOMS,DIM)  # atomic positions
        self.FRC        = array2(MAX_ATOMS,DIM)  # atomic forces
        self.cor        = 0                      # cartesian(0) or direct(1) coor.
        self.sel        = 0                      # non-selective(0) or selective dynam.(1)
        self.selective  = array2(MAX_ATOMS,DIM)  # selective dynamics tags
        self.REC        = array2(DIM,DIM)        # reciprocal lattice vectors

# ======================================================================================
# Routines
# ======================================================================================

# K-mesh generation routine
def mesh(struc,m):
    d=[0.0,0.0,0.0]
    r=[0.0,0.0,0.0]
    k=[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
    R=0.0
    K=0.0
    mesh=[0.0,0.0,0.0]
    M=[0,0,0]

    car_dir(struc)

    for i in range(0,DIM):
        d[i]=math.sqrt(struc.LAT[i][0]*struc.LAT[i][0]+struc.LAT[i][1]*struc.LAT[i][1]+struc.LAT[i][2]*struc.LAT[i][2]);
        
    k[0][0] =  struc.LAT[1][1]*struc.LAT[2][2]-struc.LAT[1][2]*struc.LAT[2][1];
    k[0][1] = -struc.LAT[1][0]*struc.LAT[2][2]+struc.LAT[1][2]*struc.LAT[2][0];
    k[0][2] =  struc.LAT[1][0]*struc.LAT[2][1]-struc.LAT[1][1]*struc.LAT[2][0];
    k[1][0] =  struc.LAT[2][1]*struc.LAT[0][2]-struc.LAT[2][2]*struc.LAT[0][1];
    k[1][1] = -struc.LAT[2][0]*struc.LAT[0][2]+struc.LAT[2][2]*struc.LAT[0][0];
    k[1][2] =  struc.LAT[2][0]*struc.LAT[0][1]-struc.LAT[2][1]*struc.LAT[0][0];
    k[2][0] =  struc.LAT[0][1]*struc.LAT[1][2]-struc.LAT[0][2]*struc.LAT[1][1];
    k[2][1] = -struc.LAT[0][0]*struc.LAT[1][2]+struc.LAT[0][2]*struc.LAT[1][0];
    k[2][2] =  struc.LAT[0][0]*struc.LAT[1][1]-struc.LAT[0][1]*struc.LAT[1][0];
  
    for i in range(0,DIM):
        for j in range(0,DIM):
            k[i][j] /= volume(struc)

    for i in range(0,DIM):
        r[i]=math.sqrt(k[i][0]*k[i][0]+k[i][1]*k[i][1]+k[i][2]*k[i][2]);

    R=1.0
    for i in range(0,DIM):
        R *= r[i];

    COS=0.0
    for i in range(0,DIM):
        COS += (struc.LAT[0][i]*struc.LAT[1][i]/(d[0]*d[1]));

    K=float(m)/float(struc.Natoms);
    
    for i in range(0,DIM):
        mesh[i]=math.pow(K*R*R,0.333333333333333333)*r[i]/R;
    print(str(mesh[0])+"  "+str(mesh[1])+"  "+str(mesh[2]))
    for i in range(0,DIM):
        if int(round(mesh[i]))%2 == 0:
            M[i]=int(round(mesh[i]))
        else:
            M[i]=int(round(mesh[i]))+1
            
    for i in range(0,DIM):
        if m == 0 or M[i] == 0:
            M[i] = 1

    f=open("KPOINTS","w")
    f.write("KPOINTS file by MAISE % d\n" % m)
    f.write("0\n")

    if (abs(d[0]-d[1]) < 1e-6 and abs(COS+0.5) < 1e-6) or M[0]+M[1]+M[2]==3:
        f.write("G\n")
    else:
        f.write("M\n")
    f.write(" % 5d % 5d % 5d\n" % (M[0],M[1],M[2]))
    f.write(" % 5d % 5d % 5d\n" % (0,0,0))
    f.close();

#----------------------------------------------------------------------------------

# Convert POSCAR from direct (fractional) to Cartesian coordinates
def car_dir(struc):
    x=[0.0,0.0,0.0]
    if struc.cor == 1:
        return

    reciprocal(struc)

    for i in range(0,struc.Natoms):
        for q in range(0,DIM):
            x[q] = 0.0
            for j in range(0,DIM):
                x[q] += (struc.POS[i][j]*struc.REC[q][j])
        for q in range(0,DIM):
            while x[q] < 1.0:
                x[q] += 1.0
            while x[q] >= 1.0:
                x[q] -= 1.0
            struc.POS[i][q] = x[q]
    struc.cor=1

#----------------------------------------------------------------------------------

# Read a POSCAR file
def read_poscar(fname):
    struc=poscar()
    struc.name=nthline(fname,0)
    struc.scale=float(nthline(fname,1))
    for i in range(2,5):
        tmp=nthline(fname,i)
        for j in range(0,DIM):
            struc.LAT[i-2][j]=float(tmp.split()[j])

    tmp=nthline(fname,5)
    try:
        val=int(tmp.split()[0])
        return 0 # error: no species!
    except ValueError:
        struc.Ntypes=int(len(tmp.split()))
        for i in range(0,struc.Ntypes):
            struc.SPC.append(tmp.split()[i])
    
    tmp=nthline(fname,6)
    for i in range(0,struc.Ntypes):
        struc.Ntypeatoms.append(int(tmp.split()[i]))
        struc.Natoms+=int(tmp.split()[i])

    tmp=nthline(fname,7) # is there Selective dynamics?
    if tmp[:1] == "S" or tmp[:1] == "s":
        struc.sel=1
        tmp=nthline(fname,8)
    else:
        struc.sel=0
    if tmp[:1] == "d" or tmp[:1] == "D":
        struc.cor=1
    else:
        struc.cor=0

    for i in range(0,struc.Natoms):
        tmp=nthline(fname,8+struc.sel+i)
        for j in range(0,DIM):
            struc.POS[i][j]=float(tmp.split()[j])
        if struc.sel == 1:
            for j in range(0,DIM):
                struc.selective[i][j]=tmp.split()[j+3]
    return struc

#----------------------------------------------------------------------------------                

# Return volume of a 3D cell
def volume(C):
    v = float(C.LAT[0][0] * C.LAT[1][1] * C.LAT[2][2] + C.LAT[0][2] * C.LAT[1][0] * C.LAT[2][1] + C.LAT[0][1] * C.LAT[1][2] * C.LAT[2][0] - C.LAT[0][0] * C.LAT[1][2] * C.LAT[2][1] - C.LAT[0][1] * C.LAT[1][0] * C.LAT[2][2] - C.LAT[0][2] 
* C.LAT[1][1] * C.LAT[2][0])
    return v

#----------------------------------------------------------------------------------

# Return the nth line of a file
def nthline(fname,n):
    with open(fname) as f:
        line=f.read().split("\n") [n]
    return line

#----------------------------------------------------------------------------------

# Define a 2-D array
def array2(x,y):
    return [[0 for j in range(y)]for i in range(x)]


# ======================================================================================
# Main task
# ======================================================================================

# Read a POSCAR file with k-mesh density in command line; export the KPOINTS file
s=read_poscar("POSCAR")
mesh(s,int(sys.argv[1]))
