#================================================================#
#                                                                #
#                       maise-net  TOOLBOX                       #
#                                                                #
#                    Python interface to maise                   #
#                  for automated data generation                 #
#                                                                #
#     Various utility routines (math, memory, unit cell, etc)    #
#                                                                #
#                             -----                              #
#                                                                #
#         With any use of this script, please cite:              #
#                                                                #
#         'MAISE: Construction of neural network interatomic     #
#             models and evolutionary structure optimization'    #
#                                                                #
#         arXiv:2005.12131                                       #
#                                                                #
#                             -----                              #
#                                                                #
#                   Questions and bug reports:                   #
#                                                                #
#        Samad  Hajinazar         hajinazar@binghamton.edu       #
#        Alexey Kolmogorov        kolmogorov@binghamton.edu      #
#                                                                #
#================================================================#

import os
import shutil
import errno
import math
from glob import glob
from subprocess import Popen, PIPE
from maisenet_defs import *

# ======================================================================================           
# Unit cell and crystal structure related routines 
# ======================================================================================         

def smpl_kmesh(struc, m, path = ".", ndim = 3): # K-mesh generation routine; default ndim=3; for ndim=2 K_z=1

    d    = [0.0,0.0,0.0]
    r    = [0.0,0.0,0.0]
    M    = [0,0,0]
    k    = mkarray_2D(3, 3, var = 0.0)
    R    = 0.0
    K    = 0.0
    mesh = [0.0,0.0,0.0]
    coor = struc.COOR

    if coor == 0:
        carttofrac(struc)
        struc.COOR = 1

    for i in range(0,DIM):
        d[i] = math.sqrt(struc.LATT[i][0]*struc.LATT[i][0]+struc.LATT[i][1]*struc.LATT[i][1]+struc.LATT[i][2]*struc.LATT[i][2])
        
    k[0][0] =  struc.LATT[1][1]*struc.LATT[2][2]-struc.LATT[1][2]*struc.LATT[2][1]
    k[0][1] = -struc.LATT[1][0]*struc.LATT[2][2]+struc.LATT[1][2]*struc.LATT[2][0]
    k[0][2] =  struc.LATT[1][0]*struc.LATT[2][1]-struc.LATT[1][1]*struc.LATT[2][0]
    k[1][0] =  struc.LATT[2][1]*struc.LATT[0][2]-struc.LATT[2][2]*struc.LATT[0][1]
    k[1][1] = -struc.LATT[2][0]*struc.LATT[0][2]+struc.LATT[2][2]*struc.LATT[0][0]
    k[1][2] =  struc.LATT[2][0]*struc.LATT[0][1]-struc.LATT[2][1]*struc.LATT[0][0]
    k[2][0] =  struc.LATT[0][1]*struc.LATT[1][2]-struc.LATT[0][2]*struc.LATT[1][1]
    k[2][1] = -struc.LATT[0][0]*struc.LATT[1][2]+struc.LATT[0][2]*struc.LATT[1][0]
    k[2][2] =  struc.LATT[0][0]*struc.LATT[1][1]-struc.LATT[0][1]*struc.LATT[1][0]
  
    for i in range(0,DIM):
        for j in range(0,DIM):
            k[i][j] /= volume_str(struc)

    for i in range(0,DIM):
        r[i] = math.sqrt(k[i][0]*k[i][0]+k[i][1]*k[i][1]+k[i][2]*k[i][2])

    R = 1.0
    for i in range(0,DIM):
        R *= r[i]

    COS = 0.0
    for i in range(0,DIM):
        COS += (struc.LATT[0][i]*struc.LATT[1][i]/(d[0]*d[1]))

    K = float(m)/float(struc.NATM)
    
    for i in range(0,DIM):
        mesh[i] = math.pow(K*R*R,0.333333333333333333)*r[i]/R

    for i in range(0,DIM):
        if int(round(mesh[i]))%2 == 0:
            M[i] = int(round(mesh[i]))
        else:
            M[i] = int(round(mesh[i]))+1
            
    for i in range(0,DIM):
        if m == 0 or M[i] == 0:
            M[i] = 1

    f = open(path+"/"+"KPOINTS","w")
    f.write("KPOINTS file by MAISE % d\n" % m)
    f.write("0\n")
    if (abs(d[0]-d[1]) < 1e-6 and abs(COS+0.5) < 1e-6) or M[0]+M[1]+M[2] == 3:
        f.write("G\n")
    else:
        f.write("M\n")
    if ndim == 3:
        f.write(" % 5d % 5d % 5d\n" % (M[0],M[1],M[2]))
    if ndim == 2:
        f.write(" % 5d % 5d % 5d\n" % (M[0],M[1],1))
    f.write(" % 5d % 5d % 5d\n" % (0,0,0))
    f.close()

    if coor == 0:
        fractocart(struc)
        struc.COOR = 0

#----------------------------------------------------------------------------------                

def reciprocal(struc): # Calculate the reciprocal lattice vectors
    vectr_prod(struc.LATT[0], struc.LATT[1], struc.RECI[2])
    vectr_prod(struc.LATT[1], struc.LATT[2], struc.RECI[0])
    vectr_prod(struc.LATT[2], struc.LATT[0], struc.RECI[1])
    t = float(1.0 / cross_prod(struc.LATT[0], struc.LATT[1], struc.LATT[2]))
    for i in range(0, DIM):
        for q in range(0, DIM):
            struc.RECI[i][q] *= t

#----------------------------------------------------------------------------------                

def fractocart(struc): # Convert POSCAR from direct (fractional) to Cartesian coordinates
    x = [0.0, 0.0, 0.0]
    if struc.COOR == 0:
        print("Already in Cartesian coordinates!")
        return
    for i in range(0, struc.NATM):
        for q in range(0, DIM):
            x[q] = 0.0
            for j in range(0, DIM):
                x[q] += struc.POSI[i][j] * struc.LATT[j][q]

        for q in range(0, DIM):
            struc.POSI[i][q] = x[q]
    struc.COOR = 0

#----------------------------------------------------------------------------------                

def carttofrac(struc): # Convert POSCAR from Cartesian to direct (fractional) coordinates
    x = [0.0, 0.0, 0.0]
    if struc.COOR == 1:
        print("Already in direct (fractional) coordinates!")
        return
    reciprocal(struc)
    for i in range(0, struc.NATM):
        for q in range(0, DIM):
            x[q] = 0.0
            for j in range(0, DIM):
                x[q] += struc.POSI[i][j] * struc.RECI[q][j]

        for q in range(0, DIM):
            while x[q] < 0.0:
                x[q] += 1.0

            while x[q] > 1.0:
                x[q] -= 1.0

            struc.POSI[i][q] = x[q]
    struc.COOR = 1

#----------------------------------------------------------------------------------                

def volume_str(struc): # Return volume of a 3D cell
    return float(struc.LATT[0][0] * struc.LATT[1][1] * struc.LATT[2][2] + struc.LATT[0][2] * struc.LATT[1][0] * struc.LATT[2][1] + struc.LATT[0][1] * struc.LATT[1][2] * struc.LATT[2][0] - struc.LATT[0][0] * struc.LATT[1][2] * struc.LATT[2][1] - struc.LATT[0][1] * struc.LATT[1][0] * struc.LATT[2][2] - struc.LATT[0][2] * struc.LATT[1][1] * struc.LATT[2][0])

#----------------------------------------------------------------------------------                

def atomic_mas(i): # Return atomic mass for Z=i; for 1>Z>95 returns 0.0
    if i in Latom_mass:
        return Latom_mass[i]
    else:
        return 0.0

#----------------------------------------------------------------------------------                

def atomtosymb(i): # Return atomic symbol for Z=i; for 1>Z>95 reutrns "--"
    if i in Latom_symb:
        return Latom_symb[i]
    else:
        return "--"

#----------------------------------------------------------------------------------                

def symbtoatom(s): # Retrun atomic number (Z) for atom symbol s; for 1>Z>95 returns 0
    if s in Lsymb_atom:
        return Lsymb_atom[s]
    else:
        return 0

#----------------------------------------------------------------------------------                

def atomic_rad(i): # Retrun atomic radius for Z=i; for 1>Z>95 returns 0.0
    if i in Latom_radi:
        return Latom_radi[i]
    else:
        return 0.0

# ======================================================================================           
# Memory allocation routines
# ======================================================================================         

def mkarray_2D(x, y, var = 0): # Define a 2-D array (x fastest index) of "var"; if no "var" array filled with 0
    return [ [ var for j in range(y) ] for i in range(x) ]

#----------------------------------------------------------------------------------                

def mkarray_3D(x, y, z, var = 0): # Define a 3-D array (x fastest index) of "var"; if no "var" array filled with 0
    return [ [ [ var for k in range(z) ] for j in range(y) ] for i in range(x) ]

# ======================================================================================           
# String and file handling routines
# ======================================================================================         

def maise_vers(addr): # find and return the version of maise executable at add/ as a float number: ##.##
    td  = "tmpmaisever"
    cwd = os.getcwd()

    if doesntexst(addr+"/maise"):
        return -1.0
    mkdir_tree(td)
    duplicates(addr+"/maise", td)
    os.chdir(td)
    os.system("./maise > out 2>&1")
    a,b,c = findfststr("out", "version")
    if a < 1:
        os.chdir(cwd)
        removes_fd(td)
        return 0.0
    s = c.split()
    if s[2][0] != "m":
        os.chdir(cwd)
        removes_fd(td)
        return 0.0
    v = s[2].split(".")
    if v[0] != "maise" or len(v) != 4 or not (is_integer(v[1]) and is_integer(v[2]) and is_integer(v[3])):
        os.chdir(cwd)
        removes_fd(td)
        return 0.0

    os.chdir(cwd)
    removes_fd(td)
    return float(v[1])*10+float(v[2])+float(v[3])/100.0

#----------------------------------------------------------------------------------                

def find_betwn(s, first, last): # Retruns sub-string in (s) between two markers (first,last) with num of words in that
    if first not in s:
        first = ""
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        if last not in s:
            return s[start:]
        else:
            return ""

#----------------------------------------------------------------------------------                

def file_nmlns(fname): # Return total number of lines for a file; -1 if does not exist
    if doesntexst(fname):
        return -1

    with open(fname) as f:
        for i, l in enumerate(f):
            pass

    return i + 1

#----------------------------------------------------------------------------------                

def file_nthln(fname, n): # Return the nth line of a file (counting from 0); "" if does not exist
    if doesntexst(fname):
        return ""
    
    with open(fname) as f:
        line = f.read().split("\n")[n]

    return line

#----------------------------------------------------------------------------------                

def findallstr(fname, string): # All occurr. of string in file: (num of occur.,line numbers,lines); num=-1 no file; num=0 no occur.
    if doesntexst(fname):
        return (-1,-1,"")

    n = 0
    l = []
    s = []

    with open(fname) as f:
        for i, line in enumerate(f, 0):
            if string in line:
                l.append(i)
                n += 1

    for i in range(0,n):
        s.append(file_nthln(fname,l[i]))

    return (n, l, s)

#----------------------------------------------------------------------------------                

def findlststr(fname, string): # Last occur. of a string in file: (j=1,line number, line); j=-1 no file; j=0 no occur.
    if doesntexst(fname):
        return (-1,-1,"")

    n = 0
    l = -1
    s = ""

    with open(fname) as f:
        for i, line in enumerate(f, 0):
            if string in line:
                n = 1
                l = i

    if n > 0:
        s = file_nthln(fname,l)

    return (n, l, s)

#----------------------------------------------------------------------------------                

def findfststr(fname, string): # Fisrt occur. of a string in file: (j=1,line number, line); j=-1 no file; j=0 no occur.
    if doesntexst(fname):
        return (-1,-1,"")

    n = 0
    l = -1
    s = ""

    with open(fname) as f:
        for i, line in enumerate(f, 0):
            if string in line and n == 0:
                n = 1
                l = i

    if n > 0:
        s = file_nthln(fname,l)

    return (n, l, s)

#----------------------------------------------------------------------------------                

def addstrfile(fname, string): # Adds a single line to the end of a file; file will be created if does not exist
    with open(fname, "a") as f:
        f.write("%s\n" % string)

#----------------------------------------------------------------------------------                

def replacenln(fname, ln, str1): # Replace nth line of file "fname" with "str1"
    if doesntexst(fname):
        return -1
    n = file_nmlns(fname)
    if n < (ln+1):
        return -2
    lines = open(fname, 'r').readlines()
    lines[ln] = str1+"\n"
    f = open(fname, 'w')
    f.writelines(lines)
    f.close()
    return 1

#----------------------------------------------------------------------------------                

def replacestr(fname, str1, str2): # Replace all instances of the str1 with str2 in fname file; return -1 if file does not exist
    if doesntexst(fname):
        return -1
    with open(fname, "r") as f:
        filedata = f.read()
    filedata = filedata.replace(str1, str2)
    with open(fname, "w") as f:
        f.write(filedata)
    return 1

#----------------------------------------------------------------------------------                

def eliminates(fname, s): # Eliminate all lines that contain string "s" from file "fname"; return -1 if file does not exist
    if doesntexst(fname):
        return -1
    f = open(fname,"r")
    O = f.read().split("\n")
    f.close()
    f = open(fname,"w")
    for i in range(0,len(O)-1):
        if not s in O[i]:
            f.write("%s\n" % O[i])
    f.close()
    return 1

#----------------------------------------------------------------------------------                

def doesntexst(name): # returns TRUE if file/directory DOES NOT EXIST!
    if not os.path.exists(name):
        return True
    else:
        return False

#----------------------------------------------------------------------------------                

def does__exst(name): # returns TRUE if file/directory EXISTS!
    if os.path.exists(name):
        return True
    else:
        return False

#----------------------------------------------------------------------------------                

def renames_fd(src,dst): # renames file/directory f0 to f1 if that exists! return 0 if file names are the same
    if src == dst:
        return 0
    if doesntexst(src):
        return -1
    os.rename(src,dst)
    return 1

#----------------------------------------------------------------------------------                

def duplicates(src, dst): # Copy any [list of] file or directory (with wildcards)

#    n = len(src.split())
#    if n > 1 and doesntexst(dst):
#        mkdir_tree(dst)
#    for i in range(0,n):
#        fil = src.split()[i]
#        lst = glob(fil)
#        for path in lst:
#            if does__exst(path):
#                os.system("cp -r "+path+" "+dst)
    if doesntexst(src):
        return -1
    try:
        shutil.copytree(src, dst)
    except OSError as exc:
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else:
            raise
    return 1

#----------------------------------------------------------------------------------                

def removes_fd(name): # removes a [list of] file or directory (no matter full/write-protected/wildcards in names)
    n = len(name.split())
    for i in range(0,n):
        fil = name.split()[i]
        lst = glob(fil)
        for path in lst:
            if does__exst(path):
                if os.path.isdir(path):
                    shutil.rmtree(path)
                if os.path.isfile(path):
                    os.remove(path)

#----------------------------------------------------------------------------------                

def mkdir_tree(path): # Makes directory tree if it does not exist; nothing happens if it exists!
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

#----------------------------------------------------------------------------------                

def searchfile(address, fname): # Searches the "address/..." for "fname" (with wildcards); results are sorted: (#,dirs,files)

    alldir  = glob(address)
    subdirs = []
    for n in range(0,len(alldir)):
        for i,j,k in os.walk(alldir[n]):
            subdirs.append(i)

    subdirs.sort()

    z = []
    for i in range(0,len(subdirs)):
        if len(glob(subdirs[i]+"/"+fname)) > 0:
            z += glob(subdirs[i]+"/"+fname)

    files = []
    direc = []
    for i in range(0,len(z)):
        a = z[i].split("/")
        files.append(a[len(a)-1])
        s = ""
        for j in range(0,len(a)-1):
            s += (a[j]+"/")
        direc.append(s)

    return (len(direc), direc, files)

# ======================================================================================           
# Math routines
# ======================================================================================         

def is_integer(s): # Check if input string is an integer number: True/False
    return s.isdigit()

#----------------------------------------------------------------------------------                

def is___float(s): # check if input string is a float number: True/False
    return s.replace('.', '', 1).isdigit()

#----------------------------------------------------------------------------------                

def vectr_prod(a, b, c): # Vector product of 2 vectors: c=axb
    c[0] = a[1] * b[2] - a[2] * b[1]
    c[1] = a[2] * b[0] - a[0] * b[2]
    c[2] = a[0] * b[1] - a[1] * b[0]


#----------------------------------------------------------------------------------                

def cross_prod(a, b, c): # Cross product of 3 vectors: = axb.c
    return a[0] * b[1] * c[2] + a[1] * b[2] * c[0] + a[2] * b[0] * c[1] - a[0] * b[2] * c[1] - a[1] * b[0] * c[2] - a[2] * b[1] * c[0]

# ======================================================================================           
# System tools
# ======================================================================================         

def is_install(name): # Check if a program is intalled on OS; returns logical output (TRUE or FALSE)
    p = Popen(['/usr/bin/which', name], stdout=PIPE, stderr=PIPE)
    p.communicate()
    return p.returncode == 0
