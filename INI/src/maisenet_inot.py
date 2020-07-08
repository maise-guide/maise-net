#================================================================#
#                                                                #
#                       maise-net  TOOLBOX                       #
#                                                                #
#                    Python interface to maise                   #
#                  for automated data generation                 #
#                                                                #
#                    Input and output routines                   #
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
import datetime
import time
from maisenet_defs  import *
from maisenet_util  import *

#----------------------------------------------------------------------------------    

def import_stp(fname): # Read the data generation setup
    struc = setup()
    with open(fname) as f:
        for i, line in enumerate(f, 0):
            s    = find_betwn(line,"","#").split()
            if len(s) > 1:
                flag  = s[0]
                multi = s[1:]    # input source for variables which are list
                singl = multi[0] # input source for variables with single value
                # int vars
                if "WAIT" == flag:
                    struc.WAIT = int(singl)
                if "TMIN" == flag:
                    struc.TMIN = int(singl)
                if "TMIN" == flag:
                    struc.TMIN = int(singl)
                if "TITR" == flag:
                    struc.TITR = int(singl)
                if "TNPP" == flag:
                    struc.TNPP = int(singl)
                if "JOBT" == flag:
                    struc.JOBT = int(singl)
                if "MAXJ" == flag:
                    struc.MAXJ = int(singl)
                if "SMER" == flag:
                    struc.SMER = int(singl)
                if "smer" == flag:
                    struc.smer = int(singl)
                if "AEOS" == flag:
                    struc.AEOS = int(singl)
                if "ATST" == flag:
                    struc.ATST = int(singl)
                if "SITR" == flag:
                    struc.SITR = int(singl)
                if "ITER" == flag:
                    struc.ITER = int(singl)
                if "NITR" == flag:
                    struc.NITR = int(singl)
                if "DATA" == flag:
                    struc.DATA = int(singl)
                if "SORT" == flag:
                    struc.SORT = int(singl)
                if "SMLR" == flag:
                    struc.SMLR = int(singl)
                if "NPAR" == flag:
                    struc.NPAR = int(singl)
                if "RPAR" == flag:
                    struc.RPAR = int(singl)
                if "NSYM" == flag:
                    struc.NSYM = int(singl)
                if "RCUT" == flag:
                    struc.RCUT = int(singl)
                if "NDIM" == flag:
                    struc.NDIM = int(singl)
                if "LSYM" == flag:
                    struc.LSYM = int(singl)
                if "nsw0" == flag:
                    struc.nsw0 = int(singl)
                if "kdns" == flag:
                    struc.kdns = int(singl)
                if "KDNS" == flag:
                    struc.KDNS = int(singl)
                if "EXTR" == flag:
                    struc.EXTR = int(singl)
                if "QUET" == flag:
                    struc.QUET = int(singl)
                if "ERNI" == flag:
                    struc.ERNI = int(singl)
                if "NGEN" == flag:
                    struc.NGEN = int(singl)
                if "ngen" == flag:
                    struc.ngen = int(singl)
                if "SORB" == flag:
                    struc.SORB = int(singl)
                # float vars
                if "SYMP" == flag:
                    struc.SYMP = float(singl)
                if "SIGM" == flag:
                    struc.SIGM = float(singl)
                if "sigm" == flag:
                    struc.sigm = float(singl)
                if "LREG" == flag:
                    struc.LREG = float(singl)
                if "FMRK" == flag:
                    struc.FMRK = float(singl)
                if "LBOX" == flag:
                    struc.LBOX = float(singl)
                if "TETR" == flag:
                    struc.TETR = float(singl)
                if "PLNT" == flag:
                    struc.PLNT = float(singl)
                if "PACK" == flag:
                    struc.PACK = float(singl)
                if "BLOB" == flag:
                    struc.BLOB = float(singl)
                if "MATE" == flag:
                    struc.MATE = float(singl)
                if "SWAP" == flag:
                    struc.SWAP = float(singl)
                if "RUBE" == flag:
                    struc.RUBE = float(singl)
                if "REFL" == flag:
                    struc.REFL = float(singl)
                if "INVS" == flag:
                    struc.INVS = float(singl)
                if "CHOP" == flag:
                    struc.CHOP = float(singl)
                if "MUTE" == flag:
                    struc.MUTE = float(singl)
                if "MCRS" == flag:
                    struc.MCRS = float(singl)
                if "SCRS" == flag:
                    struc.SCRS = float(singl)
                if "LCRS" == flag:
                    struc.LCRS = float(singl)
                if "ACRS" == flag:
                    struc.ACRS = float(singl)
                if "SDST" == flag:
                    struc.SDST = float(singl)
                if "LDST" == flag:
                    struc.LDST = float(singl)
                if "ADST" == flag:
                    struc.ADST = float(singl)
                if "ELPS" == flag:
                    struc.ELPS = float(singl)
                if "ecut" == flag:
                    struc.ecut = float(singl)
                if "ECUT" == flag:
                    struc.ECUT = float(singl)
                # string vars
                if "MAIL" == flag:
                    struc.MAIL = singl
                if "prec" == flag:
                    struc.prec = singl
                if "PREC" == flag:
                    struc.PREC = singl
                # int lists
                if "TNGN" == flag:
                    struc.TNGN = list(map(int,multi))
                if "TSTA" == flag:
                    struc.TSTA = list(map(int,multi))
                if "TSTB" == flag:
                    struc.TSTB = list(map(int,multi))
                if "TSTC" == flag:
                    struc.TSTC = list(map(int,multi))
                if "NNNU" == flag:
                    struc.NNNU = list(map(int,multi))
                if "NNGT" == flag:
                    struc.NNGT = list(map(int,multi))
                if "TSPC" == flag:
                    struc.TSPC = list(map(int,multi))
                if "aspc" == flag:
                    struc.aspc = list(map(int,multi))
                if "bspc" == flag:
                    struc.bspc = list(map(int,multi))
                if "cspc" == flag:
                    struc.cspc = list(map(int,multi))
                if "npop" == flag:
                    struc.npop = list(map(int,multi))
                if "mitr" == flag:
                    struc.mitr = list(map(int,multi))
                if "tefs" == flag:
                    struc.tefs = list(map(int,multi))
                if "ASPC" == flag:
                    struc.ASPC = list(map(int,multi))
                if "BSPC" == flag:
                    struc.BSPC = list(map(int,multi))
                if "CSPC" == flag:
                    struc.CSPC = list(map(int,multi))
                if "NPOP" == flag:
                    struc.NPOP = list(map(int,multi))
                if "MITR" == flag:
                    struc.MITR = list(map(int,multi))
                if "TEFS" == flag:
                    struc.TEFS = list(map(int,multi))
                # float lists
                if "TGPA" == flag:
                    struc.TGPA = list(map(float,multi))
                if "TPWT" == flag:
                    struc.TPWT = list(map(float,multi))
                if "PGPA" == flag:
                    struc.PGPA = list(map(float,multi))
                if "PWGT" == flag:
                    struc.PWGT = list(map(float,multi))

    return struc

#----------------------------------------------------------------------------------    

def reread_stp(fileout, setup, fname): # Re-reads setup file for data generation; quits for JOBT=89; pause for JOBT=88 
    s = import_stp(fname)
    if s.JOBT%10 == JPSE:
        export_out(fileout,"Note: pausing ...",color = 3)
        while s.JOBT%10 == JPSE:
            time.sleep(5)
            s = import_stp(fname)
        export_out(fileout,"Note: resuming ...",color = 3)
    if s.JOBT%10 == JQTE:
        export_out(fileout,"Note: terminating at user's request!",color = 3);exit()
    setup.MAXJ = s.MAXJ
    setup.WAIT = s.WAIT

#----------------------------------------------------------------------------------                

def export_out(fileout, message, time = 1, color = 0): # Write to the output file
    f = open(fileout, "a")
    if time == 1:
        if color > 0:
            f.write("% 20s > %s % s %s \n\n" % (datetime.datetime.now().strftime("%m/%d/%Y %H:%M"),COLOR_CODES[color],message,COLOR_CODES[0]))
            print("% 20s > %s % s %s \n" % (datetime.datetime.now().strftime("%m/%d/%Y %H:%M"),COLOR_CODES[color],message,COLOR_CODES[0]))
        else:
            f.write(" % 20s > % s \n\n" % (datetime.datetime.now().strftime("%m/%d/%Y %H:%M"),message))
            print(" % 20s > % s \n" % (datetime.datetime.now().strftime("%m/%d/%Y %H:%M"),message))
    else:
        if color > 0:
            f.write(COLOR_CODES[color]+" "+message +" "+COLOR_CODES[0]+ "\n")
            print(COLOR_CODES[color]+" "+message +" "+COLOR_CODES[0])
        else:
            f.write(message + "\n")
            print(message)
    f.close()

#----------------------------------------------------------------------------------                

def import_str(fname): # Read a POSCAR file
    struc = poscar()
    struc.NAME = file_nthln(fname, 0)
    struc.SCAL = float(file_nthln(fname, 1))
    for i in range(2, 5):
        tmp = file_nthln(fname, i)
        for j in range(0, DIM):
            struc.LATT[i - 2][j] = float(tmp.split()[j]) * struc.SCAL

    struc.SCAL /= struc.SCAL

    tmp = file_nthln(fname, 5)
    try:
        val = int(tmp.split()[0])
        return 0
    except ValueError:
        struc.NTYP = int(len(tmp.split()))
        for i in range(0, struc.NTYP):
            struc.SPCS.append(tmp.split()[i])

    tmp = file_nthln(fname, 6)
    for i in range(0, struc.NTYP):
        struc.NAET.append(int(tmp.split()[i]))
        struc.NATM += int(tmp.split()[i])

    tmp = file_nthln(fname, 7)
    if tmp[:1] == "S" or tmp[:1] == "s":
        struc.SLTV = 1
        tmp = file_nthln(fname, 8)
    else:
        struc.SLTV = 0
    if tmp[:1] == "d" or tmp[:1] == "D":
        struc.COOR = 1
    else:
        struc.COOR = 0
    for i in range(0, struc.NATM):
        tmp = file_nthln(fname, 8 + struc.SLTV + i)
        for j in range(0, DIM):
            struc.POSI[i][j] = float(tmp.split()[j])

        if struc.SLTV == 1:
            for j in range(0, DIM):
                struc.DYNM[i][j] = tmp.split()[j + 3]

    return struc

#----------------------------------------------------------------------------------                

def export_str(struc, fname): # Output the POSCAR structure
    f = open(fname, "w")
    f.write(" %s\n" % struc.NAME)
    f.write(" %lf\n" % struc.SCAL)

    for i in range(0, DIM):
        for j in range(0, DIM):
            f.write(" % 20.14lf " % struc.LATT[i][j])
        f.write("\n")

    for i in range(0, struc.NTYP):
        f.write(" %s " % struc.SPCS[i])
    f.write("\n")

    for i in range(0, struc.NTYP):
        f.write(" %d " % struc.NAET[i])
    f.write("\n")

    if struc.SLTV == 1:
        f.write("Selective dynamics\n")

    coor = struc.COOR
    if coor == 0:
        carttofrac(struc)
        struc.COOR = 1
    f.write("Direct\n")

    for i in range(0, struc.NATM):
        for j in range(0, DIM):
            f.write(" % 20.14lf " % struc.POSI[i][j])
        if struc.SLTV == 1:
            for j in range(0, DIM):
                f.write(" %s " % struc.DYNM[i][j])
        f.write("\n")

    f.close()

    if coor == 0:
        fractocart(struc)
        struc.COOR = 0

#----------------------------------------------------------------------------------

def outcar_out(fileout, struc, dir_out):

    if doesntexst(dir_out):
        export_out(fileout,"Error: target directory for OUTCAR (% s) does not exist" % (dir_file),color = 2);exit()
    nion = struc.NATM
    vole = volume_str(struc)/float(nion)
    tene = struc.ENER
    aene = tene/float(nion)
    inkB = struc.inkB
    enth = struc.ENTH/float(nion)
    pres = struc.PRES

    f = open(dir_out+"/OUTCAR","w")

    f.write("OUTCAR file generated by maise-net script!\n\n")

    f.write(" POSITION                                       TOTAL-FORCE (eV/Angst)\n")
    f.write(" -----------------------------------------------------------------------------------\n")

    coor = struc.COOR
    if coor == 1:
        fractocart(struc)
        struc.COOR = 0

    for i in range(0,struc.NATM):
        for j in range(0,DIM):
            f.write("% lf   " % struc.POSI[i][j])
        for j in range(0,DIM):
            f.write("% lf   " % struc.FORC[i][j])
        f.write("\n")
    f.write(" -----------------------------------------------------------------------------------\n")
    f.write(" total drift: 0 0 0\n\n")
    f.write(" vol-ene   % 16.8lf   % 16.8lf\n" % (vole,aene))
    f.write(" enthalpy  % 16.8lf and total % 16.8lf\n" %(enth,enth*float(nion)))
    f.write(" pressure  % 16.8lf\n\n" %(pres))
    f.write("in kB %s\n\n" % inkB)
    f.write("final total enthalpy=    % 16.8lf  total energy=     % 16.8lf      % 16.8lf     % 16.8lf\n\n" % (enth*float(nion),tene,enth,tene/float(nion)))
    f.write("Total CPU time used (sec): % 16.8lf\n" % tene)

    f.close()

    if coor == 1:
        carttofrac(struc)
        struc.COOR = 1

#----------------------------------------------------------------------------------                

def datdat_out(fileout, dir_file, dir_out): # Write (dir_out/dat.dat) file using (outcar with path: dir_file)
    if doesntexst(dir_file):
        export_out(fileout,"Error: input OUTCAR file (% s) does not exist" % (dir_file),color = 2);exit()
    if doesntexst(dir_out):
        export_out(fileout,"Error: target directory for dat.dat (% s) does not exist" % (dir_file),color = 2);exit()

    f= open(dir_file,"r")
    O=f.read().split("\n")
    f.close()

    N_nion = 0
    N_vole = 0
    N_tene = 0
    N_enth = 0
    N_inkB = 0
    N_frc1 = 0
    N_frc2 = 0

    for i in range(0,len(O)-1):
        if "NIONS" in O[i]:
            N_nion = i
        if "volume of cell" in O[i]:
            N_vole = i
        if "y=" in O[i]:
            N_tene = i
        if "in kB" in O[i]:
            N_inkB = i
        if "POSITION" in O[i]:
            N_frc1 = i
        if "drift" in O[i]:
            N_frc2 = i
        if "enth" in O[i]:
            N_enth = i
        if "PSTRESS=" in O[i]:
            N_pres = i

    nion = int(O[N_nion].split()[11])
    vole = float(O[N_vole].split()[4])/nion
    tene = float(O[N_tene].split()[6])
    aene = tene/nion
    inkB = O[N_inkB]
    pres = float(O[N_pres].split()[1])/10.0

    if N_enth > 0:
        enth = float(O[N_enth].split()[4])/nion
    else:
        enth = aene

    # Write the dat.dat file                        
    f = open(dir_out+"/dat.dat","w")
    f.write("% 16.8lf\n" % tene)
    f.write("%s\n" % inkB)
    for i in range(N_frc1,N_frc2+1):
        f.write("%s\n" % O[i])
    f.write(" vol-ene   % 16.8lf   % 16.8lf\n" % (vole,aene))
    f.write(" enthalpy  % 16.8lf\n" %(enth))
    f.write(" pressure  % 16.8lf\n" %(pres))
    f.close()

#----------------------------------------------------------------------------------                

def datdat_str(fileout, struc, dir_out): # Write (dir_out/dat.dat) file using a "struc" class

    if doesntexst(dir_out):
        export_out(fileout,"Error: target directory for dat.dat (% s) does not exist" % (dir_file),color = 2);exit()

    nion = struc.NATM
    vole = volume_str(struc)/float(nion)
    tene = struc.ENER
    aene = tene/float(nion)
    inkB = struc.inkB
    enth = struc.ENTH/float(nion)
    pres = struc.PRES

    # Write the dat.dat file                        
    f = open(dir_out+"/dat.dat","w")
    f.write("% 16.8lf\n" % tene)
    f.write("in kB %s\n" % inkB)
    f.write(" POSITION                                       TOTAL-FORCE (eV/Angst)\n")
    f.write(" -----------------------------------------------------------------------------------\n")
    
    coor = struc.COOR
    if coor == 1:
        fractocart(struc)
        struc.COOR = 0

    for i in range(0,struc.NATM):
        for j in range(0,DIM):
            f.write("% lf   " % struc.POSI[i][j])
        for j in range(0,DIM):
            f.write("% lf   " % struc.FORC[i][j])
        f.write("\n")
    f.write(" -----------------------------------------------------------------------------------\n")
    f.write(" total drift: 0 0 0\n")
    f.write(" vol-ene   % 16.8lf   % 16.8lf\n" % (vole,aene))
    f.write(" enthalpy  % 16.8lf\n" %(enth))
    f.write(" pressure  % 16.8lf\n" %(pres))
    f.close()

    if coor == 1:
        carttofrac(struc)
        struc.COOR = 1
#----------------------------------------------------------------------------------         

def readoutcar(fname, read = "all"): # returns total number of structures (N) and their info (strucs) from OUTCAR
    (n,l,s)=findallstr(fname,"VRHFIN")
    if n > 0:
        return outcar_vsp(fname, read)
    else:
        return outcar_ann(fname, read)

#----------------------------------------------------------------------------------         

def outcar_ann(fname, read): # Returns total number of structures (N) and their info (strucs) from MAISE OUTCAR
    strucs = []

    if read == "all":
        N , l0, s0 = findallstr(fname, "y=") # to read "total energy/structure"
        n1, l1, s1 = findallstr(fname, "in kB")
        n2, l2, s2 = findallstr(fname, "atoms of each species")
        n3, l3, s3 = findallstr(fname, "VOLUME")
        n4, l4, s4 = findlststr(fname, "pressure")
        n5, l5, s5 = findallstr(fname, "species              ")
        n6, l6, s6 = findallstr(fname, "POSITION")
        n7, l7, s7 = findallstr(fname, "LATTICE")

    if read == "first":
        N , a0, b0 = findfststr(fname, "y=") # to read "total energy/structure"
        n1, a1, b1 = findfststr(fname, "in kB")
        n2, a2, b2 = findfststr(fname, "atoms of each species")
        n3, a3, b3 = findfststr(fname, "VOLUME")
        n4, a4, b4 = findlststr(fname, "pressure")
        n5, a5, b5 = findfststr(fname, "species              ")
        n6, a6, b6 = findfststr(fname, "POSITION")
        n7, a7, b7 = findfststr(fname, "LATTICE")
        l0=[a0];l1=[a1];l2=[a2];l3=[a3];l4=[a4];l5=[a5];l6=[a6];l7=[a7]
        s0=[b0];s1=[b1];s2=[b2];s3=[b3];s4=[b4];s5=[b5];s6=[b6];s7=[b7]


    if read == "last":
        N , a0, b0 = findlststr(fname, "y=") # to read "total energy/structure"
        n1, a1, b1 = findlststr(fname, "in kB")
        n2, a2, b2 = findlststr(fname, "atoms of each species")
        n3, a3, b3 = findlststr(fname, "VOLUME")
        n4, a4, b4 = findlststr(fname, "pressure")
        n5, a5, b5 = findlststr(fname, "species              ")
        n6, a6, b6 = findlststr(fname, "POSITION")
        n7, a7, b7 = findlststr(fname, "LATTICE")
        l0=[a0];l1=[a1];l2=[a2];l3=[a3];l4=[a4];l5=[a5];l6=[a6];l7=[a7]
        s0=[b0];s1=[b1];s2=[b2];s3=[b3];s4=[b4];s5=[b5];s6=[b6];s7=[b7]

    for i in range(0, N):
        strucs.append(poscar())

    for i in range(0, N):
        strucs[i].ENER  = float(s0[i].split()[6])
        strucs[i].ENTH = float(s0[i].split()[4])

    for i in range(0, N):
        a = ""
        for j in range(2,8):
            a += "   "+s1[i].split()[j]
        strucs[i].inkB = a

    for i in range(0, N):
        strucs[i].NTYP = int(len(s2[0].split())) - 4
        for j in range(0, strucs[i].NTYP):
            strucs[i].NAET.append(int(s2[0].split()[4 + j]))
        a = 0
        for j in range(0, strucs[i].NTYP):
            a += strucs[i].NAET[j]
        strucs[i].NATM = a

    for i in range(0, N):
        strucs[i].VOLM = float(file_nthln(fname,l3[i]+2).split()[3]) * float(strucs[i].NATM)

    for i in range(0, N):
        strucs[i].PRES = float(s4.split()[2])

    for i in range(0, N):
        for j in range(0, int(len(s5[0].split()) - 1)):
            strucs[i].SPCS.append(s5[0].split()[1 + j])

    for i in range(0, N):
        for j in range(l6[i] + 2, l6[i] + strucs[i].NATM + 2):
            tmp = file_nthln(fname, j)
            for k in range(0, 3):
                strucs[i].POSI[j - l6[i] - 2][k] = float(tmp.split()[k])
                strucs[i].FORC[j - l6[i] - 2][k] = float(tmp.split()[k + 3])

    for i in range(0, N):
        for j in range(l7[i] + 2, l7[i] + 5):
            tmp = file_nthln(fname, j).split()
            strucs[i].LATT[j - l7[i] - 2][0] = float(tmp[0])
            strucs[i].LATT[j - l7[i] - 2][1] = float(tmp[1])
            strucs[i].LATT[j - l7[i] - 2][2] = float(tmp[2])

    return ("MAIES", N, strucs)

#----------------------------------------------------------------------------------     

def outcar_vsp(fname, read): # returns total number of structures (N) and their info (strucs) from VASP OUTCAR
    strucs = []

    if read == "all":
        N , l0, s0 = findallstr(fname,"y=") #total energy        
        n1, l1, s1 = findallstr(fname,"enth")
        n2, l2, s2 = findallstr(fname, "in kB")
        n3, l3, s3 = findallstr(fname," volume of cell")
        n4, l4, s4 = findallstr(fname,"PSTRESS=")
        n5, l5, s5 = findallstr(fname,"NIONS") #total number of atoms: NATM
        n6, l6, s6 = findallstr(fname,"   ions per type =")  #number of species NTYP
        n7, l7, s7 = findallstr(fname,"direct lattice vectors")  #direct lattice vectors
        n8, l8, s8 = findallstr(fname,"POSITION")  #positions and forces
        n9, l9, s9 = findallstr(fname,"VRHFIN") #types of species

    if read == "first":
        N , a0, b0 = findfststr(fname,"y=") #total energy        
        n1, a1, b1 = findfststr(fname,"enth")
        n2, a2, b2 = findfststr(fname, "in kB")
        n3, a3, b3 = findfststr(fname," volume of cell")
        n4, a4, b4 = findfststr(fname,"PSTRESS=")
        n5, a5, b5 = findfststr(fname,"NIONS") #total number of atoms: NATM
        n6, a6, b6 = findfststr(fname,"   ions per type =")  #number of species NTYP
        n7, a7, b7 = findfststr(fname,"direct lattice vectors")  #direct lattice vectors
        n8, a8, b8 = findfststr(fname,"POSITION")  #positions and forces
        n9, a9, b9 = findfststr(fname,"VRHFIN") #types of species
        l0=[a0];l1=[a1];l2=[a2];l3=[a3];l4=[a4];l5=[a5];l6=[a6];l7=[a7];l8=[a8];l9=[a9]
        s0=[b0];s1=[b1];s2=[b2];s3=[b3];s4=[b4];s5=[b5];s6=[b6];s7=[b7];s8=[b8];s9=[b9]
        n7 *= 2; l7 += l7; s7 += s7

    if read == "last":
        N , a0, b0 = findlststr(fname,"y=") #total energy        
        n1, a1, b1 = findlststr(fname,"enth")
        n2, a2, b2 = findlststr(fname, "in kB")
        n3, a3, b3 = findlststr(fname," volume of cell")
        n4, a4, b4 = findlststr(fname,"PSTRESS=")
        n5, a5, b5 = findlststr(fname,"NIONS") #total number of atoms: NATM
        n6, a6, b6 = findlststr(fname,"   ions per type =")  #number of species NTYP
        n7, a7, b7 = findlststr(fname,"direct lattice vectors")  #direct lattice vectors
        n8, a8, b8 = findlststr(fname,"POSITION")  #positions and forces
        n9, a9, b9 = findlststr(fname,"VRHFIN") #types of species
        l0=[a0];l1=[a1];l2=[a2];l3=[a3];l4=[a4];l5=[a5];l6=[a6];l7=[a7];l8=[a8];l9=[a9]
        s0=[b0];s1=[b1];s2=[b2];s3=[b3];s4=[b4];s5=[b5];s6=[b6];s7=[b7];s8=[b8];s9=[b9]
        n7 *= 2; l7 += l7; s7 += s7
        
    for i in range(0, N):
        strucs.append(poscar())

    for i in range(0, N):
        strucs[i].ENER=float(s0[i].split()[6])

    if n1 != N:
        for i in range(0,N):
            strucs[i].ENTH = strucs[i].ENER
    else:
        for i in range(0, N):
            strucs[i].ENTH = float(s1[i].split()[4])

    for i in range(0, N):
        a = ""
        for j in range(2,8):
            a += "   "+s2[i].split()[j]
        strucs[i].inkB = a

    for i in range(0, N):
        strucs[i].VOLM = float(s3[i].split()[4])

    if n4 == N:
        for i in range(0, N):
            strucs[i].PRES = float(s4[i].split()[1])/10.0

    for i in range(0, N):
        strucs[i].NATM = int(s5[0].split()[11])

    for i in range(0, N):  
        strucs[i].NTYP = int(len(s6[0].split())-4)
    for i in range(0, N):
        for j in range(0,strucs[i].NTYP):
            strucs[i].NAET.append(int(s6[0].split()[4+j]))

    for i in range(1, N+1):
        for j in range(l7[i]+1,l7[i]+4):
            tmp = file_nthln(fname,j)
            if (len(tmp.split()) == 6):
                for k in range(0,3):
                    strucs[i-1].LATT[j-l7[i]-1][k] = float(tmp.split()[k])
            else:
                strucs[i-1].LATT[j-l7[i]-1][0] = float(tmp[0:16])
                strucs[i-1].LATT[j-l7[i]-1][1] = float(tmp[16:29])
                strucs[i-1].LATT[j-l7[i]-1][2] = float(tmp[29:42])

    for i in range(0,N):
        for j in range(l8[i]+2,l8[i]+strucs[i].NATM+2):
            tmp = file_nthln(fname,j)
            for k in range(0,3):
                strucs[i].POSI[j-l8[i]-2][k] = float(tmp.split()[k])
                strucs[i].FORC[j-l8[i]-2][k] = float(tmp.split()[k+3])

    for i in range(0,N):
        for j in range(0,n9):
            start = s9[j].find("=")+1
            end   = s9[j].find(":",start)
            strucs[i].SPCS.append(s9[j][start:end])

    return ("VASP", N, strucs)

#----------------------------------------------------------------------------------                

def export_hdr(output, ver): # Output the header of the script
    
    cwd = os.getcwd()
    
    fileout = cwd+"/"+output

    if doesntexst(fileout):
        export_out(fileout,"      #================================================================#",time = 0,color = 3)
        export_out(fileout,"      #                                                                #",time = 0,color = 3)
        export_out(fileout,"      #                            maise-net                           #",time = 0,color = 3)
        export_out(fileout,"      #                                                                #",time = 0,color = 3)
        export_out(fileout,"      #                    Python interface to maise                   #",time = 0,color = 3)
        export_out(fileout,"      #                  for automated data generation                 #",time = 0,color = 3)
        export_out(fileout,"      #                                                                #",time = 0,color = 3)
        export_out(fileout,"      #"+int(round((64.0-len(ver))/2.0))*" "+ver+(64-int(round((64.0-len(ver))/2.0))-len(ver))*" "+"#",time = 0,color = 3)             
        export_out(fileout,"      #                                                                #",time = 0,color = 3)
        export_out(fileout,"      #                             -----                              #",time = 0,color = 3)
        export_out(fileout,"      #                                                                #",time = 0,color = 3)
        export_out(fileout,"      #         With any use of this script, please cite:              #",time = 0,color = 3)
        export_out(fileout,"      #                                                                #",time = 0,color = 3)
        export_out(fileout,"      #         'MAISE: Construction of neural network interatomic     #",time = 0,color = 3)
        export_out(fileout,"      #             models and evolutionary structure optimization'    #",time = 0,color = 3)
        export_out(fileout,"      #                                                                #",time = 0,color = 3)
        export_out(fileout,"      #         arXiv:2005.12131                                       #",time = 0,color = 3)
        export_out(fileout,"      #                                                                #",time = 0,color = 3)
        export_out(fileout,"      #                             -----                              #",time = 0,color = 3)
        export_out(fileout,"      #                                                                #",time = 0,color = 3)
        export_out(fileout,"      #                   Questions and bug reports:                   #",time = 0,color = 3)
        export_out(fileout,"      #                                                                #",time = 0,color = 3)
        export_out(fileout,"      #        Samad  Hajinazar         hajinazar@binghamton.edu       #",time = 0,color = 3)
        export_out(fileout,"      #        Alexey Kolmogorov        kolmogorov@binghamton.edu      #",time = 0,color = 3)
        export_out(fileout,"      #                                                                #",time = 0,color = 3)
        export_out(fileout,"      #================================================================#",time = 0,color = 3)
        export_out(fileout,"                                                                        ",time = 0,color = 3)

#----------------------------------------------------------------------------------                

def mkplot_err(fileout, setup_g): # Makes histogram and scatter plot for NN testing errors
    if not is_install("gnuplot"):
        export_out(fileout,"Warning: gnuplot does not exist; no plot will be generated",time = 0)
        return

    f = open("sc-gp","w")
    f.write("set terminal png size 800,600\n")
    f.write("set output '%s-ent.png'\n" % (setup_g.NAME))
    f.write("set title '%s'\n" % (os.getcwd()))
    for i in range(0,setup_g.PNUM+1):
        f.write("set linestyle  %d lt 1 lw 2 lc %d pt %d\n" % (i+1,i+1,i+1))
    f.write("set xrange [0:399]\n")
    f.write("set xlabel 'relative enthalpy (meV/atom)'\n")
    f.write("set ylabel 'testing error (meV/atom)'\n")
    f.write("set multiplot\n")
    f.write("plot ")
    for i in range(0,setup_g.PNUM):
        if does__exst("err-ent%d.dat" % i):
            f.write("'err-ent%d.dat' u 2:1 w points ls %d title '%6.2lf'," % (i,i+1,setup_g.PGPA[i]))
    f.write("\n")
    f.write("unset multiplot\n")
    f.close()
    os.system("gnuplot < sc-gp")
    removes_fd("sc-gp")

    f = open("sc-gp","w")
    f.write("set terminal png size 800,600\n")
    f.write("set output '%s-his.png'\n" %(setup_g.NAME))
    f.write("set title '%s'\n" % (os.getcwd()))
    for i in range(0,setup_g.PNUM+1):
        f.write("set linestyle  %d lt 1 lw 2 lc %d pt %d\n" % (i+1,i+1,i+1))
    f.write("width=5\n")
    f.write("hist(x,width)=width*floor(x/width)+width/2.0\n")
    f.write("set boxwidth width*0.9\n")
    f.write("set style fill transparent pattern 4\n")
    f.write("set xlabel 'testing error (meV/atom)'\n")
    f.write("set ylabel 'count'\n")
    f.write("set multiplot\n")
    f.write("plot ")
    for i in range(0,setup_g.PNUM):
        f.write("'err-ent%d.dat' u (hist($1,width)):(1.0) smooth freq w boxes ls %d title '%6.2lf'," % (i,i+1,setup_g.PGPA[i]))
    f.write("\n")
    f.write("unset multiplot\n")
    f.close()
    os.system("gnuplot < sc-gp")
    removes_fd("sc-gp")
