#================================================================#
#                                                                #
#                       maise-net  TOOLBOX                       #
#                                                                #
#                    Python interface to maise                   #
#                  for automated data generation                 #
#                                                                #
#           Main routines for data generation and test           #
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
from copy import deepcopy
from maisenet_defs  import *
from maisenet_task  import *
from maisenet_util  import *
from maisenet_inot  import *

#----------------------------------------------------------------------------------

def DATGEN_REF(RUN, output, setup_0): # REF DATA GENERATION

     if setup_0.NSPC > 1:
          return

     # type of the data which is being generated
     runt    = dref_.upper()

     # local copy of setup
     l_setup = deepcopy(setup_0)

     cwd = os.getcwd()
     fileout = cwd+"/"+output
     export_out(fileout,"\n                      ===== running in %s mode =====\n" % runt,time = 0,color = 3)

     # default values of NDIM=0, KDNS=0, ISIF=2,PSTRESS=0.0 are enforced
     l_setup.SIGM          = defsigma
     l_setup.NDIM          = 0
     l_setup.KDNS          = 0
     l_setup.ISIF          = 2
     l_setup.PGPA[0]       = 0.0
     datdir         = cwd+"/"+mdat_+"/"+dref_+"/"  # path to DFT runs
     pos            = poscar()
     pos.NATM       = 1
     pos.NTYP       = 1
     pos.NAET       = [1]
     pos.SPCS       = l_setup.SATM.split()
     pos.POSI[0]    = [0.5,0.5,0.5]
     pos.COOR        = 1
     if l_setup.LSYM == 1:
          pos.LATT = [[l_setup.LBOX,0,0],[0,l_setup.LBOX,0],[0,0,l_setup.LBOX]]
     else:
          pos.LATT = [[l_setup.LBOX*0.95,0,0],[0,l_setup.LBOX,0],[0,0,l_setup.LBOX*1.05]]

     if l_setup.ERNI == 1 and does__exst(datdir):
          export_out(fileout,"           Skipping... % s" % datdir,time = 0,color = 3)
          return

     if l_setup.LBOX == 0.0:
          export_out(fileout,"Error: LBOX value can not be zero for reference structures",color = 2);exit()
     if does__exst(datdir):
          export_out(fileout,"Error: % s directory already exists" % datdir,color = 2);exit()

     # ===== starting the task

     reread_stp(fileout,l_setup,cwd+"/setup")

     export_out(fileout,"Note: generating %s dataset" % runt ,color = 3)

     export_out(fileout,"Note: a small SIGMA of % 5.3lf will be used for the %s data generation" % (l_setup.SIGM,runt),color = 3)

     # submit DFT runs
     run = datdir+"00/"
     mkdir_tree(run)
     export_str(pos,run+"/POSCAR")
     (N,direc,fies) = searchfile(datdir,"POSCAR")
     export_out(fileout,"Note: submitting % 5d DFT runs for %s dataset at % s..." % (N,runt,datdir),color = 3)
     submit_dft(cwd,mini_,datdir,fileout,"jdft",l_setup)

     for j in range(1,rnumb):
          tag = format(j,"02d")
          duplicates(datdir+"00",datdir+tag)
     os.chdir(cwd)

     duplicates(cwd+"/"+mini_+"/"+iprs_+"/tag-ref",datdir+"/tag")

     (N,direc,files) = searchfile(datdir+dref_+"*","POSCAR.0")
     export_out(fileout,"Note: %s data generation is all done! total % 5d DFT structure are generated in % s" % (runt,N,datdir),color = 3)

#----------------------------------------------------------------------------------

def DATGEN_EOS(RUN, output, cyc, setup_0): # EOS DATA GENERATION

     SORT = setup_0.SORT

     # type of the data which is being generated
     if RUN != JBAS:
          runt = daes_.upper()
     if RUN == JBAS:
          runt = deos_.upper()

     # local copy of setup
     x_setup = deepcopy(setup_0)

     cwd = os.getcwd()
     fileout = cwd+"/"+output

     # list for radi; in case prototypes do not have species
     elelist = [29,46,47,50,79,11,12,20]

     # default values of ISIF=3,PSTRESS=0.0 are enforced
     x_setup.ISIF     = 3
     x_setup.PGPA[0]  = 0.0
     srcdir = cwd+"/"+mini_+"/"+ieos_+"/"       # path to source of the prototypes
     datdir = cwd+"/"+mdat_+"/"                   # path to DFT runs for *EOS
     min_z  = eminz
     max_z  = emaxz
     if RUN != JBAS:
          eosdir = cwd+"/"+meos_+"/"+daes_+format(cyc,"02d")+"/"  # path to the optimized prototypes
     if RUN == JBAS:
          eosdir = cwd+"/"+meos_+"/"+deos_+"00"+"/"  # path to the optimized prototypes

     if RUN == JBAS:
          export_out(fileout,"\n                      ===== running in %s mode =====\n" % runt,time = 0,color = 3)

     export_out(fileout,"Note: generating %s dataset" % runt,color = 3)

     # ===== starting the task

     # skip the DFT runs for prototypes for EVOS
     if RUN != JBAS:
          num    = 4
          (N,direc,files) = searchfile(eosdir,"POSCAR"+str(x_setup.NSPC)+"*")

     # prepare and submit DFT runs for prototypes for BASICS
     if RUN == JBAS:
          num    = enumb
          export_out(fileout,"Note: optimizing %s prototypes ..." % runt,color = 3)
          reread_stp(fileout,x_setup,cwd+"/setup")

          (N,direc,files) = searchfile(srcdir,"POSCAR"+str(x_setup.NSPC)+"*")
          if N == 0:
               export_out(fileout,"Error: no POSCAR prototype for this system exists in % s" % srcdir,color = 2);exit()

          for i in range(0,N):
               tag = files[i].replace("POSCAR"+str(x_setup.NSPC),"")
               run = eosdir+tag
               if does__exst(run) and x_setup.ERNI == 0:
                    export_out(fileout,"Error: % s directory already exists" % run,color = 2);exit()
               if does__exst(run) and x_setup.ERNI == 1:
                    export_out(fileout,"           Skipping... % s" % run,time = 0,color = 3)
                    continue
               mkdir_tree(run)
               duplicates(direc[i]+"/"+files[i],run+"/POSCAR")

               # find the atomic radius for species in prototype
               s = file_nthln(run+"/POSCAR",0)
               radi = 0.0
               if len(s.split()) >= x_setup.NSPC: # if mentioned; read them!
                    c = 0
                    for j in range(0,x_setup.NSPC):
                         if atomic_rad(symbtoatom(s.split()[j])) > 0.0:
                              c += 1
                              radi += atomic_rad(symbtoatom(s.split()[j]))
                    if c == x_setup.NSPC:
                         radi /= float(c)
               if radi == 0.0: # if not mentioned; take a guess!
                    for j in range(0,x_setup.NSPC):
                         radi = atomic_rad(elelist[j])/float(x_setup.NSPC)

               radi = x_setup.RADI / radi
               replacestr(run+"/POSCAR","VVVV",str(radi))
               replacestr(run+"/POSCAR","TTTT",x_setup.SATM)
               duplicates(cwd+"/"+mini_+"/INCAR", run)
               replacestr(run+"/INCAR","PPPP",str(0.0))
               replacestr(run+"/INCAR","AAAA",str(x_setup.PREC))
               replacestr(run+"/INCAR","GGGG",str(x_setup.ISIF))
               replacestr(run+"/INCAR","SMSM",str(x_setup.SMER))
               replacestr(run+"/INCAR","SGSG",str(x_setup.SIGM))
               if x_setup.SYMP != 0.0:
                    addstrfile(run+"/INCAR","SYMPREC="+str(x_setup.SYMP))
               if x_setup.SORB > 0:
                    addstrfile(run+"/INCAR","LSORBIT=.TRUE.")
                    addstrfile(run+"/INCAR","LMAXMIX="+str(x_setup.SORB))
               if x_setup.ECUT > 0.0:
                    l = "ENCUT="+str(x_setup.ECUT)
                    addstrfile(run+"/INCAR",l)
               pos = import_str(run+"/POSCAR")
               if tag[0:2] == "lr": # for 2D eos data; only x-y changes
                    dim = 2
                    # KPOINTS and INCAR should be copied here
                    eliminates(run+"/INCAR","SYMPREC")
                    smpl_kmesh(pos,8000,path = run, ndim = 2)

          submit_dft(cwd,mini_,eosdir,fileout,"jeos",x_setup)

     # collect prototypes, generate EOS structures, and submit DFT runs
     for i in range(0,N):
          reread_stp(fileout,x_setup,cwd+"/setup")
          tag  = files[i].replace("POSCAR"+str(x_setup.NSPC),"")
          if RUN != JBAS:
               cle =  tag[0:2]
               prs =  tag[2:4]
               ids =  tag[4:6]
               siz =  tag[6:len(tag)]

          if doesntexst(eosdir+tag+"/CONTCAR.0"):
               export_out(fileout,"Error: the file %s does not exist" % (eosdir+tag+"/CONTCAR.0"),color = 2);exit()
          if (SORT == 1 or SORT == 2 or SORT == 3) and does__exst(datdir+deos_+tag):
               if x_setup.ERNI == 0:
                    export_out(fileout,"Error: % s directory already exists" % (datdir+deos_+tag),color = 2);exit()
               else:
                    export_out(fileout,"           Skipping... % s" % (datdir+deos_+tag),time = 0,color = 3)
                    continue
          if SORT == 0 and does__exst(datdir+deos_+"/"+tag):
               if x_setup.ERNI == 0:
                    export_out(fileout,"Error: % s directory already exists" % (datdir+deos_+"/"+tag),color = 2);exit()
               else:
                    export_out(fileout,"           Skipping... % s" % (datdir+deos_+"/"+tag),time = 0,color = 3)
                    continue
          for j in range(0,num+1):
               lbl = format(j,"02d")
               if RUN != JBAS:
                    if SORT == 1:
                         run = datdir+daes_+x_setup.EHSH[int(prs)][int(ids)]+"/"+str(cle)+"/"+str(siz)+"/"+lbl
                    if SORT == 3:
                         run = datdir+daes_+str(prs)+str(ids)+"/"+str(cle)+"/"+str(siz)+"/"+lbl
                    if SORT == 2 or SORT == 0:
                         run = datdir+daes_+"/"+str(prs)+"/"+str(cle)+"/"+str(ids)+"/"+str(siz)+"/"+lbl
               if RUN == JBAS:
                    if SORT == 1 or SORT == 2 or SORT == 3:
                         run = datdir+deos_+tag+"/"+lbl
                    if SORT == 0:
                         run = datdir+deos_+"/"+tag+"/"+lbl

               mkdir_tree(run)
               scl = min_z+float(j)*float(max_z-min_z)/float(num)
               duplicates(cwd+"/"+mini_+"/INCAR", run)
               replacestr(run+"/INCAR","PPPP",str(0.0))
               replacestr(run+"/INCAR","AAAA",str(x_setup.PREC))
               replacestr(run+"/INCAR","GGGG",str(x_setup.ISIF))
               replacestr(run+"/INCAR","SMSM",str(x_setup.SMER))
               replacestr(run+"/INCAR","SGSG",str(x_setup.SIGM))
               if x_setup.SYMP != 0.0:
                    addstrfile(run+"/INCAR","SYMPREC="+str(x_setup.SYMP))
               if x_setup.SORB > 0:
                    addstrfile(run+"/INCAR","LSORBIT=.TRUE.")
                    addstrfile(run+"/INCAR","LMAXMIX="+str(x_setup.SORB))
               if x_setup.ECUT > 0.0:
                    l = "ENCUT="+str(x_setup.ECUT)
                    addstrfile(run+"/INCAR",l)
               pos = import_str(eosdir+tag+"/CONTCAR.0")

               if tag[0:2] == "lr": # for 2D eos data; only x-y changes
                    dim = 2
                    # KPOINTS and INCAR should be copied here
                    eliminates(run+"/INCAR","SYMPREC")
                    smpl_kmesh(pos,8000,path = run, ndim = 2)
               else:
                    dim = 3

               for x in range(0,dim):
                    for y in range(0,3):
                         pos.LATT[x][y] *= (scl)

               c = 0
               if scl <  0.995 or scl > 1.005:
                    export_str(pos,run+"/POSCAR")
               else:
                    if RUN == JBAS and c == 0:
                         duplicates(eosdir+tag+"/POSCAR.0",run)
                         duplicates(eosdir+tag+"/CONTCAR.0",run)
                         duplicates(eosdir+tag+"/OUTCAR.0",run)
                         duplicates(eosdir+tag+"/dat.dat",run)
                         c = c + 1
                    if RUN != JBAS:
                         export_str(pos,run+"/poscar")

          if RUN != JBAS:
               duplicates(cwd+"/"+mini_+"/"+iprs_+"/tag-aes",datdir+daes_+x_setup.EHSH[int(prs)][int(ids)]+"/tag")
               if SORT == 1 and doesntexst(datdir+daes_+x_setup.EHSH[int(prs)][int(ids)]+"/inf"):
                    s = ""
                    for j in range(0,x_setup.NSPC):
                         s += str(x_setup.SPCN[j][int(ids)])+" "
                    addstrfile(datdir+daes_+x_setup.EHSH[int(prs)][int(ids)]+"/inf","pressure % s" % (str(x_setup.PGPA[int(prs)])))
                    addstrfile(datdir+daes_+x_setup.EHSH[int(prs)][int(ids)]+"/inf","composit % s" % (s))
               if SORT == 3:
                    duplicates(cwd+"/"+mini_+"/"+iprs_+"/tag-aes",datdir+daes_+str(prs)+str(ids)+"/tag")
               if SORT == 2 or SORT == 0:
                    duplicates(cwd+"/"+mini_+"/"+iprs_+"/tag-aes",datdir+daes_+"/tag")
          if RUN == JBAS:
               if SORT == 1 or SORT == 2 or SORT == 3:
                    duplicates(cwd+"/"+mini_+"/"+iprs_+"/tag-eos",datdir+deos_+tag+"/tag")
               if SORT == 0:
                    duplicates(cwd+"/"+mini_+"/"+iprs_+"/tag-eos",datdir+deos_+"/tag")

     if RUN != JBAS:
          (N,direc,fies) = searchfile(datdir+daes_+"*","POSCAR")
          export_out(fileout,"Note: submitting % 5d DFT runs for %s dataset at % s..." % (N,runt,datdir+daes_+"*"),color = 3)
          if N > 0:
               N = submit_dft(cwd,mini_,datdir+daes_+"*",fileout,"jdft",x_setup)

     if RUN == JBAS:
          (N,direc,fies) = searchfile(datdir+deos_+"*","POSCAR")
          export_out(fileout,"Note: submitting % 5d DFT runs for %s dataset at % s..." % (N,runt,datdir+deos_+"*"),color = 3)
          if N > 0:
               N = submit_dft(cwd,mini_,datdir+deos_+"*",fileout,"jdft",x_setup)
          (N,direc,fies) = searchfile(datdir+deos_+"*","POSCAR.0")
          export_out(fileout,"Note: %s data generation is all done! total % 5d DFT structure are generated in % s" % (runt,N,datdir+deos_+"*"),color = 3)

#----------------------------------------------------------------------------------

def DATGEN_CLS(RUN, output, setup_0): # CLST DATA GENERATION

     if setup_0.NSPC > 3:
          return

     # type of the data which is being generated
     runt = dcls_.upper()

     # flag for sorting the training data
     SORT = setup_0.SORT

     # local copy of setup
     z_setup=deepcopy(setup_0)

     cwd = os.getcwd()
     fileout = cwd+"/"+output
     export_out(fileout,"\n                      ===== running in %s mode =====\n" % runt,time = 0,color = 3)

     # default values of KDNS=0,ISIF=2,PSTRESS=0.0 are enforced
     z_setup.KDNS    = 0
     z_setup.ISIF    = 2
     z_setup.PGPA[0] = 0.0
     z_setup.NDIM    = 0
     z_setup.SIGM    = defsigma
     datdir = cwd+"/"+mdat_+"/"           # path to DFT runs
     min_z  = z_setup.RADI*cminz
     max_z  = min_z*cmaxz
     num    = cnumb
     delta  = (max_z-min_z)/float(num)

     if z_setup.LBOX == 0.0:
          export_out(fileout,"Error: LBOX value can not be zero for clusters",color = 2);exit()

     # ===== starting the task

     export_out(fileout,"Note: generating %s dataset" % runt,color = 3)
     export_out(fileout,"Note: a small SIGMA of % 5.3lf will be used for the %s data generation" % (z_setup.SIGM,runt),color = 3)


     if SORT == 1 or SORT == 3 or SORT == 0:
          if does__exst(datdir+dcls_):
               if z_setup.ERNI == 0:
                    export_out(fileout,"Error: % s directory already exists" % (datdir+dcls_),color = 2);exit()
               else:
                    export_out(fileout,"           Skipping... % s" % (datdir+dcls_),time = 0,color = 3)
                    return

     if SORT == 2:
          if z_setup.NSPC < 3:
               if does__exst(datdir+dcls_+"02") or does__exst(datdir+dcls_+"04"):
                    if z_setup.ERNI == 0:
                         export_out(fileout,"Error: % s directory already exists" % (datdir+dcls_+"0*"),color = 2);exit()
                    else:
                         export_out(fileout,"           Skipping... % s" % (datdir+dcls_+"0*"),time = 0,color = 3)
                         return
          if z_setup.NSPC == 3:
               if does__exst(datdir+dcls_+"03"):
                    if z_setup.ERNI == 0:
                         export_out(fileout,"Error: % s directory already exists" % (datdir+dcls_+"0*"),color = 2);exit()
                    else:
                      export_out(fileout,"           Skipping... % s" % (datdir+dcls_+"0*"),time = 0,color = 3)
                      return

     # prepare CLST structures of 2- and 4- atoms, and submit DFT runs
     for i in range(0,num+1):
          reread_stp(fileout,z_setup,cwd+"/setup")
          tag = format(i,"02d")
          pos            = poscar()
          pos.COOR       = 1
          pos.NTYP       = z_setup.NSPC
          pos.LATT       = [[z_setup.LBOX,0,0],[0,z_setup.LBOX,0],[0,0,z_setup.LBOX]]
          comp           = (min_z+delta*float(i))/float(z_setup.LBOX)
          pos.SPCS       = z_setup.SATM.split()

          if z_setup.NSPC < 3:
               # 2-atom structures
               pos.NATM       = 2
               for j in range(0,z_setup.NSPC):
                    pos.NAET.append(pos.NATM/z_setup.NSPC)
               pos.POSI[0]    = [0.5,0.5,0.5]
               pos.POSI[1]    = [0.5,0.5,0.5+comp]
               if SORT == 1 or SORT == 3 or SORT == 0:
                    run = datdir+dcls_+"/02/"+tag
               if SORT == 2:
                    run = datdir+dcls_+"02/"+tag
               mkdir_tree(run)
               export_str(pos,run+"/POSCAR")
               # 4-atom structures
               pos.NATM       = 4
               for j in range(0,z_setup.NSPC):
                    pos.NAET.append(pos.NATM/z_setup.NSPC)
               pos.POSI[0]    = [0.5,0.5,0.5]
               pos.POSI[1]    = [0.5,0.5,0.5+comp]
               pos.POSI[2]    = [0.5,0.5+comp,0.5]
               pos.POSI[3]    = [0.5,0.5+comp,0.5+comp]
               if SORT == 1 or SORT == 3 or SORT == 0:
                    run = datdir+dcls_+"/04/"+tag
               if SORT == 2:
                    run = datdir+dcls_+"04/"+tag
               mkdir_tree(run)
               export_str(pos,run+"/POSCAR")
          if z_setup.NSPC == 3:
               # 3-atom structures
               pos.NATM       = 3
               for j in range(0,z_setup.NSPC):
                    pos.NAET.append(pos.NATM/z_setup.NSPC)
               pos.POSI[0]    = [0.5,0.5,0.5]
               pos.POSI[1]    = [0.5,0.5,0.5+comp]
               pos.POSI[2]    = [0.5,0.5+comp,0.5]
               if SORT == 1 or SORT == 3 or SORT == 0:
                    run = datdir+dcls_+"/03/"+tag
               if SORT == 2:
                    run = datdir+dcls_+"03/"+tag
               mkdir_tree(run)
               export_str(pos,run+"/POSCAR")

     if SORT == 1 or SORT == 3 or SORT == 0:
          duplicates(cwd+"/"+mini_+"/"+iprs_+"/tag-cls",datdir+dcls_+"/tag")
     if SORT == 2:
          if z_setup.NSPC < 3:
               duplicates(cwd+"/"+mini_+"/"+iprs_+"/tag-cls",datdir+dcls_+"02/tag")
               duplicates(cwd+"/"+mini_+"/"+iprs_+"/tag-cls",datdir+dcls_+"04/tag")
          if z_setup.NSPC == 3:
               duplicates(cwd+"/"+mini_+"/"+iprs_+"/tag-cls",datdir+dcls_+"03/tag")

     (N,direc,fies) = searchfile(datdir+dcls_+"*","POSCAR")
     export_out(fileout,"Note: submitting % 5d DFT runs for %s dataset at % s..." % (N,runt,datdir+dcls_+"*"),color = 3)
     submit_dft(cwd,mini_,datdir+dcls_+"*",fileout,"jdft",z_setup)

     (N,direc,files) = searchfile(datdir+dcls_+"*","POSCAR.0")
     export_out(fileout,"Note: %s data generation is all done! total % 5d DFT structure are generated in % s" % (runt,N,datdir+dcls_+"*"),color = 3)

#----------------------------------------------------------------------------------

def DATGEN_EVO(RUN, output, cyc, setup_0): # EVOS DATA GENERATION

     SORT = setup_0.SORT

     # type of the data which is being generated
     runt = devo_.upper()

     # local copy of setup
     m_setup = deepcopy(setup_0)

     CYCLES = format(cyc,"02d")

     cwd = os.getcwd()
     fileout = cwd+"/"+output
     export_out(fileout,"\n                      ===== running in %s mode =====\n" % runt,time = 0,color = 3)

     datdir = cwd+"/"+mdat_+"/"             # path to collected DFT data of the cycle for parsing
     srcdir = cwd+"/"+mevo_+"/"+CYCLES+"/"  # path to NN evolution runs
     prsdir = cwd+"/"+mprs_+"/"+CYCLES+"/"  # path to parsed data
     libdir = cwd+"/"+mlib_+"/"+CYCLES+"/"  # path to final NN model
     mindir = cwd+"/"+mmin_+"/"             # path to unique minima of NN runs
     model  = cwd+"/model"

     if does__exst(srcdir):
          export_out(fileout,"Error: % s already exists" % srcdir,color = 2);exit()
     if does__exst(prsdir):
          export_out(fileout,"Error: % s already exists" % prsdir,color = 2);exit()
     if does__exst(libdir):
          export_out(fileout,"Error: % s already exists" % libdir,color = 2);exit()

     # ===== starting the task

     # ===== CYCLE z DATA GENERATION                     
     if cyc == 0:

          if does__exst(mindir):
               export_out(fileout,"Error: % s already exists" % mindir,color = 2);exit()

          mkdir_tree(mindir)

          # prepare and submit ES runs for CYCLE 0 dataset
          export_out(fileout,"Note: generating %s dataset in cycle = % 3d" % (runt,cyc),color = 3)

          permin = []
          for i in range(0,m_setup.PNUM):
               permin.append(0)
          totmin  = 0

          for z in range(0,m_setup.PNUM):
               PRESSU = format(z,"02d")
               for i in range(0,len(m_setup.SIZ0)):
                    reread_stp(fileout,m_setup,cwd+"/setup")
                    POPULN  = format(i,"02d")
                    SRCDIR  = srcdir+"/"+PRESSU+"/"+POPULN  # always indexed as: EVO/cycle/pressure/population
                    if SORT == 1:
                         DATDIR = datdir+devo_+m_setup.ehsh[z][i]+"/"+CYCLES   # indexed using hash
                    if SORT == 3:
                         DATDIR = datdir+devo_+PRESSU+POPULN+"/"+CYCLES # indexed as: DAT/"evo"pressure"population"/cycle/
                    if SORT == 2:
                         DATDIR = datdir+devo_+PRESSU+"/"+CYCLES # indexed as: DAT/"evo"pressure/cycle/population
                    if SORT == 0:
                         DATDIR = datdir+devo_+"/"+PRESSU+"/"+CYCLES # indexed as: DAT/"evo/"pressure/cycle/population
                    if does__exst(DATDIR):
                         export_out(fileout,"Error: % s already exists" % DATDIR,color = 2);exit()
                    popsize = int(round(m_setup.npop[i]*m_setup.PWGT[z]))
                    pop     = submit_esj(cwd,mini_,DATDIR,mindir,SRCDIR,fileout,iesz_,i,z,cyc,m_setup)
                    totmin    += pop
                    permin[z] += pop

          export_out(fileout,"Note: collected % 5d ES-minima in cycle = % 3d" % (totmin,cyc),color = 3)
          for i in range(0,m_setup.PNUM):
               export_out(fileout,"        pressure % 6.2lf  ES-minima % 5d (cycle = % 3d)" %(m_setup.PGPA[i],permin[i],cyc),time = 0)

          # collect minima and finished DFT runs in DAT directory for cyc=0
          (NM,mdi,mfi) = searchfile(mindir,"P"+CYCLES+"*")
          for i in range(0,NM):
               fi  = mfi[i]
               C   = fi[1:3]
               P   = fi[3:5]
               R   = fi[5:7]
               I   = fi[7:10]
               sr  = file_nthln(mdi[i]+"/"+mfi[i],0).split()[0]
               if SORT == 1:
                    run = datdir+devo_+m_setup.ehsh[int(P)][int(R)]+"/"+C+"/"+I+"000"   # indexed using hash
                    tar = datdir+devo_+m_setup.ehsh[int(P)][int(R)]+"/tag"
               if SORT == 3:
                    run = datdir+devo_+P+R+"/"+C+"/"+I+"000"  # indexed as: DAT/"evo"pressure"populations"/cycle
                    tar = datdir+devo_+P+R+"/tag"
               if SORT == 2:
                    run = datdir+devo_+P+"/"+C+"/"+R+"/"+I+"000"  # indexed as: DAT/"evo"pressure/cycle/population
                    tar = datdir+devo_+P+"/tag"
               if SORT == 0:
                    run = datdir+devo_+"/"+P+"/"+C+"/"+R+"/"+I+"000"  # indexed as: DAT/"evo/"pressure/cycle/population
                    tar = datdir+devo_+"/tag"
               mkdir_tree(run)
               if SORT == 1 and doesntexst(datdir+devo_+m_setup.ehsh[int(P)][int(R)]+"/inf"):
                    s = ""
                    for j in range(0,m_setup.NSPC):
                         s += str(m_setup.SPC0[j][int(R)])+" "
                    addstrfile(datdir+devo_+m_setup.ehsh[int(P)][int(R)]+"/inf","pressure % s" % (str(m_setup.PGPA[int(P)])))
                    addstrfile(datdir+devo_+m_setup.ehsh[int(P)][int(R)]+"/inf","composit % s" % (s))
               duplicates(cwd+"/"+mini_+"/"+iprs_+"/tag-evo",tar)
               duplicates(mindir+"/"+fi,run+"/POSCAR.0") # this is equal to: duplicates(sr+"/CONTCAR.0",run+"/POSCAR.0")
               duplicates(sr+"/OUTCAR.0",run+"/OUTCAR.0")
               datdat_out(fileout,sr+"/OUTCAR.0",run)

          # parse the data for cyc=0
          submit_prs(cwd,mini_,datdir,prsdir,fileout,cyc,m_setup)

          # train the model for cyc=0
          submit_trn(cwd,mini_,libdir,prsdir,fileout,"jtrn",cyc,m_setup)

          if SORT == 1:
               (N,direc,files) = searchfile(datdir+devo_+"*/"+CYCLES,"POSCAR.0")
               export_out(fileout,"Note: %s data generation in cycle % d is all done! total % 5d DFT structure are generated in % s" % (runt,cyc,N,datdir+devo_+"*/"+CYCLES),color = 3)
          if SORT == 3:
               (N,direc,files) = searchfile(datdir+devo_+"*/"+CYCLES,"POSCAR.0")
               export_out(fileout,"Note: %s data generation in cycle % d is all done! total % 5d DFT structure are generated in % s" % (runt,cyc,N,datdir+devo_+"*/"+CYCLES),color = 3)
          if SORT == 2:
               (N,direc,files) = searchfile(datdir+devo_+"*/"+CYCLES,"POSCAR.0")
               export_out(fileout,"Note: %s data generation in cycle % d is all done! total % 5d DFT structure are generated in % s" % (runt,cyc,N,datdir+devo_+"*/"+CYCLES),color = 3)
          if SORT == 0:
               (N,direc,files) = searchfile(datdir+devo_+"/*/"+CYCLES,"POSCAR.0")
               export_out(fileout,"Note: %s data generation in cycle % d is all done! total % 5d DFT structure are generated in % s" % (runt,cyc,N,datdir+devo_+"/*/"+CYCLES),color = 3)

     # ===== CYCLE n DATA GENERATION                     
     if cyc > 0:
          N1 = 0
          N2 = 0
          N3 = 0
          if doesntexst(model):
               if does__exst(cwd+"/"+mlib_+"/"+format(cyc-1, "02d")+"/model"):
                    duplicates(cwd+"/"+mlib_+"/"+format(cyc-1, "02d")+"/model",cwd+"/model")
               else:
                    export_out(fileout,"Error: the neural network model does not exist at % s or % s " % (model,mlib_+"/"+format(cyc-1,"02d")),color = 2);exit()

          if doesntexst(mindir):
               export_out(fileout,"Error: the %s directory does not exist" % (mindir),color = 2);exit()

          permin = []
          for i in range(0,m_setup.PNUM):
               permin.append(0)
          totmin  = 0

          # prepare and submit ES runs for CYCLE 1+ dataset
          export_out(fileout,"Note: generating %s dataset in cycle = % 3d" % (runt,cyc),color = 3)

          for z in range(0,m_setup.PNUM):
               PRESSU = format(z,"02d")
               for i in range(0,len(m_setup.SIZN)):
                    reread_stp(fileout,m_setup,cwd+"/setup")
                    POPULN  = format(i,"02d")
                    SRCDIR  = srcdir+"/"+PRESSU+"/"+POPULN  # always indexed as: EVO/cycle/pressure/population
                    if SORT == 1:
                         DATDIR = datdir+devo_+m_setup.EHSH[z][i]+"/"+CYCLES # indexed by hash
                    if SORT == 3:
                         DATDIR = datdir+devo_+PRESSU+POPULN+"/"+CYCLES # indexed as: DAT/"evo"pressure"population"/cycle
                    if SORT == 2:
                         DATDIR = datdir+devo_+PRESSU+"/"+CYCLES # indexed as: DAT/"evo"pressure/cycle/population
                    if SORT == 0:
                         DATDIR = datdir+devo_+"/"+PRESSU+"/"+CYCLES # indexed as: DAT/"evo"pressure/cycle/population
                    if does__exst(DATDIR):
                         export_out(fileout,"Error: % s already exists" % DATDIR,color = 2);exit()
                    popsize = int(round(m_setup.NPOP[i]*m_setup.PWGT[z]))
                    pop     = submit_esj(cwd,mini_,DATDIR,mindir,SRCDIR,fileout,iesn_,i,z,cyc,m_setup)
                    totmin    += pop
                    permin[z] += pop

          export_out(fileout,"Note: collected % 5d ES-minima in cycle = % 3d" % (totmin,cyc),color = 3)
          for i in range(0,m_setup.PNUM):
               export_out(fileout,"        pressure % 6.2lf  ES-minima % 5d (cycle = % 3d)" %(m_setup.PGPA[i],permin[i],cyc),time = 0)

          # collect the minima and submit DFT runs in DAT directory for CYCLES 1+
          (NM,mdi,mfi) = searchfile(mindir,"P"+CYCLES+"*")
          for i in range(0,NM):
               fi  = mfi[i]
               C   = fi[1:3]
               P   = fi[3:5]
               R   = fi[5:7]
               I   = fi[7:10]
               sr  = file_nthln(mdi[i]+"/"+mfi[i],0).split()[0]
               if SORT == 1:
                    run = datdir+devo_+m_setup.EHSH[int(P)][int(R)]+"/"+C+"/"+I+"000"  # indexed as: DAT/"evo"pressure"populations"/cycle
                    tar = datdir+devo_+m_setup.EHSH[int(P)][int(R)]+"/tag"
               if SORT == 3:
                    run = datdir+devo_+P+R+"/"+C+"/"+I+"000"  # indexed as: DAT/"evo"pressure"populations"/cycle
                    tar = datdir+devo_+P+R+"/tag"
               if SORT == 2:
                    run = datdir+devo_+P+"/"+C+"/"+R+"/"+I+"000"  # indexed as: DAT/"evo"pressure/cycle/population
                    tar = datdir+devo_+P+"/tag"
               if SORT == 0:
                    run = datdir+devo_+"/"+P+"/"+C+"/"+R+"/"+I+"000"  # indexed as: DAT/"evo/"pressure/cycle/population
                    tar = datdir+devo_+"/tag"
               mkdir_tree(run)
               if SORT == 1 and doesntexst(datdir+devo_+m_setup.EHSH[int(P)][int(R)]+"/inf"):
                    s = ""
                    for j in range(0,m_setup.NSPC):
                         s += str(m_setup.SPCN[j][int(R)])+" "
                    addstrfile(datdir+devo_+m_setup.EHSH[int(P)][int(R)]+"/inf","pressure % s" % (str(m_setup.PGPA[int(P)])))
                    addstrfile(datdir+devo_+m_setup.EHSH[int(P)][int(R)]+"/inf","composit % s" % (s))
               duplicates(cwd+"/"+mini_+"/"+iprs_+"/tag-evo",tar)
               duplicates(mindir+"/"+fi,run+"/POSCAR")
               duplicates(cwd+"/"+mini_+"/INCAR", run)
               replacestr(run+"/INCAR","PPPP",str(m_setup.PGPA[int(P)]*10.0))
               replacestr(run+"/INCAR","AAAA",str(m_setup.PREC))
               replacestr(run+"/INCAR","GGGG",str(m_setup.ISIF))
               replacestr(run+"/INCAR","SMSM",str(m_setup.SMER))
               replacestr(run+"/INCAR","SGSG",str(m_setup.SIGM))
               if m_setup.SYMP != 0.0:
                    addstrfile(run+"/INCAR","SYMPREC="+str(m_setup.SYMP))
               if m_setup.SORB > 0:
                    addstrfile(run+"/INCAR","LSORBIT=.TRUE.")
                    addstrfile(run+"/INCAR","LMAXMIX="+str(m_setup.SORB))
               if m_setup.ECUT > 0.0:
                    l = "ENCUT="+str(m_setup.ECUT)
                    addstrfile(run+"/INCAR",l)

          if SORT == 1:
               export_out(fileout,"Note: submitting % 5d DFT runs for minima of %s dataset at % s..." % (NM,runt,datdir+devo_+"*/"+CYCLES),color = 3)
          if SORT == 3:
               export_out(fileout,"Note: submitting % 5d DFT runs for minima of %s dataset at % s..." % (NM,runt,datdir+devo_+"*/"+CYCLES),color = 3)
          if SORT == 2:
               export_out(fileout,"Note: submitting % 5d DFT runs for minima of %s dataset at % s..." % (NM,runt,datdir+devo_+"*/"+CYCLES),color = 3)
          if SORT == 0:
               export_out(fileout,"Note: submitting % 5d DFT runs for minima of %s dataset at % s..." % (NM,runt,datdir+devo_+"/*/"+CYCLES),color = 3)

          submit_dft(cwd,mini_,datdir+devo_+"*",fileout,"jdft",m_setup)

          # weight for sampling relaxation path; default 1.0 for each minima
          permin = []
          persam = []
          for i in range(0,m_setup.PNUM):
               permin.append(0)
               persam.append(0)

          mwegt = [] 
          (NM,mdi,mfi) = searchfile(cwd+"/"+mmin_,"P"+CYCLES+"*")
          norm = 1.0/float(NM)
          for i in range(0,NM):
               mwegt.append(norm)
          MINIMA = dict(zip(mfi,mwegt)) # list of minima and their weights

          if m_setup.SMLR == 0: # no check is performed
               export_out(fileout,"Note: no similarity check is performed: all % 5d minima are marked as unique" % (minima1+minima2),color = 3)

          if m_setup.SMLR == 1: # DFT relax and assign weights accordingly
              norm = similr_chk(cwd,mini_,mindir,datdir,fileout,cyc,"jdft",NM,MINIMA,m_setup)

          minima = 0
          for i in range(0,NM):
               if MINIMA[list(MINIMA.keys())[i]] > 0.0:
                    minima += 1
          if (minima) >= m_setup.DATA or minima == 0:
               export_out(fileout,"Warning: number of minima is zero or larger than desired number of structures, no relaxation path will be sampled!",time = 0);
               WEGT = 0
          else:
               WEGT = m_setup.DATA - minima
          trajmin = 0
          export_out(fileout,"Note: sampling relaxation path (for target % 5d structures) in cycle = % 3d ..." % (WEGT,cyc),color = 3)
          for i in range(0,NM):
               fi  = list(MINIMA.keys())[i]
               C   = fi[1:3]
               P   = fi[3:5]
               R   = fi[5:7]
               I   = fi[7:10]
               S   = file_nthln(mindir+"/"+fi,0).split()[0]
               wegt    = int(round(WEGT * MINIMA[fi]))
               trajtmp = 0
               if wegt > 0:
                    reread_stp(fileout,m_setup,cwd+"/setup")
                    (o,N,POSCARS) = readoutcar(S+"/OUTCAR.0")
                    delta = int(round((N-1)/(wegt+1)))
                    if wegt >= N or delta == 0: # if there is no enough struc to sample; do all
                         wegt  = N-1
                         delta = 1
                    for k in range(1,wegt+1):
                         if SORT == 1:
                              run = datdir+devo_+m_setup.EHSH[int(P)][int(R)]+"/"+C+"/"+I+format(N-k*delta,"03d")
                         if SORT == 3:
                              run = datdir+devo_+P+R+"/"+C+"/"+I+format(N-k*delta,"03d")
                         if SORT == 2:
                              run = datdir+devo_+P+"/"+C+"/"+R+"/"+I+format(N-k*delta,"03d")
                         if SORT == 0:
                              run = datdir+devo_+"/"+P+"/"+C+"/"+R+"/"+I+format(N-k*delta,"03d")
                         mkdir_tree(run)
                         export_str(POSCARS[N-k*delta],run+"/POSCAR")
                         duplicates(cwd+"/"+mini_+"/INCAR", run)
                         replacestr(run+"/INCAR","PPPP",str(m_setup.PGPA[int(P)]*10.0))
                         replacestr(run+"/INCAR","AAAA",str(m_setup.PREC))
                         replacestr(run+"/INCAR","GGGG",str(m_setup.ISIF))
                         replacestr(run+"/INCAR","SMSM",str(m_setup.SMER))
                         replacestr(run+"/INCAR","SGSG",str(m_setup.SIGM))
                         if m_setup.SYMP != 0.0:
                              addstrfile(run+"/INCAR","SYMPREC="+str(m_setup.SYMP))
                         if m_setup.SORB > 0:
                              addstrfile(run+"/INCAR","LSORBIT=.TRUE.")
                              addstrfile(run+"/INCAR","LMAXMIX="+str(m_setup.SORB))
                         if m_setup.ECUT > 0.0:
                              l = "ENCUT="+str(m_setup.ECUT)
                              addstrfile(run+"/INCAR",l)
                         trajtmp += 1
                         trajmin += 1
                    permin[int(fi[3:5])] += 1
                    persam[int(fi[3:5])] += trajtmp
               export_out(fileout,"          sampled % 5d structures from %s weight % 5.3lf (%s)" % (trajtmp,fi,MINIMA[fi],S),time = 0)

          export_out(fileout,"Note: collected % 5d relaxation path structures from % 5d unique-minima in cycle = % 3d" % (trajmin,minima,cyc),color = 3)
          for i in range(0,m_setup.PNUM):
               export_out(fileout,"        pressure % 6.2lf  unique-minima % 5d  sampled % 5d (cycle = % 3d)" %(m_setup.PGPA[i],permin[i],persam[i],cyc),time = 0)

          os.chdir(cwd)

          # submit DFT runs
          wait = m_setup.WAIT
          if m_setup.WAIT == 1:
               m_setup.WAIT = 9 # submit_dft will do the "wait" only if setup.WAIT == 9 (but for user it should be set to 1)!
          if SORT == 1:
               (N,direc,fies) = searchfile(datdir+devo_+"*/"+CYCLES,"POSCAR")
               export_out(fileout,"Note: submitting % 5d DFT runs for relaxation path of %s dataset at % s..." % (N,runt,datdir+devo_+"*/"+CYCLES),color = 3)
               submit_dft(cwd,mini_,datdir+devo_+"*/"+CYCLES,fileout,"jdft",m_setup)
               (N1,direc,files) = searchfile(datdir+devo_+"*/"+CYCLES,"POSCAR.0")
               export_out(fileout,"Note: %s data generation in cycle % d is all done! % 5d DFT structures are generated in % s" % (runt,cyc,N1,datdir+devo_+"*/"+CYCLES),color = 3)
          if SORT == 3:
               (N,direc,fies) = searchfile(datdir+devo_+"*/"+CYCLES,"POSCAR")
               export_out(fileout,"Note: submitting % 5d DFT runs for relaxation path of %s dataset at % s..." % (N,runt,datdir+devo_+"*/"+CYCLES),color = 3)
               submit_dft(cwd,mini_,datdir+devo_+"*/"+CYCLES,fileout,"jdft",m_setup)
               (N1,direc,files) = searchfile(datdir+devo_+"*/"+CYCLES,"POSCAR.0")
               export_out(fileout,"Note: %s data generation in cycle % d is all done! % 5d DFT structures are generated in % s" % (runt,cyc,N1,datdir+devo_+"*/"+CYCLES),color = 3)
          if SORT == 2:
               (N,direc,fies) = searchfile(datdir+devo_+"*/"+CYCLES,"POSCAR")
               export_out(fileout,"Note: submitting % 5d DFT runs for relaxation path of %s dataset at % s..." % (N,runt,datdir+devo_+"*/"+CYCLES),color = 3)
               submit_dft(cwd,mini_,datdir+devo_+"*/"+CYCLES,fileout,"jdft",m_setup)
               (N1,direc,files) = searchfile(datdir+devo_+"*/"+CYCLES,"POSCAR.0")
               export_out(fileout,"Note: %s data generation in cycle % d is all done! % 5d DFT structures are generated in % s" % (runt,cyc,N1,datdir+devo_+"*/"+CYCLES),color = 3)
          if SORT == 0:
               (N,direc,fies) = searchfile(datdir+devo_+"/*/"+CYCLES,"POSCAR")
               export_out(fileout,"Note: submitting % 5d DFT runs for relaxation path of %s dataset at % s..." % (N,runt,datdir+devo_+"/*/"+CYCLES),color = 3)
               submit_dft(cwd,mini_,datdir+devo_+"/*/"+CYCLES,fileout,"jdft",m_setup)
               (N1,direc,files) = searchfile(datdir+devo_+"/*/"+CYCLES,"POSCAR.0")
               export_out(fileout,"Note: %s data generation in cycle % d is all done! % 5d DFT structures are generated in % s" % (runt,cyc,N1,datdir+devo_+"/*/"+CYCLES),color = 3)
          m_setup.WAIT = wait
          # submit AEOS runs for unique minima of the current cycle (if any)
          if m_setup.AEOS == 1:
               (NM,mdi,mfi) = searchfile(cwd+"/"+mmin_,"P"+CYCLES+"*")
               for j in range(0,NM):
                    reread_stp(fileout,m_setup,cwd+"/setup")
                    fi   = list(MINIMA.keys())[j]
                    C    = fi[1:3]
                    P    = fi[3:5]
                    R    = fi[5:7]
                    I    = fi[7:10]
                    if SORT == 1:
                         run = datdir+devo_+m_setup.EHSH[int(P)][int(R)]+"/"+C+"/"+I+"00r"
                    if SORT == 3:
                         run = datdir+devo_+P+R+"/"+C+"/"+I+"00r"
                    if SORT == 2:
                         run = datdir+devo_+P+"/"+C+"/"+R+"/"+I+"00r"
                    if SORT == 0:
                         run = datdir+devo_+"/"+P+"/"+C+"/"+R+"/"+I+"00r"
                    if abs(norm-MINIMA[fi]) < 1.0e-3:
                         eosdir = cwd+"/"+meos_+"/"+daes_+format(cyc,"02d")+"/"+C+P+R+I+"/"
                         mkdir_tree(eosdir)
                         duplicates(run+"/CONTCAR.0",eosdir+"/POSCAR"+str(m_setup.NSPC)+C+P+R+I)
                         duplicates(run+"/CONTCAR.0",eosdir)
                         duplicates(run+"/POSCAR.0", eosdir)
                         duplicates(run+"/OUTCAR.0", eosdir)
                         duplicates(run+"/dat.dat",  eosdir)
               DATGEN_EOS(RUN,output,cyc,m_setup)
               if SORT == 1:
                    (N2,direc,files) = searchfile(datdir+daes_+"*/"+CYCLES,"POSCAR.0")
                    export_out(fileout,"Note: %s data generation in cycle % d is all done! % 5d DFT structures are generated in % s" % (daes_.upper(),cyc,N2,datdir+daes_+"*/"+CYCLES),color = 3)
               if SORT == 3:
                    (N2,direc,files) = searchfile(datdir+daes_+"*/"+CYCLES,"POSCAR.0")
                    export_out(fileout,"Note: %s data generation in cycle % d is all done! % 5d DFT structures are generated in % s" % (daes_.upper(),cyc,N2,datdir+daes_+"*/"+CYCLES),color = 3)
               if SORT == 2 or SORT == 0:
                    (N2,direc,files) = searchfile(datdir+daes_+"/*/"+CYCLES,"POSCAR.0")
                    export_out(fileout,"Note: %s data generation in cycle % d is all done! % 5d DFT structures are generated in % s" % (daes_.upper(),cyc,N2,datdir+daes_+"/*/"+CYCLES),color = 3)

          # parsing the new set of the structures
          submit_prs(cwd,mini_,datdir,prsdir,fileout,cyc,m_setup)
    
          # training NN model with the parsed data
          submit_trn(cwd,mini_,libdir,prsdir,fileout,"jtrn",cyc,m_setup)

          # TEST data generation (if any)
          if m_setup.ATST == 1: # collect TEST data if ATST run!
               reread_stp(fileout,m_setup,cwd+"/setup")
               DATGEN_TST(RUN,output,cyc,1,m_setup)
               if does__exst(mtst_+"/"+CYCLES+"/"+dtst_):
                    if SORT == 1 or SORT == 2 or SORT == 0 or SORT == 3:
                         (N,direc,files) = searchfile(mtst_+"/"+CYCLES+"/"+dtst_,"POSCAR.0")
                         for i in range(0,N):
                              s = direc[i].replace(mtst_+"/"+CYCLES+"/"+dtst_+"/","")
                              j = s[4:10]  # these lines depend on the SORT scheme!
                              c = s[11:13]
                              d = s[14:17]
                              duplicates(mtst_+"/"+CYCLES+"/"+dtst_+"/"+dtst_+j+"/"+c+"/"+d,datdir+dtst_+j+"/"+c+"/"+d)
                              duplicates(mini_+"/"+iprs_+"/tag-tst",datdir+dtst_+j+"/tag")

               if SORT == 1 or SORT == 2 or SORT == 3:
                    (N3,direc,files) = searchfile(datdir+dtst_+"*/"+CYCLES,"POSCAR.0")
                    export_out(fileout,"Note: %s data generation in cycle % d is all done! % 5d DFT structures are generated in % s" % (dtst_.upper(),cyc,N3,datdir+dtst_+"*/"+CYCLES),color = 3)
               if SORT == 0:
                    (N3,direc,files) = searchfile(datdir+dtst_+"/*/"+CYCLES,"POSCAR.0")
                    export_out(fileout,"Note: %s data generation in cycle % d is all done! % 5d DFT structures are generated in % s" % (dtst_.upper(),cyc,N3,datdir+dtst_+"/*/"+CYCLES),color = 3)

          export_out(fileout,"Note: data generations in cycle % d is completed with total % 5d DFT structures (% 5d %s, % 5d %s, % 5d %s)" % (cyc,N1+N2+N3,N1,runt,N2,daes_.upper(),N3,dtst_.upper()),color = 3)

#----------------------------------------------------------------------------------

def DATGEN_TST(RUN, output, cyc, cal, setup_0): # TEST THE NN MODEL

     if setup_0.NSPC > 3:
          return
     # type of the data which is being generated
     runt = dtst_.upper()

     SORT = setup_0.SORT

     # local copy of setup
     t_setup = deepcopy(setup_0)

     cwd = os.getcwd()
     fileout = cwd+"/"+output

     export_out(fileout,"\n                      ===== running in %s mode =====\n" % runt ,time = 0,color = 3)

     CYCLES  = format(cyc,"02d")

     basdir  = cwd+"/"+mtst_+"/"+CYCLES+"/"
     testout = basdir+"results.dat"
     datdir  = basdir+"dat/"              # path to DFT   runs
     srcdir  = basdir+"evo/"
     mindir  = basdir+"min/"              # path to minima of NN-ES
     model   = cwd+"/model"

     if does__exst(basdir):
          export_out(fileout,"Error: % s already exists" % basdir,color = 2);exit()
     if doesntexst(model):
          export_out(fileout,"Error: % s does not exist" % model,color = 2);exit()

     # initialize the variables according to the number of species
     t_setup.SPCN  = mkarray_2D(t_setup.NSPC,len(t_setup.ASPC))
     for j in range(0,t_setup.NSPC):
          for i in range(0,len(t_setup.ASPC)):
               t_setup.SPCN[j][i] = t_setup.SPCT[j][i]
     t_setup.SIZN = []
     for j in range(0,len(t_setup.ASPC)):
          t_setup.SIZN.append(t_setup.SIZT[j])
     t_setup.NPOP = []
     for j in range(0,len(t_setup.SIZN)):
          t_setup.NPOP.append(t_setup.TNPP)
     t_setup.NSW0  = 0
     t_setup.PGPA  = t_setup.TGPA[:]
     t_setup.PWGT  = t_setup.TPWT[:]
     t_setup.ITER  = t_setup.TITR
     t_setup.PNUM  = len(t_setup.TGPA)

     # ===== starting the task

     # prepare and submit ES runs for TEST run
     export_out(fileout,"Note: generating %s dataset in cycle % 3d ..." % (runt,cyc),color = 3)
     totmin  = 0
     mkdir_tree(mindir)
     for z in range(0,t_setup.PNUM):
          PRESSU = format(z,"02d")
          for i in range(0,len(t_setup.SIZN)):
               t_setup.NGEN  = t_setup.TNGN[i]
               reread_stp(fileout,t_setup,cwd+"/setup")
               POPULN = format(i,"02d")
               SRCDIR = srcdir+PRESSU+"/"+POPULN  # always indexed as: EVO/cycle/pressure/population
               if SORT == 1:
                    DATDIR = datdir+devo_+t_setup.THSH[z][i]+"/"+CYCLES   # indexed using hash
               if SORT == 3:
                    DATDIR = datdir+devo_+PRESSU+POPULN+"/"+CYCLES # indexed as: DAT/"evo"pressure"population"/cycle
               if SORT == 2:
                    DATDIR = datdir+devo_+PRESSU+"/"+CYCLES # indexed as: DAT/"evo"pressure/cycle/population
               if SORT == 0:
                    DATDIR = datdir+devo_+"/"+PRESSU+"/"+CYCLES # indexed as: DAT/"evo"pressure/cycle/population
               pop    = submit_esj(cwd,mini_,DATDIR,mindir,SRCDIR,fileout,iesn_,i,z,cyc,t_setup)
               totmin += pop

     export_out(fileout,"Note: collected % 4d minima in %s run" % (totmin,runt),color = 3)
          
     # submit DFT runs for "t_setup.TMIN" lowest NN-minima
     (NM,mdi,mfi) = searchfile(mindir,"P*")
     export_out(fileout,"Note: submitting % 5d DFT runs for selected NN minima of %s run at % s..." % (NM,runt,datdir),color = 3)
     for z in range(0,t_setup.PNUM):
          PRESSU = format(z,"02d")
          for i in range(0,len(t_setup.SIZN)):
               reread_stp(fileout,t_setup,cwd+"/setup")
               POPULN = format(i,"02d")
               if SORT == 1:
                    DATDIR = datdir+devo_+t_setup.THSH[z][i]+"/"+CYCLES   # indexed using hash
               if SORT == 3:
                    DATDIR = datdir+devo_+PRESSU+POPULN+"/"+CYCLES # indexed as: DAT/"evo"pressure"population"/cycle
               if SORT == 2:
                    DATDIR = datdir+devo_+PRESSU+"/"+CYCLES # indexed as: DAT/"evo"pressure/cycle/population
               if SORT == 0:
                    DATDIR = datdir+devo_+"/"+PRESSU+"/"+CYCLES # indexed as: DAT/"evo"pressure/cycle/population
               tag    = PRESSU+POPULN
               for j in range(0,t_setup.TMIN):
                    m = format(j,"03d")
                    if does__exst(mindir+"P"+CYCLES+tag+m):
                         run = DATDIR+"/"+m+"000"
                         mkdir_tree(run)
                         # try to symmetrize structures before the DFT runs!
                         duplicates(mindir+"P"+CYCLES+tag+m,"POSCAR")
                         os.system("./"+mini_+"/maise -spg 0.1 1000 |tail -n 1 > tmpspg")
                         s = file_nthln("tmpspg",0).split()
                         if does__exst("CONV") and len(s) == 5 and is___float(s[4]) and float(s[4]) >= 0.99:
                              duplicates("CONV",run+"/POSCAR") # spg solver worked; maise >= 2.7
                         else:
                              duplicates(mindir+"P"+CYCLES+tag+m,run+"/POSCAR")  # spg solver failed
                         removes_fd("POSCAR CONV PRIM str.cif tmpspg")
                         duplicates(cwd+"/"+mini_+"/INCAR", run)
                         replacestr(run+"/INCAR","PPPP",str(t_setup.PGPA[int(z)]*10.0))
                         replacestr(run+"/INCAR","AAAA",str(t_setup.PREC))
                         replacestr(run+"/INCAR","GGGG",str(t_setup.ISIF))
                         replacestr(run+"/INCAR","SMSM",str(t_setup.SMER))
                         replacestr(run+"/INCAR","SGSG",str(t_setup.SIGM))
                         if t_setup.SYMP != 0.0:
                              addstrfile(run+"/INCAR","SYMPREC="+str(t_setup.SYMP))
                         if t_setup.SORB > 0:
                              addstrfile(run+"/INCAR","LSORBIT=.TRUE.")
                              addstrfile(run+"/INCAR","LMAXMIX="+str(t_setup.SORB))
                         if t_setup.ECUT > 0.0:
                              l = "ENCUT="+str(t_setup.ECUT)
                              addstrfile(run+"/INCAR",l)
     submit_dft(cwd,mini_,datdir,fileout,"jeos",t_setup)

     # summarize and write the output
     export_out(fileout,"        N    P     ANN_symmetry     ANN_vol     ANN_H         DFT_H        DFT_vol  DFT_symmetry    ANN_path                                 DFT_path",time = 0)
     export_out(fileout,"-----------------------------------------------------------------------------------------------------------------------------------------------------",time = 0)
     export_out(testout,"        N    P     ANN_symmetry     ANN_vol     ANN_H         DFT_H        DFT_vol  DFT_symmetry    ANN_path                                 DFT_path",time = 0)
     export_out(testout,"-----------------------------------------------------------------------------------------------------------------------------------------------------",time = 0)

     for z in range(0,t_setup.PNUM):
          for i in range(0,len(t_setup.SIZN)):
               spg1N  = []
               spg2N  = []
               spg1D  = []
               spg2D  = []
               clrsN  = []
               clrsD  = []
               enerN  = []
               voleN  = []
               enerD  = []
               voleD  = []
               ranks  = []
               adrsN  = []
               adrsD  = []
               PRESSU = format(z,"02d")
               POPULN = format(i,"02d")
               if SORT == 1:
                    DATDIR = datdir+devo_+t_setup.THSH[z][i]+"/"+CYCLES   # indexed using hash
               if SORT == 3:
                    DATDIR = datdir+devo_+PRESSU+POPULN+"/"+CYCLES # indexed as: DAT/"evo"pressure"population"/cycle
               if SORT == 2:
                    DATDIR = datdir+devo_+PRESSU+"/"+CYCLES # indexed as: DAT/"evo"pressure/cycle/population
               if SORT == 0:
                    DATDIR = datdir+devo_+"/"+PRESSU+"/"+CYCLES # indexed as: DAT/"evo"pressure/cycle/population
               tag    = PRESSU+POPULN
               for j in range(0,t_setup.TMIN):
                    m = format(j,"03d")
                    if does__exst(mindir+"P"+CYCLES+tag+m):
                         duplicates(mindir+"P"+CYCLES+tag+m,"POSCAR")
                         os.system("./"+mini_+"/maise -spg 0.1 1000 |tail -n 1 > tmpspg")
                         if file_nthln("tmpspg",0) == "":
                              spg1N.append(" ")
                              spg2N.append(" ")
                         else:
                              spg1N.append(file_nthln("tmpspg",0).split()[0])
                              spg2N.append(file_nthln("tmpspg",0).split()[1])
                         removes_fd("POSCAR CONV PRIM str.cif tmpspg")
                         src = file_nthln(mindir+"P"+CYCLES+tag+m,0).split()[0]
                         (k,l,s) = findlststr(src+"/OUTCAR.0","y=") # to read enthalpy/atom
                         enerN.append(float(s.split()[7])) # for MAISE outcar: enthalpy/atom
                         pos = import_str(mindir+"P"+CYCLES+tag+m)
                         voleN.append(float(volume_str(pos)/float(pos.NATM)))
                         clrsN.append(j)
                         adrsN.append(src.replace(mtst_+"/"+CYCLES+"/",""))
                    run = DATDIR+"/"+m+"000"
                    if does__exst(run+"/OUTCAR.0"):
                         duplicates(run+"/CONTCAR.0","POSCAR")
                         os.system("./"+mini_+"/maise -spg 0.1 1000 |tail -n 1 > tmpspg")
                         if file_nthln("tmpspg",0) == "":
                              spg1D.append(" ")
                              spg2D.append(" ")
                         else:
                              spg1D.append(file_nthln("tmpspg",0).split()[0])
                              spg2D.append(file_nthln("tmpspg",0).split()[1])
                         removes_fd("POSCAR CONV PRIM str.cif tmpspg")
                         pos = import_str(run+"/CONTCAR.0")
                         voleD.append(float(volume_str(pos)/float(pos.NATM)))
                         if t_setup.PGPA[z] == 0.0:
                              (k,l,s) = findlststr(run+"/OUTCAR.0","y=") # to read energy/structure
                              enerD.append(float(s.split()[6])/float(pos.NATM))
                         else:
                              (k,l,s) = findlststr(run+"/OUTCAR.0","enth") # to read enthalpy/structure
                              enerD.append(float(s.split()[4])/float(pos.NATM))
                    else:
                         spg1D.append("--")
                         spg2D.append("--")
                         voleD.append(9999.0)
                         enerD.append(9999.0)
                    clrsD.append(j)
                    adrsD.append(run.replace(basdir,""))

               ranks = sorted(range(len(enerD)), key=lambda k: enerD[k])

               composition  = ""
               for k in range(0,t_setup.NSPC-1):
                    composition += str(t_setup.SPCN[k][i])+":"
               composition += str(t_setup.SPCN[t_setup.NSPC-1][i])

               for j in range(0,len(enerN)):
                    r  = ranks[j]
                    cN = j+1 # color for ANN structure
                    cD = r+1 # color for DFT structure
                    export_out(fileout," % 8s  % 5.1lf  % 5s  % 5s  % 10.3lf %s % 10.4lf %s  %s % 10.4lf %s  % 10.3lf  % 5s  % 5s    % 25s  % 25s" % (composition,t_setup.PGPA[z],spg1N[j],spg2N[j],voleN[j],COLOR_CODES[cN],enerN[j],COLOR_CODES[0],COLOR_CODES[cD],enerD[r],COLOR_CODES[0],voleD[r],spg1D[r],spg2D[r],adrsN[j],adrsD[r]),time = 0)
                    export_out(testout," % 8s  % 5.1lf  % 5s  % 5s  % 10.3lf %s % 10.4lf %s  %s % 10.4lf %s  % 10.3lf  % 5s  % 5s    % 25s  % 25s" % (composition,t_setup.PGPA[z],spg1N[j],spg2N[j],voleN[j],COLOR_CODES[cN],enerN[j],COLOR_CODES[0],COLOR_CODES[cD],enerD[r],COLOR_CODES[0],voleD[r],spg1D[r],spg2D[r],adrsN[j],adrsD[r]),time = 0)

               export_out(fileout,"-------------------------------------------------------------------------------------------",time = 0)
               export_out(testout,"-------------------------------------------------------------------------------------------",time = 0)

               # collect the DFT energy of NN minima in test/ directory
               for j in range(0,len(enerN)):
                    r  = ranks[j]
                    if does__exst(basdir+adrsD[r]+"/OUTCAR.0") and does__exst(basdir+adrsD[r]+"/POSCAR.0"):
                         atstdir = cwd+"/"+mtst_+"/"+CYCLES+"/"+dtst_+"/"+dtst_+t_setup.THSH[int(PRESSU)][int(POPULN)]+"/"+CYCLES+"/"+format(j,"03d")
                         (o,tmpn,tmps) = readoutcar(basdir+adrsD[r]+"/OUTCAR.0", read = "first")
                         if tmpn > 0:
                              mkdir_tree(atstdir)
                              datdat_str(fileout,tmps[0],atstdir)
                              duplicates(basdir+adrsD[r]+"/POSCAR.0",atstdir+"/POSCAR.0")
                              duplicates(cwd+"/"+mini_+"/"+iprs_+"/tag-tst",cwd+"/"+mtst_+"/"+CYCLES+"/"+dtst_+"/"+dtst_+t_setup.THSH[int(PRESSU)][int(POPULN)]+"/tag")
                              if doesntexst(cwd+"/"+mtst_+"/"+CYCLES+"/"+dtst_+"/"+dtst_+t_setup.THSH[int(PRESSU)][int(POPULN)]+"/inf"):
                                   s = ""
                                   for j in range(0,t_setup.NSPC):
                                        s += str(t_setup.SPCT[j][int(POPULN)])+" "
                                   addstrfile(cwd+"/"+mtst_+"/"+CYCLES+"/"+dtst_+"/"+dtst_+t_setup.THSH[int(PRESSU)][int(POPULN)]+"/inf","pressure % s" % (str(t_setup.TGPA[int(PRESSU)])))
                                   addstrfile(cwd+"/"+mtst_+"/"+CYCLES+"/"+dtst_+"/"+dtst_+t_setup.THSH[int(PRESSU)][int(POPULN)]+"/inf","composit % s" % (s))

     export_out(fileout,"Note: %s run for cycle % d is completed" % (runt,cyc),color = 3)

#----------------------------------------------------------------------------------
