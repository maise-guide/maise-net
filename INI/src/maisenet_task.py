#================================================================#
#                                                                #
#                       maise-net  TOOLBOX                       #
#                                                                #
#                    Python interface to maise                   #
#                  for automated data generation                 #
#                                                                #
#       Routines for various data generation related tasks       #
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
import time
from maisenet_defs  import *
from maisenet_util  import *
from maisenet_inot  import *

#----------------------------------------------------------------------------------          

def submit_esj(cwd, inidir, datdir, mindir, run, fileout, estype, i, z, cyc, setup_0):

     duplicates(cwd+"/"+inidir+"/"+estype,run)
     duplicates(cwd+"/"+inidir+"/maise",run)

     tag = format(z,"02d")+format(i,"02d")  # pressur+population: for MIN/P* name tags and DAT/evo*

     # for cycle z run
     if estype == iesz_:
          popsize = int(round(setup_0.npop[i]*setup_0.PWGT[z]))
          totrun  = len(setup_0.SIZ0)*setup_0.PNUM
          indrun  = len(setup_0.SIZ0)*z+i+1
          duplicates(cwd+"/"+inidir+"/jesz",run+"/INI/g0")
          duplicates(cwd+"/"+inidir+"/INCAR",run+"/INI/INCAR0")
          duplicates(cwd+"/"+inidir+"/maisenet_kmsh.py",run+"/INI")
          duplicates(cwd+"/"+inidir+"/POTCAR",run+"/INI")
          replacestr(run+"/INI/INCAR0","PPPP",str(setup_0.PGPA[z]*10.0))
          replacestr(run+"/INI/INCAR0","GGGG",str(setup_0.ISIF))
          replacestr(run+"/INI/INCAR0","SMSM",str(setup_0.SMER))
          replacestr(run+"/INI/INCAR0","SGSG",str(setup_0.SIGM))
          replacestr(run+"/INI/INCAR0","AAAA",str(setup_0.prec))
          replacestr(run+"/INI/g0","WWWW",str(setup_0.nsw0))
          replacestr(run+"/INI/g0","KKKK",str(setup_0.kdns))
          if not setup_0.SYMP == 0.0:
               addstrfile(run+"/INI/INCAR0","SYMPREC="+str(setup_0.SYMP))
          if setup_0.SORB > 0:
               addstrfile(run+"/INI/INCAR0","LSORBIT=.TRUE.")
               addstrfile(run+"/INI/INCAR0","LMAXMIX="+str(setup_0.SORB))
          if setup_0.ecut > 0.0:
               s = "ENCUT="+str(setup_0.ecut)
               addstrfile(run+"/INI/INCAR0",s)
          if setup_0.KDNS != setup_0.kdns:
               s = "python maisenet_kmsh.py "+str(setup_0.KDNS)
               replacestr(run+"/INI/g0","#kmesh",s)
          if setup_0.PREC != setup_0.prec:
               s = "sed -i '/^PREC/d' INCAR"
               replacestr(run+"/INI/g0","#prec0",s)
               s = "echo 'PREC="+str(setup_0.PREC)+"' >> INCAR"
               replacestr(run+"/INI/g0","#prec1",s)
          if setup_0.ECUT != setup_0.ecut:
               s = "sed -i '/^ENCUT/d' INCAR"
               replacestr(run+"/INI/g0","#encut0",s)
               if setup_0.ECUT > 0.0:
                    s = "echo 'ENCUT="+str(setup_0.ECUT)+"' >> INCAR"
                    replacestr(run+"/INI/g0","#encut1",s)
          duplicates(run+"/setupz",run+"/setup")
          s = " "
          for j in range(0,setup_0.NSPC):
               s += str(setup_0.SPC0[j][i])+" "
          replacestr(run+"/setup","BBBB",s)

     # for cycle n run
     if estype == iesn_:
          popsize = int(round(setup_0.NPOP[i]*setup_0.PWGT[z]))
          totrun  = len(setup_0.SIZN)*setup_0.PNUM
          indrun  = len(setup_0.SIZN)*z+i+1
          duplicates(cwd+"/"+inidir+"/jesn",run+"/INI/g0")
          duplicates(cwd+"/"+inidir+"/maise",run+"/INI/NNET/maise-nnet")
          duplicates(cwd+"/model",run+"/INI/NNET")
          replacestr(run+"/INI/NNET/setup","PPPP",str(setup_0.PGPA[z]))
          replacestr(run+"/INI/NNET/setup","QQQQ",str(setup_0.NPAR))
          replacestr(run+"/INI/NNET/setup","GGGG",str(setup_0.ISIF))
          replacestr(run+"/INI/NNET/setup","DDDD",str(setup_0.NDIM))
          replacestr(run+"/INI/NNET/setup","OOOO",str(setup_0.ITER))
          duplicates(run+"/setupn",run+"/setup")
          s = " "
          for j in range(0,setup_0.NSPC):
               s += str(setup_0.SPCN[j][i])+" "
          replacestr(run+"/setup","BBBB",s)

     replacestr(run+"/setup","VV01",str(setup_0.TETR))
     replacestr(run+"/setup","VV02",str(setup_0.PLNT))
     replacestr(run+"/setup","VV03",str(setup_0.PACK))
     replacestr(run+"/setup","VV04",str(setup_0.BLOB))
     replacestr(run+"/setup","VV05",str(setup_0.MATE))
     replacestr(run+"/setup","VV06",str(setup_0.SWAP))
     replacestr(run+"/setup","VV07",str(setup_0.RUBE))
     replacestr(run+"/setup","VV08",str(setup_0.REFL))
     replacestr(run+"/setup","VV09",str(setup_0.INVS))
     replacestr(run+"/setup","VV10",str(setup_0.CHOP))
     replacestr(run+"/setup","VV11",str(setup_0.MUTE))
     replacestr(run+"/setup","VV12",str(setup_0.MCRS))
     replacestr(run+"/setup","VV13",str(setup_0.SCRS))
     replacestr(run+"/setup","VV14",str(setup_0.LCRS))
     replacestr(run+"/setup","VV15",str(setup_0.ACRS))
     replacestr(run+"/setup","VV16",str(setup_0.SDST))
     replacestr(run+"/setup","VV17",str(setup_0.LDST))
     replacestr(run+"/setup","VV18",str(setup_0.ADST))
     replacestr(run+"/setup","VV19",str(setup_0.ELPS))
     replacestr(run+"/setup","DDDD",str(setup_0.NDIM))
     replacestr(run+"/setup","HHHH",str(setup_0.QUET))
     replacestr(run+"/setup","ZZZZ",str(setup_0.NSPC))
     replacestr(run+"/setup","XXXX",str(setup_0.LBOX))
     replacestr(run+"/setup","TTTT",str(setup_0.ZATM))
     replacestr(run+"/setup","PPPP",str(setup_0.PGPA[z]))
     replacestr(run+"/setup","MMMM",str(popsize))
     if estype == iesz_:
          replacestr(run+"/setup","GGGG",str(setup_0.ngen))
     if estype == iesn_:
          replacestr(run+"/setup","GGGG",str(setup_0.NGEN))

     os.chdir(run)

     print("Update: submitted ES job ( % 3d / % 3d ) at % s" % (indrun,totrun,run))
     os.system("./maise 2>/dev/null")

     # collect minima for cycle z; no analysis by MAISE; all minima will be collected!
     if estype == iesz_:
          mkdir_tree("POOL")
          pop = 0
          (n,l,s) = findallstr("log.out","The energy diversity is too small")
          if does__exst("ebest.dat") or n > 0:
               for j in range(0,setup_0.ngen):
                    for k in range(popsize,popsize*2):
                         m = format(pop, "03d")
                         r = "EVOS/G"+format(j, "03d")+"/M"+format(k, "03d")
                         if checks_dft(r+"/OUTCAR.1") < 0:
                              export_out(fileout,"Warning: DFT run was not successful in % s" % (r),time = 0)
                              renames_fd(r+"/POSCAR.1",r+"/POSCAR.error")
                              renames_fd(r+"/CONTCAR.1",r+"/CONTCAR.error")
                              renames_fd(r+"/OUTCAR.1",r+"/OUTCAR.error")
                         else:
                              des     = datdir+"/"+m+"000"
                              mkdir_tree(des)
                              (o,p,s) = readoutcar(r+"/OUTCAR.1", read = "last")
                              nam     = run.replace(cwd+"/","")+"/"+r+" "+str(s[0].ENTH)+" "+des.replace(cwd+"/","")
                              duplicates(r+"/CONTCAR.1","POOL/POSCAR"+m)
                              replacenln("POOL/POSCAR"+m,0,nam)
                              duplicates("POOL/POSCAR"+m,mindir+"/P"+format(cyc,"02d")+tag+m)
                              pop    += 1
               export_out(fileout,"          produced % 5d ES-minima from % s" % (pop,run),time = 0)
          else:
               export_out(fileout,"Warning: ES was not successful in % s" % (run),time = 0)

     # collect minima for cycle n; POOL is generated through MAISE's analysis
     if estype == iesn_:
          pop = 0
          (n,l,s) = findallstr("log.out","The energy diversity is too small")
          if does__exst("ebest.dat") or n > 0:
               replacestr("setup","JOBT  10","JOBT  13")
               os.system("./maise")
               if doesntexst("POOL"):
                    export_out(fileout,"Warning: ES was not successful in % s" % (run),time = 0)
               else:
                    (n,di,fi) = searchfile("POOL/","POSCAR*")
                    for j in range(0,n):
                         m    = format(j,"03d")
                         des  = datdir+"/"+m+"000"
                         s    = file_nthln("POOL/POSCAR"+m,0)
                         nam  = run.replace(cwd+"/","")+"/"+s.split()[0]+" "+s.split()[1]+" "+des.replace(cwd+"/","")
                         replacenln("POOL/POSCAR"+m,0,nam)
                         duplicates("POOL/POSCAR"+m,mindir+"/P"+format(cyc,"02d")+tag+m)
                         pop += 1
                    export_out(fileout,"          produced % 5d ES-minima from % s" % (pop,run),time = 0)
          else:
               export_out(fileout,"Warning: ES was not successful in % s" % (run),time = 0)

     os.chdir(cwd)

     return pop

#----------------------------------------------------------------------------------          
def checks_dft(fname): # Checks an OUTCAR file to see if run is properly finished and converged

     if doesntexst(fname):
          return -1 # file does not exist
     (n,p,l) = findlststr(fname,"al CPU")
     if n == 0:
          return -2 # not finished
     (n,p,l) = findlststr(fname,"Itera")
     if n > 0:
          ele = int(find_betwn(l,"(",")"))
          if ele == 60:
               return -3 # not converges
          else:
               return 1 # run is OK!
     else:
          return -2 # not finished/started

#----------------------------------------------------------------------------------          
def submit_dft(cwd, inidir, datdir, fileout, jobfile, setup_0): # Create a list of collected directories and submit DFT runs

    isub = 0
    dsub = [] # keeps the list of directories that are submitted here
    MAXJ = reread_stp(fileout,cwd+"/setup")

    # number of new jobs to be submitted
    (N,direc,files) = searchfile(datdir,"POSCAR")
    nsub = [0  for i in range(N)] # will tag directories that are submitted here with 1 value
    jsub = ["" for i in range(N)] # will hold the jobid for those submitted here
    for j in range(0,N):
        # check how many of those submitted are still running
        irun = 0
        for i in range(0,isub):
            if doesntexst(dsub[i]+"/OUTCAR.0"):
                irun += 1
        while irun >= MAXJ:
            irun = 0
            for i in range(0,isub):
                if doesntexst(dsub[i]+"/OUTCAR.0"):
                    irun += 1
            time.sleep(5)
            MAXJ = reread_stp(fileout,cwd+"/setup")

        if (doesntexst(direc[j]+"/OUTCAR.0")) and (doesntexst(direc[j]+"/OUTCAR.error")):
            os.chdir(direc[j])
            scrfile = "j"+format(isub+1,"05d")
            duplicates(cwd+"/"+inidir+"/POTCAR"  , ".")
            duplicates(cwd+"/"+inidir+"/maisenet_kmsh.py" , ".")
            duplicates(cwd+"/"+inidir+"/"+jobfile, scrfile)
            if does__exst("KPOINTS"):
                 eliminates(scrfile,"KKKK")
            else:
                 replacestr(scrfile,"KKKK",str(setup_0.KDNS))
            if doesntexst("INCAR"):
                duplicates(cwd+"/"+inidir+"/INCAR"   , ".")
                if not setup_0.NSW0 == 0:
                    replacestr("INCAR","NSW=0","NSW=%s" % str(setup_0.NSW0))
                replacestr("INCAR","PPPP",str(setup_0.PGPA[0]*10.0))
                replacestr("INCAR","AAAA",str(setup_0.PREC))
                replacestr("INCAR","GGGG",str(setup_0.ISIF))
                replacestr("INCAR","SMSM",str(setup_0.SMER))
                replacestr("INCAR","SGSG",str(setup_0.SIGM))
                if not setup_0.SYMP == 0.0:
                     addstrfile("INCAR","SYMPREC="+str(setup_0.SYMP))
                if setup_0.SORB > 0:
                     addstrfile("INCAR","LSORBIT=.TRUE.")
                     addstrfile("INCAR","LMAXMIX="+str(setup_0.SORB))
                if setup_0.ECUT > 0.0:
                     s = "ENCUT="+str(setup_0.ECUT)
                     addstrfile("INCAR",s)
            os.system(setup_0.CLUS[0]+" "+scrfile+" > tmpa;date >> tmpa")
            s = file_nthln("tmpa",0)
            t = file_nthln("tmpa",1)
            if setup_0.QUET == 0:
                jobid = s.split(".")[0]
            elif setup_0.QUET == 1:
                jobid = s.split()[len(s.split())-1]
            elif setup_0.QUET == 2:
                jobid = ""
            addstrfile("id",jobid)
            addstrfile("id",t)
            removes_fd("tmpa")
            print("Update: submitted DFT job % s id % s ( % 5d / % 5d ) at % s" % (scrfile,jobid,isub+1,N,direc[j]))
            nsub[j] = 1
            jsub[j] = jobid
            dsub.append(direc[j])
            isub += 1
            os.chdir(cwd)
        MAXJ = reread_stp(fileout,cwd+"/setup")

    if isub == 0:
        export_out(fileout,"Note: no DFT job is submitted!",color = 2)
    else:
        # make sure all DFT runs are finished!
        jobs = 0
        fail = []
        print("\n-------------------- waiting for DFT jobs to finish --------------------")
        while jobs < N:
            MAXJ = reread_stp(fileout,cwd+"/setup")
            time.sleep(30)
            jobs = 0
            for j in range(0,N):
                 if does__exst(direc[j]+"/OUTCAR.0"):
                      jobs += 1
                 else:
                      if nsub[j] == 1 and (not j in fail):
                           if not jobrunning(jsub[j],setup_0.QUET,setup_0.CLUS):
                                fail.append(j)
                                print("Warning: the DFT job is not (finished or running) in % s" % direc[j])
        print("Update: DFT jobs are all finished!\n")

    # check for unfinished DFT runs; clean up the unfinished directories; write dat.dat files
    (N,direc,files) = searchfile(datdir,"OUTCAR.0")
    tot = 0
    for j in range(0,N):
         if doesntexst(direc[j]+"/dat.dat"):
              if checks_dft(direc[j]+"/OUTCAR.0") < 0:  
                   export_out(fileout,"Warning: DFT run was not successful in % s" % (direc[j]),time = 0)
                   renames_fd(direc[j]+"/POSCAR.0",direc[j]+"/POSCAR.error")
                   renames_fd(direc[j]+"/CONTCAR.0",direc[j]+"/CONTCAR.error")
                   renames_fd(direc[j]+"/OUTCAR.0",direc[j]+"/OUTCAR.error")
              else:
                   datdat_out(fileout,direc[j]+"/OUTCAR.0",direc[j])
                   tot += 1

    time.sleep(5)
    return tot

#----------------------------------------------------------------------------------                
def jobrunning(jobid, QUET, CLUS): # check if a job is still running on cluster; returns logical value

     c = 0
     # for slurm queue
     if QUET == 1:
          os.system("%s -j %s > tmpz 2>&1" % (CLUS[2],jobid))
     # for torque queue
     if QUET == 0:
          os.system("%s % s   > tmpz 2>&1" % (CLUS[2],jobid))

     f= open("tmpz","r")
     O=f.read().split("\n")
     f.close()

     for i in range(0,len(O)-1):
          if jobid == O[i].split()[0]:
               c = 1

     removes_fd("tmpz")

     if c == 1:
          return True
     else:
          return False

#----------------------------------------------------------------------------------                

def submit_prs(cwd, inidir, datdir, prsdir, fileout, cyc, setup_0):                     # Parsing the collected data

    MAXJ = reread_stp(fileout,cwd+"/setup")

    export_out(fileout,"Note: parsing structures for cycle = % 3d in % s ..." % (cyc,prsdir),color = 2)

    duplicates(cwd+"/"+inidir+"/"+iprs_,prsdir)
    duplicates(cwd+"/"+inidir+"/maise",prsdir)
    os.chdir(prsdir)
    duplicates(str(setup_0.RCUT)+str(setup_0.NSYM)+".dat","basis")
    duplicates("setupp","setup")
    replacestr("setup","TTTT",str(setup_0.ZATM))
    replacestr("setup","ZZZZ",str(setup_0.NSPC))
    replacestr("setup","YYYY",str(setup_0.NSYM))
    replacestr("setup","CCCC",str(setup_0.NCMP))
    replacestr("setup","IIII",datdir)
    replacestr("setup","DDDD",str(setup_0.NDIM))
    replacestr("setup","RRRR",str(setup_0.FMRK))
    replacestr("setup","QQQQ",str(setup_0.NPAR))
    os.system("./maise 2>/dev/null")
    if does__exst("stamp.dat"):
        os.chdir(cwd)
    else:
        export_out(fileout,"Error: parsing did not finish in % s" % (prsdir),color = 1);exit()

#----------------------------------------------------------------------------------                

def submit_trn(cwd, inidir, libdir, prsdir, fileout, jobfile, cyc, setup_0):      # Training NN model using the parsed data

     MAXJ = reread_stp(fileout,cwd+"/setup")

     duplicates(cwd+"/"+inidir+"/"+itrn_,libdir)
     os.chdir(libdir)

     # set the training type if subsystem models are present in INI/
     if setup_0.NSPC == 2:
          model0 = cwd+"/"+inidir+"/"+atomtosymb(setup_0.TSPC[0])+".dat"
          model1 = cwd+"/"+inidir+"/"+atomtosymb(setup_0.TSPC[1])+".dat"
          if does__exst(model0) and does__exst(model1):
               duplicates(model0,libdir)
               duplicates(model1,libdir)
               setup_0.STRT = 1
          else:
               setup_0.STRT = 0
     if setup_0.NSPC == 3:
          model0 = cwd+"/"+inidir+"/"+atomtosymb(setup_0.TSPC[0])+atomtosymb(setup_0.TSPC[1])+".dat"
          model1 = cwd+"/"+inidir+"/"+atomtosymb(setup_0.TSPC[0])+atomtosymb(setup_0.TSPC[2])+".dat"
          model2 = cwd+"/"+inidir+"/"+atomtosymb(setup_0.TSPC[1])+atomtosymb(setup_0.TSPC[2])+".dat"
          if does__exst(model0) and does__exst(model1) and does__exst(model2):
               duplicates(model0,libdir)
               duplicates(model1,libdir)
               duplicates(model2,libdir)
               setup_0.STRT = 1
          else:
               setup_0.STRT = 0
               
     if setup_0.STRT == 0:
          export_out(fileout,"Note: full training of neural network model for cycle = % 3d in % s ..." % (cyc,libdir),color = 2)
     else:
          export_out(fileout,"Note: stratified training of neural network model for cycle = % 3d in % s ..." % (cyc,libdir),color = 2)

     if cyc == 0:
          t_mitr = setup_0.mitr
          t_tefs = setup_0.tefs
          t_extr = 1
     else:
          t_mitr = setup_0.MITR
          t_tefs = setup_0.TEFS
          t_extr = setup_0.EXTR

     ii = -1 # to set the starting point for training
     # if t_extr is negative (last cycle): no initial E-only training!
     if cyc > 0 and cyc == setup_0.NITR and t_extr < 0:
          for i in range(0,len(t_tefs)):
               if t_tefs[i] == 0:
                    ii = i
     t_extr = abs(t_extr)

     for i in range(ii+1,len(t_mitr)):
          duplicates(cwd+"/"+inidir+"/"+jobfile, libdir)
          duplicates(cwd+"/"+inidir+"/maise",    libdir)
          MAXJ = reread_stp(fileout,cwd+"/setup")
          duplicates("setupt","setup")
          if setup_0.STRT == 1:
               replacestr("setup","JOBT  40","JOBT  41")
          replacestr("setup","EEEE",str(t_tefs[i]))
          replacestr("setup","TTTT",str(setup_0.ZATM))
          replacestr("setup","YYYY",str(setup_0.NSYM))
          replacestr("setup","ZZZZ",str(setup_0.NSPC))
          replacestr("setup","CCCC",str(setup_0.NCMP))
          replacestr("setup","IIII",prsdir)
          replacestr("setup","DDDD",str(setup_0.NDIM))
          replacestr("setup","LLLL",str(setup_0.LREG))
          replacestr("setup","AAAA",str(len(setup_0.NNNU)))
          s1 = " "
          s2 = " "
          for j in range(0,len(setup_0.NNNU)):
               s1 += str(setup_0.NNNU[j])+" "
               s2 += str(setup_0.NNGT[j])+" "
          replacestr("setup","UUUU",s1)
          replacestr("setup","GGGG",s2)
          if cyc > 0 and cyc == setup_0.NITR and t_tefs[i] == 1 and t_extr > 0:
               replacestr("setup","SSSS",str(t_mitr[i]*t_extr))
          else:
               replacestr("setup","SSSS",str(t_mitr[i]))
          if cyc == 0:
               replacestr("setup","NTRN -90","NTRN -100")
               replacestr("setup","NTST -10","NTST  0  ")
          # submitting the trainin job
          os.system(setup_0.CLUS[0]+" "+jobfile+" > tmpa;date >> tmpa")
          s = file_nthln("tmpa",0)
          t = file_nthln("tmpa",1)
          if setup_0.QUET == 0:
               jobid = s.split(".")[0]
          elif setup_0.QUET == 1:
               jobid = s.split()[len(s.split())-1]
          elif setup_0.QUET == 2:
               jobid = ""
          addstrfile("id",jobid)
          addstrfile("id",t)
          removes_fd("tmpa")
          print("Update: submitted NN training job % s id % s ( % 4d / % 4d ) for cycle % d" % (jobfile,jobid,i+1,len(t_mitr),cyc))
          c = 0
          while c == 0:
               time.sleep(5)
               if does__exst("stamp"):
                    c = 1
          removes_fd("stamp")

          # Analysis of NN testing erros and plotting them (if gnuplot exists!)
          mdlerr_chk(fileout,setup_0)

          nndir = "nn"+format(i, "02d")
          mkdir_tree(nndir)
          duplicates("model",nndir)
          duplicates("err-ene.dat",nndir)
          duplicates("err-out.dat",nndir)
          duplicates("err-max.dat",nndir)
          for j in range(0,setup_0.PNUM):
               duplicates("err-ent%d.dat" % j,nndir)
          duplicates("err-frc.dat",nndir)
          duplicates("%s-ent.png" % (setup_0.NAME),nndir)
          duplicates("%s-his.png" % (setup_0.NAME),nndir)
          duplicates("model",cwd)
          if i < (len(t_mitr)-1):
               removes_fd("err-*.dat *.png")

     if setup_0.STRT == 0:
          replacestr("setup","JOBT  40","JOBT  50")
     else:
          replacestr("setup","JOBT  41","JOBT  50")
     os.system("./maise")

     os.chdir(cwd)
     MAXJ = reread_stp(fileout,cwd+"/setup")

#----------------------------------------------------------------------------------                

def similr_chk(cwd, inidir, mindir, datdir, fileout, cyc, jobfile, NM, MINIMA, setup_0): # Check for similarities based on DFT-relaxed minima

     export_out(fileout, "Note: checking data for similarities in cycle = % d" % cyc,color = 2)

     SORT = setup_0.SORT

     nrej = 0

     kdns = setup_0.KDNS
     nsw0 = setup_0.NSW0
     setup_0.NSW0 = 31
     setup_0.KDNS = setup_0.kdns

     for z in range(0,setup_0.PNUM):
          MAXJ = reread_stp(fileout,cwd+"/setup")
          for i in range(0,len(setup_0.SIZN)):
               C   = format(cyc,"02d")
               P   = format(z,"02d")
               R   = format(i,"02d")
               (fls,pat,nam) = searchfile(mindir,"P"+C+P+R+"*")     
               for j in range(0,fls):
                    I   = nam[j][7:10]
                    if SORT == 1:
                         adr = datdir+devo_+P+R+"/"+C+"/"+I
                    if SORT == 2:
                         adr = datdir+devo_+P+"/"+C+"/"+R+"/"+I
                    if SORT == 0:
                         adr = datdir+devo_+"/"+P+"/"+C+"/"+R+"/"+I
                    if does__exst(adr+"000/OUTCAR.0"):
                         run = adr+"00t"
                         mkdir_tree(run)
                         duplicates(adr+"000/CONTCAR.0",run+"/POSCAR")
                         duplicates(cwd+"/"+inidir+"/INCAR", run)
                         replacestr(run+"/INCAR","NSW=0","NSW=%s" % str(31))
                         replacestr(run+"/INCAR","PPPP",str(setup_0.PGPA[int(P)]*10.0))
                         replacestr(run+"/INCAR","AAAA",str(setup_0.prec))
                         replacestr(run+"/INCAR","GGGG",str(setup_0.ISIF))
                         replacestr(run+"/INCAR","SMSM",str(setup_0.SMER))
                         replacestr(run+"/INCAR","SGSG",str(setup_0.SIGM))
                         if not setup_0.SYMP == 0.0:
                              addstrfile(run+"/INCAR","SYMPREC="+str(setup_0.SYMP))
                         if setup_0.SORB > 0:
                              addstrfile(run+"/INCAR","LSORBIT=.TRUE.")
                              addstrfile(run+"/INCAR","LMAXMIX="+str(setup_0.SORB))
                         if setup_0.ecut > 0.0:
                              s = "ENCUT="+str(setup_0.ecut)
                              addstrfile(run+"/INCAR",s)
     submit_dft(cwd,inidir,datdir,fileout,jobfile,setup_0)

     setup_0.NSW0 = nsw0
     setup_0.KDNS = kdns
     for z in range(0,setup_0.PNUM):
          MAXJ = reread_stp(fileout,cwd+"/setup")
          for i in range(0,len(setup_0.SIZN)):
               C   = format(cyc,"02d")
               P   = format(z,"02d")
               R   = format(i,"02d")
               (fls,pat,nam) = searchfile(mindir,"P"+C+P+R+"*")     
               for j in range(0,fls):
                    I   = nam[j][7:10]
                    if SORT == 1:
                         adr = datdir+devo_+P+R+"/"+C+"/"+I
                    if SORT == 2:
                         adr = datdir+devo_+P+"/"+C+"/"+R+"/"+I
                    if SORT == 0:
                         adr = datdir+devo_+"/"+P+"/"+C+"/"+R+"/"+I
                    if does__exst(adr+"00t/OUTCAR.0"):
                         run = adr+"00r"
                         mkdir_tree(run)
                         duplicates(adr+"00t/CONTCAR.0",run+"/POSCAR")
                         duplicates(cwd+"/"+inidir+"/INCAR", run)
                         replacestr(run+"/INCAR","PPPP",str(setup_0.PGPA[int(P)]*10.0))
                         replacestr(run+"/INCAR","AAAA",str(setup_0.PREC))
                         replacestr(run+"/INCAR","GGGG",str(setup_0.ISIF))
                         replacestr(run+"/INCAR","SMSM",str(setup_0.SMER))
                         replacestr(run+"/INCAR","SGSG",str(setup_0.SIGM))
                         if not setup_0.SYMP == 0.0:
                              addstrfile(run+"/INCAR","SYMPREC="+str(setup_0.SYMP))
                         if setup_0.SORB > 0:
                              addstrfile(run+"/INCAR","LSORBIT=.TRUE.")
                              addstrfile(run+"/INCAR","LMAXMIX="+str(setup_0.SORB))
                         if setup_0.ECUT > 0.0:
                              s = "ENCUT="+str(setup_0.ECUT)
                              addstrfile(run+"/INCAR",s)
                         renames_fd(adr+"00t/CONTCAR.0",adr+"00t/CONTCAR.t")
                         renames_fd(adr+"00t/POSCAR.0", adr+"00t/POSCAR.t")
                         renames_fd(adr+"00t/OUTCAR.0", adr+"00t/OUTCAR.t")
     submit_dft(cwd,inidir,datdir,fileout,jobfile,setup_0)

     # assign the weights for minima
     for z in range(0,setup_0.PNUM):
          MAXJ = reread_stp(fileout,cwd+"/setup")
          for i in range(0,len(setup_0.SIZN)):
               C   = format(cyc,"02d")
               P   = format(z,"02d")
               R   = format(i,"02d")
               ENE = []
               emin = 100000.0
               emax = -100000.0
               (fls,pat,nam) = searchfile(mindir,"P"+C+P+R+"*")
               for j in range(0,fls):
                    I   = nam[j][7:10]
                    if SORT == 1:
                         run = datdir+devo_+P+R+"/"+C+"/"+I+"00r"
                    if SORT == 2:
                         run = datdir+devo_+P+"/"+C+"/"+R+"/"+I+"00r"
                    if SORT == 0:
                         run = datdir+devo_+"/"+P+"/"+C+"/"+R+"/"+I+"00r"
                    if doesntexst(run+"/OUTCAR.0"):
                         ENE.append(1e+8)
                    else:
                         pos = import_str(run+"/CONTCAR.0")
                         (n,l,s) = findlststr(run+"/OUTCAR.0","enth")
                         if n > 0:
                              e = float(s.split()[4])/float(pos.NATM)
                         else:
                              (n,l,s) = findlststr(run+"/OUTCAR.0","y=")
                              e = float(s.split()[6])/float(pos.NATM)
                         if e < emin:
                              emin = e
                         if e > emax:
                              emax = e
                         ENE.append(e)
               d = 1.1*abs(emax-emin) # 0.1 is to give chance to highest energy!
               if d < 1e-8:
                    d = 1.0
               t = -1
               for j in range(0,fls):
                    I   = nam[j][7:10]
                    fi  = "P"+C+P+R+I
                    t  += 1
                    if ENE[t] == 1e+8:
                         MINIMA[fi] = 0.0
                         nrej += 1
                    else:
                         MINIMA[fi] = 1.0 - abs(ENE[t]-emin)/d

     # check against previous cycles' DFT-relaxed; if any
     if cyc > 1:
          ene2 = []
          pth2 = []
          vol2 = []
          tot2 = 0
          for i in range(1,cyc):
               MAXJ = reread_stp(fileout,cwd+"/setup")
               if SORT == 1:
                    (N,di,fi) = searchfile(datdir+"/"+devo_+"*/"+format(i,"02d")+"/*r/","OUTCAR.0")
               if SORT == 2:
                    (N,di,fi) = searchfile(datdir+"/"+devo_+"*/"+format(i,"02d")+"/*/*r/","OUTCAR.0")
               if SORT == 0:
                    (N,di,fi) = searchfile(datdir+"/"+devo_+"/*/"+format(i,"02d")+"/*/*r/","OUTCAR.0")
               for j in range(0,N):
                    pos = import_str(di[j]+"/CONTCAR.0")
                    (k,l,s) = findlststr(di[j]+"/OUTCAR.0","y=") # to read energy/structure
                    ene2.append(float(s.split()[6])/float(pos.NATM))
                    vol2.append(float(volume_str(pos))/float(pos.NATM))
                    pth2.append(di[j])
                    tot2 += 1
          if tot2 > 0:
               duplicates(cwd+"/"+inidir+"/maise",cwd+"/")
               rdfs = 1
               for j in range(0,NM):
                    MAXJ = reread_stp(fileout,cwd+"/setup")
                    fi   = list(MINIMA.keys())[j]
                    C    = fi[1:3]
                    P    = fi[3:5]
                    R    = fi[5:7]
                    I    = fi[7:10]
                    if SORT == 1:
                         run = datdir+devo_+P+R+"/"+C+"/"+I+"00r"
                    if SORT == 2:
                         run = datdir+devo_+P+"/"+C+"/"+R+"/"+I+"00r"
                    if SORT == 0:
                         run = datdir+devo_+"/"+P+"/"+C+"/"+R+"/"+I+"00r"
                    c   = 0
                    if MINIMA[fi] > 0.0:
                         pos  = import_str(run+"/CONTCAR.0")
                         (k,l,s) = findlststr(run+"/OUTCAR.0","y=") # to read energy/structure                
                         ene1 = float(s.split()[6])/float(pos.NATM)
                         vol1 = float(volume_str(pos)/float(pos.NATM))
                         for k in range(0,tot2):
                              if c == 0:
                                   if abs(vol1 - vol2[k]) < defdlvol and abs(ene1 - ene2[k]) < defdlene:
                                        if rdfs == 1:
                                             export_out(fileout, "          List of RDF checks performed for (dftwegt) in cycle = % 4d" % cyc,time = 0)
                                        duplicates(run+"/CONTCAR.0", "POSCAR0")
                                        duplicates(pth2[k]+"/CONTCAR.0", "POSCAR1")
                                        os.system("./maise -cxc 1000 |tail -n 1 > tmprdf")
                                        val = file_nthln("tmprdf", 0)
                                        export_out(fileout, "               (% 4d) % 20s      with % 20s   rdf = % lf" % (rdfs,run+"/CONTCAR.0",pth2[k]+"/CONTCAR.0",float(val)),time = 0)
                                        rdfs += 1
                                        if float(val.split()[0]) > defdlrdf:
                                             c += 1
                         if c > 0:
                              MINIMA[fi] = 0.0 # mark structure with zero weight (rejected)
                              nrej += 1

          removes_fd("maise tmprdf CONTCAR POSCAR0 POSCAR1 RDF*dat")

     # normalize the weights                  
     totwgt = 0.0
     for i in range(0,NM):
          totwgt += MINIMA[list(MINIMA.keys())[i]]
     if totwgt > 0.0:
          for i in range(0,NM):
               MINIMA[list(MINIMA.keys())[i]] /= totwgt

     export_out(fileout, "Note: similarity check is finished: % 4d structures rejected from % 4d (cycle = % 4d)" % (nrej, NM, cyc),color = 2)
     MAXJ = reread_stp(fileout,cwd+"/setup")

#----------------------------------------------------------------------------------                

def initit_stp(RUN, inidir, output, setup_0): # Check initial directories; and initiate the required variables for runs

     cwd = os.getcwd()

     fileout = cwd+"/"+output

     # check the directories to begin with
     if doesntexst(cwd+"/"+inidir):
          export_out(fileout,"Error: %s directory does not exist" % inidir,color = 1);exit()
     if doesntexst(cwd+"/"+inidir+"/maise"):
          export_out(fileout,"Error: %s does not exist" % (inidir+"/maise"),color = 1);exit()
     if doesntexst(cwd+"/"+inidir+"/POTCAR"):
          export_out(fileout,"Error: %s does not exist" % (inidir+"/POTCAR"),color = 1);exit()

     # initializing variables and arrays with setup values
     if setup_0.nsw0 < 0:
          setup_0.nsw0 = setup_0.NSW0
     if setup_0.ecut < 0.0:
          setup_0.ecut = setup_0.ECUT
     if setup_0.prec == " ":
          setup_0.prec = setup_0.PREC
     if setup_0.kdns < 0:
          setup_0.kdns = setup_0.KDNS
     if setup_0.smer < 0:
          setup_0.smer = setup_0.SMER
     if setup_0.sigm < 0:
          setup_0.sigm = setup_0.SIGM

     # Queue type & CLUS [0] submission command; [1] deletion command; [2] status
     if setup_0.QUET > 2:
          export_out(fileout,"Error: Queue option % d is not recognized!" % setup_0.QUET,color = 1);exit()
     if   setup_0.QUET == 0:
          setup_0.CLUS[0] = "qsub"
          setup_0.CLUS[1] = "qdel"
          setup_0.CLUS[2] = "qstat"
     elif setup_0.QUET == 1:
          setup_0.CLUS[0] = "sbatch"
          setup_0.CLUS[1] = "scancel"
          setup_0.CLUS[2] = "squeue"
     elif setup_0.QUET == 2:
          setup_0.CLUS[0] = "bsub <"
          setup_0.CLUS[1] = ""
          setup_0.CLUS[2] = ""
          
     if (not type(setup_0.TSPC) is list):
          export_out(fileout,"Error: TSPC is not defined or it is not a [list]",color = 1);exit()

     # total number of species in the setup file
     setup_0.NSPC = len(setup_0.TSPC)

     # list of species with their Z value
     for i in range(0,setup_0.NSPC):
          setup_0.ZATM += str(setup_0.TSPC[i])+" "
          setup_0.SATM += atomtosymb(setup_0.TSPC[i])+" "
          setup_0.NAME += atomtosymb(setup_0.TSPC[i])

     # average atomic radius of elements in the system     
     for i in range(0,setup_0.NSPC):
          setup_0.RADI += float(atomic_rad(setup_0.TSPC[i]))/float(setup_0.NSPC)

     # check the POTCAR file for species
     (n,l,s) = findallstr(cwd+"/"+inidir+"/POTCAR","VRHFIN")
     if n != setup_0.NSPC:
         export_out(fileout,"Error: number of elements in TSPC list (% 3d) and %s/POTCAR (% 3d) do not match" %(setup_0.NSPC,inidir,n),color = 1);exit()
     for j in range(0,n):
         if not atomtosymb(setup_0.TSPC[j]) in s[j]:
             export_out(fileout,"Error: element number % 3d (z = % 3d) in the TSPC list does not exist in %s/POTCAR file" % (j+1,setup_0.TSPC[j],inidir),color = 1);exit()

     # adjust the ISIF for EVOS according to dimensionality
     if setup_0.NDIM == 3:
          setup_0.ISIF = 3
     if setup_0.NDIM == 0:
          setup_0.ISIF = 2
          if setup_0.LBOX == 0.0:
               export_out(fileout,"Error: LBOX can not be zero when NDIM=0",color = 1);exit()

     # initiate TEST run variables
     if setup_0.TMIN == 0:
          setup_0.TMIN = tstmins
     if setup_0.TITR == 0:
          setup_0.TITR = tstiter
     if setup_0.TNGN == 0:
          setup_0.TNGN = tstngen
     if setup_0.TNPP == 0:
          setup_0.TNPP = tstnpop
     if len(setup_0.TGPA) == 0 or (not type(setup_0.TGPA) is list):
          setup_0.TGPA = tstpgpa[:]
     if len(setup_0.TPWT) == 0 or (not type(setup_0.TPWT) is list):
          setup_0.TPWT = tstpwgt[:]
     if len(setup_0.TSTA) == 0 or (not type(setup_0.TSTA) is list):
          if setup_0.NSPC == 1:
               setup_0.TSTA = tst1ASPC[:]
          if setup_0.NSPC == 2:
               setup_0.TSTA = tst2ASPC[:]
          if setup_0.NSPC == 3:
               setup_0.TSTA = tst3ASPC[:]
     if len(setup_0.TSTB) == 0 or (not type(setup_0.TSTB) is list):
          if setup_0.NSPC == 2:
               setup_0.TSTB = tst2BSPC[:]
          if setup_0.NSPC == 3:
               setup_0.TSTB = tst3BSPC[:]
     if len(setup_0.TSTC) == 0 or (not type(setup_0.TSTC) is list):
          if setup_0.NSPC == 3:
               setup_0.TSTC = tst3CSPC[:]

     # number of pressure values for EVOS runs
     setup_0.PNUM = len(setup_0.PGPA)

     # training and parsing settings (for EVOS run modes)
     if RUN == JEVO: 
          G2 = 0; G4 = 0
          f = cwd+"/"+inidir+"/"+iprs_+"/"+str(setup_0.RCUT)+format(setup_0.NSYM,"02d")+".dat"
          for i in range(5,5+setup_0.NSYM):
               s = int(file_nthln(f,i).split()[1])
               if s == 2:
                    G2 += 1
               if s == 4:
                    G4 += 1
          c = 0
          for i in range(setup_0.NSPC,0,-1):
               c += i
          setup_0.NCMP = (setup_0.NSPC*G2+c*G4)
          
          if len(setup_0.NNNU) != len(setup_0.NNGT):
               export_out(fileout,"Error: length of neuron/layer list does not match that of activation functions list",color = 1);exit()

     # for EVOS (cycle = 0 ) run modes
     if RUN == JEVO and setup_0.SITR == 0: 
          if len(setup_0.mitr) != len(setup_0.tefs):
               export_out(fileout,"Error: mitr and tefs dimensions %d and %d are not correct" % (len(setup_0.mitr),len(setup_0.tefs)),color = 1);exit()

          if (not type(setup_0.npop) is list):
               export_out(fileout,"Error: npop is not defined or it is not a [list]",color = 1);exit()

          setup_0.SPC0 = mkarray_2D(setup_0.NSPC,len(setup_0.npop))
          if setup_0.NSPC >= 1:
               if (not type(setup_0.aspc) is list):
                    export_out(fileout,"Error: aspc is not defined or it is not a [list]",color = 1);exit()
               if len(setup_0.aspc) != len(setup_0.npop):
                    export_out(fileout,"Error: aspc and npop do not match",color = 1);exit()
               for i in range(0,len(setup_0.npop)):
                    setup_0.SPC0[0][i] = setup_0.aspc[i]
          if setup_0.NSPC >= 2:
               if (not type(setup_0.bspc) is list):
                    export_out(fileout,"Error: bspc is not defined or it is not a [list]",color = 1);exit()
               if len(setup_0.bspc) != len(setup_0.npop):
                    export_out(fileout,"Error: bspc and npop do not match",color = 1);exit()
               for i in range(0,len(setup_0.npop)):
                    setup_0.SPC0[1][i] = setup_0.bspc[i]
          if setup_0.NSPC == 3:
               if (not type(setup_0.cspc) is list):
                    export_out(fileout,"Error: cspc is not defined or it is not a [list]",color = 1);exit()
               if len(setup_0.cspc) != len(setup_0.npop):
                    export_out(fileout,"Error: cspc and npop do not match",color = 1);exit()
               for i in range(0,len(setup_0.npop)):
                    setup_0.SPC0[2][i] = setup_0.cspc[i]

          setup_0.SIZ0 = []
          for j in range(0,len(setup_0.npop)):
               c = 0
               for i in range(0,setup_0.NSPC):
                    c += setup_0.SPC0[i][j] 
               setup_0.SIZ0.append(c)

     # for EVOS (cycle = 1+ ) run modes
     if RUN == JEVO and setup_0.NITR > 0: 
          if len(setup_0.MITR) != len(setup_0.TEFS):
               export_out(fileout,"Error: MITR and TEFS dimensions %d and %d are not correct" % (len(setup_0.MITR),len(setup_0.TEFS)),color = 1);exit()

          if (not type(setup_0.NPOP) is list):
               export_out(fileout,"Error: NPOP is not defined or it is not a [list]",color = 1);exit()

          setup_0.SPCN = mkarray_2D(setup_0.NSPC,len(setup_0.NPOP))
          if setup_0.NSPC >= 1:
               if (not type(setup_0.ASPC) is list):
                    export_out(fileout,"Error: ASPC is not defined or it is not a [list]",color = 1);exit()
               if len(setup_0.ASPC) != len(setup_0.NPOP):
                    export_out(fileout,"Error: ASPC and NPOP do not match",color = 1);exit()
               for i in range(0,len(setup_0.NPOP)):
                    setup_0.SPCN[0][i] = setup_0.ASPC[i]
          if setup_0.NSPC >= 2:
               if (not type(setup_0.BSPC) is list):
                    export_out(fileout,"Error: BSPC is not defined or it is not a [list]",color = 1);exit()
               if len(setup_0.BSPC) != len(setup_0.NPOP):
                    export_out(fileout,"Error: BSPC and NPOP do not match",color = 1);exit()
               for i in range(0,len(setup_0.NPOP)):
                    setup_0.SPCN[1][i] = setup_0.BSPC[i]
          if setup_0.NSPC == 3:
               if (not type(setup_0.CSPC) is list):
                    export_out(fileout,"Error: CSPC is not defined or it is not a [list]",color = 1);exit()
               if len(setup_0.CSPC) != len(setup_0.NPOP):
                    export_out(fileout,"Error: CSPC and NPOP do not match",color = 1);exit()
               for i in range(0,len(setup_0.NPOP)):
                    setup_0.SPCN[2][i] = setup_0.CSPC[i]

          setup_0.SIZN = []
          for j in range(0,len(setup_0.NPOP)):
               c = 0
               for i in range(0,setup_0.NSPC):
                    c += setup_0.SPCN[i][j] 
               setup_0.SIZN.append(c)

#----------------------------------------------------------------------------------                

def mdlerr_chk(fileout, setup_0): # Analysis of testing errors and "relative enthalpy-errors"

     SORT = setup_0.SORT

     if doesntexst("err-ene.dat") or doesntexst("model"):
          export_out(fileout,"Error: err-ene.dat or model file do(es) not exist",color = 1);exit()
     if does__exst("setup"):
          (k,l,s) = findlststr("setup","DATA")
          prsdir  = s.split()[1]
     else:
          prsdir  = ""

     # find and collect structures with test error >= train_error
     Ntot = file_nmlns("err-ene.dat")
     if Ntot > 0:
          (k,l,s) = findlststr("model","train energy error")
          ERROR = float(s.split()[5])
          (k,l,s) = findlststr("model","test  energy data")
          Ntst = int(s.split()[5])
          if Ntst == 0:
               print("No test data for analysis!")
               return
          pminim = []
          pcount = []
          pindex = []
          for i in range(0,setup_0.PNUM):
               pminim.append(100000.0)
               pcount.append(0)
               pindex.append(i)
          efiles = []
          errors = []
          addres = []
          pgpaes = []
          dftene = []
          enthal = []
          volatm = []
          for i in range(Ntot-Ntst,Ntot):
               s = file_nthln("err-ene.dat",i)
               if len(s.split()) == 6:
                    route = s.split()[5]
               elif does__exst(prsdir+"/"+s.split()[1]):
                    route = file_nthln(prsdir+"/"+s.split()[1],2)
               if doesntexst(route):
                    route = "no_parsed_path_exists"

               if does__exst(route+"/POSCAR.0"):
                    pos = import_str(route+"/POSCAR.0")
                    volatm.append(float(volume_str(pos)/float(pos.NATM)))
                    if does__exst(route+"/dat.dat"):
                         n,p,o = findlststr(route+"/dat.dat","pres")
                         if n > 0:
                              for j in range(0,setup_0.PNUM):
                                   if float(o.split()[1]) == setup_0.PGPA[j]:
                                        P = str(j)
                    elif devo_ in route:
                         s = route.split(devo_,1)[1][0:3]
                         if s[0] == "/":
                              P = s[1:3]
                         else:
                              P = s[0:2]
                    elif setup_0.PNUM == 1:
                         P = "0"
                    else:
                         continue
                    try:
                         num = int(P)
                    except:
                         continue

                    if int(P) >= setup_0.PNUM:
                         continue
                         
                    addres.append(route)
                    efiles.append(s.split()[1])
                    errors.append(float(s.split()[4]))
                    dftene.append(float(s.split()[2]))
                    enthal.append(dftene[len(dftene)-1]+setup_0.PGPA[int(P)]*volatm[len(volatm)-1]/(160.21766208))
               else:
                    P = "0"
                    volatm.append(0.0)
                    addres.append(route)
                    efiles.append(s.split()[1])
                    errors.append(float(s.split()[4]))
                    dftene.append(float(s.split()[2]))
                    enthal.append(dftene[len(dftene)-1])

               pgpaes.append(int(P))
               pcount[int(P)] += 1
               
               if enthal[len(enthal)-1] < pminim[int(P)]:
                    pminim[int(P)] = enthal[len(enthal)-1]

          pindex = [[0.0 for i in range(max(pcount))] for j in range(setup_0.PNUM)]
          for i in range(0,setup_0.PNUM):
               k = 0
               for j in range(0,len(enthal)):
                    if pgpaes[j] == i:
                         pindex[i][k] = j
                         k += 1

          for i in range(0,setup_0.PNUM):
               if pcount[i] > 0:
                    f = open("tmp1234321","w")
                    for j in range(0,pcount[i]):
                         ind = pindex[i][j]
                         f.write(" % 10.2lf  % 10.2lf  % s\n" % (1000.0*errors[ind],1000.0*(enthal[ind]-pminim[i]),addres[ind]))
                    f.close()
                    os.system("sort -g -k2,2 tmp1234321 > err-ent%d.dat" % i)
                    removes_fd("tmp1234321")

          f      = open("err-max.dat","w")
          n      = 1
          abserr = [abs(k) for k in errors]
          srt    = sorted(range(len(abserr)),key  = lambda k: abserr[k], reverse = True)
          f.write("=====  num     efile       tst_err    ratio    structure\n")
          for j in range(0,setup_0.PNUM):
               f.write("===== % 6.2lf GPa\n" % setup_0.PGPA[j])
               for i in range(0,len(efiles)):
                    ind = srt[i]
                    if abs(errors[ind]) >= 1.0*ERROR and pgpaes[ind] == j:
                         f.write("      % 04d  % 10s  % 10.2lf ( % 5.3lf )  % s\n" % (n,efiles[ind],1000.0*errors[ind],abs(errors[ind])/ERROR,addres[ind]))
                         n += 1
          f.close()
               
          mkplot_err(fileout,setup_0)
