#================================================================#
#                                                                #
#                          maise-net                             #
#                                                                #
#                  Python interface to maise                     #
#                for automated data generation                   #
#                                                                #
#                                                                #
#                 version 1.0.05   07/08/2020                    # 
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

#================================================================#
#                                                                #
#                      HOW-TO and NOTES                          #
#                                                                #
#  1) This wrapper is developed for automated dataset generation #
#     with MAISE and VASP codes for neural network potentials.   #
#     It generates: single atom reference data (REFS), 2-4 atoms #    
#     clusters (CLST), equation of state for high-symmetry       #
#     structures (EOSZ), and evolutionary datasets of bulk or    # 
#     nanoparticle structures with a set of arbitrary sizes as   #
#     defined in "setup" (EVOS). Optionally, EOSN (EOS data for  #
#     the unique minima in EVOS data) and TEST (DFT relaxation   #
#     for NN-based minima found in the test evolutionary search) #
#     edata can be generated, as well.                           #      
#                                                                #
# 2)  maise-net works with any out-of-the-box version of the     #
#     Python (2.X/3.X) without any non-standard module.          #
#                                                                #
# 3)  Local optimizations and energy calculations are performed  #
#     with VASP package. VASP should be already installed on     #
#     user's system with "vasp" executable in the system path.   #
#                                                                #
# 4)  "maise" executive and proper "POTCAR" for the system       #
#     should be provided in the local INI/ directory.            #
#                                                                #
# 5)  Following files are job related and cluster specific;      #
#     their header should be manually adjusted in the indicated  #
#     part of the scripts, while the type of the cluster is      #
#     defined in the "setup" file with QUET flag:                #
#                                                                #
#     INI/jesz     for EVO run in cycle 0                        #   
#     INI/jesn     for EVO run in cycles 1+                      #    
#     INI/jdft     for DFT runs                                  #
#     INI/jtrn     for training jobs                             #
#     INI/jeos     for DFT runs in EOS and TST runs              #
#                                                                #
#================================================================#

import os.path
import sys
sys.path.insert(0, "INI/src")
from maisenet_defs  import *
from maisenet_main  import *
from maisenet_inot  import *
from maisenet_task  import initit_stp

#================================================================#
#                        LOAD THE SETUP                          #
#================================================================#

# Output file name
output = mout_

export_hdr(output,"version 1.0.05   07/08/2020")

if not os.path.exists("setup"):
     export_out(output,"Error: setup file does not exist",color = 1);exit()
     
setup = import_stp("setup")

if setup.JOBT//10 != 8:
     export_out(output,"Error: JOBT is not data generation",color = 1);exit()

RUN = setup.JOBT%10

if RUN == JQTE:
     export_out(output,"Note: terminating at user's request!",color = 2);exit()

#================================================================#
#                 INITIALIZE THE RUN PARAMETERS                  #
#================================================================#

initit_stp(RUN,mini_,output,setup)

#================================================================#
#                      REF DATA GENERATION                       #
#================================================================#

if (RUN == JBAS):

     DATGEN_REF(RUN,output,setup)

#================================================================#
#                      EOS DATA GENERATION                       #
#================================================================#

if (RUN == JBAS):

     DATGEN_EOS(RUN,output,0,setup)

#================================================================#
#                      CLS DATA GENERATION                       #
#================================================================#

if (RUN == JBAS):

     DATGEN_CLS(RUN,output,setup)

#================================================================#
#                      EVO DATA GENERATION                       #
#================================================================#

if (RUN == JEVO):

     if setup.SITR > setup.NITR:
          exit()

     for cyc in range(setup.SITR,setup.NITR+1):
          DATGEN_EVO(RUN,output,cyc,setup)

#================================================================#
#                       TST THE NN MODEL                         #
#================================================================#

if (RUN == JTST):

     DATGEN_TST(RUN,output,setup.NITR,0,setup)

#================================================================#
#                            THE END!                            #
#================================================================#
