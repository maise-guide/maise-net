#================================================================#
#                                                                #
#                       maise-net  TOOLBOX                       #
#                                                                #
#                    Python interface to maise                   #
#                  for automated data generation                 #
#                                                                #
#           Constants, parameters, and data structures           #
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

# ======================================================================================        
# Global file storage conventions (actual directory names will be set by address book)
# ======================================================================================   
#
#    INI/src contains the source code (used in maisenet.py for starting the code)
#    INI/eos contains the EOS prototypes
#    INI/trn contains training files
#    INI/prs contains parsing files
#    INI/esz contain evolutionary files for cycle 0
#    INI/esn contain evolutionary files for cycle 1+
#
#    EVO directory is organized as: EVO/cycle/pressure/population/run*** with:
#             run000 -> DFT static run of minima itself
#             run00t -> DFT optimizaion of minia during similarity check
#             run00r -> DFT static run of the relaxed minia during similarity check
#             run### -> sampled trajectory for minima
#
#    MIN/P* files are indexed as: P_cycle_pressure_population_run
#
#    DAT directory contains the EVOS data as: evo_pressure-population/cycle/run* (SROT=1)
#    DAT directory contains the REFS data as: refs/00..##
#    DAT directory contains the CLST data as: clst/0_size/00..##
#    DAT directory contains the EOSZ data as: eosz_prototypename/00..##
#    DAT directory contains the EOSN data as: eosn_pressure-population/cycle/run/00..##
#    DAT directory contains the TEST data as: test_pressure-population/cycle/minimum
#
#    LIB directory contains models as: cycle/model and cycle/nn_training_round
#
#    PRS directory contains parsed data as: cycle/
#
#    TST directory contains MIN, EVO, and DAT as well

# ======================================================================================        
# Address book (address definitions are global)
# ======================================================================================   

# main output file
mout_ = "output.dat"
# main (root) directories
mini_ = "INI"
mevo_ = "EVO"
mdat_ = "DAT"
meos_ = "EOS"
mtst_ = "TST"
mlib_ = "LIB"
mprs_ = "PRS"
mmin_ = "MIN"
# directories to store DFT data: mdat_/"name"
devo_ = "evos"
deos_ = "eosz"
daes_ = "eosn"
dcls_ = "clst"
dref_ = "refs"
dtst_ = "test"
# main sources in "mini_" ~ INI/
ieos_ = "eos"
iesn_ = "esn"
iesz_ = "esz"
iprs_ = "prs"
itrn_ = "trn"
isrc_ = "src"

# ======================================================================================        
# Job type definitions (used in the code globally)
# ======================================================================================   

# job-type definitions
JBAS  = 0                         # generate BASIC data
JEVO  = 1                         # generate EVOS data
JTST  = 7                         # TEST run
JPSE  = 8                         # pause
JQTE  = 9                         # quit

# ======================================================================================        
# Constants and parameters
# ======================================================================================   

# constants
MAX_ATOMS = 200                   # max number of atoms of the code (used for POSCAR structures)
DIM       = 3                     # max dimensionality of the code

# hard-coded parameters for DFT and similarity check
defsigma = 0.01                   # sigma value for DFT run for REFS and CLST dat
defdlvol = 0.1                    # similarity measure for vol/atom in A^3/atom
defdlene = 0.001                  # similarity measure for ene/atom in eV/atom 
defdlrdf = 0.98                   # similarity measure for maise RDF analysis

# REF data parameters
rnumb  = 10                       # number of REFS copies

# EOS data parameters
eminz  = 0.8                      # minimum lattice parameter
emaxz  = 1.2                      # maximum lattice parameter
enumb  = 40                       # number of structures per prototype

# CLS data parameters
cminz  = 1.0                      # fraction of equil. distance
cmaxz  = 2.0                      # fraction of minz
cnumb  = 50                       # number of cluster per size

# TST run parameters
tstmins  = 10                              # max number of minima for DFT run
tstiter  = 500                             # NN optimization steps in NN-ES
tstnpop  = 16                              # population size for NN-ES
tstpgpa  = [0.0,10.0,30.0]                 # pressure values for NN-ES
tstpwgt  = [1.0, 1.0, 1.0]                 # weights for pressure values
tst1ASPC = [2,4,8]                         # atom/cell for elemental search
tstngen1 = [25,25,25]                      # number of generations for NN-ES
tst2ASPC = [1,1,2,3,1,3,2,3,5,4,2,6]             # atom/cell for binary search (element 1)
tst2BSPC = [1,2,1,1,3,2,3,5,3,4,6,2]             # atom/cell for binary search (element 2)
tstngen2 = [25,25,25,25,25,25,25,25,25,25,25,25] # number of generations for NN-ES
tst3ASPC = [1,1,2,1,2,1,2]                 # atom/cell for ternary search (element 1)
tst3BSPC = [1,1,1,2,2,2,1]                 # atom/cell for ternary search (element 2)
tst3CSPC = [1,2,1,1,1,2,2]                 # atom/cell for ternary search (element 3)
tstngen3 = [25,25,25,25,25,25,25]          # number of generations for NN-ES

# ======================================================================================           
# Physical constants and values
# ======================================================================================           

# Atomic covalent radii (from doi://10.1039/b801115j)
Latom_radi = {1:0.31,2:0.28,3:1.28,4:0.96,5:0.84,6:0.76,7:0.71,8:0.66,9:0.57,10:0.58,11:1.66,12:1.41,13:1.21,14:1.11,15:1.07,16:1.05,17:1.02,18:1.06,19:2.03,20:1.76,21:1.70,22:1.60,23:1.53,24:1.39,25:1.39,26:1.32,27:1.26,28:1.24,29:1.32,30:1.22,31:1.22,32:1.20,33:1.19,34:1.20,35:1.20,36:1.16,37:2.20,38:1.95,39:1.90,40:1.75,41:1.64,42:1.54,43:1.47,44:1.46,45:1.42,46:1.39,47:1.45,48:1.44,49:1.42,50:1.39,51:1.39,52:1.38,53:1.39,54:1.40,55:2.44,56:2.15,57:2.07,58:2.04,59:2.03,60:2.01,61:1.99,62:1.98,63:1.98,64:1.96,65:1.94,66:1.92,67:1.92,68:1.89,69:1.90,70:1.87,71:1.87,72:1.75,73:1.70,74:1.62,75:1.51,76:1.44,77:1.41,78:1.36,79:1.36,80:1.32,81:1.45,82:1.46,83:1.48,84:1.40,85:1.50,86:1.50,87:2.60,88:2.21,89:2.15,90:2.06,91:2.00,92:1.96,93:1.90,94:1.87,95:1.80}

# Atomic mass values (in atomic mass units; from IUPAC 2009-2017; those without tenths are approximations)
Latom_mass = {1:1.008,2:4.0026,3:6.94,4:9.0122,5:10.81,6:12.011,7:14.007,8:15.999,9:18.998,10:20.180,11:22.990,12:24.305,13:26.982,14:28.085,15:30.974,16:32.06,17:35.45,18:39.95,19:39.098,20:40.078,21:44.956,22:47.867,23:50.942,24:51.996,25:54.938,26:55.845,27:58.933,28:58.693,29:63.546,30:65.38,31:69.723,32:72.63,33:74.922,34:78.971,35:79.904,36:83.798,37:85.468,38:87.62,39:88.906,40:91.224,41:92.906,42:95.95,43:98,44:101.07,45:102.91,46:106.42,47:107.87,48:112.41,49:114.82,50:118.71,51:121.76,52:127.6,53:126.90,54:131.29,55:132.91,56:137.33,57:138.91,58:140.12,59:140.91,60:144.24,61:145,62:150.36,63:151.96,64:157.25,65:158.93,66:162.5,67:164.93,68:167.26,69:168.93,70:173.05,71:174.97,72:178.49,73:180.95,74:183.84,75:186.21,76:190.23,77:192.22,78:195.08,79:196.97,80:200.59,81:204.38,82:207.2,83:208.98,84:209,85:210,86:222,87:223,88:226,89:227,90:232.04,91:231.04,92:238.03,93:237,94:244,95:243}

# table of atomic numbr (Z) to element symbols
Latom_symb = {1:"H",2:"He",3:"Li",4:"Be",5:"B",6:"C",7:"N",8:"O",9:"F",10:"Ne",11:"Na",12:"Mg",13:"Al",14:"Si",15:"P",16:"S",17:"Cl",18:"Ar",19:"K",20:"Ca",21:"Sc",22:"Ti",23:"V",24:"Cr",25:"Mn",26:"Fe",27:"Co",28:"Ni",29:"Cu",30:"Zn",31:"Ga",32:"Ge",33:"As",34:"Se",35:"Br",36:"Kr",37:"Rb",38:"Sr",39:"Y",40:"Zr",41:"Nb",42:"Mo",43:"Tc",44:"Ru",45:"Rh",46:"Pd",47:"Ag",48:"Cd",49:"In",50:"Sn",51:"Sb",52:"Te",53:"I",54:"Xe",55:"Cs",56:"Ba",57:"La",58:"Ce",59:"Pr",60:"Nd",61:"Pm",62:"Sm",63:"Eu",64:"Gd",65:"Tb",66:"Dy",67:"Ho",68:"Er",69:"Tm",70:"Yb",71:"Lu",72:"Hf",73:"Ta",74:"W",75:"Re",76:"Os",77:"Ir",78:"Pt",79:"Au",80:"Hg",81:"Tl",82:"Pb",83:"Bi",84:"Po",85:"At",86:"Rn",87:"Fr",88:"Ra",89:"Ac",90:"Th",91:"Pa",92:"U",93:"Np",94:"Pu",95:"Am"}

# table of element symbols to atomic number (Z)
Lsymb_atom = {"H":1,"h":1,"He":2,"he":2,"Li":3,"li":3,"Be":4,"be":4,"B":5,"b":5,"C":6,"c":6,"N":7,"n":7,"O":8,"o":8,"F":9,"f":9,"Ne":10,"ne":10,"Na":11,"na":11,"Mg":12,"mg":12,"Al":13,"al":13,"Si":14,"si":14,"P":15,"p":15,"S":16,"s":16,"Cl":17,"cl":17,"Ar":18,"ar":18,"K":19,"k":19,"Ca":20,"ca":20,"Sc":21,"sc":21,"Ti":22,"ti":22,"V":23,"v":23,"Cr":24,"cr":24,"Mn":25,"mn":25,"Fe":26,"fe":26,"Co":27,"co":27,"Ni":28,"ni":28,"Cu":29,"cu":29,"Zn":30,"zn":30,"Ga":31,"ga":31,"Ge":32,"ge":32,"As":33,"as":33,"Se":34,"se":34,"Br":35,"br":35,"Kr":36,"kr":36,"Rb":37,"rb":37,"Sr":38,"sr":38,"Y":39,"y":39,"Zr":40,"zr":40,"Nb":41,"nb":41,"Mo":42,"mo":42,"Tc":43,"tc":43,"Ru":44,"ru":44,"Rh":45,"rh":45,"Pd":46,"pd":46,"Ag":47,"ag":47,"Cd":48,"cd":48,"In":49,"in":49,"Sn":50,"sn":50,"Sb":51,"sb":51,"Te":52,"te":52,"I":53,"i":53,"Xe":54,"xe":54,"Cs":55,"cs":55,"Ba":56,"ba":56,"La":57,"la":57,"Ce":58,"ce":58,"Pr":59,"pr":59,"Nd":60,"nd":60,"Pm":61,"pm":61,"Sm":62,"sm":62,"Eu":63,"eu":63,"Gd":64,"gd":64,"Tb":65,"tb":65,"Dy":66,"dy":66,"Ho":67,"ho":67,"Er":68,"er":68,"Tm":69,"tm":69,"Yb":70,"yb":70,"Lu":71,"lu":71,"Hf":72,"hf":72,"Ta":73,"ta":73,"W":74,"w":74,"Re":75,"re":75,"Os":76,"os":76,"Ir":77,"ir":77,"Pt":78,"pt":78,"Au":79,"au":79,"Hg":80,"hg":80,"Tl":81,"tl":81,"Pb":82,"pb":82,"Bi":83,"bi":83,"Po":84,"po":84,"At":85,"at":85,"Rn":86,"rn":86,"Fr":87,"fr":87,"Ra":88,"ra":88,"Ac":89,"ac":89,"Th":90,"th":90,"Pa":91,"pa":91,"U":92,"u":92,"Np":93,"np":93,"Pu":94,"pu":94,"Am":95,"am":95}

# ======================================================================================           
# Color codes (will be used in the code globally)
# ======================================================================================           

# ordered as: red; green; yellow; blue; magenta; cyan
COLOR_CODES = {0:"\033[0m",1:"\033[31m",2:"\033[32m",3:"\033[33m",4:"\033[34m",5:"\033[35m",6:"\033[36m",7:"\033[37m",8:"\033[31;1m",9:"\033[32;1m",10:"\033[33;1m",11:"\033[34;1m",12:"\033[35;1m",13:"\033[36;1m",14:"\033[37;1m"}

# ======================================================================================           
# POSCAR data class
# ======================================================================================           

class poscar:
    def __init__(self):
        self.NAME       = "maise"                                                 # name of structure                               
        self.SCAL       = 1.0                                                     # scaling factor                                  
        self.ENER       = 0.0                                                     # structure energy                                
        self.ENTH       = 0.0                                                     # structure enthalpy
        self.NATM       = 0                                                       # number of atoms                                 
        self.NTYP       = 0                                                       # number of species                               
        self.NAET       = []                                                      # number of atoms of each species                 
        self.SPCS       = []                                                      # symbol of species                               
        self.LATT       = [[0 for col in range(DIM)] for row in range(DIM)]       # lattice parameters                              
        self.POSI       = [[0 for col in range(DIM)] for row in range(MAX_ATOMS)] # atomic positions                                
        self.FORC       = [[0 for col in range(DIM)] for row in range(MAX_ATOMS)] # atomic forces                                   
        self.COOR       = 0                                                       # cartesian(0) or direct(1) coor.                 
        self.SLTV       = 0                                                       # non-selective(0) or selective dynam.(1)         
        self.DYNM       = [[0 for col in range(DIM)] for row in range(MAX_ATOMS)] # selective dynamics tags                         
        self.RECI       = [[0 for col in range(DIM)] for row in range(DIM)]       # reciprocal lattice vectors   
        self.VOLM       = 0.0                                                     # structure volume (A^3) from OUTCAR
        self.PRES       = 0.0                                                     # structure pressure from OUTCAR
        self.inkB       = ""                                                      # "in kB" values from OUTCAR            

# ======================================================================================           
# SETUP data class
# ======================================================================================           

class setup:
    def __init__(self):
        #--- those which are being read from SETUP
        #-----------------------------------------
        self.JOBT       = 0                     # type of the job
        self.TSPC       = []                    # list of atomic number (Z) of the elements in the system
        self.QUET       = 1                     # (0) torque (1) slurm (2) IBM-lsf
        self.SITR       = 0                     # start from this cycle (0 for complete run)
        self.NITR       = 0                     # total number of cycles (1,....)
        self.DATA       = 500                   # desired number of structures per cycle (1+)
        self.NNNU       = [10,10]               # neuron/hidden_layer for NN model (by default 2 hidden layers are used)
        self.NNGT       = [1,1]                 # activation function for hidden layer neurons (1 = hypertangent)
        self.LREG       = 1e-8                  # regularization
        self.STRT       = 0                     # (0) full or (1) stratified training of NN model
        self.SMLR       = 1                     # structure similarity check (0) none (1) unique minima + their relaxation path
        self.AEOS       = 1                     # (1) Active-EOS data generation from unique minima; (0) non Active-EOS
        self.ATST       = 1                     # (1) Active-TST data generation from NN MODEL TEST; (0) non Active-TST
        self.SORT       = 1                     # (1) pressure-based and (0) cycle-based data collection
        self.NPAR       = 8                     # number of cores used for parsing and maise runs
        self.RPAR       = 1                     # number of cores for NN relaxation jobs
        self.NSYM       = 51                    # number of symmetry function components
        self.RCUT       = 1                     # cut-off radius
        self.FMRK       = 0.5                   # fraction of atoms selected for force training
        self.NDIM       = 3                     # dimensionality of the unit cell
        self.LBOX       = 0.0                   # unit cell size: non-zero only for NDIM=0
        self.LSYM       = 1                     # (1) symmetric box for REFS data (0) non-symmetric box
        # evolutionary operations
        self.TETR       = 0.00                  # random using tetris
        self.PLNT       = 0.00                  # seeded
        self.PACK       = 0.00                  # biased
        self.BLOB       = 0.00                  # random using blob shape
        self.MATE       = 0.70                  # crossover using two halves
        self.SWAP       = 0.00                  # crossover using core-shell
        self.RUBE       = 0.00                  # Rubik's cube operation
        self.REFL       = 0.00                  # symmetrization via reflection
        self.INVS       = 0.00                  # symmetrization via inversion
        self.CHOP       = 0.00                  # chop to make facets
        self.MUTE       = 0.30                  # distortion
        self.MCRS       = 0.50                  # mutation rate in crossover
        self.SCRS       = 0.10                  # crossover:  swapping rate
        self.LCRS       = 0.10                  # crossover:  mutation strength for lattice vectors
        self.ACRS       = 0.10                  # crossover:  mutation strength for atomic positions
        self.SDST       = 0.10                  # distortion: swapping rate
        self.LDST       = 0.15                  # distortion: mutation strength for lattice vectors
        self.ADST       = 0.15                  # distortion: mutation strength for atomic positions
        self.ELPS       = 0.15                  # random:     nanoparticle ellipticity
        # search parameters
        self.aspc       = []                    # number of atoms for the first element (0 cycle)
        self.bspc       = []                    # number of atoms for the second element (0 cycle)
        self.cspc       = []                    # number of atoms for the third element (0 cycle)
        self.npop       = []                    # number of structures for each population (=sum aspc+bspc+cspc)
        # DFT settings for pre-static-run
        self.nsw0       = -1                    # ionic steps in DFT run
        self.ecut       = -1.0                  # energy cut-off for DFT (0 = POTCAR default)
        self.prec       = " "                   # precision of the DFT run (norm,acc) in pre-static-run
        self.kdns       = -1                    # k-mesh density for DFT runs in pre-static-run
        self.smer       = -1                    # vasp smearing in pre-static-run
        self.sigm       = -0.1                  # sigma for smearing in pre-static-run
         # model training
        self.mitr       = []                    # number of interations/training steps (0 cycle)
        self.tefs       = []                    # type of training at each step: E/E+F (0 cycle)
        # search parameters
        self.ASPC       = []                    # number of atoms for the first element (1+ cycle)
        self.BSPC       = []                    # number of atoms for the second element (1+ cycle)
        self.CSPC       = []                    # number of atoms for the third element (1+ cycle)
        self.NPOP       = []                    # number of structures for each population (1+ cycle)
        # NN-search relaxation
        self.ITER       = 25                    # neural network model relaxation steps
        # DFT run parameters for final static runs
        self.MAXJ       = 10                    # max number of DFT jobs to be submitted at once
        self.ECUT       = 0.0                   # energy cut-off for DFT (0=default)
        self.PREC       = "acc"                 # precision of the DFT run (norm,acc)
        self.KDNS       = 4000                  # k-mesh density for DFT runs
        self.SMER       = 1                     # vasp smearing                                           
        self.SIGM       = 0.1                   # sigma for smearing         
        self.PGPA       = [0.0,10.0,30.0,50.0]  # pressure values for data generation
        self.PWGT       = [1.0,0.5,0.25,0.25]   # weight for each pressure value (X population size)
        self.SORB       = 0                     # (0) no spin-orbit; non-zero = LMAXMIX value!
        self.SYMP       = 1.0E-8                # the SYMPREC flag in VASP; with zero value this will be truned off!
        # model training
        self.MITR       = []                    # number of interations/training steps (1+ cycle)
        self.TEFS       = []                    # type of training at each step (1+ cycle)
        self.EXTR       = 1                     # extended training: times number of steps for the last cycle
        self.ERNI       = 0                     # resume flag; for Erneso!
        # variables for TST run; have defaults defined above!
        self.TSTA       = []                    # list of number of atoms species A
        self.TSTB       = []                    # list of number of atoms species B
        self.TSTC       = []                    # list of number of atoms species C
        self.TGPA       = []                    # list of pressures for ES run
        self.TPWT       = []                    # list of weight for pressures in ES runs
        self.TMIN       = 0                     # number of minima per run for DFT
        self.TITR       = 0                     # number of NN-relaxation steps
        self.TNGN       = []                    # list of number of generations per run
        self.TNPP       = 0                     # number of members per population (will be resclaed with weights)
        # internals: set and used inside the code
        #----------------------------------------
        self.NSPC       = 0                     # total number of species in the run
        self.ISIF       = 0                     # VASP relaxation scheme (initiated according to NDIM)
        self.PNUM       = 0                     # number of "pressure values" for data generation
        self.SPC0       = []                    # array of number of atoms per species/population (0 cycle)
        self.SIZ0       = []                    # array of "total atom number" for populations (0 cycle)
        self.SPCN       = []                    # array of number of atoms per species/population (1+ cycle)
        self.SIZN       = []                    # array of "total atom number" for populations (1+ cycle)
        self.NAME       = ""                    # continous string containing symbols of all species of the run  
        self.SATM       = " "                   # string containing all species of the run by their symbols
        self.ZATM       = " "                   # string containing all species of the run by their Z number
        self.RADI       = 0.0                   # average atomic radius of all species in the run
        self.CLUS       = ["","",""]            # cluster commands: submit job, delete job, list of jobs
        self.NGEN       = 1                     # number of generations for ES (cycle 1+)
        self.ngen       = 2                     # number of generations for ES (cycle 0)
        self.NSW0       = 0                     # final high-accuracy DFT runs are static

