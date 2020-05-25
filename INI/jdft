#!/bin/bash
#SBATCH -p long,main
#SBATCH -n 16
#SBATCH -N 1-1
#SBATCH --mem-per-cpu=MaxMemPerCPU

cd $SLURM_SUBMIT_DIR

# path to VASP executable
EXE=/opt/ohpc/pub/apps/VASP/vasp.5.4.4/bin/vasp_std


#============================ DO NOT CHANGE FROM HERE


python maisenet_kmsh.py KKKK

date        >> id
mpirun $EXE  > out
date        >> id
mv POSCAR    POSCAR.0
mv CONTCAR   CONTCAR.0
mv OUTCAR    OUTCAR.0
mv OSZICAR   OSZICAR.0

rm IB* EI* PC* XD* CH* WA* vasp* DO* maisenet_kmsh.py
