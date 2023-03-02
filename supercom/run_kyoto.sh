#!/bin/bash
#QSUB -q gr10310b
#QSUB -ug gr10310
#QSUB -W 24:00
#QSUB -A p=160:t=1:c=1:m=3413M
#QSUB -M masaichi2425@gmail.com
#QSUB -m be
module load anaconda3/5.3.1
source activate Masa 
make -f kaityo_makefile
mpirun -np 160 ./cps joblist.sh
