#!/bin/bash
#SBATCH -N1
##SBATCH --time 00:05:00
#SBATCH --time 02:00:00
#SBATCH --exclusive
#SBATCH --export="NONE"
##SBATCH -c 64
#SBATCH -c 128
##SBATCH --partition ndl

set -x

export DR_HOOK=1
export DR_HOOK_OPT=prof
export DR_HOOK_IGNORE_SIGNALS=-1
export DR_NVTX=1

ulimit -s unlimited
export OMP_STACK_SIZE=32G

cd $SLURM_SUBMIT_DIR

\rm *.grb

export SLURM_EXPORT_ENV=ALL
export MPIAUTOCONFIG=mpiauto.PGI.conf
#export MPIAUTOCONFIG=mpiauto.DDT.conf
#export OMP_NUM_THREADS=64
export OMP_NUM_THREADS=128

#/opt/softs/mpiauto/mpiauto --prefix-command ./nvprof.sh --nouse-slurm-mpi\
# --verbose -np 8 --wrap --wrap-stdeo -- ./spcm.x --case ../t0031l015-008mpi --write-grib-1 --write-grib-2 --stat-gp
#/opt/softs/bin/ja
#exit

#/opt/softs/mpiauto/mpiauto --prefix-command ./nvprof.sh --nouse-slurm-mpi\
# --verbose -np 8 --wrap --wrap-stdeo -- ./spcm.x --case ../t0107l070-008mpi --stat-gp
#/opt/softs/bin/ja
#exit

#/opt/softs/mpiauto/mpiauto --prefix-command ./nvprof.sh --nouse-slurm-mpi\
# --verbose -np 8 -openmp 1 --wrap --wrap-stdeo -- ./spcm.x --case ../t0149l105-008mpi --write-grib-1 --write-grib-2 --stat-gp --stat-sp
#/opt/softs/bin/ja
#exit
#type nvprof

#/opt/softs/mpiauto/mpiauto \
#  --prefix-command ./nvprof.sh --nouse-slurm-mpi \
#  --verbose -np 1 --wrap --wrap-stdeo -- \
#  ./spcm.x --case ../t0031l015-001mpi --stat-gp
#
#/opt/softs/bin/ja
#exit

/opt/softs/mpiauto/mpiauto \
   --nouse-slurm-mpi \
  --verbose -np 1 --wrap --wrap-stdeo -- \
  ./spcm.x --case ../t0499l105-001mpi --stat-gp

#/opt/softs/mpiauto/mpiauto \
#  --prefix-command ./nsys.sh --nouse-slurm-mpi \
#  --verbose -np 1 --wrap --wrap-stdeo -- \
#  ./spcm.x --case ../t0031l015-001mpi --stat-gp


#/opt/softs/mpiauto/mpiauto --nouse-slurm-mpi --verbose -np 8 -openmp 1 --wrap --wrap-stdeo -- ./spcm.x --case t0149l105-008mpi --write-grib-1 --write-grib-2 --stat-gp --stat-sp

/opt/softs/bin/ja
