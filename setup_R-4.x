module unload R
module load gcc
RPATH=/juno/work/bic/socci/opt/common/CentOS_7/R/R-4.1.0
export PATH=$RPATH/bin:$PATH
R_VERSION=$(R --version | head -1 | awk '{print $3}')
export R_LIBS=/home/socci/lib/R/CentOS${MAJOROS}/$R_VERSION
