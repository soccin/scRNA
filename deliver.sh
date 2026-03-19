#/bin/bash

workuser_tag() {
  local result
  result=$(pwd | sed -n 's|.*Work/Users/\([^/]*\)/\([^/]*\)/\([^/]*\).*|\1/\2/\3|Ip')
  if [[ -n "$result" ]]; then
      echo "$result" | tr '[:upper:]' '[:lower:]' | sed 's|proj_b|Proj_B|;s|proj_|Proj_|'
  fi
}

SDIR="$( cd "$( dirname "$0" )" && pwd )"

DEFAULT_RESDIR=/ifs/rtsia01/bic/results/$(workuser_tag)/r_001

USE_DEFAULT=0
if [ "$1" = "-d" ]; then
    USE_DEFAULT=1
    shift
fi

if [ "$USE_DEFAULT" = "1" ]; then
    RESDIR=$DEFAULT_RESDIR
    mkdir -p "$RESDIR"
elif [ "$#" = "1" ]; then
    RESDIR=$1
else
    echo
    echo "usage: deliver.sh [-d] [RESDIR]"
    echo
    echo "    -d          use default RESDIR and create it"
    echo "    RESDIR=$DEFAULT_RESDIR"
    echo
    exit
fi

echo \$RESDIR=$(realpath $RESDIR)

mkdir -p $RESDIR/seurat
mkdir -p $RESDIR/seurat/docs

rsync -avP --exclude="._*" results/stage* $RESDIR/seurat
rsync -avP --exclude="._*" *_TblClusterMarkers_*xlsx $RESDIR/seurat
rsync -avP --exclude="._*" *_ClusterMarkerPathways_*xlsx $RESDIR/seurat
cp scRNA/docs/*pdf $RESDIR/seurat/docs
