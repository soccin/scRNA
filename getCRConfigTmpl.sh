#!/bin/bash

CELLRANGER=$1

OUTS=$(find $CELLRANGER -type d | egrep "/outs$" | xargs realpath | sort -V)

echo "inputs:"
for oi in $OUTS; do
    si=$(echo $oi | sed 's/.outs$//' | tr '/' '\n' | tail -1)
    echo "- dir:" $oi
    echo "  sid:" $si
done
echo "genome: [mm10|human?]"
