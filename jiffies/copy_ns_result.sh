#!/bin/bash

# $PWD/results/tip16_0.5/tips16_0.5_001.0400/nodeSupport/nodeSupportForBaseTree.result

files="../null/tips*_*_*.0400/nodeSupport/*.result"
regex="(tips[0-9]+)_([0-9].[0-9])_([0-9]+).0400/nodeSupport/nodeSupportFor(Base|Conc|BaseConc)Tree.result"
for f in $files
do
    [[ $f =~ $regex ]]
    tips="${BASH_REMATCH[1]}"
    sym="${BASH_REMATCH[2]}"
    set_id="${BASH_REMATCH[3]}"
    tree="${BASH_REMATCH[4]}" 
    cp $f "../nodeSupport/${tips}_${sym}_0400_${set_id}_${tree}Tree.result"
done
