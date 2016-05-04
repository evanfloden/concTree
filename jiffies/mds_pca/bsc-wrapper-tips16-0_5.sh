#!/bin/bash
#BSUB -oo %J.tips16_0.5.out
#BSUB -eo %J.tips16_0.5.err
#BSUB -J tips16_0.5
#BSUB -q debug
#BSUB -W 1:00
#BSUB -x
#BSUB -n 8
mds_alignments.sh /home/bic22/bic22913/projects/concMSA-nf/jiffies/tips16_0.5_0400_001.fa 
