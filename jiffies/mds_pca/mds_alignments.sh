#!/bin/bash
HHLIB="/home/bic22/bic22913/projects/concMSA-nf/lib/hhsuite-2.0.16-linux-x86_64"
PERL5LIB="$PERL5LIB:/home/bic22/bic22913/projects/concMSA-nf/lib/hhsuite-2.0.16-linux-x86_64/scripts:/home/bic22/bic22913/projects/concMSA-nf/lib64/perl5"

rm -rf guidance2 a3m base_alignment.aln alternativeMSA hhm_db tmp

perl ~/software/guidance.v2.01/www/Guidance/guidance.pl \
        --seqFile $1\
        --msaProgram CLUSTALW \
        --seqType aa \
        --outDir $PWD/guidance2 \
        --clustalw clustalw2

    cp guidance2/MSA.CLUSTALW.aln.With_Names base_alignment.aln

    mkdir $PWD/alternativeMSA
    tar -zxf guidance2/MSA.CLUSTALW.Guidance2_AlternativeMSA.tar.gz -C alternativeMSA
perl hhsearch_mds.pl /home/bic22/bic22913/projects/concMSA-nf/jiffies/alternativeMSA /home/bic22/bic22913/projects/concMSA-nf/jiffies/a3m 20
