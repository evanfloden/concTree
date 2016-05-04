tips=(16 32 64)
syms=('0.5' '1.0' '2.0')

for tip in ${tips[@]}
do 
  for sym in ${syms[@]}
  do
    rm -rf tips${tip}_${sym}

    mkdir tips${tip}_${sym}

    for i in {001..100}
    do

     esl-reformat afa /users/cn/jchang/projects/2010-03_ConM-Coffee/results/2013-03-26_bench-simulation/data/tips${tip}/asymmetric_${sym}/${i}/seqs/${i}.0400.fa > tips${tip}_${sym}/tips${tip}_${sym}_${i}.0400.afa

     t_coffee -other_pg seq_reformat -in /users/cn/jchang/projects/2010-03_ConM-Coffee/results/2013-03-26_bench-simulation/data/tips${tip}/asymmetric_${sym}/${i}/seqs/${i}.0400.fa -output fasta_seq > tips${tip}_${sym}/tips${tip}_${sym}_${i}.0400.fa


    done
  done
done

