#!/usr/bin/env perl

use strict;
use warnings;

# Take as input a directory containing nodesupport reslts with the following filename format:
# tips16_0.5_0400_001_BaseConcTree.result

my %da;
my %depth;
my $datasetID;
my ($tips, $sym, $len, $set, $tree);

my @fileArray;
my $dir=$ARGV[0];
opendir (DIR, $dir) or die $!;
while (my $file = readdir(DIR)) {
  push @fileArray, $file;
}
close(DIR);

foreach my $resultFile (@fileArray) {
    #print "Results File $resultFile\n";
    #print "Result File: $resultFile\n";
    if ($resultFile=~/(tips\d+)_(\d.\d)_(\d+)_(\d+)_(Base|Conc|BaseConc)Tree.result/) {
        $datasetID="$1_$2";
	#print "DatasetID $datasetID\n";
        my $tips = $1;
	my $sym  = $2;
	my $len  = $3;
	my $set  = $4;
	my $tree = $5;

        my @bootstrapMethodList;
        my $node;
        
        open (F, "$dir/$resultFile");
        #print "Open $resultFile\n";
        while (<F>) {
            my $line=$_;
            # #REPLICATE:   2 partialBootstrapTreesConcat.nwk   15
	    if ( ($line=~/\#REPLICATE:/)) {
                my $bootstrapMethod;
                # weightedBootstrapTreesConcat.nwk
                #print "LINE=$line\n";
		$line=~/(bootstrap|paramastrap)sCat.txt/;
                #print "LINE = $line\n METHOD = $1\n\n";
                $bootstrapMethod=$1;
                @bootstrapMethodList=(@bootstrapMethodList,$bootstrapMethod);
            }
            elsif ( !($line=~/\#/)) {
                $node++;
                my @v=split ('\s+', $line);
                my $count=3;
                foreach my $bootstrapMethod (@bootstrapMethodList) {

                    #             0            1        2          $count < 3,4,5 etc..>
                    #           <list>      <depth> <correct> <support from replicate X..>
                    # @v = 1100000000000000    2        1     1.0000 0.9333 0.5333 1.0000

                    $depth{$datasetID}{$v[1]}++;
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'list'}           =$v[0];
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'depth'}          =$v[1];
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'treeCorrect'}    =$v[2]; #Is the tree correct: 0 or 1
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'tree'}{'pp'}                   +=$v[2]; #Number nodes correct in tree
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'tree'}{'nt'}++;
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'bs'}          =$v[$count]; #BS support for the node
		    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'tree'}{bs}                +=$v[$count]; #Sum BS supports in tree    
		    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'tree'}{nnode}++;            
                    #print "Dataset $datasetID\nTree $tree\nBootstrap $bootstrapMethod\nSet $set\nNode $node\n\n";
                    $count++;
	        
                    # Avg Number of Correct Nodes
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'tree'}{'acc'}=$da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'tree'}{'pp'}/$da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'tree'}{'nt'};

                    # Avg Accuracy Score of Node 
                    # If Tree is 0 (incorrect) and BS Method is 0.1, score is 0.9.  
                    # If Tree is 1 (correct) and BS Method is 0.8, score is 0.8.              
                    my $score;
                    if ($da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'treeCorrect'}) {
                        $score = $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'bs'};
                    }
                    else {
                        $score = 1 - $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'bs'};
                    }
                              
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'score'} = $score;               

                }
            }
        }
        close(F);   
    }
}

my %ma;
my %T;
my %bacc;
my $nexp;
foreach my $exp (sort (keys (%da))) {
    $nexp++;

    #print "EXP = $exp\n";

    foreach my $tree (sort (keys (%{$da{$exp}}))) {
        my $j=1;
        my $y;
        my @bootstrapMethods = (sort (keys (%{$da{$exp}{$tree}})));
        my $bootstrapMethod = shift(@bootstrapMethods);

        #print "TREE = $tree\n";

        #print "BS_METHOD = $bootstrapMethod\n";
   
        foreach my $set (sort (keys %{$da{$exp}{$tree}{$bootstrapMethod}})) {

            #print "SET = $set\n"; 
 
            my @ll=keys(%{$da{$exp}{$tree}{$bootstrapMethod}});
            my $nset=$#ll+1;
            $y->[$j][1]=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'bs'};
            $y->[$j][2]=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'acc'};
            $j++;

            $ma{$exp}{$tree}{$bootstrapMethod}{'acc'}+=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'acc'}/$nset;
            $ma{$exp}{$tree}{$bootstrapMethod}{'bs'}+=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'bs'}/$nset;
            
            $ma{all}{$tree}{$bootstrapMethod}{'acc'}+=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'acc'}/$nset;
            $ma{all}{$tree}{$bootstrapMethod}{'bs'}+=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'bs'}/$nset;

            #print "ACC = $da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'acc'}\n";
            
        }
    }
}


my @treeList=('Base');#;,'Conc', 'BaseConc');
my @methodsList=('paramastrap', 'bootstrap');
my @tipsList=('tips16','tips32');#,'tips64');
my @symList=('0.5','1.0','2.0');

print "Table 1.1: Average Topological Accuracy\n";
foreach my $tree (@treeList) {
    print "$tree ";
    foreach my $tips (@tipsList) {
        foreach my $sym (@symList) {
            my $exp="$tips\_$sym";
            my @bootstrapMethods = (sort (keys (%{$da{$exp}{$tree}})));
            my $bootstrapMethod = shift(@bootstrapMethods);
            #print "EXP = $exp\n";
            #print "TREE = $tree\n";
            #print "BS_METHOD = $bootstrapMethod\n";
            printf "%.3f ", $ma{$exp}{$tree}{$bootstrapMethod}{acc};           
         }
    }
    printf "%.3f ",$ma{all}{$tree}{'bootstrap'}{acc}/$nexp;
    print "\n";
}



## ROC Data


my $filename = "ROC.data";
open my $fh, '>', $filename or die;

foreach my $method (@methodsList) {
  print $fh "$method.labels,";
  foreach my $tips (@tipsList) {
    foreach my $sym (@symList) {
      my $exp="$tips\_$sym";
      foreach my $tree (@treeList) {
          foreach my $set (sort (keys %{$da{$exp}{$tree}{$method}})) { 
            foreach my $node (keys %{$da{$exp}{$tree}{$method}{$set}{node}}) {
                  print $fh "$da{$exp}{$tree}{$method}{$set}{'node'}{$node}{'treeCorrect'},";
          }
        }
      }
    }
  }
  print $fh "\b \b\n";
}


foreach my $method (@methodsList){
  print $fh "$method.predictions,";
  foreach my $tips (@tipsList) {
    foreach my $sym (@symList) {
      my $exp="$tips\_$sym";
      foreach my $tree (@treeList) {
          foreach my $set (sort (keys %{$da{$exp}{$tree}{$method}})) {
              foreach my $node (keys %{$da{$exp}{$tree}{$method}{$set}{node}}) {
                  print $fh "$da{$exp}{$tree}{$method}{$set}{'node'}{$node}{'bs'},";
          }
        }
      }
    }
  }
  print $fh "\b \b\n";
}


my $filename_mcc = "MCC.data";
open my $fh_bmcc, '>', $filename_mcc or die;

foreach my $method (@methodsList){
  foreach my $tips (@tipsList) {
    foreach my $sym (@symList) {
      my $exp="$tips\_$sym";
      foreach my $tree (@treeList) {
          my $mcc_string = '';
          $mcc_string = $mcc_string."${method}\_${exp},";
          foreach my $set (sort (keys %{$da{$exp}{$tree}{$method}})) {
              my $y;
              my $j=0;
              foreach my $node (keys %{$da{$exp}{$tree}{$method}{$set}{node}}) {
                 $y->[$j][0]=$da{$exp}{$tree}{$method}{$set}{'node'}{$node}{'treeCorrect'};
                 $y->[$j][1]=$da{$exp}{$tree}{$method}{$set}{'node'}{$node}{'bs'};
                  
                 $j++;
              }
              my $mcc=list2bmcc_alternative($y);
              $mcc_string = $mcc_string."$mcc,";
          }
          chop($mcc_string);
          $mcc_string = $mcc_string."\n";
          print $fh_bmcc $mcc_string;
       }
    }
  }
} 

# Best Mathews Correlation Coefficent

sub list2bmcc
  {
     my ($x)=@_;

     my $n= scalar(@{$x});
     my ($bmcc,$pbs,$bbs);
     my ($tp,$tn,$fp,$fn);
     my @array;
     for (my $i=0; $i<$n; $i++)
       {
         $array[$i][0]=0+$x->[$i][0];
         $array[$i][1]=0+$x->[$i][1];

         print "$i\t$array[$i][0]\n";
         print "$i\t$array[$i][1]\n\n";

         $tp+=0+$array[$i][0];
         $fp+=1-$array[$i][0];
       }
     @array=sort { $a->[1] <=> $b->[1]} @array;


     $bmcc=-999;
     $pbs=-1;
     $bbs=-1;
     for (my $i=0; $i<$n; $i++)
       {
         my $bs=$array[$i][1];
         my $pp=$array[$i][0];
         my $pn=1-$pp;
         if ($pbs!=-1 && $bs!=$pbs)
           {
             my $matN=($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn);
             my $mcc=($matN==0)?-999:(($tp*$tn)-($fp*$fn))/sqrt(($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn));
             if ($mcc >$bmcc)
               {
                 $bmcc=$mcc; $bbs=$pbs;
               }
           }
         $tp-=$pp;
         $tn+=$pn;
         $fp-=$pn;
         $fn+=$pp;
         $pbs=$bs;
       }
     my $matN=($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn);
     my $mcc=($matN==0)?-999:(($tp*$tn)-($fp*$fn))/sqrt(($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn));
     if ($mcc >$bmcc){$bmcc=$mcc; $bbs=$pbs;}
     return $bmcc;
   }

sub list2bmcc_alternative {

     my ($x)=@_;

     my $n= scalar(@{$x});
     my @array;
     for (my $i=0; $i<$n; $i++)
       {
         $array[$i][0]=0+$x->[$i][0];
         $array[$i][1]=0+$x->[$i][1];

         #print "$i\t$array[$i][0]\n";
         #print "$i\t$array[$i][1]\n\n";

       }
     @array=sort { $a->[1] <=> $b->[1]} @array;
     
     my $bmcc=-999;
     my $mcc=0;
     my ($tp,$tn,$fp,$fn);
     for (my $k=0; $k<=$n; $k++) { # for every level we draw the line
       $mcc=0;
       $tp=0;$tn=0;$fp=0;$fn=0;
       for (my $j=0; $j<$n; $j++) { # for every node we assign value
         if ($array[$j][0] == 0 && $j >= $k ) {$fp=$fp+1}
         elsif ($array[$j][0] == 0 && $j < $k ) {$tn=$tn+1} 
         elsif ($array[$j][0] == 1 && $j >= $k ) {$tp=$tp+1} 
         elsif ($array[$j][0] == 1 && $j < $k ) {$fn=$fn+1} 
         else {print "Error |$array[$j][0]| $k $j\n"; die}
       }

     print "$k = tp $tp, tn $tn, fp $fp, fn $fn \n";


     if ( ($tn == 0 && $fp == 0 && $fn ==0) || ($tp == 0 && $fp == 0 && $fn ==0 ))  {
       $mcc=1;
     }

     elsif ( ($tp+$fp)==0 || ($tp+$fn)==0 || ($tn+$fp) == 0 || ($tn+$fn) ==0 ) {
       $mcc=0;
     }

     else { 
       $mcc = (($tp*$tn)-($fp*$fn))/sqrt(($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn));
     }  
    
     if ($mcc >$bmcc){$bmcc=$mcc;}

     }      
     print "BMCC = $bmcc\n\n";
     return $bmcc;
}


## Accuracy 
