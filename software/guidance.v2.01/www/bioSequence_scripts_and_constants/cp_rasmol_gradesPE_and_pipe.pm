#!/usr/bin/perl -w

package cp_rasmol_gradesPE_and_pipe;

use Fcntl ':flock'; # import LOCK_* constants

#use strict;
#************************************************************************************************

my %consurf_rasmol_colors=("1"=>"[16,200,209]","2"=>"[140,255,255]","3"=>"[215,255,255]","4"=>"[234,255,255]",
                   "5"=>"[255,255,255]","6"=>"[252,237,244]","7"=>"[250,201,222]","8"=>"[240,125,171]",
                   "9"=>"[160,37,96]","10"=>"[255,255,150]");

my %consurf_colors_chimera=("0"=>" 255 255 150","1"=>" 16 200 209","2"=>" 140 255 255","3"=>" 215 255 255","4"=>" 234 255 255","5"=>" 255 255 255","6"=>" 252 237 244","7"=>" 250 201 222","8"=>" 240 125 171","9"=>" 160 37 96");
my $bayesInterval=3;
my %ColorScale = (0 => 9, 1 => 8, 2 => 7, 3 => 6, 4 => 5, 5 => 4, 6 => 3, 7 => 2, 8 => 1);
my %tr_aa;
    $tr_aa{LYS}="K";$tr_aa{ARG}="R";$tr_aa{HIS}="H";$tr_aa{ASP}="D";$tr_aa{GLU}="E";
    $tr_aa{TYR}="Y";$tr_aa{TRP}="W";$tr_aa{SER}="S";$tr_aa{THR}="T";$tr_aa{PHE}="F";
    $tr_aa{LEU}="L";$tr_aa{ILE}="I";$tr_aa{MET}="M";$tr_aa{CYS}="C";$tr_aa{ASN}="N";
    $tr_aa{GLN}="Q";$tr_aa{ALA}="A";$tr_aa{VAL}="V";$tr_aa{PRO}="P";$tr_aa{GLY}="G";

my $CHIMERA_SAVING_FIGURE_LINK="http://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/print.html";
#############
#subroutines#
#############

#************************************************************************************************
# read the rate4site output into an array.
# calculates layers according to the max and min grades from the output
# $Output : a pointer to array ; assigns each position in the MSA (array index) with a grade (arracy value)

sub assign_colors_according_to_r4s_layers{	
    
	my ($rate4site_filename, $Output) = @_;
    
    unless (open RATE4SITE, $rate4site_filename){
        return ("assign_colors_according_to_r4s_layers : can't open $rate4site_filename","PANIC");}
    my $line;
    my $i = 0;
    while (<RATE4SITE>) {
        $line = $_;
        chomp $line;
		# baysean
        if ($line =~ /^\s+(\d+)\s+(\w)\s+(\S+)\s+\[\s*(\S+),\s*(\S+)\]\s+\S+\s+(\d+)\/(\d+)/){
                $Output->[$i]{POS} = $1;
                $Output->[$i]{SEQ} = $2;
                $Output->[$i]{GRADE} = $3;
                $Output->[$i]{INTERVALLOW} = $4;
                $Output->[$i]{INTERVALHIGH} = $5;
                $Output->[$i]{MSA_NUM}=$6;
                $Output->[$i]{MSA_DENUM}=$7;
                $i++;
        }
		# Maximum likelihood
        elsif($line=~m/(\d+)\s+(\w)\s+(\S+)\s+(\d+)\/(\d+)/){
            $Output->[$i]{POS} = $1;
            $Output->[$i]{SEQ} = $2;
            $Output->[$i]{GRADE} = $3;
            $Output->[$i]{INTERVALLOW} = $3;
            $Output->[$i]{INTERVALHIGH} = $3;
            $Output->[$i]{MSA_NUM}=$4;
            $Output->[$i]{MSA_DENUM}=$5;
            $i++;
        }
    }
    close RATE4SITE;
	
    my $element;
    my $max_cons = $Output->[0]{GRADE};
    my $ConsColorUnity; #unity of conservation to be colored
    foreach $element (@$Output){
        if ($$element{GRADE} < $max_cons) {$max_cons = $$element{GRADE};}
    }
    $ConsColorUnity = $max_cons / 4.5 * -1; 
    if ($max_cons !~ /^\-/){$ConsColorUnity = $max_cons;}        
        
#calculates the grades for each color

    my $NoLayers = 9;
    my @ColorLayers;
    for (my $i = 0; $i <= $NoLayers; $i++) {
        $ColorLayers[$i] = $max_cons + ($i * $ConsColorUnity);
    }
#gives the color to the interval

    my $Count = 0;

    foreach my $element (@$Output) {
        
        for (my $i = 0; $i <= $#ColorLayers; $i++){
            if( ($i==$#ColorLayers) and !exists $$element{INTERVALLOWCOLOR}){
                $$element{INTERVALLOWCOLOR} = 8;
            }
            elsif ($$element{INTERVALLOW} >= $ColorLayers[$i] and $$element{INTERVALLOW} < $ColorLayers[$i + 1]) {
                $$element{INTERVALLOWCOLOR} =$i;
            } 
            elsif ( ($$element{INTERVALLOW} < $ColorLayers[$i]) and !exists $$element{INTERVALLOWCOLOR}){
                $$element{INTERVALLOWCOLOR} = 0;
            } 
            if (($i == $#ColorLayers)  and !exists $$element{INTERVALHIGHCOLOR}){
                $$element{INTERVALHIGHCOLOR} = 8;
            }
            elsif ($$element{INTERVALHIGH} >= $ColorLayers[$i] and $$element{INTERVALHIGH} < $ColorLayers[$i + 1]) {
                $$element{INTERVALHIGHCOLOR} =$i;
            } 
            elsif ( ($$element{INTERVALHIGH} < $ColorLayers[$i]) and !exists $$element{INTERVALHIGHCOLOR}){
                $$element{INTERVALHIGHCOLOR} = 0;
            }
        } # END FOR
    } # END FOREACH
	
#give the color for each position based on the grades	
	# match the colors to the grades
    foreach my $element (@$Output){ 
        for (my $i = 0; $i <= $#ColorLayers; $i++) {
            if ($i == $#ColorLayers) {
                $$element{COLOR} = $ColorScale{$i-1};
            }
            elsif ($$element{GRADE} >= $ColorLayers[$i] && $$element{GRADE} < $ColorLayers[$i + 1]) {
                $$element{COLOR} = $ColorScale{$i};         
                last;
            }            
        }
        if ((($$element{INTERVALHIGHCOLOR}-$$element{INTERVALLOWCOLOR})>$bayesInterval ) or ($$element{MSA_NUM} <= 5)){
            $$element{ISD}=1;			
        }
        else{$$element{ISD}=0;}		
    }
	return ("OK");
} 
#************************************************************************************************
#printing the the ConSurf gradesPE file
sub create_gradesPE_ConSurf{
	my ($Output,$ref_match, $ref_residue_freq, $no_isd_residue_color, $isd_residue_color, $gradesPE_file) = @_;
    my ($seq3d_grades_isd, $seq3d_grades);
    # open file
    unless (open PE, ">$gradesPE_file" ){
        return ("create_gradesPE_ConSurf : can't open '$gradesPE_file'","PANIC");}
    print PE "\t Amino Acid Conservation Scores\n";
    print PE "\t===============================\n\n";
    print PE "- POS: The position of the AA in the SEQRES derived sequence.\n";
    print PE "- SEQ: The SEQRES derived sequence in one letter code.\n";
    print PE "- 3LATOM: The ATOM derived sequence in three letter code, including the AA's positions as they appear in the PDB file and the chain identifier.\n";
    print PE "- SCORE: The normalized conservation scores.\n";
    print PE "- COLOR: The color scale representing the conservation scores (9 - conserved, 1 - variable).\n";
    print PE "- CONFIDENCE INTERVAL: When using the bayesian method for calculating rates, a confidence interval is assigned to each of the inferred evolutionary conservation scores.\n"; 
    print PE "- CONFIDENCE INTERVAL COLORS: When using the bayesian method for calculating rates. The color scale representing the lower and upper bounds of the confidence interval.\n"; 
    print PE "- MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position.\n";
    print PE "- RESIDUE VARIETY: The residues variety at each position of the multiple sequence alignment.\n\n";
    print PE " POS\t SEQ\t    3LATOM\tSCORE\t\tCOLOR\tCONFIDENCE INTERVAL\tCONFIDENCE INTERVAL COLORS\tMSA DATA\tRESIDUE VARIETY\n";
    print PE "    \t    \t        \t(normalized)\t        \t               \n";
    foreach my $elem (@$Output){
		my $pos = $$elem{POS};
		my $var = "";
		my $atom_3L=$$ref_match{$pos};
		my $score = $$elem{COLOR};
        
		printf (PE "%4d", $pos);
        printf (PE "\t%4s", "$$elem{SEQ}");	        
        printf (PE "\t%10s", $atom_3L);
        printf (PE "\t%6.3f", "$$elem{GRADE}");
        if($$elem{ISD}==1){
            printf (PE "\t\t%3d", "$$elem{COLOR}");
            printf (PE "%1s", "*");
        }
        else{printf (PE "\t\t%3d", "$$elem{COLOR}");}
        printf (PE "\t%6.3f", "$$elem{INTERVALLOW}");
        printf (PE "%1s", ",");
        printf (PE "%6.3f", "$$elem{INTERVALHIGH}");
        printf (PE "\t\t\t%5d", "$ColorScale{$$elem{INTERVALLOWCOLOR}}");
        printf (PE "%1s", ",");
        printf (PE "%1d\t\t", "$ColorScale{$$elem{INTERVALHIGHCOLOR}}");
        printf (PE "\t%8s", "$$elem{MSA_NUM}\/$$elem{MSA_DENUM}");
        for my $_aa (keys %{$ref_residue_freq->{($pos)}}){
            $var.= "$_aa,";
        }
        chop($var) if ($var =~ /,$/);
        print PE "\t$var\n";		
		# the amino-acid in that position, must be part of the residue variety in this column
		if ($var !~ /$$elem{SEQ}/i){
			close PE;
			return ("create_gradesPE_ConSurf : in position $pos, the amino-acid ".$$elem{SEQ}." does not match the residue variety: $var.","PANIC");}
		#printing the residue to the rasmol script
        #assigning grades to $seq3d strings
        if($atom_3L !~ m/\-/){
            $atom_3L =~ m/(.+):/;
            $atom_3L = $1;
            if($score=~ m/(\d)/){
                my $color = $1;
                push @{ $no_isd_residue_color->[$color] }, $atom_3L;
                #if($score=~ m/\*/){
				if ($$elem{ISD}==1){
                    push @{ $isd_residue_color->[10] }, $atom_3L;
                    $seq3d_grades_isd.="0";
                }
                else{
                    push @{ $isd_residue_color->[$color] }, $atom_3L;                        
                    $seq3d_grades_isd.="$color";
                }
                $seq3d_grades.="$color";
            }
        }
        else{
            $seq3d_grades_isd.=".";
            $seq3d_grades.=".";
        }
    }
    print PE "\n\n*Below the confidence cut-off - The calculations for this site were performed on less than 6 non-gaped homologue sequences,\nor the confidence interval for the estimated score is equal to- or larger than- 4 color grades.\n";
    close(PE);
	return ("OK",$seq3d_grades_isd, $seq3d_grades);
}
#************************************************************************************************
#printing the the ConSurf gradesPE file
sub create_gradesPE_ConSeq{
	my ($Output, $ref_residue_freq, $ref_Solv_Acc_Pred, $gradesPE_file) = @_;
    # open file
    unless (open PE, ">$gradesPE_file" ){
        return ("create_gradesPE_Conseq : can't open '$gradesPE_file'","PANIC");}
    print PE "\t Amino Acid Conservation Scores\n";
    print PE "\t===============================\n\n";
    print PE "- POS: The position of the AA in the SEQRES derived sequence.\n";
    print PE "- SEQ: The SEQRES derived sequence in one letter code.\n";
    print PE "- SCORE: The normalized conservation scores.\n";
    print PE "- COLOR: The color scale representing the conservation scores (9 - conserved, 1 - variable).\n";
    print PE "- CONFIDENCE INTERVAL: When using the bayesian method for calculating rates, a confidence interval is assigned to each of the inferred evolutionary conservation scores.\n"; 
    print PE "- CONFIDENCE INTERVAL COLORS: When using the bayesian method for calculating rates. The color scale representing the lower and upper bounds of the confidence interval.\n"; 
    print PE "- B/E: Burried (b) or Exposed (e) residue.\n"; 
    print PE "- FUNCTION: functional (f) or structural (s) residue (f - highly conserved and exposed, s - highly conserved and burried).\n"; 
    print PE "- MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position.\n";
    print PE "- RESIDUE VARIETY: The residues variety at each position of the multiple sequence alignment.\n\n";
    print PE " POS\t SEQ\tSCORE\t\tCOLOR\tCONFIDENCE INTERVAL\tCONFIDENCE INTERVAL COLORS\tB\/E\tFUNCTION\tMSA DATA\tRESIDUE VARIETY\n";
    print PE "    \t    \t(normalized)\t        \t               \n";
    foreach my $elem (@$Output){
		my $pos = $$elem{POS};
		my $var = "";
		my $Solv_Acc_Pred=$$ref_Solv_Acc_Pred{$pos};
		my $score = $$elem{COLOR};
        
		printf (PE "%4d", $pos);
        printf (PE "\t%4s", "$$elem{SEQ}");	        
        
        printf (PE "\t%6.3f", "$$elem{GRADE}");
        if($$elem{ISD}==1){
            printf (PE "\t\t%3d", "$$elem{COLOR}");
            printf (PE "%1s", "*");
        }
        else{printf (PE "\t\t%3d", "$$elem{COLOR}");}
        printf (PE "\t%6.3f", "$$elem{INTERVALLOW}");
        printf (PE "%1s", ",");
        printf (PE "%6.3f", "$$elem{INTERVALHIGH}");
        printf (PE "\t\t\t%5d", "$ColorScale{$$elem{INTERVALLOWCOLOR}}");
        printf (PE "%1s", ",");
        printf (PE "%1d\t\t", "$ColorScale{$$elem{INTERVALHIGHCOLOR}}");
        printf (PE "\t%3s", $Solv_Acc_Pred);
	## FUNCT/STRUCT COL
	if ($Solv_Acc_Pred eq "e"){
           if ($$elem{COLOR} == 9 || $$elem{COLOR} == 8 ){
	       printf (PE "\t%8s", "f");
          }
           else {
               printf (PE "\t%8s", " ");
          }
       }
       else {
          if ($$elem{COLOR} == 9){
               printf (PE "\t%8s", "s");
          }
           else {
               printf (PE "\t%8s", " ");
          }
       }
	
	
	printf (PE "\t%8s", "$$elem{MSA_NUM}\/$$elem{MSA_DENUM}");
        for my $_aa (keys %{$ref_residue_freq->{($pos)}}){
            $var.= "$_aa,";
        }
        chop($var) if ($var =~ /,$/);
        print PE "\t$var\n";		
		# the amino-acid in that position, must be part of the residue variety in this column
		if ($var !~ /$$elem{SEQ}/i){
			close PE;
			return ("create_gradesPE_ConSeq : in position $pos, the amino-acid ".$$elem{SEQ}." does not match the residue variety: $var.","PANIC");}
		#printing the residue to the rasmol script
        #assigning grades to $seq3d strings
    }
    print PE "\n\n*Below the confidence cut-off - The calculations for this site were performed on less than 6 non-gaped homologue sequences,\nor the confidence interval for the estimated score is equal to- or larger than- 4 color grades.\n";
    close(PE);
	return ("OK");
}
####################################################################
#printing the the ConSurf gradesPE file for Nucleic Acid
####################################################################
sub create_gradesPE_ConSeq_Nuc{
	my ($Output, $ref_residue_freq, $gradesPE_file) = @_;
    # open file
    unless (open PE, ">$gradesPE_file" ){
        return ("create_gradesPE_Conseq_Nuc : can't open '$gradesPE_file'","PANIC");}
    print PE "\t Nucleic Acid Conservation Scores\n";
    print PE "\t===============================\n\n";
    print PE "- POS: The position of the Nucleic Acid on the sequence.\n";
    print PE "- SEQ: The Nucleic Acid\n";
    print PE "- SCORE: The normalized conservation scores.\n";
    print PE "- COLOR: The color scale representing the conservation scores (9 - conserved, 1 - variable).\n";
    print PE "- CONFIDENCE INTERVAL: When using the bayesian method for calculating rates, a confidence interval is assigned to each of the inferred evolutionary conservation scores.\n"; 
    print PE "- CONFIDENCE INTERVAL COLORS: When using the bayesian method for calculating rates. The color scale representing the lower and upper bounds of the confidence interval.\n"; 
    print PE "- MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position.\n";
    print PE "- RESIDUE VARIETY: The residues variety at each position of the multiple sequence alignment.\n\n";
    print PE " POS\t SEQ\tSCORE\t\tCOLOR\tCONFIDENCE INTERVAL\tCONFIDENCE INTERVAL COLORS\tMSA DATA\tRESIDUE VARIETY\n";
    print PE "    \t    \t(normalized)\t        \t               \n";
    foreach my $elem (@$Output){
		my $pos = $$elem{POS};
		my $var = "";
		my $score = $$elem{COLOR};
        
		printf (PE "%4d", $pos);
        printf (PE "\t%4s", "$$elem{SEQ}");	        
        
        printf (PE "\t%6.3f", "$$elem{GRADE}");
        if($$elem{ISD}==1){
            printf (PE "\t\t%3d", "$$elem{COLOR}");
            printf (PE "%1s", "*");
        }
        else{printf (PE "\t\t%3d", "$$elem{COLOR}");}
        printf (PE "\t%6.3f", "$$elem{INTERVALLOW}");
        printf (PE "%1s", ",");
        printf (PE "%6.3f", "$$elem{INTERVALHIGH}");
        printf (PE "\t\t\t%5d", "$ColorScale{$$elem{INTERVALLOWCOLOR}}");
        printf (PE "%1s", ",");
        printf (PE "%1d\t\t", "$ColorScale{$$elem{INTERVALHIGHCOLOR}}");
	
	printf (PE "\t%8s", "$$elem{MSA_NUM}\/$$elem{MSA_DENUM}");
        for my $_aa (keys %{$ref_residue_freq->{($pos)}}){
            $var.= "$_aa,";
        }
        chop($var) if ($var =~ /,$/);
        print PE "\t$var\n";		
		# the nucleotide in that position, must be part of the residue variety in this column
		if ($var !~ /$$elem{SEQ}/i){
			close PE;
			return ("create_gradesPE_ConSeq : in position $pos, the nucleotide ".$$elem{SEQ}." does not match the residue variety: $var.","PANIC");}
		#printing the residue to the rasmol script
        #assigning grades to $seq3d strings
    }
    print PE "\n\n*Below the confidence cut-off - The calculations for this site were performed on less than 6 non-gaped homologue sequences,\nor the confidence interval for the estimated score is equal to- or larger than- 4 color grades.\n";
    close(PE);
	return ("OK");
}
#************************************************************************************************
#matches the position in the seqres/msa sequence to the position in the pdb
sub match_seqres_pdb{
	
	my ($seqres_atom_aln, $atom_position, $chain, $ref_fas2pdb) = @_;
	my (@seqres, @atoms, $length_of_seqres, $length_of_atom);
	my $pdbseq="";
	my $query_seq="";
	
	#creating arrays containing each sequences
	unless(open ALN, $seqres_atom_aln){return ("match_seqres_pdb : Could not open the file $seqres_atom_aln for reading $!", "PANIC");}
	flock ALN, LOCK_EX;
	while(<ALN>){
		if($_=~ m/^ATOM_\S+\s+(\S+)/){
			$pdbseq=$pdbseq.$1;
		}
		elsif($_=~ m/^SEQRES_\S+\s+(\S+)/){
			$query_seq=$query_seq.$1;
		}
		elsif($_=~ m/^MSA_\S+\s+(\S+)/){
			$query_seq=$query_seq.$1;
		}
		elsif($_=~ m/^QUERY\s+(\S+)/){
			$query_seq=$query_seq.$1;
		}
	}
	flock ALN,LOCK_UN;
	close(ALN);
	#arrays that conatin the sequences from the pairwise for each sequneces (including gaps).
	#the sequnces start from position 1 in the array!
	$pdbseq="X".$pdbseq;
	$query_seq="X".$query_seq;
	@atoms=split("",$pdbseq);
	@seqres=split("",$query_seq);
	
	#creating the hash that matches the position in the ATOM fasta to its position
	#in the pdb file and also the fasta ATOM position to the correct residue
	unless(open MATCH, $atom_position){return ("match_seqres_pdb : Could not open the file $atom_position for reading $!", "PANIC");}
	
	flock MATCH, LOCK_EX;
	my %match_ATOM;
	my %res_ATOM;
	while(<MATCH>){
		$_=~ m/(\w{3})\s+(\d+)\s+(\S+)/;
		my $res=$1;
		my $fas_atom=$2;
		my $pdb_atom=$3;
		$match_ATOM{$fas_atom}=$res.$pdb_atom.":".$chain;
	}
	flock MATCH,LOCK_UN;
	close MATCH;
    #creating a hash in which the key is the position in the aln (i.e the
    #position in the @seqres) and the value is the correct position in the fas.
    my $pos_count=0;
    my %aln_fas_seqres;
    for(my $n=1;$n<@seqres+0;$n++){
#        if($seqres[$n] !~ m/\-/){ #HAIM
        if(($seqres[$n] !~ m/\-/) and ($seqres[$n] !~ m/X/)){ 
            $pos_count++;
            $aln_fas_seqres{$n}=$pos_count;
        }
    }
    $length_of_seqres = $pos_count;
    #creating a hash in which the key is the position in the aln (i.e the
    #position in the @atoms) and the value is the correct position in the fas.
    $pos_count=0;
    my %aln_fas_atoms;
    for(my $n=1;$n<@atoms+0;$n++){
#        if($atoms[$n] !~ m/-/){ #HAIM
	  if(($atoms[$n] !~ m/-/) and ($atoms[$n] !~ m/X/)){
            $pos_count++;
            $aln_fas_atoms{$n}=$pos_count;
        }
    }	
    $length_of_atom = $pos_count;    
    for(my $i=1;$i<@seqres+0;$i++){
		my $fas_pos=$aln_fas_seqres{$i};
        if($seqres[$i]!~m/-/){           
#            if($atoms[$i] eq "-"){ #HAIM
            if(($atoms[$i] eq "-") or ($atoms[$i] eq "X")){
                $ref_fas2pdb->{$fas_pos}="-";
            }
            else{
                 my $match=$match_ATOM{$aln_fas_atoms{$i}};
                $ref_fas2pdb->{$fas_pos}=$match;
            }
        }
    }
    return ("OK",$length_of_seqres, $length_of_atom);
}
#************************************************************************************************

#adds the 3LATOM colomn to the new gradesPE file. also creates a rasmol script.
# its input is a reference to a hash where key is the position in the SEQRES sequence, and the value is that 3LATOM according to the PDB sequence
sub add_pdb_gradesPE{
    my ($pre_gradesPE, $outPE, $ref_match, $no_isd_residue_color, $isd_residue_color, $ref_residue_freq)=@_;
	my ($seq3d_grades_isd, $seq3d_grades, $aa_position);
    my ($grade, $score, $high, $low, $gradel, $gradeh, $msa_m, $msa_d, $var);
    my $seq = "";
    unless (open GRADESPE, $pre_gradesPE){return ("add_pdb_gradesPE : Could not open the file $pre_gradesPE for reading $!", "PANIC");}
    unless(open PE, ">$outPE"){return ("add_pdb_gradesPE : Could not open the file $outPE for writing $!", "PANIC");}
    while(<GRADESPE>){
        if($_!~ m/\s+\d+\s+\w\s+atom_res/){
            print PE "$_";    
            next;
        }
        else{
            $_=~m/^\s+(\d+)\s+/;
            my $pos=$1;
            my $atom_3L=$$ref_match{$pos};
			$atom_3L=~ m/(\w{3})/;
			my $aa=$1;			
            $_=~ m/\s+(\d+)\s+(\w)\s+atom_res\s+(\S+)\s+(\S+)\s+(\S+),\s*(\S+)\s+(\S+),\s*(\S+)\s+(\d+)\/(\d+)/;
            $pos=$1; $seq = $2;
            $grade=$3; $score=$4; $high=$5;$low=$6;$gradel=$7;
            $gradeh=$8; $msa_m=$9; $msa_d=$10; $var="";
                
            printf (PE "%4d", "$pos");
            printf (PE "\t%4s", "$seq");
            printf (PE "\t%10s", "$atom_3L");
            printf (PE "\t%6.3f", "$grade");
            if($score=~m/\*/){
                printf (PE "\t\t%3d", "$score");
                printf (PE "%1s", "\*");
            }
            else{printf (PE "\t\t%3d", "$score");}
            printf (PE "\t%6.3f", "$high");
            printf (PE "%1s", ",");
            printf (PE "%6.3f", "$low");
            printf (PE "\t\t\t%5d", "$gradel");
            printf (PE "%1s", ",");
            printf (PE "%1d\t\t", "$gradeh");
            printf (PE "\t%8s", "$msa_m\/$msa_d");
            for my $_aa (keys %{$ref_residue_freq->{$pos}}){
                $var.= "$_aa,";
            }
            
            chop($var) if ($var =~ /,$/);
            print PE "\t$var\n";
            #printing the residue to the rasmol script
            #assigning grades to $seq3d strings
            if($atom_3L !~ m/\-/){
                $atom_3L =~ m/(.+):/;
                $atom_3L = $1;
                if($score=~ m/(\d)/){
                    my $color = $1;
                    push @{ $no_isd_residue_color->[$color] }, $atom_3L;
                    if($score=~ m/\*/){
                        push @{ $isd_residue_color->[10] }, $atom_3L;
                        $seq3d_grades_isd.="0";
                    }
                    else{
                        push @{ $isd_residue_color->[$color] }, $atom_3L;                        
                        $seq3d_grades_isd.="$color";
                    }
                    $seq3d_grades.="$color";
                }
            }
            else{
                $seq3d_grades_isd.=".";
                $seq3d_grades.=".";
            }
        }
    }
    close PE;
	close GRADESPE;
	return ("OK",$seq3d_grades_isd, $seq3d_grades);
}
#************************************************************************************************
sub print_rasmol{
    my ($out_file,$rasmol_isd,$ref_colors_array,$chain, $proteopedia) = @_;
    # print out new format of rasmol
    unless (open OUT, ">".$out_file) {return ("print_rasmol : Could not open the file $out_file for writing $!", "PANIC");}
    print OUT "select all\ncolor [200,200,200]\n\n" if $proteopedia ne "yes";
    for (my $i=$#{$ref_colors_array}+1; $i>0;$i--){
        next if ($i==10 and $rasmol_isd eq "no");        
        if ($#{$ref_colors_array->[$i]}>-1){
            print OUT print_selected(\@{ $ref_colors_array->[$i] }, "no");
            print OUT "\nselect selected and :$chain\n";
            print OUT "color $consurf_rasmol_colors{$i}\nspacefill\n";
            print OUT "define CON".($i)." selected\n\n";            
        }
    }
#    print OUT
    close OUT;
	return ("OK");
}
#************************************************************************************************
sub print_selected{
    my $arr_ref = shift;
    my $print_for_pipe = shift;
    my @all_aa = @$arr_ref;
    my $string = "";
    my $total_text = "";
    if ($print_for_pipe eq "yes"){
        $string = "! select ";}
    else{
        $string = "select ";}
    my $total_length = length $string;
	my $aa_num=@all_aa;
	print "all_AA ($aa_num):",join(",",@all_aa),"\n";
	if ($aa_num>0)
	{
		foreach (@all_aa){        
			my $_aa = $_;
			$total_length+=length $_aa;
			#print "aa is: $_aa total is: $total_length\n";
			if ($total_length >80){
				if ($string =~ m/, $/){
					chop $string; chop $string; 
            }
				$total_text .= $string."\n";            
				if ($print_for_pipe eq "yes"){
					$string="! select selected or $_aa, ";}
				else{
					$string="select selected or $_aa, ";}            
            $total_length = length $string;
			}
			else{
				$string .="$_aa, ";
				$total_length+=2; # for the ',' and ' '
			}
			$_aa = "";
		}
	}
	else
	{
		if ($print_for_pipe eq "yes"){
			$string="! select none";}
		else{
			$string="select none";}
		$total_text.=$string;
	}
    if ($string =~ m/, $/){
        chop $string; chop $string;
        $total_text.= $string;
    }
	print "$total_text\n";
    return $total_text;
}
#************************************************************************************************
#input: $input_pdb_file, $ref_return_arr
#output: ("OK"/"ERR",$header_line,$title_lines,$compnd_lines)

sub extract_data_from_pdb{
	my $input_pdb_file = shift;
	my @return_arr = ();
	unless (open PDB, $input_pdb_file){
		$return_arr[0] = "extract_data_from_pdb : Could not open the file $input_pdb_file $!";
		$return_arr[1] = "PANIC";
	}	
	else{
		flock PDB, 2;
		while (<PDB>){
			if(/^HEADER/)   {
				$_ =~ s/\s+$//;
				$return_arr[1] = $_;
			}
			elsif (/^TITLE\s+/){
				$_ =~ s/\s+$//;
				$_ =~ /^TITLE\s+\d*\s(.*)\s*/;
				$return_arr[2].=$1." ";
			}
			 elsif (/^COMPND\s+/) {
				$_ =~ s/\s+$//;
				$_ =~ /^COMPND\s+\d*\s(.*)/;
				$return_arr[3].=$1." ";
			}
			elsif (/^SOURCE/ || /^KEYWDS/ || /^AUTHOR/ || /^SEQRES/ || /^ATOM/) {# no nead to go over all the pdb
				last;
			}
		}
		flock PDB, 8;
		close PDB;
		$return_arr[0] = "OK";
	}
	return @return_arr;
}
#************************************************************************************************
# take a string aaaaaaa and returns it in this format: ! "aaa" +\n! "aa";\n
sub design_string_for_pipe{
    my $string_to_format = shift;
    my $part = $string_to_format;
    my $newPart = "";
    while (length $part > 73){
        $newPart .= "! \"".(substr $part, 0, 73)."\" +\n";
        $part = substr $part, 73;
    }
    $newPart .= "! \"".$part."\" ;";
    return $newPart;
}
#************************************************************************************************
sub design_string_with_spaces_for_pipe{    
    my $part = shift;    
    my @words = split /\s+/, $part;
    my $newPart = "! \"".$words[0];
    $part = "";
    for (my $i=1; $i<@words; $i++){
        # if adding another word to the string will yeild a too long string - we cut it.
        if ((length($words[$i]) +1+ length($newPart)) > 76){
            $part .= $newPart."\" +\n" ;
            $newPart="! \"".$words[$i];
        }
        else{
            $newPart.=" $words[$i]";
        }
    }
    $part .= $newPart."\" ;";
    return $part;
}
#************************************************************************************************
# input : 1 - a path to a directory where sub-directories hold every chain,
#  2 - a reference to hash that will store chains and their identical chains as a string
# please note: if chain A and B are identical, $ref_identical_chains->{"A"} will hold only "B" and $ref_identical_chains->{"B"} will hold only "A"
sub find_identical_chains{
    my ($all_identical,$ref_identical_chains) = @_;
    my %chain_seq = ();    
    my $seq = "";	
    
    opendir PROTEIN, $protein_path;
    while (my $chain = readdir PROTEIN){
        next if ($chain !~ /^[A-Za-z]$/);
        my $seq_file = $protein_path.$chain."/seq.fas";
        if (-e $seq_file){
            open SEQ, $seq_file;
            <SEQ>;
            $seq = <SEQ>;
            close SEQ;
            next if ($seq !~ /\w+/);
            if (exists $chain_seq{$seq}){
                if (length $chain_seq{$seq} > 1){
                    my @chains = split //, $chain_seq{$seq};
                    foreach (@chains){
                        $ref_identical_chains->{$chain}.=$_;
                        $ref_identical_chains->{$_}.=$chain;
                    }                    
                }
                else{
                    $ref_identical_chains->{$chain}.=$chain_seq{$seq};
                    $ref_identical_chains->{$chain_seq{$seq}} .= $chain;
                }
                $chain_seq{$seq}.=$chain;                
            }
            else{
                $chain_seq{$seq} = $chain;
            }
        }
    }
    closedir PROTEIN;  
}
#************************************************************************************************
# the routine creates 2 hashes, according to information it reads from a MSA
# for each position in the query sequence :
# 1. it collects all the residues that are aligned to this positions. it also counts the number of time each residue appeared in that position in the MSA
# 2. it counts the total number of residues which aligned to this position
sub read_residue_variety
{
	my ($msa_filename, $msa_ref_sequence, $ref_residue_frequency, $ref_position_totalAA) = @_;
	my $line;
	my $position = 1;
	my $num_of_seqs = 0;
	my $flag=0;
	my $data_length;
	my @elements_to_remove = ();
	my $ret_seq = "";

	# open msa file
	unless (open MSA, $msa_filename) {return ("Could not open the file $msa_filename for reading $!", "PANIC");}
	flock MSA, LOCK_EX;
	# read header lines
	$line = <MSA>;
	$line = <MSA>;
	while ($line = <MSA>)
	{
	   if ($line =~ /^(.+) +([A-Za-z-]+)$/)
	    {
			$num_of_seqs++ if ($flag==0); #count the number of sequences, only in the first block
            if($flag==1){
                #inc position and read next (blank) line
                $position += $data_length;
                $flag=2;
            }
	    	my $seq_name = $1;
	        my $seq_data = $2;
	        $seq_name =~ s/\s+$//; # remove spaces
	        $data_length = length($seq_data);
	        for (my $position_offset = 0;$position_offset < length($seq_data);$position_offset++)
	        {
	            # get residue at this position
	            my $current_aa = substr($seq_data,$position_offset,1);
	            
	            # if this is the query seq and no aa is found we need to save
	            # this position in order to erase it in the end
				if ($seq_name eq $msa_ref_sequence and $current_aa !~ /^[ABCDEFGHIKLMNPQRSTVWY]$/){
						push @elements_to_remove,$position+$position_offset;
						$ret_seq .= "in seq $seq_name ignored $current_aa in position ".($position+$position_offset)." " if ($current_aa !~ /^[-Xx]$/);
				}				
				elsif ($current_aa !~ /^-$/){
					# save residue value
					if (!exists $ref_residue_frequency->{$position+$position_offset}){
						$ref_residue_frequency->{$position+$position_offset}->{$current_aa} = 1;
					}
					elsif (!exists $ref_residue_frequency->{$position+$position_offset}->{$current_aa}){
						$ref_residue_frequency->{$position+$position_offset}->{$current_aa} = 1;
					}
					else{                
						$ref_residue_frequency->{$position+$position_offset}->{$current_aa}++;
					}
					$ref_position_totalAA->{$position+$position_offset}++;
				}
	        }
	    }
	    elsif($line =~ /^\s+/){
			$flag=1;
	    }
	}
	flock MSA, LOCK_UN;
	close MSA;
	
	# remove all positions where the query seq was skipped
	my $index = pop(@elements_to_remove);
	while (defined($index))
	{
        residues_shift_left($index, $ref_residue_frequency);
        residues_shift_left($index, $ref_position_totalAA);
		$index = pop(@elements_to_remove);
	}
	# If there was an unexpected character in the query sequence in the MSA, report it and exit
	if ($ret_seq ne "")
		{return ($ret_seq, 'PANIC');}
	else {return ("OK", $num_of_seqs);}
}
#************************************************************************************************

sub residues_shift_left
{
    my ($start_position, $ref_residue_frequency) = @_;
    my $index = $start_position;
    
    while (exists($ref_residue_frequency->{$index+1}))
    {
        $ref_residue_frequency->{$index} = $ref_residue_frequency->{$index+1};
        $index++;
    }
    delete($ref_residue_frequency->{$index});    
}
#************************************************************************************************
# for each position, calculate the % for each residue which is found in the MSA in that position
# the output is written in "Cvs" format, meaning there are ',' signs for each tab carachter
sub print_precentage{
    my ($ref_residue_freq,$ref_position_totalAA, $out_file, $ref_ConSurfGrades) = @_;  # HAIM ADD ref_ConSurfGrades as optional input  
    my @aa_arr = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","OTHER");
    my ($val, $aa_found, $total, @aa_in_position, $other, $other_val);
    unless (open OUT, ">".$out_file) {return ("print_precentage : Could not open the file $out_file for writing $!", 'PANIC');}
	print OUT "\"The table details the residue variety in % for each position in the query sequence.\"\n\"Each column shows the % for that amino-acid, found in position ('pos') in the MSA.\"\n\"In case there are residues which are not a standard amino-acid in the MSA, they are represented under column 'OTHER'\"\n\n";
    print OUT "pos";
    print OUT ",$_" foreach (@aa_arr);
    print OUT ",MAX AA,ConSurf Grade\n";
    # order the lines according to AA position in the sequence
    for my $position (sort {$a <=> $b} keys %$ref_residue_freq){
        $total=0; $aa_found = "no"; $other="";$other_val=0;
        my $max_precent=0;my $max_AA="";
	print OUT $position;
        # the total number of animo acids found in the MSA for that position
        $total += $ref_position_totalAA->{$position};
        #For each position , sort the amino acids variety
        @aa_in_position = sort keys %{$ref_residue_freq->{$position}};
        my $j=0;
        $aa = $aa_in_position[$j];
        $val = ($ref_residue_freq->{$position}->{$aa})/$total*100;
        # in order to print a table, we go by the sorted array
        my $i=0; 
	while ($i<@aa_arr){
	    # if there are non standart aa, we calculate its total value seperately and add it under column OTHER
	    while ($aa !~ /^[KRHDEYWSTFLIMCNQAVPG]$/i and $j < @aa_in_position){
		$other_val+=$val;
		$j+=1;
		if ($j < @aa_in_position){
		    $aa = $aa_in_position[$j];
		    $val = ($ref_residue_freq->{$position}->{$aa})/$total*100;
		}
	    }
	    if ($aa_arr[$i] eq 'OTHER' and $other_val!=0){
		$other = $other_val;
		print OUT ",$other";
	    }			
#            elsif ($aa_arr[$i] eq $aa){				
	    elsif ($aa_arr[$i] eq uc ($aa)){
		if ($val>$max_precent)
		{
		    $max_precent=$val;
		    $max_AA=$aa;
		}
		elsif ($val==$max_precent)
		{
		    $max_AA.=$aa;
		}	
		print OUT ",".sprintf("%.3f",$val);
		$j+=1;
		if ($j < @aa_in_position){
		    $aa = $aa_in_position[$j];
		    $val = ($ref_residue_freq->{$position}->{$aa})/$total*100;
		    
		}
	    }			
	    else{
		print OUT ",";
	    }
	    $i++;
	}
        print OUT ",$max_AA ".sprintf("%.3f",$max_precent);
	if (defined $ref_ConSurfGrades->[($position-1)]) # HAIM ADD 2/6/2013
	{
	    my $GradesPE_elem=$ref_ConSurfGrades->[($position-1)];
	    print OUT ",$$GradesPE_elem{COLOR}";
	    if($$GradesPE_elem{ISD}==1){
		print OUT "*";
	    }
	}
	print OUT "\n";
	
    }
    close OUT;
    return ("OK");
}
#************************************************************************************************
# for each position, calculate the % for each nucleotide which is found in the MSA in that position
# the output is written in "Cvs" format, meaning there are ',' signs for each tab carachter
# SHOULD BE REWIRITTEN - THE LAST IS NOT PRINTED
sub print_precentage_nuc{
    my ($ref_nucleotide_freq,$ref_position_totalNuc, $out_file) = @_;    
    my @nuc_arr = ("A","C","G","T","U","OTHER");
    my ($val, $nuc_found, $total, @nuc_in_position, $other, $other_val);
    unless (open OUT, ">".$out_file) {return ("print_precentage_nuc : Could not open the file $out_file for writing $!", 'PANIC');}
	print OUT "\"The table details the nucleic acid variety in % for each position in the query sequence.\"\n\"Each column shows the % for that nucleic-acid, found in position ('pos') in the MSA.\"\n\"In case there are nucleotides  which are not a standard nucleic-acid in the MSA, they are represented under column 'OTHER'\"\n\n";
    print OUT "pos";
    print OUT ",$_" foreach (@nuc_arr);
    print OUT "\n";
    # order the lines according to AA position in the sequence
    for my $position (sort {$a <=> $b} keys %$ref_nucleotide_freq){
        $total=0; $nuc_found = "no"; $other="";$other_val=0;
        print OUT $position;
        # the total number of nucleic acids found in the MSA for that position
        $total += $ref_position_totalNuc->{$position};
        #For each position , sort the nucleic acids variety
        @nuc_in_position = sort keys %{$ref_nucleotide_freq->{$position}};
	print "POS:$position\tNUC_IN_POS:",join(",",@nuc_in_position),"\n";
        my $j=0;
        $nuc = $nuc_in_position[$j];
	print "NUC:$nuc\n";
	$val = ($ref_nucleotide_freq->{$position}->{$nuc})/$total*100;
#	print "VAL:$val\n";
        # in order to print a table, we go by the sorted array
        my $i=0; 
		while ($i<@nuc_arr){
			# if there are non standart aa, we calculate its total value seperately and add it under column OTHER
			while ($nuc !~ /^[ACTGU]$/i and $j < @nuc_in_position){
				$other_val+=$val;
				$j+=1;
				if ($j < @nuc_in_position){
					$nuc = $nuc_in_position[$j];
					$val = ($ref_nucleotide_freq->{$position}->{$nuc})/$total*100;
				}
			}
			if ($nuc_arr[$i] eq 'OTHER' and $other_val!=0){
				$other = $other_val;
				print OUT ",$other";
			}			
			elsif ($nuc_arr[$i] eq uc($nuc)){				
			  print OUT ",$val";
			  $j+=1;
			  if ($j < @nuc_in_position){
			    $nuc = $nuc_in_position[$j];
			    $val = ($ref_nucleotide_freq->{$position}->{$nuc})/$total*100;
			    
			  }
			}			
			else{
			  print OUT ",";
			}
			print "POS\t$position\tI:$i\t$nuc_arr[$i]\tVAL:$val\n";
			$i++;
		      }
        print OUT "\n";
      }
    close OUT;
    return ("OK");
  }
#************************************************************************************************
# for each position, calculate the % for each residue which is found in the MSA in that position
# the values in the line might not sum to exactly 100%, as we round the value to be written with
# 2 digits after the '.'
sub print_precentage_text{
    my ($ref_residue_freq,$ref_position_totalAA, $out_file) = @_;    
    my @aa_arr = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","OTHER");
    my ($val, $aa_found, $total, @aa_in_position, $other, $other_val);
    unless (open OUT, ">".$out_file) {return ("Could not open the file $out_file for writing $!", 'PANIC');}
	print OUT "The table details the residue variety in % for each position in the query sequence.\nEach column shows the % for that amino-acid, found in position (\"pos\") in the MSA.\nIn case there are residues which are not a standard amino-acid, they are represented under column \"OTHER\"\n\n";
    print OUT "pos  |  ";
    print OUT "   -$_-  " foreach (@aa_arr);
    print OUT "\n" ;
    print OUT "---------" foreach (@aa_arr);
    print OUT "\n";
    # order the lines according to AA position in the sequence
    for my $position (sort {$a <=> $b} keys %$ref_residue_freq){
        $total=0; $aa_found = "no"; $other="";$other_val=0;
        printf OUT '%4s', $position;
        print OUT " | ";
        # the total number of animo acids found in the MSA for that position
        $total += $ref_position_totalAA->{$position};
        #For each position , sort the amino acids variety
        @aa_in_position = sort keys %{$ref_residue_freq->{$position}};
        my $j=0;
        $aa = $aa_in_position[$j];
        $val = ($ref_residue_freq->{$position}->{$aa})/$total*100;
        # in order to print a table, we go by the sorted array
        my $i=0; 
		while ($i<@aa_arr){
			# if there are non standart aa, we calculate its total value seperately and add it under column OTHER
			while ($aa !~ /^[KRHDEYWSTFLIMCNQAVPG]$/i and $j < @aa_in_position){
				$other_val+=$val;
				$j+=1;
				if ($j < @aa_in_position){
					$aa = $aa_in_position[$j];
					$val = ($ref_residue_freq->{$position}->{$aa})/$total*100;
				}
			}
			if ($aa_arr[$i] eq 'OTHER' and $other_val!=0){
				$other = sprintf '%*2$.2f', $other_val;
				print OUT "     ".$other;
			}			
#            elsif ($aa_arr[$i] eq $aa){				
	     elsif ($aa_arr[$i] eq uc($aa)){		
                printf OUT '%*2$.2f', $val, 8;
                $j+=1;
                if ($j < @aa_in_position){
                    $aa = $aa_in_position[$j];
                    $val = ($ref_residue_freq->{$position}->{$aa})/$total*100;
					
                }
            }			
            else{
                print OUT "        ";
            }
			$i++;
        }
        print OUT "\n";
    }
    close OUT;
	return ("OK");
}
#************************************************************************************************
sub round{
    my $in = shift;
    if ($in =~ /(\d+)\.(\d)(\d)(\d)\d*/){
        $x = $1; $y = $2; $z = $3; $w = $4;
        if (less_than(5,$w)){return $x.".".$y.$z;}
        else{
            if(less_than(9,$z)){return $x.".".$y.($z+1);}
            else{
                if(less_than(9,$y)){return $x.".".($y+1)."0";}
                else{return ($x+1).".00";}
            }
        }
    }
    else{
        return $in;
    }
}
sub less_than{
    my $_val = shift;
    my $in = shift;
    if ($in>=0 and $in<$_val) {return 1;}
    else{return 0;}
}
#************************************************************************************************
# creating part of the pipe file, which contains all the non-unique information.
# each chain will use this file to construct the final pdb_pipe file, to be viewed with FGiJ
sub create_part_of_pipe{	
    my ($pipe_file, $num_of_seq_in_msa, $IN_db_file, $seq3d_grades_isd, $seq3d_grades, $length_of_seqres, $length_of_atom, $ref_isd_residue_color, $ref_no_isd_residue_color) = @_;
	
    #############################################
    # collect data for printing the "pipe" file #
    #############################################
    # database where sequences collected from:
    unless (open DB, $IN_db_file){return ("cannot open the file $IN_db_file for reading $!", 'PANIC');}
    <DB> =~ /Database: (.+)/;
    if ($1 =~ /sprot/i or $1 =~ /swissprot/i){
        $db = "SWISS-PROT";}
    else{
        $db = "UNIPROT";}        
    close DB;    
    
    # design the seq3d to be printted out to the pipe file
    $seq3d_grades_isd = rasmol_gradesPE_and_pipe::design_string_for_pipe($seq3d_grades_isd);
    $seq3d_grades = rasmol_gradesPE_and_pipe::design_string_for_pipe($seq3d_grades);
    
    # creating the frequencies array which corresponds the number of residues in each grade
    my ($consurf_grade_freqs_isd, $consurf_grade_freqs) = freq_array($ref_isd_residue_color, $ref_no_isd_residue_color);
    
    ##########################
    # write to the pipe file #
    ##########################
    unless (open PIPE, ">".$pipe_file){return ("cannot open the file $pipe_file for writing $!", 'PANIC');}
    
    print PIPE <<EndOfPipe;
! consurf_psi_blast_e_value = 0.001;
! consurf_psi_blast_database = "$db";
! consurf_psi_blast_iterations = 3;
! consurf_max_seqs = 300;
! consurf_alignment = "MUSCLE";
! consurf_method = "bayesian";
! consurf_substitution_model =  "JTT";
!
! consurf_seqres_length = $length_of_seqres;
! consurf_atom_seq_length = $length_of_atom;
! consurf_unique_seqs = $num_of_seq_in_msa;
! consurf_grade_freqs_isd = $consurf_grade_freqs_isd;
! consurf_grade_freqs = $consurf_grade_freqs;
!
! seq3d_grades_isd =
$seq3d_grades_isd
!
! seq3d_grades = 
$seq3d_grades
!
!
!! ====== CONTROL PANEL OPTIONS SECTION ======
!js.init
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! pipe_title_enlarged = false;
! pipe_background_color = "white";
!
!! Specify the custom consurf control panel
!!
! pipe_cp1 = "consurf/consurf.htm";
!
!! If you want the frontispiece to be reset every time you enter this
!! page, use false. If this is a one-page presentation (no contents)
!! and you want to be able to return from QuickViews without resetting
!! the view, use true.
!!
! frontispiece_conditional_on_return = true;
!
!! Open the command input slot/message box to 30% of window height.
!!
! pipe_show_commands = true;
! pipe_show_commands_pct = 30;
!
!! Don't show the PiPE presentation controls in the lower left frame.
!!
! pipe_hide_controls = true;
!
!! Hide development viewing mode links at the bottom of the control panel.
!!
! pipe_tech_info = false; 
!
!! pipe_start_spinning = true; // default is PE's Preference setting.
!! top.nonStopSpin = true; // default: spinning stops after 3 min.
!!
!! ====== COLORS SECTION ======
!!
!color color_carbon C8C8C8
!color color_sulfur FFC832
!
!! Ten ConSurf color grades follow:
!!
!color color_grade0 FFFF96 insufficient data yellow
!color color_grade1 10C8D1 turquoise variable
!color color_grade2 8CFFFF
!color color_grade3 D7FFFF
!color color_grade4 EAFFFF
!color color_grade5 FFFFFF
!color color_grade6 FCEDF4
!color color_grade7 FAC9DE
!color color_grade8 F07DAB
!color color_grade9 A02560 burgundy conserved
!
!
!! ====== SCRIPTS SECTION ======
!!----------------------------------------
!!
!spt #name=select_and_chain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!!----------------------------------------
!!
!spt \#name=view01
! \@spt consurf_view_isd
!
!!----------------------------------------
!!
!spt \#name=hide_all
! restrict none
! ssbonds off
! hbonds off
! dots off
! list \* delete
!
!!----------------------------------------
!! common_spt uses CPK carbon gray (or phosphorus yellow) for backbones.
!!
!spt \#name=common_spt
! \@spt hide_all
! select all
! color [xC8C8C8] \# rasmol/chime carbon gray
! select nucleic
! color [xFFA500] \# phosphorus orange
! select hetero
! color cpk
! select not hetero
! backbone 0.4
! javascript top.water=0
! 
! ssbonds 0.3
! set ssbonds backbone
! color ssbonds \@color_sulfur
! 
! select hetero and not water
! spacefill 0.45
! wireframe 0.15
! dots 50
! 
! select protein
! center selected
! 
!!----------------------------------------
!!
!spt \#name=consurf_view_isd
! \@spt common_spt
! \@for \$=0, 9
! \@spt select_isd_grade\$
! \@spt select_and_chain
! color \@color_grade\$
! spacefill
! \@endfor
! zoom 115
!
!!----------------------------------------
EndOfPipe
    
	my $lineToPrint = "";
    for (my $i = 9; $i>0; $i--){
        print PIPE "!!\n!spt \#name=select_isd_grade$i\n!\n";
		$lineToPrint = cp_rasmol_gradesPE_and_pipe::print_selected(\@{$ref_isd_residue_color->[$i]}, "yes");
        print PIPE $lineToPrint."\n" if ($lineToPrint =~ /select/);
        print PIPE "!\n!\n!!----------------------------------------\n";
		$lineToPrint = "";
    }
    print PIPE "!!\n!spt \#name=select_isd_grade0\n";
	$lineToPrint = cp_rasmol_gradesPE_and_pipe::print_selected(\@{$ref_isd_residue_color->[10]}, "yes");
	print PIPE $lineToPrint."\n" if ($lineToPrint =~ /select/);
    print PIPE "!\n!\n!!----------------------------------------\n";
    for (my $i = 9; $i>0; $i--){
        print PIPE "!!\n!spt \#name=select_grade$i\n!\n";
		$lineToPrint = cp_rasmol_gradesPE_and_pipe::print_selected(\@{$ref_no_isd_residue_color->[$i]}, "yes");
        print PIPE $lineToPrint."\n" if ($lineToPrint =~ /select/);
        print PIPE "!\n!\n!!----------------------------------------\n";
		$lineToPrint = "";
    }
    print PIPE "!!\n!spt \#name=select_grade0\n! select none\n!!\n";
    print PIPE "!! ====== END OF CONSURF PiPE BLOCK ======\n";    
    close PIPE;
	return ("OK", $db);
}
#************************************************************************************************
# design the frequencies array
sub freq_array{
	my ($isd_residue_color, $no_isd_residue_color) = @_;
	my ($consurf_grade_freqs_isd, $consurf_grade_freqs);
# the insufficient data should be the first in this array
    $consurf_grade_freqs_isd = "Array(";
    $consurf_grade_freqs_isd .=$#{$isd_residue_color->[10]}+1;
    for (my $i=1; $i<10; $i++){
        $consurf_grade_freqs_isd .=",";
        $consurf_grade_freqs_isd .= $#{$isd_residue_color->[$i]}+1;    
    }
    $consurf_grade_freqs_isd .=")";
    $consurf_grade_freqs = "Array(0";
    print "\n";
    for (my $i=1; $i<10; $i++){
        $consurf_grade_freqs .=",";
        $consurf_grade_freqs .= $#{$no_isd_residue_color->[$i]}+1;    
    }
    $consurf_grade_freqs .=")";
	return ($consurf_grade_freqs_isd, $consurf_grade_freqs);
}

#************************************************************************************************
# creating part of the pipe file, which contains all the non-unique information.
# each chain will use this file to construct the final pdb_pipe file, to be viewed with FGiJ
sub create_part_of_pipe_new{	
    my ($pipe_file, $unique_seqs, $db, $seq3d_grades_isd, $seq3d_grades, $length_of_seqres, $length_of_atom, $ref_isd_residue_color, $ref_no_isd_residue_color, $E_score,$iterations,$max_num_homol,$MSAprogram,$algorithm,$matrix) = @_;

	
#    #############################################
#    # collect data for printing the "pipe" file #
#    #############################################
#    # database where sequences collected from:
#    unless (open DB, $IN_db_file){return ("cannot open the file $IN_db_file for reading $!", 'PANIC');}
#    <DB> =~ /Database: (.+)/;
#    if ($1 =~ /sprot/i or $1 =~ /swissprot/i){
#        $db = "SWISS-PROT";}
#    else{
#        $db = "UNIPROT";}        
#    close DB;    
#    
    # design the seq3d to be printted out to the pipe file
    $seq3d_grades_isd = rasmol_gradesPE_and_pipe::design_string_for_pipe($seq3d_grades_isd);
    $seq3d_grades = rasmol_gradesPE_and_pipe::design_string_for_pipe($seq3d_grades);
    
    # creating the frequencies array which corresponds the number of residues in each grade
    my ($consurf_grade_freqs_isd, $consurf_grade_freqs) = freq_array($ref_isd_residue_color, $ref_no_isd_residue_color);
    
    # Taking Care of Strings
    if ($max_num_homol eq "all"){$max_num_homol="\"$max_num_homol\""};
    	
    ##########################
    # write to the pipe file #
    ##########################
    unless (open PIPE, ">".$pipe_file){return ("cannot open the file $pipe_file for writing $!", 'PANIC');}
    
    print PIPE <<EndOfPipe;
! consurf_psi_blast_e_value = $E_score;
! consurf_psi_blast_database = "$db";
! consurf_psi_blast_iterations = $iterations;
! consurf_max_seqs = $max_num_homol;
! consurf_alignment = "$MSAprogram";
! consurf_method = "$algorithm";
! consurf_substitution_model =  "$matrix";
!
! consurf_seqres_length = $length_of_seqres;
! consurf_atom_seq_length = $length_of_atom;
! consurf_unique_seqs = $unique_seqs;
! consurf_grade_freqs_isd = $consurf_grade_freqs_isd;
! consurf_grade_freqs = $consurf_grade_freqs;
!
! seq3d_grades_isd =
$seq3d_grades_isd
!
! seq3d_grades = 
$seq3d_grades
!
!
!! ====== CONTROL PANEL OPTIONS SECTION ======
!js.init
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! pipe_title_enlarged = false;
! pipe_background_color = "white";
!
!! Specify the custom consurf control panel
!!
! pipe_cp1 = "consurf/consurf.htm";
!
!! If you want the frontispiece to be reset every time you enter this
!! page, use false. If this is a one-page presentation (no contents)
!! and you want to be able to return from QuickViews without resetting
!! the view, use true.
!!
! frontispiece_conditional_on_return = true;
!
!! Open the command input slot/message box to 30% of window height.
!!
! pipe_show_commands = true;
! pipe_show_commands_pct = 30;
!
!! Don't show the PiPE presentation controls in the lower left frame.
!!
! pipe_hide_controls = true;
!
!! Hide development viewing mode links at the bottom of the control panel.
!!
! pipe_tech_info = false; 
!
!! pipe_start_spinning = true; // default is PE's Preference setting.
!! top.nonStopSpin = true; // default: spinning stops after 3 min.
!!
!! ====== COLORS SECTION ======
!!
!color color_carbon C8C8C8
!color color_sulfur FFC832
!
!! Ten ConSurf color grades follow:
!!
!color color_grade0 FFFF96 insufficient data yellow
!color color_grade1 10C8D1 turquoise variable
!color color_grade2 8CFFFF
!color color_grade3 D7FFFF
!color color_grade4 EAFFFF
!color color_grade5 FFFFFF
!color color_grade6 FCEDF4
!color color_grade7 FAC9DE
!color color_grade8 F07DAB
!color color_grade9 A02560 burgundy conserved
!
!
!! ====== SCRIPTS SECTION ======
!!----------------------------------------
!!
!spt #name=select_and_chain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!!----------------------------------------
!!
!spt \#name=view01
! \@spt consurf_view_isd
!
!!----------------------------------------
!!
!spt \#name=hide_all
! restrict none
! ssbonds off
! hbonds off
! dots off
! list \* delete
!
!!----------------------------------------
!! common_spt uses CPK carbon gray (or phosphorus yellow) for backbones.
!!
!spt \#name=common_spt
! \@spt hide_all
! select all
! color [xC8C8C8] \# rasmol/chime carbon gray
! select nucleic
! color [xFFA500] \# phosphorus orange
! select hetero
! color cpk
! select not hetero
! backbone 0.4
! javascript top.water=0
! 
! ssbonds 0.3
! set ssbonds backbone
! color ssbonds \@color_sulfur
! 
! select hetero and not water
! spacefill 0.45
! wireframe 0.15
! dots 50
! 
! select protein
! center selected
! 
!!----------------------------------------
!!
!spt \#name=consurf_view_isd
! \@spt common_spt
! \@for \$=0, 9
! \@spt select_isd_grade\$
! \@spt select_and_chain
! color \@color_grade\$
! spacefill
! \@endfor
! zoom 115
!
!!----------------------------------------
EndOfPipe
    
	my $lineToPrint = "";
    for (my $i = 9; $i>0; $i--){
        print PIPE "!!\n!spt \#name=select_isd_grade$i\n!\n";
		$lineToPrint = cp_rasmol_gradesPE_and_pipe::print_selected(\@{$ref_isd_residue_color->[$i]}, "yes");
        print PIPE $lineToPrint."\n" if ($lineToPrint =~ /select/);
        print PIPE "!\n!\n!!----------------------------------------\n";
		$lineToPrint = "";
    }
    print PIPE "!!\n!spt \#name=select_isd_grade0\n";
	$lineToPrint = cp_rasmol_gradesPE_and_pipe::print_selected(\@{$ref_isd_residue_color->[10]}, "yes");
	print PIPE $lineToPrint."\n" if ($lineToPrint =~ /select/);
    print PIPE "!\n!\n!!----------------------------------------\n";
    for (my $i = 9; $i>0; $i--){
        print PIPE "!!\n!spt \#name=select_grade$i\n!\n";
		$lineToPrint = cp_rasmol_gradesPE_and_pipe::print_selected(\@{$ref_no_isd_residue_color->[$i]}, "yes");
        print PIPE $lineToPrint."\n" if ($lineToPrint =~ /select/);
        print PIPE "!\n!\n!!----------------------------------------\n";
		$lineToPrint = "";
    }
    print PIPE "!!\n!spt \#name=select_grade0\n! select none\n!!\n";
    print PIPE "!! ====== END OF CONSURF PiPE BLOCK ======\n";    
    close PIPE;
	return ("OK");
}
#************************************************************************************************
# Create the pipe file for FGiJ
sub create_consurf_pipe_new{
	my ($results_dir, $IN_pdb_id_capital, $chain, $ref_header_title, $final_pipe_file, $identical_chains,$partOfPipe,$current_dir,$run_number,$msa_filename,$query_name_in_msa,$tree_filename,$submission_time,$completion_time,$run_date) = @_;
	
        if (!defined $submission_time){$submission_time="";} # OPTIONAL ARGUMENT
	if (!defined $tree_filename){$tree_filename="";}     #OPTIONAL ARGUMENT
	if (!defined $run_date){$run_date="";} #OPTIONAL ARGUMENT
	if (!defined $completion_time) {$completion_time="";} #OPTIONAL ARGUMENT
        if (!defined $query_name_in_msa) {$query_name_in_msa="";} #OPTIONAL ARGUMENT
	$chain=uc($chain);
        if ($chain eq "NONE"){$chain="";}
	my $pdb_dir = $results_dir.$IN_pdb_id_capital."/";
	my ($header_line, $title_line, $identical_chains_line, $ref_outputs,$cmd);
	
	# read info from the pdb file
	my @ans = @$ref_header_title;
	if (exists $ans[1]){
	    $header_line = $ans[1];
	}
	if (exists $ans[2]){
	    $title_line = $ans[2];
	}
	elsif (exists $ans[3]){
	    $title_line = $ans[3];
	}
	if ($title_line eq ""){
		$title_line = "! \"No title or compound description was found in the PDB file\";";
	}
	else{
		$title_line = design_string_with_spaces_for_pipe($title_line);}
	
	if ($current_dir !~ /(\/)$/){$current_dir=$current_dir."/";}	
	#$current_dir = $pdb_dir.$chain."/";
	
	# design the identical chains line
	$identical_chains_line = "! consurf_identical_chains = \"$identical_chains\";";
	
	# in case there is a source dir - we determine the var $query_name_in_msa
	if(-e $current_dir."_sourcedir"){
		unless (open SOURCEDIR, $current_dir."_sourcedir"){return ("create_consurf_pipe : cannot open ".$current_dir."_sourcedir $!", 'PANIC');}
		my $line = <SOURCEDIR>;
		if ( $line =~ /(\d[\d\w]{3})\/(\w)/){
			$query_name_in_msa = $1.$2;
			chomp($line);
		}
		close SOURCEDIR;	
	}		
	#### if we wish to read from source
	## look for the file partOfPipe
	#if (-e $current_dir."partOfPipe"){
	#	$dir_to_copy_from = $current_dir;
	#}
	#elsif(-e $current_dir."_sourcedir"){
	#	unless (open SOURCEDIR, $current_dir."_sourcedir"){return ("create_consurf_pipe : cannot open ".$current_dir."_sourcedir $!", 'PANIC');}
	#	my $line = <SOURCEDIR>;
	#	if ( $line =~ /(\d[\d\w]{3})\/(\w)/){
	#		$query_name_in_msa = $1.$2;
	#		chomp($line);
	#		$dir_to_copy_from = $results_dir.$line."/";
	#	}
	#	close SOURCEDIR;	
	#}
	#### if we wish to read from source
	if ($query_name_in_msa eq ""){
	    $query_name_in_msa = $IN_pdb_id_capital.$chain; # should be the reference sequence name in the MSA
	}
	##########################
	# write to the pipe file #
	##########################
#	unless (open PIPE_PART, $current_dir."partOfPipe"){return("create_consurf_pipe : cannot open $current_dir"."partOfPipe for reading $!",'PANIC');}
	unless (open PIPE_PART, $partOfPipe){return("create_consurf_pipe : cannot open $partOfPipe for reading $!",'PANIC');}
	unless (open PIPE, ">".$final_pipe_file){return("create_consurf_pipe : cannot open ".$final_pipe_file." for writing $!", 'PANIC');}
	if ($header_line ne ""){
		print PIPE $header_line."\n";
	}
	else{
		print PIPE "HEADER                                 [THIS LINE ADDED FOR JMOL COMPATIBILITY]\n";
	}
	print PIPE <<EndOfPipe;	
!! ====== IDENTIFICATION SECTION ======
!js.init
! consurf_version = "3.0";
! consurf_run_number = "$run_number";
! consurf_run_date = "$run_date";
! consurf_run_submission_time = "$submission_time";
! consurf_run_completion_time = "$completion_time";
!
! consurf_pdb_id = "$IN_pdb_id_capital";
! consurf_chain = "$chain";
$identical_chains_line
! consurf_msa_filename = "$msa_filename";
! consurf_msa_query_seq_name = "$query_name_in_msa";
! consurf_tree_filename = "$tree_filename";
!
EndOfPipe
	my $titleFlag = 0;
	while(<PIPE_PART>){
		if ($_  =~ /^~~~+/){
			if( $titleFlag==0){
				print PIPE "! pipe_title = \"<i>ConSurf View:</i> $IN_pdb_id_capital chain $chain.\"\n";
				print PIPE "!! pipe_subtitle is from TITLE else COMPND\n";
				print PIPE "!!\n";
				print PIPE "! pipe_subtitle =\n";
				print PIPE $title_line."\n";
				$titleFlag=1;
			}
			else{
			  if ($chain ne "")
			    {
			      print PIPE "! select selected and :$chain\n";
			    }
			  else
			    {
			       print PIPE "! select selected and protein\n";
			    }
			}
		}
		else{
			print PIPE $_;}
	}        
	close PIPE_PART;
	close PIPE;
	unlink $partOfPipe;
	return ("OK");
}

#************************************************************************************************
# create the file to be shown using FGiJ. read the pdb file and concat header pipe to it.
sub add_pdb_data_to_pipe{
	my ($pdb_file,$pipe_file) = @_;
	
	unless (open (PIPE,">>$pipe_file")) {return ("add_pdb_data_to_pipe: cannot open $pipe_file for writing $!");}
	unless (open (PDB_FILE, $pdb_file)) {return ("add_pdb_data_to_pipe: cannot open the $pdb_file for writing $!");}
	
	while (my $line=<PDB_FILE>)
	{
		if ($line!~/^HEADER/){
                        print PIPE $line;
                    }
        }
	close (PIPE);
	close (PDB_FILE);
	return ("OK");
}
#************************************************************************************************
sub create_consurf_pipe{
	my ($results_dir, $IN_pdb_id_capital, $chain, $ref_header_title, $final_pipe_file, $identical_chains) = @_;
	
	my $pdb_dir = $results_dir.$IN_pdb_id_capital."/";
	my ($header_line, $title_line, $identical_chains_line, $current_dir, $ref_outputs,$cmd);
	my $query_name_in_msa = "";
	# read info from the pdb file
	my @ans = @$ref_header_title;
	if (exists $ans[1]){
	    $header_line = $ans[1];
	}
	if (exists $ans[2]){
	    $title_line = $ans[2];
	}
	elsif (exists $ans[3]){
	    $title_line = $ans[3];
	}
	if ($title_line eq ""){
		$title_line = "! \"No title or compound description was found in the PDB file\";";
	}
	else{
		$title_line = design_string_with_spaces_for_pipe($title_line);}
	
	$current_dir = $pdb_dir.$chain."/";
	
	# design the identical chains line
	$identical_chains_line = "! consurf_identical_chains = \"$identical_chains\";";
	
	# in case there is a source dir - we determine the var $query_name_in_msa
	if(-e $current_dir."_sourcedir"){
		unless (open SOURCEDIR, $current_dir."_sourcedir"){return ("create_consurf_pipe : cannot open ".$current_dir."_sourcedir $!", 'PANIC');}
		my $line = <SOURCEDIR>;
		if ( $line =~ /(\d[\d\w]{3})\/(\w)/){
			$query_name_in_msa = $1.$2;
			chomp($line);
		}
		close SOURCEDIR;	
	}		
	#### if we wish to read from source
	## look for the file partOfPipe
	#if (-e $current_dir."partOfPipe"){
	#	$dir_to_copy_from = $current_dir;
	#}
	#elsif(-e $current_dir."_sourcedir"){
	#	unless (open SOURCEDIR, $current_dir."_sourcedir"){return ("create_consurf_pipe : cannot open ".$current_dir."_sourcedir $!", 'PANIC');}
	#	my $line = <SOURCEDIR>;
	#	if ( $line =~ /(\d[\d\w]{3})\/(\w)/){
	#		$query_name_in_msa = $1.$2;
	#		chomp($line);
	#		$dir_to_copy_from = $results_dir.$line."/";
	#	}
	#	close SOURCEDIR;	
	#}
	#### if we wish to read from source
	if ($query_name_in_msa eq ""){
	    $query_name_in_msa = $IN_pdb_id_capital.$chain; # should be the reference sequence name in the MSA
	}
	##########################
	# write to the pipe file #
	##########################
	unless (open PIPE_PART, $current_dir."partOfPipe"){return("create_consurf_pipe : cannot open $current_dir"."partOfPipe for reading $!",'PANIC');}
	unless (open PIPE, ">".$current_dir.$final_pipe_file){return("create_consurf_pipe : cannot open ".$current_dir.$final_pipe_file." for writing $!", 'PANIC');}
	if ($header_line ne ""){
		print PIPE $header_line."\n";
	}
	else{
		print PIPE "HEADER                                 [THIS LINE ADDED FOR JMOL COMPATIBILITY]\n";
	}
	print PIPE <<EndOfPipe;	
!! ====== IDENTIFICATION SECTION ======
!js.init
! consurf_version = "3.0";
! consurf_run_number = "ConSurf_DB/$IN_pdb_id_capital/$chain";
! consurf_run_date = "";
! consurf_run_submission_time = "";
! consurf_run_completion_time = "";
!
! consurf_pdb_id = "$IN_pdb_id_capital";
! consurf_chain = "$chain";
$identical_chains_line
! consurf_msa_filename = "msa.aln";
! consurf_msa_query_seq_name = "$query_name_in_msa";
! consurf_tree_filename = "";
!
EndOfPipe
	my $titleFlag = 0;
	while(<PIPE_PART>){
		if ($_  =~ /^~~~+/){
			if( $titleFlag==0){
				print PIPE "! pipe_title = \"<i>ConSurf View:</i> $IN_pdb_id_capital chain $chain.\"\n";
				print PIPE "!! pipe_subtitle is from TITLE else COMPND\n";
				print PIPE "!!\n";
				print PIPE "! pipe_subtitle =\n";
				print PIPE $title_line."\n";
				$titleFlag=1;
			}
			else{
				print PIPE "! select selected and :$chain\n";
			}
		}
		else{
			print PIPE $_;}
	}        
	close PIPE_PART;
	close PIPE;
	unlink $current_dir."partOfPipe";
	return ("OK");
}

#************************************************************************************************
# Go over the PDB file and find all the chains that are identical to the given chain
sub find_identical_chains_on_PDB_File{
	my $PDBfile=shift;
	my $chain=shift;
	unless (open PDBFILE , "$PDBfile") {return "find_identical_chains_on_PDB_File : cannot open $PDBfile for reading $!\n";}
        my $IdenticalChainsLine = $chain;
        my %chainsHash=();
	while (<PDBFILE>)
	{
	      chomp($_);
    	      if (/^SEQRES\s+\d+\s+([A-Z0-9])\s+\d+\s+([A-Z\s]+)\s+/i)
		{
              	my $key = $1;
              	my  $line = $2;
              	$line =~ s/\s+//g;
              	$chainsHash{$key} .= $line;
		}
	}
	close PDBFILE;
	foreach my $key (keys  %chainsHash)
	{
    		if ($chainsHash{$key} eq $chainsHash{$chain} && ($key ne $chain))
		{
        	$IdenticalChainsLine .= $key;
    		}
	}
	return ("OK",$IdenticalChainsLine);
}

#************************************************************************************************
# create a file with all the ATOM records from a pdb.
# a list of residues, serial number and their position according to the pdb
# all the outpus are written to a hash (which is an input to the routine), so in order to extract output information:
# $ref_to_return_hash->{ERROR} : error
# $ref_to_return_hash->{NMR} : will be "yes" if it is NMR
# $ref_to_return_hash->{INFO} : information regarding nmr model
# $ref_to_return_hash->{WARNING} : in case a non-standard residue was found, will give its description
# $ref_to_return_hash->{ATOM_AA_SEQ} : will contain the ATOM sequence in 1 letter amino-acids
sub create_atom_position_file{
	my ($pdb_filename ,$atom_position_filename, $chain_id, $ref_to_return_hash) = @_;
	my $nmr_model = "no";
	my $info_ret = "";
	my $warn_ret = "";
	my ($residue,$pos,$last_pos);
	my $first=1;
	my $fas="";
	my $pdb_line =0;
	if (-z $pdb_filename or !-e $pdb_filename){
		$ref_to_return_hash->{ERROR} = "rasmol_gradesPE_and_pipe::create_atom_position_file : the file $pdb_filename does not exist ot is empty";
		return;
	}
	unless (open PDB, $pdb_filename){
		$ref_to_return_hash->{ERROR} = "rasmol_gradesPE_and_pipe::create_atom_position_file : could not open the file $pdb_filename for reading $!";
		return;
	}
	unless (open CORR, ">$atom_position_filename"){
		$ref_to_return_hash->{ERROR} = "rasmol_gradesPE_and_pipe::create_atom_position_file : could not open the file $atom_position_filename for writing $!";
		return;
	}	
    #going over the PDB file and extracting the chain sequence.
    while(<PDB>){
		$pdb_line++;
        # in case of NMR, we only take the first model's info
        if (/^MODEL\s+1\s*$/){
            $nmr_model = "yes";
            $ref_to_return_hash->{INFO} = "found NMR model 1.";
        }
        elsif (/^MODEL\s+(\d+)\s*$/ and $nmr_model eq "yes"){
            $ref_to_return_hash->{INFO} .= " stop reading PDB at NMR model $1";
            last;
        }
        #reaches an ATOM line, extracts the residue and prints out its position.
        #elsif($_=~ m/^ATOM.................$chain_id/){
		elsif($_=~ m/^ATOM/ and (substr ($_,21,1) eq $chain_id)){
			$residue=substr ($_,17,3);
            $pos=substr ($_,22,5);;
            if(exists($tr_aa{$residue})){
				if($first ==1){
					$fas=$fas.$tr_aa{$residue};   
                    $last_pos=$pos;
                    $fas_pos=1;
					print CORR "$residue\t$fas_pos\t$pos\n";
				}
				$first=0;
                if($pos ne $last_pos){
			if ($pos=$last_Pos+1)
			{
                    $fas=$fas.$tr_aa{$residue};   
                    $fas_pos++;
					print CORR "$residue\t$fas_pos\t$pos\n";
			}
			else # HAIM ADD FOR DISORDER REGION
			{
			$fas=$fas."X";   
                     $fas_pos++;
			}
                }
                $last_pos=$pos;
            }
			else{
				$ref_to_return_hash->{WARNING}.= "line $pdb_line : residue $residue is not legal\n";
			}
		}
	}
	close PDB;
	close CORR;
	if (!-e $atom_position_filename or -z $atom_position_filename){
		$ref_to_return_hash->{ERROR} = "rasmol_gradesPE_and_pipe::create_atom_position_file : the file $atom_position_filename was not created or of size 0";
		return;
	}
	$ref_to_return_hash->{ATOM_AA_SEQ} = $fas;
	$ref_to_return_hash->{NMR} = "yes" if ($nmr_model eq "yes");
}
#************************************************************************************************
#############################################################
# Description:  replace the tempFactor column in the PDB file
#               with another column. (for example: grades)
# Arguments: 1. PDB file
#            2. chain
#            3. reference to hash
#               key - the residue sequence number
#               value - the field to insert instead of the tempFactor
#############################################################
sub replace_tempFactor {

    my $PdbFile = shift;
    my $chain = shift;
    my $HashRef = shift;
    my $Out=shift;
    
    my $PDBchain;
    my $ResNum;
    my $iCode;
    my $ret = "nothing"; 
    
    # read the PDB file to an array
    open READPDB, "<$PdbFile";
    my @PDB = <READPDB>;
    close READPDB;

    # write the PDB file and replace the tempFactor column
    # with the new one.
    open WRITEPDB, ">$Out";
    
    foreach my $line (@PDB) {
	if ($line =~ /^ATOM/){	 
	    $PDBchain = substr ($line, 21, 1);
            $ResNum = substr ($line, 22, 4);
	    $ResNum =~ s/\s*//g;
	    $iCode = substr ($line, 26, 1);
            if ($iCode =~ /^\s*$/){
                $iCode = '';
            }
            if ($PDBchain eq $chain){
                my $ResNumiCode=$ResNum . $iCode;
                $$HashRef{$ResNumiCode} =~ s/\s*$//;
		while (length($$HashRef{$ResNumiCode})<6)
		{
			$$HashRef{$ResNumiCode}=" ".$$HashRef{$ResNumiCode};
		}
	    	if (length($line)>60 && length($line)>=66)# to make sure that the temprature factor is written in the line
                    {substr ($line, 60, 6) = $$HashRef{$ResNumiCode};}
                else{
                    $ret = "no_TF";
                }
                print WRITEPDB $line;
            }
            else { # Not The Requested Chain
	        if (length($line)>60 && length($line)>=66)# to make sure that the temprature factor is written in the line
                    {substr ($line, 60, 6) = "      ";}
                else{
                    $ret = "no_TF";
                }
                print WRITEPDB $line;
            }
	}
	else {
	    print WRITEPDB $line;
	}
    }
    close WRITEPDB;
    return $ret;
}
#************************************************************************************************
# Creates The ATOM section with ConSurf grades instead of the TempFactor column
sub ReplaceTempFactConSurfScore
{
    my $query_chain = shift;
    my $in = shift;		#path to <file>.ent
    my $PE = shift;		#path to <file>.gradesPE
    my $out1 = shift;	#path to the output file
    my $out2 = shift;	#path to the insuficient data output file
    my ($out, $line, $residue, $insufficient, $tempFactor_isd);
    my @grades;
    my %gradesPE_info_with0 = ();
    my %gradesPE_info = ();
	
    
    $insufficient = &read_ConSurf_gradesPE($PE, \%gradesPE_info, \%gradesPE_info_with0);
    ###################################################
    #$in should be the name of the nest output pdb file
    ##################################################### 
  
    unless  (open (PDB, "<$in")){
        return "could not open file '$in' $!\n";
    }
    else{
        #open (OUT, ">$out") or die "could not open file"; 
        
        my $i = -1;
        $_ = 'A B C D E F G H I J K L M N O P Q R S T U V W X Y Z'; my @j = split; # to define letter that will be the chain, to split in case the input has more than one pdb. not really necessary now.
        my $flag=0;
        
        unless (open (OUTP, ">$out1")){return "could not open the file '$out1' $!\n";}
        # If * was found in gradesPE file, we create another output, where the * grades are replaces with 0
        if ($insufficient eq "yes"){
            unless (open OUT_ISD, ">".$out2) {return "could not open the file '$out2' $!\n";}
        }
        while ($line = <PDB>){
            $_=$line;
            if (m/^\s*$/){next;}
            #if (m/^MODEL/) {$i++; print OUT $_; }
            
            elsif (m/^ATOM/){
                my $atom = substr($line, 0, 6); 
                my $atom_num = substr($line, 6, 5);
                my $atom_name = substr($line, 11, 5);
                my $altLoc = substr($line, 16, 1);
                my $res_name = substr($line, 17, 4);
                my $chain = substr($line, 21, 1);
                my $res_num = substr($line, 22, 4); 
                my $iCode = substr($line, 26, 4);
                my $x = substr($line, 30, 8);
                my $y = substr($line, 38, 8);
                my $z = substr($line, 46, 8);
                my $occupancy = substr($line, 54, 6);
                my $tempFactor = substr($line, 60, 12);
                #my $segID = substr($line, 72, 4);
                #my $element = substr($line, 76, 2);
                
                $tempFactor_isd = $tempFactor;
                if ($chain eq $query_chain or ($chain =~ m/\s+/ and $query_chain =~ m/none/i)){
                    $residue = $res_num;
                    $residue =~ s/\s//g;
                    if (exists $gradesPE_info{$residue}){
                        # the TF is updated with the number from gradesPE
                        $tempFactor = "     $gradesPE_info{$residue}      "; 
                    }
                    if ($insufficient eq "yes" && exists $gradesPE_info_with0{$residue}){
                        $tempFactor_isd = "     $gradesPE_info_with0{$residue}      "; 
                    }                        
                }
		else
		{
			$tempFactor = "            ";
			$tempFactor_isd = "            ";
		}
                print OUTP "$atom$atom_num$atom_name$altLoc$res_name$chain$res_num$iCode$x$y$z$occupancy$tempFactor\n";
                print OUT_ISD "$atom$atom_num$atom_name$altLoc$res_name$chain$res_num$iCode$x$y$z$occupancy$tempFactor_isd\n" if ($insufficient eq "yes" );
            }                
        }#while               
    }
    close OUTP;
    close OUT_ISD;
    #!system ("rm -f $out");
    return "OK";
}
#************************************************************************************************
sub read_ConSurf_gradesPE{
# the routine matches each position in the gradesPE file its ConSurf grade. In case there was a grade mark with *, we put it in a seperate hash with the grade 0.
# the routine returns "yes" if a * was found and "no" otherwise
    my $gradesPE_file = shift;
    my $gradesPE_hash_ref = shift;
    my $gradesPE_0_hash_ref = shift;
	my $insufficient = "no";
    
    open GRADES, $gradesPE_file;
    while (<GRADES>){
        if (/^\s*\d+\s+\w/ ){
		    my @grades=split;            
		    $grades[2] =~ s/[a-z\:]//gi;
            if ($grades[4] =~/\d\*?/){
				# if it is insufficient color - we change its grade to 0, which will be read as light yellow
				if ($grades[4] =~/(\d)\*/){
					$gradesPE_hash_ref->{$grades[2]} = $grades[4];
                    $gradesPE_0_hash_ref->{$grades[2]} = 0;
                    $insufficient = "yes";
				}
				else{
					$gradesPE_hash_ref->{$grades[2]} = $grades[4];
                    $gradesPE_0_hash_ref->{$grades[2]} = $grades[4];
				}
            }
        }
    }
    close GRADES; 
    return $insufficient;
}
#************************************************************************************************
sub read_Rate4Site_gradesPE{
# the routine matches each position in the gradesPE file its Rate4Site grade. 
    my $gradesPE_file = shift;
    my $gradesPE_hash_ref = shift;
    
    open GRADES, $gradesPE_file;
    while (<GRADES>){
        if (/^\s*\d+\s+\w/ ){
		    my @grades=split;            
		    $grades[2] =~ s/[a-z\:]//gi;
	            $gradesPE_hash_ref->{$grades[2]} = $grades[3];
	}
    }
    close GRADES; 
}
#************************************************************************************************
# Chimera Outputs Functions
#************************************************************************************************
sub color_with_chimera{
# For chimera scripts, in order to color the positions in the MSA only according to query sequence and ignore positions where there are gaps, the MSA and the R4S files should be read simultaneously. 
# the rorutine will only create 'insufficient data' scripts, if insufficient data was calculated from r4s grades.
# input: MSA - CLUSTALW format only!!
	#input:
	my $dir_path = shift; # where input files are located
	my $query_name = shift;
	my $msa_file = shift;
	my $rate4site_filename = shift;
	# files that will be created:
	my $chimera_scf = shift;
	my $headers =  shift;
	my $chimera_scf_isd = shift;
	my $headers_isd = shift;	
	my $out_path=shift; # Optional, if not given equals to the $dir_path
	if ($dir_path !~/(\/)$/) {$dir_path=$dir_path."/";}	
	if (!defined $out_path or $out_path eq ""){$out_path=$dir_path;}
	if ($out_path !~/(\/)$/) {$out_path=$out_path."/";}
	my @r4s_position_grade = ();
	my %aln_position_grade = ();
	my %aln_position_grade_isd = ();
	my $isd = "no";	
	
	my @ans = assign_colors_according_to_r4s_layers($dir_path.$rate4site_filename, \@r4s_position_grade);
	if ($ans[0] ne "OK") {die $ans[0];}
	
	my $ref_r4s_position_grade = \@r4s_position_grade;
	
	
	#open MSA file, read only the query sequence lines. assign each position with increasing number. if it is not a '-' : read the grade from the array @r4s_position_grade and put it in the hash %aln_position_grade.
	#if there is insufficient data: save it to a saparate hash
	
	my $position_in_msa = 0;
	my $position_in_r4s = 0;
	my $FOUND_QUERY=0; #Found the Query exact name
	unless (open MSA, $dir_path.$msa_file) {return ("cannot open $msa_file $!");}
	while(<MSA>){
		chomp;
		if (/^(.+)\s+([A-Za-z-]+)$/){
			my $seq_name = $1;
			my $seq = $2;
			$seq_name =~ s/\s+$//;
			#next if ($query_name !~ $seq_name);
			next if ($query_name ne $seq_name); #FIND QUERY EXACT NAME
			# we are in the query name
			$FOUND_QUERY=1;
			$seq =~ s/\s+$//; # remove spaces
			my $seq_length = length($seq);
			#print "LENGTH:$seq_length\tSEQ:*$seq*\n";
			my @sequence_block = split //, $seq;
			foreach(@sequence_block){
				#print "POS IN R4S:$position_in_r4s\t*$_*\n";
				unless(/-/){
					#print "ISD:".$ref_r4s_position_grade->[$position_in_r4s]{ISD}." ";
					if ($ref_r4s_position_grade->[$position_in_r4s]{ISD} == 1){						
						$aln_position_grade_isd{$position_in_msa} =0;
						$isd = "yes" if ($isd eq "no");						
					}
					else{
						$aln_position_grade_isd{$position_in_msa} = $ref_r4s_position_grade->[$position_in_r4s]{COLOR};
					}
					$aln_position_grade{$position_in_msa} = $ref_r4s_position_grade->[$position_in_r4s]{COLOR};
					$position_in_r4s++;
				}
				$position_in_msa++;
			}
		}
	}
	close MSA;
	if (!$FOUND_QUERY) # DONT FIND QUERY EXACT NAME THUS LOOP AGAIN AND LOOK FOR LINE CONTAIN THE QUERY NAME (FOR TRUNCATED NAMES)
	{
		unless (open MSA, $dir_path.$msa_file) {return ("cannot open $msa_file $!");}
		while(<MSA>){
			chomp;
			if (/^(.+)\s+([A-Za-z-]+)$/){
				my $seq_name = $1;
				my $seq = $2;
				$seq_name =~ s/\s+$//;
				next if ($query_name !~ $seq_name);
				# we are in the query name
				$seq =~ s/\s+$//; # remove spaces
				my $seq_length = length($seq);
				#print "LENGTH:$seq_length\tSEQ:*$seq*\n";
				my @sequence_block = split //, $seq;
				foreach(@sequence_block){
					#print "POS IN R4S:$position_in_r4s\t*$_*\n";
					unless(/-/){
						#print "ISD:".$ref_r4s_position_grade->[$position_in_r4s]{ISD}." ";
						if ($ref_r4s_position_grade->[$position_in_r4s]{ISD} == 1){						
							$aln_position_grade_isd{$position_in_msa} =0;
							$isd = "yes" if ($isd eq "no");						
						}
						else{
							$aln_position_grade_isd{$position_in_msa} = $ref_r4s_position_grade->[$position_in_r4s]{COLOR};
						}
						$aln_position_grade{$position_in_msa} = $ref_r4s_position_grade->[$position_in_r4s]{COLOR};
						$position_in_r4s++;
					}
					$position_in_msa++;
				}
			}
		}
		close MSA;
	}
	
	&print_chimera_scf_and_histogram($out_path.$chimera_scf, $out_path.$headers, \%aln_position_grade, $position_in_msa);
	&print_chimera_scf_and_histogram($out_path.$chimera_scf_isd,$out_path.$headers_isd, \%aln_position_grade_isd, $position_in_msa) if ($isd eq "yes");
	return ("OK",$isd);	
}
#************************************************************************************************
sub print_chimera_scf_and_histogram{
# FORMERLY: sub print_chimera_scripts{
    my $chimera_scf = shift;
    my $chimera_headers = shift;
    my $ref_aln_position_grade = shift;
    my $position_in_msa = shift;
    
    open (SCF, ">".$chimera_scf) || die "Cna't Open $chimera_scf $!";
    open HEADER, ">".$chimera_headers;
    print HEADER "name: Conservation scores\nstyle: character\n";
    my $j=0;
    for(my $i=0; $i<=$position_in_msa; $i++){
        $j=$i+1;
        if ((exists $ref_aln_position_grade->{$i}) and ($ref_aln_position_grade->{$i} ne "")){        
            print SCF "$i $i 0 0 $consurf_colors_chimera{$ref_aln_position_grade->{$i}}\n";
            print HEADER "\t$j\t$ref_aln_position_grade->{$i}\tblack\n";
        }
    }
    close SCF;
    
    print HEADER "#  ConSurf data shown as histogram\nname: ConSurf histogram\nstyle: numeric\n";
    $j=0;
    my $histo;
    for(my $i=0; $i<=$position_in_msa; $i++){
        $j=$i+1;
        if ((exists $ref_aln_position_grade->{$i}) and ($ref_aln_position_grade->{$i} ne "")){
            if ($ref_aln_position_grade->{$i} == 9){$histo = "1.0";}
            else{$histo = ".".($ref_aln_position_grade->{$i}+1);}
            
            print HEADER "\t$j\t$histo\tblack\n";
        }
    }
    close HEADER;
}
#************************************************************************************************
sub create_chimera_image_script{
	my $script_file = shift;  # full path to open the chimerax script from
	my $pdb_file_name = shift;  # should end with '.pdb'
	my $pdb_path = shift;  # the www path to the working directory
	if ($pdb_path !~/(\/)$/) {$pdb_path.="/";}
	$pdb_path.=$pdb_file_name;
	
	unless (open CHIMERAX, ">".$script_file) {return "Failed to open $script_file $!";}
	print CHIMERAX <<EndOfScript;
<?xml version="1.0"?>
  <ChimeraPuppet type="std_webdata">
<web_files>
<file name="$pdb_file_name" format="text" loc="$pdb_path"/>
</web_files>
<commands>
  <mid_cmd>colordef CONS10 1.00 1.00 0.59</mid_cmd>
  <mid_cmd>colordef CONS9 0.63 0.15 0.38</mid_cmd>
  <mid_cmd>colordef CONS8 0.94 0.49 0.67</mid_cmd>
  <mid_cmd>colordef CONS7 0.98 0.79 0.87</mid_cmd>
  <mid_cmd>colordef CONS6 0.99 0.93 0.96</mid_cmd>
  <mid_cmd>colordef CONS5 1.00 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS4 0.92 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS3 0.84 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS2 0.55 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS1 0.06 0.78 0.82</mid_cmd>
  <mid_cmd>color CONS10 @/bfactor=0</mid_cmd>
  <mid_cmd>color CONS9 @/bfactor=9 </mid_cmd>
  <mid_cmd>color CONS8 @/bfactor=8</mid_cmd>
  <mid_cmd>color CONS7 @/bfactor=7</mid_cmd>
  <mid_cmd>color CONS6 @/bfactor=6</mid_cmd>
  <mid_cmd>color CONS5 @/bfactor=5</mid_cmd>
  <mid_cmd>color CONS4 @/bfactor=4</mid_cmd>
  <mid_cmd>color CONS3 @/bfactor=3</mid_cmd>
  <mid_cmd>color CONS2 @/bfactor=2</mid_cmd>
  <mid_cmd>color CONS1 @/bfactor=1</mid_cmd>
  <mid_cmd>preset apply pub 3;repr cpk;show;focus;color red ligand</mid_cmd>
</commands>
</ChimeraPuppet>
	
EndOfScript
	close CHIMERAX;
	chmod 0755, $script_file;
	return "OK";
}

#************************************************************************************************
sub create_chimera_script_align_tree{
	my ($www_dir, $files_dir, $chimerax_file, $msa_file, $tree_file, $scf_file, $hdr_file) = @_;
	
	if ($www_dir !~/(\/)$/){$www_dir=$www_dir."/";}
	if ($files_dir !~/(\/)$/){$files_dir=$files_dir."/";}
	
	open CHIMERAX, ">".$files_dir.$chimerax_file;
	print CHIMERAX <<EndOfFile;
<?xml version="1.0"?>
  <ChimeraPuppet type="std_webdata">
<web_files>
EndOfFile
	if ($msa_file ne ""){
		print CHIMERAX <<EndOfFile;
<file  name=\"$msa_file\" format=\"text\" loc=\"$www_dir$msa_file\"/>
</web_files>
<commands>

<!-- the following 3 lines locate the Multalign Viewer instance
	that was created by opening the alignment file, and stores a reference
	to the instance as the variable 'mav' -->
  <py_cmd>from MultAlignViewer.MAViewer import MAViewer</py_cmd>
  <py_cmd>from chimera.extension import manager</py_cmd>
  <py_cmd>mav = [inst for inst in manager.instances if isinstance(inst, MAViewer)][-1]</py_cmd>
<!-- read in/show Consurf tree -->

<!-- hide initial headers, show phylogeny tree, load the coloring file,
	and make residue letters black.
	This uses two possible sets of calls:  one for the 1.2540 release
	and one for later versions that uses a better API.
	The 'if' condition guarantees that the code will work no
	matter what version the user has -->
<py_cmd>
if hasattr(mav, 'loadScfFile'):
	mav.hideHeaders(mav.headers(shownOnly=True))
	mav.usePhylogenyFile("$www_dir$tree_file", askReorder=False)
	mav.loadScfFile("$www_dir$scf_file")
	mav.useColoringFile(None)
else:
	mav.hideHeaders(mav.seqCanvas.headerDisplayOrder())
	mav.usePhylogenyFile("$www_dir$tree_file")
	mav.regionBrowser.loadScfFile("$www_dir$scf_file")
	from MultAlignViewer.prefs import RC_BLACK
	mav.seqCanvas.setColorFunc(RC_BLACK)
</py_cmd>
<!-- read in/show Consurf headers -->
	<py_cmd>mav.readHeaderFile(\"$www_dir$hdr_file\")</py_cmd>
EndOfFile
	}
	# if there is no MSA, print the relevant colors in order to color the molecule
	else{
		print CHIMERAX <<EndOfColors;
</web_files>
<commands>
  <mid_cmd>colordef CONS10 1.00 1.00 0.59</mid_cmd>
  <mid_cmd>colordef CONS9 0.63 0.15 0.38</mid_cmd>
  <mid_cmd>colordef CONS8 0.94 0.49 0.67</mid_cmd>
  <mid_cmd>colordef CONS7 0.98 0.79 0.87</mid_cmd>
  <mid_cmd>colordef CONS6 0.99 0.93 0.96</mid_cmd>
  <mid_cmd>colordef CONS5 1.00 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS4 0.92 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS3 0.84 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS2 0.55 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS1 0.06 0.78 0.82</mid_cmd>
  <mid_cmd>color CONS10 @/bfactor=0</mid_cmd>
  <mid_cmd>color CONS9 @/bfactor=9 </mid_cmd>
  <mid_cmd>color CONS8 @/bfactor=8</mid_cmd>
  <mid_cmd>color CONS7 @/bfactor=7</mid_cmd>
  <mid_cmd>color CONS6 @/bfactor=6</mid_cmd>
  <mid_cmd>color CONS5 @/bfactor=5</mid_cmd>
  <mid_cmd>color CONS4 @/bfactor=4</mid_cmd>
  <mid_cmd>color CONS3 @/bfactor=3</mid_cmd>
  <mid_cmd>color CONS2 @/bfactor=2</mid_cmd>
  <mid_cmd>color CONS1 @/bfactor=1</mid_cmd>
  <mid_cmd>preset apply pub 3;repr cpk;show;focus;color red ligand</mid_cmd>
EndOfColors
	}
	print CHIMERAX "\n  </commands>\n</ChimeraPuppet>";
	close CHIMERAX;
	chmod 0755, $chimerax_file;
}

sub create_chimera_script{
	my ($pdb_file_name, $www_dir, $files_dir, $chimerax_file, $msa_file, $tree_file, $scf_file, $hdr_file) = @_;
	
	if ($www_dir !~/(\/)$/){$www_dir=$www_dir."/";}
	if ($files_dir !~/(\/)$/){$files_dir=$files_dir."/";}
	
	open CHIMERAX, ">".$files_dir.$chimerax_file;
	print CHIMERAX <<EndOfFile;
<?xml version="1.0"?>
  <ChimeraPuppet type="std_webdata">
<web_files>
<file  name="$pdb_file_name" format="text" loc="$www_dir$pdb_file_name"/>
EndOfFile
	if ($msa_file ne ""){
		print CHIMERAX <<EndOfFile;
<file  name=\"$msa_file\" format=\"text\" loc=\"$www_dir$msa_file\"/>
</web_files>
<commands>
  <mid_cmd>preset apply pub 3;repr cpk;show;focus;color red ligand</mid_cmd>

<!-- the following 3 lines locate the Multalign Viewer instance
	that was created by opening the alignment file, and stores a reference
	to the instance as the variable 'mav' -->
  <py_cmd>from MultAlignViewer.MAViewer import MAViewer</py_cmd>
  <py_cmd>from chimera.extension import manager</py_cmd>
  <py_cmd>mav = [inst for inst in manager.instances if isinstance(inst, MAViewer)][-1]</py_cmd>
<!-- read in/show Consurf tree -->

<!-- hide initial headers, show phylogeny tree, load the coloring file,
	and make residue letters black.
	This uses two possible sets of calls:  one for the 1.2540 release
	and one for later versions that uses a better API.
	The 'if' condition guarantees that the code will work no
	matter what version the user has -->
<py_cmd>
if hasattr(mav, 'loadScfFile'):
	mav.hideHeaders(mav.headers(shownOnly=True))
	mav.usePhylogenyFile("$www_dir$tree_file", askReorder=False)
	mav.loadScfFile("$www_dir$scf_file")
	mav.useColoringFile(None)
else:
	mav.hideHeaders(mav.seqCanvas.headerDisplayOrder())
	mav.usePhylogenyFile("$www_dir$tree_file")
	mav.regionBrowser.loadScfFile("$www_dir$scf_file")
	from MultAlignViewer.prefs import RC_BLACK
	mav.seqCanvas.setColorFunc(RC_BLACK)
</py_cmd>
<!-- read in/show Consurf headers -->
	<py_cmd>mav.readHeaderFile(\"$www_dir$hdr_file\")</py_cmd>
<!-- show chains other than the one in the alignment as gray ribbon -->
<py_cmd>
import chimera
m = chimera.openModels.list(modelTypes=[chimera.Molecule])[0]
for seq in m.sequences():
	if hasattr(seq.residues[0], "mavConSurfHistogram"):
		continue
	from chimera.colorTable import getColorByName
	gray = getColorByName("gray")
	for r in seq.residues:
		for a in r.atoms:
			a.display = False
		r.ribbonDisplay = True
		r.ribbonDrawMode = chimera.Residue.Ribbon_Round
		r.ribbonColor = gray
</py_cmd>
EndOfFile
	}
	# if there is no MSA, print the relevant colors in order to color the molecule
	else{
		print CHIMERAX <<EndOfColors;
</web_files>
<commands>
  <mid_cmd>colordef CONS10 1.00 1.00 0.59</mid_cmd>
  <mid_cmd>colordef CONS9 0.63 0.15 0.38</mid_cmd>
  <mid_cmd>colordef CONS8 0.94 0.49 0.67</mid_cmd>
  <mid_cmd>colordef CONS7 0.98 0.79 0.87</mid_cmd>
  <mid_cmd>colordef CONS6 0.99 0.93 0.96</mid_cmd>
  <mid_cmd>colordef CONS5 1.00 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS4 0.92 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS3 0.84 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS2 0.55 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS1 0.06 0.78 0.82</mid_cmd>
  <mid_cmd>color CONS10 @/bfactor=0</mid_cmd>
  <mid_cmd>color CONS9 @/bfactor=9 </mid_cmd>
  <mid_cmd>color CONS8 @/bfactor=8</mid_cmd>
  <mid_cmd>color CONS7 @/bfactor=7</mid_cmd>
  <mid_cmd>color CONS6 @/bfactor=6</mid_cmd>
  <mid_cmd>color CONS5 @/bfactor=5</mid_cmd>
  <mid_cmd>color CONS4 @/bfactor=4</mid_cmd>
  <mid_cmd>color CONS3 @/bfactor=3</mid_cmd>
  <mid_cmd>color CONS2 @/bfactor=2</mid_cmd>
  <mid_cmd>color CONS1 @/bfactor=1</mid_cmd>
  <mid_cmd>preset apply pub 3;repr cpk;show;focus;color red ligand</mid_cmd>
EndOfColors
	}
	print CHIMERAX "\n  </commands>\n</ChimeraPuppet>";
	close CHIMERAX;
	chmod 0755, $chimerax_file;
}
#************************************************************************************************
sub create_pymol_page{
        my $PyMol_Instruction_Page=shift;
	my $PDB_ID=shift;
	my $Chain_ID=shift;
	
	my $Conf_Link=shift;
	my $PyMol_ConSurf_Script=shift;
	
	my $PDB_With_ConSurf_Scores=shift;
	my $PDB_With_ConSurf_Scores_Showing_ISD=shift;
	
	
	unless (open (PYMOL_PAGE,">$PyMol_Instruction_Page")) {return "cp_rasmol_gradesPE_and_pipe::create_pymol_page: Can't open PyMol Output '$PyMol_Instruction_Page' $!\n";}
	
	print PYMOL_PAGE <<EndOfHTML;
<html><title>ConSurf Image using PyMOL for $PDB_ID Chain: $Chain_ID</title><br>
<body>
<?
include (\"/var/www/html/ConSurf/php/templates/definitions.tpl\");
include (\"/var/www/html/ConSurf/php/templates/output_header_no_refresh.tpl\");
?>	
<center> <h3>Create a high resolution PyMol figure for $PDB_ID Chain: $Chain_ID</h3></center>
It is recommended to create a new directory where the attached files will be saved.<br>
Please note: For all PyMOL commands, you should type absolute paths for the files you are about to download, from their new location (on your computer).<br><br>
EndOfHTML
if ($PDB_With_ConSurf_Scores_Showing_ISD ne ""){
	print PYMOL_PAGE "<br> You may view your protein colored according to conservation scores with a unique color for positions of low confidence cutoff. If this is your preferred option, please choose option 1(a). Otherwise - choose option 1(b).<br>\n";
	print PYMOL_PAGE "Click <a href = \"$Conf_Link\">HERE</a> for more information regarding confidence cutoff.\n";
	print PYMOL_PAGE "<br><br><ol><li>Download the <b>PDB_FILE</b> updated with ConSurf\'s colors.<br>\n";
	print PYMOL_PAGE "(a) <a href= \"$PDB_With_ConSurf_Scores_Showing_ISD\">PDB_FILE</a> showing Insufficient Data<br>\n";
	print PYMOL_PAGE "(b) <a href= \"$PDB_With_ConSurf_Scores\">PDB_FILE</a> hiding Insufficient Data\n";
	print PYMOL_PAGE "</li>\n";
}
else{
	print PYMOL_PAGE "<br><br><ol><li>Download the <b>PDB_FILE</b> updated with ConSurf\'s colors.<br>\n";
	print PYMOL_PAGE "<a href= \"$PDB_With_ConSurf_Scores\">PDB_FILE</a> hiding Insufficient Data\n";
}
print PYMOL_PAGE <<EndOfHTML_2;	
<li>Download the file <a href= \"$PyMol_ConSurf_Script\" target =\"PyMOL_py\"><b>consurf_new.py</b></a>.</li>
<li>Start PyMOL.</li>
<li>A. Load the <font color=red> Consurf's modified pdb file (not the original pdb file that used as an input) </b> </font><br>In the PyMOL viewer window type:<br>
<font face=Courier New>PyMOL>\"load <b> PDB_FILE</b>\"
(no quotes) and hit return.<br></font>
</li>
B. Run the script to define ConSurf\'s color; Type:<br>
<font face=Courier New>PyMOL>\"run <b> consurf_new.py</b>\" 
      (no quotes) and hit return.<br> </font>
	 <br><font color=red>
	 Note: Please don't forget to write the path of the script file unless the script file is located at the pymol's root directory.
	 <br></font>
	 
</font>
</ol></body></html>

<?
    include(\"/var/www/html/ConSurf/php/templates/footer.tpl\");
?>
EndOfHTML_2
return ("OK");
}
#************************************************************************************************
sub create_chimera_page{
# Create the html describing how to create Chimera Hige resolution images	
	my $chimera_instructions_file = shift;
	my $conf_link= shift;
	my $chimera_consurf_commands = shift;
	
	
	my $chimerax_script = shift;
	my $PDB_atoms_ConSurf = shift;
	
	my $isd_chimerax_script =shift;
	my $isd_PDB_atoms_ConSurf = shift;
	
	my $url=shift;
	
	if (!defined $url) {$url="";}
	unless (open CHIMERA_HTML, ">".$chimera_instructions_file) {die "canot open $chimera_instructions_file $!\n";}
	print CHIMERA_HTML <<EndOfHTML;
	
<html><title>ConSurf Image using Chimera</title><br />
<body>
<?
include (\"/var/www/html/ConSurf/php/templates/definitions.tpl\");
include (\"/var/www/html/ConSurf/php/templates/output_header_no_refresh.tpl\");
?>
	<center><h3>Create a high resolution Chimera figure</h3></center>
EndOfHTML
	if ($isd_PDB_atoms_ConSurf ne ""){
		print CHIMERA_HTML <<EndOfSection;
	You may view your protein colored according to conservation scores with a unique color for positions of low confidence cutoff. If this is your preferred option, please choose option 1(a). Otherwise - choose option 1(b).
    <br>
    Click <a href = "$conf_link">HERE</a> for more information regarding confidence cutoff.
    <br /><br />		

    	<b><font color="green">You can open the molecule with Chimera on a click!</font></b>
    <ol>
	    <li>(a) <a href = "$url/$isd_chimerax_script" type="application/x-chimerax">Click here</a> for PDB showing Insufficient Data<br />
	    (b) <a href = "$url/$chimerax_script" type="application/x-chimerax">Click here</a> for PDB hiding Insufficient Data</li>
    </ol>
    		
EndOfSection
	}
	else{
		print CHIMERA_HTML "<b><font color=\"green\">You can open the molecule with Chimera <a href = \"$chimerax_script\" type=\"application/x-chimerax\">on a click!</font></b></a><br /><br />";
	}
	
	print CHIMERA_HTML"    <b><font color=\"green\">Alternatively, you can download the files and view it in a later stage.</font></b><br/>\n
	It is recommended to create a new directory where the attached files will be saved.<br>\n\t<ol>";

	if ($isd_PDB_atoms_ConSurf ne ""){

		print CHIMERA_HTML "
		<li>Download the  <b>PDB_FILE</b> updated with ConSurf's colors.<br>\n
		(a) <a href=$isd_PDB_atoms_ConSurf>PDB_FILE</a> showing Insufficient Data<br>\n
		(b) <a href=$PDB_atoms_ConSurf>PDB_FILE</a> hiding Insufficient Data
		</li>\n";
	}
	else{
		print CHIMERA_HTML "<li>Download the <a href=$PDB_atoms_ConSurf><b>PDB_FILE</b></a> updated with ConSurf's colors.</li>\n";
	}
	print CHIMERA_HTML <<EndOfHTML;
		<li>Download the <a href="$chimera_consurf_commands"><b>coloring script</b></a> for chimera.</li>
		<li>Start Chimera program.</li>
		<li>Load the PDB_FILE from the File/Open menu.</li>
		<li>Load the coloring script from the File/Open menu.</li>		
	</ol>
	<br />
	Once the molecule is loaded and coloured, you may save a high quality figure following <a href=$CHIMERA_SAVING_FIGURE_LINK>these instructions</a>.<br />
<?
    include(\"/var/www/html/ConSurf/php/templates/footer.tpl\");
?>
</body></html>

EndOfHTML
	close CHIMERA_HTML;
	return ("OK");
}
#************************************************************************************************
#************************************************************************************************

1;
