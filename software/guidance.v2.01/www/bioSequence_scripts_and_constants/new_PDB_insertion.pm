package new_PDB_insertion;

#********************************************* #
# PDB package 
# -----------
# contains functions connected to 
# read, write and modify pdb files format
#
# This package returns the array full with data
# of the pdb file in it.
#
# All pdb data starts from @array[1]
#********************************************* #

use strict;
use Vectors;

BEGIN {  }


END {  }


my %AAdata = ( 
                "ALA" => { 
                    "name"    => "A", 
                    "area" => "222" 
                    }, 
                "ARG" => { 
                    "name"    => "R", 
                    "area" => "366" 
                    }, 
                "ASN" => { 
                    "name"    => "N", 
                    "area" => "274" 
                    }, 
                "ASP" => { 
                    "name"    => "D", 
                    "area" => "267" 
                    }, 
                "CYS" => { 
                    "name"    => "C", 
                    "area" => "249" 
                    }, 
                "GLN" => { 
                    "name"    => "Q", 
                    "area" => "300" 
                    }, 
                "GLU" => { 
                    "name"    => "E", 
                    "area" => "294" 
                    }, 
                "GLY" => { 
                    "name"    => "G", 
                    "area" => "193" 
                    }, 
                "HIS" => { 
                    "name"    => "H", 
                    "area" => "309" 
                    }, 
                "ILE" => { 
                    "name"    => "I", 
                    "area" => "301" 
                    }, 
                "LEU" => { 
                    "name"    => "L", 
                    "area" => "301" 
                    }, 
                "LYS" => { 
                    "name"    => "K", 
                    "area" => "336" 
                    }, 
                "MET" => { 
                    "name"    => "M", 
                    "area" => "314" 
                    }, 
                "MSE" => {  #SELENIOMETHIONINE
                    "name"    => "M", 
                    "area" => "314" 
                    }, 
                "PHE" => { 
                    "name"    => "F", 
                    "area" => "329" 
                    }, 
                "PRO" => { 
                    "name"    => "P", 
                        "area" => "257" 
                        }, 
                "SER" => { 
                    "name"    => "S", 
                    "area" => "235" 
                    }, 
                "THR" => { 
                        "name"    => "T", 
                        "area" => "260" 
                        }, 
                "TRP" => { 
                    "name"    => "W", 
                    "area" => "372" 
                    }, 
                "TYR" => { 
                    "name"    => "Y", 
                    "area" => "349" 
                    }, 
                "VAL" => { 
                    "name"    => "V", 
                    "area" => "274" 
                    }, 
                
                ); 


################################################
##############################################

#********************************************* #
sub pdb_to_fasta{
    my $file = shift;  #file to read
    my $chain = shift; #chain to read
    my $pdb = shift;
    my $mode = shift; # 'msa' or 'fasta'
    my $run_name = shift;
    my $OutputHtml = shift;
    my $pdb_data = shift;
    my $pdb_data_insertion = shift;
    my $insertionsListFile = shift;
    my $out_file_name = shift;   # REMARK: add the path to the file, from the general file.
    # REMARK: should be passed as argument
    #my $insertionsListFile = 'insertions.list';
    my $atom_length_file = shift;  # "atom.length"
    my $seqres_length_file = shift;  # "seqres.length"
    my $aa_seq_file = shift; # was $NamePref. ".seq" , should be renamed to "protein_aa.seq"
    my $seqres_pdb_strings_file = shift; # was $NamePref. ".fasta2"
    my $pdb_in_fasta_file = shift; # was $NamePref. ".pdbfasta"
    my $pdb_title_file = shift; # was $NamePref. ".title"
    
    
    my @sequence = (); # to hold the SEQRES sequence
    my @residues = (); # to hold the ATOM sequence
    my $seq_string = "";  # the SEQRES string
    my $atom_string = ""; # the ATOM string
    my ($seq_length, $atom_length, $title, @ans);  
    ########## insertions ##########
    my  $insertionsList = ',';
    ###################################
    
    open PDBFILE, "<$file";
    my @orig_pdb = <PDBFILE>;
    close PDBFILE;
    
    # REMARK: should be sent as file handle
    #open LOG, ">>$OutLogFile";
    #print LOG "entered pdb_to_fasta.pl\n";
    
    my $printwarning = 0; # in order to print the warning only once
    ### read the SEQRES and ATOM sequences (use the PDB module)
    # in case of MSA, there is no need to read the SEQRES
    @sequence = read_SEQRES_pdb($file) if ($mode eq "fasta");
    @residues = read_complex($file, $title);
    my $res_arr;
    $chain = " " if ($chain eq "NONE"); 
    open (PDBDATA, "> " .$pdb_data);   # to PE 
    open (PDBDATAINSERTION, "> " .$pdb_data_insertion);   # to PE         
    # check the ATOM sequence
    @ans = check_sequence("ATOM", \@residues);
    close PDBDATAINSERTION;
    close PDBDATA;
    
    if ($ans[0] eq "ok"){
        $atom_string = $ans[1];
        $atom_length = $ans[2];
    }
    else{
        open ERROR_OUT, $out_file_name;
        print ERROR_OUT "HTML: $ans[1]\n";
        print ERROR_OUT "LOG: ".$ans[2];
        close ERROR_OUT;
        exit;
    }
    
    open (INSERTIONSLISTFILE, ">$insertionsListFile");
    print INSERTIONSLISTFILE "$insertionsList";
    close INSERTIONSLISTFILE;
    #############
    #for PE New Version
    open (ATOMLENGTH, ">".$atom_length_file);   # to PE 
    print ATOMLENGTH $atom_length;
    close ATOMLENGTH;

    my @sequencePE;
    my $seq_string_PE;
    my $seq_length_PE;
    
    @sequencePE = &read_SEQRES_pdb($file);
    ($seq_string_PE, $seq_length_PE) = check_sequence("SEQRES", \@sequencePE);
    
    open (SEQRESLENGTH, ">".$seqres_length_file);
    print SEQRESLENGTH $seq_length_PE;
    close SEQRESLENGTH;

    ##############
    # check the SEQRES sequence
    if ($mode eq "fasta"){    
        ($seq_string, $seq_length) = check_sequence("SEQRES", \@sequence);
    
        # in case there is no SEQRES
        if ($seq_string eq "NO_SEQRES"){
            $seq_string = $atom_string;
            $seq_length = $atom_length;
        }
    }

    ### print the sequences in the related files
    
    if ($mode eq "fasta"){
    
        # the SEQRES file
        open SEQ, ">" .$aa_seq_file;
        print SEQ ">" .$pdb . "_" . $chain . "  " . $seq_length . "\n" . $seq_string . "\n" ;
        close SEQ;
    
        # the SEQRES and ATOM file
        open OUT, ">" .$seqres_pdb_strings_file;
        print OUT ">SEQRES_$chain\n";
        print OUT "$seq_string\n\n";
        print OUT ">PDB_$chain\n";
        print OUT "$atom_string\n";
        close OUT;
    }
    
    # the ATOM file
    open PDBFASTA, ">" .$pdb_in_fasta_file;
    print PDBFASTA ">PDB_$chain\n";
    print PDBFASTA "$atom_string\n";
    close PDBFASTA;
    
    # the title file
    open TITLE, ">" .$pdb_title_file;   
    print TITLE $title;
    print TITLE "\nCHAIN: $chain\n";
    close TITLE;

#********************************************* #
# sub read_sequence_pdb($) 
# ------------------
# reads a pdb sequence file containg
# into @residues array, on a 1LETTER format
# This array will contain an entry per each 
# amino acid, containing all the relevant data
#********************************************* #
sub read_SEQRES_pdb($) {
 
    my $file = $_[0];   #the file to be opened and read
    my $chain; #the chain
    my $ResList; #list of residue names.

    my @names3L = (); #here are written the names in 3L   
    my $resCount = 0; 
    my $AA3L;
    my @residues = ();
    my @sequence = (); #main data to be returned

    open (PDB, $file) || die "cannot open file = $file";
    
    while (<PDB>) {
        
        if (/^SEQRES/){
            
            $chain = substr ($_, 11, 1);
            $ResList = substr ($_, 19, 51);    
            @names3L = split (/\s+/, $ResList);
      
            #*************************
            #make data structure like in read_complex
 
            foreach $AA3L (@names3L) {
                
                $sequence[$resCount]{name1L} = $AAdata{$AA3L}{name};   #residue name 1L
                $sequence[$resCount]{name3L} = $AA3L; #residue 3L name
                $sequence[$resCount]{chain} = $chain;  #chain        
                $resCount++;
            }
        }
    }
    close (PDB);    
  
    return (@sequence); #returns entire array with data
}

##############################################
##############################################
#********************************************* #
# sub read_complex($) 
# ------------------
# reads a pdb file containg
# a pdb complex (more than one protein chain
# into @residues array.
# This array will contain an entry per each 
# amino acid, containing all the relevant data
#********************************************* #

sub read_complex($$) {
    
    my $resCount = 0; 
    my @residues = ();
    my $i;
    my $file = $_[0]; #the file to be opened and read
    my $Res3L; #residue before 3letter validation (AARG/ BARG)
    my $chain;
    my $number;
    my $iCode; # insertion code
    my $insCode='';
    my $atomNo;
    my $chainInsCode = '';
    $residues[0]{number} = -300; #to initialize this scalar

    open (PDB, $file) || die "cannot open file = $file";    
    
    while (<PDB>) {        
        if (/^ATOM/){
  ############# insertion begin ###########
	    $atomNo .= substr ($_, 7, 5);
	############# insertion end ###########
            $Res3L = substr ($_, 17, 3);
            $chain = substr ($_, 21, 1);
            $number = substr ($_, 22, 4);
	    $number =~  s/\s+//g;
	    $iCode = substr ($_, 26, 2);
	 print "AAAAAAACCDD   \$Res3L = $Res3L,  \$chain = $chain, \$number = $number$iCode \n"; 
	    if ($iCode !~ /[A-Z]/i){
		$chainInsCode = '';
	    }
########################$$$$$$$$$$$$$$$$$$###########################
	    $number =~ s/\s+//;
	    if (($iCode =~ /[A-Z]/i) && ($chainInsCode !~ /$iCode/))
	    {
		$number =~  s/\s+//g;		
	        $insCode .= $Res3L . "," . $chain . "," . $number . "," . $iCode  . ";";	
		$chainInsCode =  $insCode;
		$chainInsCode =~ s/.*,(.)/$1/;
	    }
            ($number) =~ s/\s+//g;
  ############# insertion begin ###########
	    ($atomNo) =~ s/\s+//g;
             $atomNo .=','; 
	############# insertion end ###########
	    $iCode =~ s/\s*//g;
            my $numberAndInsertion = $number . $iCode; 
	    if ($numberAndInsertion ne  ($residues[$resCount - 1]{number})) { #create residue hash
        ############# insertion begin ###########
                $residues[$resCount-1]{atomNo} = $atomNo;
                $residues[$resCount-1]{atomNo} =~ /^(\d+),.*\,(\d+)\,$/;
                my $leftAtomNo = $1-1; #?
                my $rightAtomNo = $2-1; #?
                $residues[$resCount-1]{atomNo} = '(atomno >= ' . $leftAtomNo . ' and atomno <= ' . $rightAtomNo . ')';
	############# insertion end ###########
           
        ##############insertion begin #######################
		if ($iCode =~ /[A-Z]/i){    
                    $insCode =~ s/\;$//;
		    $residues[$resCount]{name1L} = $AAdata{$Res3L}{name};   #residue name 1L           
		    $residues[$resCount]{name3L} = $Res3L;   #residue name 3L
		    $residues[$resCount]{chain} = $chain;  #residue chain
		    $residues[$resCount]{number} = $number . $iCode;
                    $resCount++;
		    $insCode='';
		    $atomNo ='';
	        }
                else
		{
		    $residues[$resCount]{name1L} = $AAdata{$Res3L}{name};   #residue name 1L           
		    $residues[$resCount]{name3L} = $Res3L;   #residue name 3L
                    $residues[$resCount]{chain} = $chain;  #residue chain
		    $residues[$resCount]{number} = $number; #residue number
		    $resCount++;
		    $atomNo ='';
		}
                ############### insertion end #####################
            }
        }
        # READ TITLES
        elsif (/^TITLE(.+)/){
            $_[1] .= $1;
        }
        # for NMR files - read only one model
        elsif (/^ENDMDL/){
            last;
        }
    }
    close (PDB);    
    return (@residues); #returns entire array with data
}


#********************************************* #
# sub de_reference($) 
# ------------------
# goes over all pairs of residues on\@array
# of a pdb complex (more than one protein chain
# into @residues array).
# This sub returns the 
#********************************************* #
sub de_reference ($) {
    my $first_res = $_[0];
}

#********************************************* #
# sub all_pairs($) 
# ------------------
# goes over all pairs of residues on\@array
# of a pdb complex (more than one protein chain
# into @residues array).
# This sub returns the 
#********************************************* #

sub all_pairs ($) {

    my $array = $_[0];
    my $first = 1;
    my $second = 1;
    my @pair;

    for ($first = 1; $first < @$array; $first++){

        for ($second = $first + 1; $second < @$array; $second++) {
            Vectors::distance (@$array[$first]->{atoms}{CA}, @$array[$second]->{atoms}{CA});
        }
    }
}
#********************************************* #
# sub print_pdb_1L (@) 
# ------------------
# goes over residues on @array{name1L}
# of a pdb complex and prints its 1L fasta format
#********************************************* #

sub print_pdb_1L (@) {

    my @residues = @_;
    my $AA;
}

#############################################################
# Description:  replace the tempFactor column in the PDB file
#               with another column. (for example: grades)
# Arguments: 1. PDB file
#            2. chain
#            3. reference to hash
#               key - the residue sequence number
#               value - the field to insert instead of the tempFactor
#############################################################
sub insert_tempFactor {

    my $PdbFile = shift;
    my $chain = shift;
    my $HashRef = shift;

    my $PDBchain;
    my $ResNum;
    my $iCode;
    
    # read the PDB file to an array
    open READPDB, "<$PdbFile";
    my @PDB = <READPDB>;
    close READPDB;

    # write the PDB file and replace the tempFactor column
    # with the new one.
    open WRITEPDB, ">$PdbFile";
    
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
                substr ($line, 60, 6) = $$HashRef{$ResNumiCode};
                print WRITEPDB $line;
            }
            else {
                print WRITEPDB $line;
            }
	}
	else {
	    print WRITEPDB $line;
	}
    }
    close WRITEPDB;
}
################################################################
# check if the SEQRES or ATOM sequence exists and contain standard amino-acids, 
# and write the file prefix.pdbdata
# returns $atom_string, $atom_length
sub check_sequence {

    my $field = shift; # 'ATOM' or 'SEQRES'
    my $array = shift;
    my $OutputHtml  = shift;
    my $chain_link = shift;

    my $seq = "";
    my $anyAA = 0;  # to check if there is a sequence at all
    my $CountAA = 0; # to count the AAs of the given chain
    my $StandardAA = 0;
    my %uniqueInsertion;
    my $html_error;
 
    foreach my $AA (@{$array}) {
	$anyAA++ unless ($$AA{name3L} =~ /^\s*$/);
        #for chain $chain only
        if ($$AA{chain} eq $chain) { 	    
            $CountAA++;	    
	    if ($$AA{name1L} ne ""){ # a standard AA
		$seq .= $$AA{name1L}; #AA to print
		$StandardAA++;
		if ($field eq "ATOM"){
###############insertion begin #####################
		    print PDBDATA "$$AA{name3L}$$AA{number}:$$AA{chain}\n";#PDBDATA was opened by the function who called the check_sequence procedure
		    my $uniqueIns;
		    if ($$AA{number} =~ /[A-Z]$/){
			my $insertionNumber = $$AA{number};
			$insertionNumber =~ s/[A-Z]$//;
			$insertionsList .= "$insertionNumber\," if ($insertionsList !~ /,$insertionNumber\,/);		
		    }	
		    print PDBDATAINSERTION "$$AA{name3L}$$AA{number}:$$AA{chain},$$AA{name3L}$$AA{number},$$AA{atomNo}\n";			    
############### insertion end ######################		   
		}
	    }
        }
    }
    ### There is no sequence at all
    if ($anyAA == 0){
	# no SEQRES
	if ($field eq "SEQRES"){
	    $seq = "NO_SEQRES";            
            if ($printwarning == 0){
                open OUTPUT, ">>$OutputHtml";
                print OUTPUT "\n<p><ul><li>Warning: There is no SEQRES derived information in the PDB file. The calculation will be based on the ATOM derived sequence. If this sequence is incomplete, we recommend to run ConSurf again using an external multiple sequence alignment file, which is based on the complete protein sequence.</li></ul></p>\n";
                close OUTPUT;
                $printwarning =1;
            }
	    return "ok", $seq, $StandardAA;
	}
	# no ATOM
	else {
            $html_error = "<b>There is no $field derived information in the PDB file.<br>Please refer to the OVERVIEW for detailed information about the PDB format.</b></p>\n";
            my $log_error = "\npdb_to_fasta_insertion.pl: exit - There is no $field field.\n";

            return "error", $html_error, $log_error;
	}
    }

    ### The given chain is not correct
    if ($CountAA == 0){

	# the given chain was "none" 
	if ($chain eq " "){
	    if ($pdb eq "FILE"){  # an uploaded PDB file
		$html_error ="<b>The <a href=$chain_link>chain identifier</a> column in the $field field of the <a href=$file target=PDB>PDB file</a> is not empty.<br>Please check your file or run ConSurf again with a specific chain identifier.</b></p>\n";
		
	    }
	    else {
		$html_error ="<b>The <a href=$chain_link>chain identifier</a> column in the <a href=$file target=PDB>PDB file</a> is not empty.<br>Please run ConSurf again with a specific chain identifier (A, B, C, etc.)</b></p>\n";
	    }
	}	
	else {
	    if ($pdb eq "FILE"){  # an uploaded PDB file

		$html_error = "<b>Chain \"$chain\" does not exist in the $field field of the <a href=$file target=PDB>PDB file</a>.<br> Please check your file or run ConSurf again with another chain. If there is no <a href=$chain_link>chain identifier</a>, set the 'Chain Identifier' field to \"none\".</b></p>\n";
	    }
	    else {
		$html_error = "<b>Chain \"$chain\" does not exist in the <a href=$file target=PDB>PDB file</a>.<br> Please run ConSurf again with another chain. If there is no <a href=$chain_link>chain identifier</a>, set the 'Chain Identifier' field to \"none\".</b></p>\n";
	    }
	}
	my $log_error = "\npdb_to_fasta_insertion.pl: exit - the chain is not correct\n";
        return "error", $html_error, $log_error;
    }

    ### The chain exist but it is not a protein
    if ($StandardAA == 0){
	if ($pdb eq "FILE"){  # an uploaded PDB file
	    $html_error = "<b>Chain \"$chain\" in the $field field of the <a href=$file target=PDB>PDB file</a> does not contain any standard amino acids.<br>Please check your file or run ConSurf again with another chain.</b></p>\n";
	}
	else {
	    $html_error = "<b>Chain \"$chain\" in the $field field of the <a href=$file target=PDB>PDB file</a> does not contain any standard amino acids.<br>Please run ConSurf again with another chain.</b></p>\n";
	}	
	my $log_error = "\npdb_to_fasta_insertion.pl: exit - the chain is not a protein.\n";
        return "error", $html_error, $log_error;
    }
    return "ok", $seq, $StandardAA;
}
################################################################
}

1;






