#!/usr/bin/env perl

use Getopt::Long;
use List::Util 'shuffle';
use Bio::AlignIO;
use Bio::Align::Utilities qw(:all);
use Bio::LocatableSeq;

my ($help, $intersect, $random, $size, $replicate, $output, @files, @alns);
my $MAX_ALNS=2000;

parse_options();
preprocess();
parse_alignment_file();
if(!defined $intersect)
{
	add_empty_seq_4_union();
}
concatenate_output();

sub concatenate_output
{
	$con_aln = cat(@alns);
	if(!defined $output)
	{
	  $out = Bio::AlignIO->new(-fh => \*STDOUT);
	}
	else
	{
	  $format = "fasta";
	  if($output =~ m/[^\.]*.phylip/ )
	  {
	    $format = "phylip";
	  }  
	  $out = Bio::AlignIO->new(-file => ">$output" , '-format' => $format);
	}	  
	$con_aln->set_displayname_flat(); #Makes all the sequences be displayed as just their name, not name/start-end
	$con_aln->uppercase(); #Sets all the sequences to uppercase
	$out->write_aln($con_aln);
}

sub union
{
	my ($arr_ref1, $arr_ref2) = @_;
	my %union = ();
	foreach(@{$arr_ref1},@{$arr_ref2})
	{
		$union{$_}=1;
	}	
	@{$arr_ref1} = keys %union;
}

sub parse_alignment_file
{
	for $file (@files)
	{
 	  	$str = Bio::AlignIO->new('-file' => $file);
  		$aln = $str->next_aln();
		push(@alns, $aln);
	}
}

sub add_empty_seq_4_union
{
	@aln_seq_arr = ();
	@output_name_array = ();
	foreach $alnObj (@alns)
	{
		@seq_names = ();
		foreach my $seqObj ($alnObj->each_seq)
		{
			push(@seq_names, $seqObj->display_id);
		}
		push @aln_seq_arr, [ @seq_names ];

		if (scalar(@output_name_array) < 1)
		{
			@output_name_array = @seq_names;
		}
		else
		{
			union(\@output_name_array, \@seq_names);
		}
	}
	for($i = 0; $i < scalar(@aln_seq_arr); $i++)
	{
		my %hash=map{$_=>1} @{$aln_seq_arr[$i]};
		my @miss_arrary=grep(!defined $hash{$_}, @output_name_array);
		if(scalar(@miss_arrary) > 0)
		{
			add_seq_2aln(\$alns[$i], @miss_arrary);
		}
	}
}

sub add_seq_2aln
{
  	$empty_character = "X";
	my ($alnObj_ref, @add_seq_arr) = @_;
	
  	$miss_seq = $empty_character x $$alnObj_ref->length;
	foreach (@add_seq_arr)
	{
		 $seq = new Bio::LocatableSeq(-seq => $miss_seq, -id  => $_);
  		 $$alnObj_ref->add_seq($seq);
	}
}

sub usage
{
#-- prints usage if no command line parameters are passed or there is an unknown
#   parameter or help option is passed
  	print "Unknown option: @_\n" if ( @_ );
	print "\nUnion concatenate alignments according to inputting order (default model),";
	print "\nBioPerl is need for running concatenate.pl";
	print "\n\nUSAGE:";
	print "\nconcatenate.pl [options] --aln alignment1 alignment2 ... (at least two alignments) --out file";
	print "\n\t--intersect: only concatenate sequences happening in all alignments";
	print "\n\t--random: concatenate order is random";
	print "\n\t--size SIZE: only concatenate x alignments from inputting, exclusive with replicate mode";
	print "\n\t--replicate NUM: replicate one alignment in x times, exclusive with size mode";
	print "\n";
	print "\nrunning example and output message";
	print "\n-------------------------------------------------------------------------------------";
	print "\n../bin/concatanate.pl --aln Dmel_56.phylip_aln Dana_34.phylip_aln Dere_43.phylip_aln";
	print "\nHow to concatenate:";
	print "\n\tsize: 3 alignments";
        print "\n\talignments:";
        print "\n\t\tDmel_56.phylip_aln";
        print "\n\t\tDana_34.phylip_aln";
        print "\n\t\tDere_43.phylip_aln";
	print "\n-------------------------------------------------------------------------------------";
	print "\n";
	exit;
}

sub parse_options
{
	#no input options
	if((!@ARGV) or	!GetOptions('help|?'=>\$help,'intersect'=>\$intersect,'random'=>\$random,'size=i'=>\$size,'replicate=i'=>\$replicate,"aln=s{1,$MAX_ALNS}"=>\@files,'out=s'=>\$output) 
	or defined $help 
	or ($#ARGV ne -1)) #after parse options, there are remaining options
	{
		print "\n[ERROR] unkonwn options = @ARGV\n";
		usage();
	}
	else
	{
		if ((defined $size) and (defined $replicate)){
			print "\n[ERROR] size and replicate modes are exclusive\n";
			usage();
		}
		else{
			if ((!defined $size) and (!defined $replicate)){
				$size = scalar(@files);
			}
			elsif(defined $size){
				if($size < 2){
					print "\n[ERROR] specific size = $size, should be larger than 2\n";
					exit;
				}
				elsif ($size > scalar(@files)){
					print "\n[ERROR]size = $size > # of alignments = ".scalar(@files)."\n";
					exit;
				}
			}
			elsif(defined $replicate){
				if(scalar(@files) > 1){
					print "\n[ERROR]replicate, # of alignments = ".scalar(@files)." > 1\n";
					exit;
				}
			}
		}
	}
}

sub preprocess
{
	@files = shuffle(@files)  if (defined $random);
	if(defined $size){
		@files = @files[0..($size-1)] if($size < scalar(@files));
	}
	elsif(defined $replicate){
		$single_file=$files[0];
		@files = ($files[0]) x $replicate;
	}
	print "\nHow to concatenate:";
	print "\n\tintersect: only concatenate sequences happen in all alignments" if(defined $intersect);
	print "\n\torder: random" if (defined $random);
	print "\n\tsize: $size alignments" if (defined $size);
	print "\n\treplicate: $replicate times" if (defined $replicate);
	print "\n\talignments:";
	foreach (@files)
	{
		print "\n\t\t$_";
	}
	print "\n\toutput: $output";
	print " \n";
}
