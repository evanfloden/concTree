#!/usr/bin/perl

use Env;
use strict;
use warnings;
use Algorithm::Cluster qw/kcluster/;

my $fa_dir  = $ARGV[0];
my $a3m_dir = $ARGV[1];
my $numClusters = $ARGV[2];

system("mkdir -p $a3m_dir");
system("mkdir -p hhm_db");

opendir(FA_DIR, $fa_dir) or die $!;
while (my $fname = readdir(FA_DIR)) {
    if ($fname =~ /(\S+).fasta/) {
        system("$HHLIB/scripts/reformat.pl fas a3m $fa_dir/$fname $a3m_dir/$1.a3m.temp");
        system("sed '1s/.*/>$1/' $a3m_dir/$1.a3m.temp > $a3m_dir/$1.a3m");
        system("rm $a3m_dir/$1.a3m.temp");
    }
}
closedir(FA_DIR);

system("$HHLIB/scripts/hhblitsdb.pl -o hhm_db/databaseA -ia3m $a3m_dir");
                                              
opendir(A3M_DIR, $a3m_dir);
while (my $fname = readdir(A3M_DIR)) {
    if ($fname =~ /.a3m/) {
        system("$HHLIB/bin/hhsearch -glob -alt 1 -i $a3m_dir/$fname -d hhm_db/databaseA_hhm_db -id 100");
    }
}
closedir(A3M_DIR);

my %results;

opendir(A3M_DIR1, $a3m_dir);
while (my $fname = readdir(A3M_DIR1)) {
    next unless (-f "$a3m_dir/$fname");
    next unless ($fname =~ m/\.hhr$/);
    open my $result_file, '<', "$a3m_dir/$fname" or die "Couldn't open $fname: $!";;
    my $query;
    my @hits;
    while (my $line=<$result_file>) {
        if ($line =~ /^Query\s+(\S+)/) {
            $query=$1;
        }
        if ($line =~ /^(\s+)?\d+\s+(\S+)\s+\d+.\d+\s+\S+\s+\S+\s+(\d+.\d+)/) {
            my $score = 1/$3; 
            my @hit = [$2,$score];
            push @hits, @hit;
        }
    }
    my @hits_sorted = sort { $a->[0] cmp $b->[0] } @hits;
    $results{$query}=\@hits_sorted;
}
closedir(A3M_DIR1);

my @query_names;
my @query_data;
my @mask;
my @weight;
my $i=0;
my $j=0;
my $size = keys %results;
foreach my $key (sort (keys %results)) {
    $query_names[$i]       =   $key;
    my $data_array_ref     = $results{$key};
    my @data_array         = @$data_array_ref;
    my @query_array;
    for my $k (0..$size-1) {
        push @query_array, $data_array[$k][1];
    }    
    $query_data[$i]=\@query_array;

    $mask[$i]          = [(1)x$size];
    $i++;
}

# Define the paramaters for kcluster
my %params = (
    nclusters=> $numClusters,
    transpose=>            0,
    npass    =>          100,
    method   =>          'a',
    dist     =>          'e',
    data     => \@query_data,
    mask     =>       \@mask,
    weight   =>     \@weight
);

my ($clusters, $error, $found) = kcluster(%params);
my %query_by_rowid;
$i=0;
$query_by_rowid{$i++} = $_, foreach(@query_names);

my %query_by_cluster;
$i=0;
foreach(@{$clusters}) {
    push @{$query_by_cluster{$_}}, $query_by_rowid{$i++};
}

# Print out the first alignment name from the set of alignments
my $repClusters='repCluster.txt';

open(my $file, '>', $repClusters);
for ($i = 0; $i < $params{"nclusters"}; $i++) {
    print $file "$query_by_cluster{$i}[0].fasta\n";
}

exit 0;
