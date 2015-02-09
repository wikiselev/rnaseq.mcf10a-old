#!/usr/bin/perl
use warnings;
use strict;

chdir('../../miRNA_raw_data/') or die "$!";

my @sam_files = <*sam>;

my %counts;

for my $index (0..$#sam_files) {
    my $sam_file = $sam_files[$index];

    warn "Processing $sam_file\n";

    open (IN,$sam_file) or die $!;

    while (<IN>) {
	next if (/^@/);

	my (undef,undef,$mirna) = split(/\t/);
	next if ($mirna eq '*');

	$counts{$mirna}->[$index]++;
    }
}

open (OUT,'>','summarised_mirna_counts.txt') or die $!;

print OUT join("\t",("miRNA",@sam_files)),"\n";

foreach my $mirna (keys %counts) {

    my @counts = @{$counts{$mirna}};

    for (0..$#counts) {
	$counts[$_] = 0 unless ($counts[$_]);
    }

    for ($#counts+1..$#sam_files) {
	$counts[$_] = 0;
    }

    print OUT join("\t",($mirna,@counts)),"\n";
} 


close OUT or die $!;