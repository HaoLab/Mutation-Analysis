#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(dirname basename);
use Parallel::ForkManager;
use Data::Dumper;

unless (@ARGV==1) {
	print "perl\t$0\t<sh>\n";
	exit 1;
}
my $pm = new Parallel::ForkManager(5);



open IN,$ARGV[0] or die $!;
while(<IN>){
	chomp;
	$pm->start and  next;
	lumpy($_);
	$pm->finish;
}

$pm->wait_all_children;



sub lumpy{
	my($bam) = shift @_;
	system($bam);
}







