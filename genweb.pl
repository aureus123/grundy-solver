#!/usr/bin/perl -w
use strict;

# Web graph generator.
#
# (c) 2008 Daniel E. Severin.

print "genweb - Web graph generator\n\n";

sub manual {
	print "Usage: ./genweb.pl n p\n";
	print "Must be n >= 5, p >= 1, n >= 2p.\n";
	exit 1;
}

if ($#ARGV < 1) { manual; }

my $n = $ARGV[0];
my $p = $ARGV[1];
if ($n < 5 || $p < 1 || $n < 2*$p) { manual; }

my $file = "W_" . $n . "_" . $p;
my $file_grafo = $file . ".graph";

open OUT, '>', 'edges.temp' or die "Error creating file: $!";

my $v1;
my $v2;
my $edges = 0;
for ($v1 = 0; $v1 <= $n-2; $v1++) {
	for ($v2 = $v1 + 1; $v2 <= $n-1; $v2++) {
		my $dif = abs($v1-$v2);
		if (($dif > 0 && $dif <= $p) || ($dif >= $n - $p)) {
			print OUT "$v1\n$v2\n";
			$edges++;
		}
	}
}
close OUT;

open IN, '<', 'edges.temp' or die "Error reading file: $!";
print "Output file: $file_grafo\n";
open OUT, '>', $file_grafo or die "Error creating file: $!";

print OUT "$n:$edges\n";
my $e;
for ($e = 1; $e <= $edges; $e++) {
	$v1 = <IN>;
	chomp $v1;
	$v2 = <IN>;
	chomp $v2;
	print OUT "$v1,$v2\n";
}

close OUT;
close IN;

unlink 'edges.temp';

print "W: |V| = $n, |E| = $edges, p = $p\n";
print "Ready!\n";
exit 0;
