#!/usr/bin/env perl

use strict;
use v5.10;
$/ = ">";
while (<>){
#    say "one";
    chomp;
    m/^\S+\n(\S+)\n/ms;
    my $firstLine = $1;
    $firstLine =~ s/-//g;
    say $firstLine if length $firstLine >= 26;
}
