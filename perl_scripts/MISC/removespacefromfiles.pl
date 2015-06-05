#!/usr/bin/perl -w 
my @l = <*png>;
foreach my $i (@l){
  my $orig = $i ;
  $i =~ s/\s*//g;
  system ("mv \"$orig\" $i");
}
