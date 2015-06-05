#!/bin/csh 
if($#argv != 5  ) then 
  echo "exec rawcount outfile cutoff idx1 idx2"
  exit 
endif 

set rawcount=$1 
set outfile=$2
set cutoff=$3
set idx1=$4
set idx2=$5
extractforR.pl -in $rawcount -outf $outfile -cutoff $cutoff -idx $idx1 -idx $idx2
echo  "wilcox.test(X, Y, correct = F)" >> $outfile
R --no-save < $outfile > & ! XXXX
