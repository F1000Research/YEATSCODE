#!/bin/csh 
if($#argv != 2  ) then 
  echo "exec> rawcount outfile trs"
  exit 
endif 

set rawcount=$1 
set outfile=$2
extractSinglelineforR.pl -in $rawcount -outf $outfile 
#echo  "wilcox.test(X, Y, correct = F)" >> $outfile

#R --no-save < $outfile > & ! QUANTILE/$trs.quant
