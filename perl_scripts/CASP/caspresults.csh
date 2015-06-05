#!/bin/csh -f


setenv APBSDIR APBSDIR/
setenv PDBDIR $cwd

ln -s ~/pd.* .
if(! -e oooo) then 
    escapist.pl -outf oooo -con $CONFIGGRP -lis list -score pd.CB.score
    distResidues.pl -outf xxxx -con $CONFIGGRP -lis list
endif

sort.pl -in oooo -idx 5 -out results
extractindexfromfile.pl -in results -idx 21 -out ttt
cat ttt | head -1 > ! results.escapist
set i=`cat results.escapist`
$SRC/ALIGN/specfpr.pl -infil results -li list -idx 5 -outf qqqqqq -tag aaa


grep $i xxxx 
cat spec.5.aaa
