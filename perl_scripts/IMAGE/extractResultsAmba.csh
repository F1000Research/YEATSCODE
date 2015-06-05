#!/bin/csh 
if($#argv != 4  ) then
  echo "Usage : specify resolution "
    exit
	endif

set thresh = $2
set ignoredist = $3
set randomsampling = $4


ff "out.peri" > ! list
$SRC/IMAGE/post3DprocessAmba.pl -li list -outf oooo -resolution $1  -thresh $thresh -ignoredist $ignoredist -randomsampling $randomsampling
frequencyDistribution.pl -out out.all -inf mindist -max 7 -del 1
frequencyDistributionABS.pl -out out.all.abs -inf mindist -max 7 -del 1
$SRC/MISC/makecumu.pl -in out.all


#grep upper list > !  list.upper
#$SRC/IMAGE/post3DprocessAmba.pl -li list.upper -outf oooo
#frequencyDistribution.pl -out out.upper -inf mindist -max 7 -del 1

