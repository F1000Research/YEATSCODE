#!/bin/csh 

if($#argv != 2  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set i = $1 
extractoneBWA.pl -outf iii -inf $2  -tr $i > ! $i.tex
\rm -f KK.tex 
ln -s $i.tex KK.tex





 

