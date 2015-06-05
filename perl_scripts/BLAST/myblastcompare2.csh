#!/bin/csh -f
if($#argv != 2  ) then 
  echo "Usage : $#argv ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
 
$SRC/BLAST/myblastcompare2.pl -p1 $1 -p2 $2 

