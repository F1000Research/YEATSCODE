#!/bin/csh -f
if($#argv != 4  ) then 
  echo "Usage : $#argv ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 

#mkdir -p ~/junk/TMPFASTADIR
#\cp -f $1 $2 ~/junk/TMPFASTADIR

 
$SRC/BLAST/myblastcompare2.pl -p1 $1 -p2 $2 -filesgiven -outf $3 -what $4

