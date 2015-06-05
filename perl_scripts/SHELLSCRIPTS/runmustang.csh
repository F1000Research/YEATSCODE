#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 

#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1
set outfile = $list.mustang

createMustangInput.pl -li $list

touch $list.in

~/DATA/mustang/MUSTANG_v3.2.1/bin/mustang-3.2.1 -F fasta -f $outfile

pymol.2oneprotein.pl -out $outfile.p1m -inf $list.in -p results.pdb -formustang

#pymol $outfile.p1m 
