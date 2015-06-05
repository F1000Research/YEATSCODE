#!/bin/csh 

if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1 

 
newfile.csh list.pdb
parsePDB_SEQRES.pl  -inf $list -dir RRRR -out list.pdb -or 1 -silent

ls $PDBDIR/ > ! list.5
replacestring.pl -in list.5 -whic "..pdb" -with "" -out ooo -same
TW list.pdb list.5


mv ofhinAbutnotinB list.2getthistime

