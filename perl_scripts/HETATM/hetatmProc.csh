#!/bin/csh -f

if($#argv != 3  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 

set pdb=$1
set hetatm=$2
set ispolar=$3


hetatmProc.pl -outf ooo -l l -het $hetatm -ispolar $ispolar -p $pdb

