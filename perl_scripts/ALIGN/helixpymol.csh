#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 

echo $1 
helixwheel.pl -out ooooooo -con $CONFIGGRP -aa ~/aalist -prote $1
unlink in
pymol.2oneprotein.csh $1 & 
sleep 1 
kill %%
cat ooo.p1m pymol.in > ! XXX.p1m 
pymol XXX.p1m 
