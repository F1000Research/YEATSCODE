#!/bin/csh 

if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 

echo install.packages\(\"$1\", repos=\"http://cran.r-project.org\"\) > ! xxxxxx
sudo R --no-save < xxxxxx
unlink xxxxxx

 

