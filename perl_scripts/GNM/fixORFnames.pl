#!/bin/csh 

if($#argv != 2  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1 
set DB = $2 

newfile.csh $DB

mkdir TMP
\rm TMP/*

foreach i (`cat $list`)

    echo ">$i" > ! TMP/$i.ALL.1.fasta
    cat FASTADIR_ORFLATEST/$i.ALL.1.fasta | grep -v ">" >> TMP/$i.ALL.1.fasta
    cat TMP/$i.ALL.1.fasta >> $DB
    \cp -f TMP/$i.ALL.1.fasta FASTADIR_ORFLATEST/$i.ALL.1.fasta
    
end 


pruneSameSequenceFromMadeFasta.pl -outf $list -inf $DB
 

newfile.csh $DB
foreach i (`cat $list`)

    cat FASTADIR_ORFLATEST/$i.ALL.1.fasta >> $DB
    
end 

