#!/bin/csh -f

if($#argv != 4) then 
  echo "Usage : <id> <inputfile> <dist> "
  exit 
endif 

set i = $1
set inputfile = $2
set num = $3
set tag = $4


$SRC/SHELLSCRIPTS/createCLASPinput.csh  $1 $2 $num $num 
ann2simpleinput.pl -out ANNOTATE.$num/$i.in -in ANNOTATE.$num/$i.outconf.annotated

$SRC/MISC/createPremoninput.pl -out $i.$num.$tag.premon.in -con $CONFIGGRP -li ANNOTATE.$num/$i.in -pr $i

