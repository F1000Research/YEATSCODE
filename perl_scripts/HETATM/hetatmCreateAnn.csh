#!/bin/csh 

set i=$1 
preparePremonConfigs.pl -outf lll -lis $i.4.1.table.in -pr $i -pol 1 -size 4 -anndir ANNOTATE.4
\cp -f $i.4.1.clasp.in ANNOTATE.4/
createPremoninput.pl -out ANNOTATE.4/$i.4.1.premon.in -con $CONFIGGRP -li ANNOTATE.4/$i.4.1.clasp.in -pr $i
