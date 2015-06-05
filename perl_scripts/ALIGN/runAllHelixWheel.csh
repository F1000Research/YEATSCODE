#!/bin/csh -f


setenv PDBDIR $cwd

set list = allnames

newfile.csh allvalues

foreach i (`cat $list`)
echo -n "$i  "
helixwheel.pl -out ooooooo -con $CONFIGGRP -aa ~/aalist -prote $i > & ! /dev/null
end 

echo 
