#!/bin/csh -f
mv pot0-PE0.dx.atompot pot0.dx.atompot
if(-e pot0.dx.atompot) then 
    ln -s  pot0.dx.atompot pot1.dx.atompot
endif 
#mv pot1-PE0.dx.atompot pot1.dx.atompot
