#!/bin/csh -f
alias sp 'set path = ( \!* $path )'

setenv BIOPERLHOME /home/sandeepc/Bio/Code/perl_scripts/perl_scripts/
setenv PRISM $BIOPERLHOME/PRISM

setenv PDB2PQR /home/sandeepc/pdb2pqr/



setenv PDBDIR /home/data/pdbsnamd/
setenv APBSDIR  /home/data/apbsnamd/


setenv PDBDIR /home/data/pdbsdecoy/

setenv PDBDIRDECOY /home/data/pdbsdecoy/
setenv PDBDIR /home/data/pdbs/
#setenv PDBDIR /home/data/pdbs/
setenv PDBDIR ~/PDBMANAGER/PDBDIR/
setenv APBSDIR ~/PDBMANAGER/apbs/
setenv FASTADIR ~/PDBMANAGER/FASTADIR/
setenv HELIXDIR ~/PDBMANAGER/HELIXDIR/
setenv DSSP ~/PDBMANAGER/DSSP/
setenv FPOCKET ~/PDBMANAGER/FPOCKET/
setenv BLASTDB ~/PDBMANAGER/BLASTDB/
setenv BLASTOUT ~/PDBMANAGER/BLASTOUT/

setenv UNIPROT  /home/data/uniprot/
setenv PREMONITION  /home/sandeepc/DATA/CSA/PREMONITION/

setenv CONFIGGRP /home/sandeepc/DATA/data/config.grp 
setenv MAPPINGFILE /home/sandeepc/DATA/data/mapping.txt

setenv ANNDIR  /home/sandeepc/DATA/CSA/ANNOTATE.4/
setenv CACHEPDB  $SRC
setenv MATCH3D  $SRC
setenv RESULTDIR  $SRC



sp $BIOPERLHOME
sp $BIOPERLHOME/APBS
sp $BIOPERLHOME/DELPHI
sp $BIOPERLHOME/EXTERNALTOOLS
sp $BIOPERLHOME/SHELLSCRIPTS
sp $BIOPERLHOME/BIOPERL
sp $BIOPERLHOME/PDBSEQRES
sp $BIOPERLHOME/MISC
sp $BIOPERLHOME/ALIGN
sp $BIOPERLHOME/CSA
sp $BIOPERLHOME/WEB
sp $BIOPERLHOME/IMAGE
sp $BIOPERLHOME/PRIMER
sp $BIOPERLHOME/DECAAF
sp $BIOPERLHOME/BLASE
sp $BIOPERLHOME/BLAST
sp $BIOPERLHOME/NAMD
sp $BIOPERLHOME/FRAGALWEB
sp $BIOPERLHOME/QP
sp $BIOPERLHOME/HETATM
sp $BIOPERLHOME/ENM
sp $BIOPERLHOME/CASP
sp $BIOPERLHOME/DB
sp $BIOPERLHOME/GNM
sp $BIOPERLHOME/RSCRIPTS
sp $BIOPERLHOME/CODON
sp $PRISM
setenv PERLLIB "/home/sandeepc/Downloads/BioPerl-1.6.924:/home/sandeepc/Bio/Code/perl_scripts/perl_scripts/:/home/sandeepc/Bio/Code/perl_scripts/perl_scripts/PRISM:$BIOPERLHOME/PRIMER/:/home/sandeepc/Bio/Code/perl_scripts/perl_scripts/BIOPERL:$BIOPERLHOME/QP:$PERLLIB"



sp /home/sandeepc/DownloadedTools/apbs/
# for multivalue - 
sp /home/sandeepc/DownloadedTools/apbs/apbs-1.2.1-source/tools/mesh/
# for pdb2pqr.py
sp /home/sandeepc/pdb2pqr
