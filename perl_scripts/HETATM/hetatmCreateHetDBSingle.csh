#!/bin/csh -f

set hetname = ` cat name` 
mkdir -p ANNOTATE.4

if(! -e "list.50.het") then 
     writeFasta.csh list.all 
     checkIdentity.pl -out list.50 -simi 50 -needleout oooooooooo -arg ~/needle.arg -list list.all
     hetatmPruneListforHet.pl -outf cmds2convert -het2 $PDBDIR/HET2PDB -pdb ~/pdb_seqres.txt -lis list.50 -size 4
endif 


apbs.csh list.50.het
apbs.csh list.50.nohet

echo "Doing for polar"
foreach i (` cat list.50.het `)
	if(! -e ANNOTATE.4/$i.in) then 
        hetatmProc.pl -outf ooo  -het $hetname  -p $i -ispolar 1 
	    \cp -f $i.4.clasp.in ANNOTATE.4/$i.in
	endif
end

source cmds2convert polar

exit 

echo "Doing for non polar"
foreach i (` cat list.50.het `)
    hetatmProc.pl -outf ooo  -het $hetname  -p $i -ispolar 0
	\cp -f $i.4.clasp.in ANNOTATE.4/$i.in
end

source cmds2convert polar
