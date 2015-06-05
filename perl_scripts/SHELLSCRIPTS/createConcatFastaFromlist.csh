#!/bin/csh -f


newfile.csh MSA
newfile.csh MSA.list
foreach i ($argv)
    
	echo $i >> MSA.list
	if(! -e $FASTADIR/$i.ALL.1.fasta) then
	echo $FASTADIR/$i.ALL.1.fasta does not exist 
	exit
	endif 
	cat $FASTADIR/$i.ALL.1.fasta  >> MSA	
end 

replacestring.pl -with "" -whi ";" -outf kkkkkk.tex -in MSA -same
#clustalwfromonefile.csh MSA
mafft.csh MSA
seaview MSA.aln & 
