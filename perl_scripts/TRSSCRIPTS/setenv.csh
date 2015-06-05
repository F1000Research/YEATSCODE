setenv FASTADIR $cwd/../FASTADIR/
#setenv BLASTOUT $cwd/BLASTOUT_ALLNT/
setenv BLASTOUT $cwd/BLASTOUT/

#setenv BLASTDB /media/sandeepc/USB30FD/NT
setenv BLASTDB /media/sandeepc/USB30FD/PLANTS/ftp.ncbi.nlm.nih.gov/genomes/PLANTS/
setenv BLASTDB /media/sandeepc/USB30FD/NR
setenv BLASTDB $cwd/../

mkdir -p $FASTADIR 
mkdir -p $BLASTOUT 
mkdir -p $BLASTDB 
