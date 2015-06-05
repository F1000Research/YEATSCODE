#procRaf.pl -in DifExpProteinsTable.csv -out smalltable -control 10 -control 11 -control 12 -treated 13 -treated 14 -treated 15


newfile.csh annotations.txt
foreach i (`cat $1`)
	if(-e $BLASTOUT/$i.blast.nt) then
       ~/parsewebblast.pl -in $BLASTOUT/$i.blast.nt -outf annotations.txt
	endif 
end

#extractindexfromfile.pl -in report.csv -out smalltable.csv -id 0 -idx 3 -idx 4 -idx 5 -idx 6 -idx 7 -idx 8 -sep "   "


#$SRC/MISC/needleFastalistSaveParse.pl -in datafile -outf common -cutof 0.00000000001


#$SRC/MISC/runProbislistall2all.pl -outf UUU -lis list.pdbs -tmp ttt


#../procSame.pl -inf annotations.txt.sort -out out.cnt





pruneSameSequenceFromMadeFasta.pl -outf list.transcriptome.merged -infi ../ORFbestFasta.fasta -ignore warning.notanno

mergeEquivalent.pl -mapp list.transcriptome.clean.merged.mapping -outf ../raw_bwa_counts_by_trs.txt.normalized.merged -inf ../raw_bwa_counts_by_trs.txt.normalized -ignoref warning.notanno

mergeEquivalent.pl -mapp list.transcriptome.clean.merged.overlap.mapping.f -infile ../raw_bwa_counts_by_trs.txt.normalized.merged -out ../raw_bwa_counts_by_trs.txt.normalized.merged.overlap -ignoref warning.notanno -mergenames

mergeEquivalent.pl -mapp data.TRS2TRS.2.groups.sort -infile ../raw_bwa_counts_by_trs.txt.normalized.merged.overlap -out ../raw_bwa_counts_by_trs.txt.normalized.merged.overlap.final -ignoref warning.notanno

../preprocessBWA.pl -inf ../raw_bwa_counts_by_trs.txt.normalized.merged.overlap.final -ignoref warning.notanno -outf trs36k

sort.pl -idx 3 -in trs36k -rev -out trs36k.mean -cutoff 100
sort.pl -idx 5 -in trs36k.mean -rev -out trs36k.mean.sd


../extractoneBWA.pl -outf iii -inf ../raw_bwa_counts_by_trs.txt.normalized -tr $i > ! $i.tex


../normalizeRawCountsUsingMeanSD.pl -igno warning.notanno -map meanSDall.names -inf ../raw_bwa_counts_by_trs.txt.normalized.merged.overlap.final

./preprocessBWA.pl -inf ../raw_bwa_counts_by_trs.txt.normalized.merged.overlap.final.usingMeanSD -ignoref warning.notanno -outf trs36k.new

../findPWDistance.pl -ignore warning.notanno -inf ../raw_bwa_counts_by_trs.txt.normalized.merged.overlap.final.usingMeanSD




~/parseblastReflection.pl -outf data.TRS2TRS -inf BLASTOUT_TRS2TRS/$i.blast.merged
~/parseblastReflection.pl -outf data.TRS2TRS.1 -inf BLASTOUT_TRS2TRS/$i.blast.out
~/parseblastReflection.pl -outf data.TRS2TRS.2 -inf BLASTOUT_TRS2TRS.2/$i.blast.out




../parseTRS2TRS.pl -outf data.TRS2TRS.2.mapping -in data.TRS2TRS.2 -tag BLASTOUT_TRS2TRS.2
UNIQ hasmorethanone
wc hasmorethanone.uniq onlyone  ## this is same as the length
groupBasedonCutoff.pl -outf data.TRS2TRS.2.groups -inf data.TRS2TRS.2.mapping -cutoff 0.001 -dir 0
cat onlyone data.TRS2TRS.2.groups >  ! list.transcriptome.finalgenes
sort.pl -idx 3 -in data.TRS2TRS.2.groups
../processMappingForSamePrefix.pl -in data.TRS2TRS.2.groups.sort -out jjjjj





../joinTRS.pl -inf ../DB/DB.list.transcriptome.merged.TRS2TRS -outf data.TRS2TRS.OVERLAP7 
../joinTRSPostProcess.pl -inf data.TRS2TRS.OVERLAP7 -size 7 -anno annotations.txt

The above 2 are combined for by:
source ../joinTRS.csh ../DB/DB.list.transcriptome.merged.TRS2TRS 
(the result of the above is in a dir called MERGED.1 - this can be done iteratively, the second iteration in the DIR MERGED)
cat MERGED.1/*include > ! INC
cat MERGED.1/*exclude > ! EXC


groupBasedonCutoff.pl -outf XXXX -inf list.transcriptome.clean.merged.overlap.mapping -cutoff 0.001 -dir 0
sort.pl -idx 3 -in XXXX
../chainChainedTRS.pl -outf list.transcriptome.clean.merged.overlap.mapping.final -in list.transcriptome.clean.merged.overlap.mapping -list XXXX.sort 

../findRevInOrfs.pl -outf rev.TRS -list list.transcriptome.clean


extractslicefromfasta.pl -out kkkk -inf SKDH.super790.fastarev.fasta -sta 42000 -end 46795



extractAlignFromBlastout.pl -in jjjj -outf xxxx -query ../../FASTADIR/C16577_G1_I1.ALL.1.fasta -sub jcf7180001222206.ALL.1.fasta.comp.fasta -tag C16577_G1_I1
extractAlignFromBlastout.pl -in BLASTOUT_2WGS/C45998_G1_I1.blast.nt -outf xxxx -query ../../FASTADIR/C45998_G1_I1.ALL.1.fasta -sub super778.ALL.1.fasta.comp.fasta -tag C45998_G1_I1


myblastcompare2Fastafiles.csh ../FASTADIR/$i.ALL.1.fasta super778.ALL.1.fasta.comp.fasta jjjj
extractAlignFromBlastout.pl -in jjjj -outf xxxx -query ../FASTADIR/C45998_G1_I1.ALL.1.fasta -sub super778.ALL.1.fasta.comp.fasta -tag C45998_G1_I1 -doint 0


extractsinglefastafromfile.pl -inf ../../DB/wgs.5d.scafSeq -wh jcf7180001222202


~/generateScripts.pl -outf kkkk -in SKDH.TRS.ALL.anno




../chainChainedTRS.pl -outf list.transcriptome.clean.merged.overlap.mapping.f -in list.transcriptome.clean.merged.overlap.mapping -list XXXX.sort -app list.transcriptome.clean.merged.overlap

../normalizeRawCounts.pl -igno warning.notanno -out iii -map map.TRS2length -inf ../raw_bwa_counts_by_trs.txt



../parseWebBlastOuts.pl -outf jjj -tr $i -cutoff 0.00001 -chooseone

../../iterativeRunFindGene.csh

processReannotate.pl -out oooo -inf infointrons.csv

