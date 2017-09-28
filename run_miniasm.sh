#miniasm
name='nine_exp_2d_poretools'
gunzip $name.fasta.gz
head $name.fasta
#minimap -Sw5 -L100 -m0 -t8 $name.fasta $name.fasta | gzip -1 > $name.paf.gz
#miniasm -f $name.fasta $name.paf.gz > $name.miniasm.gfa
#gzip $name.fasta
#gfa to fasta
echo $name.miniasm.gfa
#grep ^S $name.miniasm.gfa | cut -f2,3 | sed -e 's/\t/\n/g' | sed -e 's/utg/>utg/g' > $name.miniasm.fasta
#bioawk -c fastx '{ print $name, length($seq) }' < $name.miniasm.fasta

