fasta='nine_exp_2d_poretools'
cut -f1,2 *.sam | grep -v '^@' |  grep $'\t'4  | grep '2d' | cut -f 1 > nonhost_list.txt
bioawk -cfastx 'BEGIN{while((getline k <"nonhost_list.txt")>0)i[k]=1}{if(i[$name])print ">"$name" "$comment"\n"$seq}' ${fasta}.fasta > ${fasta}_nohost.fasta
