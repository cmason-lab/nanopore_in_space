fi="nine_exp_2d_poretools"
#grep -v '^@' $fi.sam | grep NC_000913.3-Escherichia_coli | cut -f 1 | cut -d':' -f 2 | grep '2d' > ecoli_read_list.txt

#bioawk -cfastx 'BEGIN{while((getline k <"ecoli_read_list.txt")>0)i[k]=1}{if(i[$name])print "@"$name"\n"$seq"\n+\n"$qual}' post_sep_flight.fastq > post_sep_flight_2d_ecoli.fastq
bioawk -cfastx 'BEGIN{while((getline k <"ecoli_read_list.txt")>0)i[k]=1}{if(i[$name])print ">"$name" "$comment"\n"$seq}' ../${fi}.fasta > ${fi}_ecoli.fasta
