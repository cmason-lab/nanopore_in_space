ref='/zenodotus/masonlab/epitranscriptome_scratch/abm237/REFS/fastas/mm10_hg38.fa'
for fi in "$@"; do 
   name=$(echo $fi | cut -d'.' -f 1)
   echo $name
   graphmap align -r $ref -d $fi -o $name.sam
done
