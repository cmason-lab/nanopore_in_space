dir='/home/abm237/zeno_abm237/runs_to_date'

for type in qualdist; do #stats
for exp in 20160826; do
 for loc in flight; do
   echo $exp $type
   #ls $dir/${exp}_nasa/$loc/downloads | wc -l
   #ls /$dir/${exp}_nasa/flight/ | head -10
   #poretools $type --type 2D $dir/${exp}_nasa/$loc/downloads > ${exp}_${loc}.$type
   poretools $type $dir/${exp}_nasa/$loc/downloads > ${exp}_${loc}.$type
 done
done
done

for type in qualdist; do #stats
for exp in 20160903 20160907 20160913 20161018 20161025 20161126 20170112; do
 for loc in ground; do
   echo $exp $type
   #ls /$dir/${exp}_nasa/flight/downloads/ | head -10
   for qual in pass fail; do
      #ls $dir/${exp}_nasa/$loc/downloads/$qual | wc -l
      #poretools $type --type 2D $dir/${exp}_nasa/$loc/downloads/$qual >> ${exp}_${loc}.$type
      poretools $type $dir/${exp}_nasa/$loc/downloads/$qual >> ${exp}_${loc}.$type
   done
 done
done
done
