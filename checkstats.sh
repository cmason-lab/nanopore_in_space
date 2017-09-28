bash poretools_extract4stats.sh
for fi in *.stats; do bp=$(grep "total base pairs" $fi | awk -F '\t' '{sum+=$2}END{print sum}'); run=$(echo $fi | cut -d'.' -f1); echo $run $bp  >> bp_per_run.txt; done
