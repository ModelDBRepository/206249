for ws in 1 5 10 20 30 60 90 120 180 240; do
	(
	for run in 0 1 2 3 4 5 6 7 8 9 ; do
		./lamodel -P 2  -T $ws -S 1980$run -s "weak2_${ws}_${run}"  
		./lamodel -P 2  -T $ws -S 1980$run -s "weak2L_${ws}_${run}"  -L
	done
	) &
done



