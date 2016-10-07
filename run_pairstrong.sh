for ws in  90; do
	(
	for run in 0 1 2 3 4 5 6 7 8 9 ; do
		./lamodel -P 2  -T $ws -S 1980$run -s "pairstrongAlt_0_${ws}_${run}"  &
		./lamodel -P 2  -T $ws -S 1980$run -s "pairstrongAlt_${ws}_0_${run}"  &

		./lamodel -P 2  -T $ws -S 1980$run -s "pairstrongAltL_0_${ws}_${run}"  -L &
		./lamodel -P 2  -T $ws -S 1980$run -s "pairstrongAltL_${ws}_0_${run}"  -L &

		#./lamodel -P 2 -w 1 -T $ws -S 1980$run -s "weakstrongLN_0_${ws}_${run}"  -L -n
		#./lamodel -P 2 -w 2 -T $ws -S 1980$run -s "weakstrongLN_${ws}_0_${run}"  -L -n

		#./lamodel -P 2 -w 1 -T $ws -S 1980$run -s "weakstrongN_0_${ws}_${run}"  -n
		#./lamodel -P 2 -w 2 -T $ws -S 1980$run -s "weakstrongN_${ws}_0_${run}"  -n
	done
	) &
done



