

function la_run {
	qsub -v "LAPARAMS=$LAPARAMS" submit_lamodel.sh
	#./lamodel $LAPARAMS  
	#sleep 1
}


function run_weaks {
		LAPARAMS="$WEAKS -P 3  -T $ws -S 191$run -s three${SUFF}_${ws}_${run}  "
		la_run 

		LAPARAMS="$WEAKS -P 3  -n -T $ws -S 191$run -s three${SUFF}N_${ws}_${run}  "
		la_run 

		LAPARAMS="$WEAKS -P 3  -T $ws -S 191$run -s three${SUFF}G_${ws}_${run}  -G"
		la_run 

		LAPARAMS="$WEAKS -P 3  -T $ws -S 191$run -s three${SUFF}L_${ws}_${run}  -L"
		la_run 

}

for ws in 60  120 180 300  1440;  do
	for run in  {0..10}; do

	WEAKS="-w 1 -w 2 -w 3"
	SUFF="www"
	run_weaks 

	WEAKS="-w 1  -w 3"
	SUFF="wsw"
	run_weaks 

	WEAKS="-w 2 "
	SUFF="sws"
	run_weaks 

		#LAPARAMS=" -P 2  -T $ws -S 1980$run -s "three2N_${ws}_${run}"   -n"
	#	qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
	#	sleep 1
	#	LAPARAMS=" -P 2  -T $ws -S 1980$run -s "three2NL_${ws}_${run}"  -n -L "
	#	qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
	#	sleep 1

		#LAPARAMS=" -P 2  -T $ws -S 1980$run -s "three2Alt_${ws}_${run}" -C "
		#qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		#LAPARAMS="-P 2  -T $ws -S 1980$run -s "three2LAlt_${ws}_${run}"  -C -L"
		#qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		#LAPARAMS=" -P 2  -T $ws -S 1980$run -s "three2NAlt_${ws}_${run}"  -C -n"
		#qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		#LAPARAMS=" -P 2  -T $ws -S 1980$run -s "three2NLAlt_${ws}_${run}" -C -n -L "
		#qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		#sleep 1
	done
done



