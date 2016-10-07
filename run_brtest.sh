
for ws in 10 20 30 40; do
	for run in  {0..10}; do

		LAPARAMS=" -T 180 -J -B $ws -P 2   -S 191$run -s brtest_${ws}_${run}  "
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh

		LAPARAMS=" -T 180 -J -B $ws -P 2   -S 191$run -s brtestL_${ws}_${run}  -L"
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh

		LAPARAMS=" -T 180 -J -B $ws -P 2   -S 191$run -s brtestG_${ws}_${run}  -G"
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh

		#LAPARAMS=" -P 2  -T $ws -S 1980$run -s "strong2N_${ws}_${run}"   -n"
	#	qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
	#	sleep 1
	#	LAPARAMS=" -P 2  -T $ws -S 1980$run -s "strong2NL_${ws}_${run}"  -n -L "
	#	qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
	#	sleep 1
		#LAPARAMS=" -P 2  -T $ws -S 1980$run -s "strong2Alt_${ws}_${run}" -C "
		#qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		#LAPARAMS="-P 2  -T $ws -S 1980$run -s "strong2LAlt_${ws}_${run}"  -C -L"
		#qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		#LAPARAMS=" -P 2  -T $ws -S 1980$run -s "strong2NAlt_${ws}_${run}"  -C -n"
		#qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		#LAPARAMS=" -P 2  -T $ws -S 1980$run -s "strong2NLAlt_${ws}_${run}" -C -n -L "
		#qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		#sleep 1

	done
done



