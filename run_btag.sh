function la_run {
	#qsub -v "LAPARAMS=$LAPARAMS" submit_lamodel.sh
	#sleep 1
	./lamodel $LAPARAMS  
	#sleep 1
}




COND=weakstrong

for ws in 60 120 180 300 1440; do
	for run in {0..9}; do
	

		LAPARAMS="-P 2  -w 1 -n  -T $ws -S 1980$run -s ${COND}N_0_${ws}_${run}  "
		la_run

		LAPARAMS="-P 2 -w 2   -n  -T $ws -S 1980$run -s ${COND}N_${ws}_0_${run}  "
		la_run

		LAPARAMS="-P 2 -w 1  -n  -T $ws -S 1980$run -s ${COND}NL_0_${ws}_${run}  -L  "
		la_run

		LAPARAMS="-P 2 -w 2  -n   -T $ws -S 1980$run -s ${COND}NL_${ws}_0_${run}  -L   "
		la_run

		LAPARAMS="-P 2 -w 1  -n  -T $ws -S 1980$run -s ${COND}NG_0_${ws}_${run}  -G  "
		la_run

		LAPARAMS="-P 2 -w 2  -n   -T $ws -S 1980$run -s ${COND}NG_${ws}_0_${run}  -G   "
		la_run




		LAPARAMS="-P 2  -w 1   -T $ws -S 1980$run -s ${COND}_0_${ws}_${run}  "
		la_run

		LAPARAMS="-P 2 -w 2     -T $ws -S 1980$run -s ${COND}_${ws}_0_${run}  "
		la_run

		LAPARAMS="-P 2 -w 1    -T $ws -S 1980$run -s ${COND}L_0_${ws}_${run}  -L  "
		la_run

		LAPARAMS="-P 2 -w 2     -T $ws -S 1980$run -s ${COND}L_${ws}_0_${run}  -L   "
		la_run

		LAPARAMS="-P 2 -w 1    -T $ws -S 1980$run -s ${COND}G_0_${ws}_${run}  -G  "
		la_run

		LAPARAMS="-P 2 -w 2     -T $ws -S 1980$run -s ${COND}G_${ws}_0_${run}  -G   "
		la_run




	done
done


