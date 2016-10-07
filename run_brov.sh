function la_run {
	qsub -v "LAPARAMS=$LAPARAMS" submit_lamodel.sh
	#echo ./lamodel $LAPARAMS 
	#./lamodel $LAPARAMS &
}

for BROV in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do

	for run in 0 1 2 3 4 5 6 7 8 9 ; do
		LAPARAMS=" -P 1 -O $BROV -T 180 -J -S 1980$run -s brov${BROV}_${run} "
		la_run

		LAPARAMS=" -P 1 -O $BROV -T 180 -J -S 1980$run -s brovL${BROV}_${run}  -L "
		la_run

		LAPARAMS=" -P 1 -O $BROV -T 180 -J -S 1980$run -s brovG${BROV}_${run}  -G "
		la_run




		#sleep 1
		#LAPARAMS=" -P 1 -T 180 -J -S 1980$run -s sparseAlt_${run} -C"
		#la_run
		#
		#sleep 1
		#LAPARAMS=" -P 1 -T 180 -J -S 1980$run -s sparseLAlt_${run}  -L -C"
		# la_run
		#
		#sleep 1
	done
done


