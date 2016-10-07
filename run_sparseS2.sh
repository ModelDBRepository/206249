function la_run {
	qsub -v "LAPARAMS=$LAPARAMS" submit_lamodel.sh
	#./lamodel $LAPARAMS &
}


for run in {0..9}  ; do

	LAPARAMS=" -P 1 -T 180 -J -S 1980$run -s sparseS2_${run} "
	la_run

	LAPARAMS=" -P 1 -T 180 -J -S 1980$run -s sparseS2L_${run}  -L "
	la_run

	LAPARAMS=" -P 1 -T 180 -J -S 1980$run -s sparseS2G_${run}  -G "
	la_run


	#LAPARAMS=" -P 1 -T 180 -n -J -S 1980$run -s sparseS2N_${run} "
	#la_run

	#LAPARAMS=" -P 1 -T 180 -n -J -S 1980$run -s sparseNL_${run}  -L "
	#la_run
done


