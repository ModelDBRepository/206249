function la_run {
	#qsub -v "LAPARAMS=$LAPARAMS" submit_lamodel.sh
	./lamodel $LAPARAMS 
}


for run in {0..9}  ; do

	LAPARAMS=" -P 1 -T 180 -J -S 1980$run -s sparse_${run} "
	la_run

	LAPARAMS=" -P 1 -T 180 -J -S 1980$run -s sparseL_${run}  -L "
	la_run

	LAPARAMS=" -P 1 -T 180 -J -S 1980$run -s sparseG_${run}  -G "
	la_run


	LAPARAMS=" -P 1 -T 180 -n -J -S 1980$run -s sparseN_${run} "
	la_run

	#LAPARAMS=" -P 1 -T 180 -n -J -S 1980$run -s sparseNL_${run}  -L "
	#la_run
done


