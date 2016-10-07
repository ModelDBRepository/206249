function la_run {
	qsub -v "LAPARAMS=$LAPARAMS  -o dendSpikeThresh=3.0  -o connectivityParam=2.363" submit_lamodel.sh
	#./lamodel $LAPARAMS &
}


for run in {0..9}  ; do

	LAPARAMS=" -P 1 -T 180 -J -S 1980$run -s Nsparse_${run} "
	la_run

	LAPARAMS=" -P 1 -T 180 -J -S 1980$run -s NsparseL_${run}  -L "
	la_run

	LAPARAMS=" -P 1 -T 180 -J -S 1980$run -s NsparseG_${run}  -G "
	la_run


	LAPARAMS=" -P 1 -T 180 -n -J -S 1980$run -s NsparseN_${run} "
	la_run

	#LAPARAMS=" -P 1 -T 180 -n -J -S 1980$run -s sparseNL_${run}  -L "
	#la_run
done


