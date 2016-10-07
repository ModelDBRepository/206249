
for s in  0 1 2 3 4 ; do
(
	for j in 1 2 3 4 5 6 ; do



		LAPARAMS=" -P $j -T 1440 -S 1980$s  -s repeated_${j}_${s} -R  "
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		echo $LAPARAMS

		LAPARAMS=" -P $j -T 1440 -S 1980$s  -s repeatedL_${j}_${s} -R  -L  "
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		echo $LAPARAMS

		LAPARAMS=" -P $j -T 1440 -S 1980$s  -s repeatedG_${j}_${s} -R  -G  "
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		echo $LAPARAMS


		LAPARAMS=" -P $j -T 1440 -S 1980$s  -s repeatedN_${j}_${s} -R  -n "
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		echo $LAPARAMS

		LAPARAMS=" -P $j -T 1440 -S 1980$s  -s repeatedNL_${j}_${s} -R  -L  -n "
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		echo $LAPARAMS

		LAPARAMS=" -P $j -T 1440 -S 1980$s  -s repeatedNG_${j}_${s} -R -G -n  "
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		echo $LAPARAMS



		LAPARAMS=" -P $j -T 1440 -S 1980$s  -s repeatedU_${j}_${s} -R -U "
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		echo $LAPARAMS

		LAPARAMS=" -P $j -T 1440 -S 1980$s  -s repeatedUL_${j}_${s} -R -U -L  "
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		echo $LAPARAMS

		LAPARAMS=" -P $j -T 1440 -S 1980$s  -s repeatedUG_${j}_${s} -R -U -G  "
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		echo $LAPARAMS


		LAPARAMS=" -P $j -T 1440 -S 1980$s  -s repeatedUN_${j}_${s} -R -U -n "
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		echo $LAPARAMS

		LAPARAMS=" -P $j -T 1440 -S 1980$s  -s repeatedUNL_${j}_${s} -R -U -L  -n "
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		echo $LAPARAMS

		LAPARAMS=" -P $j -T 1440 -S 1980$s  -s repeatedUNG_${j}_${s} -R -U -G -n  "
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		echo $LAPARAMS

		#LAPARAMS=" -P $j -T 1440 -S 1980$s  -s repeatedUC_${j}_${s}  "
	done
) 
done

