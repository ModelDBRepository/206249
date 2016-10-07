function la_run {
	#qsub -v "LAPARAMS=$LAPARAMS" submit_lamodel.sh
	./lamodel $LAPARAMS 
}



for s in 0 1 2 3 4 5 6 7 8 9; do
	for i in 60 300 1440; do

		#LAPARAMS="-P 10 -T $i -S 19$s  -s multi_${i}_${s}"
		#qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh

		LAPARAMS="-P 10   -T $i -S 19$s  -L   -s multiL_${i}_${s}"
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh

		LAPARAMS="-P 10   -T $i -S 19$s  -G   -s multiG_${i}_${s}"
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh

		LAPARAMS="-P 10  -T $i -S 19$s  -G  -n  -s multiGN_${i}_${s}"
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh

		LAPARAMS="-P 10 -T $i -S 19$s  -n -L  -s multiLN_${i}_${s}" 
		qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh

		#LAPARAMS=" -U -T $i -S 19$s -H 1 -n  -s multiUHN_${i}_${s}"
		#qsub -v "LAPARAMS=${LAPARAMS}" submit_lamodel.sh
		#sleep 1


	done
done

