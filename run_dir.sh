

function la_run {
	qsub -v "LAPARAMS=$LAPARAMS" submit_lamodel.sh
	#./lamodel $LAPARAMS  
	#sleep 1
}


for ws in 60 120 300 1440 2880 ;  do
	for run in  {0..5}; do
		
		LAPARAMS=" -o protocol=2  -T ${ws} -S 191$run -s dir1_${ws}_${run}  "
		la_run 
		LAPARAMS=" -o protocol=3  -T ${ws} -S 191$run -s dir2_${ws}_${run}  "
		la_run 

		LAPARAMS=" -o protocol=2  -T ${ws} -S 191$run -s dir1L_${ws}_${run} -L "
		la_run 
		LAPARAMS=" -o protocol=3  -T ${ws} -S 191$run -s dir2L_${ws}_${run} -L "
		la_run 

		LAPARAMS=" -o protocol=2  -T ${ws} -S 191$run -s dir1G_${ws}_${run} -G "
		la_run 
		LAPARAMS=" -o protocol=3  -T ${ws} -S 191$run -s dir2G_${ws}_${run} -G "
		la_run 


	done
done



