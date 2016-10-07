

function la_run {
	#qsub -v "LAPARAMS=$LAPARAMS" submit_lamodel.sh
	./lamodel $LAPARAMS  
	#sleep 1
}


#altparams="-o connectivityParam=2.363 -o dendSpikeThresh=3.0 -o homeostasisTimeParam=3.0"
altparams=""

for ws in 60 120 300 1440 2880 ;  do
	for run in  {0..10}; do
		
		LAPARAMS=" -P 2  -T $ws -S 191$run -s strong2_${ws}_${run}  ${altparams}"
		#la_run 

		LAPARAMS=" -P 2  -T $ws -S 191$run -s strong2G_${ws}_${run}  -G   ${altparams}"
		#la_run 

		LAPARAMS=" -P 2  -T $ws -S 191$run -s strong2L_${ws}_${run}  -L   ${altparams}"
		#la_run 



		LAPARAMS=" -P 2  -T $ws -S 191$run -s strong2N_${ws}_${run}  -n ${altparams} "
		la_run 

		LAPARAMS=" -P 2  -T $ws -S 191$run -s strong2NL_${ws}_${run} -L  -n  ${altparams}"
		la_run 

		LAPARAMS=" -P 2  -T $ws -S 191$run -s strong2NG_${ws}_${run} -G  -n ${altparams} "
		la_run 

	done
done



