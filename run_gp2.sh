function la_run {
	qsub -v "LAPARAMS=$LAPARAMS" submit_lamodel.sh
	#echo  not running ./lamodel $LAPARAMS 
	#./lamodel $LAPARAMS &
	
}

function ws_run {
	qsub -v "LAPARAMS=$LAPARAMS" submit_lamodel.sh
	#echo ./lamodel $LAPARAMS 
	#./lamodel $LAPARAMS &
}




for  CREBTimeParam in 0.5 1.0 1.5; do
	for  connectivityParam in 0.5 1.0 1.5; do
		for  globalPRPThresh in 0.5 1.0 1.5; do
			for  localPRPThresh in 0.5 1.0 1.5; do
				for  dendSpikeThresh in 0.5 1.0 1.5; do
					for  initWeight in 0.5 1.0 1.5; do
for run in 0 1 2 3 4;
do
			LAPARAMS=" -P 1 -T 180 -S 1980$run -s SENS_${CREBTimeParam}_${connectivityParam}_${globalPRPThresh}_${localPRPThresh}_${dendSpikeThresh}_${initWeight}_${run}  -o CREBTimeParam=${CREBTimeParam} -o connectivityParam=${connectivityParam} -o globalPRPThresh=${globalPRPThresh} -o localPRPThresh=${localPRPThresh} -o dendSpikeThresh=${dendSpikeThresh} -o initWeight=${initWeight}"
			la_run
			LAPARAMS=" -P 2 -T 60 -w 2 -S 1980$run -s SENSWS_${CREBTimeParam}_${connectivityParam}_${globalPRPThresh}_${localPRPThresh}_${dendSpikeThresh}_${initWeight}_${run}  -o CREBTimeParam=${CREBTimeParam} -o connectivityParam=${connectivityParam} -o globalPRPThresh=${globalPRPThresh} -o localPRPThresh=${localPRPThresh} -o dendSpikeThresh=${dendSpikeThresh} -o initWeight=${initWeight}"
			ws_run 
done

					done
				done
			done
		done
	done
done

