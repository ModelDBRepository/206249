function la_run {
	qsub -v "LAPARAMS=$LAPARAMS" submit_lamodel.sh
	#echo ./lamodel $LAPARAMS 
	#./lamodel $LAPARAMS &
}



for run in 0 1 2 3 4 5 6 7 8 9
do
	for ParamVal in 0.7 0.8 0.9 1.0 1.1 1.2 1.3
	do

		for ParamName in CREBTimeParam connectivityParam inhibitionParam initWeight maxWeight dendSpikeThresh globalPRPThresh localPRPThresh homeostasisTimeParam nNeuronsParam nBranchesParam
		do

			#LAPARAMS=" -P 1 -T 180 -S 1980$run -s NPERT_${ParamName}_${ParamVal}_${run} -o connectivityParam=2.363 -o dendSpikeThresh=3.0  -o ${ParamName}=${ParamVal}"

			LAPARAMS=" -P 1 -T 180 -S 1980$run -s PERT_${ParamName}_${ParamVal}_${run}   -o ${ParamName}=${ParamVal}"
			la_run

			LAPARAMS=" -P 1 -T 180 -S 1980$run -s PERTL_${ParamName}_${ParamVal}_${run}  -o ${ParamName}=${ParamVal} -L"
			la_run

			LAPARAMS=" -P 1 -T 180 -S 1980$run -s PERTG_${ParamName}_${ParamVal}_${run}  -o ${ParamName}=${ParamVal} -G"
			la_run

			LAPARAMS=" -w 2 -P 2 -T 60 -S 1980$run -s PERTWS_${ParamName}_${ParamVal}_${run}  -o ${ParamName}=${ParamVal} "
			la_run

			LAPARAMS=" -w 2 -P 2 -T 60 -S 1980$run -s PERTWSG_${ParamName}_${ParamVal}_${run}   -o ${ParamName}=${ParamVal} -G"
			la_run

			LAPARAMS=" -w 2 -P 2 -T 60 -S 1980$run -s PERTWSL_${ParamName}_${ParamVal}_${run}  -o ${ParamName}=${ParamVal} -L"
			la_run
		done
	done
done

