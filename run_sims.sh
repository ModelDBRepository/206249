for ITV in 60 120 180 300
do
	for i in 1 2 3 4 5 6 
	do
		sh job_sims.sh $ITV $i sn "" &
	done
done

