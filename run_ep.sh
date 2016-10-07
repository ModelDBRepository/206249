
for s in  80 81 82 83 84 85 86 87 88 89; do
	(
	for i in 60 180; do
		#./lamodel -T $i -S 19$s 
		#./lamodel -T $i -S 19$s  -L 
		#./lamodel -T $i -S 19$s  -n
		./lamodel -T $i -S 19$s  -L -n
	done
	) &
done

