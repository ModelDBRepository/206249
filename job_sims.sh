for i in 0 1 2 3 4 5 6 7 8 9
do
	#echo ./lamodel -S 198$1 $*
	./lamodel -S 198$i -P 2 -p 6 -i 1 -T $1  -d $2 -s "_sim$2"
	./lamodel -S 198$i -P 2 -p 6 -i 1 -T $1 -c  -d $2 -s "_sim$2"
	python stats.py N400.B20.I12.i1.P2.p6.T$1.S198$i.sn_sim$2 NO 0
	python stats.py N400.B20.I12.i1.P2.p6.T$1.S198$i.sc_sim$2 NO 0
	echo Done
done


