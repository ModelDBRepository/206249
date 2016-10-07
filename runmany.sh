
for i in 0 1 2 3 4 5 6 7 8 9
do
	./lamodel -S 198$i $* &
	echo Done
done


