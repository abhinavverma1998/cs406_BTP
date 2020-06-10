depth=( 10000 50000 80000 100000 130000 )
for i in "${depth[@]}"
do
	./runme.sh $i >> log_iterative 
	echo "" >> log_iterative
done
