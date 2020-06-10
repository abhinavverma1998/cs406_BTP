gcc=( g++-4.8 )
input=( 1025 )
type=( fixed var )
opt_flags=( O0 O1 O2 O3 Os Ofast)
echo "" > log_2

for i in "${gcc[@]}"
do
	for j in "${input[@]}"
	do
		for l in "${opt_flags[@]}"
		do	
			for k in "${type[@]}"
			do  
				echo "GCC version : " $i " Input size : " $j " Type : " $k >> log_2
				./runme.sh $i $j $k $l>> log_2
				echo "" >> log_2
			done
		done
	done
done
