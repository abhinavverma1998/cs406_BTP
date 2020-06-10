make clean
if [ "$3" == "fixed" ]; then
	cat fixed > GlobalTypes.h
elif [ "$3" == "var" ]; then
	cat var > GlobalTypes.h
fi
make GCC=$1 OPT=$4
./MCM input_medium.txt $2