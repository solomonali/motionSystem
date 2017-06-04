#!/bin/bash

echo -n "Is a new user? "
read text
if [ $text == "yes" ]
then
	rm *.net
fi
if [ ! -e ./*.net ]
then
	(cd ./in_field_training; sh generate_extraction.sh)
fi

echo ""
echo ">>>>COMPILING<<<<"

rm -f file*
make clean
make


echo ""
echo ">>>>EXECUTING imu_data<<<<"
./imu_data &
imu_process=$!

echo ""
echo ">>>>executing extract_stride_data<<<<"
./extract_features_data &
extract_process=$!	

#echo ""
#echo ">>>>executing extract_stride_data<<<<"
#./test_neural_network &
 
read -p "press anything to continue..." -n1 -s

echo "-------------"
kill $imu_process
kill $extract_process
