#!/bin/bash

#ifile=$1
#activity=$2
#if [ $3 ]; then
#		threshold=$3	
#else
#		threshold=100
#fi
#if  false 
#then	
echo ""
echo ">>>>COMPILING<<<<"
make clean
make

./imu_data

./features_ext_walking_speeds file0.csv pt_output.csv strides.csv 80 w1
./features_ext_walking_speeds file1.csv pt_output.csv strides.csv 80 w2       	
./features_ext_walking_speeds file2.csv pt_output.csv strides.csv 80 w3
nLines=$(wc -l < ws_temp.txt)
nSamples=$((($nLines)/2))
echo "$nSamples 2 3" > ws_train_set.txt
cat ws_temp.txt >> ws_train_set.txt

./features_ext_running_speeds file3.csv pt_output.csv strides.csv 80 w1
./features_ext_running_speeds file4.csv pt_output.csv strides.csv 80 w2       	
nLines=$(wc -l < rs_temp.txt)
nSamples=$((($nLines)/2))
echo "$nSamples 2 2" > rs_train_set.txt
cat rs_temp.txt >> rs_train_set.txt

./features_ext_jumping_levels file5.csv pt_output.csv strides.csv 80 w1
./features_ext_jumping_levels file6.csv pt_output.csv strides.csv 80 w2       	
nLines=$(wc -l < jl_temp.txt)
nSamples=$((($nLines)/2))
echo "$nSamples 2 2" > jl_train_set.txt
cat jl_temp.txt >> jl_train_set.txt

./features_ext_stairA_speeds file9.csv pt_output.csv strides.csv 80 w1
./features_ext_stairA_speeds file10.csv pt_output.csv strides.csv 80 w2       	
./features_ext_stairA_speeds file11.csv pt_output.csv strides.csv 80 w3
nLines=$(wc -l < sas_temp.txt)
nSamples=$((($nLines)/2))
echo "$nSamples 2 3" > sas_train_set.txt
cat sas_temp.txt >> sas_train_set.txt

./features_ext_stairD_speeds file12.csv pt_output.csv strides.csv 80 w1
./features_ext_stairD_speeds file13.csv pt_output.csv strides.csv 80 w2       	
./features_ext_stairD_speeds file14.csv pt_output.csv strides.csv 80 w3
nLines=$(wc -l < sds_temp.txt)
nSamples=$((($nLines)/2))
echo "$nSamples 2 3" > sds_train_set.txt
cat sds_temp.txt >> sds_train_set.txt

####
cat file0.csv > walking_data.csv
cat file1.csv >> walking_data.csv
cat file2.csv >> walking_data.csv

cat file3.csv > running_data.csv
cat file4.csv >> running_data.csv

cat file5.csv > jumping_data.csv
cat file6.csv >> jumping_data.csv

cat file7.csv > leftTurn_data.csv
cat file8.csv > rightTurn_data.csv


cat file9.csv > stairA_data.csv
cat file10.csv >> stairA_data.csv
cat file11.csv >> stairA_data.csv

cat file12.csv > stairD_data.csv
cat file13.csv >> stairD_data.csv
cat file14.csv >> stairD_data.csv

rm file*
cat running_data.csv > others_walking.csv
cat jumping_data.csv >> others_walking.csv
cat leftTurn_data.csv >> others_walking.csv
cat rightTurn_data.csv >> others_walking.csv
cat stairA_data.csv >> others_walking.csv
cat stairD_data.csv >> others_walking.csv

cat walking_data.csv > others_running.csv
cat jumping_data.csv >> others_running.csv
cat leftTurn_data.csv >> others_running.csv
cat rightTurn_data.csv >> others_running.csv
cat stairA_data.csv >> others_running.csv
cat stairD_data.csv >> others_running.csv

cat walking_data.csv > others_jumping.csv
cat running_data.csv >> others_jumping.csv
cat leftTurn_data.csv >> others_jumping.csv
cat rightTurn_data.csv >> others_jumping.csv
cat stairA_data.csv >> others_jumping.csv
cat stairD_data.csv >> others_jumping.csv

cat walking_data.csv > others_turnning.csv
cat jumping_data.csv >> others_turning.csv
cat running_data.csv >> others_turning.csv
cat stairA_data.csv >> others_turning.csv
cat stairD_data.csv >> others_turning.csv


cat walking_data.csv > others_stairA.csv
cat jumping_data.csv >> others_stairA.csv
cat leftTurn_data.csv >> others_stairA.csv
cat rightTurn_data.csv >> others_stairA.csv
cat running_data.csv >> others_stairA.csv
cat stairD_data.csv >> others_stairA.csv

cat walking_data.csv > others_stairD.csv
cat jumping_data.csv >> others_stairD.csv
cat leftTurn_data.csv >> others_stairD.csv
cat rightTurn_data.csv >> others_stairD.csv
cat running_data.csv >> others_stairD.csv
cat stairA_data.csv >> others_stairD.csv


#####
./features_ext_turning leftTurn_data.csv pt_output.csv strides.csv 80 w1
./features_ext_turning rightTurn_data.csv pt_output.csv strides.csv 80 w2 
./features_ext_turning others_turning.csv pt_output.csv strides.csv 80 w3  
nLines=$(wc -l < t_temp.txt)
nSamples=$((($nLines)/2))
echo "$nSamples 1 3" > t_train_set.txt
cat t_temp.txt >> t_train_set.txt


./features_ext_walking walking_data.csv pt_output.csv strides.csv 80 w1
./features_ext_walking others_walking.csv pt_output.csv strides.csv 80 w2
#./features_ext_walking running_data.csv pt_output.csv strides.csv 80 w2
#./features_ext_walking jumping_data.csv pt_output.csv strides.csv 80 w3
#./features_ext_walking leftTurn_data.csv pt_output.csv strides.csv 80 t1
#./features_ext_walking rightTurn_data.csv pt_output.csv strides.csv 80 t2
#./features_ext_walking stairA_data.csv pt_output.csv strides.csv 80 w3
#./features_ext_walking stairD_data.csv pt_output.csv strides.csv 80 w3

nLines=$(wc -l < w_temp.txt)
nSamples=$((($nLines)/2))
echo "$nSamples 3 2" > w_train_set.txt
cat w_temp.txt >> w_train_set.txt


./features_ext_running running_data.csv pt_output.csv strides.csv 80 w1
./features_ext_running others_running.csv pt_output.csv strides.csv 80 w2
nLines=$(wc -l < r_temp.txt)
nSamples=$((($nLines)/2))
echo "$nSamples 3 2" > r_train_set.txt
cat r_temp.txt >> r_train_set.txt

./features_ext_jumping jumping_data.csv pt_output.csv strides.csv 80 w1
./features_ext_jumping others_jumping.csv pt_output.csv strides.csv 80 w2      
nLines=$(wc -l < j_temp.txt)
nSamples=$((($nLines)/2))
echo "$nSamples 4 2" > j_train_set.txt
cat j_temp.txt >> j_train_set.txt

./features_ext_stairA stairA_data.csv pt_output.csv strides.csv 80 w1
./features_ext_stairA others_stairA.csv pt_output.csv strides.csv 80 w2       	
nLines=$(wc -l < sa_temp.txt)
nSamples=$((($nLines)/2))
echo "$nSamples 4 2" > sa_train_set.txt
cat sa_temp.txt >> sa_train_set.txt

./features_ext_stairD stairD_data.csv pt_output.csv strides.csv 80 w1
./features_ext_stairD others_stairD.csv pt_output.csv strides.csv 80 w2       	
nLines=$(wc -l < sd_temp.txt)
nSamples=$((($nLines)/2))
echo "$nSamples 4 2" > sd_train_set.txt
cat sd_temp.txt >> sd_train_set.txt
#fi
./w_train_neural_net
./ws_train_neural_net
./r_train_neural_net
./rs_train_neural_net
./j_train_neural_net
./jl_train_neural_net
./t_train_neural_net
./sa_train_neural_net
./sas_train_neural_net
./sd_train_neural_net
./sds_train_neural_net








#echo ""
#echo "Number of lines in train_set.txt"
#wc -l train_set.txt




