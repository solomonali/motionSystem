CC=gcc
LDFLAGS=-lmraa -lfann
CFLAGS=-Wall
#SOURCES= test_neural_network.c train_neural_net.c
EXECUTABLES=$(SOURCES:.c=)


all: imu_data features_ext_walking features_ext_walking_speeds features_ext_running features_ext_running_speeds features_ext_jumping features_ext_jumping_levels features_ext_turning features_ext_stairA features_ext_stairA_speeds features_ext_stairD features_ext_stairD_speeds w_train_neural_net ws_train_neural_net r_train_neural_net rs_train_neural_net j_train_neural_net jl_train_neural_net t_train_neural_net sa_train_neural_net sas_train_neural_net sd_train_neural_net sds_train_neural_net


imu_data: imu_data.c
	gcc -lmraa -lm -o imu_data imu_data.c LSM9DS0.c

features_ext_walking: features_ext_walking.c
	gcc -o features_ext_walking features_ext_walking.c -lm
features_ext_walking_speeds: features_ext_walking_speeds.c
	gcc -o features_ext_walking_speeds features_ext_walking_speeds.c -lm

features_ext_running: features_ext_running.c
	gcc -o features_ext_running features_ext_running.c -lm

features_ext_running_speeds: features_ext_running_speeds.c
	gcc -o features_ext_running_speeds features_ext_running_speeds.c -lm

features_ext_jumping: features_ext_jumping.c
	gcc -o features_ext_jumping features_ext_jumping.c -lm

features_ext_jumping_levels: features_ext_jumping_levels.c
	gcc -o features_ext_jumping_levels features_ext_jumping_levels.c -lm

features_ext_turning: features_ext_turning.c
	gcc -o features_ext_turning features_ext_turning.c -lm

features_ext_stairA: features_ext_stairA.c
	gcc -o features_ext_stairA features_ext_stairA.c -lm

features_ext_stairA_speeds: features_ext_stairA_speeds.c
	gcc -o features_ext_stairA_speeds features_ext_stairA_speeds.c -lm

features_ext_stairD: features_ext_stairD.c
	gcc -o features_ext_stairD features_ext_stairD.c -lm

features_ext_stairD_speeds: features_ext_stairD_speeds.c
	gcc -o features_ext_stairD_speeds features_ext_stairD_speeds.c -lm

w_train_neural_net: ../stage3/w_train_neural_net.c
	$(CC) $(CFLAGS) -o w_train_neural_net ../stage3/w_train_neural_net.c $(LDFLAGS)

ws_train_neural_net: ../stage3/ws_train_neural_net.c
	$(CC) $(CFLAGS) -o ws_train_neural_net ../stage3/ws_train_neural_net.c $(LDFLAGS)

r_train_neural_net: ../stage3/r_train_neural_net.c
	$(CC) $(CFLAGS) -o r_train_neural_net ../stage3/r_train_neural_net.c $(LDFLAGS)

rs_train_neural_net: ../stage3/rs_train_neural_net.c
	$(CC) $(CFLAGS) -o rs_train_neural_net ../stage3/rs_train_neural_net.c $(LDFLAGS)

j_train_neural_net: ../stage3/j_train_neural_net.c
	$(CC) $(CFLAGS) -o j_train_neural_net ../stage3/j_train_neural_net.c $(LDFLAGS)

jl_train_neural_net: ../stage3/jl_train_neural_net.c
	$(CC) $(CFLAGS) -o jl_train_neural_net ../stage3/jl_train_neural_net.c $(LDFLAGS)

t_train_neural_net: ../stage3/t_train_neural_net.c
	$(CC) $(CFLAGS) -o t_train_neural_net ../stage3/t_train_neural_net.c $(LDFLAGS)

sa_train_neural_net: ../stage3/sa_train_neural_net.c
	$(CC) $(CFLAGS) -o sa_train_neural_net ../stage3/sa_train_neural_net.c $(LDFLAGS)

sas_train_neural_net: ../stage3/sas_train_neural_net.c
	$(CC) $(CFLAGS) -o sas_train_neural_net ../stage3/sas_train_neural_net.c $(LDFLAGS)

sd_train_neural_net: ../stage3/sd_train_neural_net.c
	$(CC) $(CFLAGS) -o sd_train_neural_net ../stage3/sd_train_neural_net.c $(LDFLAGS)

sds_train_neural_net: ../stage3/sds_train_neural_net.c
	$(CC) $(CFLAGS) -o sds_train_neural_net ../stage3/sds_train_neural_net.c $(LDFLAGS)


clean:
	rm -rf imu_data *~
	rm -rf features_ext_walking *~
	rm -rf features_ext_walking_speeds *~
	rm -rf features_ext_running *~
	rm -rf features_ext_running_speeds *~
	rm -rf features_ext_jumping *~
	rm -rf features_ext_jumping_speeds *~
	rm -rf features_ext_turning *~
	rm -rf features_ext_stairA *~
	rm -rf features_ext_stairA_speeds *~
	rm -rf features_ext_stairD *~
	rm -rf features_ext_stairA_speeds *~
	rm -f train_neural_net
	rm -f *~
	rm -f ../*.net
	rm *.txt *.csv


