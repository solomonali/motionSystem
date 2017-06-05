#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <mraa/i2c.h>
#include <sys/time.h>
#include "LSM9DS0.h"
#include <sys/file.h>

#define MILLION 1000000.0f

sig_atomic_t volatile run_flag = 1;

void do_when_interrupted(int sig) 
{
	if (sig == SIGINT)
		run_flag = 0;
}

double parse_tv(struct timeval *tv)
{
	return (double) tv->tv_sec + (double) tv->tv_usec/MILLION;
}

int main() {
	mraa_i2c_context 	accel, gyro;
	accel_scale_t 		a_scale = A_SCALE_8G;
	gyro_scale_t 		g_scale = G_SCALE_500DPS;
	
	data_t 			accel_data, gyro_data;
	float 			a_res, g_res;

	char			file_name[10];
	FILE			*fp;
	struct timeval 		tv_before, tv_after,tv_start;


	
/*	int cntr = 0;
	char number[10];
	char file[10] = "file";
	char extension[10]=".csv";
*/	char *file[]= {"file0.csv","file1.csv","file2.csv","file3.csv","file4.csv","file5.csv",
		       "file6.csv","file7.csv","file8.csv","file9.csv","file10.csv","file11.csv",
		       "file12.csv","file13.csv","file14.csv"};
	double startTime = 0;
	int fd;
	int success_flag;
	int motion=0;
	char *motionArr[]={"Walking:speed 1","Walking:speed 2","Walking:speed 3"
	,"Running:speed 1","Running:speed 2",
	"Jumping: level 1","Jumping:level 2",
	"90 degree left turn","90 degree right turn",
	"Stair Ascent:speed 1","Stair Ascent:speed 2","Stair Ascent:speed 3",
	"Stair Descent:speed 1","Stair Descent:speed 2","Stair Descent:speed 3",
	};
				

	while(motion < 15)
	{
		signal(SIGINT, &do_when_interrupted);
		if(!run_flag){break;}
//		if(cntr==20){cntr=0;}
//		sprintf(number,"%d",cntr);
//		strcpy(file_name,file);
//		strcat(file_name,number);
//		strcat(file_name,extension);
		fp = fopen(file[motion], "w");
		if (fp == NULL) {
			fprintf(stderr, "Failed to open file \'%s\'. Exiting.\n",
					file_name);
		exit(EXIT_FAILURE);
		}	
//		cntr++;
//		fd = fileno(fp);

//		flock(fd,LOCK_EX);
		//initialize sensors, set scale, and calculate resolution.
		accel = accel_init();
		set_accel_scale(accel, a_scale);	
		a_res = calc_accel_res(a_scale);

		gyro = gyro_init();
		set_gyro_scale(gyro, g_scale);
		g_res = calc_gyro_res(g_scale);

/*		fprintf(fp, "%s,%s,%s,%s,%s,%s,%s,%s\n",
				"timestamp_before",
				"timestamp_after",
				"accel_x", "accel_y", "accel_z",
				"gyro_x", "gyro_y", "gryo_z"
		       );
*/
		printf("Please do %s\n", motionArr[motion]);
		printf("please press the return key when you are ready...\n");

				

		do {
			success_flag = getchar();
		} while (success_flag != '\n');
		
		printf("go\n");


		//Read the sensor data and print them.
		gettimeofday(&tv_start,NULL);
		startTime = tv_start.tv_sec;
		while (((tv_after.tv_sec)-startTime)<6){
			if(!run_flag){break;}
			gettimeofday(&tv_before, NULL);
			accel_data = read_accel(accel, a_res);
			gyro_data = read_gyro(gyro, g_res);
			gettimeofday(&tv_after, NULL);
			fprintf(fp, "%20.10lf,%20.10lf,",
				       	parse_tv(&tv_before),
					parse_tv(&tv_after)
			       );
			fprintf(fp, "%8.4lf,%8.4lf,%8.4lf,",
					accel_data.x,
					accel_data.y,
					accel_data.z
			       );
			fprintf(fp, "%8.4lf,%8.4lf,%8.4lf\n",
					gyro_data.x,
					gyro_data.y,
					gyro_data.z
			      );
		//	usleep(100);
		}	
		fclose(fp);
		motion++;
		printf("stop\n\n");
	
	}
	printf("Done!\n");

	return 0;	
}

