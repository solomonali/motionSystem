/* for file and terminal I/O */
#include <stdio.h>
/* for string manip */
#include <string.h>
/* for exit() */
#include <stdlib.h>
/* for fabsf() */
#include <math.h>

#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
#include <sys/file.h>
#include "floatfann.h"
#include <string.h>

#define BUFF_SIZE 1024

void clear_buffer(float *arr, float val, int n); 
float calculate_mean(float *arr, int start, int end);
void calculate_Max_Min_Range(float *arr,int start,int end,
		float *max, float *min, float *range);
void calculate_Statistics (float *arr, int start, int end, float mean, 
		float *MAD, float *variance, float *std,
		float *skewness, float *kurtosis);
int find_peaks_and_troughs(float *arr,int n_samples,float E,float *P, float *T,
		 	int *n_P, int *n_T);
int stride_extraction(int ind_thr,int n_P, float *P_i,float *T_i,
		 float *S_imax, float *S_imin);
void featureExtraction(int n_S,float *S_imax,float **features,float *axis);

sig_atomic_t volatile run_flag = 1;

void do_when_interrupted(int sig)
{
	if(sig == SIGINT)
		run_flag = 0;
}

int main()
{
	/* Generic variables */
	int i,ind_thr,rv;
	/* Variables for reading file line by line */
//	char ifile_name[100];
	FILE *fp;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int N_SAMPLES;


	/* Variables for storing the data and storing the return values */
	double *t,*t_b, *t_a; 
	/* variables for data collected from input file */
	float *x_ac, *y_ac, *z_ac, *x_gy, *y_gy, *z_gy; 	
	float pk_threshold = 100.00;	// pk-threshold value
       	/* Variables for peak-trough detection */	
	float *P_i; 	// indicies of each peak found by peak detection
	float *T_i; 	// indicies of each trough found by trough detection
	float *S_imax,*S_imin; 	// indicies of the start of each stride
	int n_P; 	// number of peaks
	int n_T; 	// number of troughs
	int n_S = 1000; 	// number of strides

//	char number[10];
	int icntr = 0;
	char *ifile_name[] = {"file0.csv","file1.csv","file2.csv","file3.csv",
		"file4.csv","file5.csv","file6.csv","file7.csv","file8.csv",
		"file9.csv","file10.csv"};
//	char iextension[10] = ".csv";
	int fd;

	/* general variables */   
       	int ii,j,k,start,end;
  	float max, offset=0;


    	fann_type *calc_out;
    	fann_type input[10];

	/* feature extraction variables */ 
	float mean_zgy[n_S],kurtosis_zgy[n_S],holder[n_S],skewness_zgy[n_S],max_zgy[n_S],
	      mean_xac[n_S],min_xac[n_S],variance_xac[n_S],mean_xac3[n_S],min_xac4[n_S],
	      mean_xgy[n_S],max_xgy[n_S],
	      ratio_yac[n_S],max_yac[n_S],min_yac[n_S],mean_yac[n_S],
	      skewness_yac[n_S],variance_yac[n_S],range_yac[n_S],mean_yac3[n_S];
	//segmentation variables
	int subIndx = 0;
	int subSeg =4;
	int nft = 9;
	float xftrs[1000][nft][subSeg];

	/* motions variables */ 
	int w_nfeatures,w_nOutputs,w_result,
	    ws_nfeatures,ws_nOutputs,ws_result,
	    r_nfeatures,r_nOutputs,r_result,
	    rs_nfeatures,rs_nOutputs,rs_result,
	    j_nfeatures,j_nOutputs,j_result,
	    jl_nfeatures,jl_nOutputs,jl_result,
	    t_nfeatures,t_nOutputs,t_result,
	    sa_nfeatures,sa_nOutputs,sa_result,
	    sas_nfeatures,sas_nOutputs,sas_result,
	    sd_nfeatures,sd_nOutputs,sd_result,
	    sds_nfeatures,sds_nOutputs,sds_result;
	w_nfeatures=2;w_nOutputs=2;w_result=0;
	ws_nfeatures=2,ws_nOutputs=2,ws_result=0,
	r_nfeatures=2;r_nOutputs=2;r_result=0;
        rs_nfeatures=2,rs_nOutputs=2,rs_result=0,
	j_nfeatures=2;j_nOutputs=2;j_result=0;
        jl_nfeatures=2,jl_nOutputs=2,jl_result=0,
        t_nfeatures=1;t_nOutputs=3;t_result=0;
        sa_nfeatures=4;sa_nOutputs=2;sa_result=0;
        sas_nfeatures=2,sas_nOutputs=2,sas_result=0,
        sd_nfeatures=4;sd_nOutputs=2;sd_result=0,
        sds_nfeatures=2,sds_nOutputs=2,sds_result=0;


    	struct fann *w_ann,
		    *ws_ann,
		    *r_ann,
		    *rs_ann,
		    *j_ann,
		    *jl_ann,
		    *t_ann,
		    *sa_ann,
		    *sas_ann,
		    *sd_ann,
		    *sds_ann;
    	w_ann = fann_create_from_file("w_TEST.net");
	ws_ann = fann_create_from_file("./ws_TEST.net");
  	r_ann = fann_create_from_file("r_TEST.net");
	rs_ann = fann_create_from_file("./rs_TEST.net");
	j_ann = fann_create_from_file("j_TEST.net");
	jl_ann = fann_create_from_file("./jl_TEST.net");
	t_ann = fann_create_from_file("t_TEST.net");
	sa_ann = fann_create_from_file("sa_TEST.net");
	sas_ann = fann_create_from_file("./sas_TEST.net");
	sd_ann = fann_create_from_file("./sd_TEST.net");
	sds_ann = fann_create_from_file("../sds_TEST.net");

	float *w_features[]={mean_yac,variance_xac};
	float *ws_features[]={mean_xac,max_xgy};
	float *r_features[]={mean_yac3,mean_xac3};
	float *rs_features[]={min_xac,max_xgy};
	float *j_features[]={mean_xac3,min_xac4};
	float *jl_features[]={max_zgy,max_xgy};
	float *t_features[]={mean_xgy};
	float *sa_features[]={ratio_yac,max_yac,skewness_yac,skewness_zgy};
	float *sas_features[]={range_yac,variance_yac};
	float *sd_features[]={mean_zgy,skewness_yac,mean_yac,variance_yac};
	float *sds_features[]={range_yac,variance_yac};
	
	////begin of reading the files/////	
	while(1)
	{
		signal(SIGINT, &do_when_interrupted);
		if(!run_flag){ break;}
		if(icntr==10){icntr=0;}
//		sprintf(number,"%d",icntr);
//		strcpy(ifile_name,ifile);
//		strcat(ifile_name,number);
//		strcat(ifile_name,iextension);
		fp = fopen(ifile_name[icntr], "r");
		if (fp==NULL) {	
			fprintf(stderr, "Failed to open file \'%s\'. Exiting.\n"
					,ifile_name[icntr]);
			exit(EXIT_FAILURE);
		}
		fd = fileno(fp);
		flock(fd,LOCK_EX);
		//printf("Attempting to read from file \'%s\'.\n", ifile_name);

		/* count the number of lines in the file */
		read = getline(&line, &len, fp); //discard header of file
		N_SAMPLES = 0;
		while ((read = getline(&line, &len, fp)) != -1) {
			N_SAMPLES++;
		}
	
		/* go back to the start of the file so that the data can be read */
		rewind(fp);
		read = getline(&line, &len, fp); //discard header of file

		/* start reading the data from the file into the data structures */
		i = 0;
		t =   (double *) malloc(sizeof(double) * N_SAMPLES);
		t_b = (double *) malloc(sizeof(double) * N_SAMPLES);
		t_a = (double *) malloc(sizeof(double) * N_SAMPLES);
		x_ac = (float *) malloc(sizeof(float) * N_SAMPLES);
		y_ac = (float *) malloc(sizeof(float) * N_SAMPLES);
		z_ac = (float *) malloc(sizeof(float) * N_SAMPLES);
		x_gy = (float *) malloc(sizeof(float) * N_SAMPLES);
		y_gy = (float *) malloc(sizeof(float) * N_SAMPLES);
		z_gy = (float *) malloc(sizeof(float) * N_SAMPLES);
	
		while ((read = getline(&line, &len, fp)) != -1) {
			/* parse the data */
			rv = sscanf(line, "%lf,%lf,%f,%f,%f,%f,%f,%f\n",
			&t_b[i],&t_a[i],&x_ac[i],&y_ac[i],&z_ac[i],&x_gy[i],
		       	&y_gy[i], &z_gy[i]);
			if (rv != 8) {
				fprintf(stderr,
						"%s %d \'%s\'. %s.\n",
						"Failed to read line",
						i,
						line,
						"Exiting"
				       );
				exit(EXIT_FAILURE);
			}
			i++;
		}
		fclose(fp);	
		remove(ifile_name[icntr]); 

		for ( ii = 0; ii < N_SAMPLES; ii++)
		{
			t[ii] = (t_b[ii] + t_a[ii])/2.0;
	
		}
	
		/* 
		 * From selected thresholds, 
		 * find indicies of peaks
		 * find indicies of troughs
		 */
	
		ind_thr	= 300;
		P_i = (float *) malloc(sizeof(float) * N_SAMPLES);
		T_i = (float *) malloc(sizeof(float) * N_SAMPLES);
		rv = find_peaks_and_troughs(z_gy,N_SAMPLES,pk_threshold,P_i, 
				T_i, &n_P, &n_T);
		if (rv < 0) {
			fprintf(stderr, "find_peaks_and_troughs failed\n");
			exit(EXIT_FAILURE);
		}

		S_imax = (float *) malloc(sizeof(float) * BUFF_SIZE); // P
		S_imin = (float *) malloc(sizeof(float) * BUFF_SIZE); //T 

	        n_S = stride_extraction(ind_thr,n_P,P_i,T_i,S_imax,S_imin);	

        	//////extracting features///////////////////
		/*x_ac*/
		for(k=0;k<n_S;k++)
		{
			if( S_imax[k]>250){ offset = 250;}
			else{ continue;}
			start =abs((int)S_imax[k]-offset); //shift to the valley
			end = (int)S_imax[k]+offset;
			mean_xac[k] = calculate_mean(x_ac,start,end);
			calculate_Max_Min_Range(x_ac,start,end,
					(holder+k),(min_xac+k),
					(holder+k));
			calculate_Statistics (x_ac,start,end,mean_xac[k],
				(holder+k),(variance_xac+k),(holder+k),
				(holder+k),(holder+k));
		}	
		for(k=0;k<n_S;k++)
		{
			if( S_imax[k]>250){ offset = 250;}
			else {continue;} 
			start = abs((int)S_imax[k]-offset); //shift to  valley
			end = (int)S_imax[k]+offset;
			subIndx = (end-start)/subSeg;
			for(j=0;j<subSeg;j++)
			{
				end = start + subIndx;
		    		xftrs[k][0][j] = calculate_mean(x_ac,start,end);
        	    		calculate_Max_Min_Range(x_ac,start,end,
					&xftrs[k][1][j],&xftrs[k][2][j],
					&xftrs[k][3][j]);
		       	    	calculate_Statistics(x_ac,start,end,
					xftrs[k][0][j],&xftrs[k][4][j],
					&xftrs[k][5][j],&xftrs[k][6][j],
					&xftrs[k][7][j],&xftrs[k][8][j]);
				start = end;
		    
			}
			mean_xac3[k]= xftrs[k][0][2];         //mean subseg 3
		        min_xac4[k] = xftrs[k][2][3];	

		}
		
		/*x_gy*/
		for(k=0;k<n_S;k++)
		{
			if( S_imax[k]>250){ offset = 250;}
			else{ continue;}
			start =abs((int)S_imax[k]-offset); //shift to the valley
			end = (int)S_imax[k]+offset;
			mean_xgy[k] = calculate_mean(x_gy,start,end);
			calculate_Max_Min_Range(x_gy,start,end,
					(max_xgy+k),(holder+k),
					(holder+k));
		}	
		/*y_ac*/
		for(k=0;k<n_S;k++)
		{
			if( S_imax[k]>250){ offset = 250;}
			else{ continue;}
			start =abs((int)S_imax[k]-offset); //shift to the valley
			end = (int)S_imax[k]+offset;
			mean_yac[k] = calculate_mean(y_ac,start,end);
			calculate_Max_Min_Range(y_ac,start,end,
					(max_yac+k),(min_yac+k),
					(holder+k));
			calculate_Statistics (y_ac,start,end,mean_yac[k],
				(holder+k),(variance_yac+k),(holder+k),
				(skewness_yac+k),(holder+k));
			ratio_yac[k] = max_yac[k]/min_yac[k];
		}
	
		for(k=0;k<n_S;k++)
		{
			if( S_imax[k]>250){ offset = 250;}
			else {continue;} 
			start = abs((int)S_imax[k]-offset); //shift to  valley
			end = (int)S_imax[k]+offset;
			subIndx = (end-start)/subSeg;
			for(j=0;j<subSeg;j++)
			{	
				end = start + subIndx;
		    		xftrs[k][0][j] = calculate_mean(y_ac,start,end);
        	    		calculate_Max_Min_Range(y_ac,start,end,
					&xftrs[k][1][j],&xftrs[k][2][j],
					&xftrs[k][3][j]);
       	    	    		calculate_Statistics(y_ac,start,end,
					xftrs[k][0][j],&xftrs[k][4][j],
					&xftrs[k][5][j],&xftrs[k][6][j],
					&xftrs[k][7][j],&xftrs[k][8][j]);
				start = end;
		    
			}
			mean_yac3[k]= xftrs[k][0][2];         //mean subseg 3 

		}


		/*z_gy*/
		for(k=0;k<n_S;k++)
		{
			if( S_imax[k]>250){ offset = 250;}
			else{ continue;}
			start =abs((int)S_imax[k]-offset); //shift to  valley
			end = (int)S_imax[k]+offset;
			mean_zgy[k] = calculate_mean(z_gy,start,end);
			calculate_Max_Min_Range(z_gy,start,end,
					(max_zgy+k),(holder+k),
					(holder+k));
			calculate_Statistics (z_gy,start,end,mean_zgy[k],
				(holder+k),(holder+k),(holder+k),
				(skewness_zgy+k),(kurtosis_zgy+k));
		}


		if (n_S==0)
		{
			printf("standing\n");
		}	
		for (i = 0; i < (n_S); i++) {
			//walking network
			for(j=0;j<w_nfeatures;j++)
			{
		    		input[j]=(float)(w_features[j][i]/100);
			}		
			calc_out = fann_run(w_ann, input);

			max = -100;
		        for (ii = 0; ii < w_nOutputs; ii++) {
       			    if (calc_out[ii] > max) {
               			max = calc_out[ii];
	       	       	 	w_result = ii;
        	    	    }
			}
			//walking speeds network
			for(j=0;j<ws_nfeatures;j++)
			{
		    		input[j]=(float)(ws_features[j][i]/100);
			}		
			calc_out = fann_run(ws_ann, input);

			max = -100;
		        for (ii = 0; ii < ws_nOutputs; ii++) {
       			    if (calc_out[ii] > max) {
               			max = calc_out[ii];
	       	       	 	ws_result = ii;
        	    	    }
			}
			//running network
			for(j=0;j<r_nfeatures;j++)
			{
		    		input[j]=(float)(r_features[j][i]/100);
			}		
			calc_out = fann_run(r_ann, input);

			max = -100;
		        for (ii = 0; ii < r_nOutputs; ii++) {
       			    if (calc_out[ii] > max) {
               			max = calc_out[ii];
	       	       	 	r_result = ii;
        	    	    }
			}
			//running speeds network
			for(j=0;j<rs_nfeatures;j++)
			{
		    		input[j]=(float)(rs_features[j][i]/100);
			}		
			calc_out = fann_run(rs_ann, input);

			max = -100;
		        for (ii = 0; ii < rs_nOutputs; ii++) {
       			    if (calc_out[ii] > max) {
               			max = calc_out[ii];
	       	       	 	rs_result = ii;
        	    	    }
			}
			//jumping network
			for(j=0;j<j_nfeatures;j++)
			{
		    		input[j]=(float)(j_features[j][i]/100);
			}		
			calc_out = fann_run(j_ann, input);

			max = -100;
		        for (ii = 0; ii < j_nOutputs; ii++) {
       			    if (calc_out[ii] > max) {
               			max = calc_out[ii];
	       	       	 	j_result = ii;
        	    	    }
			}
			//jumping levels network
			for(j=0;j<jl_nfeatures;j++)
			{
		    		input[j]=(float)(jl_features[j][i]/100);
			}		
			calc_out = fann_run(jl_ann, input);

			max = -100;
		        for (ii = 0; ii < jl_nOutputs; ii++) {
       			    if (calc_out[ii] > max) {
               			max = calc_out[ii];
	       	       	 	jl_result = ii;
        	    	    }
			}

			//turning network
			for(j=0;j<t_nfeatures;j++)
			{
		    		input[j]=(float)(t_features[j][i]/100);
			}		
			calc_out = fann_run(t_ann, input);

			max = -100;
		        for (ii = 0; ii < t_nOutputs; ii++) {
       			    if (calc_out[ii] > max) {
               			max = calc_out[ii];
	       	       	 	t_result = ii;
        	    	    }
			}
			//stair Ascent network
			for(j=0;j<sa_nfeatures;j++)
			{
		    		input[j]=(float)(sa_features[j][i]/100);
			}		
			calc_out = fann_run(sa_ann, input);

			max = -100;
		        for (ii = 0; ii < sa_nOutputs; ii++) {
       			    if (calc_out[ii] > max) {
               			max = calc_out[ii];
	       	       	 	sa_result = ii;
        	    	    }
			}
			//stair Ascent speeds network
			for(j=0;j<sas_nfeatures;j++)
			{
		    		input[j]=(float)(sas_features[j][i]/100);
			}		
			calc_out = fann_run(sas_ann, input);

			max = -100;
		        for (ii = 0; ii < sas_nOutputs; ii++) {
       			    if (calc_out[ii] > max) {
               			max = calc_out[ii];
	       	       	 	sas_result = ii;
        	    	    }
			}

			//stair Descent network
			for(j=0;j<sd_nfeatures;j++)
			{
		    		input[j]=(float)(sd_features[j][i]/100);
			}		
			calc_out = fann_run(sd_ann, input);

			max = -100;
		        for (ii = 0; ii < sd_nOutputs; ii++) {
       			    if (calc_out[ii] > max) {
               			max = calc_out[ii];
	       	       	 	sd_result = ii;
        	    	    }
			}
			//stair Descent speeds network
			for(j=0;j<sds_nfeatures;j++)
			{
		    		input[j]=(float)(sds_features[j][i]/100);
			}		
			calc_out = fann_run(sds_ann, input);

			max = -100;
		        for (ii = 0; ii < sds_nOutputs; ii++) {
       			    if (calc_out[ii] > max) {
               			max = calc_out[ii];
	       	       	 	sds_result = ii;
        	    	    }
			}
			
			switch(w_result){
				case 0:
					printf("walking\n");
					break;
				case 1:
					printf("others\n");
					break;
				default:
					printf("error motion\n");
			}
			switch(r_result){
				case 0:
					printf("running\n");
					break;
				case 1:
					printf("others\n");
					break;
				default:
					printf("error motion\n");
			}
			switch(j_result){
				case 0:
					printf("jumping\n");
					break;
				case 1:
					printf("others\n");
					break;
				default:
					printf("error motion\n");
			}	
			switch(t_result){
				case 0:
					printf("left turn\n");
					break;
				case 1:
					printf("right\n");
					break;
				case 2:
					printf("others\n");
					break;

				default:
					printf("error motion\n");
			}	
			switch(sa_result){
				case 0:
					printf("stair ascent\n");
					break;
				case 1:
					printf("others\n");
					break;
				default:
					printf("error motion\n");
			}	
			switch(sd_result){
				case 0:
					printf("stair descent\n");
					break;
				case 1:
					printf("others\n");
					break;
				default:
					printf("error motion\n");
			}
			printf("\n");	

		}
		icntr++;
	}

	
	printf("extract_stride_data completed successfuly. Exiting.\n");
   	fann_destroy(w_ann);
        fann_destroy(ws_ann);
  	fann_destroy(r_ann);
  	fann_destroy(rs_ann);
  	fann_destroy(j_ann);
  	fann_destroy(jl_ann);
  	fann_destroy(t_ann);
  	fann_destroy(sa_ann);
  	fann_destroy(sas_ann);
  	fann_destroy(sd_ann);
  	fann_destroy(sds_ann);


/*	free(t); free(t_a); free(t_b); free(x_ac); free(y_ac); free(z_ac);
	free(x_gy); free(y_gy); free(z_gy); free(S_imax); free(S_imin); 
	free(P_i); free(T_i); 
*/	
	
	return 0;


}
///////////////////////////
void featureExtraction(int n_S,float *S_imax,float **features,float *axis)
{	
	float *mean,*max,*min,*range,*MAD,*variance,
	      *std,*skewness,*kurtosis;
        mean = features[0];           max = features[1];
	min = features[2];            range = features[3];
	MAD = features[4];            variance = features[5];
	std = features[6];            skewness = features[7];
	kurtosis = features[8];

	int start,end;
	int k, offset;
	for(k=0;k<n_S;k++)
	{
		if((k+1)!=n_S){	offset = (S_imax[k+1]-S_imax[k])/2;}
		else{ offset = 200;}

		start = S_imax[k]-offset; //shift from the pick to the valley
		end = S_imax[k]+offset;
		mean[k] = calculate_mean(axis,start,end);
        	calculate_Max_Min_Range(axis,start,end,(max+k),(min+k),
				(range+k));
        	calculate_Statistics (axis,start,end,mean[k],(MAD+k),
				(variance+k),(std+k),(skewness+k),(kurtosis+k));
	}

}


void clear_buffer(float *arr, float val, int n) 
{
	int i;
	for (i = 0; i < n; i++) {
		arr[i] = val;
	}
}

/*
 * Caculates mean of first <n> samples in <*arr>
 */
float calculate_mean(float *arr, int start, int end)
{
	float total;
	int i, n;

	n = end - start;
	total = 0.0f;
	for (i = 0; i < n; i++) {
		total += arr[start + i];
	}

	return total/(n);
}

void calculate_Max_Min_Range(float *arr,int start,int end,
		float *max, float *min, float *range)
{

	int i, n;

	n = end - start;
	
	*max = arr[start];
	*min = arr[start];
	
	for (i = 0; i < n; i++) {
		if(arr[start+i]> *max)
		{
			*max = arr[start+i];
		}
		if(arr[start+i]< *min)
		{
			*min = arr[start+i];
		}
	}
	*range = (*max-*min);
}


void calculate_Statistics (float *arr, int start, int end, float mean, 
		float *MAD, float *variance, float *std,
		float *skewness, float *kurtosis)
{	
	float total1, total2, total3, total4, holder; 
	total1 = 0.0f;	total2 = 0.0f; total3 = 0.0f;
	total4 = 0.0f;  holder = 0.0f;
	int n = end - start;
	int i;
	for (i = 0; i < n; i++) {
		holder = arr[start + i]-mean;
		total1 += abs(holder);
		total2 += holder*holder;
		total3 += holder*holder*holder;
		total4 += holder*holder*holder*holder;	
	}

	*MAD = total1/(n);
	*variance = total2/(n);
	*std = sqrt(*variance);
	holder = *std;
	*skewness = (total3/(n))/(holder*holder*holder);
	*kurtosis = (total4/(n))/(holder*holder*holder*holder);
}


int find_peaks_and_troughs(
		float *arr, 	// signal 
		int n_samples, 	// number of samples present in the signal
		float E, 	// threshold for peak detection
		// arrays that will store the indicies of the located
		// peaks and troughs
		float *P, float *T,
		// number of peaks (n_P) and number of troughs (n_T)
		// found in the data set *arr
		int *n_P, int *n_T
		)
{
	int a, b, i, d, _n_P, _n_T;

	i = -1; d = 0; a = 0; b = 0;
	_n_P = 0; _n_T = 0;

	clear_buffer(P, 0.0f, n_samples);
	clear_buffer(T, 0.0f, n_samples);

	while (i != n_samples) {
		i++;
		if (d == 0) {
			if (arr[a] >= (arr[i] + E)) {
				d = 2;
			} else if (arr[i] >= (arr[b] + E)) {
				d = 1;
			}
			if (arr[a] <= arr[i]) {
				a = i;
			} else if (arr[i] <= arr[b]) {
				b = i;
			}
		} else if (d == 1) {
			if (arr[a] <= arr[i]) {
				a = i;
			} else if (arr[a] >= (arr[i] + E)) {
				/*
				 * Peak has been detected.
				 * Add index at detected peak
				 * to array of peak indicies
				 * increment count of peak indicies
				 */
				P[_n_P] = a;
				_n_P++;
				b = i;
				d = 2;
			}
		} else if (d == 2) {
			if (arr[i] <= arr[b]) {
				b = i;
			} else if (arr[i] >= (arr[b] + E)) {
				/*
				 * Trough has been detected.
				 * Add index at detected trough
				 * to array of trough indicies
				 * increment count of trough indicies
				 */
				T[_n_T] = b;
				_n_T++;
				a = i;
				d = 1;
			}
		}
	}

	(*n_P) = _n_P;
	(*n_T) = _n_T;
	return 0;
}

int stride_extraction(int ind_thr,int n_P, float *P_i,float *T_i,
		 float *S_imax, float *S_imin)
{
	int a , b , c, n_S;
	a = 0; b = 1; n_S = 0; c = 0; 
	while( b < n_P)
        {	if((P_i[b] - P_i[a]) > ind_thr)
		{
			S_imax[n_S] = P_i[a];
			S_imin[n_S] = T_i[a];
			n_S++;
			a=a+c+1;              
			b=b+1;
			c=0;
			/*checking for boundary conditions: if the last stride*/ 
			if((b == n_P) && ((P_i[b-1] - P_i[a-1])>300)) //P
		        {
				S_imax[n_S] = P_i[b-1];	//P
				S_imin[n_S] = T_i[b-1];  //T
				n_S++;
			}	
		}
		else{
			b++;
			c++;
                        /*checking for boundary conditions: if the last stride*/
			if((b == n_P) && ((P_i[a] - P_i[a-1])>300))  //P
			{		
				S_imax[n_S] = P_i[a];   //P
                        	S_imin[n_S] = T_i[a];  //T
				n_S++;
			}	
		}
	
	}
	return n_S;
}






