/* for file and terminal I/O */
#include <stdio.h>
/* for string manip */
#include <string.h>
/* for exit() */
#include <stdlib.h>
/* for fabsf() */
#include <math.h>

#define BUFF_SIZE 1024
/*
 * sets first <n> values in <*arr> to <val>
 */
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
void peaks_throughs_file(FILE *fp,char* ofile_pt_name, int n_P,float *P_i,
		int n_T,float *T_i, double *t,float *axis);
void stride_file(FILE *fp,char* ofile_st_name, int n_S,float *S_imax,
		float *S_imin, double *t,float *axis);
void training_file(FILE *fp,char* train_file_name, int n_S,float *S_imax,
		float **features,int activityCode,double *t,int nfeatures,int
		nOutputs);
void featureExtraction(int n_S,float *S_imax,float **features,float *axis);
void training_file2(FILE *fp,char* train_file_name, int n_S,float *S_imax,
		int activityCode,double *t,int nfeatures,int subSeg,
		float Features[n_S][nfeatures][subSeg]);

int main(int argc, char **argv)
{
	/* Generic variables */
	int i, idx, ind_thr;
	int rv;
	/* Variables for reading file line by line */
	char *ifile_name, *ofile_pt_name, *ofile_st_name;
	FILE *fp;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int N_SAMPLES;

	char *activity;

	/* Variables for storing the data and storing the return values */
	double *t,*t_b, *t_a; 
	/* variables for data collected from input file */
	float *x_ac, *y_ac, *z_ac, *x_gy, *y_gy, *z_gy; 	
	float pk_threshold;	// pk-threshold value
       	/* Variables for peak-trough detection */	
	float *P_i; 	// indicies of each peak found by peak detection
	float *T_i; 	// indicies of each trough found by trough detection
	float *S_imax,*S_imin; 	// indicies of the start of each stride
	int n_P; 	// number of peaks
	int n_T; 	// number of troughs
	int n_S = 0; 	// number of strides
	
	/*
	 * Check if the user entered the correct command line arguments
	 * Usage: 
	 * ./extract_stride_data <ifile_name> <output_peaks> <output_strides>
	 * 				<threshold_value_float>
	 * Or 
	 * ./extract_stride_data
	 */
	if (argc != 6) {
		ifile_name = (char *) malloc(sizeof(char) * BUFF_SIZE);
		memset(ifile_name, 0, BUFF_SIZE);
		snprintf(ifile_name, 
				BUFF_SIZE, 
				"data.csv"
			);
		ofile_pt_name = (char *) malloc(sizeof(char) * BUFF_SIZE);
		memset(ofile_pt_name, 0, BUFF_SIZE);
		snprintf(ofile_pt_name, BUFF_SIZE, "pt_output.csv");
		ofile_st_name = (char *) malloc(sizeof(char) * BUFF_SIZE);
		memset(ofile_st_name, 0, BUFF_SIZE);
		snprintf(ofile_st_name, BUFF_SIZE, "strides.csv");
		pk_threshold = 100.0;
	} else {
		ifile_name = argv[1];
		ofile_pt_name = argv[2];
		ofile_st_name = argv[3];
		pk_threshold = atof(argv[4]);
		activity = argv[5];
	}

	int activityCode = 0 ;
	char act = *activity;	
	int number = atoi ((activity+1));
	printf("%d\n",number);
	switch (act)
	{
		case 'w':
			activityCode = -1 + number;
			break;
		case 't':
			activityCode = 2 + number;
			break;
		case 'j':
			activityCode = 4 + number;
			break;
		case 'u':
			activityCode = 6 + number;
			break;
		case 'd':
			activityCode = 9 + number;
			break;
                case 'r':
			activityCode = 12 + number;
			break;
		default :
			printf("invalid activity code \n");
	}

	printf("Arguments used:\n\t%s=%s\n\t%s=%s\n\t%s=%s\n\t%s=%f\n",
			"ifile_name", ifile_name,
			"ofile_peak_trough_name", ofile_pt_name,
			"ofile_stride_name", ofile_st_name,
			"peak_threshold", pk_threshold
	      );

	/* open the input file */
	printf("Attempting to read from file \'%s\'.\n", ifile_name);
	fp = fopen(ifile_name, "r");
	if (fp == NULL) {
		fprintf(stderr, 
				"Failed to read from file \'%s\'.\n", 
				ifile_name
		       );
		exit(EXIT_FAILURE);
	}

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
		rv = sscanf(line, "%lf,%lf,%f,%f,%f,%f,%f,%f\n", &t_b[i], &t_a[i], &x_ac[i], &y_ac[i],&z_ac[i],&x_gy[i], &y_gy[i], &z_gy[i]);
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
        
	int ii;
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
	rv = find_peaks_and_troughs(z_gy,N_SAMPLES,pk_threshold,P_i, T_i, 
			&n_P, &n_T);
	if (rv < 0) {
		fprintf(stderr, "find_peaks_and_troughs failed\n");
		exit(EXIT_FAILURE);
	}

	S_imax = (float *) malloc(sizeof(float) * BUFF_SIZE); // P
	S_imin = (float *) malloc(sizeof(float) * BUFF_SIZE); //T 

        n_S = stride_extraction(ind_thr,n_P,P_i,T_i,S_imax,S_imin);	
	peaks_throughs_file(fp,ofile_pt_name,n_P,P_i,n_T,T_i,t,z_gy);
	stride_file(fp,ofile_st_name,n_S,S_imax,S_imin,t,z_gy);

	/*extracting features*/
	char * train_file_name = "jl_temp.txt";
	char * analysis_file_name = "jl_analysis.csv";
	int start,end;
	int j,k;
	float offset=0;
	float holder[n_S];
	int nfeatures;
	int nOutputs;

 	/*x_ac*/
	float mean_xac3[n_S],min_xac4[n_S],range_xac4[n_S],variance_xac4[n_S];
	int subIndx = 0;
	int subSeg =4;
	int nft = 9;
	float xftrs[n_S][nft][subSeg];
	//printf("n_S: %d\n",n_S);
	for(k=0;k<n_S;k++)
	{
		//if((k+1)!=n_S){	offset = (S_imax[k+1]-S_imax[k])/2;}
		if( S_imax[k]>250){ offset = 250;}
		else {continue;} 
		//printf("Indx: %f\n",S_imax[k]);
		//printf("offset: %f\n",offset);
		start = abs((int)(S_imax[k]-offset)); //shift from the pick to the valley
		end = (int)S_imax[k]+offset;
		subIndx = (end-start)/subSeg;
		//printf("subIndx: %d\n",subIndx);
		//printf("start %d\n",start);
		for(j=0;j<subSeg;j++)
		{
			end = start + subIndx;
		        //printf("segEnd %d\n",end);

		    	xftrs[k][0][j] = calculate_mean(x_ac,start,end);
        	    	calculate_Max_Min_Range(x_ac,start,end,&xftrs[k][1][j],
				    &xftrs[k][2][j],&xftrs[k][3][j]);
        	    	calculate_Statistics(x_ac,start,end,xftrs[k][0][j],
				&xftrs[k][4][j],&xftrs[k][5][j],&xftrs[k][6][j],
				&xftrs[k][7][j],&xftrs[k][8][j]);
			start = end;
		    
		}
		mean_xac3[k]= xftrs[k][0][2];         //mean subseg 3 
		min_xac4[k] = xftrs[k][2][3];         //min  subseg 4
		range_xac4[k]= xftrs[k][3][3];        //range subseg 4 
		variance_xac4[k]=xftrs[k][5][3];      //var   subseg 4

		//printf("k: %d\n",k);

	}

	training_file2(fp,analysis_file_name,n_S,S_imax,activityCode,t, nft,
		subSeg,	xftrs);

	/*y_ac*/
	float mean_yac[n_S],max_yac[n_S],min_yac[n_S],variance_yac[n_S],
	      skewness_yac[n_S],ratio_yac[n_S],
	      mean_yac3[n_S],max_yac4[n_S],range_yac4[n_S],kurtosis_yac2[n_S],
	      variance_yac4[n_S];

	
	for(k=0;k<n_S;k++)
	{
		//if((k+1)!=n_S){	offset =(S_imax[k+1]-S_imax[k])/2;}
		if( S_imax[k]>250){ offset = 250;}
		else { continue;}

		start = abs((int) S_imax[k]-offset); //shift from the pick to the valley
		end = (int)S_imax[k]+offset;
		mean_yac[k] = calculate_mean(y_ac,start,end);
        	calculate_Max_Min_Range(y_ac,start,end,(max_yac+k),(min_yac+k),
				(holder+k));
        	calculate_Statistics (y_ac,start,end,mean_yac[k],(holder+k),
	          (variance_yac+k),(holder+k),(skewness_yac+k),(holder+k));
		
		ratio_yac[k] = max_yac[k]/min_yac[k];
	}
/////////////seg
	for(k=0;k<n_S;k++)
	{
	        //if((k+1)!=n_S){offset = (S_imax[k+1]-S_imax[k])/2;}
		if( S_imax[k]>250){ offset = 250;}
		else {continue;} 
		//printf("Indx: %f\n",S_imax[k]);
		//printf("offset: %f\n",offset);
		start = abs((int)(S_imax[k]-offset)); //shift from the pick to the valley
		end = (int)S_imax[k]+offset;
		subIndx = (end-start)/subSeg;
		//printf("subIndx: %d\n",subIndx);
		//printf("start %d\n",start);
		for(j=0;j<subSeg;j++)
		{
			end = start + subIndx;
		        //printf("segEnd %d\n",end);

		    	xftrs[k][0][j] = calculate_mean(y_ac,start,end);
        	    	calculate_Max_Min_Range(y_ac,start,end,&xftrs[k][1][j],
				    &xftrs[k][2][j],&xftrs[k][3][j]);
        	    	calculate_Statistics(y_ac,start,end,xftrs[k][0][j],
				&xftrs[k][4][j],&xftrs[k][5][j],&xftrs[k][6][j],
				&xftrs[k][7][j],&xftrs[k][8][j]);
			start = end;
		    
		}
		mean_yac3[k] = xftrs[k][0][2];         //mean subseg 3 
		max_yac4[n_S] = xftrs[k][1][3];
		range_yac4[n_S] = xftrs[k][3][3];
		kurtosis_yac2[n_S] = xftrs[k][8][1];
	      	variance_yac4[n_S] = xftrs[k][5][3];

//	printf("k: %d\n",k);

	}

//	training_file2(fp,analysis_file_name,n_S,S_imax,activityCode,t, nft,
//		subSeg,	xftrs);


	/*z_ac*/
	float mean_zac[n_S],max_zac[n_S],min_zac[n_S],variance_zac[n_S];

	for(k=0;k<n_S;k++)
	{
		//if((k+1)!=n_S){offset = (S_imax[k+1]-S_imax[k])/2;}
		if( S_imax[k]>250){ offset = 250;}
		start = abs((int)S_imax[k]-offset); //shift from the pick to the valley
		end = (int)S_imax[k]+offset;
		mean_zac[k] = calculate_mean(z_ac,start,end);
        	calculate_Max_Min_Range(z_ac,start,end,(max_zac+k),(min_zac+k),
				(holder+k));
        	calculate_Statistics (z_ac,start,end,mean_zac[k],(holder+k),
			(variance_zac+k),(holder+k),(holder+k),(holder+k));
	}
/////////////seg
	for(k=0;k<n_S;k++)
	{
	        //if((k+1)!=n_S){offset = (S_imax[k+1]-S_imax[k])/2;}
		if( S_imax[k]>250){ offset = 250;}
		else {continue;} 
		//printf("Indx: %f\n",S_imax[k]);
		//printf("offset: %f\n",offset);
		start = abs((int)(S_imax[k]-offset)); //shift from the pick to the valley
		end = (int)S_imax[k]+offset;
		subIndx = (end-start)/subSeg;
		//printf("subIndx: %d\n",subIndx);
		//printf("start %d\n",start);
		for(j=0;j<subSeg;j++)
		{
			end = start + subIndx;
		        //printf("segEnd %d\n",end);

		    	xftrs[k][0][j] = calculate_mean(z_ac,start,end);
        	    	calculate_Max_Min_Range(z_ac,start,end,&xftrs[k][1][j],
				    &xftrs[k][2][j],&xftrs[k][3][j]);
        	    	calculate_Statistics(z_ac,start,end,xftrs[k][0][j],
				&xftrs[k][4][j],&xftrs[k][5][j],&xftrs[k][6][j],
				&xftrs[k][7][j],&xftrs[k][8][j]);
			start = end;
		    
		}
	//	printf("k: %d\n",k);

	}

//	training_file2(fp,analysis_file_name,n_S,S_imax,activityCode,t, nft,
//		subSeg,	xftrs);


	/*x_gy*/
	float mean_xgy[n_S],kurtosis_xgy[n_S],max_xgy[n_S],min_xgy[n_S];
	for(k=0;k<n_S;k++)
	{
		//if((k+1)!=n_S){offset = (S_imax[k+1]-S_imax[k])/2;}
		if( S_imax[k]>250){ offset = 250;}
		else{ continue;}
		start = abs((int)S_imax[k]-offset); //shift from the pick to the valley
		end = (int)S_imax[k]+offset;
		mean_xgy[k] = calculate_mean(x_gy,start,end);
 	      calculate_Max_Min_Range(x_gy,start,end,(max_xgy+k),(min_xgy+k),
			(holder+k));
       		calculate_Statistics (x_gy,start,end,mean_xgy[k],(holder+k),
	         (holder+k),(holder+k),(holder+k),(kurtosis_xgy+k));
	}

/////////////seg
	for(k=0;k<n_S;k++)
	{
	        //if((k+1)!=n_S){offset = (S_imax[k+1]-S_imax[k])/2;}
		if( S_imax[k]>250){ offset = 250;}
		else {continue;} 
		//printf("Indx: %f\n",S_imax[k]);
		//printf("offset: %f\n",offset);
		start = abs((int)(S_imax[k]-offset)); //shift from the pick to the valley
		end = (int)S_imax[k]+offset;
		subIndx = (end-start)/subSeg;
		//printf("subIndx: %d\n",subIndx);
		//printf("start %d\n",start);
		for(j=0;j<subSeg;j++)
		{
			end = start + subIndx;
		        //printf("segEnd %d\n",end);

		    	xftrs[k][0][j] = calculate_mean(x_gy,start,end);
        	    	calculate_Max_Min_Range(x_gy,start,end,&xftrs[k][1][j],
				    &xftrs[k][2][j],&xftrs[k][3][j]);
        	    	calculate_Statistics(x_gy,start,end,xftrs[k][0][j],
				&xftrs[k][4][j],&xftrs[k][5][j],&xftrs[k][6][j],
				&xftrs[k][7][j],&xftrs[k][8][j]);
			start = end;
		    
		}
	//	printf("k: %d\n",k);

	}

//	training_file2(fp,analysis_file_name,n_S,S_imax,activityCode,t, nft,
//		subSeg,	xftrs);



	/*y_gy*/
/////////////seg
	for(k=0;k<n_S;k++)
	{
	        //if((k+1)!=n_S){offset = (S_imax[k+1]-S_imax[k])/2;}
		if( S_imax[k]>250){ offset = 250;}
		else {continue;} 
		//printf("Indx: %f\n",S_imax[k]);
		//printf("offset: %f\n",offset);
		start = abs((int)(S_imax[k]-offset)); //shift from the pick to the valley
		end = (int)S_imax[k]+offset;
		subIndx = (end-start)/subSeg;
		//printf("subIndx: %d\n",subIndx);
		//printf("start %d\n",start);
		for(j=0;j<subSeg;j++)
		{
			end = start + subIndx;
		        //printf("segEnd %d\n",end);

		    	xftrs[k][0][j] = calculate_mean(y_gy,start,end);
        	    	calculate_Max_Min_Range(y_gy,start,end,&xftrs[k][1][j],
				    &xftrs[k][2][j],&xftrs[k][3][j]);
        	    	calculate_Statistics(y_gy,start,end,xftrs[k][0][j],
				&xftrs[k][4][j],&xftrs[k][5][j],&xftrs[k][6][j],
				&xftrs[k][7][j],&xftrs[k][8][j]);
			start = end;
		    
		}
	//	printf("k: %d\n",k);

	}

//	training_file2(fp,analysis_file_name,n_S,S_imax,activityCode,t, nft,
//		subSeg,	xftrs);



	/*z_gy*/
	float mean_zgy[n_S],max_zgy[n_S],min_zgy[n_S],MAD_zgy[n_S],
	      skewness_zgy[n_S],kurtosis_zgy[n_S];
	
	for(k=0;k<n_S;k++)
	{
		//if((k+1)!=n_S){offset =(S_imax[k+1]-S_imax[k])/2;}
		if( S_imax[k]>250){ offset = 250;}
		else { continue;}
		start = abs((int)S_imax[k]-offset); //shift from the pick to the valley
		end = (int)S_imax[k]+offset;
		mean_zgy[k] = calculate_mean(z_gy,start,end);
 	      	calculate_Max_Min_Range(z_ac,start,end,(max_zgy+k),(min_zgy+k),
				(holder+k));
       		calculate_Statistics (z_gy,start,end,mean_zgy[k],(MAD_zgy+k),
	         (holder+k),(holder+k),(skewness_zgy+k),(kurtosis_zgy+k));
	}

/////////////seg
	for(k=0;k<n_S;k++)
	{
	        //if((k+1)!=n_S){offset = (S_imax[k+1]-S_imax[k])/2;}
		if( S_imax[k]>250){ offset = 250;}
		else {continue;} 
		//printf("Indx: %f\n",S_imax[k]);
		//printf("offset: %f\n",offset);
		start = abs((int)(S_imax[k]-offset)); //shift from the pick to the valley
		end = (int)S_imax[k]+offset;
		subIndx = (end-start)/subSeg;
		//printf("subIndx: %d\n",subIndx);
		//printf("start %d\n",start);
		for(j=0;j<subSeg;j++)
		{
			end = start + subIndx;
		        //printf("segEnd %d\n",end);

		    	xftrs[k][0][j] = calculate_mean(z_gy,start,end);
        	    	calculate_Max_Min_Range(z_gy,start,end,&xftrs[k][1][j],
				    &xftrs[k][2][j],&xftrs[k][3][j]);
        	    	calculate_Statistics(z_gy,start,end,xftrs[k][0][j],
				&xftrs[k][4][j],&xftrs[k][5][j],&xftrs[k][6][j],
				&xftrs[k][7][j],&xftrs[k][8][j]);
			start = end;
		    
		}
	//	printf("k: %d\n",k);

	}

//	training_file2(fp,analysis_file_name,n_S,S_imax,activityCode,t, nft,
//		subSeg,	xftrs);




	float *features[]={max_xgy,max_zgy};


	nfeatures = 2;
	nOutputs = 2;
	training_file(fp,train_file_name,n_S,S_imax,features,activityCode,t,
			nfeatures,nOutputs);



	printf("extract_stride_data completed successfuly. Exiting.\n");
/*	free(t); free(t_a); free(t_b); free(x_ac); free(y_ac); free(z_ac);
	free(x_gy); free(y_gy); free(z_gy); free(S_imax); free(S_imin); 
	free(ifile_name); free(ofile_st_name); free(ofile_pt_name); 
	free(P_i); free(T_i); 
*/	return 0;
}

////////////////////////////

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

void training_file(FILE *fp,char* train_file_name, int n_S,float *S_imax,
	float **Features,int activityCode,double *t,int nfeatures,int nOutputs)
{
	double period;
	int i,j,k,idx_max,idx_next;
	
	int arrCode[nOutputs];
	
	for(i=0; i<nOutputs; i++)
	{
		arrCode[i] = -1;
	}
	arrCode[activityCode] = 1;

        /* open the training file for the neutal network */
	printf("Attempting to write to file \'%s\'.\n", train_file_name);
	fp = fopen(train_file_name, "a");
	if (fp == NULL) {
		fprintf(stderr, 
				"Failed to write to file \'%s\'.\n", 
				train_file_name
		       );
		exit(EXIT_FAILURE);
	}
        
	idx_max = 0;
	idx_next = 0;
//	fprintf(fp,"start\n");
//	

	for (i = 0; i < (n_S); i++) {
		//idx_max = (int) S_imax[i];
		//idx_next = (int) S_imax[i+1];
		//((i+1)!=n_S)? period = t[idx_next]- t[idx_max]: period;
		//fprintf(fp, "%lf",period/10.0);
		
		for(j=0;j<nfeatures;j++)
		{
			fprintf(fp,"%f ",Features[j][i]/100);
		}	
		fprintf(fp, "\n");
		for(k=0;k<nOutputs;k++)
		{
			fprintf(fp,"%d ",arrCode[k]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void training_file2(FILE *fp,char* train_file_name, int n_S,float *S_imax,
		int activityCode,double *t,int nfeatures,int subSeg,
		float Features[n_S][nfeatures][subSeg])
{
	double period;
	int i,j,k,idx_max,idx_next;
	int nOutputs = 5;
	int arrCode[nOutputs];
	
	for(i=0; i<nOutputs; i++)
	{
		arrCode[i] = -1;
	}
	arrCode[activityCode] = 1;

        /* open the training file for the neutal network */
	printf("Attempting to write to file \'%s\'.\n", train_file_name);
	fp = fopen(train_file_name, "a");
	if (fp == NULL) {
		fprintf(stderr, 
				"Failed to write to file \'%s\'.\n", 
				train_file_name
		       );
		exit(EXIT_FAILURE);
	}
        
	idx_max = 0;
	idx_next = 0;
	fprintf(fp,"start\n");

	for (i = 0; i < (n_S); i++) {
		//idx_max = (int) S_imax[i];
		//idx_next = (int) S_imax[i+1];
		//((i+1)!=n_S)? period = t[idx_next]- t[idx_max]: period;
		//fprintf(fp, "%lf",period/10.0);
		
		for(j=0;j<nfeatures;j++)
		{
			for(k=0;k<subSeg;k++)
			{
				fprintf(fp,"%f,",Features[i][j][k]);
		
			}	
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
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

	return total/((float) n);
}

void calculate_Max_Min_Range(float *arr,int start,int end,
		float *max, float *min, float *range)
{
	float total;
	int i, n;

	n = end - start;
	total = 0.0f;
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

	*MAD = total1/((float) n);
	*variance = total2/((float) n);
	*std = sqrt(*variance);
	holder = *std;
	*skewness = (total3/((float) n))/(holder*holder*holder);
	*kurtosis = (total4/((float) n))/(holder*holder*holder*holder);
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
	while( b < n_P)                     //P
	{
		if((P_i[b] - P_i[a]) > ind_thr)   //P
		{
			S_imax[n_S] = P_i[a];         //P
		        S_imin[n_S] = T_i[a];	  // T
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

void peaks_throughs_file(FILE *fp,char* ofile_pt_name, int n_P,float *P_i,
		int n_T,float *T_i, double *t,float *axis)
{
	int i, idx;
	/* open the output file to write the peak and trough data */
	printf("Attempting to write to file \'%s\'.\n", ofile_pt_name);
	fp = fopen(ofile_pt_name, "w");
	if (fp == NULL) {
		fprintf(stderr, 
				"Failed to write to file \'%s\'.\n", 
				ofile_pt_name
		       );
		exit(EXIT_FAILURE);
	}

	fprintf(fp, "P_i,P_t,P_v,T_i,T_t,T_v\n");            
	for (i = 0; i < n_P || i < n_T; i++) {
		/* Only peak data if there is peak data to write */
		if (i < n_P) {
			idx = (int) P_i[i];
			fprintf(fp, "%d,%20.10lf,%lf,",
					idx,
					t[idx],
					axis[idx]               //x
			       );
		} else {
			fprintf(fp, ",,,");
		}
		/* Only trough data if there is trough data to write */
		if (i < n_T) {
			idx = (int) T_i[i];
			fprintf(fp, "%d,%20.10lf,%lf\n",
					idx,
					t[idx],
					axis[idx]                //x
			       );
		} else {
			fprintf(fp, ",,\n");
		}
	}
	fclose(fp);
}

void stride_file(FILE *fp,char* ofile_st_name, int n_S,float *S_imax,
		float *S_imin, double *t,float *axis)
{
	double period;
	int i,idx_max,idx_min,idx_next;
        /* open the output file to write the stride data */
	printf("Attempting to write to file \'%s\'.\n", ofile_st_name);
	fp = fopen(ofile_st_name, "w");
	if (fp == NULL) {
		fprintf(stderr, 
				"Failed to write to file \'%s\'.\n", 
				ofile_st_name
		       );
		exit(EXIT_FAILURE);
	}
       	fprintf(fp, "S_i,S_t,S_max,S_i,S_t,S_min,period\n");          
	for (i = 0; i < n_S; i++) {
		idx_max = (int) S_imax[i];
		idx_min = (int) S_imin[i];
		idx_next = (int) S_imax[i+1];
		((i+1)!=n_S)? period = t[idx_next]- t[idx_max]: period;
		fprintf(fp, "%d,%20.10lf,%f,%d,%20.10lf,%f,%lf\n",
				idx_max,
				t[idx_max],
				axis[idx_max],		
				idx_min,
				t[idx_min],
                                axis[idx_min],  
				period
		       );
	}
	fclose(fp);
}



