}
int main(int argc, char **argv)
{
	

	/*extracting features*/
	char * train_file_name = "train_set.txt";
	int start,end;
	int j,k;
	float offset=0;
	float holder[n_S];
	int nfeatures = 19;
        /*x_ac*/
	/*y_ac*/
	float mean_yac[n_S],max_yac[n_S],min_yac[n_S],variance_yac[n_S],
	      skewness_yac[n_S],ratio_yac[n_S];

	
	for(k=0;k<n_S;k++)
	{
		if((k+1)!=n_S){	offset = (S_imax[k+1]-S_imax[k])/2;}
		else if( S_imax[k]>200){ offset = 200;}

		start = (int) S_imax[k]-offset; //shift from the pick to the valley
		end = (int)S_imax[k]+offset;
		mean_yac[k] = calculate_mean(y_ac,start,end);
        	calculate_Max_Min_Range(y_ac,start,end,(max_yac+k),(min_yac+k),
				(holder+k));
        	calculate_Statistics (y_ac,start,end,mean_yac[k],(holder+k),
	          (variance_yac+k),(holder+k),(skewness_yac+k),(holder+k));
		
		ratio_yac[k] = max_yac[k]/min_yac[k];
	}




	/*z_ac*/
	float mean_zac[n_S],max_zac[n_S],min_zac[n_S],variance_zac[n_S];

	for(k=0;k<n_S;k++)
	{
		if((k+1)!=n_S){	offset = (S_imax[k+1]-S_imax[k])/2;}
		else if( S_imax[k]>200){ offset = 200;}
		start = (int)S_imax[k]-offset; //shift from the pick to the valley
		end = (int)S_imax[k]+offset;
		mean_zac[k] = calculate_mean(z_ac,start,end);
        	calculate_Max_Min_Range(z_ac,start,end,(max_zac+k),(min_zac+k),
				(holder+k));
        	calculate_Statistics (z_ac,start,end,mean_zac[k],(holder+k),
			(variance_zac+k),(holder+k),(holder+k),(holder+k));
	}

	/*x_gy*/
	float mean_xgy[n_S],kurtosis_xgy[n_S];
	for(k=0;k<n_S;k++)
	{
		if((k+1)!=n_S){	offset = (S_imax[k+1]-S_imax[k])/2;}
		else if( S_imax[k]>200){ offset = 200;}
		start = (int)S_imax[k]-offset; //shift from the pick to the valley
		end = (int)S_imax[k]+offset;
		mean_xgy[k] = calculate_mean(x_gy,start,end);
 	      //calculate_Max_Min_Range(x_ac,start,end,(max_xgy+k),(min_xgy+k),
		//		(holder+k));
       		calculate_Statistics (x_gy,start,end,mean_xgy[k],(holder+k),
	         (holder+k),(holder+k),(holder+k),(kurtosis_xgy+k));
	}
	/*y_gy*/

	/*z_gy*/
	float mean_zgy[n_S],max_zgy[n_S],min_zgy[n_S],MAD_zgy[n_S],
	      skewness_zgy[n_S],kurtosis_zgy[n_S];

	for(k=0;k<n_S;k++)
	{
		if((k+1)!=n_S){	offset = (S_imax[k+1]-S_imax[k])/2;}
		else if( S_imax[k]>200){ offset = 200;}
		start = (int)S_imax[k]-offset; //shift from the pick to the valley
		end = (int)S_imax[k]+offset;
		mean_zgy[k] = calculate_mean(z_gy,start,end);
 	      	calculate_Max_Min_Range(z_ac,start,end,(max_zgy+k),(min_zgy+k),
				(holder+k));
       		calculate_Statistics (z_gy,start,end,mean_zgy[k],(MAD_zgy+k),
	         (holder+k),(holder+k),(skewness_zgy+k),(kurtosis_zgy+k));
	}



	float *features[]={
		mean_yac,max_yac,min_yac,ratio_yac,variance_yac,skewness_yac,
	        min_zac,kurtosis_xgy,
	        max_zgy,min_zgy,MAD_zgy,kurtosis_zgy};

//	float **features[]={features_xac,features_yac};

	nfeatures = 12;
	training_file(fp,train_file_name,n_S,S_imax,features,activityCode,t,
			nfeatures);



	printf("extract_stride_data completed successfuly. Exiting.\n");
	return 0;
}

////////////////////////////



void training_file(FILE *fp,char* train_file_name, int n_S,float *S_imax,
		float **Features,int activityCode,double *t,int nfeatures)
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
//	fprintf(fp,"start\n");

	for (i = 0; i < (n_S); i++) {
		//idx_max = (int) S_imax[i];
		//idx_next = (int) S_imax[i+1];
		//((i+1)!=n_S)? period = t[idx_next]- t[idx_max]: period;
		//fprintf(fp, "%lf",period/10.0);
		
		for(j=0;j<nfeatures;j++)
		{
			fprintf(fp,"%f ",Features[j][i]/10);
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
		float **Features,int activityCode,double *t,int nfeatures)
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
//	fprintf(fp,"start\n");

	for (i = 0; i < (n_S); i++) {
		//idx_max = (int) S_imax[i];
		//idx_next = (int) S_imax[i+1];
		//((i+1)!=n_S)? period = t[idx_next]- t[idx_max]: period;
		//fprintf(fp, "%lf",period/10.0);
		
		for(j=0;j<nfeatures;j++)
		{
			fprintf(fp,"%f ",Features[j][i]/10);
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


