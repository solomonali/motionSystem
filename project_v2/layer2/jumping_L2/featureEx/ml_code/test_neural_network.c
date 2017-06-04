#include <unistd.h>

#include <stdio.h>
#include "floatfann.h"

#include <string.h>
#include <stdlib.h>

int main()
{
    int i,k,nf=1;
    int motion;
    int nOutputs = 2; 
    float max;
    fann_type *calc_out;
    fann_type input[nf];
    struct fann *ann;
   
    
    ann = fann_create_from_file("TEST.net");

    FILE *fp;
    char *file_name = "test_data.txt";
    float f[nf];
    fp = fopen(file_name,"r");
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    
    read = getline(&line, &len, fp); //discard the frist line


    while ((read = getline(&line, &len, fp)) != -1) {
//        sscanf(line,"%f %f %f\n",&a,&b,&c);
        sscanf(line,"%f\n",&f[0]);
     
      	for(k=0;k<nf;k++)
	{
 		input[k] = (float) f[k];
	}

        calc_out = fann_run(ann, input);
        max = -100;
        for (i = 0; i < nOutputs; i++) {
            if (calc_out[i] > max) {
                max = calc_out[i];
                motion = i;
            }
        }

	switch(motion){
	    case 0:
		printf("running 1\n");
		break;
            case 1:
		printf("running 2\n");
		break;
	    case 2:
		printf("jumping\n");
		break;
	    case 3:
		printf("stair ascent\n");
		break;
	    case 4:
		printf("stair descent\n");
		break;
	    default:
		printf("error in motion\n");
		
	}
	sleep(1);
        read = getline(&line, &len, fp); //to discard every other line
    }
    fclose(fp);
    fann_destroy(ann);
    return 0;
}
