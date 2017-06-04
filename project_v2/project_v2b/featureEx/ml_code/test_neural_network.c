#include <unistd.h>

#include <stdio.h>
#include "floatfann.h"

#include <string.h>
#include <stdlib.h>

int main()
{
    int i;
    int motion;
   
    float max;
    fann_type *calc_out;
    fann_type input[12];
    struct fann *ann;
   
    
    ann = fann_create_from_file("TEST.net");

    FILE *fp;
    char *file_name = "test_data.txt";
    float a,b,c,d,e,f,g,h,l,m,n,o;//,p,q;
    fp = fopen(file_name,"r");
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    
    read = getline(&line, &len, fp); //discard the frist line


    while ((read = getline(&line, &len, fp)) != -1) {
//        sscanf(line,"%f %f %f\n",&a,&b,&c);
        sscanf(line,"%f %f %f %f %f %f %f %f %f %f %f %f\n",&a,&b,&c,&d,&e
			,&f,&g,&h,&l,&m,&n,&o);
     
        input[0] = (float) a;
        input[1] = (float) b;
        input[2] = (float) c;
	input[3] = (float) d;
        input[4] = (float) e;
        input[5] = (float) f;
	input[6] = (float) g;
        input[7] = (float) h;
        input[8] = (float) l;
        input[9] = (float) m;
  	input[10] = (float) n;
   	input[11] = (float) o;
//   	input[12] = (float) p;
//	input[13] = (float) q;


        calc_out = fann_run(ann, input);
        max = -100;
        for (i = 0; i < 5; i++) {
            if (calc_out[i] > max) {
                max = calc_out[i];
                motion = i;
            }
        }

	switch(motion){
	    case 0:
		printf("walking\n");
		break;
            case 1:
		printf("running\n");
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
