// Author: Hang-Hyun Jo (h2jo23@gmail.com) and Tibebe Birhanu (tibebebirhanu@gmail.com)
// Description: main functions

void get_ready(char *folder, char *filename, long *timings, long *iet_sequence, long **tree, long num_event, long isTreeRead){
    long num_iet;

    /* read timings
    timings = vector_long(0, num_event_max);
    num_event = read_timings(timings, num_event_max, filename, folder);
    */

    num_iet = num_event - 1;

    // get interevent time sequence
    get_iet_sequence_from_timings(timings, num_event, iet_sequence);

    // get tree
    if(isTreeRead == 0){ // detect tree after obtaining inter-event time sequence
        get_tree(iet_sequence, num_iet, tree);
        print_tree(folder, filename, tree, num_iet);
        printf("get and print tree done\n");
    }
    else{ // read tree from the file
        read_tree(folder, filename, tree, num_iet);
        printf("read tree done\n");
    }
}


void burst_analysis_option(char *folder, char *filename, long ens, long *timings, long *iet_sequence, long **tree, long num_event, char *options){
    long i, j, k, num_iet, train_Dts[6], num_Dt, iet_max;
    double autocorrel_params[3], logbin_params[3];
    char filename_ens[500];

    num_iet = num_event - 1;
    iet_max = find_max(iet_sequence, num_iet);
    //if(ens >= 0) sprintf(filename_ens, "%s_ens%ld", filename, ens);
    //else sprintf(filename_ens, "%s", filename);

    // kernel 2D
    if(options[1] == '1'){
        printf("get kernel2D\n");
        logbin_params[0] = 5; logbin_params[1] = 0.18;
        get_kernel2DLog(folder, filename, tree, num_iet, logbin_params);
    }

// USING iet_sequence
    // iet distribution
    if(options[2] == '1'){
        printf("get iet distr\n");
        logbin_params[0] = 2; // binStart
        logbin_params[1] = 0.2; // binSize
        logbin_params[2] = 10; // xc = crossover of x (?)
        get_iet_distr(folder, filename, iet_sequence, num_iet, logbin_params);
    }

// train size distribution
    if(options[3] == '1'){
        printf("get train distr\n");
        train_Dts[0] = 600;    train_Dts[1] = 721;
        train_Dts[2] = 751;  train_Dts[3] = 781;
        train_Dts[4] = 901; train_Dts[5] = 931;
        logbin_params[0] = 1; // binStart
        logbin_params[1] = 0.2; // binSize
        logbin_params[2] = 10; // xc = crossover of x (?)
        //num_Dt = 3; // Wikipedia & JUNEC
        //num_Dt = 2; // Twitter
        num_Dt = 6; // heartbeat
        get_train_distr(folder, filename, iet_sequence, num_iet, train_Dts, num_Dt, logbin_params);
    }
// USING tree
    // memory coeff of trains, LC, suscept
    if(options[4] == '1'){
        printf("get percolation\n");
        logbin_params[0] = 5; logbin_params[1] = 0.18;
        get_percolation(folder, filename, tree, num_iet, logbin_params);
    }
    
//calculate the memory coefficient between iets
// memory coeff of iet
    if(options[5] == '1'){
        printf("get memory iet\n");
        get_memory_iet(folder, filename, "", iet_sequence, num_iet, 0); // 1 = ensemble
    }  
    
      
//burstiness calculation for bursty time series
    if(options[6] = '1'){
       printf("get burstiness \n");
       
        num_Dt = iet_max;
       logbin_params[0] = 5; logbin_params[1] = 0.18;
       get_burstiness_burst(folder, filename, iet_sequence, num_iet, num_Dt, logbin_params);
    }  
 
 // memory coefficients b/n bursts based on iet sequence
 if(options[7] == '1'){
       printf("get memory_burst\n");
       
       num_Dt = iet_max;
       logbin_params[0] = 5; logbin_params[1] = 0.18;
       get_memory_burst_model(folder, filename, iet_sequence, num_iet, num_Dt, logbin_params);
   }     
}

void summarize_ensemble_option(char *folder, char *filename, char *options, long ens){
     long i, j, k;
     long train_Dts[6];
     char prefix[100];
     
     if(options[3] == '1'){
        printf("get train distr\n");
        train_Dts[0] = 700;    train_Dts[1] = 800;
        train_Dts[2] = 900;  train_Dts[3] = 1000;
        train_Dts[4] = 1100; train_Dts[5] = 1200;
        for(i = 0; i < 6; i++){
            sprintf(prefix, "train_distr_Dt%ld", train_Dts[i]);
            printf("%s\n", prefix);
            summarize_distributions(folder, filename, prefix, ens);
        
        }
     }  
}

void  find_average(double **curves, long index, long ens, double *avg, double *std){
      long i, j, k, counter;
      double x, y, avg0, std0;
      avg0 = std0 = 0;
      counter = 0;
      for(i = 1; i <= ens; i++){
          x = curves[index][i];
          if(x){
          avg0 += x;
          std0 += x * x;
          counter++;
          }
      }
      if(counter){
        avg0 /= (double)counter;
        std0 = sqrt(std0/(double)counter - avg0 * avg0)/sqrt(counter);
      }
      else{
          avg0 = 0;
          std0 = 0;   
      }
      *avg = avg0;
      *std = std0;
}


void summarize_kernelDiag(char *folder, char *filename, long ens){
     long i, j, k, x_min, x_max;
     double **curves, **curves2, x, y, avg, std;
     char input[2000], input2[2000], output[2000], output2[500], prefix[500], prefix2[500];
     FILE *file_in, *file_in2, *file_out, *file_out2;
     
     x_min = 0;
     x_max = 0;
     strcpy(prefix, "diagonal_kernel");
     strcpy(prefix2, "diagonal_kernel2");
     for(i = 1; i <= ens; i++){
        sprintf(input, "%s%s_%s%ld.txt", folder, prefix, filename, i);
        k = count_line(input);
        if(k > x_max)x_max = k;
     }
     printf("maximum count line = %ld\n", k);
     
     curves = matrix_double(x_min, x_max, 0, ens);
     for(i = x_min; i < x_max;i++){
        for(j = 0; j <= ens; j++){
            curves[i][j] = 0;
        }
     }
     
    curves2 = matrix_double(x_min, x_max, 0, ens);
     for(i = x_min; i < x_max;i++){
        for(j = 0; j <= ens; j++){
            curves2[i][j] = 0;
        }
     }
    
     for(i = 1; i <= ens; i++){
        sprintf(input, "%s%s_%s%ld.txt", folder, prefix, filename, i); 
        file_in = fopen(input, "r");
        j = 0;
        while(fscanf(file_in, "%lf %lf", &x, &y) && !feof(file_in)){
             curves[j][0] = x;
             curves[j][i] = y;
             j++;
        }
       fclose(file_in);
     }
     
     //for the opposite diagonal
     for(i = 1; i <= ens; i++){
        sprintf(input2, "%s%s_%s%ld.txt", folder, prefix2, filename, i); 
        file_in2 = fopen(input2, "r");
        j = 0;
        while(fscanf(file_in2, "%lf %lf", &x, &y) && !feof(file_in2)){
             curves2[j][0] = x;
             curves2[j][i] = y;
             j++;
        }
       fclose(file_in2);
     }
      
     sprintf(output, "%s%s_%s_summary.txt", folder, prefix, filename);
     file_out = fopen(output, "w");
     for(i = 0; i < k; i++){
        find_average(curves, i, ens, &avg, &std);
        if(avg){
        fprintf(file_out, "%lf %.10lf %.10lf\n", curves[i][0], avg, std);
        }
     }
     fclose(file_out);
     
     //for the opposite digonal 
     sprintf(output2, "%s%s_%s_summary.txt", folder, prefix2, filename);
     file_out2 = fopen(output2, "w");
     for(i = 0; i < k; i++){
        find_average(curves2, i, ens, &avg, &std);
        if(avg){
        fprintf(file_out2, "%lf %.10lf %.10lf\n", curves2[i][0], avg, std);
        }
     }
     fclose(file_out2);
     
     free_matrix_double(curves, x_min, x_max, 0, ens);
     free_matrix_double(curves2, x_min, x_max, 0, ens);
}


void summarize_burstiness(char *folder, char *filename, long ens){
     long i, j, k;
     double x, y, *burstiness;
     char prefix[200], output[2000], output2[2000];
     FILE *file_out, *file_out2;
     
     burstiness = vector_double(1, ens);
     
     sprintf(prefix, "burstiness");
     for(i = 1; i <= ens; i++){
        sprintf(output, "%s%s_%s%ld.txt", folder, prefix, filename, i);
        file_out = fopen(output, "r");
        while(fscanf(file_out,"%lf", &x) == 1){
              burstiness[i] = x;
              printf("%lf\n", burstiness[i]);
        }
        fclose(file_out);    
     }
     sprintf(output2, "%s%s_%s_summary.txt", folder, prefix, filename);
     file_out2 = fopen(output2, "w");
     for (i = 1; i <= ens; i++){
          fprintf(file_out2, "%lf\n", burstiness[i]);
     }
     fclose(file_out2);
     
     
     
     free_vector_double(burstiness, 1, ens);

}

void summarize_memory(char *folder, char *filename, long ens){
     long i, j, k;
     double x, y, *memory;
     char prefix[200], output[2000], output2[2000];
     FILE *file_out, *file_out2;
     
     memory = vector_double(1, ens);
     
     sprintf(prefix, "memory_iet");
     for(i = 1; i <= ens; i++){
        sprintf(output, "%s%s_%s%ld.txt", folder, prefix, filename, i);
        file_out = fopen(output, "r");
        while(fscanf(file_out,"%lf", &x) == 1){
              memory[i] = x;
              printf("%lf\n", memory[i]);
        }
        fclose(file_out);    
     }
     sprintf(output2, "%s%s_%s_summary.txt", folder, prefix, filename);
     file_out2 = fopen(output2, "w");
     for (i = 1; i <= ens; i++){
          fprintf(file_out2, "%lf\n", memory[i]);
     }
     fclose(file_out2);
     
     
     
     free_vector_double(memory, 1, ens);

}

void summarize_dt_interval(char *folder, char *filename, long ens){
     long i, j, k;
     double x, y, *dt_interval;
     char prefix[200], input[2000], output[2000];
     FILE *file_in, *file_out;
     
     dt_interval = vector_double(1, ens);
     
     sprintf(prefix, "Dt_interval");
     for(i = 1; i <= ens; i++){
        sprintf(input, "%s%s_%s%ld.txt", folder, prefix, filename, i);
        file_in = fopen(input, "r");
        while(fscanf(file_in,"%lf", &x) == 1){
              dt_interval[i] = x;
              printf("%lf\n", dt_interval[i]);
        }
        fclose(file_in);    
     }
     sprintf(output, "%s%s_%s_summary.txt", folder, prefix, filename);
     file_out = fopen(output, "w");
     for (i = 1; i <= ens; i++){
          fprintf(file_out, "%lf %lf\n",3., dt_interval[i]);
     }
     fclose(file_out);
     
     
     
     free_vector_double(dt_interval, 1, ens);

}

//summarize the smoothed timescale interval
void summarize_sdt_interval(char *folder, char *filename, long ens){
     long i, j, k;
     double x, y, *dt_interval;
     char prefix[200], input[2000], output[2000];
     FILE *file_in, *file_out;
     
     dt_interval = vector_double(1, ens);
     
     sprintf(prefix, "smoothed_Dt_interval");
     for(i = 1; i <= ens; i++){
        sprintf(input, "%s%s_%s%ld.txt", folder, prefix, filename, i);
        file_in = fopen(input, "r");
        while(fscanf(file_in,"%lf", &x) == 1){
              dt_interval[i] = x;
              printf("%lf\n", dt_interval[i]);
        }
        fclose(file_in);    
     }
     sprintf(output, "%s%s_%s_summary.txt", folder, prefix, filename);
     file_out = fopen(output, "w");
     for (i = 1; i <= ens; i++){
          fprintf(file_out, "%.1lf %lf\n",3., dt_interval[i]);
     }
     fclose(file_out);
     
     
     
     free_vector_double(dt_interval, 1, ens);

}

void summarize_option(char *folder, char *filename, char *options, long ens){
     long i, j, k;
     
     // summarize the diagonal kernel 
     if(options[1] = '1'){
        summarize_kernelDiag(folder, filename, ens);
     }
     //collect the values of burstiness parameter value
     if(options[2] = '1'){
        summarize_burstiness(folder, filename, ens);
     }
     
     //collect the values of memory coefficient
     if(options[5] = '1'){
       summarize_memory(folder, filename, ens); 
     }
     
     if(options[6] = '1'){
       summarize_dt_interval(folder, filename, ens);
       summarize_sdt_interval(folder, filename, ens);
     }
}

void burst_analysis_option_ens(char *folder, char *filename, long ens, long *timings, long *iet_sequence, long **tree, long num_event, long num_event_max, char *options, long isTreeRead){

     long i, j, k, num_iet;
     char filename_ens[2000];
    
         for(i = 1; i <= ens; i++){
            sprintf(filename_ens, "%s%ld", filename,i);
            k = read_timings(folder, filename_ens, timings, num_event_max);
            printf("number of event %ld\n", k);
            get_ready(folder, filename_ens, timings, iet_sequence, tree, k, isTreeRead);
            burst_analysis_option(folder, filename_ens, ens, timings, iet_sequence, tree, k, options);    
         }
      
         summarize_option(folder, filename, options, ens);
}    
