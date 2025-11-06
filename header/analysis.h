// Author: Hang-Hyun Jo (h2jo23@gmail.com) and Tibebe Birhanu (tibebebirhanu@gmail.com)
// Description: functions for some basic operations

// get distribution
#include <stdbool.h>
void get_distr_long(long *sequence, long num_item, long *distr, long x_min, long x_max){
    long i, x;

    for(x = x_min; x <= x_max; x ++) distr[x] = 0;

    for(i = 0; i < num_item; i ++){
        x = sequence[i];
        distr[x] ++;
    }
}

void get_distr_double(double *sequence, long num_item, double resolution, long *distr, long x_min, long x_max){ // x is in [x_min, x_max]
    long i, x;

    for(x = x_min; x <= x_max; x ++) distr[x] = 0;

    for(i = 0; i < num_item; i ++){
        x = (long)(sequence[i] / resolution + 0.5);
        if(x < x_min || x > x_max) printf("out of range\n");
        distr[x] ++;
    }
}

long find_max(long *sequence, long num_item){
    long i, x, x_max;

    x_max = -1;
    for(i = 0; i < num_item; i ++){
        if(sequence[i] > x_max)
            x_max = sequence[i];
    }

    return x_max;
}
//find the maximum value fro 2D array
double find_max_2D(double **sequence, long num_Dt){
    long i, x;
    double  x_max;

    x_max = -1;
    for(i = 0; i < num_Dt; i ++){
        if(sequence[i][1] > x_max)
            x_max = sequence[i][1];
    }

    return x_max;
}
//find the corresponding values from the maximum value(the second column)
double find_max_2D_2(double **sequence, long num_Dt){
    long i, x;
    double  x_max, y_value;

    x_max = -1;
    for(i = 0; i < num_Dt; i ++){
        // to find the first maximum: exit the for loop when you get the first max
        if(sequence[i][1] < x_max)break;
        if(sequence[i][1] > x_max){
            x_max = sequence[i][1];
            y_value = sequence[i][0];
            printf("max_burst%lf max_dt%lf\n", x_max, y_value);
            }
    }

    return y_value;
}


long find_min(long *sequence, long num_item){
    long i, x, x_min;

    x_min = 1e8;
    for(i = 0; i < num_item; i ++){
        if(sequence[i] < x_min)
            x_min = sequence[i];
    }

    return x_min;
}

// for event sequence analysis
long read_timings(char *folder, char *filename, long *timings, long n_max){
    long i, timing, timing0; 
    char output[500]; 
    FILE *file_in; 

    for(i = 0; i < n_max; i ++) timings[i] = 0;

    sprintf(output, "%stimings_%s.txt", folder, filename);
    file_in = fopen(output, "r"); 
    i = 0;
    timing0 = -1;
    while(fscanf(file_in, "%ld", &timing) && !feof(file_in)){
        if(timing == timing0) continue;
        timings[i] = timing;
        i ++;
        if(i == n_max) printf("data size exceeding the limit in read_timings()\n");
        timing0 = timing;
    }
    fclose(file_in); 

    return i;
}

// exchange one form to another
long get_iet_sequence_from_timings(long *timings, long n, long *iet_sequence){
    long i, timing0, timing, iet, iet_max;

    for(i = 0; i < n - 1; i ++) iet_sequence[i] = 0;

    iet_max = -1;
    timing0 = -1;
    for(i = 0; i < n; i ++){
        timing = timings[i];
        if(timing0 >= 0){
            iet = timing - timing0;
            iet_sequence[i - 1] = iet;
            if(iet > iet_max) iet_max = iet;
        }
        timing0 = timing;
    }

    return iet_max;
}

void get_timing_sequence_from_iet_sequence(long *iet_sequence, long num_iet, long *timing_sequence){
    long i, timing;

    for(i = 0; i <= num_iet; i ++) timing_sequence[i] = 0;

    timing = 0;
    for(i = 0; i < num_iet; i ++){
        timing += iet_sequence[i];
        timing_sequence[i + 1] = timing;
    }
}

void get_iet_sequence_from_tree_direct(long **tree, long *iet_sequence, long num_iet){
    long i, loc;

    for(i = 0; i < num_iet; i ++){
        loc = tree[i][0];
        iet_sequence[loc] = tree[i][1];
    }
}

void get_iet_sequence_from_tree(long **tree, long root_rank, long root_loc, long *iet_sequence, long num_iet){
    long rank_left, rank_right, train0, train1, loc;

    rank_left = tree[root_rank][2];
    rank_right = tree[root_rank][3];

    if(rank_left < 0 && rank_right < 0) return;

    if(rank_left >= 0){
        train1 = tree[rank_left][5]; // train size of right branch of "rank_left"
        loc = root_loc - train1;
        iet_sequence[loc] = tree[rank_left][1];
        get_iet_sequence_from_tree(tree, rank_left, loc, iet_sequence, num_iet);
    }
    if(rank_right >= 0){
        train0 = tree[rank_right][4]; // train size of left branch of "rank_right"
        loc = root_loc + train0;
        iet_sequence[loc] = tree[rank_right][1];
        get_iet_sequence_from_tree(tree, rank_right, loc, iet_sequence, num_iet);
    }
}
// iet distribution
void get_iet_distr(char *folder, char *filename, long *iet_sequence, long num_iet, double *logbin_params){
    long i, iet_min, iet_max, *iet_distr, index, numBin, num_event, mode_iet, mode_iet2;
    double ** distr, num, y, xc, binSize, max_distr, max_distr2;
    double avg, std, t1, t2, r, A_n;
    char prefix[500], prefix2[500], output[2000],output2[2000], output3[2000], output4[2000];
    FILE *file_out, *file_out2, *file_out3, *file_out4;
    
    num_event = num_iet + 1;
    t1 = sqrt(num_event + 1);
    t2 = sqrt(num_event - 1);
    xc = 0;
    binSize = 30;
    iet_min = 0;
    iet_max = find_max(iet_sequence, num_iet);

    iet_distr = vector_long(iet_min, iet_max);
    get_distr_long(iet_sequence, num_iet, iet_distr, iet_min, iet_max);

    // print
    sprintf(prefix, "iet_distr");
    sprintf(prefix2, "burstiness");
    //print_distr_noTotal(folder, filename, prefix, iet_distr, iet_min, iet_max);
    distr = matrix_double(0, iet_max, 0, 1);
    for(i = 0; i <= iet_max; i ++) distr[i][0] = distr[i][1] = 0;

    index = 0;
    num = 0;
    for(i = 0; i <= iet_max; i ++){
        y = iet_distr[i];
       if(y > 0){
            distr[index][0] = i;
            distr[index][1] = y;
            num += y;
            index ++;
        }
    }
    
    //calculate the median in iet distribution or from the iet sequence
    max_distr = -1;
    for(i = 0; i < index; i++){
       if(distr[i][1] > max_distr){
          max_distr = distr[i][1];
          mode_iet = distr[1][0];
       }
    
    }
    sprintf(output3, "%smode_without_bin_%s.txt", folder, filename);
    file_out3 = fopen(output3, "w");
    fprintf(file_out3, "%lf %ld", max_distr, mode_iet);
    fclose(file_out3);
    
	for(i = 0; i < index; i ++) distr[i][1] /= num;
	avg = 0;
	std = 0;
	for(i = 0; i < index; i++){
	   avg += distr[i][0] * distr[i][1];
	   std += distr[i][0] * distr[i][0] * distr[i][1];
	}
        std -= avg*avg;
        r = sqrt(std)/avg;
        A_n = (t1 * r - t2)/((t1 - 2)*r + t2);
        sprintf(output2, "%s%s_%s.txt", folder, prefix2, filename);
        file_out2 = fopen(output2,"w");
        fprintf(file_out2, "%lf", A_n);
        fclose(file_out2);
	numBin = LinBin_number(distr, 0, index-1, xc, binSize);
	//numBin = LogBin_period_mix_nx(distr, 0, index - 1, binStart, xc, binSize);

	sprintf(output, "%s%s_%s.txt", folder, prefix, filename);
    file_out = fopen(output, "w");
    for(i = 0; i < numBin; i ++){
        fprintf(file_out, "%.15lf %.15lf\n", distr[i][0], distr[i][1]);
    }
    fclose(file_out);
    max_distr2 = -1;
    for(i = 0; i < numBin; i++){
       if(distr[i][1] > max_distr2){
          max_distr2 = distr[i][1];
          mode_iet2 = distr[i][0];
       }
    
    }
    sprintf(output4, "%smode_after_bin_%s.txt", folder, filename);
    file_out4 = fopen(output4, "w");
    fprintf(file_out4, "%lf %ld", max_distr2, mode_iet2);
    fclose(file_out4);
	free_matrix_double(distr, 0, iet_max, 0, 1);

    // print logbinned distr
    get_logbin_period_noAvg(folder, filename, prefix, iet_distr, iet_max, logbin_params);

    free_vector_long(iet_distr, iet_min, iet_max);
}

// train size distributions
long get_train_sequence_given_Dt(long *iet_sequence, long num_iet, long *train_sequence, long Dt){
    long i, iet, train, num_train, num_event;

    num_event = num_iet + 1;

    for(i = 0; i < num_event; i ++) train_sequence[i] = 0;

    train = 1;
    num_train = 0;
    for(i = 0; i < num_iet; i ++){
        iet = iet_sequence[i];
        if(iet <= Dt) train ++;
        else{
            train_sequence[num_train] = train;
            num_train ++;
            train = 1;
        }
    }

    return num_train;
}

void get_train_distr(char *folder, char *filename, long *iet_sequence, long num_iet, long *Dts, long num_Dt, double *logbin_params){
    long i, Dt, num_train, num_event, train_min, train_max;
    long *train_sequence, *train_distr;
    char prefix[500];

    num_event = num_iet + 1;
    
    train_sequence = vector_long(0, num_event);
    train_distr = vector_long(0, num_event);

    for(i = 0; i < num_Dt; i ++){
        Dt = Dts[i];
        printf("Dt=%ld\n", Dt);
        num_train = get_train_sequence_given_Dt(iet_sequence, num_iet, train_sequence, Dt);
        //train_min = find_min(train_sequence, num_train);
        train_min = 0;
        train_max = find_max(train_sequence, num_train);
        get_distr_long(train_sequence, num_train, train_distr, train_min, train_max);

        // print
        sprintf(prefix, "train_distr_Dt%ld", Dt);
        print_distr_noTotal(folder, filename, prefix, train_distr, train_min, train_max);
        
        // print logbinned distr
        get_logbin_period(folder, filename, prefix, train_distr, train_max, logbin_params);
    }

    free_vector_long(train_sequence, 0, num_event);
    free_vector_long(train_distr, 0, num_event);
}


// memory coeff of iet (but also for any sequence)
double calculate_memory_coeff(long *iet_sequence, long num_iet){ 
    long i, iet; 
    double iet_sum, iet2_sum, iet0, iet1, mean1, mean2, std1, std2, ietiet, memory;

    iet_sum = iet2_sum = ietiet = 0;
    for(i = 0; i < num_iet; i ++){
        iet = iet_sequence[i];
        iet_sum += iet;
        iet2_sum += iet * iet;
        if(i > 0) ietiet += iet * iet_sequence[i - 1];
    }
    iet0 = iet_sequence[0];
    iet1 = iet_sequence[num_iet - 1];

    mean1 = (iet_sum - iet1) / (double)(num_iet - 1);
    mean2 = (iet_sum - iet0) / (double)(num_iet - 1);
    std1 = sqrt((iet2_sum - iet1 * iet1) / (double)(num_iet - 1) - mean1 * mean1);
    std2 = sqrt((iet2_sum - iet0 * iet0) / (double)(num_iet - 1) - mean2 * mean2);
    ietiet /= (double)(num_iet - 1); //ietiet /= (double)(n - 2); <- wrong?

    memory = (ietiet - mean1 * mean2) / std1 / std2;

    return memory;
}

// memory coeff of iet (but also for any sequence)
void get_memory_iet(char *folder, char *filename, char *rrm, long *iet_sequence, long num_iet, long option){ 
    long i, iet; 
    double memory;
    char prefix[500], output[2000];
    FILE *file_out;

    memory = calculate_memory_coeff(iet_sequence, num_iet);

    sprintf(prefix, "memory_iet");

    if(option == 0){ // original value
        sprintf(output, "%s%s_%s.txt", folder, prefix, filename);
        file_out = fopen(output, "w");
        fprintf(file_out, "%lf\n", memory);
        fclose(file_out);
    }
    else if(option == 1){ // ensemble
        sprintf(output, "%s%s_%s%s_ens.txt", folder, prefix, filename, rrm);
        file_out = fopen(output, "a"); // add not write
        fprintf(file_out, "%lf\n", memory);
        fclose(file_out);
    }
}

double measure_burst(long num_train, long *train_sequence){
    long i, j, train, next_train;
    double train_sum, train2_sum, train0, train1, traintrain, mean, std, r, bursti;
    
        train_sum = train2_sum = traintrain = 0;
        
        for(j = 0; j < num_train; j++){
        train = train_sequence[j];
        train_sum += train;
        train2_sum += train*train; 
        //printf("j %ld train %ld \n", j, train);
        }
        mean = (train_sum) / (double)(num_train);
        std = sqrt((train2_sum) / (double)(num_train) - mean * mean);
        r = std/mean;
        bursti = (sqrt(num_train + 1)*r - sqrt(num_train - 1))/((sqrt(num_train + 1) - 2.0)*r + sqrt(num_train - 1));
        
       return bursti ;
}

double measure_memory(long num_train, long *train_sequence){
    long i, j, train, next_train;
    double train_sum, train2_sum, train0, train1, traintrain, mean1, mean2, std1, std2, mem;
    
        train_sum = train2_sum = traintrain = 0;
       
        for(j = 0; j < num_train; j++){
            train = train_sequence[j];
            next_train = train_sequence[j+1];
            train_sum += train;
            train2_sum += train*train; 
            traintrain += train * next_train;
           // printf("j %ld train %ld next_train %ld\n", j, train, next_train);
        }
        
        train0 = train_sequence[0];
        train1 = train_sequence[num_train-1];
        mean1 = (train_sum - train1) / (double)(num_train - 1);
        mean2 = (train_sum - train0) / (double)(num_train - 1);
        std1 = sqrt((train2_sum - train1 * train1) / (double)(num_train - 1) - mean1 * mean1);
        std2 = sqrt((train2_sum - train0 * train0) / (double)(num_train - 1) - mean2 * mean2);
        traintrain /= (double)(num_train - 1);
        mem = (traintrain - mean1 * mean2) / std1 / std2;
        
       return mem ;
}

double calculate_min(double c_min, double **burstiness, long k) {
    long i;
    double x_min = -1000;
    double t_min = -1; // Initialize t_min to an invalid value
    for (i = 0; i < k; i++) {
        if (burstiness[i][1] > x_min) {
            x_min = burstiness[i][1];
        }
        if (x_min > c_min) {
            t_min = burstiness[i][0]; // Record the time scale
            printf("The min c and timescale min: %lf %lf\n", x_min, t_min);
            break; // Exit after finding the first valid minimum
        }
    }
    return t_min;
}
double calculate_max(double c_max, double **burstiness, long k) {
    long i;
    double x_max = -1000;
    double t_max = -1; // Initialize t_max to an invalid value
    for (i = 0; i < k; i++) {
        if (burstiness[i][1] > x_max) {
            x_max = burstiness[i][1]; // Correctly update x_max
        }
        if (x_max > c_max) {
            t_max = burstiness[i - 1][0]; // Record the time scale of the last valid maximum
            x_max = burstiness[i - 1][1]; // Update x_max to the last valid maximum
            printf("The max c and timescale max: %lf %lf\n", x_max, t_max);
            break; // Exit after finding the last valid maximum
        }
    }
    return t_max;
}

//moving average
void moving_average_smooth(double *data, double *smoothed, long n, long window_size){
     long i, j, k, half_window, count;
     double sum; 
     
     half_window = window_size/2;
     for(i = 0; i < n; i++){
        sum = 0;
        count = 0;
        for(j = i - half_window; j <= i + half_window; j++){
           if(j >= 0 && j < n){
              sum += data[j];
              count++;
           }
        }   
        
       smoothed[i] = sum/count;

      }  
}

void get_burstiness_burst(char *folder, char *filename, long *iet_sequence, long num_iet, long num_Dt, double *logbin_params) {
    long i, j, k, Dt, num_train, num_event, train_min, train_max, train, next_train, count, window_size;
    long iet_min, m2, Dt_m2, m3, Dt_m3;
    double train_sum, train2_sum, train0, train1, traintrain, mean1, mean2, std1, std2, burst, max_burst, Dt_max;
    double **burstiness, burst_m2, burst_m3, c_min, c_max, dt_min, dt_max, dt_interval, *x_values, *y_values, *smoothed_y, **xy_values;
    double sc_min, sc_max, sdt_min, sdt_max, sdt_interval;
    long *train_sequence;
    char prefix[500], prefix2[500], prefix3[500], prefix4[500], output1[2000], output2[2000], output3[2000], output4[2000];
    FILE *file1_out, *file2_out, *file3_out, *file4_out;
    num_event = num_iet + 1;

    train_sequence = vector_long(0, num_event);
    burstiness = matrix_double(0, num_Dt, 0, 2);
    xy_values = matrix_double(0, num_Dt, 0, 1);
    x_values = vector_double(0, num_Dt);
    y_values = vector_double(0, num_Dt);
    
    for(i = 0; i <= num_Dt; i++) burstiness[i][0] = burstiness[i][1] = burstiness[i][2] = 0;
    for(i = 0; i <= num_Dt; i++) xy_values[i][0] = xy_values[i][1] = 0;
    for(i = 0; i <= num_Dt; i++) x_values[i] = y_values[i] = 0;
    
    iet_min = find_min(iet_sequence, num_iet);
    printf("minimum iet %ld \n", iet_min);
    max_burst = 0; 
    k = 0;   

    // Main loop to calculate burstiness for increments of 10
    for(i = iet_min - 10; i < num_Dt; i += 10) {
        Dt = i;
        num_train = get_train_sequence_given_Dt(iet_sequence, num_iet, train_sequence, Dt);
        if (num_train >= 2) {
            burst = measure_burst(num_train, train_sequence);
            burstiness[k][0] = Dt;
            burstiness[k][1] = burst;
            burstiness[k][2] = num_train;
            k++;
        }
    }
    
    // Write the unique values to the output file
    strcpy(prefix, "burstiness_burst");
    sprintf(output1, "%s%s_%s.txt", folder, prefix, filename); 
    file1_out = fopen(output1, "w");
    count = 0;
    for(i = 0; i < k; i++) {
        if (burstiness[i][0] && burstiness[i][1])
        x_values[count] = burstiness[i][0];
        y_values[count] = burstiness[i][1];
        fprintf(file1_out, "%lf %lf %lf\n", burstiness[i][0], burstiness[i][1], burstiness[i][2]);
        count++;
    }
    fclose(file1_out);
    
    c_min = -0.8;
    c_max = 0.5;
    dt_min = calculate_min(c_min, burstiness, k);
    dt_max = calculate_max(c_max, burstiness, k);
    dt_interval = dt_max - dt_min;
    strcpy(prefix2, "Dt_interval");
    sprintf(output2, "%s%s_%s.txt", folder, prefix2, filename);
    file2_out = fopen(output2, "w");
    fprintf(file2_out, "%lf", dt_interval);
    
    fclose(file2_out);
    
    //smoothing the curve and calculate the timescale interval
    window_size = 9;
    
     smoothed_y = vector_double(0, count);
    
      
     moving_average_smooth(y_values, smoothed_y, count, window_size);
      
     strcpy(prefix3, "smoothed_burstiness_burst");
     sprintf(output3, "%s%s_%s.txt", folder, prefix3, filename);
     file3_out = fopen(output3, "w");
       for(i = 0; i < count; i++){
          xy_values[i][0] = x_values[i];
          xy_values[i][1] = smoothed_y[i];
          fprintf(file3_out, "%.10lf %.10lf\n", x_values[i], smoothed_y[i]);        
      
       }
    
       fclose(file3_out);
     
     
    sc_min = -0.8;
    sc_max = 0.45;
    sdt_min = calculate_min(sc_min, xy_values, count);
    sdt_max = calculate_max(sc_max, xy_values, count);
    sdt_interval = sdt_max - sdt_min;
    strcpy(prefix4, "smoothed_Dt_interval");
    sprintf(output4, "%s%s_%s.txt", folder, prefix4, filename);
    file4_out = fopen(output4, "w");
    fprintf(file4_out, "%lf", sdt_interval);
    
    fclose(file4_out);
     
     
    free_vector_long(train_sequence, 0, num_event);
    free_matrix_double(burstiness, 0, num_Dt, 0, 2);
    free_matrix_double(xy_values, 0, num_Dt, 0, 1);
    free_vector_double(x_values, 0, num_Dt);
    free_vector_double(y_values, 0, num_Dt);
    free_vector_double(smoothed_y, 0, count);
}


// memory coefficients b/n bursts based on iet sequence
void get_memory_burst_model(char *folder, char *filename, long *iet_sequence, long num_iet, long num_Dt, double *logbin_params){
    long i, j, k, Dt, num_train, num_event, train_min, train_max, train, next_train, m3, iet_min;
    double train_sum, train2_sum, train0, train1, traintrain, mean1, mean2, std1, std2, memory, memory_m3;
    double **memories;
    long *train_sequence;
    char prefix[500], output1[2000];
    FILE *file1_out;
    num_event = num_iet + 1;

    train_sequence = vector_long(0, num_event);
    memories = matrix_double(0, num_Dt, 0, 1);
    
    for(i = 0; i <= num_Dt; i++)memories[i][0] = memories[i][1] = 0;
         k = 0;
         iet_min = find_min(iet_sequence, num_iet);
    for(i = iet_min - 10; i < num_Dt; i += 10){
        Dt = i;//Dts[i];
        //printf("Dt=%ld\n", Dt);
        num_train = get_train_sequence_given_Dt(iet_sequence, num_iet, train_sequence, Dt);
        printf("Dt %ld num train %ld \n", Dt, num_train);
              if (num_train >= 3){
                 memory = measure_memory(num_train, train_sequence);
       
                 memories[k][0] = Dt;
                 memories[k][1] = memory;
                 k++;
                 printf("memory is %lf\n", memory);
              }
      
      }
      
      
    strcpy(prefix, "burst_memory");
    sprintf(output1, "%s%s_%s.txt", folder, prefix, filename); 
    file1_out = fopen(output1, "w");
    for(i = 0; i < k; i ++){
        fprintf(file1_out, "%lf %lf\n", memories[i][0], memories[i][1]);
    }
      fclose(file1_out);
    //get_logbin_number_direct(folder, filename, prefix, memories, num_Dt, logbin_params);
  
    free_vector_long(train_sequence, 0, num_event);
    free_matrix_double(memories, 0, num_Dt, 0, 1);
   
} 
