void find_avg_std(double **curves, long index, long ENS, double *avg, double *std){
    long i, ens, index_x;
    double x, avg0, std0;

    avg0 = std0 = 0;
    for(ens = 0; ens < ENS; ens ++){
        x = curves[index][ens + 2];
        avg0 += x;
        std0 += x * x;
    }
    avg0 /= (double)ENS;
    std0 = sqrt(std0 / (double)ENS - avg0 * avg0);

    *avg = avg0;
    *std = std0;
}


// summarize distributions
void summarize_distributions(char *folder, char *filename, char *prefix, long ENS){
    long i, j, k, x_min, x_max, ens, x, y;
    long *distr;
    double logbin_params[3];
    char input[500], output[500];
    FILE *file_in, *file_out;

    x_min = 0;
    x_max = 1000000;
    for(ens = 1; ens <= ENS; ens ++){
        sprintf(input, "%s%s_%s%ld.txt", folder, prefix, filename, ens);
        k = count_line(input);
        if(k > x_max) x_max = k;
    }
    printf("max line count=%ld\n", k);
    distr = vector_long(x_min, x_max);

    for(i = x_min; i <= x_max; i ++){
        distr[i] = 0;
    }

    for(ens = 1; ens <= ENS; ens ++){
        sprintf(input, "%s%s_%s%ld.txt", folder, prefix, filename, ens);
        file_in = fopen(input, "r");
        while(fscanf(file_in, "%ld %ld", &x, &y) && !feof(file_in)){
            distr[x] += y;
        }
        fclose(file_in);
    }

    sprintf(output, "%s%s_%s_summary.txt", folder, prefix, filename);
    file_out = fopen(output, "w");
    for(i = x_min; i <= x_max; i ++){
        y = distr[i];
        if(y) fprintf(file_out, "%ld %ld\n", i, y);
    }
    fclose(file_out);

    // print logbinned distr
    logbin_params[0] = 1; // binStart
    logbin_params[1] = 0.2; // binSize
    logbin_params[2] = 10; // xc = crossover of x (?)
    get_logbin_period(folder, filename, prefix, distr, x_max, logbin_params);

    free_vector_long(distr, x_min, x_max);
}

