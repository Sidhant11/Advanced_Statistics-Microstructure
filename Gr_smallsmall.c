/* program for finding the average pair velocity in a many body suspension */

# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <time.h>
# define cy 302 // Number of cycles averaging throguh
# define g 2000 // Total number of particles
# define g1 5 //Number of big particles
# define s cy*g


main() {
    int i, j, k, m, co;
    float re;
    int rd, td, pd;
    int st = 200;         //Start cycle to behin averaging
    int td_s;              /* tangential g(r)=g(-r) bin population */
    int nd, nt;
    float norm_f;          /* Normalization factor in each cycle*/
    float volume_f;
    float xm, ym, zm;
    float counter, p_counter;
    float d_r, d_t, alfa;  /* Gr oriented*/
    float rad_c;
    double pi;
    float rm, vm;
    float er1, er2, er3, et1, et2, et3, ep1, ep2, ep3;
    float erm, etm, epm;
    float r1, r2, r3;
    float vx, vy, vz;
    float vr, vt, vp;
    float rp1, rp2, rp3, gammadot;
    float elapsed;
    clock_t start, end;
    float phi, tetha;
    float phi_s, tetha_s;  /* Symmetry imposing angles */
    float cr1, cr2;
    double **M;
    float N[g][5] = {{0.0}};
    float **HP;
    float **HPG;
    float **MP;
    FILE *infile_ptr;
    FILE *outfile_ptr;
    
    start = clock();
    
    
    pi= acos(-1.0);
    d_r = 0.01;  //Radial resolution
    nd = 60;    
    nt = 180;    //Maximum angle considered - Mirror statistics employed
    rad_c = 8.0;  //Maximum normalized distance considered
    d_t = (2*pi)/nt;
    alfa = 1.1;
    rp1 = 1.0;
    rp2 = 1.0;
    rp3 = 0.5*(rp1 + rp2);
    gammadot = 1.2312e-04;
    xm = 33.7937;
    ym = 33.7937;
    zm = 33.7937;
    cr1 = (pi/2.0)-(12.5*pi/180.0);
    cr2 = (pi/2.0)+(12.5*pi/180.0);
    
    
    /* Allocate M and Input M */
    M = (double **) calloc(s, sizeof(double *));
    for (i = 0; i < s; i++)
        M[i] = (double *)calloc(8, sizeof(double));
    
    infile_ptr = fopen("T.txt","r");
    
    if (infile_ptr == NULL)
    {   printf("\n can not open the file");
        exit(1); }
    
    for (i = 0; i < s; i++)
        fscanf(infile_ptr, "%le %le %le %le %le %le %le %le", &M[i][0], &M[i][1], &M[i][2], &M[i][3], &M[i][4], &M[i][5], &M[i][6], &M[i][7]);
    
    fclose(infile_ptr);
    MP = (float **) calloc(nd, sizeof(float *));           /* Loop histogram matrix */
    for (i = 0; i < nd; i++)
        MP[i] = (float *)calloc(nt, sizeof(float));
    
    HP = (float **) calloc(nd, sizeof(float *));          /* middle histogram matrix */
    for (i = 0; i < nd; i++)
        HP[i] = (float *)calloc(nt, sizeof(float));
    
    HPG = (float **) calloc(nd, sizeof(float *));          /* Grand histogram matrix */
    for (i = 0; i < nd; i++)
        HPG[i] = (float *)calloc(nt, sizeof(float));
    /*
    for (i = 0; i < s; i++)
    {
        
        M[i][2] = (1.0/rp2)*M[i][2];
        M[i][3] = (1.0/rp2)*M[i][3];
        M[i][4] = (1.0/rp2)*M[i][4];
    }
    
    xm /= rp2;
    ym /= rp2;
    zm /= rp2;
    */
    for (co = st; co < cy; co ++)
    {                /* #M */
        
        //printf("\n the cycle number is = %d ", co);
        
        
        /* store the chunk of data needed for a cylce repeat */
        for (i = 0; i < g; i++)
        {
            N[i][0] = M[i+(co)*g][0];
            N[i][1] = M[i+(co)*g][1];
            N[i][2] = M[i+(co)*g][2];
            N[i][4] = M[i+(co)*g][3];
            N[i][3] = M[i+(co)*g][4];
        }
        
        /*Loop over a configuration*/
        
        for (i = 0; i < g ; i++ )
        {                /* #1 */
            
            if (N[i][1]==1.0)
            {   /* #2 */
                p_counter += 1.0;
                
                for (j = 0; j < g; j++)
                {              /* #3 */
                    if (i != j && N[j][1]==1.0)
                    {                             /* #4 */
                        
                        r1 = N[j][2] - N[i][2];
                        r2 = N[j][3] - N[i][3];
                        r3 = N[j][4] - N[i][4];
                        
                        /*periodic interaction correction*/
                        
                        if (fabs(r1) > xm/2.0)
                        {
                            if (r1 > 0.0)
                            {
                                r1 = -1*(xm - fabs(r1));
                            }
                            else
                            {
                                r1 = xm - fabs(r1);
                            }
                            
                        }
                        
                        if (fabs(r3) > zm/2.0)
                        {
                            if (r3 > 0.0)
                            {
                                r3 = -1.0*(zm - fabs(r3));
                            }
                            else
                            {
                                r3 = zm - fabs(r3);
                            }
                        }
                        
                        if (fabs(r2) > ym/2.0)
                        {
                            if (r2 > 0.0)
                            {
                                r2 = -1.0*(ym - fabs(r2));
                            }
                            else
                            {
                                r2 = ym - fabs(r2);
                            }
                        }

                        rm = sqrt(r1*r1 + r2*r2 + r3*r3);
                        
                        if (rm <= rad_c)
                        {                            /* #5 */
                          //  printf("\n interparticle distance %f ", rm);
                            counter += 1.0;   /* It accumulates the total number of particles for finding density */
                            
                            phi = atan2(r2,r1);
                            tetha = acos(r3/rm);
                            
                            if (tetha > cr1 && tetha < cr2)
                            {  /* #6 */
                                /* Imposing g(r) = g(-r) */
                                phi_s = atan2(-1.0*r2,-1.0*r1);
                                tetha_s = acos(-1.0*r3/rm);   /* Not necessary */
                                
                                /* Bin construction and population*/
                                
                                if (rm <= 2.0)
                                {
                                    rm = 2.0001;
                                }
                                
                                rd = (int)ceil(log(1.0 - (1.0 - alfa)*(rm - 2.0)/d_r)/log(alfa));  /*check in test code !*/
                                //printf("\n interparticle distance %d ", rd);
                                //re = 2 + (exp(rd)*d_r*alfa-1)/(alfa - 1);
                                
                                //                        re = (1.0 - exp(log(alfa)*rd))*d_r/(1-alfa)+2.0;
                                //printf("\n %d  \t %f \t %f", rd, rm, re);
                                
                                if (phi >= 0.0)
                                {
                                    td = (int) ceil(phi/d_t);
                                }
                                
                                if (phi < 0.0)
                                {
                                    td = (int) ceil((phi + 2.0*pi)/d_t);
                                }
                                
                                /* For imposed symmetry pair */
                                
                                if (phi_s >= 0.0)
                                {
                                    td_s = (int) ceil(phi_s/d_t);
                                }
                                
                                if (phi_s < 0.0)
                                {
                                    td_s = (int) ceil((phi_s + 2.0*pi)/d_t);
                                }
                                
                                
                                /*Populate matrices*/
                                MP[rd-1][td-1] += 1.0;
                                //printf("\n interparticle distance %f ", MP[rd-1][td-1]);
                                /* Imposed symmetry bin population */
                                MP[rd-1][td_s-1] += 1.0;
                                
                                /*end of matrix population step*/
                                
                            }  /* #6 */
                            
                        }/* #5 */
                        
                    }/* #4 */
                }                     /* #3 */
                
                /****** Bin population and Normalization step *******/
                
                volume_f = (4.0/3.0)*pi*(rad_c*rad_c*rad_c);
                norm_f = 2.0*counter/volume_f;  /* 2 multiplies the previous norm for imposed symmetry consideration */
                counter = 0.0;             /* reset the counter */
                
                for (k = 0; k < nd; k++)
                {
                    for (m = 0; m < nt; m++)
                    {
                        
                        HP[k][m] += (MP[k][m]/norm_f);
                       //printf("\n interparticle distance %f ", HP[k][m]);
                        MP[k][m] = 0.0;
                    }
                    
                }
                
                
                
            }                           /* #2 */
            
            
        }                                /* #1 */
        
        for (k = 0; k < nd; k++)
        {
            for (m = 0; m < nt; m++)
            {
                HP[k][m] /= p_counter;   /* Devide by the number of legitimate particles because you do it that many times*/
                HPG[k][m] += HP[k][m];
                //printf("\n HPG %f ", HPG[k][m]);
                HP[k][m] = 0.0;
            }
           // printf("\n HPG %f ", HPG[k][m]);
        }
        //printf("\n HPG %f ", HPG[k-1][m-1]);
        p_counter = 0.0;
    }        /* Main loop ended */
    
    printf("\n I am at writing section");
    /*
    for (k = 0; k < nd; k++)
    {
        for (m = 0; m < nt; m++)
        {

            printf("\n HPG %f ", HPG[k][m]);
        }
    }
    // Writing the output
    
    
    outfile_ptr = fopen("HRs7.txt","w");
    if (outfile_ptr == NULL)
    {
        printf("\n Cannot open file HRs7.txt");
        exit(1);
    }
    
    for (i = 0; i < nd; i++)
    {
        for (j = 0; j < nt; j++)
        {
            
            fprintf(outfile_ptr, "%d \t %d \t %e \n ", i, j, HPG[i][j]);

            fprintf(outfile_ptr, " %e ", HPG[i][j]);
            
            
            if (j == nt - 1 )
            {
                fprintf(outfile_ptr, "\n");
            }
            
        }
    }
    fclose(outfile_ptr);
    */
    outfile_ptr = fopen("HRs6.txt","w");
    if (outfile_ptr == NULL)
    {
        printf("\n Cannot open file HRs4.txt");
        exit(1);
    }
    
    for (i = 0; i < nd; i++)
    {
        for (j = 0; j < nt; j++)
        {
            
         //   printf("\n hpg real %f ", HPG[i][j]);
            fprintf(outfile_ptr, "%d \t %d \t %e \n ", i, j, HPG[i][j]);
            
        }
    }
    fclose(outfile_ptr);
    
    for (i = 0; i < s; i ++)
    {
        free((void *) M[i]);
    }
    
    free((void *) M);
    
    
    for (i = 0; i < nd; i++)
    {
        free((void *) MP[i]);
        free((void *) HP[i]);
        free((void *) HPG[i]);
    }
    
    free((void *) MP);
    free((void *) HP);
    free((void *) HPG);
    
    
    end = clock();
    elapsed = ((double) (end - start))/ CLOCKS_PER_SEC;
    printf("\n the elapsed time is: %-1.12f", elapsed);
    
    
}


