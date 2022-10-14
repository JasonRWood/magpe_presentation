#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <time.h>
#include <cmath>
#include <unistd.h>
#include <iomanip> 
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <wrapper.h>
using namespace std;

float TINY = 1e-3;
float EPS = 1e-6;

int trait_space_length = 101;

namespace solvers{
    
    Quick_solver::Quick_solver() {};
    Quick_solver::Quick_solver(float beta_max, float alpha_max, float b, float q, float d, float gamma, int seed){
        this->beta_max = beta_max;
        this->alpha_max = alpha_max;
        this->b = b;
        this->q = q;
        this->d = d;
        this->gamma = gamma;
        this->seed = seed;

    };
    
    
    Quick_solver::~Quick_solver() {};

    
    float Quick_solver::fastmax(float a, float b){
        if (a > b){
            return a;
        }
        else{
            return b;
        }
    }

    float Quick_solver::fastmin(float a, float b){
        if (a < b){
            return a;
        }
        else{
            return b;
        }
    }

            // float fastmin(float a, float b){
            //     return (a<b)?a:b;
            // }

    int Quick_solver::dynamics(float* dydt, float* y, int num_strains, int* alpha_inds, float b, float q, float d, float* beta, float* alpha, float gamma){
        
        int i;
        float Isum = 0.0, N = 0.0;
        float S, I[num_strains];
        float Sdot, Idot[num_strains];
        float infection_force_para[num_strains], para_infection_force = 0.0;
        
        S = y[0];

        for(i=0;i<num_strains;i++){
            I[i] = y[i + 1];
            Isum = Isum + I[i];
            infection_force_para[i] = beta[alpha_inds[i]]*I[i];
            para_infection_force = para_infection_force + infection_force_para[i];
        }

        N = Isum + S;

        Sdot =  b*(1 - q*N)*N - (para_infection_force + d)*S + gamma*Isum;
        
        for(i = 0;i<num_strains;i++){
            Idot[i] = (beta[alpha_inds[i]]*S  - (d + alpha[alpha_inds[i]] + gamma))*I[i];
        }

        dydt[0] = Sdot;

        for(i = 0;i<num_strains;i++){

            dydt[i + 1] = Idot[i];
        }
        return 0;
    }

    int Quick_solver::RK45(float* y_out, float* y, float* y_err, int num_strains, int* alpha_inds, float b, float q, float d, float* beta, float* alpha, float gamma, float* h){

        float* y_1 = new float[num_strains + 1];
        float* y_2 = new float[num_strains + 1];
        float* y_3 = new float[num_strains + 1];
        float* y_4 = new float[num_strains + 1];
        // float y_1[num_strains+1], y_2[num_strains+1], y_3[num_strains+1], y_4[num_strains+1];
        int i;
        
        for (i = 0; i<(num_strains + 1); i++){
            y_1[i] = 0.0;
            y_2[i] = 0.0;
            y_3[i] = 0.0;
            y_4[i] = 0.0;
        }

        Quick_solver::dynamics(y_1, y, num_strains, alpha_inds, b, q, d, beta, alpha, gamma);

        for (i = 0; i<(num_strains + 1); i++){
            y_out[i] = y[i] + (y_1[i]*h[0]);
        }

        for (i = 0; i<(num_strains + 1); i++){
            if(isnan(y_1[i])){
                std::cout << "Ind " << i << " has created a NaN after 1 step\n";
            }
        }

        Quick_solver::dynamics(y_2, y_out, num_strains, alpha_inds, b, q, d, beta, alpha, gamma);

        for (i = 0; i<(num_strains + 1); i++){
            y_out[i] = y[i] + (y_2[i]*h[0])/2;
        }

        for (i = 0; i<(num_strains + 1); i++){
            if(isnan(y_2[i])){
                std::cout << "Ind " << i << " has created a NaN after 2 steps\n";
            }
        }
        Quick_solver::dynamics(y_3, y_out, num_strains, alpha_inds, b, q, d, beta, alpha, gamma);

        for (i = 0; i<(num_strains + 1); i++){
            y_out[i] = y[i] + (y_3[i]*h[0])/2;
        }
        
        for (i = 0; i<(num_strains + 1); i++){
            if(isnan(y_3[i])){
                std::cout << "Ind " << i << " has created a NaN after 3 steps\n";
            }
        }

        Quick_solver::dynamics(y_4, y_out, num_strains, alpha_inds, b, q, d, beta, alpha, gamma);

        for (i = 0; i<(num_strains + 1); i++){
            y_out[i] = y[i] + h[0]*((y_1[i] + 2*(y_2[i] + y_3[i]) + y_4[i])/6);
            y_err[i]= h[0]*((y_1[i] + 2*(y_2[i] + y_3[i]) + y_4[i])/6);
        }

        for (i = 0; i<(num_strains + 1); i++){
            if(isnan(y_4[i])){
                std::cout << "Ind " << i << " has created a NaN after 4 steps\n";
            }

            if(isnan(y_err[i])){
                std::cout << "error calculater " << i << " has returned a NaN\n";
            }
        }
        delete y_1;
        delete y_2;
        delete y_3;
        delete y_4;
        return 0;
    }


    int Quick_solver::perform_RK45_step(float* y, float* t, int num_strains, int* alpha_inds, float b, float q, float d, float* beta, float* alpha, float gamma, float* h, float* discrepancy, int* counter, int* num_plus, int* num_down){

        float y_out[num_strains+1], y_err[num_strains + 1];
        bool flag = true;
        int i;
        
            while ((h[0] >= 1e-10) && (h[0] <= 1e8) && (flag))
            {  
                for (i=0;i<(num_strains + 1);i++){
                    y_out[i] = 0.0;
                }

                discrepancy[0] = 0.0;

                Quick_solver::RK45(y_out,y,y_err,num_strains,alpha_inds, b,q,d,beta,alpha,gamma,h);
                float ysum = 0.0;
                for (i=0;i<(num_strains+1);i++){
                    ysum += y[i];
                }
                for (i=0;i<(num_strains+1);i++){
                    discrepancy[0] = fastmax(discrepancy[0],abs((y_out[i] - y[i])/(y[i] + EPS)));
                }

                if (discrepancy[0]/(1e-2) < 1){
                    break;
                }
                discrepancy[0] = fastmax(0.0000000001, discrepancy[0]);

                h[0] = h[0]*pow(discrepancy[0]/(1e-2), -0.2);
                counter[0]++;
            

                if(counter[0] >= 10){
                    flag = false;
                }
            }
            for (i = 0; i<(num_strains + 1); i++){
                y[i] = y_out[i];
            }
            t[0] = t[0] + h[0];
        return 0;
    }

    int Quick_solver::rkck(float *y, float *dydt, float *yout, float *yerr, float h, int num_strains, int* alpha_inds, float b, float q, float d, float* beta, float* alpha, float gamma)
    {
        int i;
        float* y_1 = new float[num_strains + 1];
        float* y_2 = new float[num_strains + 1];
        float* y_3 = new float[num_strains + 1];
        float* y_4 = new float[num_strains + 1];
        float* y_5 = new float[num_strains + 1];
        float* y_temp = new float[num_strains + 1];
        // float y_temp[num_strains + 1];
        static float b21=0.2,b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42=-0.9,b43=1.2,
        b51=-11.0/54.0,b52=2.5,b53=-70.0/27.0,b54=35.0/27.0,b61=1631.0/55296,
        b62=175.0/512.0,b63=575.0/13824.0,b64=44275.0/110592,b65=253.0/4096.0,
        c1=37.0/378.0,c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,dc5=-277.00/14336;
        float dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,dc4=c4-13525.0/55296.0,
        dc6=c6-0.25;
        
        for(i=0;i<(num_strains + 1);i++){
            y_temp[i] = y[i] + b21*h*dydt[i];
        }
        Quick_solver::dynamics(y_1, y_temp, num_strains, alpha_inds, b, q, d, beta, alpha, gamma);

        for(i=0;i<(num_strains + 1);i++){
            y_temp[i] = y[i]+h*(b31*dydt[i]+b32*y_1[i]);
        }
        Quick_solver::dynamics(y_2, y_temp, num_strains, alpha_inds, b, q, d, beta, alpha, gamma);

        for(i=0;i<(num_strains + 1);i++){
            y_temp[i]= y[i]+h*(b41*dydt[i]+b42*y_1[i]+b43*y_2[i]);
        }

        Quick_solver::dynamics(y_3, y_temp, num_strains, alpha_inds, b, q, d, beta, alpha, gamma);
        
        for(i=0;i<(num_strains + 1);i++){
            y_temp[i]= y[i]+h*(b51*dydt[i]+b52*y_2[i]+b53*y_2[i]+b54*y_3[i]);
        }

        Quick_solver::dynamics(y_4, y_temp, num_strains, alpha_inds, b, q, d, beta, alpha, gamma);
        
        for(i=0;i<(num_strains + 1);i++){
            y_temp[i]= y[i]+h*(b61*dydt[i]+b62*y_1[i]+b63*y_2[i]+b64*y_3[i]+b65*y_4[i]);
        }

        Quick_solver::dynamics(y_5, y_temp, num_strains, alpha_inds, b, q, d, beta, alpha, gamma);
        
        for(i=0;i<(num_strains + 1);i++){
            yout[i]= y[i]+h*(c1*dydt[i]+c3*y_2[i]+c4*y_3[i]+c6*y_5[i]);
            yerr[i]= h*(dc1*dydt[i]+dc3*y_2[i]+dc4*y_3[i]+dc5*y_4[i]+dc6*y_5[i]);
        }

        delete y_1;
        delete y_2;
        delete y_3;
        delete y_4;
        delete y_5;
        delete y_temp;
        return 0;
    }

    int Quick_solver::rkqs(float *y, float *dydt, float *h,float *hnext,float *yscale, int num_strains, int* alpha_inds, float b, float q, float d, float* beta, float* alpha, float gamma, int* counter, int* num_up, int* num_down, float* discrepancy, int* int_tracker, float* t)
    {
        // float* y_temp = new float[num_strains + 1];
        // float* yerr = new float[num_strains + 1];
        float y_temp[num_strains + 1], yerr[num_strains + 1];
        float htemp,errmax;
        int i;//, counter_internal;
        
        htemp= *h;
        
        Quick_solver::rkck(y, dydt, y_temp, yerr, *h, num_strains, alpha_inds, b, q, d, beta, alpha, gamma);
        *hnext= *h;
        
        // for(i=0;i<(num_strains + 1);i++){
        //     y_sum += y[i];
        // }
        // counter_internal = 0;
        for(;;)
        {
            Quick_solver::rkck(y, dydt, y_temp, yerr, *h, num_strains, alpha_inds, b, q, d, beta, alpha, gamma);
            errmax= 0.0;
            
            for(i=0;i<(num_strains + 1);i++){
                if (errmax < abs(yerr[i]/yscale[i])){
                    errmax = abs(yerr[i]/yscale[i]);
                    int_tracker[0] = i;
                }
                
            }
            errmax/= 1e-4;
            discrepancy[0] = errmax;
            if(errmax<=1.0){ 
                break;
            }
            
            htemp= 0.9*(*h)*pow(errmax,-0.25);
            
            *h= (*h>=0.0 ? Quick_solver::fastmax(htemp,0.1*(*h)) : Quick_solver::fastmin(htemp,0.1*(*h)));
            counter[0]++;
        }
        if(errmax > 1.89E-4) {
            *hnext= 0.9*(*h)*pow(errmax,-0.2);
        } 
        else {
            *hnext= 5.0*(*h);
            
        }
        
        for (i = 0; i<(num_strains + 1); i++){
            
            y[i] = y_temp[i];
            
        }
        return 0;
    }

    void Quick_solver::alpha_evo_only(float beta_max, float alpha_max, float b, float q, float d, float gamma, int seed, int alpha_init){
        
        time_t start,end;
        time (&start);
        

        srand(seed);
        
        int tmax = 1000, evo_steps = 1000, i, j;
        float beta[trait_space_length*trait_space_length], alpha[trait_space_length];
        float t[1], tcheck[1], Hsum;
        int num_strains = 1, alpha_inds[trait_space_length], ind,num_strains_2;
        int  num_strains_cleaned;

        float* y = new float[trait_space_length+1];
        float* y_check = new float[trait_space_length + 1];
        float* dydt = new float[trait_space_length + 1];
        float* yScale = new float[trait_space_length + 1];
        
        float* y_max = new float[trait_space_length + 1];
        float* y_min = new float[trait_space_length + 1];
        float* y_max_next = new float[trait_space_length + 1];
        float* y_min_next = new float[trait_space_length + 1];

        float discrep_check;
        int alpha_inds_cleaned[trait_space_length],counter=0;
        int sigma_inds_cleaned[trait_space_length];
        float I_temp[trait_space_length], H_temp[trait_space_length],total_density[trait_space_length];
        float cum_density[trait_space_length], cum_props[trait_space_length];
        int alpha_inds_2[trait_space_length];
        int sigma_inds_2[trait_space_length], int_tracker[1];
        float* y_temp = new float[trait_space_length+1];
        
        float h[1], hnext[1], discrepancy[1];

        int num_poss_outputs = 30000;
        int output_alpha_inds[num_poss_outputs], output_sigma_inds[num_poss_outputs], output_counter = 0, evo_step_tracker[num_poss_outputs];
        float output_parasite_density[num_poss_outputs], output_hyperparasite_density[num_poss_outputs], output_host_density[num_poss_outputs];
        bool output_hyperparasite_truths[num_poss_outputs];
        
        ind = 0;
        ofstream myfile, tracker_file;

        for(i=0;i<trait_space_length;i++){
            alpha[i] = (alpha_max*(i))/(trait_space_length - 1);
        }

        for(i=0;i<trait_space_length;i++){
            beta[i] = beta_max*sqrt(alpha[i]/alpha_max);
        }

        ofstream beta_values, alpha_values;
        beta_values.open("../data/beta_vals.csv");
        alpha_values.open("../data/alpha_vals.csv");

        for (i = 0; i<trait_space_length;i++){
            beta_values << i << ",";
            alpha_values << i << ",";
        }
        beta_values << "\n";
        alpha_values << "\n";

        

        for(i=0;i<trait_space_length;i++){
            beta_values << beta[i] << ",";
        }

        beta_values << "\n";
        alpha_values << "\n";

        for (i = 0; i<trait_space_length;i++){
            alpha_values << sqrt(alpha[i]/alpha_max) << ",";
        }

        alpha_values << "\n";

        beta_values.close();
        alpha_values.close();

        y[0] = (4.0)/(b/q);
        y[1] = (4.0)/(b/q);

        alpha_inds[0] = alpha_init;
        
        num_strains = 1;
        for (int evo_counter = 0; evo_counter < evo_steps; evo_counter++){
            
            t[0] = 0.0;
            h[0] = 1e-1;
            hnext[0] = 1e-1;

            float y_check_sum = 0.0;
            
            tcheck[0] = 10.0;
            for(i=0;i<(num_strains + 1);i++){
                y_check[i] = y[i];
                y_max[i] = y[i];
                y_min[i] = y[i];
                y_max_next[i] = y[i];
                y_min_next[i] = y[i];
                y_check_sum += y[i];
            }
            int step_count = 0;
            int check_count = 0;
            while (t[0] <= tmax){   
                int counter[1],num_plus[1],num_down[1];
                
                counter[0] = 0;
                num_plus[0] = 0;
                num_down[0] = 0;

                dynamics(dydt, y, num_strains, alpha_inds, b, q, d, beta, alpha, gamma);

                /* Adjust the step size to maintain accuracy */
                for (i=0; i<(num_strains + 1); i++){
                    yScale[i]=fabs(y[i])+fabs(dydt[i]*(*h))+EPS;
                    // grad_max = Quick_solver::fastmax(grad_max, abs(dydt[i]));
                }
                
                rkqs(y, dydt, h,hnext,yScale, num_strains, alpha_inds, b, q, d, beta, alpha, gamma,counter,num_plus,num_down,discrepancy, int_tracker, t);
                
                if (y[0] <= TINY){
                    y[0] = 0.0;
                }

                for(i=0;i<(num_strains);i++){
                    if (y[i+1] <= TINY){
                        y[i+1] = 0.0;
                    }
                }

                t[0]+=h[0];
                h[0] = hnext[0];
                
                for(i=0;i<(num_strains + 1);i++){
                    y_max[i] = fastmax(y[i], y_max[i]);
                    y_min[i] = fastmin(y[i], y_min[i]);
                }
                step_count += 1;
                
                if ((((t[0] - tcheck[0]) >= 10.0) && (step_count - check_count) >= 200) && (y[0] != 0.0)){
                    discrep_check = 0.0;
                    check_count = step_count;
                    
                    for(i=0;i<(num_strains + 1);i++){

                        if (y[i] >= TINY){
                            if (discrep_check < fastmax(abs(y_max[i] - y_max_next[i]),abs(y_min[i] - y_min_next[i]))){
                                
                                discrep_check = fastmax(abs(y_max[i] - y_max_next[i]),abs(y_min[i] - y_min_next[i]));
                            }
                        }
                    }

                    if (discrep_check >= 1e-3){

                        y_check_sum = 0.0;
                        tcheck[0] = t[0];

                        for(i=0;i<(num_strains + 1);i++){
                            y_check[i] = y[i];
                            y_max[i] = y[i];
                            y_min[i] = y[i];
                            y_max_next[i] = y[i];
                            y_min_next[i] = y[i];
                        }
                    }
                    else{
                        t[0] = tmax + 1.0;
                    }
                }
            }

            if (y[0] <= TINY){
                std::cout << "Healthy hosts extinct, breaking" << endl;
                break;
            }
            
            num_strains_cleaned = 0;
            for(i=0;i<num_strains;i++){
                if (y[i + 1] > TINY){
                    num_strains_cleaned = num_strains_cleaned + 1;
                }
            }

            if (num_strains_cleaned == 0){
                cout << "All parasites are dead\n";
                break;
            }

            float rand_stored1, rand_stored2;

            counter = 0;

            Hsum = 0.0;

            for (i=0;i<num_strains;i++){
                if(y[i + 1]  > TINY){
                    alpha_inds_cleaned[counter] = alpha_inds[i];
                    I_temp[counter] = y[i +1];
                    total_density[counter] = I_temp[counter];
                    counter++;
                }
            }
            

            cum_density[0] = total_density[0];
            if (num_strains_cleaned > 1){
                for(i=1;i<num_strains_cleaned;i++){
                    cum_density[i] = total_density[i] + cum_density[i - 1];
                }
            }

            for(i=0;i<num_strains_cleaned;i++){
                cum_props[i] = cum_density[i]/cum_density[num_strains_cleaned-1];
            }

            rand_stored1 = (float(rand())/float((RAND_MAX)));
            rand_stored2 = (float(rand())/float((RAND_MAX)));

            for(i=0;i<num_strains_cleaned;i++){
                if(cum_props[i] >= rand_stored1){
                    ind = i;
                    break;
                }
            }

            bool flag = true;
            int increment_ind;
            num_strains_2 = num_strains_cleaned;
            
            if (((rand_stored2 >=0.5) && (alpha_inds_cleaned[ind] != (trait_space_length-1))) || (alpha_inds_cleaned[ind] == 0)){
                for(i=0;i<num_strains_cleaned;i++){
                    if (alpha_inds_cleaned[i] == (alpha_inds_cleaned[ind] + 1)){
                        flag = false;
                        increment_ind = i;
                        break;
                    }
                }
            }
            else{
                for(i=0;i<num_strains_cleaned;i++){
                    if (alpha_inds_cleaned[i] == (alpha_inds_cleaned[ind] - 1)){
                        flag = false;
                        increment_ind = i;
                        break;
                    }
                }
            }
            
            if(flag){
                num_strains_2++;
            }

            y_temp[0] = y[0];
            for(i = 0; i<num_strains_2;i++){
                y_temp[i+1] = 0.0;
            }
            for(i=0;i<num_strains_cleaned;i++){
                alpha_inds_2[i] = alpha_inds_cleaned[i];
                y_temp[i+1] = I_temp[i];
            }
            num_strains = num_strains_2;

            for(i = 0; i<(num_strains+1); i++){
                y[i] = 0.0;
            }

            y[0] = y_temp[0];
        
            for(i=0;i<num_strains;i++){
                alpha_inds[i] = alpha_inds_2[i];
                y[i+1] = y_temp[i+1];
            }

            if (((rand_stored2 >=0.5) && (alpha_inds_cleaned[ind] != trait_space_length-1)) || (alpha_inds_cleaned[ind] == 0)){
                if(flag){
                    alpha_inds[num_strains-1] = alpha_inds_2[ind] + 1;
                    y[num_strains] = (y_temp[ind+1]/100);
                }
                else{
                    y[1+increment_ind] = y[1+increment_ind] + (y_temp[ind+1]/100);
                }
            }
            else{
                if(flag){
                    alpha_inds[num_strains-1] = alpha_inds_2[ind] - 1;
                    y[num_strains] = (y_temp[ind+1]/100);
                }
                else{
                    y[1+increment_ind] = y[1+increment_ind] + (y_temp[ind+1]/100);
                }
            }
            
            for(i=0;i<num_strains;i++){
                output_alpha_inds[output_counter] = alpha_inds[i];
                output_parasite_density[output_counter] = y[i+1];
                evo_step_tracker[output_counter] = evo_counter + 1;
                output_host_density[output_counter] = y[0];
                output_counter++;
            }

            if (evo_counter%100 == 0){
                std::cout << evo_counter << "\r";
            }
        }
        

        string seed_str;
        seed_str = std::to_string(seed);

        myfile.open("../data/evo_sims/data_set"+seed_str+".csv");
        myfile << "Trait_index_1,Density_of_Hosts,Density_of_parasite,Evolutionary_step,alpha_val,beta_val\n";
        
        for(i=0;i<output_counter;i++){
            myfile << output_alpha_inds[i] << "," << output_host_density[i] << "," << output_parasite_density[i] <<  ","<< evo_step_tracker[i] << "," << alpha[output_alpha_inds[i]] <<  "," << beta[output_alpha_inds[i]] << "\n";
        }

        myfile.close();
        
        delete y_check;
        delete y;
        delete yScale;
        delete dydt;
        delete y_temp;
        time (&end);
        float dif = difftime (end,start);
        printf ("Elasped time is %.2lf seconds.\n", dif ); 
        return;
    }

   
};
