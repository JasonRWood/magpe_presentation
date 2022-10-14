    #ifndef WRAPPER_H
    #define WRAPPER_H

    namespace solvers{

        class Quick_solver{

            public:
                float beta_max, alpha_max, sigma_max, b, q, d, rho, eta, gamma, lambda, c1, c2, hyper;
                int seed;
                Quick_solver();
                Quick_solver(float beta_max, float alpha_max, float b, float q, float d, float gamma, int seed);
                ~Quick_solver();
                float fastmax(float a, float b);
                float fastmin(float a, float b);
                int dynamics(float* dydt, float* y, int num_strains, int* alpha_inds, float b, float q, float d, float* beta, float* alpha, float gamma);
                int RK45(float* y_out, float* y, float* y_err, int num_strains, int* alpha_inds, float b, float q, float d, float* beta, float* alpha, float gamma, float* h);
                int perform_RK45_step(float* y, float* t, int num_strains, int* alpha_inds, float b, float q, float d, float* beta, float* alpha, float gamma, float* h, float* discrepancy, int* counter, int* num_plus, int* num_down);
                int rkck(float *y, float *dydt, float *yout, float *yerr, float h, int num_strains, int* alpha_inds, float b, float q, float d, float* beta, float* alpha, float gamma);
                int rkqs(float *y, float *dydt, float *h,float *hnext,float *yscale, int num_strains, int* alpha_inds,float b, float q, float d, float* beta, float* alpha, float gamma, int* counter, int* num_up, int* num_down, float* discrepancy, int* int_tracker, float* t);
                void alpha_evo_only(float beta_max, float alpha_max, float b, float q, float d, float gamma, int seed, int alpha_init);
       };
    }


    #endif