#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

void expectedValues(double *mu,
                    double *k,
                    int *ymax,
                    int *n,
                    double *Es2,
                    double *Es2y,
                    double *Es1s2,
                    double *Es3) {
    
    
    int r, s;
    double temp1, temp2, temp3, temp4, dens;
    double tempa, tempb, tempc;
    for(r = 0; r < *n; r++) {
        temp1 = 0.0;
        temp2 = 0.0;
        temp3 = 0.0;
        temp4 = 0.0;
        tempa = 0.0;
        tempb = 0.0;
        tempc = 0.0;
        for(s = 0; s <= *ymax; s++) {
            if (s > 1 ) {
                tempa +=   (s-1) / (*k * (s-1) + 1);
                tempb +=   R_pow_di((s-1), 2) / R_pow_di((*k * (s-1) + 1), 2);
                tempc +=   R_pow_di((s-1), 3) / R_pow_di((*k * (s-1) + 1), 3);
            }
            dens =  dnbinom_mu(s, 1/k[0],  mu[r], 0);
            temp1 += dens*tempa*tempb;
            temp2 += dens*tempb;
            temp3 += dens*tempc;
            temp4 += dens*tempb*s;
        }
        Es2[r] = temp2;
        Es2y[r] = temp4;
        Es3[r] = temp3;
        Es1s2[r] = temp1;
    }
}


