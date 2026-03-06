#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <assert.h>

#define N 5
#define RayNum 100000


double vars[N];
double prevVars[N];
double chi;
double h;
double b,d;
double errtolmin, errtolmax;
double starting_r, kerr_starting_r;
int zcheck, zinit = 0;
double Mdot=0., eta, rho;
double kn_isco;
double g_tt_source;

double rCylindCoord[RayNum], emisDelta[RayNum], incDelta[RayNum];//, rmean[RayNum], intens0[RayNum], phflag[RayNum], meanphflag[RayNum];
double rSurf[RayNum], thetaSurf[RayNum];
double rBeforeHit[RayNum], thetaBeforeHit[RayNum], rAfterHit[RayNum];

double Q, g[4][4], dg_dr[4][4];
int imax = 100;


// #include "eqkerr.h"
#include "eqkn.h"
#include "rk45.h"
#include "kn.h"

using namespace std;

int main(int argc, char* argv[])
{
    clock_t begin = clock();


    double chi2, r2;// alpha, beta;
    double g_tt, g_tp, g_rr, g_thth, g_pp;
    double g_tt_r, g_tp_r, g_pp_r, rHit, rHit_prev;
    double E;
    double kt, kphi;
    double z, zPrev;
    int check;
    char filename[imax];
    double thetaAfterHitj, thetaBeforeHitj, rAfterHitj, rBeforeHitj, thetaDisk, dthdrhau;
    FILE *out;

    double height, delta, stepDelta;
    double I, A, LorentzF, Omega, g_lp;
    double disk_height;
    // Intensity, proper area, Lorentz factor, orbital velocity of disk, redshift factor, disk's radial thickness 

    double Emis_delta[2], III[imax-1], IIInr[imax-1], ED[imax-1], HD[imax-1];
    double rDisk[imax], rDiskNew[imax],  RayCount[imax], unorm, denom2;
    double frac,sign;
    double gamma=2, redshift;
    int i, j;
    int arr_index = 0;

    M = 1.;

    a = atof(argv[1]);      /* spin parameter */
    height = atof(argv[2]);   /* height of the source */
    // kerr_starting_r = atof(argv[3]);
    //Mdot = atof(argv[3]);
    Q= atof(argv[3]);   


    chi = a;

    snprintf(filename, sizeof(filename), "data/lp_%f_%f_%f.dat", a, height, Q);
    printf("%s  \n", filename);
    out = fopen(filename,"w");

    chi2 = chi*chi;

    h = -1.0;
    errtolmin = 1.0e-6;
    errtolmax = 1.0e-8;
    d = 2000.0;
    check = 0;

    kn_isco = find_isco();
    //fprintf(out, "%f %f %f\n", chi, height, johannsen_isco);
    // printf(" a = %f, h = %f, a13=%f, a22=%f, a52=%f, eps3=%f\n", a, height, a13, a22, a52, eps3);

    // starting_r = johannsen_isco * kerr_starting_r / kerr_rms(a);
    starting_r = 0.999 * kn_isco;
    eta = 1.0 - specific_energy(kn_isco);
    // printf(" isco1=%f, kerr=%f\n", kn_isco, kerr_rms(chi));


    for (i = 0; i < imax - 1; i++) {
        rDisk[i] = pow((double)i/(imax-2),3) * (1000. - kn_isco) + kn_isco;
        RayCount[i] = 0;
    }
    rDisk[imax-1] = 1050;

    //  main loop
    Emis_delta[1] = M_PI - 1e-8;
    Emis_delta[0] = 1E-8;
    delta = Emis_delta[0];
    stepDelta = (Emis_delta[1] - Emis_delta[0]) / ((double)RayNum - 1.);
    delta -= stepDelta;

    // for (i = 0; i < RayNum; i++) {
    j = 0; 
    for (i = 0; delta < Emis_delta[1]; i++, delta += stepDelta) {
        // printf("%d\n", i);

        // delta += stepDelta;

        zinit = 0;  //  additional parameter for runge-kutta 

        //  Initial values
        vars[0] = height;
        vars[1] = 1E-8;//asin(0.7747764 / height);
        vars[4] = 0.0;

        r2 = vars[0]*vars[0];

        denom2 = sqrt(r2 - 2*vars[0] + a*a);

        z = vars[0]*cos(vars[1]);
        zPrev = 10000.;

        metric(vars[0], vars[1], g);
        metric_rderivatives(vars[0], vars[1], dg_dr);

        g_thth = g[2][2];
        g_tt = g[0][0];
        g_tt_source = g[0][0];
        g_tp = g[0][3];
        g_rr = g[1][1];
        g_pp = g[3][3];

        vars[2] = cos(delta)/sqrt(g_rr);// denom2;//(-sin(vars[1]) * sin(delta) + cos(vars[1]) * cos(delta));
        vars[3] = sin(delta)/ sqrt(g_thth);//vars[3] = sin(delta)/ denom2;

        kphi = 0.;

        E = sqrt(g_tp*g_tp*kphi*kphi - g_tt*(g_rr*vars[2]*vars[2]+g_thth*vars[3]*vars[3]+g_pp*kphi*kphi));
        kt = -(g_tp*kphi+E)/g_tt;
        b = -(g_pp*kphi+g_tp*kt)/(g_tt*kt+g_tp*kphi);

        unorm = -(g_rr * pow(vars[2], 2) + g_thth * pow(vars[3], 2) + g_pp * pow(kphi, 2) + 2.0 * g_tp * kt * kphi) / (g_tt * pow(kt, 2));

        vars[2] /= E;
        vars[3] /= E;


        zcheck = 0;   // additional parameter for runge-kutta 
        rho = vars[0]*sin(fabs(vars[1]));
        // printf(" rho=%f, z=%f\n", rho, z);
        if(rho > kn_isco)
            disk_height = 3.0*Mdot*(1.0 - sqrt(kn_isco/rho))/eta;
        else
            disk_height = 0.0;


        for(; z > disk_height && vars[0] <= d && vars[0] > 1.0 + sqrt(1.0 - a*a - Q*Q) + 1.0e-4; ) {
            cache();
            rk();

            rho = vars[0]*sin(fabs(vars[1]));
            if(rho > kn_isco)
                disk_height = 3.0 * Mdot * (1.0 - sqrt(kn_isco/rho))/eta;
            else
                disk_height = 0.0;

            zPrev = z;
            z = vars[0]*cos(vars[1]);

            if (z <= disk_height && fabs(z - zPrev) > 1E-8) {
                zcheck = 1;
                zinit = 1;
                rcache();
                z = zPrev;
            }
        }


        if (vars[0] > d)
            delta = M_PI; // Finish the loop
        else if (vars[0] > 0.95 * starting_r) {

            rHit = find_radii(vars, prevVars);
            vars[1] = (vars[1]);
            emisDelta[arr_index] = delta;
            incDelta[arr_index] = cal_delta_inc(delta, vars[1], height, rHit);
            rCylindCoord[arr_index] = fabs(rHit*sin(vars[1]));
            rSurf[arr_index] = fabs(rHit);
            thetaSurf[arr_index] = (vars[1]);
            // phflag[arr_index] = 1

            rBeforeHit[arr_index] = fabs(prevVars[0]);
            thetaBeforeHit[arr_index] = (prevVars[1]);
            rAfterHit[arr_index] = fabs(vars[0]);

            arr_index++;
            zcheck = 0;     //  additional parameter for runge-kutta

            
            if (rDisk[j] <= rHit && rHit < rDisk[j+1])
            {
                RayCount[j] = sin(delta)*stepDelta/fabs(rHit - rHit_prev);
                rDiskNew[j] = rHit;
                j = j+1;
            }


            rHit_prev = rHit;

//            if (rDisk[0] <= rHit && rHit < rDisk[imax-1]) {
//                j = 0;
//                while(!(rDisk[j] <= rHit && rHit < rDisk[j+1])) j++;                
//                RayCount[j] += 1.0;
//            }
        }

    }  // i mainloop 



    j=0;
    for(i=0; i<imax-1; i++) {


        while(!(rCylindCoord[j] < rDisk[i] && rDisk[i] < rCylindCoord[j+1])) j++;

        assert(rCylindCoord[j] < rDisk[i] && rDisk[i] < rCylindCoord[j+1]);
        
        frac = (rCylindCoord[j+1] - rDisk[i]) / (rCylindCoord[j+1] - rCylindCoord[j]);
        thetaDisk = 0.5 * (1-frac) * (thetaSurf[j] + thetaBeforeHit[j]) + 0.5 * frac * (thetaSurf[j+1] + thetaBeforeHit[j+1]);

        thetaAfterHitj = (1-frac) * thetaSurf[j] + frac * thetaSurf[j+1];
        thetaBeforeHitj = (1-frac)*thetaBeforeHit[j] + frac*thetaBeforeHit[j+1];
        rAfterHitj = (1-frac)*rAfterHit[j] + frac*rAfterHit[j+1];
        rBeforeHitj = (1-frac)*rBeforeHit[j] + frac*rBeforeHit[j+1];

        denom2 = rAfterHitj * sin(thetaAfterHitj) - rBeforeHitj * sin(thetaBeforeHitj);

        if (denom2 != 0)
            dthdrhau = (thetaAfterHitj - thetaBeforeHitj) / denom2;


        ED[i] = (1-frac) * emisDelta[j] + frac * emisDelta[j+1];
        HD[i] = (1-frac) * incDelta[j] + frac * incDelta[j+1];

        if (ED[i] > M_PI/2.0) HD[i] *= (-1.0);

        metric(rDiskNew[i], thetaDisk, g);
        metric_rderivatives(rDiskNew[i], thetaDisk, dg_dr);

        if (Mdot == 0) A = 2 * M_PI * sqrt(g[1][1] * g[3][3]);
        else A = 2 * M_PI * sqrt(g[1][1] + g[2][2] * dthdrhau * dthdrhau * sin(thetaDisk) * sin(thetaDisk)) * sqrt(g[3][3] / sin(thetaDisk) / sin(thetaDisk));

        //A *= (rDisk[i+1] - rDisk[i]); A => dA/dr
        Omega = (-dg_dr[0][3] + sqrt(dg_dr[0][3]*dg_dr[0][3] - dg_dr[0][0]*dg_dr[3][3])) / dg_dr[3][3];                 
        LorentzF = pow(pow(Omega + g[0][3] / g[3][3], 2) * g[3][3] * g[3][3] / (g[0][0]*g[3][3] - g[0][3]*g[0][3]) + 1.0, -0.5);
        redshift = sqrt(g_tt_source / (g[0][0] + 2 * g[0][3] * Omega + g[3][3] * Omega * Omega));

        III[i] = RayCount[i] / (A * LorentzF);
        // III[i] *= pow(redshift, gamma);


    }


    for (i = 0; i<imax-1; i++) {

        if (isnan(III[i])) III[i] = 0.;
        // fprintf(out, "%f %f %f %e %f %f\n", a, height, rDisk[i], III[i], ED[i], HD[i]);
        //fprintf(out, "%f %e %f %f\n", rDisk[i], III[i], ED[i], HD[i]);
        fprintf(out, "%f %e\n", rDiskNew[i], III[i]);

    }




    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf(" Execution time is %f s\n", time_spent);

    fprintf(out, "\n");
    fclose(out);
    // fclose(fout);
    return 0;


}   

//////////////////////////////////////////////////
/////////////////  END OF MAIN ///////////////////
//////////////////////////////////////////////////


