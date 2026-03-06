#ifndef max
#define max(a,b) (( (a) >= (b)) ? (a) : (b))
#endif

double a1 = 1.0/4.0;
double b1 = 3.0/32.0;
double b2 = 9.0/32.0;
double c1 = 1932.0/2197.0;
double c2 = -7200.0/2197.0;
double c3 = 7296.0/2197.0;
double d1 = 439.0/216.0;
double d2 = -8.0;
double d3 = 3680.0/513.0;
double d4 = -845.0/4104.0;
double e1 = -8.0/27.0;
double e2 = 2.0;
double e3 = -3544.0/2565.0;
double e4 = 1859.0/4104.0;
double e5 = -11.0/40.0;
double x1 = 25.0/216.0;
double x2 = 0.0;
double x3 = 1408.0/2565.0;
double x4 = 2197.0/4104.0;
double x5 = -1.0/5.0;
double z1 = 16.0/135.0;
double z2 = 0.0;
double z3 = 6656.0/12825.0;
double z4 = 28561.0/56430.0;
double z5 = -9.0/50.0;
double z6 = 2.0/55.0;

int zch;


void rk();

void rk()
{
	double vars_temp[N], diffs[N], k1[N], k2[N], k3[N], k4[N], k5[N], k6[N];
	double y[N], z[N];
	double err;
	int i, check = 0;
  if (zinit == 0) zch = 0;

	
	equations(vars, diffs);
	for(i = 0; i < N; i++)
	{
		k1[i] = h*diffs[i];
		vars_temp[i] = vars[i] + a1*k1[i];
	}

	equations(vars_temp, diffs);
	for(i = 0; i < N; i++)
	{
		k2[i] = h*diffs[i];
    	vars_temp[i] = vars[i] + b1*k1[i] + b2*k2[i];
  }

	equations(vars_temp, diffs);
	for(i = 0; i < N; i++)
	{
   		k3[i] = h*diffs[i];
    	vars_temp[i] = vars[i] + c1*k1[i] + c2*k2[i] + c3*k3[i];
  }

  	equations(vars_temp, diffs);
  	for(i = 0; i < N; i++)
  	{
   		k4[i] = h*diffs[i];
    	vars_temp[i] = vars[i] + d1*k1[i] + d2*k2[i] + d3*k3[i] + d4*k4[i];
  	}
  	
  	equations(vars_temp, diffs);
  	for(i = 0; i < N; i++)
  	{
   		k5[i] = h*diffs[i];
    	vars_temp[i] = vars[i] + e1*k1[i] + e2*k2[i] + e3*k3[i] + e4*k4[i] + e5*k5[i];
  	}
  	
  	equations(vars_temp, diffs);
  	for(i = 0; i < N; i++)
   		k6[i] = h*diffs[i];
  	
  	for(i = 0; i < N; i++)
  	{
	  	y[i] = vars[i] + x1*k1[i] + x2*k2[i] + x3*k3[i] + x4*k4[i] + x5*k5[i];
	  	z[i] = vars[i] + z1*k1[i] + z2*k2[i] + z3*k3[i] + z4*k4[i] + z5*k5[i] + z6*k6[i];
	  	
	  	err = fabs((y[i]-z[i])/max(vars[i], y[i]));
	  	
	  	if(err > errtolmin || zcheck) {
          check = 1;
          if (zcheck) {
              zch = 1;
              zcheck = 0;
          }
        }
  		else if(err < errtolmax && check != 1 && zch != 1)
	  		check = -1;
  	}
    
  	if(check == 1)
  	{
	  	h /= 2.0;
      // printf(" ---  h = %f\n", h);
	  	rk();
  	}
  	else if(check == -1)
  	{
  		// if (fabs(h) < 0.01)
	  	h *= 2.0;
	  	for(i = 0; i < N; i++)
  			vars[i] = y[i];
	  }
  	else
  		for(i = 0; i < N; i++)
  			vars[i] = y[i];

       // printf(" %f %f %f %f\n", vars[0], vars[1], vars[2], vars[3]);


}