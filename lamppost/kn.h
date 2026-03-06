// #include "metric.h"

const  double Pi  = 3.141592653589793;
const  double Piby2 = 1.5707963267948966192;
double a, M;



void cache();
void rcache();

double R(double vars[5]);

double kerr_rms(double a);
double find_hitting_angle(double nvars[5], double ovars[5]);
double find_radii(double nvars[5], double ovars[5]);
double interp_lin_1d(double ifac_r, double rlo, double rhi);
double cal_q(double delta, double h);
double cal_delta_inc(double delta, double theta, double h, double r);
double cal_kt(double r, double theta);
double cal_gi(double delta, double theta, double h, double r);
double r0_grid(int i);
double r_mean(int i);
double I0(int i);
double II(double r);
double flat_I(double r,  double h);
void sortd(double ds[]);



// extern double a13;

void metric(double r, double th, double g[][4]);
void metric_rderivatives(double r, double th, double dg[][4]);
void metric_r2derivatives(double r, double th, double dg2[][4]);
void uppermetric(double r, double th, double rth[]);
double find_isco();

// void find_isco( double spin,  double spin2,  double eps3,  double a13,  double a22,  double a52,  double z1,  double& isco);



void metric(double r, double th, double g[][4])
{
double spin = chi;
double t3,  t6, t7, t8, t9, t10, t11, t12, t13, t17, t2, t5;
double gtt, grr, gthth, gpp, gtp;
  t2 = pow(r,2);
  t3 = pow(spin,2);
  t10 = pow(Q,2);
  t11 = -2*r;
  t5 = cos(th);
  t6 = pow(t5,2);
  t7 = t3*t6;
  t8 = t2 + t7;
  t9 = 1/t8;
  t12 = sin(th);
  t13 = pow(t12,2);
  t17 = t10 + t11 + t2 + t3;

  gtt = -((t10 + t11 + t2 + t3 - t13*t3)*t9);
  grr = t8/t17;
  gthth = t8;
  gpp = t13*(-(t13*t17*t3) + pow(t2 + t3,2))*t9;
  gtp = spin*(t10 + t11)*t13*t9;
  g[0][0] = gtt;
  g[0][3] = gtp;
  g[1][1] = grr;
  g[2][2] = gthth;
  g[3][0] = g[0][3];
  g[3][3] = gpp;
}

void metric_rderivatives(double r, double th, double dg[][4])
{
  double spin = chi;
  double  t2, t3, t5, t6, t7, t8, t9, t10, t11, t13, t14, t15, t16, t17, t22;
	double dgttdr, dgtpdr, dgppdr;
	t5 = pow(r,2);
	t6 = pow(spin,2);
	t7 = cos(th);
	t8 = pow(t7,2);
	t9 = t6*t8;
	t10 = t5 + t9;
	t11 = 1/t10;
	t16 = sin(th);
	t17 = pow(t16,2);
	t2 = 2*r;
	t3 = -2 + t2;
	t13 = pow(t10,-2);
	t22 = t5 + t6;
	t14 = pow(Q,2);
	t15 = -2*r;

	dgttdr = -(t11*t3) + 2*r*t13*(t14 + t15 + t5 + t6 - t17*t6);
	dgppdr = t11*t17*(4*r*t22 - t17*t3*t6) - 2*r*t13*t17*(pow(t22,2) - t17*t6*(t14 + t15 + t5 + t6));
	dgtpdr = -2*spin*t11*t17 - 2*r*spin*t13*(t14 + t15)*t17;

  dg[0][0] = dgttdr;
  dg[0][3] = dgtpdr;
  dg[3][0] = dg[0][3];
  dg[3][3] = dgppdr;
}

void metric_r2derivatives(double r, double th, double dg2[][4])
{
  double spin = chi;
	double  dgttdr2, dgtpdr2, dgppdr2;
	double t2, t3, t5, t6, t7, t8, t9, t10, t11, t13, t15, t16, t17, t18, t19, t20, t21, t26, t35, t36, t37, t38, t42;

	t5 = pow(r,2);
	t6 = pow(spin,2);
	t7 = cos(th);
	t8 = pow(t7,2);
	t9 = t6*t8;
	t10 = t5 + t9;
	t11 = pow(t10,-2);
	t16 = pow(Q,2);
	t17 = -2*r;
	t18 = sin(th);
	t19 = pow(t18,2);
	t20 = -(t19*t6);
	t21 = t16 + t17 + t20 + t5 + t6;
	t13 = 1/t10;
	t26 = t5 + t6;
	t2 = 2*r;
	t3 = -2 + t2;
	t15 = pow(t10,-3);
	t35 = pow(t26,2);
	t36 = t16 + t17 + t5 + t6;
	t37 = -(t19*t36*t6);
	t38 = t35 + t37;
	t42 = t16 + t17;
	
	dgttdr2 = -2*t13 + 2*t11*t21 + 4*r*t11*t3 - 8*t15*t21*t5;
	dgppdr2 = -2*t11*t19*t38 + 8*t15*t19*t38*t5 + t13*t19*(4*t26 + 8*t5 - 2*t19*t6) - 4*r*t11*t19*(4*r*t26 - t19*t3*t6);
	dgtpdr2 = 8*r*spin*t11*t19 - 2*spin*t11*t19*t42 + 8*spin*t15*t19*t42*t5;


  dg2[0][0] = dgttdr2;
  dg2[0][3] = dgtpdr2;
  dg2[3][0] = dg2[0][3];
  dg2[3][3] = dgppdr2;
}


void uppermetric(double r, double th, double rth[])
{
  double spin = chi;
	double t2, t3, t5, t6, t8, t9, t10, t11, t12, t13;
	double gurr, guthth;

	t5 = pow(r,2);
	t6 = pow(spin,2);
	t8 = cos(th);
	t9 = pow(t8,2);
	t10 = t6*t9;
	t11 = t10 + t5;
	t12 = 1/t11;
	
	gurr = t12*(pow(Q,2) - 2*r + t5 + t6);
	guthth = t12;

  rth[0] = gurr;
  rth[1] = guthth;
}


/* Calculate ISCO radius */
double find_isco()
{
    double d2Veff, d2Veff2, E_var, Lz_var, Omega_var, denom, den, den_r, den_rr, num, num_r, num_rr;
    double m[4][4],dmdr[4][4],dmdr2[4][4];
    double l, j, r, d2Veff_last = 0, d2Veff_last2 = 1000, rin;
    double jmin, jmax, lmin, lmax, rmin, factor, rstep;

    rstep = 0.1;
    rmin = 1.1;
    jmax = 15.0;
    jmin = rmin;
    factor = 1.0e2;

    for(j=jmax; j>jmin && j>rmin; j-=rstep/factor)
    {
        r = j;

    /* Calculate 2nd-derivative of the effective potential - d2Veff */
        metric(r, Pi/2., m);
        metric_rderivatives(r, Pi/2., dmdr);
        metric_r2derivatives(r, Pi/2., dmdr2);

        Omega_var = (-dmdr[0][3] + sqrt(dmdr[0][3]*dmdr[0][3] - dmdr[0][0]*dmdr[3][3])) / dmdr[3][3];
        denom = sqrt(-(m[0][0] + 2.0*m[0][3]*Omega_var + m[3][3]*Omega_var*Omega_var));
        E_var = -(m[0][0] + m[0][3]*Omega_var) / denom;
        Lz_var = (m[0][3] + m[3][3]*Omega_var) / denom;

        den = m[0][3]*m[0][3]-m[0][0]*m[3][3];
        den_r = 2.0*m[0][3]*dmdr[0][3]-dmdr[0][0]*m[3][3]-m[0][0]*dmdr[3][3];
        den_rr = 2.0*dmdr[0][3]*dmdr[0][3]+2.0*m[0][3]*dmdr2[0][3]-dmdr2[0][0]*m[3][3]-2.0*dmdr[0][0]*dmdr[3][3]-m[0][0]*dmdr2[3][3];
        num = E_var*E_var*m[3][3]+2.0*E_var*Lz_var*m[0][3]+Lz_var*Lz_var*m[0][0];
        num_r = E_var*E_var*dmdr[3][3]+2.0*E_var*Lz_var*dmdr[0][3]+Lz_var*Lz_var*dmdr[0][0];
        num_rr = E_var*E_var*dmdr2[3][3]+2.0*E_var*Lz_var*dmdr2[0][3]+Lz_var*Lz_var*dmdr2[0][0];

        d2Veff = (num_rr + (-2.0*num_r*den_r - num*den_rr + 2.0*num*den_r*den_r/den) / den) / den;

    /* Search for when d2Veff flips sign */
        if( (d2Veff < 0.0 && d2Veff_last > 0.0) || (d2Veff > 0.0 && d2Veff_last < 0.0) )
        {
            rin = j;
            lmax = rin + rstep/factor;
            lmin = rin - rstep/factor;
      jmin = lmin;
            factor *= 100.0;
            d2Veff_last2 = 1000;

      /* Zoom in on where d2Veff flips sign to increase ISCO accuracy */
            while(factor < 1.0e10)
            {
                for(l=lmax; l>lmin && l>rmin; l-=rstep/factor)
                {
                    r = l;

          /* Calculate 2nd-derivative of the effective potential - d2Veff */
                    metric(r, Pi/2., m);
                    metric_rderivatives(r, Pi/2., dmdr);
                    metric_r2derivatives(r, Pi/2., dmdr2);

                    Omega_var = (-dmdr[0][3] + sqrt(dmdr[0][3]*dmdr[0][3] - dmdr[0][0]*dmdr[3][3])) / dmdr[3][3];
                    denom = sqrt(-(m[0][0] + 2.0*m[0][3]*Omega_var + m[3][3]*Omega_var*Omega_var));
                    E_var = -(m[0][0] + m[0][3]*Omega_var) / denom;
                    Lz_var = (m[0][3] + m[3][3]*Omega_var) / denom;

                    den = m[0][3]*m[0][3]-m[0][0]*m[3][3];
                    den_r = 2.0*m[0][3]*dmdr[0][3]-dmdr[0][0]*m[3][3]-m[0][0]*dmdr[3][3];
                    den_rr = 2.0*dmdr[0][3]*dmdr[0][3]+2.0*m[0][3]*dmdr2[0][3]-dmdr2[0][0]*m[3][3]-2.0*dmdr[0][0]*dmdr[3][3]-m[0][0]*dmdr2[3][3];
                    num = E_var*E_var*m[3][3]+2.0*E_var*Lz_var*m[0][3]+Lz_var*Lz_var*m[0][0];
                    num_r = E_var*E_var*dmdr[3][3]+2.0*E_var*Lz_var*dmdr[0][3]+Lz_var*Lz_var*dmdr[0][0];
                    num_rr = E_var*E_var*dmdr2[3][3]+2.0*E_var*Lz_var*dmdr2[0][3]+Lz_var*Lz_var*dmdr2[0][0];

                    d2Veff2 = (num_rr + (-2.0*num_r*den_r - num*den_rr + 2.0*num*den_r*den_r/den) / den) / den;

                    if(fabs(d2Veff2)<fabs(d2Veff_last2))
                        rin = l;

                    d2Veff_last2 = d2Veff2;
                }

                lmax = rin + rstep/factor;
                lmin = rin - rstep/factor;
                factor *= 100.0;
                d2Veff_last2 = 1000;
            }

            factor = 1.0e2;
        }

        d2Veff_last = d2Veff;
    }

    // isco = rin;
    return rin;
}






double flat_I(double r,  double h) {
  return h / pow(r*r + h*h, 3./2.) / (2 * M_PI);
  // return 1. / pow(pow(r/h, 2.) + 1., 3./2.) / (2. * M_PI * pow(h, 2.));
}


void sortd(double ds[]){
  int i, j;
  double tmp;
  for (i=0;i<3;i++) 
    for (j=i+1; j<3; j++)
      if (ds[i] > ds[j]) {
        tmp = ds[i];
        ds[i] = ds[j];
        ds[j] = tmp;
      }
  
}


double cal_q(double delta, double height) {
  double a = chi;
  // double q = sqrt(a*a - 2.0*h + h*h) * sin(delta) / sqrt(a*a + h*h);
  double q = (a*a + height*height) * sin(delta);
  q *= q;
  q /= (height*height - 2*height + a*a);
  q -= a*a;
  q = sqrt(fabs(q));
  // q = sqrt((h*h - 2*h + a*a) / (a*a + h*h)) * sin(delta);
  return q;
}

void redshift(double r, double th, double ktkp, double& gg);

double cal_delta_inc(double delta, double th, double height, double r) {

  // th = M_PI/2.0;
  th = fabs(th);
  double theta = th;
  double rder,rder1,rder2,thder,thder1,thder2,beta,t1,c1,c2,g[4][4],dg[4][4], rth[2];
  double angle;
  double gg;
  double Zr, Zth, K, cs1, ss1, ss3, r3, gurr, guthth;
  
  metric(height, 0.0, g);
  g_tt_source = g[0][0];

  metric(r, th, g);
  metric_rderivatives(r, th, dg_dr);
  uppermetric(r, th, rth);

  double Omega = (-dg_dr[0][3] + sqrt(dg_dr[0][3]*dg_dr[0][3] - dg_dr[0][0]*dg_dr[3][3])) / dg_dr[3][3];

  gg = sqrt(g_tt_source / (g[0][0] + 2 * g[0][3] * Omega + g[3][3] * Omega * Omega));
  if (isnan(gg)) {
    printf(" gg is %f, %f, %f, %f, %f\n", gg, -g_tt_source, (-g[0][0] - 2 * g[0][3] * Omega - g[3][3] * Omega * Omega), -g_tt_source / (-g[0][0] - 2 * g[0][3] * Omega - g[3][3] * Omega * Omega), Omega);
  }
  K = 3. / eta * Mdot;
  cs1 = cos(th);
  ss1 = sin(th);
  ss3 = ss1*ss1*ss1;
  r3 = r*r*r;
  // Zr = 0.5 * K * sqrt(johannsen_isco / (r3 * ss1)) - cs1;
  Zr = cs1 - 0.5 * K * sqrt(kn_isco / (r3 * ss1));
  Zth = - 0.5 * K * cs1 * sqrt(kn_isco / (r * ss3)) - r*ss1;

  gurr = rth[0];
  guthth = rth[1];

  double vr = (vars[2] + prevVars[2]) / 2.0;
  double vth = (vars[3] + prevVars[3]) / 2.0;
  vr = prevVars[2];
  vth = prevVars[3];

  // printf(" vr=%f, vth=%f\n", vr, vth);

  angle = 1. / sqrt(gurr * Zr * Zr + guthth * Zth * Zth) * (- gurr * g[1][1] * Zr * vr - guthth * g[2][2] * Zth * vth);
  angle *= sqrt(-(g[0][0] + 2 * g[0][3] * Omega + g[3][3] * Omega * Omega));
  angle = acos(angle);

  // printf(" %f %f \n", angle, acos(cal_q(delta, height) / (r * sqrt(-1 / (g[0][0] + 2 * g[0][3] * Omega + g[3][3] * Omega * Omega)))));
  

  if (isnan(angle)) 
    printf(" angle, %f, %f %f\n", gg, gurr * Zr * Zr + guthth * Zth * Zth, Zr * vr + Zth * vth);

  return angle;  
}

// double cal_gi(double delta, double theta, double h, double r) {
//   double a=chi;
//   return cal_kt(r, theta)*sqrt(h*h-2*h+a*a)/sqrt(h*h + a*a);
// }

double cal_kt(double r, double theta) {
  double g_tt, g_pp, g_tp, denom, kt, a = chi, M=1.0;
  double r2, a2, s2, c2;
  r2 = r*r;
  a2 = a*a;
  s2 = sin(theta)*sin(theta);
  c2 = cos(theta)*cos(theta);
  g_tt = -1.0 + (2.0*M*r)/(r2 + a2*c2);
  // g_tt = -1. + (2.*M*r)/(pow(r,2) + pow(a,2)*pow(cos(theta),2));
  // g_pp = pow(sin(theta),2)*(pow(a,2) + pow(r,2) + (2*pow(a,2)*M*r*pow(sin(theta),2))/(pow(r,2) + pow(a,2)*pow(cos(vars[1]),2)));
  g_pp = s2*(a2 + r2 + (2.0*a2*M*r*s2)/(r2 + a2*c2));
  g_tp = (-2.0*a*M*r*s2)/(r2 + a2*c2);
  // g_tp = (-2.0*a*M*r*pow(sin(theta),2))/(pow(r,2) + pow(a,2)*pow(cos(theta),2));  

  
  denom = (g_tt*g_pp-g_tp*g_tp);
  kt = -(g_pp+b*g_tp)/denom;
  // printf(" gtt=%f, gpp=%f, gtp=%f ", g_tt, g_pp, g_tp);
  return kt;
}


void cache() {
  int i = 0;
  for (;i<N;i++)
    prevVars[i] = vars[i];
}

void rcache() {
  int i = 0;
  for (;i<N;i++)
    vars[i] = prevVars[i];
}

double R(double vars[5]) {
  double x, y, z;
  
  x = vars[0]*sin(vars[1])*cos(vars[4]);
  y = vars[0]*sin(vars[1])*sin(vars[4]);
  z = vars[0]*cos(vars[1]);

  return sqrt(x*x + y*y + z*z);
}

double interp_lin_1d(double ifac_r, double rlo, double rhi){
  return ifac_r*rhi + (1.0-ifac_r)*rlo;
}


double find_radii(double nvars[5], double ovars[5]) {
  return (ovars[0]+nvars[0])/2.0;//sqrt(x * x + y * y);
}

void redshift(double r, double th, double ktkp, double& gg)
{
  double Omega;
  double uet;
  double g00,g03,g33,rderg00,rderg03,rderg33;
  double spin = chi;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32;

	t2 = pow(r,2);
	t3 = pow(spin,2);
	t5 = cos(th);
	t6 = pow(t5,2);
	t7 = t3*t6;
	t8 = t2 + t7;
	t9 = 1/t8;
	t12 = sin(th);
	t13 = pow(t12,2);
	t10 = pow(Q,2);
	t11 = -2*r;
	t14 = -(t13*t3);
	t15 = t10 + t11 + t14 + t2 + t3;
	t17 = t2 + t3;
	t25 = 2*r;
	t26 = -2 + t25;
	t28 = pow(t8,-2);
	t18 = pow(t17,2);
	t19 = t10 + t11 + t2 + t3;
	t20 = -(t13*t19*t3);
	t21 = t18 + t20;
	t23 = t10 + t11;

	g00 = -(t15*t9);
	g33 = t13*t21*t9;
	g03 = spin*t13*t23*t9;

	rderg00 = 2*r*t15*t28 - t26*t9;
	rderg33 = -2*r*t13*t21*t28 + t13*(4*r*t17 - t13*t26*t3)*t9;
	rderg03= -2*r*spin*t13*t23*t28 - 2*spin*t13*t9;

  Omega  = (-rderg03 + sqrt(rderg03*rderg03 - rderg00*rderg33))/rderg33;
    
  uet = sqrt(-g00 - 2.*g03*Omega - g33*Omega*Omega);
    
  gg = uet/(1. - ktkp*Omega);
}

double specific_energy(double r)
{
  double Omega, SE;
  double g00,g03,g33,rderg00,rderg03,rderg33;
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32;
  double spin = chi;

	t2 = pow(r,-2);
	t6 = pow(r,2);
	t9 = pow(spin,2);
	t3 = pow(Q,2);
	t5 = -2*r;
	t7 = t3 + t5 + t6;
	t18 = 2*r;
	t19 = -2 + t18;
	t10 = t6 + t9;
	t21 = pow(r,-3);
	t11 = pow(t10,2);
	t12 = t3 + t5 + t6 + t9;
	t13 = -(t12*t9);
	t14 = t11 + t13;
	t16 = t3 + t5;

	g00 = -(t2*t7);
	g33 = t14*t2;
	g03 = spin*t16*t2;

	rderg00= -(t19*t2) + 2*t21*t7;
	rderg33  = -2*t14*t21 + t2*(4*r*t10 - t19*t9);
	rderg03 = -2*spin*t2 - 2*spin*t16*t21;


  Omega  = (-rderg03 + sqrt(rderg03*rderg03 - rderg00*rderg33))/rderg33;
    
  SE = -1.0*(g00 + Omega*g03)/sqrt(-1.0*g00 - 2.0*Omega*g03 - Omega*Omega*g33);
    
  return SE;
}

