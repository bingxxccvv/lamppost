void equations(double varst[], double diffs[]);

void equations(double varst[], double diffs[])
{
	double r, th;
	double kt, kphi;
	double kt2, kr2, kth2, kp2, ktp, krth;
	double ch_rtt, ch_rtp, ch_rrr, ch_rrth, ch_rthth, ch_rpp;
	double ch_thtt, ch_thtp, ch_thrr, ch_thrth, ch_ththth, ch_thpp;
	double g_tt, g_tp, g_pp, gurr, guthth;
	double dgttdr, dgttdth, dgtpdr, dgtpdth, dgrrdr, dgrrdth, dgththdr, dgththdth, dgppdr, dgppdth;
	double hgurr, hguthth;
	double denom;
	double spin = chi;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52;
	double unorm, g_rr, g_thth;
	spin = chi;
	r = varst[0];
	th = varst[1];

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
	t19 = t10 + t11 + t2 + t3;
	t14 = -(t13*t3);
	t15 = t10 + t11 + t14 + t2 + t3;
	t26 = 2*r;
	t27 = -2 + t26;
	t17 = t2 + t3;
	t29 = pow(t8,-2);
	t18 = pow(t17,2);
	t20 = -(t13*t19*t3);
	t21 = t18 + t20;
	t23 = t10 + t11;
	t32 = 1/t19;
	t51 = pow(t12,3);

	g_tt = -(t15*t9);
	g_pp = t13*t21*t9;
	g_tp = spin*t13*t23*t9;

	gurr = t19*t9;
	guthth = t9;

	dgttdr = 2*r*t15*t29 - t27*t9;
	dgrrdr = 2*r*t32 - (t27*t8)/pow(t19,2);
	dgththdr = t26;
	dgppdr = -2*r*t13*t21*t29 + t13*(4*r*t17 - t13*t27*t3)*t9;
	dgtpdr = -2*r*spin*t13*t23*t29 - 2*spin*t13*t9;

	dgttdth = -2*t12*t15*t29*t3*t5 + 2*t12*t3*t5*t9;
	dgrrdth = -2*t12*t3*t32*t5;
	dgththdth = -2*t12*t3*t5;
	dgppdth = 2*t21*t29*t3*t5*t51 + 2*t12*t21*t5*t9 - 2*t19*t3*t5*t51*t9;
	dgtpdth = 2*pow(spin,3)*t23*t29*t5*t51 + 2*spin*t12*t23*t5*t9;

	g_thth = hgurr; 
	g_rr = hguthth;

	hgurr = 0.5*gurr;
	hguthth = 0.5*guthth;
	
	ch_rtt = -hgurr*dgttdr;
	ch_rtp = -hgurr*dgtpdr;
	ch_rrr = hgurr*dgrrdr;
	ch_rrth = hgurr*dgrrdth;
	ch_rthth = -hgurr*dgththdr;
	ch_rpp = -hgurr*dgppdr;
	ch_thtt = -hguthth*dgttdth;
	ch_thtp = -hguthth*dgtpdth;
	ch_thrr = -hguthth*dgrrdth;
	ch_thrth = hguthth*dgththdr;
	ch_ththth = hguthth*dgththdth;
	ch_thpp = -hguthth*dgppdth;
	
	//printf("%Le %Le\n", r, th);
	//printf("%Le %Le %Le %Le %Le %Le\n", ch_rtt, ch_rtp, ch_rrr, ch_rrth, ch_rthth, ch_rpp);
	//printf("%Le %Le %Le %Le %Le %Le\n", ch_thtt, ch_thtp, ch_thrr, ch_thrth, ch_ththth, ch_thpp);
	
	denom = (g_tt*g_pp-g_tp*g_tp);
	
	kt = -(g_pp+b*g_tp)/denom;
	kphi = (g_tp+b*g_tt)/denom;

	diffs[0] = varst[2];
	
	diffs[1] = varst[3];
	
	kt2 = kt*kt;
	kr2 = varst[2]*varst[2];
	kth2 = varst[3]*varst[3];
	kp2 = kphi*kphi;
	ktp = kt*kphi;
	krth = varst[2]*varst[3];


	// unorm = -(g_rr * kr2 + g_thth * kth2 + g_pp * pow(kphi, 2) + 2.0 * g_tp * kt * kphi) / (g_tt * pow(kt, 2));
	// printf(" unorm=%f\n", unorm);



	diffs[2] = -(ch_rtt*kt2+ch_rrr*kr2+ch_rthth*kth2+ch_rpp*kp2+2.0*(ch_rtp*ktp+ch_rrth*krth));
				
	diffs[3] = -(ch_thtt*kt2+ch_thrr*kr2+ch_ththth*kth2+ch_thpp*kp2+2.0*(ch_thtp*ktp+ch_thrth*krth));
}	