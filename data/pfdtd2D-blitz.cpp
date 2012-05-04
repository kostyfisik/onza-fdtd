#include <ctime>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <blitz/array.h>

bool pendry = false;
// bool pendry = true;
bool zhao = false;
// bool zhao = true;
bool cai = false;
// bool cai = true;
bool belov = false;
// bool belov = true;
bool tretyakov = true;
// bool tretyakov = false;

using namespace blitz;
BZ_DECLARE_STENCIL2(averX,B,A)
    A = (B(0,0) + B(0,1) + B(-1,0) + B(-1,1)) / 4.0;
BZ_END_STENCIL
BZ_DECLARE_STENCIL2(averY,B,A)
    A = (B(0,0) + B(1,0) + B(0,-1) + B(1,-1)) / 4.0;
BZ_END_STENCIL
#define PI M_PI           // from cmath
#define NPML 10           // Depth of PML region in # of cells
#define IE 400            // Width
#define JE 500            // Height
#define NMAX 1200          // Total time steps
void print2DToFile(string fileName, Array<double,2> inArray);
double getEpsTheta(double dist, double innerR, double outerR);
double getEpsR(double dist, double innerR, double outerR);
double getMuZ(double dist, double innerR, double outerR);
double getMuA(double dist, double innerR, double outerR);

int main(int argc, char *argv[])
{
    std::cout<<std::endl<<"*** Start... ***"<<std::endl;
    double c_0=2.9979e8,         // Speed of light in free space
	mu_0=4.0*PI*1.0e-7,   // Permeability of free space
	eps_0=8.8542e-12;     // Permittivity of free space

    // double dx = 4.0e-3;           // Space increment of square lattice
    double dx = 0.55e-3;           // Space increment of square lattice
    // double dt = dx/c_0/sqrt(2);   // Time step
    double dt = dx/c_0/sqrt(2);   // Time step

    int R1 = round(0.03/dx);   // Inner R of cloak
    int R2 = round(0.11/dx);   // Outer R

    double freq=8.49e9;        // Source freq
    // freq = 8.4e9           //Belov freq
    // omega_0 = 8e9          //Belov freq
    int is = IE/2;   // Location of plane-wave source
    int js = 30;     // Location of plane-wave source
    clock_t t1 = clock();
    int  m, n, l, i, j,	ic, jc, N, M, ST=300, ne=0;
    double omega, tau, delay, x, g,
	sigmax, rhomax, re, rm, dist, eps_r, eps_theta,
	mu_z, sinx, siny, cosx, cosy, tan_l, omega_px, gamma_px, 
	omega_py, gamma_py, omega_pz, gamma_pz, A;
    
    omega = 2.0*PI*freq;

    Array<double,1> sig(NPML), rho(NPML);
    Array<double,1> ca(NPML), cb(NPML), da(NPML), db(NPML);
    Array<double,1> source(NMAX), err(NMAX);
    //***********************************************************************
    // Initialise matrices for field components
    //***********************************************************************
    int ie = IE;     // # of grid cells in x-direction
    int je = JE;     // # of grid cells in y-direction
    Array<double, 2> Dx(ie,je), Dx_h1(ie,je), Dx_h2(ie,je);
    Array<double, 2> Dx_a(ie,je), Dx_a_h1(ie,je), Dx_a_h2(ie,je);
    Array<double, 2> Dy(ie,je), Dy_h1(ie,je), Dy_h2(ie,je);
    Array<double, 2> Dy_a(ie,je), Dy_a_h1(ie,je), Dy_a_h2(ie,je);
    Array<double, 2> caDx(ie,je), cbDx(ie,je);
    Array<double, 2> caDy(ie,je), cbDy(ie,je);
    Array<double, 2> Ex(ie,je), Ex_h1(ie,je), Ex_h2(ie,je), Ex_r(ie,je), Ex_i(ie,je);
    Array<double, 2> Ey(ie,je), Ey_h1(ie,je), Ey_h2(ie,je), Ey_r(ie,je), Ey_i(ie,je);
    Array<double, 2> Hz(ie,je), Hz_h1(ie,je), Hz_h2(ie,je), Hz_p(ie,je), Hz_r(ie,je), Hz_i(ie,je);
    Array<double, 2> Bz(ie,je), Bz_h1(ie,je), Bz_h2(ie,je);
    Array<double, 2> Bzx(ie,je), Bzy(ie,je);
    Array<double, 2> daBzx(ie,je), dbBzx(ie,je);
    Array<double, 2> daBzy(ie,je), dbBzy(ie,je);
    // Coefficients for cloak
    Array<double, 2> a0x(ie,je), a1x(ie,je), a2x(ie,je);
    Array<double, 2> b0xx(ie,je), b1xx(ie,je), b2xx(ie,je);
    Array<double, 2> b0xy(ie,je), b1xy(ie,je), b2xy(ie,je);
    Array<double, 2> a0y(ie,je), a1y(ie,je), a2y(ie,je);
    Array<double, 2> b0yy(ie,je), b1yy(ie,je), b2yy(ie,je);
    Array<double, 2> b0yx(ie,je), b1yx(ie,je), b2yx(ie,je);
    Array<double, 2> a0z(ie,je), a1z(ie,je), a2z(ie,je);
    Array<double, 2> b0z(ie,je), b1z(ie,je), b2z(ie,je);

    Dx = 0.0; Dx_h1 = 0.0; Dx_h2 = 0.0;
    Dx_a = 0.0; Dx_a_h1 = 0.0; Dx_a_h2 = 0.0;
    Dy = 0.0; Dy_h1 = 0.0; Dy_h2 = 0.0;
    Dy_a = 0.0; Dy_a_h1 = 0.0; Dy_a_h2 = 0.0;
    caDx = 1.0; cbDx = dt/eps_0/dx;
    caDy = 1.0; cbDy = dt/eps_0/dx;
    Ex = 0.0; Ex_h1 = 0.0; Ex_h2 = 0.0; Ex_r = 0.0; Ex_i = 0.0;
    Ey = 0.0; Ey_h1 = 0.0; Ey_h2 = 0.0; Ey_r = 0.0; Ey_i = 0.0;
    Hz = 0.0; Hz_h1 = 0.0; Hz_h2 = 0.0; Hz_p = 1.0; Hz_r = 0.0; Hz_i = 0.0;
    Bz = 0.0; Bz_h1 = 0.0; Bz_h2 = 0.0;
    Bzx = 0.0; Bzy = 0.0;
    daBzx = 1.0; dbBzx = dt/mu_0/dx;
    daBzy = 1.0; dbBzy = dt/mu_0/dx;
    // Coefficients for cloak
    a0x = 1.0; a1x = 0.0; a2x = 0.0;
    b0xx = 1.0; b1xx = 0.0; b2xx = 0.0;
    b0xy = 0.0; b1xy = 0.0; b2xy = 0.0;
    a0y = 1.0; a1y = 0.0; a2y = 0.0;
    b0yy = 1.0; b1yy = 0.0; b2yy = 0.0;
    b0yx = 0.0; b1yx = 0.0; b2yx = 0.0;
    a0z = 1.0; a1z = 0.0; a2z = 0.0;
    b0z = 1.0; b1z = 0.0; b2z = 0.0;

    int firstI = 0;    int lastI = IE-1;
    int firstJ = 0;    int lastJ = JE-1;
    Range woLastI = Range(firstI, lastI-1);
    Range woLastJ = Range(firstJ, lastJ-1);
    Range woEdgeI = Range(firstI+1, lastI-1);
    Range woEdgeJ = Range(firstJ+1, lastJ-1);

    tau = 1.0/freq;
    delay = 3.0*tau/dt;
    N = round(tau/dt);
    M = round(c_0/freq/dx);

    /* Source excitation */
    for(n=0;n<NMAX;n++) {
        if(n<ST*N) {
            x = 1.0 - 1.0*(ST*N-n)/(ST*N);
            g = 10.0*pow3(x) - 15.0*pow4(x) + 6.0*pow5(x);
            source(n) = g * sin(2*PI*freq*n*dt);
        }
        else {
            source(n) = sin(2*PI*freq*n*dt);
        }
//        source(n) = exp(-((n-delay)*dt/tau)^2) * sin(2*PI*freq*n*dt);
    }

    //***********************************************************************
    // Setup Berenger's PML material constants and coefficients
    //***********************************************************************
    sigmax = -3.0*eps_0*c_0*log(1.0e-5)/(2.0*dx*NPML);
    rhomax = sigmax*(mu_0/eps_0);
    for(m=0;m<NPML;m++) {
        sig(m) = sigmax*pow2((m+0.5)/(NPML+0.5));
        rho(m) = rhomax*pow2((m+1)/(NPML+0.5));
        re = sig(m)*dt/eps_0;
        rm = rho(m)*dt/mu_0;
        ca(m) = exp(-re);
        cb(m) = -(exp(-re)-1.0)/sig(m)/dx;
        da(m) = exp(-rm);
        db(m) = -(exp(-rm)-1.0)/rho(m)/dx;
    }
    //***********************************************************************
    // Initialize all of the matrices for the Berenger's PML
    //***********************************************************************
    int ip = IE - 1 - NPML;   //PML border
    int jp = JE - 1 - NPML;   //PML border
    for(i=0;i<IE-1;i++) {     //Left & right
	for(j=0;j<NPML;j++) {
	    m = NPML - 1 - j;
	    daBzy(i,j) = da(m);
	    dbBzy(i,j) = db(m);
	}
	for(j=1;j<NPML+1;j++) {
	    m = NPML - j;
	    caDx(i,j) = ca(m);
	    cbDx(i,j) = cb(m);
	}
	for(j=jp;j<JE-1;j++) {
	    m = j - jp;
	    daBzy(i,j) = da(m);
	    dbBzy(i,j) = db(m);
	}
	for(j=jp;j<JE-1;j++) {
	    m = j - jp;
	    caDx(i,j) = ca(m);
	    cbDx(i,j) = cb(m);
	}
    }
    for(j=0;j<JE-1;j++) { // Front & back
	for(i=0;i<NPML;i++) {
	    m = NPML - 1 - i;
	    daBzx(i,j) = da(m);
	    dbBzx(i,j) = db(m);
	}
	for(i=ip;i<IE-1;i++) {
	    m = i - ip;
	    daBzx(i,j) = da(m);
	    dbBzx(i,j) = db(m);
	}
	for(i=1;i<NPML+1;i++) {
	    m = NPML - i;
	    caDy(i,j) = ca(m);
	    cbDy(i,j) = cb(m);
	}
	for(i=ip;i<IE-1;i++) {
	    m = i - ip;
	    caDy(i,j) = ca(m);
	    cbDy(i,j) = cb(m);
	}
    }
    //***********************************************************************
    // Define cloak structure
    //***********************************************************************
    ic = IE/2;
    jc = JE/2;
    tan_l = 0.0;
    for(i=0;i<IE-1;i++) { // woLastI
	for(j=0;j<JE-1;j++) { //woLastJ
	    // X components
	    dist = sqrt(pow2(i+0.5-ic) + pow2(j-jc));
	    if(dist>R1 && dist<=R2) {   
		eps_theta = getEpsTheta(dist, R1, R2);
		eps_r = getEpsR(dist, R1, R2);
		gamma_px = tan_l*omega*eps_r/(1.0-eps_r);
		omega_px = sqrt(4.0*(eps_r-1.0)*(cos(omega*dt)-1.0)
				/(pow2(dt)*(cos(omega*dt)+1.0)));
		sinx = 1.0*(j-jc)/dist;
		cosx = 1.0*(i+0.5-ic)/dist;
		a0x(i,j) = eps_theta*(1.0/pow2(dt) + gamma_px/2.0/dt 
				      + pow2(omega_px)/4.0);
		a1x(i,j) = eps_theta*(-2.0/pow2(dt) + pow2(omega_px)/2.0);
		a2x(i,j) = eps_theta*(1.0/pow2(dt) - gamma_px/2.0/dt 
				      + pow2(omega_px)/4.0);
		b0xx(i,j) = pow2(sinx)*(1.0/pow2(dt) + gamma_px/2.0/dt + pow2(omega_px)/4.0) 
		    + eps_theta*pow2(cosx)*(1.0/pow2(dt) + gamma_px/2.0/dt);
		b1xx(i,j) = pow2(sinx)*(-2.0/pow2(dt) + pow2(omega_px)/2.0) 
		    - 2.0*eps_theta*pow2(cosx)/pow2(dt);
		b2xx(i,j) = pow2(sinx)*(1.0/pow2(dt) - gamma_px/2.0/dt + pow2(omega_px)/4.0) 
		    + eps_theta*pow2(cosx)*(1.0/pow2(dt) - gamma_px/2.0/dt);
		b0xy(i,j) = sinx*cosx*(eps_theta*(1.0/pow2(dt) + gamma_px/2.0/dt) 
				       - (1.0/pow2(dt) + gamma_px/2.0/dt + pow2(omega_px)/4.0));
		b1xy(i,j) = sinx*cosx*(eps_theta*(-2.0/pow2(dt)) 
				       - (-2.0/pow2(dt) + pow2(omega_px)/2.0));
		b2xy(i,j) = sinx*cosx*(eps_theta*(1.0/pow2(dt) - gamma_px/2.0/dt)
				       - (1.0/pow2(dt) - gamma_px/2.0/dt + pow2(omega_px)/4.0));
	    }
	    else if(dist<=R1) {
		cbDx(i,j) = 0.0;
	    }
	    // Y components
	    dist = sqrt(pow2(i-ic) + pow2(j+0.5-jc));
	    if(dist>R1 && dist<=R2) {
		eps_theta = getEpsTheta(dist, R1, R2);
		eps_r = getEpsR(dist, R1, R2);
		gamma_py = tan_l*omega*eps_r/(1.0-eps_r);
		omega_py = sqrt(4.0*(eps_r-1.0)*(cos(omega*dt)-1.0)
				/(pow2(dt)*(cos(omega*dt)+1.0)));
		siny = 1.0*(j+0.5-jc)/dist;
		cosy = 1.0*(i-ic)/dist;
		a0y(i,j) = eps_theta*(1.0/pow2(dt) + gamma_py/2.0/dt + pow2(omega_py)/4.0);
		a1y(i,j) = eps_theta*(-2.0/pow2(dt) + pow2(omega_py)/2.0);
		a2y(i,j) = eps_theta*(1.0/pow2(dt) - gamma_py/2.0/dt + pow2(omega_py)/4.0);
		b0yy(i,j) = pow2(cosy)*(1.0/pow2(dt) + gamma_py/2.0/dt + pow2(omega_py)/4.0) 
		    + eps_theta*pow2(siny)*(1.0/pow2(dt) + gamma_py/2.0/dt);
		b1yy(i,j) = pow2(cosy)*(-2.0/pow2(dt) + pow2(omega_py)/2.0)
		    - 2.0*eps_theta*pow2(siny)/pow2(dt);
		b2yy(i,j) = pow2(cosy)*(1.0/pow2(dt) - gamma_py/2.0/dt + pow2(omega_py)/4.0) 
		    + eps_theta*pow2(siny)*(1.0/pow2(dt) - gamma_py/2.0/dt);
		b0yx(i,j) = siny*cosy*(eps_theta*(1.0/pow2(dt) + gamma_py/2.0/dt) 
				       - (1.0/pow2(dt) + gamma_py/2.0/dt + pow2(omega_py)/4.0));
		b1yx(i,j) = siny*cosy*(eps_theta*(-2.0/pow2(dt)) 
				       - (-2.0/pow2(dt) + pow2(omega_py)/2.0));
		b2yx(i,j) = siny*cosy*(eps_theta*(1.0/pow2(dt) - gamma_py/2.0/dt) 
				       - (1.0/pow2(dt) - gamma_py/2.0/dt + pow2(omega_py)/4.0));
	    }
	    else if(dist<=R1) {
		cbDy(i,j) = 0.0;
	    }
	    // Z components
	    dist = sqrt(pow2(i-ic) + pow2(j-jc));
	    if(dist>R1 && dist<=R2) {
		A = getMuA(dist, R1, R2);
		mu_z = getMuZ(dist, R1, R2);
		if (tretyakov)
		    a0z(i,j) = mu_z;
		else {
		    gamma_pz = 0.0;
		    omega_pz = sqrt(4.0*(mu_z-1.0)*(cos(omega*dt)-1.0)/(pow2(dt)*(cos(omega*dt)+1.0)));
		    a0z(i,j) = A*(1.0/pow2(dt) + gamma_pz/2.0/dt + pow2(omega_pz)/4.0);
		    a1z(i,j) = A*(-2.0/pow2(dt) + pow2(omega_pz)/2.0);
		    a2z(i,j) = A*(1.0/pow2(dt) - gamma_pz/2.0/dt + pow2(omega_pz)/4.0);
		    b0z(i,j) = 1.0/pow2(dt) + gamma_pz/2.0/dt;
		    b1z(i,j) = -2.0/pow2(dt);
		    b2z(i,j) = 1.0/pow2(dt) - gamma_pz/2.0/dt;
		}
	    }
	}
    }
    //***********************************************************************
    // Begin time-stepping loop
    //***********************************************************************
    std::cout<<"Total n = "<<NMAX<<" , running :"<<std::endl;
    for(n=0;n<NMAX;n++) {
	if (n % 20 == 0) std::cout<<n<<" "<<std::flush;
	/* Save past field*/
	Dx_h2 = Dx; cycleArrays(Dx_h2, Dx_h1, Dx);  Dy_h2 = Dy; cycleArrays(Dy_h2, Dy_h1, Dy);
	Ex_h2 = Ex; cycleArrays(Ex_h2, Ex_h1, Ex);  Ey_h2 = Ey; cycleArrays(Ey_h2, Ey_h1, Ey);
	Bz_h2 = Bz; cycleArrays(Bz_h2, Bz_h1, Bz);  Hz_h2 = Hz; cycleArrays(Hz_h2, Hz_h1, Hz);
	/* Update D field */
	Dx(woLastI,woEdgeJ)= caDx(woLastI,woEdgeJ)*Dx(woLastI,woEdgeJ) 
	    + cbDx(woLastI,woEdgeJ)*(Hz(woLastI,woEdgeJ) - Hz(woLastI,woEdgeJ-1));
	Dy(woEdgeI,woLastJ) = caDy(woEdgeI,woLastJ)*Dy(woEdgeI,woLastJ) 
	    + cbDy(woEdgeI,woLastJ)*(Hz(woEdgeI-1,woLastJ) - Hz(woEdgeI,woLastJ));
	/* Averaged D field */
	applyStencil(averX(), Dx, Dx_a);
	applyStencil(averX(), Dx_h1, Dx_a_h1);
	applyStencil(averX(), Dx_h2, Dx_a_h2);
	applyStencil(averY(), Dy, Dy_a);
	applyStencil(averY(), Dy_h1, Dy_a_h1);
	applyStencil(averY(), Dy_h2, Dy_a_h2);
	
	/* Update E field */
	Ex(woLastI,woEdgeJ) = (b0xx(woLastI,woEdgeJ)*Dx(woLastI,woEdgeJ) 
			       + b1xx(woLastI,woEdgeJ)*Dx_h1(woLastI,woEdgeJ) 
			       + b2xx(woLastI,woEdgeJ)*Dx_h2(woLastI,woEdgeJ)
			       + b0xy(woLastI,woEdgeJ)*Dy_a(woLastI,woEdgeJ) 
			       + b1xy(woLastI,woEdgeJ)*Dy_a_h1(woLastI,woEdgeJ) 
			       + b2xy(woLastI,woEdgeJ)*Dy_a_h2(woLastI,woEdgeJ)
			       - (a1x(woLastI,woEdgeJ)*Ex_h1(woLastI,woEdgeJ) 
				  + a2x(woLastI,woEdgeJ)*Ex_h2(woLastI,woEdgeJ)))
	    /a0x(woLastI,woEdgeJ);
	Ey(woEdgeI,woLastJ) = (b0yy(woEdgeI,woLastJ)*Dy(woEdgeI,woLastJ) 
			       + b1yy(woEdgeI,woLastJ)*Dy_h1(woEdgeI,woLastJ) 
			       + b2yy(woEdgeI,woLastJ)*Dy_h2(woEdgeI,woLastJ)
			       + b0yx(woEdgeI,woLastJ)*Dx_a(woEdgeI,woLastJ) 
			       + b1yx(woEdgeI,woLastJ)*Dx_a_h1(woEdgeI,woLastJ) 
			       + b2yx(woEdgeI,woLastJ)*Dx_a_h2(woEdgeI,woLastJ)
			       - (a1y(woEdgeI,woLastJ)*Ey_h1(woEdgeI,woLastJ) 
				  + a2y(woEdgeI,woLastJ)*Ey_h2(woEdgeI,woLastJ)))
	    /a0y(woEdgeI,woLastJ);
	/* Update B & H field */
	Bzx(woLastI,woLastJ) = daBzx(woLastI,woLastJ)*Bzx(woLastI,woLastJ) 
	    + dbBzx(woLastI,woLastJ)*(Ey(woLastI,woLastJ) - Ey(woLastI+1,woLastJ));
	Bzy(woLastI,woLastJ) = daBzy(woLastI,woLastJ)*Bzy(woLastI,woLastJ) 
	    + dbBzy(woLastI,woLastJ)*(Ex(woLastI,woLastJ+1) - Ex(woLastI,woLastJ));
	Bz = Bzx + Bzy;
	Hz = (b0z*Bz + b1z*Bz_h1 + b2z*Bz_h2 - (a1z*Hz_h1 + a2z*Hz_h2))/ a0z;
	/* Source excitation */
	for (int i =0; i<ie; ++i){
	    Hz(i,js) = Hz(i,js) + source(n);
	}
	/* Calculate real & imaginary parts */
	if(n>NMAX-N) {
	    double cosX = cos(2*PI*freq*dt*n); double sinX = sin(2*PI*freq*dt*n);
	    Ex_r = Ex_r + Ex*cosX; Ex_i = Ex_i - Ex*sinX;  Ey_r = Ey_r + Ey*cosX;
	    Ey_i = Ey_i - Ey*sinX; Hz_r = Hz_r + Hz*cosX;  Hz_i = Hz_i - Hz*sinX;
	}
	/* Calculate errors for convergence */
	if(n % N==0) {
	    if(max(abs(Hz_p))==0) { err(ne) = 100.0; }
	    else { err(ne) = max(abs(abs(Hz) - abs(Hz_p))) / max(Hz_p) * 100.0; }
	    ne++;
	    Hz_p(woLastI,woLastJ) = Hz(woLastI,woLastJ);
	}
    }     // End of time-stepping loop 
    Ex_r = Ex_r*2/N; Ex_i = Ex_i*2/N;
    Ey_r = Ey_r*2/N; Ey_i = Ey_i*2/N;
    Hz_r = Hz_r*2/N; Hz_i = Hz_i*2/N;
    /* Save data to file */
    print2DToFile("Hz_src.field2D", Hz);
    print2DToFile("Hz_r_src.field2D", Hz_r);
    print2DToFile("Hz_i_src.field2D", Hz_i);
    print2DToFile("Ex_src.field2D", Ex);
    print2DToFile("Ex_r_src.field2D", Ex_r);
    print2DToFile("Ex_i_src.field2D", Ex_i);
    print2DToFile("Ey_src.field2D", Ey);
    print2DToFile("Ey_r_src.field2D", Ey_r);
    print2DToFile("Ey_i_src.field2D", Ey_i);

    FILE *fp;
    fp = fopen("err_src.field1D", "w");
    for(n=0;n<round(NMAX/N);n++) {
	fprintf(fp, "%.4f\n", err(n));
    }
    fclose(fp);

    clock_t t2 = clock();
    std::cout<<std::endl<<"*** Finish! ***"<<std::endl;
    printf("%.2lf seconds of processing\n", (t2-t1)/(double)CLOCKS_PER_SEC);
    return 0;
}
void print2DToFile(string fileName, Array<double,2> inArray){
    FILE *fp;
    double maxx = abs(max(inArray));
    double minn = abs(min(inArray));
    double big = maxx>minn ? maxx:minn;
    inArray(0,0) = big;
    inArray(0,1) = -big;
    fp = fopen(fileName.c_str(), "w");
    for(int j=inArray.lbound(secondDim);j<=inArray.ubound(secondDim); ++j) {
	for(int i=inArray.lbound(firstDim);i<=inArray.ubound(firstDim); ++i) {
	    fprintf(fp, "%.8f ", inArray(i,j));
	}
	fprintf(fp, "\n");
    }
    fclose(fp);
}
double getEpsR(double dist, double innerR, double outerR){
    if (pendry) return 1.0*(dist-innerR)/dist;
    if (cai) return 1.0*pow2(outerR/(outerR-innerR))*pow2((dist-innerR)/dist);
    if (tretyakov) return 1.0*(outerR/(outerR-innerR))*pow2((dist-innerR)/dist);
    return 1.0;
};
// double getEpsPhi(double dist, double innerR, double outerR){
//     if (pendry) return 1.0*dist/(dist-innerR);
//     if (tretyakov) return 1.0*outerR/(outerR-innerR);
//     return 1.0;
// };
double getEpsTheta(double dist, double innerR, double outerR){
    if (pendry) return 1.0*dist/(dist-innerR);
    if (cai) return 1.0*pow2(outerR/(outerR-innerR));
    if (belov||tretyakov) return 1.0*outerR/(outerR-innerR);
    return 1.0;
};

double getMuZ(double dist, double innerR, double outerR){
    if (zhao) return 1.0*(outerR/(outerR-innerR))*(dist-innerR)/dist;
    if (pendry) return 1.0*pow2(outerR/(outerR-innerR))*(dist-innerR)/dist;
    if (cai) return 1.0;
    if (tretyakov) return 1.0*outerR/(outerR-innerR);
    return 1.0;
};
double getMuA(double dist, double innerR, double outerR){
    return 1.0*outerR/(outerR-innerR);
};
