/* Function to evaluate the Huygens principle by pure integration.
 * 
 *    [E] = window_integ(window,trans,rec,lamda,[sum,intType,tol])
 *
 * Procedure will perform an integration accordint to Huygens principle
 * via the window front defined by "window" with a transmitter at 
 * coordinates "trans" and a receiver at "rec". The integral will be 
 * evaluated for all wave length in the vector lambda.
 * The function evaluated will be for each lambda_0
 * the scalar Kirchhoff diffraction.
 *
 * Please note: A normalistion of E on the free space loss 
 *              at the carrier frequency is NOT done!
 *
 * Input: window = [y_min; y_max; z_min; z_max] as definition of the 
 *                 aperture in space with x=0. Each coloumn defines a 
 *                 different window whose result will be superimposed
 *                 to the other ones. WINDOWS ARE COLOUMNWISE !!!
 *        trans  = [x,y,z] as definition of the transmitter position.
 *        rec    = [x,y,z] as definition of the receiver position. Multiple
 *                 receiver positions are handled as one receiver position per
 *                 coloumn.
 *        lambda = Wave length \lambda_0..\lambda_{N-1} the integral shall be evaluated on.
 *        sum = If sum==1 (default), then E will be sum of all window contributions,
 *              while for sum!=1, E is separated for each window.
 *        intType = Type of diffraction calculation
 *                   1 = Scalar Kirchhoff
 *                   2 = Exponentials only
 *                   3 = Advanced Fresnel, taking the cross elements also into account
 *                   4 = Pure Fresnel for checking Matlab code (does not work
 *                        multiple lambdas!!)
 *                   5 = Maggi-Rubinowicz boundary diffraction wave
 *                        (equivalent to Kirchhoff, but should be faster)
 *                   6 = Maggi-Rubinowicz according to Albani [1] (default)
 *                        (equivalent to Kirchhoff, but should be faster
 *                         and have a better numerical stability compared to 5)
 *                   7 = mmMagic approximate model (no numerical integration).
 *                   8 = Proposed update for the Fresnel integral method (no numerical integration).
 *                   9 = Fresnel diffraction approximation (equivalent to 4 but without numerical integration).
 *        tol = Tolerance to be used for the integration (default is 1e-6).
 * Output E = Electrical field as vector for each wave length given as 
 *            complex value.
 *             Matrix will be [receiver_positions x lambda x windows] if sum==0.
 *             Matrix will be [receiver_positions x lambda] if sum!=0.
 *
 * [1] M. Albani, "Boundary Diffracted Wave and Incremental Geometrical Optics:
 *     A Numerically Efficient and Physically Appealing Line-Integral Representation
 *     of Radiation Integrals. Aperture Scalar Case", IEEE Trans. on Ant. and Prop.,
 *     vol. 59, no. 2, Feb. 2011
 *
 * Copyright 2016 Thomas Jost, German Aerospace Center (DLR)
 *
 *  This file is part of the satellite-to-indoor channel simulator.
 * 
 *     The satellite-to-indoor channel simulator is free software: 
 *     you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 2 of the License, or
 *     (at your option) any later version.
 * 
 *     The satellite-to-indoor channel simulator 
 *     is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with the program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#include <stdio.h>
#include <complex>
#include <vector>
#include <math.h>
#include "mex.h"

#define VERSION 1.4

#define PI 3.14159265358979323846264338327950288419716939937510
#define DOUBLE_EPS 2.22044604925031e-16   // taking eps from Matlab
#define TOLERANCE 1e-6
#define TRUE 1

#define KIRCHHOFF 1
#define EXP 2
#define ADV_FRESNEL 3
#define FRESNEL 4
#define RUBI 5
#define ALBANI 6
#define MM_MAGIC 7
#define UP_FRESNEL 8
#define FRESNEL_NO_INT 9

#define INT_NULL 0
#define INT_2D 1
#define INT_BOUND 2

#define YDIM 1  // internal for Maggi-Rubinowicz to indicate integration direction
#define ZDIM 2

#define SQR(x) ((x)*(x))
#define SIGN(x) (((x)>=0) ? 1.0:-1.0)
 
// Simple swapping function for double
inline void swap(double &a, double &b)
{
 double c = a;
 a = b;
 b = c;
}

// sinc routine as sin(x)/x
inline double sinc(double x)
{
 double tol = 1e3*DOUBLE_EPS; // tolerance, where Maclaurin series should be used
 if (fabs(x)<=tol)
  return 1-1/6*x*x+1/120*x*x*x*x; // use Maclaurin series expansion around x==0
 else
  return sin(x)/x;
}


/* --------------- Routines to calculate the Fresnel integral --------------- */

 /* Using:
 * [1] K. D. Mielenz, "Computation of Fresnel Integrals II",
 *     Journal of Research of the NIST, vol. 105, pp. 589-590, 2000
 * [2] J. Boersma, "Computation of Fresnel integrals", Math. Computation,
 *     vol. 14, p. 380, 1960 
 */

#define NR_COEFFS 12
 /* coefficients according to [1], Table 1. */
const double fn[NR_COEFFS] = {0.318309844, 9.34626e-8, -0.09676631, 
                              0.000606222, 0.325539361, 0.325206461, 
                              -7.450551455, 32.20380908, -78.8035274, 
                              118.5343352, -102.4339798, 39.06207702};
const double gn[NR_COEFFS] = {0.0, 0.101321519, -4.07292e-5, -0.152068115,
		                            -0.046292605, 1.622793598, -5.199186089, 
                              7.477942354, -0.695291507, -15.10996796, 
                              22.28401942, -10.89968491};                       


/* Function to calculate the real and the imaginary part
 * of the Fresnel coefficient according to the interpolation
 * of Mielenz [1] in eq. (2d) using the coefficients in Table 1.
 *
 * Input: x = input for the Fresnel integral
 * Output: re = real part of the calculated Fresnel integral.
 *         im = imaginary part of the calculated Fresnel integral.
 *
 */
void mielenzInterpolation(double x, double& re, double& im)
{
 double sign = 1.0;
 if (fabs(x)<=1.6)
  mexErrMsgTxt("Milenz's interpolation only valid for |x|>1.6!");
 if (x<0.0)     /* check for negative value of x */
  sign = -1.0;
 x = fabs(x);   /* calculate using absolute value, i.e. positive */
 double f=fn[0]/x, g=gn[0]/x, x2n = 1.0;   /* init for n=0 */
 for (unsigned int n=1; n<NR_COEFFS; n++)   /* loop to calcualte f(x) and g(x) according to eq. (2d), [1] */
 {
  x2n = x2n*x*x;   /* x^{2n} */
  f = f+fn[n]/(x2n*x);   /* eq. (2d), [1] */
  g = g+gn[n]/(x2n*x);   /* eq. (2d), [1] */
 }
 re = sign*(0.5+f*sin(PI*x*x/2.0)-g*cos(PI*x*x/2.0));  /* eq. (1a), [1] */
 im = sign*(0.5-f*cos(PI*x*x/2.0)-g*sin(PI*x*x/2.0));  /* eq. (1b), [1] */ 
}


/* Function to calculate the real and the imaginary part
 * of the Fresnel coefficient according to the interpolation
 * of Boersma [2] adopted to the notation of [1] based on 
 * a Taylor series expansion.
 *
 * Input: x = input for the Fresnel integral
 * Output: re = real part of the calculated Fresnel integral.
 *         im = imaginary part of the calculated Fresnel integral.
 *
 */
void boersmaInterpolation(double x, double& re, double& im)
{
 if (fabs(x)>1.6)
  mexErrMsgTxt("Boersma's interpolation only valid for |x|<=1.6!");
 unsigned int maxIter = 11;    /* maximum number of iterations */
 double c = 1.0;         /* coefficient c0, eq. (3a) in [1] */
 double s = PI/6.0;      /* coefficient s0, eq. (3b) in [1] */
 double x4n = 1.0;
 re = c*x;                 /* c0*x^1 for n=0 */
 im = s*x*x*x;           /* s0*x^3 for n=0 */
 for (unsigned int n=1; n<=maxIter; n++)    
 {
  x4n = x4n*x*x*x*x;    /* x^{4n} */
  c = -PI*PI*(4.0*(n-1)+1.0)*c/(4.0*(2.0*(n-1)+1.0)*(2.0*(n-1)+2.0)*(4.0*(n-1)+5.0));  /* calculate c_{n} from c{n-1} according to eq. (3a) in [1] */
  s = -PI*PI*(4.0*(n-1)+3.0)*s/(4.0*(2.0*(n-1)+2.0)*(2.0*(n-1)+3.0)*(4.0*(n-1)+7.0));  /* calculate s_{n} from s{n-1} according to eq. (3b) in [1] */  
  re = re+c*x4n*x;        /* c(x) = \sum_{0}^{\infty} c_n x^{4n+1}   */
  im = im+s*x4n*x*x*x;    /* s(x) = \sum_{0}^{\infty} s_n x^{4n+3}   */
 } 
}


/* Function to calculate the Fresnel integral using the Boersma and Mielenz
 * interpolations. Retunr values are given as real part in re and imaginary 
 * part in im.
 */
void fcs(double x, double& re, double& im)
{
 if (fabs(x)<=1.6)
  boersmaInterpolation(x, re, im);  /* taking Taylor series expansion method */
 else
  mielenzInterpolation(x, re, im);
}


/* --------------- End of the Fresnel integral --------------- */


using namespace std;

struct Param;  // forward declaration for struct Param

/* --------------- Class Adapt --------------- */
// Implementation of class Adapt from Adap.h (Numerical Recipies in C++)
class Adapt 
{
 public:
	 double TOL,toler;
  bool terminate,out_of_tolerance;
  Adapt();
	 Adapt(double tol);
	 void integrate(void func(double,vector<complex<double> >&,void*), vector<complex<double> > &ret, 
                 const double a, const double b, void* param);
 private:
  void* par;   // parameter struct for function evaluation
	 void adaptlob(void func(double,vector<complex<double> >&,void*), vector<complex<double> > &ret,
                const double a, const double b, const vector<complex<double> > &fa, 
                const vector<complex<double> > &fb, const double is);
};

// Empty contructor, seting tolerance to standard value of 1e-6
Adapt::Adapt() : TOL(1e-6),terminate(true),out_of_tolerance(false),par(NULL)
{
	const double EPS=2.22044604925031e-016; // eps taken from Matlab
	if (TOL < 10.0*EPS)
		TOL=10.0*EPS;
}

// Constructor with tolerance as parameter
Adapt::Adapt(double tol) : TOL(tol),terminate(true),out_of_tolerance(false),par(NULL)
{
	const double EPS=2.22044604925031e-016; // eps taken from Matlab
	if (TOL < 10.0*EPS)
		TOL=10.0*EPS;
}

/* Function to perform the integration.
 * Input : func = function pointer to the expression to be integrated.
 *         ret = vector of complex double as return value.
 *         a,b = real valued integration boundaries.
 *         par = pointer to additional parameter for func.
 */
void Adapt::integrate(void func(double,vector<complex<double> >&,void*), vector<complex<double> > &ret, 
                      const double a, const double b, void* param)
{
 const double x[13] = {-1.0,-0.942882415695480,-sqrt(2.0/3.0),-0.641853342345781,-1.0/sqrt(5.0),
                       -0.236383199662150,0.0,0.236383199662150,1.0/sqrt(5.0),0.641853342345781,
                       sqrt(2.0/3.0),0.942882415695480,1.0};  // x-values for 13-point Kronrod rule with x in [-1,1]
 const unsigned int x4GL[4] = {0,4,8,12};   // points used by 4-point Gauss-Lobatto rule
 const unsigned int x7Kron[7] = {0,2,4,6,8,10,12};  // points used by 7-point Kronrod rule
 const double w4GL[4] = {1.0/12.0, 5.0/12.0, 5.0/12.0, 1.0/12.0}; // weights for 4-point Gauss-Lobatto rule
 const double w7Kron[7] = {77.0/2940.0, 432.0/2940.0, 625.0/2940.0, 672.0/2940.0, 625.0/2940.0,
                           432.0/2940.0, 77.0/2940.0};   // weights for 7-point Kronrod rule
 const double w13Kron[13] = {+0.0158271919734802/2.0, 
                             +0.0942738402188500/2.0,
                             +0.155071987336585/2.0,
                             +0.188821573960182/2.0,
                             +0.199773405226859/2.0,
                             +0.224926465333340/2.0,
                             +0.242611071901408/2.0,
                             +0.224926465333340/2.0,
                             +0.199773405226859/2.0,
                             +0.188821573960182/2.0,
                             +0.155071987336585/2.0,
                             +0.0942738402188500/2.0,
                             +0.0158271919734802/2.0};   // weights for 13-point Kronrod rule
 vector<complex<double> > y[13];  // variable for function evaluations
 unsigned int vec_size=0;
 double m = (b+a)/2.0, h = (b-a)/2.0, wh = (b-a);  // pre-factors etc.
 vector<complex<double> > i4GL,i7Kron,i13Kron;         // init integral values
 double erri1=0.0,erri2=0.0,r,tol=TOL,is=0.0;
 par = param;     // save parameter struct for later use in adaptlob
	for (unsigned int i=0; i<13; i++)
		func(m+h*x[i],y[i],par);            // fill up the array with vectors from function evaluation
 vec_size = (unsigned int) y[0].size();   // size of the vectors
 i4GL.resize(vec_size);   // resize vectors to correct size
 i7Kron.resize(vec_size);
 i13Kron.resize(vec_size);
 for (unsigned int l=0; l<vec_size; l++)  // go through each dimension of the output vector
 {
  i4GL[l] = 0.0;
  for (unsigned int i=0; i<4; i++)
   i4GL[l] += wh*w4GL[i]*y[x4GL[i]][l];        // get integral for 4-point Gauss-Lobatto rule
  i7Kron[l] = 0.0;
  for (unsigned int i=0; i<7; i++)
   i7Kron[l] += wh*w7Kron[i]*y[x7Kron[i]][l];  // get integral for 7-point Kronrod rule
  i13Kron[l] = 0.0;
  for (unsigned int i=0; i<13; i++)
   i13Kron[l] += wh*w13Kron[i]*y[i][l];        // get integral for 13-point Kronrod rule
  erri1 += abs(i7Kron[l]-i13Kron[l]);    // calculate as the norm(x,inf)=sum(abs(x))
  erri2 += abs(i4GL[l]-i13Kron[l]);  // calculate as the norm(x,inf)=sum(abs(x))
  is += abs(i13Kron[l]);        // calculate as the norm(x,inf)=sum(abs(x))
 } 
 r = (erri2!=0.0) ? (erri1/erri2) : (1.0);
 toler = ((r>0.0) && (r<1.0)) ? (TOL/r) : (TOL); // create tolerance value
 //is = is*toler/DOUBLE_EPS;   // first estimate for the integral
 if (is==0.0)
  is = wh;
 adaptlob(func,ret,a,b,y[0],y[12],is);   // perform integration
}

/* Internal function to perform the integration.
 * Input : func = function pointer to the expression to be integrated.
 *         ret = vector of complex double as return value.
 *         a,b = real valued integration boundaries.
 *         fa,fb = evalueated func at points a and b.
 *         is = first error estimate from integrate() for the integral.
 */
void Adapt::adaptlob(void func(double,vector<complex<double> >&,void*), vector<complex<double> > &ret,
                     const double a, const double b, const vector<complex<double> > &fa, 
                     const vector<complex<double> > &fb, const double is)
{
 const double x7Kron[7] = {-1.0,-sqrt(2.0/3.0),-1.0/sqrt(5.0),
                            0.0,1.0/sqrt(5.0),sqrt(2.0/3.0),1.0};  // x-values for 7-point Kronrod rule with x in [-1,1]
 const unsigned int x4GL[4] = {0,2,4,6};   // points used by 4-point Gauss-Lobatto rule
 const double w4GL[4] = {1.0/12.0, 5.0/12.0, 5.0/12.0, 1.0/12.0}; // weights for 4-point Gauss-Lobatto rule
 const double w7Kron[7] = {77.0/2940.0, 432.0/2940.0, 625.0/2940.0, 672.0/2940.0, 625.0/2940.0,
                           432.0/2940.0, 77.0/2940.0};   // weights for 7-point Kronrod rule
 vector<complex<double> > y[7];
 unsigned int vec_size = (unsigned int) fa.size();
 vector<complex<double> > i4GL(vec_size), i7Kron(vec_size);
 double xval[7], err_sum=0.0;                 // abcissa values
 if (ret.size()!=vec_size)
  ret.resize(vec_size);  // resize output vector if necessary         
 double m = (b+a)/2.0, h = (b-a)/2.0, wh = (b-a);  // pre-factors etc.
 if ((fabs(h)<=DOUBLE_EPS) || (m==a) || (m==b))
 {
  for (unsigned int i=0; i<vec_size; i++)
   ret[i] = h*(fa[i]+fb[i]);   // we got a singularity most propbably, so stop here
 }
 y[0] = fa;  y[6] = fb;   // points already evaluated     
 xval[0] = a;  xval[7] = b;
 for (unsigned int i=1; i<6; i++)
 {
  xval[i] = m+h*x7Kron[i];  // x-value within [a,b]
		func(xval[i],y[i],par);   // fill up the array with vectors from function evaluation
 }
 for (unsigned int l=0; l<vec_size; l++)
 {
  i4GL[l] = 0.0;
  for (unsigned int i=0; i<4; i++)
   i4GL[l] += wh*w4GL[i]*y[x4GL[i]][l];   // get integral for 4-point Gauss-Lobatto rule
  i7Kron[l] = 0.0;
  for (unsigned int i=0; i<7; i++)
   i7Kron[l] += wh*w7Kron[i]*y[i][l];     // get integral for 7-point Kronrod rule
  err_sum += abs(i4GL[l]-i7Kron[l])/vec_size;   // calculate as the norm(x,inf)/dim=mean(abs(x)) 
 }
 if ((is+err_sum <= is+toler) | (m<=a) | (b<=m))   // check if error is small enough, i.e. less than the tolerance wanted
 {
  if ((xval[1] <= a || b <= xval[5]) && terminate) 
  {
			out_of_tolerance = true;
   terminate = false;
		}
  ret = i7Kron;    // finished
 }
	else
 {
  vector<complex<double> > dum[6];  // dummy vectors for calculation
  adaptlob(func,dum[0],a,xval[1],fa,y[1],is);
		adaptlob(func,dum[1],xval[1],xval[2],y[1],y[2],is);
		adaptlob(func,dum[2],xval[2],xval[3],y[2],y[3],is);
		adaptlob(func,dum[3],xval[3],xval[4],y[3],y[4],is);
		adaptlob(func,dum[4],xval[4],xval[5],y[4],y[5],is);
		adaptlob(func,dum[5],xval[5],b,y[5],fb,is);
  for (unsigned int i=0; i<vec_size; i++)
  {
   ret[i] = dum[0][i]+dum[1][i]+dum[2][i]+
            dum[3][i]+dum[4][i]+dum[5][i];  // add all of them together for output
  }
 }        
}


/* -------- Integrating function --------- */

// Functions implementation to be integrated

// To verify use quad(@(x)sin(x).*cos(x)+2*j*x,0,2)
//  or dblquad(@(x,y)sin(x).*cos(x).*cos(y)+j*2*x.*sin(y),0,2,-1,5) for double integration
//  or dblquad(@(x,y)sin([0:5].'*x(:).').*cos(ones(6,1)*x(:).').*cos(ones(6,1)*y(:).')+j*2*(ones(6,1)*y(:).').*sin(ones(6,1)*y(:).'),0,2,-1,5,@quadv) for vector double integration

// Parameter struct for function 
struct Param
{
 double ymin, ymax; // integration from ymin...ymax
 double zmin, zmax; // integration from zmin...zmax
 double *trans;     // transmitter positon as [x,y,z]
 double *rec;       // receiver positon as [x,y,z]
 unsigned int ndim; // dimension of array lams
 double *lams;      // array of the wave lengths
 double n;          // direction toward receiver on x-axis (either + or -1)
 double ds[3];      // direction vector for boundary (used only in RUBI)
 unsigned int integDir; // integration over YDIM or ZDIM so element of {YDIM,ZDIM}
 unsigned intKind;  // type of integration to be used (KIRCHHOFF, EXP, ...)
 
 void (*f_inner)(double,vector<complex<double> >&,void*); // pointer to inner integration function
 double y_save;     // internal saving variable 
 Adapt int_save;    // internal integration instance
};


/* KIRCHHOFF diffraction calculation
 * inner function evaluates at y_save and z */
void f_inner_kirchhoff(double z, vector<complex<double> > &ret, void* param)
{
 //printf(" in kirchhoff !\n");
 double arg = 0.0;   // dummy variable
 Param *par = (Param*) param;
   // pre-calculate d = norm([0,y,z]-trans)+norm([0,y,z]-rec)
 double dt = sqrt(par->trans[0] * par->trans[0]+
                  (par->trans[1] - par->y_save)*(par->trans[1] - par->y_save)+
                  (par->trans[2] - z)*(par->trans[2] - z));     // norm([0,y,z]-trans)
 double dr = sqrt(par->rec[0] * par->rec[0]+
                  (par->rec[1] - par->y_save)*(par->rec[1] - par->y_save)+
                  (par->rec[2] - z)*(par->rec[2] - z));         // norm([0,y,z]-rec)
 if (ret.size()!=par->ndim)  // check for output size
  ret.resize(par->ndim);  // allocate correct memory
 for (unsigned int i=0; i<par->ndim; i++)
 {
  arg = 2*PI*(dt+dr)/par->lams[i];  // argument for sin and cos for exp(.)
  // evalute 
  // E=exp(-j*2*pi*d/lam)/(4*pi*dt*dr)*n*(r2*(1+j*k*|r2|)/|r2|^2-r1*(1+j*k*|r1|)/|r1|^2)
  // with r2=vector to receiver and r1=vector to transmitter
  ret[i] = complex<double>(cos(arg)/(4*PI*(dt*dr)),-sin(arg)/(4*PI*dt*dr))*
           par->n*(par->rec[0]*complex<double>(1.0,2*PI*dr/par->lams[i])/(dr*dr)-
              par->trans[0]*complex<double>(1.0,2*PI*dt/par->lams[i])/(dt*dt));
  //printf("  at %i = %f, %f d=%f\n",i,cos(arg)/(4*PI*d),-sin(arg)/(4*PI*d),d);
 }
}


/* EXP diffraction calculation
 * inner function evaluates at y_save and z */
void f_inner_exp(double z, vector<complex<double> > &ret, void* param)
{
 //printf(" in exp !\n");
 double arg = 0.0;   // dummy variable
 Param *par = (Param*) param;
   // pre-calculate d = norm([0,y,z]-trans)+norm([0,y,z]-rec)
 double dt = sqrt(par->trans[0] * par->trans[0]+
                  (par->trans[1] - par->y_save)*(par->trans[1] - par->y_save)+
                  (par->trans[2] - z)*(par->trans[2] - z));     // norm([0,y,z]-trans)
 double dr = sqrt(par->rec[0] * par->rec[0]+
                  (par->rec[1] - par->y_save)*(par->rec[1] - par->y_save)+
                  (par->rec[2] - z)*(par->rec[2] - z));         // norm([0,y,z]-rec)
 double d = dt+dr;  // distance
 if (ret.size()!=par->ndim)  // check for output size
  ret.resize(par->ndim);  // allocate correct memory
 for (unsigned int i=0; i<par->ndim; i++)
 {
  arg = 2*PI*d/par->lams[i];  // argument for sin and cos for exp(.)
  // evalute E=exp(-j*2*pi*d/lam)
  ret[i] = complex<double>(cos(arg)/par->lams[i],-sin(arg)/par->lams[i]);  // E=exp(-j*2*pi*d/lam)/lam
 }
}


/* ADVANCED FRESNEL diffraction calculation
 * inner function evaluates at y_save and z.
 * Additional: This is not a very efficient implementation.
 * Several calculations could be calculated outside of the 
 * integral.*/
void f_inner_adv_fresnel(double za, vector<complex<double> > &ret, void* param)
{
 //printf(" in advanced Fresnel !\n");
 double k = 0.0;   // dummy variable, wavelength
 double arg = 0.0;  // dummy variable
 Param *par = (Param*) param;
 double y01 = -par->trans[0]*(par->rec[1]-par->trans[1])/(par->rec[0]-par->trans[0])+par->trans[1]; // origin y-coordinate
 double z01 = -par->trans[0]*(par->rec[2]-par->trans[2])/(par->rec[0]-par->trans[0])+par->trans[2]; // origin z-coordinate
 double p2 = sqrt(par->trans[0]*par->trans[0]+(par->trans[1]-y01)*(par->trans[1]-y01)+   // distance origin - transmitter
                 (par->trans[2]-z01)*(par->trans[2]-z01));
 double p1 = sqrt(par->rec[0]*par->rec[0]+(par->rec[1]-y01)*(par->rec[1]-y01)+           // distance origin - receiver
                 (par->rec[2]-z01)*(par->rec[2]-z01));
 double F = 2*p1*p2/(p1+p2);   // 2 * focal length (Foc)
 double yp1 = (par->rec[1]-y01)/p1; // y-value for normalised vector p1
 double zp1 = (par->rec[2]-z01)/p1; // z-value for normalised vector p1
  // now shift coordinate system also for the integration variables
 double y = par->y_save-y01;
 double z = za-z01;
  // initialise loop over lambda
 if (ret.size()!=par->ndim)  // check for output size
  ret.resize(par->ndim);  // allocate correct memory
 for (unsigned int i=0; i<par->ndim; i++)
 {
  k = 2*PI/par->lams[i];  // 2*pi/wavelength
  arg = k*(p1+p2+((1-yp1*yp1)*y*y-2*yp1*zp1*y*z+(1-zp1*zp1)*z*z)/F);
   // for testing
  //arg = k*(p1+p2+(y*y+z*z)/F);  // results into simple Fresnel
  // evalute E=exp(-j*k*(p1+p2))*exp(-j*(k/2*Foc)*((1-yp1^2)*y-2*yp1*zp1*y*z+(1-zp1^2)*z)
  //   E = exp(-j*k*(p1+p2+((1-yp1^2)*y-2*yp1*zp1*y*z+(1-zp1^2)*z)/(2*Foc))
  ret[i] = complex<double>(cos(arg)/par->lams[i],-sin(arg)/par->lams[i]);  
 }
}


/* FRESNEL diffraction calculation
 * inner function evaluates at y_save and z 
 * IMPLEMENTS THE SIMPLE FRESNEL EQUATION !! */
void f_inner_fresnel(double z, vector<complex<double> > &ret, void* param)
{
 //printf(" in fresnel !\n");
 double arg = 0.0;   // dummy variable
 Param *par = (Param*) param;
   // pre-calculate d = norm([0,y,z]-trans)+norm([0,y,z]-rec)
 double d = sqrt((par->trans[0]-par->rec[0])*(par->trans[0]-par->rec[0])+  
                  (par->trans[1]-par->rec[1])*(par->trans[1]-par->rec[1])+
                  (par->trans[2]-par->rec[2])*(par->trans[2]-par->rec[2])); // distance receiver <-> transmitter
 if (ret.size()!=par->ndim)  // check for output size
  ret.resize(par->ndim);  // allocate correct memory
 for (unsigned int i=0; i<par->ndim; i++)
 {
  // evalute E=exp(-j*2*pi*d/lam)*exp(-j*pi/2*(y^2+z^2))
  ret[i] = complex<double>(cos(PI/2*(par->y_save*par->y_save+z*z)+2*PI*d/par->lams[i]),
                          -sin(PI/2*(par->y_save*par->y_save+z*z)+2*PI*d/par->lams[i]));  // E=exp(-j*2*pi*d/lam-pi/2*(y^2+z^2))
 }
}


/* MAGGI-RUBINOWICZ diffraction calculation
 * inner function evaluates at y_save and z */
void f_inner_rubi(double z, vector<complex<double> > &ret, void* param)
{
 //printf(" in rubi !\n");
 double preCalc = 0.0;  // precalculated value before loop
 double arg = 0.0, crProd = 0.0;   // dummy variable
 Param *par = (Param*) param;
 double q[3] = {0.0,par->y_save,z};        // current point on the boundary
 if (par->integDir==YDIM)  // integration direction is via y
  swap(q[1],q[2]);  // y and z value need to be swapped
 double drho = sqrt(par->trans[0] * par->trans[0]+
                   (par->trans[1] - q[1])*(par->trans[1] - q[1])+
                   (par->trans[2] - q[2])*(par->trans[2] - q[2]));     // norm([0,y,z]-trans)
 double dr = sqrt(par->rec[0] * par->rec[0]+
                 (par->rec[1] - q[1])*(par->rec[1] - q[1])+
                 (par->rec[2] - q[2])*(par->rec[2] - q[2]));         // norm([0,y,z]-rec)
 if (ret.size()!=par->ndim)  // check for output size
  ret.resize(par->ndim);  // allocate correct memory
 // now calculate (r x rho)*ds 
 for (unsigned int i=0; i<3; i++)
 {        
  if ((par->ds[i]>DOUBLE_EPS) || (par->ds[i]<-DOUBLE_EPS)) // check for non-zero
  {
   crProd += ((q[(i+1)%3]-par->rec[(i+1)%3])*(q[(i+2)%3]-par->trans[(i+2)%3])-
              (q[(i+2)%3]-par->rec[(i+2)%3])*(q[(i+1)%3]-par->trans[(i+1)%3]))*par->ds[i]; // (r x rho)*ds
  }
 }    
  // precalculate the term  (r x rho)*ds/(norm(rho)*norm(r)*(norm(r)*norm(rho)+r.'*rho))
 preCalc = crProd/(drho*dr*(dr*drho+(q[0]-par->rec[0])*(q[0]-par->trans[0])+
                                    (q[1]-par->rec[1])*(q[1]-par->trans[1])+
                                    (q[2]-par->rec[2])*(q[2]-par->trans[2])));
 for (unsigned int i=0; i<par->ndim; i++)  // loop over all frequencies
 {
  arg = -2*PI*(drho+dr)/par->lams[i];  // argument for sin and cos for exp(.)
  // evalute 
  // E=(1/(norm(rho)*norm(r)))*exp(-1i*par.k*(norm(rho)+norm(r)))*r.'*cross(rho,par.ds)/(norm(r)*norm(rho)+r.'*rho);
  ret[i] = complex<double>(cos(arg),sin(arg))*preCalc;
//   printf("  at (%g, %g) = %g + i*%g  preCalc=%g ds=[%g,%g,%g]\n",par->y_save,z,real(ret[i]),imag(ret[i]),crProd,par->ds[0],par->ds[1],par->ds[2]);
 }
}


/* MAGGI-RUBINOWICZ diffraction calculation according to ALBANI TAP 2011
 * inner function evaluates at y_save and z */
void f_inner_albani(double z, vector<complex<double> > &ret, void* param)
{
 //printf(" in albani !\n");
 double preCalc = 0.0;  // precalculated value before loop
 double arg = 0.0, crProd = 0.0;   // dummy variable
 Param *par = (Param*) param;
 double q[3] = {0.0,par->y_save,z};   // current point on the boundary
 complex<double> j = complex<double>(0.0,1.0); // dummy j for formula
 if (par->integDir==YDIM)  // integration direction is via y
  swap(q[1],q[2]);  // y and z value need to be swapped
 //--------------------------------------------
 // Point q is on the boundary!
 // vector r is from q to receiver
 // vector rd is from transmitter to q
 // vector R is from transmitter to receiver 
 //--------------------------------------------  
 double drd = sqrt(par->trans[0] * par->trans[0]+
                  (par->trans[1] - q[1])*(par->trans[1] - q[1])+
                  (par->trans[2] - q[2])*(par->trans[2] - q[2]));     // norm([0,y,z]-trans)
 double dr = sqrt(par->rec[0] * par->rec[0]+
                 (par->rec[1] - q[1])*(par->rec[1] - q[1])+
                 (par->rec[2] - q[2])*(par->rec[2] - q[2]));          // norm([0,y,z]-rec)
 double dR = sqrt((par->trans[0]-par->rec[0])*(par->trans[0]-par->rec[0])+  
                  (par->trans[1]-par->rec[1])*(par->trans[1]-par->rec[1])+
                  (par->trans[2]-par->rec[2])*(par->trans[2]-par->rec[2])); // norm(rec-trans)
 double rxrd = (par->rec[0]-q[0])*(q[0]-par->trans[0])+
               (par->rec[1]-q[1])*(q[1]-par->trans[1])+
               (par->rec[2]-q[2])*(q[2]-par->trans[2]);   // r.'*rd without normalisation
 if (ret.size()!=par->ndim)  // check for output size
  ret.resize(par->ndim);  // allocate correct memory
 // now calculate (r x rd)*ds 
 for (unsigned int i=0; i<3; i++)
 {        
  if ((par->ds[i]>DOUBLE_EPS) || (par->ds[i]<-DOUBLE_EPS)) // check for non-zero
  {
   crProd += ((par->rec[(i+1)%3]-q[(i+1)%3])*(q[(i+2)%3]-par->trans[(i+2)%3])-
              (par->rec[(i+2)%3]-q[(i+2)%3])*(q[(i+1)%3]-par->trans[(i+1)%3]))*par->ds[i]; // (r x rd)*ds (not normalised)
  }
 }    
  // precalculate the term  (r x rho)*ds/(norm(r)*norm(rd))
 preCalc = crProd/(dr*drd); // the prefactor 1/(4*pi) is taken into account in f_boundary()
 for (unsigned int i=0; i<par->ndim; i++)  // loop over all frequencies
 {
  double k = 2*PI/par->lams[i];   // dummy wave factor
  arg = -k*(drd+dr+dR)/2.0;  // argument for sin and cos for exp(.)
  // evalute eq. (20) of the paper
  ret[i] = preCalc*(2.0*j*k*sinc(k*(dr+drd-dR)/2.0)*complex<double>(cos(arg),sin(arg))/(dr+drd+dR)+
                    complex<double>(cos(-k*dR),sin(-k*dR))/(dr*drd+rxrd));
//   printf("  at (%g, %g) = %g + i*%g  preCalc=%g ds=[%g,%g,%g]\n",par->y_save,z,real(ret[i]),imag(ret[i]),crProd,par->ds[0],par->ds[1],par->ds[2]);
 }
}




/* Function to transform the normally used coordinate system into
 * the one used within ITU-R P.526-13 where the origin is located
 * on the direct line between transmitter and receiver and the 
 * screen on the x/y plane at z=0. While the normal coordinate is 
 * defined by the screen located at x=0.
 * Input par = Struct according to Param with 
 *              par.trans as transmitter coordinates
 *              par.rec as reiver coordinates
 *              par.ymin, par.ymax, par.zmin, par.zmax to define the window
 * Output [xt,yt,zt] = transmitter coordinates
 *        [xr,yr,zr] = receiver coordinates
 *        [x1,x2,y1,y2] = window coordinates, now as [xmin, xmax, ymin, ymax]              
 */
void coorTransP526(Param &par, double *xt, double *yt, double *zt, 
                               double *xr, double *yr, double *zr, 
                   double *x1, double *x2, double *y1, double *y2)
{
 double xa_y = -par.trans[0]*(par.rec[1]-par.trans[1])/(par.rec[0]-par.trans[0])+par.trans[1];  // new origin (x=0), y-coordinate
 double xa_z = -par.trans[0]*(par.rec[2]-par.trans[2])/(par.rec[0]-par.trans[0])+par.trans[2];  // new origin (x=0), z-coordinate
 // printf("xa_y=%g, xa_z=%g \n",xa_y,xa_z);
 *x1 = par.ymin-xa_y;    *x2 = par.ymax-xa_y;
 *y1 = par.zmin-xa_z;    *y2 = par.zmax-xa_z;
 *xt = par.trans[1]-xa_y;    // transmitter position in the other coordinate system including shift for new origin
 *yt = par.trans[2]-xa_z;
 *zt = par.trans[0]*par.n;   // making sure, transmitter is at negative z-axis
 *xr = par.rec[1]-xa_y;        // receiver position in the other coordinate system including shift for new origin
 *yr = par.rec[2]-xa_z;
 *zr = par.rec[0]*par.n;
}


/* mmMagic approximate model as provided inside the project. 
 * Source code developed from the input document 3J/146 to ITU
 * August 2017 at Geneva.*/
void calcMM_Magic(vector<complex<double> > &ans, Param &par)
{
 double xt,yt,zt, xr,yr,zr, x1,x2,y1,y2;  // variables to define transmitter, receiver position and window
 double x11,x12,y21,y22,nu11,nu12,nu21,nu22;
 //complex<double> ph11,ph12,ph21,ph22,Ph;
 complex<double> ph11PhG,ph12PhG,ph21PhG,ph22PhG;        
 //printf(" in mmMagic routine !\n");
 coorTransP526(par, &xt,&yt,&zt, &xr,&yr,&zr, &x1,&x2,&y1,&y2);  // perform coordinate transformation
 //printf(" got trans = [%g,%g,%g], rec = [%g,%g,%g], wind=[%g,%g,%g,%g] \n",xt,yt,zt,xr,yr,zr,x1,x2,y1,y2);
 double phi11 = atan((y1-yr)/zr)-atan((y1-yt)/zt);  // eq. (73e-h) to calculate the diffraction angle for each edge
 double phi12 = atan((y2-yr)/zr)-atan((y2-yt)/zt);
 double phi21 = atan((x1-xr)/zr)-atan((x1-xt)/zt);
 double phi22 = atan((x2-xr)/zr)-atan((x2-xt)/zt);  
 double Dprojt_11 = sqrt(SQR(zt)+SQR(y1-yt));    // eq. (78a-h)
 double Dprojr_11 = sqrt(SQR(zr)+SQR(y1-yr));
 double Dprojt_12 = sqrt(SQR(zt)+SQR(y2-yt));
 double Dprojr_12 = sqrt(SQR(zr)+SQR(y2-yr));
 double Dprojt_21 = sqrt(SQR(zt)+SQR(x1-xt));
 double Dprojr_21 = sqrt(SQR(zr)+SQR(x1-xr));
 double Dprojt_22 = sqrt(SQR(zt)+SQR(x2-xt));
 double Dprojr_22 = sqrt(SQR(zr)+SQR(x2-xr)); 
 double rproj1 = sqrt(SQR(zr-zt)+SQR(yr-yt)); // eq. (77a)
 double rproj2 = sqrt(SQR(zr-zt)+SQR(xr-xt)); // eq. (77b)
 double r = sqrt(SQR(zr-zt)+SQR(xr-xt)+SQR(yr-yt)); // distance trans->rec r=norm(trans-rec) 
 if (fabs(SQR(Dprojr_11)-SQR(Dprojt_11))<DOUBLE_EPS)  // check for being zero
  x11 = (xt+xr)/2.0;
 else
  x11 = (xt*SQR(Dprojr_11)-xr*SQR(Dprojt_11)-Dprojt_11*Dprojr_11*(xt-xr))/(SQR(Dprojr_11)-SQR(Dprojt_11));
 if (fabs(SQR(Dprojr_12)-SQR(Dprojt_12))<DOUBLE_EPS)  // check for being zero
  x12 = (xt+xr)/2.0;
 else
  x12 = (xt*SQR(Dprojr_12)-xr*SQR(Dprojt_12)-Dprojt_12*Dprojr_12*(xt-xr))/(SQR(Dprojr_12)-SQR(Dprojt_12));
 if (fabs(SQR(Dprojr_21)-SQR(Dprojt_21))<DOUBLE_EPS)  // check for being zero
  y21 = (yt+yr)/2.0;
 else
  y21 = (yt*SQR(Dprojr_21)-yr*SQR(Dprojt_21)-Dprojt_21*Dprojr_21*(yt-yr))/(SQR(Dprojr_21)-SQR(Dprojt_21));
 if (fabs(SQR(Dprojr_22)-SQR(Dprojt_22))<DOUBLE_EPS)  // check for being zero
  y22 = (yt+yr)/2.0;
 else
  y22 = (yt*SQR(Dprojr_22)-yr*SQR(Dprojt_22)-Dprojt_22*Dprojr_22*(yt-yr))/(SQR(Dprojr_22)-SQR(Dprojt_22)); 
 double D11 = sqrt(SQR(zr)+SQR(yr-y1)+SQR(xr-x11))+sqrt(SQR(zt)+SQR(yt-y1)+SQR(xt-x11));
 double D12 = sqrt(SQR(zr)+SQR(yr-y2)+SQR(xr-x12))+sqrt(SQR(zt)+SQR(yt-y2)+SQR(xt-x12));
 double D21 = sqrt(SQR(zr)+SQR(yr-y21)+SQR(xr-x1))+sqrt(SQR(zt)+SQR(yt-y21)+SQR(xt-x1));
 double D22 = sqrt(SQR(zr)+SQR(yr-y22)+SQR(xr-x2))+sqrt(SQR(zt)+SQR(yt-y22)+SQR(xt-x2));
 // printf("x11=%g, x12=%g, y21=%g, y22=%g, D11=%g, D12=%g, D21=%g, D22=%g \n",x11,x12,y21,y22,D11,D12,D21,D22); 
//  printf("phi11=%g, phi12=%g, phi21=%g, phi22=%g \n",phi11,phi12,phi21,phi22);
//  printf("Dprojt_11=%g, Dprojr_11=%g, Dprojt_12=%g, Dprojr_12=%g, Dprojt_21=%g, Dprojr_21=%g, Dprojt_22=%g, Dprojr_22=%g\n",
//          Dprojt_11,Dprojr_11,Dprojt_12,Dprojr_12,Dprojt_21,Dprojr_21,Dprojt_22,Dprojr_22);
//  printf("rproj1=%g, rproj2=%g \n",rproj1,rproj2);
 for (unsigned int i=0; i<par.ndim; i++)  // go through each lambda, i.e. frequency
 {
//   ph11 = polar(1.0,-2.0*PI*D11/par.lams[i]);  // terms phij*Gij/Ph
//   ph12 = polar(1.0,-2.0*PI*D12/par.lams[i]);
//   ph21 = polar(1.0,-2.0*PI*D21/par.lams[i]);
//   ph22 = polar(1.0,-2.0*PI*D22/par.lams[i]); 
  //Ph = polar(1.0,-2.0*PI*r/par.lams[i]);
   // check for positive lengths as due to numerical errors it might be slightly negative
  nu11 = (Dprojt_11+Dprojr_11-rproj1>=0) ? (2.0*sqrt((Dprojt_11+Dprojr_11-rproj1)/par.lams[i])):(0);
  nu12 = (Dprojt_12+Dprojr_12-rproj1>=0) ? (2.0*sqrt((Dprojt_12+Dprojr_12-rproj1)/par.lams[i])):(0);
  nu21 = (Dprojt_21+Dprojr_21-rproj2>=0) ? (2.0*sqrt((Dprojt_21+Dprojr_21-rproj2)/par.lams[i])):(0);
  nu22 = (Dprojt_22+Dprojr_22-rproj2>=0) ? (2.0*sqrt((Dprojt_22+Dprojr_22-rproj2)/par.lams[i])):(0); 
  //printf("nu11=%g, nu12=%g, nu21=%g, nu22=%g \n",nu11,nu12,nu21,nu22);
  //printf("  dist = %g \n", Dprojt_21+Dprojr_21-rproj2);
//   double G11 = cos(phi11/2.0)*(0.5-atan(nu11*1.4)/PI);  // taking the proposed empirical parameter of 1.4
//   double G12 = cos(phi12/2.0)*(0.5-atan(nu12*1.4)/PI);
//   double G21 = cos(phi21/2.0)*(0.5-atan(nu21*1.4)/PI);
//   double G22 = cos(phi22/2.0)*(0.5-atan(nu22*1.4)/PI);
//   printf("G11=%g, G12=%g, G21=%g, G22=%g \n",G11,G12,G21,G22);
  ph11PhG = polar(cos(phi11/2.0)*(0.5-atan(nu11*1.4)/PI),-2.0*PI*(D11-r)/par.lams[i]);  // terms phij*Gij/Ph
  ph12PhG = polar(cos(phi12/2.0)*(0.5-atan(nu12*1.4)/PI),-2.0*PI*(D12-r)/par.lams[i]);
  ph21PhG = polar(cos(phi21/2.0)*(0.5-atan(nu21*1.4)/PI),-2.0*PI*(D21-r)/par.lams[i]);
  ph22PhG = polar(cos(phi22/2.0)*(0.5-atan(nu22*1.4)/PI),-2.0*PI*(D22-r)/par.lams[i]);
  //ans[i] = (SIGN(phi11)*(0.5-ph11*G11/Ph)-SIGN(phi12)*(0.5-ph12*G12/Ph))*
  //         (SIGN(phi21)*(0.5-ph21*G21/Ph)-SIGN(phi22)*(0.5-ph22*G22/Ph)); // final model normalized to LoS
  //ans[i] = (SIGN(phi11)*(0.5-ph11PhG)-SIGN(phi12)*(0.5-ph12PhG))*
  //         (SIGN(phi21)*(0.5-ph21PhG)-SIGN(phi22)*(0.5-ph22PhG)); // final model normalized to LoS
  ans[i] = (SIGN(phi11)*(0.5-ph11PhG)-SIGN(phi12)*(0.5-ph12PhG))*
           (SIGN(phi21)*(0.5-ph21PhG)-SIGN(phi22)*(0.5-ph22PhG))*polar(1.0/r,-2.0*PI*r/par.lams[i]); // final model not normalized to LoS, so including whole propagation
//   printf("ph11PhG=(%g,%g), ph12PhG=(%g,%g), ph21PhG=(%g,%g), ph22PhG=(%g,%g)\n",
//           real(ph11PhG),imag(ph11PhG),real(ph12PhG),imag(ph12PhG),real(ph21PhG),imag(ph21PhG),real(ph22PhG),imag(ph22PhG));
//   printf("abs(ph11PhG)=(%g,%g), abs(ph12PhG)=(%g,%g), abs(ph21PhG)=(%g,%g), abs(ph22PhG)=(%g,%g)\n",
//           abs(ph11PhG),arg(ph11PhG),abs(ph12PhG),arg(ph12PhG),abs(ph21PhG),arg(ph21PhG),abs(ph22PhG),arg(ph22PhG));
//   printf("ph11=%.9g, ph12=%.9g, ph21=%.9g, ph22=%.9g\n",
//           -2.0*PI*(D11-r)/par.lams[i],-2.0*PI*(D12-r)/par.lams[i],-2.0*PI*(D21-r)/par.lams[i],-2.0*PI*(D22-r)/par.lams[i]);
//   printf(" RX-TX = %g \n",r); 
  //printf(" Got final E = (%g, %g) \n",real(ans[i]),imag(ans[i]));
 }
}


/* Proposed updated Fresnel diffraction calculations as approximate model. 
 * Source code developed from the input document to ITU August 2017 at Geneva.*/
void calcUpFresnel(vector<complex<double> > &ans, Param &par)
{
 double xt,yt,zt, xr,yr,zr, x1,x2,y1,y2;  // variables to define transmitter, receiver position and window
 double Cx1,Cx2,Cy1,Cy2,Sx1,Sx2,Sy1,Sy2;
 //printf(" in updated Fresnel routine !\n");
 coorTransP526(par, &xt,&yt,&zt, &xr,&yr,&zr, &x1,&x2,&y1,&y2);  // perform coordinate transformation
 //printf(" got trans = [%g,%g,%g], rec = [%g,%g,%g], wind=[%g,%g,%g,%g] \n",xt,yt,zt,xr,yr,zr,x1,x2,y1,y2);
 double phi11 = atan((y1-yr)/zr)-atan((y1-yt)/zt);  // eq. (73e-h) to calculate the diffraction angle for each edge
 double phi12 = atan((y2-yr)/zr)-atan((y2-yt)/zt);
 double phi21 = atan((x1-xr)/zr)-atan((x1-xt)/zt);
 double phi22 = atan((x2-xr)/zr)-atan((x2-xt)/zt);
 double r = sqrt(SQR(zr-zt)+SQR(xr-xt)+SQR(yr-yt)); // distance trans->rec r=norm(trans-rec) 
 for (unsigned int i=0; i<par.ndim; i++)  // go through each lambda, i.e. frequency
 {
   // x1 edge:
  fcs(SIGN(x1)*sqrt(2.0*pow(fabs(x1),1.18)*pow(1.0/zr-1.0/zt,0.18)*pow(fabs(phi21),0.82)/par.lams[i]), Cx1, Sx1);
   // x2 edge:
  fcs(SIGN(x2)*sqrt(2.0*pow(fabs(x2),1.18)*pow(1.0/zr-1.0/zt,0.18)*pow(fabs(phi22),0.82)/par.lams[i]), Cx2, Sx2);
   // y1 edge:
  fcs(SIGN(y1)*sqrt(2.0*pow(fabs(y1),1.18)*pow(1.0/zr-1.0/zt,0.18)*pow(fabs(phi11),0.82)/par.lams[i]), Cy1, Sy1);
   // y2 edge:
  fcs(SIGN(y2)*sqrt(2.0*pow(fabs(y2),1.18)*pow(1.0/zr-1.0/zt,0.18)*pow(fabs(phi12),0.82)/par.lams[i]), Cy2, Sy2);
  //printf("Cx1=%g, Cx2=%g, Cy1=%g, Cy2=%g \n",Cx1,Cx2,Cy1,Cy2);
  //ans[i] = complex<double>(0.5*((Cx2-Cx1)*(Sy2-Sy1)+(Sx2-Sx1)*(Cy2-Cy1)), 
  //                         -0.5*((Sx2-Sx1)*(Sy2-Sy1)-(Cx2-Cx1)*(Cy2-Cy1)));  // final calculation normalized version
  ans[i] = complex<double>(0.5*((Cx2-Cx1)*(Sy2-Sy1)+(Sx2-Sx1)*(Cy2-Cy1)), 
                           -0.5*((Sx2-Sx1)*(Sy2-Sy1)-(Cx2-Cx1)*(Cy2-Cy1)))*polar(1.0/r,-2.0*PI*r/par.lams[i]);  // final calculation not normalized version
  //printf(" Got final E = (%g, %g) \n",real(ans[i]),imag(ans[i]));
 }
}


/* Fresnel diffraction calculations as approximate model. 
 * Method is equivalent to the other Fresnel except that no numerical integration is used.*/
void calcFresnelNoInt(vector<complex<double> > &ans, Param &par)
{
 double xt,yt,zt, xr,yr,zr, x1,x2,y1,y2;  // variables to define transmitter, receiver position and window
 double Cx1,Cx2,Cy1,Cy2,Sx1,Sx2,Sy1,Sy2;
 //printf(" in Fresnel no integration routine !\n");
 coorTransP526(par, &xt,&yt,&zt, &xr,&yr,&zr, &x1,&x2,&y1,&y2);  // perform coordinate transformation for simplicity
 //printf(" got trans = [%g,%g,%g], rec = [%g,%g,%g], wind=[%g,%g,%g,%g] \n",xt,yt,zt,xr,yr,zr,x1,x2,y1,y2);
 double r = sqrt(SQR(zr-zt)+SQR(xr-xt)+SQR(yr-yt)); // distance trans->rec r=norm(trans-rec) 
 double d1 = sqrt(SQR(zt)+SQR(xt)+SQR(yt)); // distance trans->origin 
 double d2 = sqrt(SQR(zr)+SQR(xr)+SQR(yr)); // distance rec->origin 
 for (unsigned int i=0; i<par.ndim; i++)  // go through each lambda, i.e. frequency
 {
   // x1 edge:
  fcs(sqrt((2.0/par.lams[i])*(1.0/d1+1.0/d2))*x1, Cx1, Sx1);
   // x2 edge:
  fcs(sqrt((2.0/par.lams[i])*(1.0/d1+1.0/d2))*x2, Cx2, Sx2);
   // y1 edge:
  fcs(sqrt((2.0/par.lams[i])*(1.0/d1+1.0/d2))*y1, Cy1, Sy1);
   // y2 edge:
  fcs(sqrt((2.0/par.lams[i])*(1.0/d1+1.0/d2))*y2, Cy2, Sy2);
  //printf("Cx1=%g, Cx2=%g, Cy1=%g, Cy2=%g \n",Cx1,Cx2,Cy1,Cy2);
  //ans[i] = complex<double>(0.5*((Cx2-Cx1)*(Sy2-Sy1)+(Sx2-Sx1)*(Cy2-Cy1)), 
  //                         0.5*((Sx2-Sx1)*(Sy2-Sy1)-(Cx2-Cx1)*(Cy2-Cy1)));  // final calculation normalized version
  ans[i] = complex<double>(0.5*((Cx2-Cx1)*(Sy2-Sy1)+(Sx2-Sx1)*(Cy2-Cy1)), 
                           -0.5*((Sx2-Sx1)*(Sy2-Sy1)-(Cx2-Cx1)*(Cy2-Cy1)))*polar(1.0/r,-2.0*PI*r/par.lams[i]);  // final calculation not normalized version
  //printf(" Got final E = (%g, %g) \n",real(ans[i]),imag(ans[i]));
 }
}



/* -------- Housekeeping functions --------- */


// outer function, do the integration on z with a given y, for 2d integration
// for methods KIRCHHOFF, EXP, ADV_FRESNEL, FRESNEL
void f_outer(double y, vector<complex<double> > &ret, void* param)
{
 Param *par = (Param*) param;
 par->y_save = y;
 par->int_save.integrate(*(par->f_inner),ret,par->zmin,par->zmax,par);   // perform integration in z-direction
//  for (unsigned int i=0; i<ret.size(); i++)
//   printf(" got ret[%i] = %f, %f\n",i,real(ret[i]),imag(ret[i]));
 if (par->int_save.out_of_tolerance)
  mexErrMsgTxt("Required tolerance not reached in f()!\n");
}


/* Function to check for the LoS component through the window
 * defined in the struct param.
 *  ret = Vector of complex<double> as return for each lambda.
 *  param = Struct of type Param for the parameters. */
void checkForLos(vector<complex<double> > &ret, void* param)
{
 Param *par = (Param*) param;  // cast param to correct type of struct
 double lam = -par->rec[0]/(par->trans[0]-par->rec[0]); // calculate lambda for line from receiver to transmitter at x=0
 double v_y = par->rec[1]+lam*(par->trans[1]-par->rec[1]);  // y-coordinate
 double v_z = par->rec[2]+lam*(par->trans[2]-par->rec[2]);  // z-coordinate
 double d = sqrt((par->trans[0]-par->rec[0])*(par->trans[0]-par->rec[0])+  
                 (par->trans[1]-par->rec[1])*(par->trans[1]-par->rec[1])+
                 (par->trans[2]-par->rec[2])*(par->trans[2]-par->rec[2])); // distance receiver <-> transmitter
 if ((v_y>=par->ymin) && (v_y<=par->ymax)   // check for y boundary
     && (v_z>=par->zmin) && (v_z<=par->zmax)) // check for z boundary
 {
  for (unsigned int i=0; i<par->ndim; i++)  // go through each lambda
   ret[i] = complex<double>(cos(-2*PI*d/par->lams[i]),sin(-2*PI*d/par->lams[i]))/d;
 }
 else
 {
  for (unsigned int i=0; i<par->ndim; i++)  // go through each lambda
   ret[i] = complex<double>(0.0,0.0);
 }
}


/* Integration method over the boundary of the window.
 *  ret = Vector of complex<double> as return for each lambda.
 *  s = Instance to the object Adapt for integration.
 *  param = Struct of type Param for the parameters. */
void f_boundary(vector<complex<double> > &ret, Adapt &s, void* param)
{
 // ret is a vector for one receiver, transmitter and one window,
 // but over all different lambdas...
 Param *par = (Param*) param;
 vector<complex<double> > dum;  // to store intermediate data
 double preFac = -1/(4*PI);     // pre-factor
 par->ds[0] = 0.0;    
  // --- direct component
 if (par->intKind==RUBI) 
  checkForLos(ret,param);  // check if we have a direct illumination (only for RUBI)
 else
 {
  for (unsigned int i=0; i<par->ndim; i++)
   ret[i] = complex<double>(0.0,0.0);   // initialise output vector for ALBANI
 }
  // --- first edge (at par->ymin along z-axis)
 par->ds[1] = par->n*0.0;  par->ds[2] = par->n*(1.0); // direction z
 par->integDir=ZDIM;  
 par->y_save = par->ymin;  
 s.integrate(*(par->f_inner),dum,par->zmin,par->zmax,par);  // perform vector integration, output in variable dum
 for (unsigned int i=0; i<par->ndim; i++)  // go through each lambda
  ret[i] += preFac*dum[i];   // add values
  // --- second edge (at zmax along y-axis)
 par->ds[1] = par->n*1.0;  par->ds[2] = par->n*(0.0); // direction y
 par->integDir=YDIM; 
 par->y_save = par->zmax;  
 s.integrate(*(par->f_inner),dum,par->ymin,par->ymax,par);  // perform vector integration, output in variable dum
 for (unsigned int i=0; i<par->ndim; i++)  // go through each lambda
  ret[i] += preFac*dum[i];   // add values
  // --- third edge (at ymax along (-z)-axis)
 par->ds[1] = par->n*0.0;  par->ds[2] = par->n*(-1.0); // direction -z
 par->integDir=ZDIM;
 par->y_save = par->ymax;  
 s.integrate(*(par->f_inner),dum,par->zmin,par->zmax,par);  // perform vector integration, output in variable dum
 for (unsigned int i=0; i<par->ndim; i++)  // go through each lambda
  ret[i] += preFac*dum[i];   // add values 
  // --- fourth edge (at ymin along (-y)-axis)
 par->ds[1] = par->n*(-1.0);  par->ds[2] = par->n*(0.0); // direction -y
 par->integDir=YDIM;
 par->y_save = par->zmin;  
 s.integrate(*(par->f_inner),dum,par->ymin,par->ymax,par);  // perform vector integration, output in variable dum
 for (unsigned int i=0; i<par->ndim; i++)  // go through each lambda
  ret[i] += preFac*dum[i];   // add values  
}


/* Function to calculate the multiplicative parameter
 * in front of the integral for KIRCHHOFF, EXP and FRESNEL. */
complex<double> getMultFac(unsigned int intKind, void *param)
{
 Param *par = (Param*) param;  // get parameter struct
 if ((intKind==KIRCHHOFF) || (intKind==RUBI) || (intKind==ALBANI) || 
     (intKind==MM_MAGIC) || (intKind==UP_FRESNEL) || (intKind==FRESNEL_NO_INT))
  return complex<double>(1.0,0.0); 
 else if(intKind==EXP)
 {
  double norm_c1 = sqrt(par->trans[0]*par->trans[0]+   // length of vector
           (par->trans[1]-(par->ymax+par->ymin)/2)*(par->trans[1]-(par->ymax+par->ymin)/2)+
           (par->trans[2]-(par->zmax+par->zmin)/2)*(par->trans[2]-(par->zmax+par->zmin)/2)); 
  double norm_c2 = sqrt(par->rec[0]*par->rec[0]+
           (par->rec[1]-(par->ymax+par->ymin)/2)*(par->rec[1]-(par->ymax+par->ymin)/2)+
           (par->rec[2]-(par->zmax+par->zmin)/2)*(par->rec[2]-(par->zmax+par->zmin)/2));
  double c1 = par->trans[0]/norm_c1;  // normalised x-coordinate
  double c2 = par->rec[0]/norm_c2;
  return complex<double>(0.0,0.5*par->n*(c2-c1)/(norm_c1*norm_c2));
 }
 else if(intKind==ADV_FRESNEL)
 {
  double norm_c1 = sqrt(par->trans[0]*par->trans[0]+   // length of vector
           (par->trans[1]-(par->ymax+par->ymin)/2)*(par->trans[1]-(par->ymax+par->ymin)/2)+
           (par->trans[2]-(par->zmax+par->zmin)/2)*(par->trans[2]-(par->zmax+par->zmin)/2)); 
  double norm_c2 = sqrt(par->rec[0]*par->rec[0]+
           (par->rec[1]-(par->ymax+par->ymin)/2)*(par->rec[1]-(par->ymax+par->ymin)/2)+
           (par->rec[2]-(par->zmax+par->zmin)/2)*(par->rec[2]-(par->zmax+par->zmin)/2));
  double c1 = par->trans[0]/norm_c1;  // normalised x-coordinate
  double c2 = par->rec[0]/norm_c2;
  return complex<double>(0.0,0.5*par->n*(c2-c1)/(norm_c1*norm_c2));
 }
 else if(intKind==FRESNEL)
 {
  double d = sqrt((par->trans[0]-par->rec[0])*(par->trans[0]-par->rec[0])+  
                  (par->trans[1]-par->rec[1])*(par->trans[1]-par->rec[1])+
                  (par->trans[2]-par->rec[2])*(par->trans[2]-par->rec[2])); // distance receiver <-> transmitter
  return complex<double>(0.0,0.5/d);  // pre-factor as j/(2*d) in front of integral
 } 
 else
  mexErrMsgTxt("Impossible to be here. Parameter integration type is wrong! intKind value out of range!");
}


/* -------- Main Matlab wrapper --------- */

// Function to output the help of the routine
void printHelp()
{
 printf(" Function to evaluate the Huygens principle by pure integration.\n");
 printf("\n");
 printf("    [E] = window_integ(window,trans,rec,lamda,[sum,intType,tol])\n");
 printf("\n");
 printf(" Procedure will perform an integration accordint to Huygens principle\n");
 printf(" via the window front defined by \"window\" with a transmitter at \n");
 printf(" coordinates \"trans\" and a receiver at \"rec\". The integral will be \n");
 printf(" evaluated for all wave length in the vector lambda.\n");
 printf(" The function evaluated will be for each lambda_0\n");
 printf(" the scalar Kirchhoff diffraction.\n");
 printf("\n");
 printf(" Please note: A normalistion of E on the free space loss \n");
 printf("              at the carrier frequency is NOT done!\n");
 printf("\n");
 printf(" Input: window = [y_min; y_max; z_min; z_max] as definition of the \n");
 printf("                 aperture in space with x=0. Each coloumn defines a \n");
 printf("                 different window whose result will be superimposed\n");
 printf("                 to the other ones.  WINDOWS ARE COLOUMNWISE !!!\n");
 printf("        trans  = [x,y,z] as definition of the transmitter position.\n");
 printf("        rec    = [x,y,z] as definition of the receiver position. Multiple\n");
 printf("                 receiver positions are handled as one receiver position per\n");
 printf("                 coloumn.\n");
 printf("        lambda = Wave length \\lambda_0..\\lambda_{N-1} the integral shall be evaluated on.\n");
 printf("        sum = If sum==1 (default), then E will be sum of all window contributions,\n");
 printf("              while for sum!=1, E is separated for each window.\n");
 printf("        tol = Tolerance to be used for the integration (default is 1e-6).\n");
 printf("        intType = Type of diffraction calculation\n");
 printf("                   1 = Scalar Kirchhoff\n");
 printf("                   2 = Exponentials only\n");
 printf("                   3 = Advanced Fresnel, taking the cross elements also into account\n");
 printf("                   4 = Pure Fresnel for checking Matlab code (does not work\n");
 printf("                        multiple lambdas!!)\n");
 printf("                   5 = Maggi-Rubinowicz boundary diffraction wave\n");
 printf("                        (equivalent to Kirchhoff, but should be faster)\n");
 printf("                   6 = Maggi-Rubinowicz according to Albani [1] (default)\n");
 printf("                        (equivalent to Kirchhoff, but should be faster\n");
 printf("                         and have a better numerical stability compared to 5)\n");
 printf("                   7 = mmMagic approximate model (no numerical integration).\n");
 printf("                   8 = Proposed update for the Fresnel integral method (no numerical integration).\n");
 printf("                   9 = Fresnel diffraction approximation (equivalent to 4 but without numerical integration).\n");
 printf(" Output E = Electrical field as vector for each wave length given as \n");
 printf("            complex value.\n");
 printf("             Matrix will be [receiver_positions x lambda x windows] if sum==0.\n");
 printf("             Matrix will be [receiver_positions x lambda] if sum!=0.\n\n");
 printf(" [1] M. Albani, \"Boundary Diffracted Wave and Incremental Geometrical Optics:\n");
 printf("     A Numerically Efficient and Physically Appealing Line-Integral Representation\n");
 printf("     of Radiation Integrals. Aperture Scalar Case\", IEEE Trans. on Ant. and Prop.,\n");
 printf("     vol. 59, no. 2, Feb. 2011\n\n");
 printf(" This is version %g from %s of file %s.\n",VERSION,__DATE__,__FILE__);
 printf("\n Copyright 2016 Thomas Jost, German Aerospace Center (DLR)\n\n");
 printf("   This file is part of the satellite-to-indoor channel simulator.\n");
 printf("  \n");
 printf("      The satellite-to-indoor channel simulator is free software: \n");
 printf("      you can redistribute it and/or modify\n");
 printf("      it under the terms of the GNU General Public License as published by\n");
 printf("      the Free Software Foundation, either version 2 of the License, or\n");
 printf("      (at your option) any later version.\n");
 printf("  \n");
 printf("      The satellite-to-indoor channel simulator \n");
 printf("      is distributed in the hope that it will be useful,\n");
 printf("      but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
 printf("      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n");
 printf("      GNU General Public License for more details.\n");
 printf("  \n");
 printf("      You should have received a copy of the GNU General Public License\n");
 printf("      along with the program.  If not, see <http://www.gnu.org/licenses/>.\n");
}


/* Function to evaluate the Huygens principle by pure integration.
 * 
 *    [E] = window_integ(window,trans,rec,lamda,[sum,intType,tol])
 *
 * Procedure will perform an integration accordint to Huygens principle
 * via the window front defined by "window" with a transmitter at 
 * coordinates "trans" and a receiver at "rec". The integral will be 
 * evaluated for all wave length in the vector lambda.
 * The function evaluated will be for each lambda_0
 * the scalar Kirchhoff diffraction.
 *
 * Please note: A normalistion of E on the free space loss 
 *              at the carrier frequency is NOT done!
 *
 * Input: window = [y_min; y_max; z_min; z_max] as definition of the 
 *                 aperture in space with x=0. Each coloumn defines a 
 *                 different window whose result will be superimposed
 *                 to the other ones. WINDOWS ARE COLOUMNWISE !!!
 *        trans  = [x,y,z] as definition of the transmitter position.
 *        rec    = [x,y,z] as definition of the receiver position. Multiple
 *                 receiver positions are handled as one receiver position per
 *                 coloumn.
 *        lambda = Wave length \lambda_0..\lambda_{N-1} the integral shall be evaluated on.
 *        sum = If sum==1 (default), then E will be sum of all window contributions,
 *              while for sum!=1, E is separated for each window.
 *        intType = Type of diffraction calculation
 *                   1 = Scalar Kirchhoff
 *                   2 = Exponentials only
 *                   3 = Advanced Fresnel, taking the cross elements also into account
 *                   4 = Pure Fresnel for checking Matlab code (does not work
 *                        multiple lambdas!!)
 *                   5 = Maggi-Rubinowicz boundary diffraction wave
 *                        (equivalent to Kirchhoff, but should be faster)
 *                   6 = Maggi-Rubinowicz according to Albani [1] (default)
 *                        (equivalent to Kirchhoff, but should be faster
 *                         and have a better numerical stability compared to 5)
 *                   7 = mmMagic approximate model (no numerical integration).
 *                   8 = Proposed update for the Fresnel integral method (no numerical integration).
 *                   9 = Fresnel diffraction approximation (equivalent to 4 but without numerical integration).
 *        tol = Tolerance to be used for the integration (default is 1e-6).
 * Output E = Electrical field as vector for each wave length given as 
 *            complex value.
 *             Matrix will be [receiver_positions x lambda x windows] if sum==0.
 *             Matrix will be [receiver_positions x lambda] if sum!=0.
 *
 * [1] M. Albani, "Boundary Diffracted Wave and Incremental Geometrical Optics:
 *     A Numerically Efficient and Physically Appealing Line-Integral Representation
 *     of Radiation Integrals. Aperture Scalar Case", IEEE Trans. on Ant. and Prop.,
 *     vol. 59, no. 2, Feb. 2011
 */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
 Param par;  // parameter struct for integrated function
 double tol=TOLERANCE;  // tolerance value for the integration
 unsigned int summing=1; 
 unsigned int intKind = ALBANI; // default value for type of numerical integration
 complex<double> multFac = complex<double>(1.0,0.0);    // multiplicative factor for integration 
 unsigned int intType;  // type of integration to be done {INT_2D or INT_BOUND} 
 
  /* ---- check arguments ---- */
 if  (nrhs==0)
 {
  printHelp();
  return;
 }
 if (nlhs>1)
 {
  printf("Function call : [E] = window_integ(window,trans,rec,lamda,[sum,intType,tol])\n");
  mexErrMsgTxt("Too many ouput arguments!");
 }
 if ((nrhs!=4) && (nrhs!=5) && (nrhs!=6) && (nrhs!=7))
 {
  printf("Function call : [E] = window_integ(window,trans,rec,lamda,[sum,intType,tol])\n");
  mexErrMsgTxt("Wrong number of input arguments!");
 }
  /* check parameters for the number of elements */
 if ((mxIsComplex(prhs[0])) || 
     (mxGetM(prhs[0])!=4))
  mexErrMsgTxt("Something wrong with the window-vector! The number of elements per row is unequal to 4!");
 if ((mxIsComplex(prhs[1])) || 
     (mxGetM(prhs[1])*mxGetN(prhs[1])!=3))
  mexErrMsgTxt("Something wrong with the transmitter vector! The number of elements is unequal to 3!");
 if ((mxIsComplex(prhs[2])) || 
     ((mxGetM(prhs[2])*mxGetN(prhs[2])!=3) && ((mxGetM(prhs[2])==1) || (mxGetN(prhs[2])==1))) ||
     ((mxGetM(prhs[2])*mxGetN(prhs[2])!=3) && (mxGetM(prhs[2])!=3)) )
  mexErrMsgTxt("Something wrong with the receiver vector! The number of elements is unequal to 3!");
 if ((mxIsComplex(prhs[3])) || 
     (mxGetM(prhs[3])*mxGetN(prhs[3])==0))
  mexErrMsgTxt("Something wrong with the wave length vector! The number of elements should be larger than 0!");

 /* ---- process arguments ---- */
 if (nrhs>=5) 
  summing = (unsigned int) mxGetScalar(prhs[4]);    // read summing bit parameter
 if (nrhs>=6) 
  intKind = (unsigned int) mxGetScalar(prhs[5]);    // reading kind of integration
 if (nrhs==7) 
  tol = mxGetScalar(prhs[6]);    // read tolerance value to be set for integration
 par.trans = mxGetPr(prhs[1]);   // get transmitter position
 par.ndim = (unsigned int) (mxGetN(prhs[3])*mxGetM(prhs[3])); // get dimension of lams vector
 par.lams = mxGetPr(prhs[3]);   // get pointer to lams vector
 par.intKind = intKind;         // take over diffraction method
 par.int_save = Adapt(tol);
 switch (intKind)
 {
  case KIRCHHOFF: {par.f_inner = &f_inner_kirchhoff; 
                   intType=INT_2D; 
                   break;}
  case EXP: {par.f_inner = &f_inner_exp; 
             intType=INT_2D; 
             break;}
  case ADV_FRESNEL: {par.f_inner = &f_inner_adv_fresnel; 
                     intType=INT_2D; 
                     break;}
  case FRESNEL: {par.f_inner = &f_inner_fresnel; 
                 intType=INT_2D; 
                 break;}
  case RUBI: {par.f_inner = &f_inner_rubi; 
              intType=INT_BOUND; 
              break;}
  case ALBANI: {par.f_inner = &f_inner_albani; 
                intType=INT_BOUND; 
                break;}
  case MM_MAGIC: {intType=INT_NULL;  /* no integration needed */
                  break;}  
  case UP_FRESNEL: {intType=INT_NULL;  /* no integration needed */
                    break;}
  case FRESNEL_NO_INT: {intType=INT_NULL;  /* no integration needed */
                        break;}					
  default : mexErrMsgTxt("Parameter integration type is wrong! Value must be 1,2,3,4,5,6,7,8 or 9!"); break;
 }
 if ((intKind==FRESNEL) && (par.ndim>1))
  mexErrMsgTxt("The Fresnel integration type is currently not able to handle more than one frequency!\n The case is only for comparison with Matlab, so please use a Matlab implementation instead!");
 
 //printf("using method %i \n",intKind);
 
//  printf("Transmitter %f, %f, %f \n",par.trans[0],par.trans[1],par.trans[2]);
//  printf("Receiver %f, %f, %f \n",par.rec[0],par.rec[1],par.rec[2]);

 /* ---- main processing ---- */
 vector<complex<double> > ans(par.ndim);  // result of individual window
 vector<complex<double> > res(par.ndim,complex<double>(0.0,0.0));  // result of individual window
 Adapt s(tol);    // instance for integration
 
 unsigned int nr_rx_pos = (unsigned int) (mxGetM(prhs[2])*mxGetN(prhs[2])/3); // number of different receiver positions
 unsigned int nr_wins = (unsigned int) mxGetN(prhs[0]);  // number of windows
// printf("nr of rx pos %i \n",nr_rx_pos);

  /* ---- create output matrix ---- */
 if (summing==TRUE)
 {
  plhs[0] = mxCreateDoubleMatrix(nr_rx_pos,par.ndim,mxCOMPLEX); // output will be summed over all window contributions
//    // output for testing
//   for (unsigned int i=0; i<nr_rx_pos*par.ndim; i++)
//    mxGetPr(plhs[0])[i] = i;  
 }
else
{
    mwSize ndims = 3;
    mwSize dims[3] = {
        static_cast<mwSize>(nr_rx_pos),
        static_cast<mwSize>(par.ndim),
        static_cast<mwSize>(nr_wins)
    };

    plhs[0] = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxCOMPLEX);

    // No mxGetComplexDoubles here, because we're using legacy separate R/I below:
    // You already write with mxGetPr(plhs[0]) and mxGetPi(plhs[0]) later.
}
  
  /* ---- calculate output ---- */
 for (unsigned int r=0; r<nr_rx_pos; r++)   // go through receiver positions
 {
  par.rec = &(mxGetPr(prhs[2])[r*3]);     // get pointer to receiver position at r
  par.n = (par.rec[0]>0) ? 1.0:-1.0;      // direction toward receiver seen from the window
//  printf(" rec at [%f,%f,%f] \n",par.rec[0],par.rec[1],par.rec[2]);
  if (summing==TRUE)   // initialise matrix before calculating the window contribution
  {
   for (unsigned int j=0; j<par.ndim; j++)
   {
    mxGetPr(plhs[0])[r+j*nr_rx_pos] = 0;  mxGetPi(plhs[0])[r+j*nr_rx_pos] = 0;   // initialise values
   }
  }  
  for (unsigned int i=0; i<nr_wins; i++)  // over all windows
  { 
   par.ymin = mxGetPr(prhs[0])[i*4+0]; // ymin value for window i
   par.ymax = mxGetPr(prhs[0])[i*4+1]; // ymax value for window i
   par.zmin = mxGetPr(prhs[0])[i*4+2]; // zmin value for window i
   par.zmax = mxGetPr(prhs[0])[i*4+3]; // zmax value for window i
   if (intKind!=KIRCHHOFF)  
    multFac = getMultFac(intKind,&par); // calculate multiplicative factor
   if (intKind==FRESNEL)  // in case of FRESNEL, we need to change the integration ranges
   {
    double y01 = -par.trans[0]*(par.rec[1]-par.trans[1])/(par.rec[0]-par.trans[0])+par.trans[1];
    double z01 = -par.trans[0]*(par.rec[2]-par.trans[2])/(par.rec[0]-par.trans[0])+par.trans[2];
    double d1 = sqrt(par.trans[0]*par.trans[0]+(par.trans[1]-y01)*(par.trans[1]-y01)+
                     (par.trans[2]-z01)*(par.trans[2]-z01));
    double d2 = sqrt(par.rec[0]*par.rec[0]+(par.rec[1]-y01)*(par.rec[1]-y01)+
                     (par.rec[2]-z01)*(par.rec[2]-z01));
    double R1 = sqrt(par.lams[0]*d1*d2/(d1+d2));  // WORKS ONLY FOR A SINGLE LAMBDA !!!!
    par.ymin = sqrt(2.0)*(par.ymin-y01)/R1;  // lower value for u, u1
    par.ymax = sqrt(2.0)*(par.ymax-y01)/R1;  // upper value for u, u2
    par.zmin = sqrt(2.0)*(par.zmin-z01)/R1;  // lower value for v, v1
    par.zmax = sqrt(2.0)*(par.zmax-z01)/R1;  // upper value for v, v2
   }
   //printf(" %i win = [%f,%f,%f,%f] \n",i,par.ymin,par.ymax,par.zmin,par.zmax);
   switch (intType)
   {
    case INT_2D: {// standard 2d integration, integrate over y
                  s.integrate(f_outer,ans,par.ymin,par.ymax,&par);  // perform vector integration, output in variable ans
                  break;}
    case INT_BOUND: {f_boundary(ans,s,&par);  // use the boundary integration method
                     break;}
    case INT_NULL: { // no integral methods
                    switch (intKind)
                    {
                     case MM_MAGIC: {calcMM_Magic(ans, par);   // ans will be an array of complex double with dimensions par.ndim equal to the number of lambdas, i.e frequencies
                                     break;}
                     case UP_FRESNEL: {calcUpFresnel(ans,par); // ans will be an array of complex double with dimensions par.ndim equal to the number of lambdas, i.e frequencies
                                       break;}
                     case FRESNEL_NO_INT: {calcFresnelNoInt(ans,par); // ans will be an array of complex double with dimensions par.ndim equal to the number of lambdas, i.e frequencies
                                           break;}
                     default: mexErrMsgTxt("Something wrong, cannot identify the method to be used!\n");
                    }
                    break;}
    default: mexErrMsgTxt("Something wrong, cannot identify the integral method to be used!\n");
   }
   if (s.out_of_tolerance)
    mexErrMsgTxt("Required tolerance not reached from outer integration!\n");
   for (unsigned int j=0; j<par.ndim; j++)   // saving for each lambda given
   {
    if (summing==TRUE)   // summing over all window contributions
    {   
     mxGetPr(plhs[0])[r+j*nr_rx_pos] += real(multFac*ans[j]);
     mxGetPi(plhs[0])[r+j*nr_rx_pos] += imag(multFac*ans[j]);
    }
    else                 // saving each window contribution separately
    {
     mxGetPr(plhs[0])[r+j*nr_rx_pos+i*nr_rx_pos*par.ndim] = real(multFac*ans[j]);
     mxGetPi(plhs[0])[r+j*nr_rx_pos+i*nr_rx_pos*par.ndim] = imag(multFac*ans[j]);
    }
   }
    //res[j] = res[j]+ans[j];  // add results for the window i to overall result
  }
 }
//  for (unsigned int i=0; i<res.size(); i++)
//   printf("%i, result is %.10f + %.10f *j \n",i,real(res[i]),imag(res[i]));   
 
 
//  /* ---- create output ---- */
//  plhs[0] = mxCreateDoubleMatrix(1,par.ndim,mxCOMPLEX); // create vector
//  for (unsigned int i=0; i<res.size(); i++)
//  {
//   mxGetPr(plhs[0])[i] = real(res[i]);  // copy values
//   mxGetPi(plhs[0])[i] = imag(res[i]);
//  }
}

