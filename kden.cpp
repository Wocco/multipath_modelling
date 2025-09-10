/* Function to estimation the pdf using a kernel estimator.
 * 
 *    [f] = kden(xi,x,[h],[lb,ub],[type])
 *    [f] = kden(xi,x,[h],[type])
 *    [f] = kden(xi,x,[lb,ub],[type])
 *
 * Function provides a kernel estimator using the data points x
 * at the x-values xi with smoothing factor h.
 *
 * Please note that in cases where the CDF shall be estimated
 * and a boundary is used, a numerical integration according to [3]
 * is used.
 * 
 * Input : xi = points to be evaluated.
 *         x = data to calculate the pdf from.
 *         h = smoothing factor (optional)
 *         [lb,ub] = lower and upper bound for bounded pdfs.
 *         type = either 'pdf' or 'cdf' for pdf of cdf calculation.
 *                type='pdf' is default.
 * Ouput : fi = density estimate at points xi.
 *
 * [1] B. W. Silverman, "Density Estimation for Statistics and Data Analysis",
 *     Chapman and Hall LtD, 1986
 * [2] M. C. Jones, "Simple boundary correction for kernel density estimation",
 *     Statistics and Computing, vol. 3, 1993, pp. 135-146
 * [3] W. Gander, W. Gautschi, "Adaptive Quadrature - Revisited",
 *     BIT Numerical Mathematics, 2000, Volume 40, Issue 1, pp. 84-101 
 *     http://www.inf.ethz.ch/personal/gander/
 * 
 *
 * For testing use:
    % test using the same PDF kernel estimation method, no boundaries
   x = 0:0.1:50;
   data = randn(1,1000)*2+30;
   fx = ones(length(data),1)*x-data(:)*ones(1,length(x));
   h = 1.06*std(data)*length(data)^(-1/5);
   y = sum(exp(-0.5*(fx./h).^2))/(length(data)*h)/sqrt(2*pi);
   fi = kden(x,data);
   plot(x,y,'.-',x,fi,'.r'); legend('Matlab','C'); grid on;
   fprintf('Difference between Matlab and C is in sum over all points : %g\n',sum(abs(y-fi)));
 *
 * or
    % test using the boundary PDF kernel estimation method
   x = 0:0.01:3;
   data = rand(1,10000)*2+0.5;
   y = zeros(1,length(x)); y((x>0.5)&(x<=2.5))=1/2;   % true pdf
   fi = kden(x,data,[0.5,2.5]);
   plot(x,y,'.-',x,fi,'.r'); legend('Matlab','C'); grid on;
 *
 * or
    % test using the CDF kernel estimation method
   x = 0:0.1:50;
   data = randn(1,1000)*2+30;
   fi = kden(x,data,'cdf');  % estimate the CDF
   plot(x,0.5*(1+erf((x-30)./sqrt(2*4))),'.-',x,fi,'.r'); legend('Matlab','C'); grid on;
 *
 * or
    % test using the CDF kernel estimation method
   x = 0:0.01:3;
   data = rand(1,10000)*2+0.5;
   y = 0.5*x-0.25;  y(x<0.5)=0; y(x>=2.5)=1;   % true cdf
   fi = kden(x,data,[0.5,2.5],'cdf');
   plot(x,y,'.-',x,fi,'.r'); legend('Matlab','C'); grid on;
 *
 * or (stats toolbox needed)
    % test using a Weibull distribution
   x = 0:0.01:3;
   data = wblrnd(1,1.2,1,10000);  % scale = 1, form = 1.2
   fi = kden(x,data,[0,inf],'cdf');
   plot(x,wblcdf(x,1,1.2),'.-',x,fi,'.r'); legend('Matlab','C'); grid on;
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
#include <math.h>
#include <string.h>
#include <mex.h> 


#define VERSION "3.0"
/* Version 1.0 was the initial version included the kernel estimator incl. the boundary conditions.
 * Version 2.0 includes the CDF estimation additionally.
 * Version 3.0 took out source code from Numerical Recipies in C++
 */

#define PI 3.14159265358979323846264338327950288419716939937510
#define DOUBLE_EPS 2.22044604925031e-16   /* taking eps from Matlab */
#define TOLERANCE 1e-4
#define TRUE 1
#define FALSE 0
#define INF 1e20 
#define parCDF "cdf"
#define parPDF "pdf"

#define SIMPSON 0
#define LOBATTO 1
#define USE_INT LOBATTO  // integration rule to be used for cdf


struct Param;  // forward declaration for struct Param, used in the integration method
void printParam(void *param);  // forward declaration for testing
//double erf(double x);  /* forward declaration */


using namespace std;

/* --------------- Class Adapt --------------- */
// Implementation of class Adapt from Adap.h (Numerical Recipies in C++)
class Adapt 
{
 public:
	 double TOL,toler;
	 bool terminate,out_of_tolerance;
  Adapt();
	 Adapt(double tol);
	 double integrate_lobatto(double func(double,void*),  
                           const double a, const double b, void* param);
  double integrate_simpson(double func(double,void*),  
                           const double a, const double b, void* param);
 private:
  void* par;   // parameter struct for function evaluation
	 double adaptlob(double func(double,void*), 
                  const double a, const double b, double fa, 
                  double fb, const double is);
  double adaptsim(double func(double,void*), 
                  const double a, const double b, double fa, double fm, 
                  double fb, const double is);
};



// Empty contructor, seting tolerance to standard value of 1e-6
Adapt::Adapt() : TOL(1e-6),terminate(true),out_of_tolerance(false),par(NULL)
{
	if (TOL < 10.0*DOUBLE_EPS)
		TOL=10.0*DOUBLE_EPS;
}

// Constructor with tolerance as parameter
Adapt::Adapt(double tol) : TOL(tol),terminate(true),out_of_tolerance(false),par(NULL)
{
	if (TOL < 10.0*DOUBLE_EPS)
		TOL=10.0*DOUBLE_EPS;
}

/* Function to perform the integration.
 * Function uses the adaptive Lobatto rule (see Gandner, Computermathematik, 1992).
 * Input : func = function pointer to the expression to be integrated.
 *         a,b = real valued integration boundaries.
 *         par = pointer to additional parameter for func.
 */
double Adapt::integrate_lobatto(double func(double,void*), const double a, const double b, void* param)
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
 double y[13];  // variables for function evaluations
 double m = (b+a)/2.0, h = (b-a)/2.0, wh = (b-a);  // pre-factors etc.
 double i4GL=0.0, i7Kron=0.0, i13Kron=0.0,is;         // init integral values
 double erri1,erri2,r;
 par = param;     // save parameter struct for later use in adaptlob
	for (unsigned int i=0; i<13; i++)
		y[i] = func(m+h*x[i],par);            // fill up the array with vectors from function evaluation
 for (unsigned int i=0; i<4; i++)
  i4GL += wh*w4GL[i]*y[x4GL[i]];        // get integral for 4-point Gauss-Lobatto rule
 for (unsigned int i=0; i<7; i++)
  i7Kron += wh*w7Kron[i]*y[x7Kron[i]];  // get integral for 7-point Kronrod rule
 for (unsigned int i=0; i<13; i++)
  i13Kron += wh*w13Kron[i]*y[i];        // get integral for 13-point Kronrod rule
 erri1 = fabs(i4GL-i13Kron);    // calculate as the norm(x,inf)=sum(abs(x))
 erri2 = fabs(i7Kron-i13Kron);  // calculate as the norm(x,inf)=sum(abs(x))
 r = (erri2!=0.0) ? (erri1/erri2) : (1.0);
 toler = ((r>0.0) && (r<1.0)) ? (TOL/r) : (TOL); // create tolerance value
 //printf("got err=%g, toler=%g \n",erri2,toler);
 if (erri2 <= toler)
  return i13Kron;    // seems to be good enough 
 is = i13Kron*toler/DOUBLE_EPS;   // first estimate for the integral
 if (i13Kron==0.0)
  is = wh;
 return adaptlob(func,a,b,y[0],y[12],is);   // perform integration
} 


/* Internal function to perform the integration.
 * Function uses the adaptive Lobatto rule (see Gandner, Computermathematik, 1992).
 *  Q = Adapt::adaptlob(func,a,b,fa,fb,is) tries to
 *   approximate the integral of func(double,void*) from a to b to
 *   an appropriate relative error. The argument 'func' is
 *   a pointer containing the pointer to the function func.  The remaining
 *   arguments are generated by Adapt::integrate_lobatto or by recursion.
 * Input : func = function pointer to the expression to be integrated.
 *         a,b = real valued integration boundaries.
 *         fa,fb = evalueated func at points a and b.
 *         is = first error estimate from integrate() for the integral.
 */
double Adapt::adaptlob(double func(double,void*),
                     const double a, const double b, double fa, 
                     double fb, const double is)
{
 const double x7Kron[7] = {-1.0,-sqrt(2.0/3.0),-1.0/sqrt(5.0),
                            0.0,1.0/sqrt(5.0),sqrt(2.0/3.0),1.0};  // x-values for 7-point Kronrod rule with x in [-1,1]
 const unsigned int x4GL[4] = {0,2,4,6};   // points used by 4-point Gauss-Lobatto rule
 const double w4GL[4] = {1.0/12.0, 5.0/12.0, 5.0/12.0, 1.0/12.0}; // weights for 4-point Gauss-Lobatto rule
 const double w7Kron[7] = {77.0/2940.0, 432.0/2940.0, 625.0/2940.0, 672.0/2940.0, 625.0/2940.0,
                           432.0/2940.0, 77.0/2940.0};   // weights for 7-point Kronrod rule
 double y[7], xval[7];   
 double m = (b+a)/2.0, h = (b-a)/2.0, wh = (b-a);  // pre-factors etc.
 double i4GL=0.0, i7Kron=0.0;                      // init integral values
 double err_sum = 0.0;
 if ((fabs(h)<=DOUBLE_EPS) || (m==a) || (m==b))
  return h*(fa+fb);   // we got a singularity most propbably, so stop here
 y[0] = fa;  y[6] = fb;   // points already evaluated                           
 xval[0] = a;  xval[7] = b;
 for (unsigned int i=1; i<6; i++)
 {
  xval[i] = m+h*x7Kron[i];    // x-value within [a,b]
		y[i] = func(xval[i],par);   // fill up the array with vectors from function evaluation
 }
 for (unsigned int i=0; i<4; i++)
  i4GL += wh*w4GL[i]*y[x4GL[i]];   // get integral for 4-point Gauss-Lobatto rule
 for (unsigned int i=0; i<7; i++)
  i7Kron += wh*w7Kron[i]*y[i];     // get integral for 7-point Kronrod rule
 err_sum = fabs(i4GL-i7Kron);   // calculate as the norm(x,inf)/dim=mean(abs(x)) 
 if ((is+err_sum <= is+toler) | (m<=a) | (b<=m))   // check if error is small enough, i.e. less than the tolerance wanted
  return i7Kron;    // finished
 else
 {
   // recursive call for all sub-intervals
  return adaptlob(func,a,xval[1],fa,y[1],is)+
		       adaptlob(func,xval[1],xval[2],y[1],y[2],is)+
		       adaptlob(func,xval[2],xval[3],y[2],y[3],is)+
		       adaptlob(func,xval[3],xval[4],y[3],y[4],is)+
		       adaptlob(func,xval[4],xval[5],y[4],y[5],is)+
		       adaptlob(func,xval[5],b,y[5],fb,is);  
 }
}

/* Function to perform the integration.
 * Function uses the adaptive Simpson rule (see Gandner, Computermathematik, 1992).
 * Input : func = function pointer to the expression to be integrated.
 *         a,b = real valued integration boundaries.
 *         par = pointer to additional parameter for func.
 */
double Adapt::integrate_simpson(double func(double,void*), const double a, const double b, void* param)
{
 const double invervSim[5]={0.9501,0.2311,0.6068,0.4860,0.8913};
 double y[3],sumyy=0.0,is; 
 unsigned int i;
 par = param;   // save parameter struct for later use in adaptsim
 y[0] = func(a,par);          // function at a
 y[1] = func((a+b)/2.0,par);  // function at (a+b)/2
 y[2] = func(b,par);          // function at b
 for (i=0; i<5; i++)
  sumyy = sumyy+func(a+invervSim[i]*(b-a),par);
 is = (b-a)/8.0*(y[0]+y[1]+y[2]+sumyy);
 if (is==0.0)
  is = b-a;
 is = is*TOL/DOUBLE_EPS;
	return adaptsim(func,a,b,y[0],y[1],y[2],is);   // perform integration
}

/* Internal function to perform the integration.
 * Function uses the adaptive Simpson rule (see Gandner, Computermathematik, 1992).
 *  Q = Adapt::adaptlob(func,a,b,fa,fb,is) tries to
 *   approximate the integral of func(double,void*) from a to b to
 *   an appropriate relative error. The argument 'func' is
 *   a pointer containing the pointer to the function func.  The remaining
 *   arguments are generated by Adapt::integrate_simpson or by recursion.
 * Input : func = function pointer to the expression to be integrated.
 *         a,b = real valued integration boundaries.
 *         fa,fb,fm = evalueated func at points a, b and (a+b)/2.
 *         is = first error estimate from integrate() for the integral.
 */
double Adapt::adaptsim(double func(double,void*),
                       const double a, const double b, double fa, double fm,
                       double fb, const double is)
{
 double m = (a+b)/2.0;
 double fml,fmr,i1,i2;
  // evaluate two more functions in the area [a,b]
 fml = func((b+3.0*a)/4.0,par);  // function evaluations
 fmr = func((3.0*b+a)/4.0,par);  // function evaluations
 i1 = (b-a)*(fa + 4.0*fm + fb)/6.0;   // Three point Simpson's rule.
 i2 = (b-a)/12.0 * (fa + 4.0 * (fml+fmr) + 2.0*fm + fb); // Five point double Simpson's rule.
 i1 = (16.0*i2 - i1)/15.0;   // One step of Romberg extrapolation.
 if (fabs(i1-i2)<=TOL)  // needed tolerance reached?
  return i1;  // output
 if (((b-a)<=DOUBLE_EPS) || (a==m) || (b==m))
 {
   // Error, Interval contains no more machine number,
   //        required tolerance may not be met.
  out_of_tolerance=true;   
	 terminate=false;
 }
  // perform another loop ...
 return adaptsim(func,a,m,fa,fml,fm,is)+
        adaptsim(func,m,b,fm,fmr,fb,is);  // output otherwise
}


/* -------- Kernel functions --------- */
/* -- Gaussian Kernel -- */

double gaussPdf(double x)
{
 return exp(-0.5*x*x)/sqrt(2.0*PI);  /* Gaussian kernel, [1], p.43 */
}

double gaussCdf(double x)
{
 return 0.5*(1+erf(x/sqrt(2.0)));  /* using the Gaussian kernel */
}

/* a0(lb,ub) for boundary option using linear combination */
double gauss_a0(double lb, double ub)
{
 return 0.5*erf(ub/sqrt(2.0))-0.5*erf(lb/sqrt(2.0));
}

/* a1(lb,ub) for boundary option using linear combination */
double gauss_a1(double lb, double ub)
{
 return -exp(-0.5*ub*ub)/sqrt(2*PI)+exp(-0.5*lb*lb)/sqrt(2.0*PI);
}

/* a2(lb,ub) for boundary option using linear combination */
double gauss_a2(double lb, double ub)
{
 return -ub*exp(-0.5*ub*ub)/sqrt(2*PI)+0.5*erf(ub/sqrt(2.0))+lb*exp(-0.5*lb*lb)/sqrt(2*PI)-0.5*erf(lb/sqrt(2.0));
}


/* -------- Helper Subroutines --------- */

/* Function to calculate the standard deviation of a vector.
 * Input : x = array of double values
 *         n = length of array x
 */
double stdev(double *x, unsigned int n)
{
 double sum=0,mean,sq_diff_sum=0;
 unsigned int i;
 if (n<=1)    /* using the unbiased formula, so no calculation for n=1 */
  return 0.0;
 for (i=0; i<n; i++) /* first loop to calculate the mean */
  sum += x[i];     
 mean = sum/n; /* calculate the mean */
 /* printf("mean = %g\n",mean); */
 sq_diff_sum=0;
 for (i=0; i<n; i++)   /* second loop for the square */
  sq_diff_sum += (x[i]-mean)*(x[i]-mean);  /* calculate sum((xi-mx)^2) */
 return sqrt(sq_diff_sum/(n-1));   /* return standard deviation, unbiased version */
/* return sqrt(sq_diff_sum/n);   // return standard deviation, biased version */
}

/* Function to calculate the PDF kernel estimation without any boundaries.
 * Input : x = x-value of the PDF to be estimated.
 *         data = data values provided.
 *         nData = length of data vector.
 *         h = smoothing factor of the estimator.
 *         kernel = kernel function to be used.
 * Ouput : Estimated PDF value.  
 */
double kernelEstPDF(double x, double *data, unsigned int nData, 
                    double h, double (*kernel)(double))
{
 unsigned int i;
 double fi=0.0; /* set up to zero for startup */
 for (i=0; i<nData; i++)  /* go through each data point for summation */
 {
  fi += kernel((x-data[i])/h);  /* perform kernel estimation */
 }
 fi /= nData*h;  /* normalise to n*h, [1], p.15, eq. (2.2a) */
 return fi;
}

/* Function to calculate the CDF kernel estimation without any boundaries.
 * Input : x = x-value of the CDF to be estimated.
 *         data = data values provided.
 *         nData = length of data vector.
 *         h = smoothing factor of the estimator.
 *         kernel = kernel function to be used.
 * Ouput : Estimated CDF value.  
 */
double kernelEstCDF(double x, double *data, unsigned int nData, 
                    double h, double (*kernel)(double))
{
 unsigned int i;
 double fi=0.0; /* set up to zero for startup */
 for (i=0; i<nData; i++)  /* go through each data point for summation */
 {
  fi += kernel((x-data[i])/h);  /* perform kernel estimation */
 }
 fi /= nData;  /* normalise to n */
 return fi;
}


// Parameter struct for bounded estimation 
struct Param
{
 double *data;                 // Pointer to sample data
 unsigned int nData;           // Length of the array data
 double h;                     // Smoothing factor
 double lb;                    // Lower bound
 double ub;                    // Upper bound
 double (*kernel)(double);     // Kernel function
 double (*a0)(double,double);  // Boundary functions
 double (*a1)(double,double); 
 double (*a2)(double,double);
};

/* for testing */
void printParam(void *param)
{
 Param *par = (Param*) param;  /* cast from void to type Param */
 printf(" data = %x \n",par->data);
 printf(" nData = %i \n",par->nData);
 printf(" h = %g \n",par->h);
 printf(" lb = %g \n",par->lb);
 printf(" ub = %g \n",par->ub);
 printf(" kernel = %x \n",par->kernel);
 printf(" a0 = %x \n",par->a0);
 printf(" a1 = %x \n",par->a1);
 printf(" a2 = %x \n",par->a2); 
}

/* Function to calculate the PDF kernel estimation with boundaries.
 * Input : x = x-value of the PDF to be estimated.
 *         param = parameter struct 
 * Ouput : Estimated PDF value.  
 */
double kernelEstBoundPDF(double x, void *param)
{
 unsigned int i;
 double z,lv,uv;        /* dummy variables */
 double fi=0.0; /* set up to zero for startup */
 Param *par = (Param*) param;  /* cast from void to type Param */
 if ((x>=par->lb) && (x<=par->ub))  /* perform only, when within boundary */
 {
  for (i=0; i<par->nData; i++)  /* go through each data point for summation */
  {
   z = (x-par->data[i])/par->h;
   lv = (par->lb-x)/par->h;    uv = ((par->ub)-x)/par->h;
   fi += (par->kernel(z))*(par->a2(lv,uv)-par->a1(-uv,-lv)*z)/(par->a2(lv,uv)*par->a0(lv,uv)-par->a1(-uv,-lv)*par->a1(-uv,-lv));  /* perform kernel estimation */
   //printf(" at %g, z=%g, lv=%g, uv=%g fi=%g \n",x,z,lv,uv,fi);
  }
  fi /= (par->nData)*(par->h); /* normalise to n*h, [1], p.15, eq. (2.2a) */
  if (fi<0)   
   fi = 0.0; /* just to secure, should not happen */
 }
 return fi;
}  


/* -------- Main Matlab wrapper --------- */

/* Function to output the help of the routine */
void printHelp()
{
 printf(" Function to calculate the kernel density estimate.\n");
 printf("\n");
 printf("    [f] = kden(xi,x,[h],[lb,ub],[type])\n");
 printf("    [f] = kden(xi,x,[h],[type])\n");
 printf("    [f] = kden(xi,x,[lb,ub],[type])\n");
 printf("\n");
 printf(" Function provides a kernel estimator using the data points x\n");
 printf(" at the x-values xi with smoothing factor h.\n");
 printf("\n"); 
 printf(" Please note that in cases where the CDF shall be estimated\n"); 
 printf(" and a boundary is used, a numerical integration according to [3]\n"); 
 printf(" is used.\n"); 
 printf("\n"); 
 printf(" Input : xi = points to be evaluated.\n");
 printf("         x = data to calculate the pdf from.\n");
 printf("         h = smoothing factor (optional)\n");
 printf("         [lb,ub] = lower and upper bound for bounded pdfs.\n");
 printf("         type = either 'pdf' or 'cdf' for pdf of cdf calculation.\n");
 printf("                type='pdf' is default.\n");
 printf(" Ouput : fi = density estimate at points xi.\n");
 printf("\n");
 printf(" [1] B. W. Silverman, \"Density Estimation for Statistics and Data Analysis\",\n");
 printf("     Chapman and Hall LtD, 1986\n");
 printf(" [2] M. C. Jones, \"Simple boundary correction for kernel density estimation\",\n");
 printf("     Statistics and Computing, vol. 3, 1993, pp. 135-146\n"); 
 printf(" [3] W. Gander, W. Gautschi, \"Adaptive Quadrature - Revisited\",\n");
 printf("     BIT Numerical Mathematics, 2000, Volume 40, Issue 1, pp. 84-101 \n");
 printf("     http://www.inf.ethz.ch/personal/gander/\n");
 printf("\n");
 printf("  For testing use:\n");
 printf("    %% test using the same PDF kernel estimation method, no boundaries\n");
 printf("   x = 0:0.1:50;\n");
 printf("   data = randn(1,1000)*2+30;\n");
 printf("   fx = ones(length(data),1)*x-data(:)*ones(1,length(x));\n");
 printf("   h = 1.06*std(data)*length(data)^(-1/5);\n");
 printf("   y = sum(exp(-0.5*(fx./h).^2))/(length(data)*h)/sqrt(2*pi);\n");
 printf("    fi = kden(x,data);\n");
 printf("   plot(x,y,'.-',x,fi,'.r'); legend('Matlab','C'); grid on;\n");
 printf("   fprintf('Difference between Matlab and C is in sum over all points : ');\n");
 printf("   sum(abs(y-fi)) \n");
 printf("\n"); 
 printf(" or\n");
 printf("   %% test using the boundary PDF kernel estimation method\n");
 printf("  x = 0:0.01:3;\n");
 printf("  data = rand(1,10000)*2+0.5;\n");
 printf("  y = zeros(1,length(x)); y((x>0.5)&(x<=2.5))=1/2;   %% true pdf\n");
 printf("  fi = kden(x,data,[0.5,2.5]);\n");
 printf("  plot(x,y,'.-',x,fi,'.r'); legend('Matlab','C'); grid on;\n");
 printf("\n");
 printf(" or\n");
 printf("   %% test using the CDF kernel estimation method\n");
 printf("  x = 0:0.1:50;\n");
 printf("  data = randn(1,1000)*2+30;\n");
 printf("  fi = kden(x,data,'cdf');  %% estimate the CDF\n");
 printf("  plot(x,0.5*(1+erf((x-30)./sqrt(2*4))),'.-',x,fi,'.r'); legend('Matlab','C'); grid on;\n");
 printf("\n");
 printf(" or\n");
 printf("   %% test using the CDF kernel estimation method\n");
 printf("  x = 0:0.01:3;\n");
 printf("  data = rand(1,10000)*2+0.5;\n");
 printf("  y = 0.5*x-0.25;  y(x<0.5)=0; y(x>=2.5)=1;   %% true cdf\n");
 printf("  fi = kden(x,data,[0.5,2.5],'cdf');\n");
 printf("  plot(x,y,'.-',x,fi,'.r'); legend('Matlab','C'); grid on;\n");
 printf("\n");
 printf(" or (stats toolbox needed)\n");
 printf("   %% test using a Weibull distribution\n");
 printf("  x = 0:0.01:3;\n");
 printf("  data = wblrnd(1,1.2,1,10000);  %% scale = 1, form = 1.2\n");
 printf("  fi = kden(x,data,[0,inf],'cdf');\n");
 printf("  plot(x,wblcdf(x,1,1.2),'.-',x,fi,'.r'); legend('Matlab','C'); grid on;\n");
 printf("\n\n");
 printf(" This is version %s from %s of file %s.\n",VERSION,__DATE__,__FILE__);
 #if USE_INT==SIMPSON
  printf("  Here the SIMPSON rule is used.\n\n");
 #elif USE_INT==LOBATTO
  printf("  Here the LOBATTO rule is used.\n\n");
 #else
  #error "The parameter USE_INT is set wrongly, don't know the integration rule to be used!"
 #endif 
 printf(" Copyright 2016 Thomas Jost, German Aerospace Center (DLR)\n\n");
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


/* --- Main Matlab Routine --- */

/* Function to estimation the pdf using a kernel estimator.
 * 
 *    [f] = kden(xi,x,[h],[lb,ub],[type])
 *    [f] = kden(xi,x,[h],[type])
 *    [f] = kden(xi,x,[lb,ub],[type])
 *
 * Function provides a kernel estimator using the data points x
 * at the x-values xi with smoothing factor h.
 * 
 * Input : xi = points to be evaluated.
 *         x = data to calculate the pdf from.
 *         h = smoothing factor (optional)
 *         [lb,ub] = lower and upper bound for bounded pdfs.
 *         type = either 'pdf' or 'cdf' for pdf of cdf calculation.
 *                type='pdf' is default.
 * Ouput : fi = density estimate at points xi.
 *
 * [1] B. W. Silverman, "Density Estimation for Statistics and Data Analysis",
 *     Chapman and Hall LtD, 1986
 * [2] M. C. Jones, "Simple boundary correction for kernel density estimation",
 *     Statistics and Computing, vol. 3, 1993, pp. 135-146
 * 
 *
 * For testing use:
    % test using the same PDF kernel estimation method, no boundaries
   x = 0:0.1:50;
   data = randn(1,1000)*2+30;
   fx = ones(length(data),1)*x-data(:)*ones(1,length(x));
   h = 1.06*std(data)*length(data)^(-1/5);
   y = sum(exp(-0.5*(fx./h).^2))/(length(data)*h)/sqrt(2*pi);
   fi = kden(x,data);
   plot(x,y,'.-',x,fi,'.r'); legend('Matlab','C'); grid on;
   fprintf('Difference between Matlab and C is in sum over all points : %g\n',sum(abs(y-fi)));
 *
 * or
    % test using the boundary PDF kernel estimation method
   x = 0:0.01:3;
   data = rand(1,10000)*2+0.5;
   y = zeros(1,length(x)); y((x>0.5)&(x<=2.5))=1/2;   % true pdf
   fi = kden(x,data,[0.5,2.5]);
   plot(x,y,'.-',x,fi,'.r'); legend('Matlab','C'); grid on;
 *
 * or
    % test using the CDF kernel estimation method
   x = 0:0.1:50;
   data = randn(1,1000)*2+30;
   fi = kden(x,data,'cdf');  % estimate the CDF
   plot(x,0.5*(1+erf((x-30)./sqrt(2*4))),'.-',x,fi,'.r'); legend('Matlab','C'); grid on;
 *
 * or
    % test using the CDF kernel estimation method
   x = 0:0.01:3;
   data = rand(1,10000)*2+0.5;
   y = 0.5*x-0.25;  y(x<0.5)=0; y(x>=2.5)=1;   % true cdf
   fi = kden(x,data,[0.5,2.5],'cdf');
   plot(x,y,'.-',x,fi,'.r'); legend('Matlab','C'); grid on;
 *
 * or (stats toolbox needed)
    % test using a Weibull distribution
   x = 0:0.01:3;
   data = wblrnd(1,1.2,1,10000);  % scale = 1, form = 1.2
   fi = kden(x,data,[0,inf],'cdf');
   plot(x,wblcdf(x,1,1.2),'.-',x,fi,'.r'); legend('Matlab','C'); grid on;
 *
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

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
 unsigned int nData=0;  /* Number of data points, so length(x) */
 unsigned int n=0;      /* Number of point to be calculated, so length(xi) */
 double *xData, *xi, *fi;    /* Pointer input/output vectors */
 double (*kernel)(double);   /* Function pointer to the kernel */ 
 double (*a0)(double,double);  /* extra function for boundary conditions */
 double (*a1)(double,double);  
 double (*a2)(double,double);  
 unsigned int l,i;      /* run variables */
 double h=0;            /* smoothing factor */
 double lb=0,ub=0;      /* lower and upper bound */
 char parTypeVal[4];    /* dummy for saving the value of paramter type */ 
 unsigned int typePDF=TRUE; /* calculated type={'pdf','cdf'}, typePDF=TRUE->type='pdf' */
 unsigned int valPar_h=0, valPar_b=0, valPar_type=0;  /* parameter value for optional input */
 Adapt integ(TOLERANCE);      /* object instance for integration */
 Param par;             /* parameter struct for bounded estimation */
  
  /* ---- check arguments ---- */
 if  (nrhs==0)
 {
  printHelp();
  return;
 }
 if (nlhs>1)
 {
  printf("Function call : [fi] = kden(xi,x,[h],[lb,ub],[type])\n");
  printf("                [fi] = kden(xi,x,[lb,ub],[type])\n");
  printf("                [fi] = kden(xi,x,[type])\n");
  mexErrMsgTxt("Too many ouput arguments!");
 }
 if (nrhs>5)
 {
  printf("Function call : [fi] = kden(xi,x,[h],[lb,ub],[type])\n");
  printf("                [fi] = kden(xi,x,[lb,ub],[type])\n");
  printf("                [fi] = kden(xi,x,[type])\n");
  mexErrMsgTxt("Wrong number of input arguments!");
 }
  /* check parameters for the number of elements */
 if ((mxIsComplex(prhs[0])) || (mxIsComplex(prhs[1])) || 
     ((nrhs>=3) && (mxIsComplex(prhs[2]))) || ((nrhs==4) && (mxIsComplex(prhs[3]))))
  mexErrMsgTxt("Some values are complex! This is invalid.");

 /* ---- start processing ---- */
 nData = (unsigned int) (mxGetN(prhs[1])*mxGetM(prhs[1])); /* Number of data points */
 n = (unsigned int) (mxGetN(prhs[0])*mxGetM(prhs[0])); /* Number of points to be evaluated */
 xData = mxGetPr(prhs[1]);   /* get pointer to data vector */
 xi = mxGetPr(prhs[0]);      /* get pointer to xi values */
  /* now check for the other parameters */
 for (i=2; i<(unsigned int) nrhs; i++)  /* go through each optional parameter */
 {
  if ((valPar_h==0) && (mxGetN(prhs[i])*mxGetM(prhs[i])==1))  /* check for paramter to be h */
   valPar_h = i;
  else if ((valPar_b==0) && (mxGetN(prhs[i])*mxGetM(prhs[i])==2))  /* check for paramter to be bounds [lb,ub] */
   valPar_b = i;
  else if ((valPar_type==0) && (mxIsChar(prhs[i])))   /* check for type, i.e. {'pdf','cdf'} */
  {
   valPar_type = i;
   mxGetString(prhs[i],parTypeVal,4);  /* read out string */
  }
  else
  {
   printf("Parameter %i is unrecognized! Please check!\n",i+1);
   mexErrMsgTxt("Wrong parameter! \n");
  }
 }
  /* now read out parameters */
 if (valPar_h!=0)
  h = mxGetScalar(prhs[valPar_h]);   /* get smoothing factor */
 else
  h = 1.06*stdev(xData,nData)*pow(nData,-0.2);  /* calculate smoothing value h according to [1], p.45 */
 if (valPar_b!=0)                   
 {
   /* set boundaries, if boundaries are not set, then ub=lb=0 */
  lb = (mxGetPr(prhs[valPar_b])[0]<mxGetPr(prhs[valPar_b])[1]) ? (mxGetPr(prhs[valPar_b])[0]):(mxGetPr(prhs[valPar_b])[1]);  /* get lower bound */
  ub = (mxGetPr(prhs[valPar_b])[0]>=mxGetPr(prhs[valPar_b])[1]) ? (mxGetPr(prhs[valPar_b])[0]):(mxGetPr(prhs[valPar_b])[1]);  /* get upper bound */
   /* now check boundaries */
  if (mxIsInf(ub))  
   ub = INF;
  if (mxIsInf(lb))
   lb = -INF;
  if ((ub==INF) && (lb==-INF))
  {
    ub=lb=0;  /* basically no boundary, so erase */
  } 
 }
 if ((valPar_type!=0) && (strncmp(parCDF,parTypeVal,3)==0))  /* check if type shall be 'cdf' */
  typePDF = FALSE;   /* set type to 'cdf' */
  /* set kernel functions */
 if ((typePDF==TRUE) || ((ub!=0) || (lb!=0)))
 {
  kernel = &gaussPdf;
 }
 else
 {
  kernel = &gaussCdf;
 }
 /* set up additional functions for boundary condition */
 par.a0 = &gauss_a0;
 par.a1 = &gauss_a1;
 par.a2 = &gauss_a2;
 
  /* now start with output */
 plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]),mxGetN(prhs[0]),mxREAL); /* output vector as [1 x n] */
 fi = mxGetPr(plhs[0]);   /* put up the pointer for simplified access */
 if ((ub==0) && (lb==0))  /* so no boundary condition */
 {
  for (l=0; l<n; l++)    /* go through each point to be calculated */
  {
   if (typePDF==TRUE)  
    fi[l] = kernelEstPDF(xi[l], xData, nData, h, kernel); /* perform kernel estimation for PDF */
   else
    fi[l] = kernelEstCDF(xi[l], xData, nData, h, kernel); /* perform kernel estimation for CDF */
  }
 }
 else
 {
   /* use the boundary version */
   /* set up parameter struct */
  par.data = xData;
  par.nData = nData;
  par.h = h;
  par.lb = lb;
  par.ub = ub;
  par.kernel = kernel;
  for (l=0; l<n; l++)    /* go through each point to be calculated */
  {
   if (typePDF==TRUE)  
   {
     /* calculate the PDF using boundary correction */
    fi[l] = kernelEstBoundPDF(xi[l], &par); /* perform kernel estimation with boundaries */
    //printf("Got at %g %g \n\n",xi[l],fi[l]);
   }
   else
   {
     /* calculate the value for the CDF using boundary correction */
    if ((xi[l]>lb) && (xi[l]<ub))
     #if USE_INT==SIMPSON
      fi[l] = integ.integrate_simpson(&kernelEstBoundPDF,lb,xi[l],&par);   // perform integration over pdf
     #elif USE_INT==LOBATTO
      fi[l] = integ.integrate_lobatto(&kernelEstBoundPDF,lb,xi[l],&par);   // perform integration over pdf
     #else
      #error "The parameter USE_INT is set wrongly, don't know the integration rule to be used!"
     #endif      
    else if (xi[l]<=lb)
     fi[l] = 0.0;
    else if (xi[l]>=ub)
     fi[l] = 1.0;
     /* check for integration */
    if (integ.out_of_tolerance)
     mexErrMsgTxt("Required tolerance not reached in f() to calculate the CDF!\n");
   }
  }  
 } 
 
 
  /* for testing of erf(), please make sure to call kden(x,randn(1,10)) */
/*  for (i=0; i<n; i++)
     mxGetPr(plhs[0])[i] = erf(xi[i]); */ /* calc erf() and return output */
}

