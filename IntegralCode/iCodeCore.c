/*
 * =====================================================================================
 *
 *       Filename:  gkc-iCode.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *
 *         Author:  Paul P. Hilscher (2012) 
 *
 * =====================================================================================
 */

#include "iCodeCore.h"
#include "SpecialMath.h"

#include <complex.h>
#include <math.h>

#include<omp.h>

#include <fenv.h>

#define M_PI 3.141592653589793238462643383279502884197169


inline  Complex PlasmaDispersion(Complex z) {
   
        static const double  SQRT_PI = 1.7724538509055160272981674833411451;
        Complex  Z;
        Complex z2 = z*z;
       

        // sometimes imaginary part is too high, due to test functions, limit
        if(cimag(z) > 25.) z = creal(z) +  25.I;

        // Use asymptotic expansions for large values |z| >> 1
        // See definition in NRL Plasma Formulary
        if(creal(z*z) > 100) {
            
            const double x = creal(z); 
            const double y = cimag(z);
            
            const double sigma = (y >  0.) ? 0 : (( y == 0) ? 1. : 2.);
            Complex z4 = z2*z2;
            Complex z6 = z4*z2;
            Z = 1.0I *   SQRT_PI * sigma * cexp(-z2) - 1./z * ( 1. + 1./(2.*z2) + 3./(4.*z4) + 15./(8.*z6));		
  
        }
        else Z = 1.0I * SQRT_PI * cexp(- z2) * cerfc(- 1.0I * z);
    return Z;
}; 

///////////////////// Based on Z. Gao, PoP, 2005, swETG in Toroidal Plasma ////////////////


Complex getH(Complex w, double k_t, double k, double k_, double *X,  double Ls,  double Ln, const int Nx,
              const double q, const double mass, const double T, const double eta) 
{

/*
   // get ballooning angle
  const double theta  = k  / (shear * k_t); 
  const double theta_ = k_ / (shear * k_t); 

   const double g_tt=  (shear +1.) * (sin(theta) - sin(theta_)) 
                        - shear * (theta * cos(theta) - theta_ * cos(theta_)) 
                        - alpha/2. * (( theta - theta_ ) - sin(theta - theta_));
 
   const double a_j = 1. - ii * 2. * eps_n * tau * g_tt / (theta - theta_) * tau_j * Z;
   const double lambda_j = (1. + a_j)/2.;

   kp2  = pow2(k_t) * (1. + pow2(shear * theta  - alpha * sin(theta )))
   kp2_ = pow2(k_t) * (1. + pow2(shear * theta_ - alpha * sin(theta_)))

   const double b_g = (kp  *  kp_)/2. * (M / (tau * pow2(Z)));
   const double b_a = (kp2 + kp2_)/2. * (M / (tau * pow2(Z)));

   // integrate over tau
   H += 
   
   R =  exp( - pow2(kx - kx_)/pow2(2.*tau*L_s) * tau_j * M_j * a_j) 
     * (w * tau_j * Z - 1. + 3./2. * eta - eta / pow2(lambda_j) * 
         * (lambda_j - b_a + b_g  * I_1 / I_0) - eta * pow2(kx - kx_) / pow(2.*tau*L_s) * tau * M_j ) * G0
   
   ////// H_00
   //0 / infty
   H += exp(-ii * w * tau) / tau * R_ ;


   return ii * sqrt(tau_i * M_j) * sign(q_j)  / (2. *L_s * lambda_j) * H;
   */
};

///////////////////// Based on Z. Gao, PoP, 2002, swITG in Toroidal Plasma ////////////////

Complex getL(Complex w, double ky, double kx, double kx_, double *X,  double Ls,  double Ln, const int Nx,
            const double q, const double mass, const double T, const double eta) 
{

  Complex R= 0.;

  const double t = 0.;

         // We can ignore  terms with kx - kx_ >> 1, as integration over x is 0 due to oscillating term std::exp( ii * (kx - kx_) * X[x]);
         // ignore modes to far away as physically no connection but numerical problems
         // if(std::abs(kx-kx_) > 5.) continue;

         const double kp2    = ky*ky + kx * kx ;
         const double kp2_   = ky*ky + kx_* kx_;
     
         const double Omega  =  - q  / mass;

         const double v_th   = sqrt(2. * T / mass);
         //const double v_th   = sqrt(T / mass);
         const double rho    = v_th / Omega ;

         const double w_star = ky * T / ( Omega * mass * Ln );
         const double w_D    = 0. ;//#- w_star * Ln / LB;
     
         const double b  = kp2  * rho*rho / 2.;
         const double b_ = kp2_ * rho*rho / 2.;
         
         // Note : b_a >= b_g for every b,b_
         const double b_a = (b + b_) / 2.; // arithmetic mean
         const double b_g = sqrt(b * b_);  // geometric  mean
    
         // need asymptoctic approximation otherwise it fails, (points corresponds where
         // errors from double precision are large than due to expansion.
         const double G0 = (b_g < 10.) ? i0(b_g) * exp(-b_a) 
                                       : exp(b_g - b_a)/sqrt(2.*M_PI*b_g) * ( 1. + 1./(8.*b_g) +  9./(128.*b_g*b_g));
         const double G1 = (b_g < 10.) ? i1(b_g) * exp(-b_a)
                                       : exp(b_g - b_a)/sqrt(2.*M_PI*b_g) * ( 1. - 3./(8.*b_g) - 15./(128.*b_g*b_g));


         //use Sinh-Tanh Rule
         const double x_i[Nx], w_i[Nx];

  
         // Integrate over x, using fixed kx, kx_
   for(int x = 0; x < Nx; x++) {
         if(X[x] == 0.) continue;
         const double k_p    = (X[x]/Ls) * ky;


         const Complex  zeta   = (w - w_D * t) / (v_th * fabs(k_p));
         const Complex  Z      = PlasmaDispersion(zeta) ;
     
              
         const Complex z = w_star/w; 
         const Complex zeta2 = zeta*zeta;

         const Complex    Lq  =  (1.   - z) * zeta * Z * G0  
                                - eta * z * (zeta2 + (zeta2 -3./2.) * zeta * Z) * G0 
                                - eta * z * zeta * Z * ((1. - b_a) * G0 + b_g * G1);

         // Note : int_0^2\pi exp(1.j * x) = 0.
         R += Lq *  cexp( 1.0I * (kx_ - kx) * X[x]);
  } 
  
  return R;

};


void setupMatrix(Complex w, double ky, double *X, double *K,  double Ls,  double Ln, int Nx, Complex *M, double dh,
                double q, double mass, double T, double eta, double n) 
{
       
      #pragma omp parallel for collapse(2)
      for(int x_k = 0; x_k < Nx; x_k++) {  for(int w_k = 0; w_k < Nx; w_k++) {
         
        const int idx = x_k*Nx + w_k;              
        
          M[idx]  +=  - n * (q*q)  / T * ( (w_k == x_k ? 1. : 0. ) + 1./(2.*M_PI) *  getL(w, ky, K[x_k], K[w_k], X, Ls, Ln, Nx, q, mass, T, eta) * dh );

       } }

      return;
}
