/*
 * =====================================================================================
 *
 *       Filename:  gkc-iCode.h
 *
typedef double _Complex Complex; 
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/16/2012 10:23:13 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef GKC_ICODE
#define GKC_ICODE


#include <complex.h>

typedef double _Complex Complex; 

void setupMatrix(Complex w, double ky, double* X, double* K,  double Ls,  double Ln, int Nx, 
                 Complex* A,  double h, double q, double mass, double T, double eta, double rho); 


#endif // GKC_ICODE
