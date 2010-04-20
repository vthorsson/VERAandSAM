#include <math.h>
#include "arraystats.h"

/* global external variables */
extern int N, M, *m;
extern float **X, **Y, *muX, *muY;
extern double *cov; 
extern int gene_index;
extern int N1, N2, *m1, *m2;
extern float **X1, **X2, **Y1, **Y2;
extern double *cov1, *cov2;
extern int gene_index_1, gene_index_2;

/* global definitions */
double Sex, Sey, Re, Sdx, Sdy, Rd; 
double Vex, Vey, Vdx, Vdy; 
double C11, C12, C22, detC ; 
double B11, B12, B22; 

#define PI 3.1415926


/********************************************************************************************/
/**                                                                                        **/
/**  Z.c - Objective function (variables beta)                                            **/
/**                                                                                        **/
/**  April 2000 Copyright (c) Vesteinn Thorsson and Trey Ideker                            **/
/********************************************************************************************/

/* X and Y must be shifted my mus */

double Z( double *beta ){

  int i,j, pisum;
  double logsum, quadsum;

  Sex = beta[1];
  Sey = beta[2]; 
  Re  = beta[3];
  Sdx = beta[4];
  Sdy = beta[5];
#ifdef USERHOD  
  Rd  = beta[6];
#else
  Rd = 0;
#endif
  
  Vex = Sex*Sex; 
  Vey = Sey*Sey;
  Vdx = Sdx*Sdx;
  Vdy = Sdy*Sdy;
  
  pisum = 0;
  logsum = 0.;
  quadsum = 0.;
  for (i=0 ; i<N ; i++ ){ /* loop over genes */

    /* construct C, the covariance matrix */    
    C11 = Vex* muX[i]*muX[i] + Vdx;
    C12 = Sex*Sey*Re* muX[i]*muY[i] + Sdx*Sdy*Rd;
    C22 = Vey* muY[i]*muY[i] + Vdy;
    detC = C11*C22 - C12*C12;

   logsum += m[i] * log(detC); /* one contribution from each sample */

    /* construct B, the inverse of C */
    B11 = C22 / detC;
    B12 = -C12 / detC;
    B22 = C11 / detC;     
    
    for (j=0 ; j<m[i]; j++ ){ /* loop over repeats */
      quadsum += B11 *  X[i][j] * X[i][j] \
	+ B22 * Y[i][j] * Y[i][j] \
	+ 2 * B12 * X[i][j] * Y[i][j];
    }

    pisum += m[i]; /* will give total number samples over all genes */

  }

  return ( log(2.*PI) * pisum + 0.5 * logsum + 0.5 * quadsum) ;
  
}


/********************************************************************************************/
/**                                                                                        **/
/**  gradZ.c - Gradient of objective function (variables beta)                            **/
/**                                                                                        **/
/**  April 2000 Copyright (c) Vesteinn Thorsson and Trey Ideker                            **/
/********************************************************************************************/

/* X and Y must be shifted my mus */


void gradZ( double *beta, double *grad ){


  int i,j;
  double temp, x, y, mu_x, mu_y;
  double temp1, temp2, temp3, temp4, temp5;
  double ddetC_dg1, ddetC_dg2, ddetC_dg3, ddetC_dg4, ddetC_dg5; 
#ifdef USERHOD
  double ddetC_dg6, temp6
#endif 

  Sex = beta[1];
  Sey = beta[2]; 
  Re  = beta[3];
  Sdx = beta[4];
  Sdy = beta[5];
#ifdef USERHOD  
  Rd  = beta[6];
#else
  Rd =0;
#endif
  
  Vex = Sex*Sex; 
  Vey = Sey*Sey;
  Vdx = Sdx*Sdx;
  Vdy = Sdy*Sdy;
  
  for (i=1; i<=NDIM; i++){
    grad[i]=0.;
  }

  for (i=0 ; i<N ; i++ ){ /* loop over genes */

    mu_x = muX[i]; 
    mu_y = muY[i];

    /* construct C, the covariance matrix of x and y */    
    C11 = Vex* mu_x*mu_x + Vdx;
    C12 = Sex*Sey*Re* mu_x*mu_y + Sdx*Sdy*Rd;
    C22 = Vey* mu_y*mu_y + Vdy;
    detC = C11*C22 - C12*C12;

    /* construct B, the inverse of C */
    B11 = C22 / detC;
    B12 = -C12 / detC;
    B22 = C11 / detC;     

    /* Compute gradient of detC
       for s = { g1, ..., g6} :
       d(detC)/ds = d(C11)/ds * C22 + d(C22)/ds * C11 - 2 C12 * d(C12)/ds */

    ddetC_dg1 = 2.*Sex* mu_x*mu_x * C22 - 2.*Sey*Re *C12 *mu_x*mu_y ;
    ddetC_dg2 = 2.*Sey* mu_y*mu_y * C11 - 2.*Sex*Re *C12 *mu_x*mu_y ;  
    ddetC_dg3 = -2.*Sex*Sey * C12 * mu_x*mu_y; 
    ddetC_dg4 = 2.*Sdx * C22 - 2*Sdy*Rd* C12 ;
    ddetC_dg5 = 2.*Sdy * C11 - 2*Sdx*Rd* C12 ;  
#ifdef USERHOD
    ddetC_dg6 = -2.*Sdx*Sdy * C12;
#endif

    /* Compute gradient of objective function 
       -1/(2 detC) * ( d(detC)/ds ( 1 - [x y]B[xy]) + x^2 d(C22)/ds +
       y^2 d(C11)/ds - 2 xy d(C12)/ds )  */

    temp1 = m[i] * ddetC_dg1;
    temp2 = m[i] * ddetC_dg2;
    temp3 = m[i] * ddetC_dg3;
    temp4 = m[i] * ddetC_dg4;
    temp5 = m[i] * ddetC_dg5;
#ifdef USERHOD
    temp6 = m[i] * ddetC_dg6;
#endif

    for (j=0 ; j<m[i]; j++ ){ /* loop over repeats */

      x=X[i][j];
      y=Y[i][j];

      temp = B11* x * x + B22 * y * y+ 2. * B12 * x * y;

      temp1 += -ddetC_dg1 * temp + 2.*Sex* mu_x*mu_x*y*y - 2.*Sey*Re* mu_x*mu_y*x*y ;
      temp2 += -ddetC_dg2 * temp + 2.*Sey* mu_y*mu_y*x*x - 2.*Sex*Re* mu_x*mu_y*x*y ;
      temp3 += -ddetC_dg3 * temp - 2.*Sex*Sey* mu_x*mu_y*x*y ;
      temp4 += -ddetC_dg4 * temp + 2.*Sdx* y*y - 2.*Sdy*Rd* x*y ;
      temp5 += -ddetC_dg5 * temp + 2.*Sdy* x*x - 2.*Sdx*Rd* x*y ; 
#ifdef USERHOD
      temp6 += -ddetC_dg6 * temp - 2.*Sdx*Sdy* x*y ; 
#endif

    }

    grad[1] += temp1/detC;
    grad[2] += temp2/detC;
    grad[3] += temp3/detC;
    grad[4] += temp4/detC;
    grad[5] += temp5/detC;
#ifdef USERHOD
    grad[6] += temp6/detC;
#endif

  }

  for(i=1;i<=NDIM;i++){
    grad[i] /= (double)2.;
  }

}

/********************************************************************************************/
/**                                                                                        **/
/**  gradZe.c - Gradient of objective function in beta[1,2,3] directions only             **/
/**             (epsilon contribution). (Returns gradient[4,5,6]=0. )                      **/
/**                                                                                        **/
/**  April 2000 Copyright (c) Vesteinn Thorsson and Trey Ideker                            **/
/********************************************************************************************/

/* X and Y must be shifted my mus */

void gradZe( double *beta, double *grad ){

  int i,j;
  double temp, temp1, temp2, temp3, x, y, mu_x, mu_y;
  double ddetC_dg1, ddetC_dg2, ddetC_dg3; 

  Sex = beta[1];
  Sey = beta[2]; 
  Re  = beta[3];
  Sdx = beta[4];
  Sdy = beta[5];
#ifdef USERHOD  
  Rd  = beta[6];
#else
  Rd =0;
#endif  

  Vex = Sex*Sex; 
  Vey = Sey*Sey;
  Vdx = Sdx*Sdx;
  Vdy = Sdy*Sdy;
  
  for (i=1; i <= NDIM; i++){
    grad[i]=0.;
  }

  for (i=0 ; i<N ; i++ ){ /* loop over genes */

    mu_x = muX[i]; 
    mu_y = muY[i];

    /* construct C, the covariance matrix of x and y */    
    C11 = Vex* mu_x*mu_x + Vdx;
    C12 = Sex*Sey*Re* mu_x*mu_y + Sdx*Sdy*Rd;
    C22 = Vey* mu_y*mu_y + Vdy;
    detC = C11*C22 - C12*C12;

    /* construct B, the inverse of C */
    B11 = C22 / detC;
    B12 = -C12 / detC;
    B22 = C11 / detC;     

    /* Compute gradient of detC
       for s = { g1, ..., g6} :
       d(detC)/ds = d(C11)/ds * C22 + d(C22)/ds * C11 - 2 C12 * d(C12)/ds */

    ddetC_dg1 = 2.*Sex* mu_x*mu_x * C22 - 2.*Sey*Re *C12 *mu_x*mu_y ;
    ddetC_dg2 = 2.*Sey* mu_y*mu_y * C11 - 2.*Sex*Re *C12 *mu_x*mu_y ;  
    ddetC_dg3 = -2.*Sex*Sey * C12 * mu_x*mu_y; 

    /* Compute gradient of objective function 
       -1/(2 detC) * ( d(detC)/ds ( 1 - [x y]B[xy]) + x^2 d(C22)/ds +
       y^2 d(C11)/ds - 2 xy d(C12)/ds )  */

    temp1 = m[i] * ddetC_dg1;
    temp2 = m[i] * ddetC_dg2;
    temp3 = m[i] * ddetC_dg3;

    for (j=0 ; j<m[i]; j++ ){ /* loop over repeats */

      x=X[i][j];
      y=Y[i][j];

      temp = B11* x * x + B22 * y * y+ 2. * B12 * x * y;

      temp1 += -ddetC_dg1 * temp + 2.*Sex* mu_x*mu_x*y*y - 2.*Sey*Re* mu_x*mu_y*x*y ;
      temp2 += -ddetC_dg2 * temp + 2.*Sey* mu_y*mu_y*x*x - 2.*Sex*Re* mu_x*mu_y*x*y ;
      temp3 += -ddetC_dg3 * temp - 2.*Sex*Sey* mu_x*mu_y*x*y ;

    }

    grad[1] += temp1/detC;
    grad[2] += temp2/detC;
    grad[3] += temp3/detC;

  }

  for(i=1;i<=3;i++){
    grad[i] /= (double)2.;
  }

}


/********************************************************************************************/
/**                                                                                        **/
/**  gradZd.c - Gradient of objective function in beta[4,5,6] directions only             **/
/**             (delta contribution). (Returns gradient[1,2,3]=0. )                        **/
/**                                                                                        **/
/**  April 2000 Copyright (c) Vesteinn Thorsson and Trey Ideker                            **/
/********************************************************************************************/

/* X and Y must be shifted my mus */

void gradZd( double *beta, double *grad ){

  int i,j;
  double temp, temp4, temp5, x, y, mu_x, mu_y;
  double ddetC_dg4, ddetC_dg5; 
#ifdef USERHOD
  double temp6, ddetC_dg6; 
#endif

  Sex = beta[1];
  Sey = beta[2]; 
  Re  = beta[3];
  Sdx = beta[4];
  Sdy = beta[5];
#ifdef USERHOD  
  Rd  = beta[6];
#else
  Rd =0;
#endif  

  Vex = Sex*Sex; 
  Vey = Sey*Sey;
  Vdx = Sdx*Sdx;
  Vdy = Sdy*Sdy;
  
  for (i=1; i <= NDIM; i++){
    grad[i]=0.;
  }

  for (i=0 ; i<N ; i++ ){ /* loop over genes */

    mu_x = muX[i]; 
    mu_y = muY[i];

    /* construct C, the covariance matrix of x and y */    
    C11 = Vex* mu_x*mu_x + Vdx;
    C12 = Sex*Sey*Re* mu_x*mu_y + Sdx*Sdy*Rd;
    C22 = Vey* mu_y*mu_y + Vdy;
    detC = C11*C22 - C12*C12;

    /* construct B, the inverse of C */
    B11 = C22 / detC;
    B12 = -C12 / detC;
    B22 = C11 / detC;     

    /* Compute gradient of detC
       for s = { g1, ..., g6} :
       d(detC)/ds = d(C11)/ds * C22 + d(C22)/ds * C11 - 2 C12 * d(C12)/ds */

    ddetC_dg4 = 2.*Sdx * C22 - 2*Sdy*Rd* C12;
    ddetC_dg5 = 2.*Sdy * C11 - 2*Sdx*Rd* C12 ;
#ifdef USERHOD  
    ddetC_dg6 = -2.*Sdx*Sdy * C12;
#endif

    /* Compute gradient of objective function 
       -1/(2 detC) * ( d(detC)/ds ( 1 - [x y]B[xy]) + x^2 d(C22)/ds +
       y^2 d(C11)/ds - 2 xy d(C12)/ds )  */

    temp4 = m[i] * ddetC_dg4;
    temp5 = m[i] * ddetC_dg5;
#ifdef USERHOD
    temp6 = m[i] * ddetC_dg6;
#endif

    for (j=0 ; j<m[i]; j++ ){ /* loop over repeats */

      x=X[i][j];
      y=Y[i][j];

      temp = B11* x * x + B22 * y * y+ 2. * B12 * x * y;

      temp4 += -ddetC_dg4 * temp + 2.*Sdx* y*y - 2.*Sdy*Rd* x*y ;
      temp5 += -ddetC_dg5 * temp + 2.*Sdy* x*x - 2.*Sdx*Rd* x*y ; 
#ifdef USERHOD
      temp6 += -ddetC_dg6 * temp - 2.*Sdx*Sdy* x*y ; 
#endif

    }

    grad[4] += temp4/detC;
    grad[5] += temp5/detC;
#ifdef USERHOD
    grad[6] += temp6/detC;
#endif

  }

  for(i=4;i<=NDIM;i++){
    grad[i] /= (double)2.;
  }


}


/********************************************************************************************/
/**                                                                                        **/
/**  Z_gene.c - Single gene objective function (variables mu_x, mu_y)                      **/
/**                                                                                        **/
/**  April 2000 Copyright (c) Vesteinn Thorsson and Trey Ideker                            **/
/********************************************************************************************/

/* X (Y) must not have muX (muY) subtraction */

double Z_gene( double* mju ){

  int j; 
  double logsum, quadsum, x, y;

  Sex = cov[1];
  Sey = cov[2]; 
  Re  = cov[3];
  Sdx = cov[4];
  Sdy = cov[5];
#ifdef USERHOD  
  Rd  = cov[6];
#else
  Rd =0;
#endif

  Vex = Sex*Sex; 
  Vey = Sey*Sey;
  Vdx = Sdx*Sdx;
  Vdy = Sdy*Sdy;
  
  /* Define value at negative mu to be equal to positive */
#ifdef MUNONNEG
  mju[1] = fabs( mju[1] ); 
  mju[2] = fabs( mju[2] );
#endif

  /* Covariance matrix C and determinant */
  C11 = Vex* mju[1] * mju[1] + Vdx;
  C12 = Sex*Sey*Re * mju[1]*mju[2] + Sdx*Sdy*Rd;
  C22 = Vey* mju[2] * mju[2] + Vdy;
  detC = C11*C22 - C12*C12;
  if ( detC <= 0 ){
    printf("Error: Determinant of correlation matrix is not positive\n");
    exit(1);
  }

  logsum = m[gene_index] * log(detC); /* one contribution from each sample */

  /* construct B, the inverse of C */
  B11 = C22 / detC;
  B12 = -C12 / detC;
  B22 = C11 / detC;     
  
  quadsum=0.;
  for (j=0 ; j<m[gene_index]; j++ ){ /* loop over repeats */
    x = X[gene_index][j] - mju[1];
    y = Y[gene_index][j] - mju[2];
    quadsum += B11*x*x + B22*y*y + 2.*B12*x*y;
  }
  
  return ( m[gene_index]*log(2.*PI) + 0.5 * logsum + 0.5 * quadsum) ;
  
}

/********************************************************************************************/
/**                                                                                        **/
/**  gradZ_gene.c - Gradient of single gene objective function (variables mu_x, mu_y)      **/
/**                                                                                        **/
/**  April 2000 Copyright (c) Vesteinn Thorsson and Trey Ideker                            **/
/********************************************************************************************/

/* X (Y) must not have muX (muY) subtraction */

void gradZ_gene( double *mju, double *grad ){

  int j, s1, s2;
  double ddetC_dg1, ddetC_dg2, temp, x, y; 
  
  Sex = cov[1];
  Sey = cov[2]; 
  Re  = cov[3];
  Sdx = cov[4];
  Sdy = cov[5];
#ifdef USERHOD  
  Rd  = cov[6];
#else
  Rd =0;
#endif
  
  Vex = Sex*Sex; 
  Vey = Sey*Sey;
  Vdx = Sdx*Sdx;
  Vdy = Sdy*Sdy;

#ifdef MUNONNEG  
  /* take absolute values of mju[1] and mju[2] but retain sign to adjust gradient */
  if ( mju[1] >= 0 ){ s1 = 1 ; } else { s1 = -1; }
  if ( mju[2] >= 0 ){ s2 = 1 ; } else { s2 = -1; }
  mju[1] =  fabs( mju[1] );
  mju[2] =  fabs( mju[2] );
#endif

  /* construct C, the covariance matrix */
  C11 = Vex* mju[1]*mju[1] + Vdx;
  C12 = Sex*Sey*Re* mju[1]*mju[2] + Sdx*Sdy*Rd;
  C22 = Vey* mju[2]*mju[2] + Vdy;
  detC = C11*C22 - C12*C12;

  /* Compute gradient of detC
     for s = {mux, muy} :
     d(detC)/ds = d(C11)/ds * C22 + d(C22)/ds * C11 - 2 C12 * d(C12)/ds  */

  ddetC_dg1 = 2.*Vex*mju[1] * C22 - 2.* C12 *Sex*Sey*Re *mju[2] ; 
  ddetC_dg2 = 2.*Vey*mju[2] * C11 - 2.* C12 *Sex*Sey*Re *mju[1] ;

  /* construct B, the inverse of C */
  B11 = C22 / detC;
  B12 = -C12 / detC;
  B22 = C11 / detC;     

 /* Compute gradient of objective function 
     1/(2 detC) * ( d(detC)/ds ( 1 - [x y]B[xy]) + x^2 d(C22)/ds +
     y^2 d(C11)/ds - 2 xy d(C12)/ds ) + 1/2 * ( B11*d(x^2)/ds+ B22 d(y^2)/ds + 2 B12 d(xy)/ds ) */

  grad[1] = m[gene_index] * ddetC_dg1/detC/((double)2.);
  grad[2] = m[gene_index] * ddetC_dg2/detC/((double)2.);

  for (j=0 ; j<m[gene_index]; j++ ){ /* loop over repeats */

    x = X[gene_index][j] - mju[1];
    y = Y[gene_index][j] - mju[2];
    
    temp = B11*x*x + B22*y*y + 2.*B12*x*y;

    grad[1] += (-ddetC_dg1 * temp + 2.*y*y*Vex*mju[1] -2.*x*y*Sex*Sey*Re*mju[2])/detC/((double)2.);
    grad[2] += (-ddetC_dg2 * temp + 2.*x*x*Vey*mju[2] -2.*x*y*Sex*Sey*Re*mju[1])/detC/((double)2.);

    grad[1] -= B11*x+B12*y;
    grad[2] -= B22*y+B12*x; 

  }

#ifdef MUNONNEG  

  /* adjust gradient to reflect redefinition */

  grad[1] = s1 * grad[1] ; 
  grad[2] = s2 * grad[2] ; 

  /* special case : in trough, define gradient to be only in direction of trough */
  /* note that with above definitions, the derivative is define as the limit from above */
  /* if cusp, we can move away from the cusp i.e. use non-zero gradient */
 
  if ( (mju[1] >= 0 && mju[1] <=  1e-8 && grad[1] > 0.) ||
      ( mju[1] <= 0 && mju[1] >= -1e-8 && grad[1] < 0.) )  grad[1] = 0 ;
  if ( (mju[2] >= 0 && mju[2] <=  1e-8 && grad[2] > 0.) ||
       (mju[2] <= 0 && mju[2] >= -1e-8 && grad[2] < 0.) )  grad[2] = 0 ;

  /* threshold was set emperically */
  /* For the conjugate gradient method, this is not an ultimate solution, as the search */
  /* direction is not always the gradient.  Re-initializing might solve this. */

#endif 

}

/********************************************************************************************/
/**                                                                                        **/
/**  hessianZ_gene.c - Hessian of single gene objective function (variables mu_x, mu_y)    **/
/**                                                                                        **/
/**  The Hessian is returned as a three-vector (H11,H22,H12=H21)                           **/
/**  June 2002  Copyright (c) Vesteinn Thorsson and Andrew Siegel                          **/
/********************************************************************************************/

/* X (Y) must not have muX (muY) subtraction */

void hessianZ_gene( double *mju, double *hessian ){

  int j, s1, s2;
  double temp, x, y; 
  double B11_1,B11_2,B22_1,B22_2,B12_1,B12_2 ;
  double C11_1, C22_2,C12_1,C12_2,detC_1, detC_2; 
  double C11_11,C22_22,C12_12,detC_11,detC_22,detC_12;
  double H11C,H22C,H12C;
  double B11_11,B11_22,B11_12,B22_11,B22_22,B22_12,B12_11,B12_22,B12_12;
  double detC_sq,detC_cub, detC_1_sq,detC_2_sq;
  double H11B,H22B,H12B;
  
  Sex = cov[1];
  Sey = cov[2]; 
  Re  = cov[3];
  Sdx = cov[4];
  Sdy = cov[5];
#ifdef USERHOD  
  Rd  = cov[6];
#else
  Rd =0;
#endif
  
  Vex = Sex*Sex; 
  Vey = Sey*Sey;
  Vdx = Sdx*Sdx;
  Vdy = Sdy*Sdy;

#ifdef MUNONNEG  
  /* take absolute values of mju[1] and mju[2] but retain sign to adjust gradient */
  if ( mju[1] >= 0 ){ s1 = 1 ; } else { s1 = -1; }
  if ( mju[2] >= 0 ){ s2 = 1 ; } else { s2 = -1; }
  mju[1] =  fabs( mju[1] );
  mju[2] =  fabs( mju[2] );
#endif

  /* construct C, the covariance matrix */
  C11 = Vex* mju[1]*mju[1] + Vdx;
  C12 = Sex*Sey*Re* mju[1]*mju[2] + Sdx*Sdy*Rd;
  C22 = Vey* mju[2]*mju[2] + Vdy;
  detC = C11*C22 - C12*C12;

  /* construct B, the inverse of C */
  B11 = C22 / detC;
  B12 = -C12 / detC;
  B22 = C11 / detC;     

  //Gradient of detC
  //  for s = {mux, muy} :
  // d(detC)/ds = d(C11)/ds * C22 + d(C22)/ds * C11 - 2 C12 * d(C12)/ds  */
  C11_1 = 2.*Vex*mju[1];
  C22_2 = 2.*Vey*mju[2];
  C12_1 = Sex*Sey*Re*mju[2];
  C12_2 = Sex*Sey*Re*mju[1];
  detC_1 = C11_1 * C22 - 2.*C12 * C12_1; 
  detC_2 = C22_2 * C11 - 2.*C12 * C12_2; 

  detC_sq = pow(detC,2.);
  detC_cub = pow(detC,3.);
  detC_1_sq = pow(detC_1,2.);
  detC_2_sq = pow(detC_2,2.);

  //Quadratic (B) form contribution
  B11_1 = - C22 / detC_sq * detC_1 ;
  B11_2 = C22_2 / detC - C22 / detC_sq * detC_2 ;
  B22_1 = C11_1 / detC - C11 / detC_sq * detC_1 ;
  B22_2 = - C11 / detC_sq * detC_2 ;
  B12_1 = - C12_1 / detC + C12 / detC_sq * detC_1 ;
  B12_2 = - C12_2 / detC + C12 / detC_sq * detC_2 ; 

  // Notation:
  // X_1 denotes partial derivative of X w.r.t mu_x 
  // X_2 denotes partial derivative of X w.r.t mu_y
  // X_s here, s can be one or two
  // X_st, X_11, X_12 examples of partial derivatives
  // C12_12: the mu_x, mu_y derivative of matrix element C12 !
  //
  // H11C : The log(detC) contribution to the Hessian element 1,1
  // H12B: The B quadratic form contrib. to Hessian element 1,2 
  //
  // The overall factor of 1/2 is left until the end
  //
  
  // detC derivatives
  // Diagonal derivatives
  // detC_ss = C22 C11_ss + C11 C22_ss + 2 C11_s C22_s - 2 (C12_s)^2 - 2 C12 C12_ss
  // Mixed derivatives 
  // detC_st = C22 C11_st + C11 C22_st + C11_s C22_t + C11_t C22_s 
  //         - 2 C12_s C12_t - 2 C12 C12_st 

  C11_11 = 2.*Vex; 
  C22_22 = 2.*Vey; 
  C12_12 = Sex*Sey*Re ;
  
  detC_11 = C22*C11_11 - 2.*pow(C12_1,2) ;
  detC_22 = C11*C22_22 - 2.*pow(C12_2,2) ;
  detC_12 = C11_1*C22_2 - 2*C12_1*C12_2 - 2*C12*C12_12 ; 

  // log(detC) contribution
  // Diagonal derivatives 
  // (log(detC))_ss = -1/(detC)^2 ( detC_s )^2 + 1/detC (detC)_ss
  // Mixed derivatives 
  // (log(detC))_st = -1/(detC)^2 detC_s detC_t + 1/detC (detC)_st
  
  H11C = -1./pow(detC,2) * pow( detC_1,2 ) + 1./detC * detC_11 ;
  H22C = -1./pow(detC,2) * pow( detC_2,2 ) + 1./detC * detC_22 ;
  H12C = -1./pow(detC,2) * detC_1 * detC_2 + 1./detC * detC_12 ;
  H11C = m[gene_index]*H11C ;
  H22C = m[gene_index]*H22C ;
  H12C = m[gene_index]*H12C ;

  // 
  // Quadratic, B, contribution
  // 
  // B matrix double derivatives
  // For B_ij = a/detC    (a=C22,C11,-C12 for ij=11,22,12)
  // Diagonal 
  // (a/detC)_ss = a_ss/detC - 2 a_s/(detC)^2 detC_s 
  //             + 2 a/(detC)^3 (detC_s)^2 - 2 a/(detC)^2 detC_ss
  // Mixed 
  // (a/detC)_st = a_st/detC - a_s/(detC)^2 detC_t - a_t/(detC)^2 detC_s
  //             + 2 a/(detC)^3 detC_s detC_t - 2 a/(detC)^2 detC_st

  
  B11_11 = 2*C22/detC_cub*detC_1_sq - C22/detC_sq*detC_11 ;
  B11_22 = C22_22/detC - 2*C22_2/detC_sq*detC_2 +
    2*C22/detC_cub*detC_2_sq - C22/detC_sq*detC_22 ;
  B11_12 = - C22_2/detC_sq*detC_1 + 
    2*C22/detC_cub*detC_1*detC_2 - C22/detC_sq*detC_12 ;
  
  B22_11 = C11_11/detC - 2*C11_1/detC_sq*detC_1+
    2*C11/detC_cub*detC_1_sq - C11/detC_sq*detC_11 ;
  B22_22 = 2*C11/detC_cub*detC_2_sq - C11/detC_sq*detC_22 ;
  B22_12 = - C11_1/detC_sq*detC_2 + 
    2*C11/detC_cub*detC_1*detC_2 - C11/detC_sq*detC_12 ;
  
  B12_11 = 2*C12_1/detC_sq*detC_1- 
    2*C12/detC_cub*detC_1_sq + C12/detC_sq*detC_11;
  B12_22 = 2*C12_2/detC_sq*detC_2- 
    2*C12/detC_cub*detC_2_sq + C12/detC_sq*detC_22;
  B12_12 = - C12_12/detC + C12_1/detC_sq*detC_2 + C12_2/detC_sq*detC_1 -
    2*C12/detC_cub*detC_1*detC_2 + C12/detC_sq*detC_12 ; 
  

  H11B=0; 
  H22B=0;
  H12B=0; 

  for ( j=0 ; j<m[gene_index] ; j++ ){

    x = X[gene_index][j] - mju[1];
    y = Y[gene_index][j] - mju[2];
	 
    H11B += B11_11*x*x + B22_11*y*y + 2.*B12_11*x*y 
      - 4*B11_1*x - 4*B12_1*y + 2*B11 ;
    
    H22B += B11_22*x*x + B22_22*y*y + 2.*B12_22*x*y 
      - 4*B22_2*y - 4*B12_2*x + 2*B22 ; 
    
    H12B += B11_12*x*x + B22_12*y*y + 2.*B12_12*x*y 
      - 2*B11_2*x - 2*B12_2*y - 2*B22_1*y - 2*B12_1*x + 2*B12 ; 
  }

  // Combine contributions and include overall 1/2 factor

  hessian[1] = ( H11C + H11B ) / 2.;  
  hessian[2] = ( H22C + H22B ) / 2.;  
  hessian[3] = ( H12C + H12B ) / 2.;  

#ifdef MUNONNEG  
  /* adjust gradient to reflect redefinition */
  hessian[1] = s1 * s1 * hessian[1] ; 
  hessian[2] = s2 * s2 * hessian[2] ; 
  hessian[3] = s1 * s2 * hessian[3] ; 
#endif 

}


/********************************************************************************************/
/**                                                                                        **/
/**  Z_gene_constrained.c - Single gene objective function                                 **/
/**                                (constrained variables mju=mu_x=mu_y)                   **/
/**                                                                                        **/
/**  April 2000 Copyright (c) Vesteinn Thorsson and Trey Ideker                            **/
/********************************************************************************************/

/* X (Y) must not have muX (muY) subtraction */

double Z_gene_constrained( double mju ){

  int j; 
  double logsum, quadsum, x, y;

  Sex = cov[1];
  Sey = cov[2]; 
  Re  = cov[3];
  Sdx = cov[4];
  Sdy = cov[5];
#ifdef USERHOD  
  Rd  = cov[6];
#else
  Rd =0;
#endif
  
  Vex = Sex*Sex; 
  Vey = Sey*Sey;
  Vdx = Sdx*Sdx;
  Vdy = Sdy*Sdy;

#ifdef MUNONNEG
  /* Take absolute value of mu for redefinition */
  mju = fabs( mju );
#endif
  
  /* Covariance matrix C and determinant */
  C11 = Vex* mju * mju + Vdx;
  C12 = Sex*Sey*Re * mju*mju + Sdx*Sdy*Rd;
  C22 = Vey* mju * mju + Vdy;
  detC = C11*C22 - C12*C12;
  if ( detC <= 0 ){
    printf("Error: Determinant of correlation matrix is not positive\n");
    exit(1);
  }

  logsum = m[gene_index] * log(detC); /* one contribution from each sample */

  /* construct B, the inverse of C */
  B11 = C22 / detC;
  B12 = -C12 / detC;
  B22 = C11 / detC;     
  
  quadsum=0.;
  for (j=0 ; j<m[gene_index]; j++ ){ /* loop over repeats */
    x = X[gene_index][j] - mju;
    y = Y[gene_index][j] - mju;
    quadsum += B11*x*x + B22*y*y + 2.*B12*x*y;
  }
  
  return ( m[gene_index]*log(2.*PI) + 0.5 * logsum + 0.5 * quadsum) ;
  
}

/********************************************************************************************/
/**                                                                                        **/
/**  Z_gene_12_xu_yc                                                                       **/
/**                Single gene objective function (variables mu_x1, mu_x2, mu_y1= mu_y2)   **/
/**                                                                                        **/
/**  Jan 2001   Copyright (c) Vesteinn Thorsson and Trey Ideker                            **/
/********************************************************************************************/

/* X (Y) must not have muX (muY) subtraction */

double Z_gene_12_xu_yc (double* mu_vec  ){

  int i;
  double mu_x1, mu_x2, mu_y, term1, term2; 
  double mu[3] ; 

  /* allocate for 1-gene data matrix, max 10 samples  */
  X = (float **) calloc(1,sizeof(float*));
  Y = (float **) calloc(1,sizeof(float*));
  X[0] = (float *) calloc( 10 ,sizeof(float));      
  Y[0] = (float *) calloc( 10 ,sizeof(float));  	 

  mu_x1 = mu_vec[1]; 
  mu_x2 = mu_vec[2]; 
  mu_y  = mu_vec[3]; 

  /* use Z_gene by first redefining cov, gene_index, *mu, X and Y, and m */

  gene_index = 0 ; /* create data array consisting of one entry */
  
  /* contribution from first data set */

  for (i=1 ; i<= NDIM ; i++ ) cov[i] = cov1[i] ; 
  mu[1] = mu_x1;
  mu[2] = mu_y ; 

  m[0] = m1[gene_index_1]; 

  for ( i=0 ; i<m1[gene_index_1] ; i++ ){
    X[0][i] = X1[gene_index_1][i];
    Y[0][i] = Y1[gene_index_1][i];
  }
  
  term1 = Z_gene ( mu ) ; 


  /* contribution from second data set */

  for (i=1 ; i<= NDIM ; i++ ) cov[i] = cov2[i] ; 
  mu[1] = mu_x2;
  mu[2] = mu_y ; 

  m[0] = m2[gene_index_2]; 
  for ( i=0 ; i<m2[gene_index_2] ; i++ ){
    X[0][i] = X2[gene_index_2][i];
    Y[0][i] = Y2[gene_index_2][i];
  }
  
  term2 = Z_gene ( mu ) ; 

  free(X[0]); free(Y[0]); free(X) ; free(Y); 

  /* printf("%12.3f %12.3f\n", term1, term2 );  */

  return (term1 + term2) ; 

}


/********************************************************************************************/
/**                                                                                        **/
/**  grad_Z_gene_12_xu_yc Gradient of :                                                    **/
/**                Single gene objective function (variables mu_x1, mu_x2, mu_y1= mu_y2)   **/
/**                                                                                        **/
/**  Jan 2001   Copyright (c) Vesteinn Thorsson and Trey Ideker                            **/
/********************************************************************************************/

/* X (Y) must not have muX (muY) subtraction */


void gradZ_gene_12_xu_yc (double* mu_vec, double *grad ){

  int i;
  double mu_x1, mu_x2, mu_y; 
  double mu[3] ; 
  float dZ1_dmu_x1, dZ1_dmu_y, dZ2_dmu_x2, dZ2_dmu_y ; 
  double *gradient; 

  /* allocate for 1-gene data matrix, max 10 samples  */
  X = (float **) calloc(1,sizeof(float*));
  Y = (float **) calloc(1,sizeof(float*));
  X[0] = (float *) calloc( 10 ,sizeof(float));      
  Y[0] = (float *) calloc( 10 ,sizeof(float));  	 
  gradient = (double*) calloc (3, sizeof (double)); /* for use with grad_Z */

  mu_x1 = mu_vec[1]; 
  mu_x2 = mu_vec[2]; 
  mu_y  = mu_vec[3]; 
  gene_index = 0 ; 

  /* compute gradient contributions from set 1 */

  for (i=1 ; i<= NDIM ; i++ ) cov[i] = cov1[i] ; 
  mu[1] = mu_x1;
  mu[2] = mu_y ; 

  m[0] = m1[gene_index_1]; 
  for ( i=0 ; i<m1[gene_index_1] ; i++ ){
    X[0][i] = X1[gene_index_1][i];
    Y[0][i] = Y1[gene_index_1][i];
  }
  
  gradZ_gene( mu, gradient ); 
  dZ1_dmu_x1 = gradient[1] ; 
  dZ1_dmu_y  = gradient[2] ; 

  /* compute gradient contributions from set 2 */

  for (i=1 ; i<= NDIM ; i++ ) cov[i] = cov2[i] ; 
  mu[1] = mu_x2;
  mu[2] = mu_y ; 

  m[0] = m2[gene_index_2]; 
  for ( i=0 ; i<m2[gene_index_2] ; i++ ){
    X[0][i] = X2[gene_index_2][i];
    Y[0][i] = Y2[gene_index_2][i];
  }
  
  gradZ_gene ( mu, gradient );
  dZ2_dmu_x2 = gradient[1] ; 
  dZ2_dmu_y  = gradient[2] ;

  free(X[0]); free(Y[0]); free(X) ; free(Y); 

  /* Combine to form total gradient */

  /* (del Z)/(del mu_x1) = (del Z1)/(del mu_x1) */  

  grad[1] = dZ1_dmu_x1 ; 

  /* (del Z)/(del mu_x2) = (del Z2)/(del mu_x2) */  

  grad[2] = dZ2_dmu_x2 ; 

  /* (del Z)/(del mu_y) = (del Z1)/(del mu_y) + (del Z2)/(del mu_y) */  

  grad[3] = dZ1_dmu_y + dZ2_dmu_y;

}

/********************************************************************************************/
/**                                                                                        **/
/**  Z_gene_12_xc_yc                                                                       **/
/**                Single gene objective function (variables mu_x1 =mu_x2, mu_y1= mu_y2)   **/
/**                                                                                        **/
/**  Jan 2001   Copyright (c) Vesteinn Thorsson and Trey Ideker                            **/
/********************************************************************************************/

/* X (Y) must not have muX (muY) subtraction */

double Z_gene_12_xc_yc (double* mu_vec ){

  int i;
  double mu_x , mu_y; 
  double mu[3] ; 
  double term1, term2; 

  /* allocate for 1-gene data matrix, max 10 samples  */
  X = (float **) calloc(1,sizeof(float*));
  Y = (float **) calloc(1,sizeof(float*));
  X[0] = (float *) calloc( 10 ,sizeof(float));      
  Y[0] = (float *) calloc( 10 ,sizeof(float));  	 

  mu_x = mu_vec[1]; 
  mu_y = mu_vec[2]; 

  /* use Z_gene by first redefining cov, gene_index, *mu, X and Y, and m */

  gene_index = 0 ; /* create data array consisting of one entry */
  
  /* contribution from first data set */

  for (i=1 ; i<= NDIM ; i++ ) cov[i] = cov1[i] ; 
  mu[1] = mu_x;
  mu[2] = mu_y; 

  m[0] = m1[gene_index_1]; 
  for ( i=0 ; i<m1[gene_index_1] ; i++ ){
    X[0][i] = X1[gene_index_1][i];
    Y[0][i] = Y1[gene_index_1][i];
  }
  
  term1 = Z_gene ( mu ) ; 

  /* contribution from second data set */

  for (i=1 ; i<= NDIM ; i++ ) cov[i] = cov2[i] ; 
  mu[1] = mu_x;
  mu[2] = mu_y; 

  m[0] = m2[gene_index_2]; 
  for ( i=0 ; i<m2[gene_index_2] ; i++ ){
    X[0][i] = X2[gene_index_2][i];
    Y[0][i] = Y2[gene_index_2][i];
  }
  
  term2 = Z_gene ( mu ) ; 

  free(X[0]); free(Y[0]); free(X) ; free(Y); 

  return (term1 + term2) ; 

}


/********************************************************************************************/
/**                                                                                        **/
/**  grad_Z_gene_12_xc_yc Gradient of :                                                    **/
/**                Single gene objective function (variables mu_x1=mu_x2, mu_y1= mu_y2)   **/
/**                                                                                        **/
/**  Jan 2001   Copyright (c) Vesteinn Thorsson and Trey Ideker                            **/
/********************************************************************************************/

/* X (Y) must not have muX (muY) subtraction */


void gradZ_gene_12_xc_yc (double* mu_vec, double *gradd ){

  int i;
  double mu_x , mu_y; 
  double mu[3] ; 
  float dZ1_dmu_x, dZ1_dmu_y, dZ2_dmu_x, dZ2_dmu_y ; 
  double *gradient; /* gradient for use with Z_gene */

  /* allocate for 1-gene data matrix, max 10 samples  */
  X = (float **) calloc(1,sizeof(float*));
  Y = (float **) calloc(1,sizeof(float*));
  X[0] = (float *) calloc( 10 ,sizeof(float));      
  Y[0] = (float *) calloc( 10 ,sizeof(float));

  gradient = (double *) calloc( 3, sizeof(double) ); 

  mu_x = mu_vec[1]; 
  mu_y = mu_vec[2]; 

  /* use Z_gene by first redefining cov, gene_index, *mu, X and Y, and m */

  gene_index = 0 ; /* create data array consisting of one entry */

  /* compute gradient contributions from set 1 */

  for (i=1 ; i<= NDIM ; i++ ) cov[i] = cov1[i] ; 
  mu[1] = mu_x ;
  mu[2] = mu_y ; 

  m[0] = m1[gene_index_1]; 
  for ( i=0 ; i<m1[gene_index_1] ; i++ ){
    X[0][i] = X1[gene_index_1][i];
    Y[0][i] = Y1[gene_index_1][i];
  }
  
  gradZ_gene( mu, gradient ); 
  dZ1_dmu_x  = gradient[1] ; 
  dZ1_dmu_y  = gradient[2] ; 

  /* compute gradient contributions from set 2 */

  for (i=1 ; i<= NDIM ; i++ ) cov[i] = cov2[i] ; 
  mu[1] = mu_x ;
  mu[2] = mu_y ; 

  m[0] = m2[gene_index_2]; 
  for ( i=0 ; i<m2[gene_index_2] ; i++ ){
    X[0][i] = X2[gene_index_2][i];
    Y[0][i] = Y2[gene_index_2][i];
  }
  
  gradZ_gene ( mu, gradient );
  dZ2_dmu_x  = gradient[1] ; 
  dZ2_dmu_y  = gradient[2] ;

  free(X[0]); free(Y[0]); free(X) ; free(Y); 

  /* Combine to form total gradient */

  /* (del Z)/(del mu_x) = (del Z1)/(del mu_x) + (del Z2)/(del mu_x) */  

  gradd[1] = dZ1_dmu_x + dZ2_dmu_x;

  /* (del Z)/(del mu_y) = (del Z1)/(del mu_y) + (del Z2)/(del mu_y) */  

  gradd[2] = dZ1_dmu_y + dZ2_dmu_y;


}
