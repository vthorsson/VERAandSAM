
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

