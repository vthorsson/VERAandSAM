/********************************************************************************************/
/**                                                                                        **/
/**  SAM.c - Given parameters of covariant error model, find                               **/
/**       significantly differentially expressed genes in two-channel expression data set  **/
/**                                                                                        **/    
/**                 evaluates lambda =                                                     **/
/**      -2*Log( (max unconstrained likelihood)/(max constrained likelihood) )             **/
/**         constrained: mu=mu_x=mu_y   unconstrained mu_x != mu_y                         **/
/**                                                                                        **/
/** Reference:                                                                             **/ 
/** T. E. Ideker, V. Thorsson, A. F. Siegel, and L. E. Hood                                **/   
/** Testing for Differentially-expressed Genes by                                          **/
/** Maximum Likelihood Analysis of Microarray Data                                         **/ 
/** Journal of Computational Biology, 2001, Vol. 7, No. 6, pps. 805-817                    **/ 
/** Implementation: Vesteinn Thorsson thorsson@systemsbiology.org                          **/
/** Version: February, 2001                                                                **/
/********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "util.h" 
#include "arraystats.h"
#include "objective.h"
#include "conjugate_gradient.h"
#include "io.h"
#define NRANSI

#define FTOL_HIGH 1.e-8  /* high fractional functional tolerance change for frprmn */
#define FTOL_LOW 1.e-10  /* low fractional functional tolerance change for frprmn */
#define LINMIN_TOL_HIGH  1.e-9 /* high tolerance for line minimization */
#define LINMIN_TOL_LOW  1.e-10 /* low tolerance for line minimization */

main (int argc,char *argv[] ){

  int i, j, k, iter, num_sig=0;
  double result, fret;
  double *p, *del; 
  char *infile, *modfile, *outfile;
  FILE *fpin, *fpmod, *fpout;
  double *mu_vec, *grad_mu;
  double mu_x, mu_y;
  double frac_change_1, frac_change_2 ; 
  double ax, bx, xx, fa, fx, fb, xmin, xmin2;
  double numerator, denominator, llr, *loglr, numerator2;
  double mu_x_from_avg, mu_y_from_avg, mu_x_from_cns, mu_y_from_cns;
  double denominator_from_avg, denominator_from_cns;

  /* output version number and compilation parameters  */
  /* printf("SAM version: %s\n", __DATE__); */
#ifdef USERHOD
  printf("rho_delta may be non-zero, ");
#else
  printf("rho_delta forced to zero, ");
#endif
#ifdef MUNONNEG
  printf("non-negative mus only.\n ");
#else
  printf("mus may be negative.\n ");
#endif

  /* initialize command line variables */
  verbose = FALSE;
  disp_iter = FALSE;        
  infile = NULL; 
  modfile = NULL;
  outfile = NULL;
  /* start_time(); start timer */

  /* parse command line, check to see if verbose and/or display iterations flag is on  */
  for (i=1; i<argc; i++) {
    if (strcmp(argv[i], "-verbose") == 0) verbose = TRUE ;
    else if (strcmp(argv[i], "-iter") == 0) disp_iter = TRUE;
    else if (infile == NULL) infile = argv[i];
    else if (modfile == NULL) modfile = argv[i];
    else if (outfile == NULL) outfile = argv[i];
    else die (5) ;
  }
  if (infile == NULL || modfile == NULL || outfile == NULL ) die(5); 

  fpin = fopen(infile, "r");
  if (fpin == NULL) die (2);
  fpmod = fopen(modfile, "r");
  if (fpmod == NULL) die (2); 
  fpout = fopen(outfile, "w");
  if (fpout == NULL) die (2);

  /* read model */
  printf ("Reading model file %s\n",  modfile ); 
  cov = (double *) calloc(NDIM+1, sizeof(double));  
  read_mod( fpmod );

  /* read data */
  printf ("Reading data file  %s\n",  infile ); 
  read_data( fpin ); 
  loglr = (double *) calloc(N, sizeof(double));

  /* allocate memory for mu variables and set them to zero */
  muX = (float *) calloc(N, sizeof(float));
  muY = (float *) calloc(N, sizeof(float));  
  for (i=0 ; i<N ; i++ ){
    muX[i] = 0.; 
    muY[i] = 0.;
  }

  /* Compute means as initial guesses for mus  */
  for(i=0 ; i<N ; i++ ){
    for (j=0; j<m[i] ; j++ ){
      muX[i] += X[i][j] ;
      muY[i] += Y[i][j] ;
    }
    muX[i] /= (float)(m[i]);
    muY[i] /= (float)(m[i]);
  }

  /****************************************************************************/
  /*********       Memory allocation, Initialization                     ******/
  /****************************************************************************/

  /* parameters of covariance matrix and vectors for single gene objective function */
  mu_vec = (double *) calloc(NDIM_MU+1, sizeof(double));
  /*grad_mu = (double *) calloc(NDIM_MU+1, sizeof(double));  */

  /*read covariance and mus (here, means) */
  
  for(i=0; i<N; i++ ){ 

    mu_x = (double) muX[i] ;
    mu_y = (double) muY[i] ;
    gene_index = i;

    if( disp_iter == TRUE )
      printf("++++++gene %d (%s,%s)+++++++\n",internal2full[i], yorf[i], gene_name[i]);
  
    /******** numerator ***********/
    LINMIN_TOL = LINMIN_TOL_LOW;
    ax = muX[i];
    xx = muY[i];

    mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,Z_gene_constrained);
    numerator=brent(ax,xx,bx,Z_gene_constrained,LINMIN_TOL,&xmin);
#ifdef MUNONNEG
    xmin = fabs( xmin );
#endif 
    if ( disp_iter ==  TRUE )
      printf("constrained mu %f, numerator %f\n", xmin, numerator); 

    /******* denominator, start search from averages   **********/
    
    mu_vec[1] = muX[i];
    mu_vec[2] = muY[i];

    if ( disp_iter == TRUE )
      printf("mu_vec before minimization %f %f\n", mu_vec[1], mu_vec[2]);
    LINMIN_TOL = LINMIN_TOL_LOW;
    frprmn(mu_vec, NDIM_MU, FTOL_LOW, &iter, &denominator, Z_gene, gradZ_gene); 
#ifdef MUNONNEG
    mu_vec[1] = fabs( mu_vec[1] );
    mu_vec[2] = fabs( mu_vec[2] );
#endif 
    if ( disp_iter == TRUE ){
      printf("mu_vec after minimization  %f %f\n", mu_vec[1], mu_vec[2]); 
      printf("denominator %f\n", denominator); 
    }

    mu_x_from_avg = mu_vec[1] ; 
    mu_y_from_avg = mu_vec[2] ; 
    denominator_from_avg = denominator ; 


    /******* denominator, start search from constrained mu  **********/
    
    mu_vec[1] = xmin;
    mu_vec[2] = xmin;

    if ( disp_iter == TRUE )
      printf("mu_vec before minimization %f %f\n", mu_vec[1], mu_vec[2]);
    LINMIN_TOL = LINMIN_TOL_LOW;
    frprmn(mu_vec, NDIM_MU, FTOL_LOW, &iter, &denominator, Z_gene, gradZ_gene); 
#ifdef MUNONNEG
    mu_vec[1] = fabs( mu_vec[1] );
    mu_vec[2] = fabs( mu_vec[2] );
#endif
    if ( disp_iter == TRUE ){
      printf("mu_vec after minimization  %f %f\n", mu_vec[1], mu_vec[2]); 
      printf("denominator %f\n", denominator); 
    }

    mu_x_from_cns = mu_vec[1] ; 
    mu_y_from_cns = mu_vec[2] ; 
    denominator_from_cns = denominator ; 

    if ( denominator_from_avg <= denominator_from_cns ){
      mu_vec[1] = mu_x_from_avg ; 
      mu_vec[2] = mu_y_from_avg ; 
      denominator = denominator_from_avg ; 
    } else { 
      mu_vec[1] = mu_x_from_cns ; 
      mu_vec[2] = mu_y_from_cns ; 
      denominator = denominator_from_cns ; 
    }

    denominator = Z_gene( mu_vec ) ; 

    /* replace revised muX, muY */
    muX[i] = mu_vec[1];
    muY[i] = mu_vec[2];

    /******** if numerator changes with revised guess, replace it ***********/
    LINMIN_TOL = LINMIN_TOL_LOW;
    ax = muX[i];
    xx = muY[i];

    mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,Z_gene_constrained);
    numerator2=brent(ax,xx,bx,Z_gene_constrained,LINMIN_TOL,&xmin2);
#ifdef MUNONNEG
    xmin2 = fabs( xmin2 );
#endif
    if ( disp_iter ==  TRUE )
      printf("constrained mu %f, numerator %f\n", xmin2, numerator); 

    if ( numerator2 <= numerator ) {
      numerator = numerator2 ; 
      xmin = xmin2 ;
    }

    if ( disp_iter == TRUE ) {
      printf("Final constrained mu %f , numerator %f\n", xmin, numerator);       
    }

    /* compute log likelihood ratio */ 
    llr = 2.*(numerator-denominator);
    if( llr < 0 ) llr = 0. ;  /* a few genes have slightly negative lambda, due to local minima. For now set to 0 */ 
    if( llr > 5000. ) llr = 5000. ;  /* single instance observed of astronomical lambda, for now set to maximum for warning */
    loglr[i] = llr;
    /* if ( llr > 3.84145882 ) num_sig++;  */
  }
  
  write_llr( fpin, fpout, loglr);

  fclose(fpmod);
  free(loglr);
  
  return 0;

}





