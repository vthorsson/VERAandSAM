/********************************************************************************************/
/**                                                                                        **/
/**  VERA_OC.c - Optimize least squares fit to paramaters of an error model                **/
/**                 distribution from one-channel expression data                          **/
/**                                                                                        **/
/**                 Error model                                                            **/
/**                 X = mu_x epsilon_x + delta_x                                           **/
/**                                                                                        **/
/**                           Normal distribution parameters                               **/
/**                 beta[1] = sigma_epsilon_x                                              **/
/**                 beta[2] = sigma_delta_x                                                **/
/**                                                                                        **/
/** Reference:                                                                             **/ 
/** T. E. Ideker, V. Thorsson, A. F. Siegel, and L. E. Hood                                **/   
/** Testing for Differentially-expressed Genes by                                          **/
/** Maximum Likelihood Analysis of Microarray Data                                         **/ 
/** Journal of Computational Biology, 2001, Vol. 7, No. 6, pp.805-817                      **/ 
/** Implementation: Vesteinn Thorsson thorsson@systemsbiology.org                          **/
/** Version: March, 2002                                                                   **/
/********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "util.h" 
#include "arraystats.h"
#include "objective.h"
#include "conjugate_gradient.h"
#include "io.h"
#define NRANSI

#define FTOL_HIGH 1.e-6 /* high fractional functional tolerance change for frprmn */
#define FTOL_LOW 1.e-10  /* low fractional functional tolerance change for frprmn */
#define LINMIN_TOL_HIGH  1.e-6 /* high tolerance for line minimization */
#define LINMIN_TOL_LOW  1.e-10 /* low tolerance for line minimization */

#define MAX_MASTER 250 /* maximum number of master loop iterations */
#define MIN_MASTER 5 /* minimum number of master loop iterations */
#define MAX_QUICK 10 /* number of high tolerance subloop iterations */
double EPS_BETA=0.0005; /* stop criterion:  abs(beta_increment) < EPS_BETA */

main ( int argc, char *argv[] ){

  int i, j, k, l, iter, internal_index;
  double fret;
  double *p, *del, *p_old, abs_dev, Z_old, Z_diff; 
  char *in_data, *in_model, *out_model; 
  FILE *fp_in_data, *fp_in_model, *fp_out_model, *fp_evolution;
  char evolution_file_name[100] ;
  double *mu_vec, *grad_mu;
  double mu_x, mu_y;
  double FTOL; 
  double mu_scalar; 

  double ax, bx, xx, fa, fx, fb, xmin ;


  /* output version number and compilation parameters  */
  printf("VERA_OC version: %s\n", __DATE__);
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
  start_mod = FALSE;
  evolution = FALSE; 
  stop_crit = FALSE;
  in_data   = NULL;
  in_model  = NULL;
  out_model = NULL;
  
  start_time(); /* start timer */

  /* parse command line */
  for (i=1; i<argc; i++) {

    if (strcmp(argv[i], "-verbose") == 0) verbose = TRUE ;
    else if (strcmp(argv[i], "-iter") == 0) disp_iter = TRUE;
    else if (strcmp(argv[i], "-init") == 0) start_mod = TRUE;
    else if (strcmp(argv[i], "-evol") == 0) evolution = TRUE;

    else if (in_model  == NULL && start_mod == TRUE ) in_model  = argv[i];

    else if (strcmp(argv[i], "-crit") == 0) stop_crit = TRUE; 
    else if (stop_crit == TRUE ){ EPS_BETA  = atof(argv[i]); stop_crit = FALSE; }

    else if (in_data   == NULL) in_data   = argv[i];
    else if (out_model == NULL) out_model = argv[i];
    else die (1) ;
  }

  if (in_data == NULL || (in_model == NULL && start_mod == TRUE) || out_model == NULL ) die(1);

  fp_in_data = fopen(in_data, "r");
  if (fp_in_data == NULL) die (2);
  if ( start_mod == TRUE ){
    fp_in_model = fopen(in_model, "r");
    if (fp_in_model == NULL) die (2);
  }
  fp_out_model = fopen(out_model, "w");
  if (fp_out_model == NULL) die (2);

  if ( evolution == TRUE ){ 
    strcpy ( evolution_file_name, in_data );
    strcat ( evolution_file_name, ".evolution" );
    fp_evolution = fopen( evolution_file_name, "w" ); 
    fprintf(fp_evolution, "iter  sigma_eps_x sigma_del_x    objective function\n");
    fflush( fp_evolution );
  }
    
  /* read in the data */
  printf ("Reading expression data (merge) file  %s\n",  in_data ); 
  read_data_OC( fp_in_data ); 

  /* allocate memory for mu variables and set them to zero */
  muX = (float *) calloc(N, sizeof(float));
  for (i=0 ; i<N ; i++ ){
    muX[i] = 0.; 
  }

  /* read in the model */
  cov = (double *) calloc(NDIM+1, sizeof(double));
  if ( start_mod == TRUE ) {
    printf ("Reading input model (initial guess) file %s\n",  in_model ); 
    read_mod( fp_in_model );
  } else {  /* if model is not read in use starting guess */
    cov[1] = 0.3 ; 
    cov[2] = 20 ; 
    /* cov[2] = 10 ; */
  }

  /* If mu values were not read in by read_mod, assign mu to means  */
  if( muX[0] == 0. && muX[N-1] == 0.  ){ /* this criterion could be improved e.g. use sum over all entries */ 
    for(i=0 ; i<N ; i++ ){
      for (j=0; j<m[i] ; j++ ){
	muX[i] += X[i][j] ;
      }
      muX[i] /= (float)(m[i]);
    }
  }

  /****************************************************************************/
  /*********       Memory allocation, Initialization                     ******/
  /****************************************************************************/

  /* beta vector and gradient indexed 1,..,6 for NumRec routines*/
  p = (double *) calloc(NDIM+1, sizeof(double));
  del  = (double *) calloc(NDIM+1 ,sizeof(double));
  p_old = (double *) calloc(NDIM+1, sizeof(double));

  /* parameters of covariance matrix and vectors for single gene objective function */

  grad_mu = (double *) calloc(NDIM_MU+1, sizeof(double));  

  for( i=1 ; i<= NDIM ; i++ ){
    p[i] = cov[i];
  }

  if ( evolution == TRUE ){ 
    /* print initial parameter evolution file */
    fprintf ( fp_evolution, "%4d ", 0);
    for ( j=1 ; j<= NDIM; j++ ){ 
      fprintf ( fp_evolution," %9.4f ", p[j]);
      fflush( fp_evolution );
    }
  } 

  /*****************************************************************************/
  /** Master loop begins here                                             ******/
  /*****************************************************************************/

  if ( disp_iter == FALSE ) printf("Iterations in progress: "); 
  fflush(stdout);

  for (k=0 ; k<MAX_MASTER ; k++){ /* master loop */

    /* display . for iteration */
    if ( disp_iter == FALSE ){
      if ( (k+1) % 5  ){
	printf(".");
	fflush(stdout);
      } else {
	printf("%d" , k+1 );
	fflush(stdout);
      }
    }

    if ( disp_iter == TRUE )
      printf("\n\n ********** Starting master loop iteration %d ************** \n", k);

    /****************************************************************************/
    /*********        Optimize beta                                       ******/
    /****************************************************************************/
    
    for( i=0 ; i<N ; i++ ){
      for (j=0; j<m[i] ; j++ ){
	X[i][j] -= muX[i] ;
      }
    }

    /* save current values of beta to later assess progress of iteration */
    for( i=1 ; i<= NDIM ; i++ ){
      p_old[i] = p[i];
    }
    Z_old = Z_OC(p);

    if ( disp_iter == TRUE )
      printf("Objective function at beginning of master it. no. %d is %14f\n", k, Z_old);   

    /* log10(tolerance) drops with each master iteration */
    LINMIN_TOL = pow(10., log10(LINMIN_TOL_HIGH) + k* (log10(LINMIN_TOL_LOW)-log10(LINMIN_TOL_HIGH)) / (double) MAX_MASTER ) ; 
    FTOL       = pow(10., log10(FTOL_HIGH) + k* (log10(FTOL_LOW)-log10(FTOL_HIGH)) / (double) MAX_MASTER );
    
    if ( disp_iter == TRUE )
      printf("\nlog10(LINMIN_TOL): %f, log10(FTOL): %f\n\n", log10(LINMIN_TOL), log10(FTOL)); 

    /* shift intensities by mus for this evaluation below */
    /*
      for( i=0 ; i<N ; i++ ){
      for (j=0; j<m[i] ; j++ ){
      X[i][j] -= muX[i] ;
      }
      }
    */

    if ( disp_iter == TRUE )
      printf("Objective function at beginning of master it. no. %d is %14f\n", k, Z_OC(p));   
    if ( evolution == TRUE ) 
      fprintf( fp_evolution, "%18.3f\n", Z_OC(p) ); 

    if ( disp_iter == TRUE )
      printf("\n+++Starting epsilon+delta optimization++++++\n");
    LINMIN_TOL = LINMIN_TOL_LOW;
    frprmn(p, NDIM, FTOL_LOW, &iter, &fret, Z_OC, gradZ_OC); 
    
    if (verbose == TRUE){
      printf("Iterations: %3d\n",iter);    
      printf("frprmn solution vector: (");
      for (j=1 ;j<=NDIM ; j++) printf("%6.4f,", p[j]);
      printf(")\n");
      printf("Func. value at solution %14f\n",fret);      
    }

    /* print result of beta optimization */
    
    if ( disp_iter == TRUE ){
      printf( " Beta: ");
      for (j=1 ;j<=NDIM ; j++) printf("%6.4f,", p[j]);
      printf("\n");
      printf( " Value of objective function: %14f\n", fret);
    }

    /* copy covariance into global location */
    for (i=1 ; i<=NDIM ; i++ ) cov[i] = p[i];   

    /****************************************************************************/
    /*********        Optimize mus                                         ******/
    /****************************************************************************/

    if ( disp_iter == TRUE )    
      printf("\n+++Starting mu optimization++++++\n");

   /* Need unshifted intensities below */
    for( i=0 ; i<N ; i++ ){
      for (j=0; j<m[i] ; j++ ){
	X[i][j] += muX[i] ;
      }
    }

    for(i=0 ; i<N; i++ ){ /* loop over genes */

      if ( disp_iter == TRUE )
	printf("mu optimization for gene %d begins \n", internal2full[i]); 
      mu_x = (double) muX[i] ;
      gene_index = i;                    
      
      /* seed unconstrained with result from above */
      mu_scalar  = mu_x;


      /********* Will have to change to single variable minimization ****/

      if ( disp_iter == TRUE )
	printf("mu_scalar before minimization %f\n", mu_scalar ); 
      
      /*frprmn(mu_scalar, NDIM_MU, FTOL, &iter, &fret, Z_gene_OC, derivativeZ_gene_OC); */

      ax = mu_x; 
      xx = mu_x+ 5.; /* who knows? */ 
      mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,Z_gene_OC);
      brent(ax,xx,bx,Z_gene_OC,LINMIN_TOL,&xmin);
      mu_scalar=xmin;

      if ( disp_iter == TRUE ){
	printf("mu_scalar after minimization  %f\n", mu_scalar );
	printf("mu derivative: %f\n", derivativeZ_gene_OC(mu_scalar) ); 
      }

#ifdef MUNONNEG
      mu_scalar = fabs( mu_scalar );
#endif

      /* replace revised muX */
      muX[i] = mu_scalar;
    } /* mu optimization ends here. X is not mu subtracted at this point*/ 


    /* correct for instances of finding negative sigma */ 
    if ( p[1] < 0 ){ 
      p[1] = -p[1];
    } 
    if ( p[2] < 0 ){
      p[2] = -p[2]; 
    } 

    /* summarize  result of beta optimization */
    if ( disp_iter == TRUE ){
      printf( " Beta after master loop iteration %d :  ",k);
      for (j=1 ;j<=NDIM ; j++) printf(" %6.4f ", p[j]);
      printf( "\n");
    }

    if ( evolution == TRUE ){ 
      /* print results of master loop iteration in parameter evolution file */
      fprintf ( fp_evolution, "%4d ", k+1);
      for ( j=1 ; j<= NDIM; j++ ){ 
	fprintf ( fp_evolution," %9.4f ", p[j]);
	fflush( fp_evolution );
      }
    }

    /* exit master loop if maximum absolute increment is below threshold */
    abs_dev = 0.;
    for ( j=1 ; j<=NDIM ; j++ ){
      if ( fabs(p[j]-p_old[j]) > abs_dev ) abs_dev = fabs(p[j]-p_old[j]);
    }

    /* shift for evaluation */
    for( i=0 ; i<N ; i++ ){
      for (j=0; j<m[i] ; j++ ){
	X[i][j] -= muX[i] ;
      }
    }

    

    Z_diff = fabs(Z_OC(p)-Z_old);

    if ( disp_iter == TRUE )
      printf("Objective function at end of master it. no. %d is %14f\n", k, Z_OC(p));   


    /* unshift */
    for( i=0 ; i<N ; i++ ){
      for (j=0; j<m[i] ; j++ ){
	X[i][j] += muX[i] ;
      }
    }

    if ( Z_diff > abs_dev ) abs_dev = Z_diff ; 
    if ( abs_dev < EPS_BETA && k > MIN_MASTER )
      break; 

  } /* close master loop */
    
  /* shift intensities by mus for this evaluation below */
  for( i=0 ; i<N ; i++ ){
    for (j=0; j<m[i] ; j++ ){
      X[i][j] -= muX[i] ;
    }
  }
  
  printf("\nOptimization completed in %d iterations.\n", k+1);
  if( k == MAX_MASTER ) printf("Warning: This is the maximum of number of allowed iterations\n");

  printf( "Final objective value %f\n\n", Z_OC(cov)); 
  if ( evolution == TRUE )
    fprintf( fp_evolution, "%18.3f\n", Z_OC(cov));

  /****************************************************************************/
  /*********   Print out final model                                     ******/
  /****************************************************************************/

  fprintf(fp_out_model,"sig_eps_x sig_del_x\n");
  fprintf(fp_out_model," %8.5f %8.5f",
	     p[1],p[2]);
  fprintf(fp_out_model,"\n");

  fclose(fp_in_data); 

  if ( start_mod == TRUE ){
  fclose(fp_in_model); 
  }

  if ( evolution == TRUE ) 
    fclose( fp_evolution );

  stdout_prn_time(); /* display elapsed time */

  return (0); 

}

















