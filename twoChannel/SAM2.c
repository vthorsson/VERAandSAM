/********************************************************************************************/
/**                                                                                        **/
/**  SAM2.c - Significance of differential gene expression for two samples related         **/
/**       by a common reference                                                            **/
/**                                                                                        **/    
/**                 evaluates lambda =                                                     **/
/**      -2*Log( (max unconstrained likelihood)/(max constrained likelihood) )             **/
/**         constrained: mu=mu_x1=mu_x2   unconstrained mu_x1 != mu_x2                     **/
/**         mu_y's constrained to be equal, both in numerator and denominator              **/
/**                                                                                        **/
/** Implementation: Vesteinn Thorsson thorsson@systemsbiology.org                          **/
/**  Version: Nov2, 2001                                                                   **/
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
#include "string.h"
#define NRANSI

#define FTOL_HIGH 1.e-8  /* high fractional functional tolerance change for frprmn */
#define FTOL_LOW 1.e-12  /* low fractional functional tolerance change for frprmn */
#define LINMIN_TOL_HIGH  1.e-9 /* high tolerance for line minimization */
#define LINMIN_TOL_LOW  1.e-10 /* low tolerance for line minimization */

/* #define Mkeep 1 keep only genes with samples Mkeep or greater */
#define MAXLENGTH 15

main (int argc,char *argv[] ){

  char res  ;
  int i, j, k, n, iter, num_sig=0, n1, n2;
  int N1, N2; /* number of genes in each data file */
  int index_2 = 0 ; /* index counter for set 2 */
  int **partner, **partner_temp, n_partner; 
  char **yorf_1, **yorf_2, **gene_name_1, **gene_name_2; /* unique ids and common gene names */ 
  double result, fret;
  double *p, *del; 
  char *infile_1, *modfile_1, *infile_2, *modfile_2, *outfile;
  FILE *fpin_1, *fpmod_1, *fpin_2, *fpmod_2, *fpout;
  double *mu_vec, *grad_mu;
  double mu_x, mu_y;
  double frac_change_1, frac_change_2 ; 
  double ax, bx, xx, fa, fx, fb, xmin, xmin2;
  double numerator, denominator, llr, *loglr, numerator2;
  double mu_x_from_avg, mu_y_from_avg, mu_x_from_cns, mu_y_from_cns;
  double denominator_from_avg, denominator_from_cns;
  float avg_x1, avg_x2, avg_y1, avg_y2;
  float tempered_log_ratio ( float, float, float, float, char* ); /* for some reason, this is not in io.h */
  float duh1, duh2;
  
  /* definitions for tempering log ratio for output */
  char return_flag ; 
  float factor;
  float t_x1, t_x2, sigma_delta_x1, sigma_delta_x2;
  int tempered_count = 0 ;
  float log_ratio_displayed ; 
  /* end definitions for tempering log ratio */

  /* output version number and compilation parameters  */
  printf("SAM2 version: %s\n", __DATE__);
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
  infile_1 = NULL; 
  modfile_1 = NULL;
  infile_2 = NULL; 
  modfile_2 = NULL;
  outfile = NULL;
  /* start_time(); start timer */

  /* parse command line, check to see if verbose and/or display iterations flag is on  */
  for (i=1; i<argc; i++) {
    if (strcmp(argv[i], "-verbose") == 0) verbose = TRUE ;
    else if (strcmp(argv[i], "-iter") == 0) disp_iter = TRUE;
    else if (infile_1 == NULL) infile_1 = argv[i];
    else if (modfile_1 == NULL) modfile_1 = argv[i];
    else if (infile_2 == NULL) infile_2 = argv[i];
    else if (modfile_2 == NULL) modfile_2 = argv[i];
    else if (outfile == NULL) outfile = argv[i];
    else die (6) ;
  }
  if (infile_1 == NULL || modfile_1 == NULL ||infile_2 == NULL || modfile_2 == NULL || outfile == NULL ) die(6); 

  fpin_1 = fopen(infile_1, "r");
  if (fpin_1 == NULL) die (2);
  fpmod_1 = fopen(modfile_1, "r");
  if (fpmod_1 == NULL) die (2); 
  fpin_2 = fopen(infile_2, "r");
  if (fpin_2 == NULL) die (2);
  fpmod_2 = fopen(modfile_2, "r");
  if (fpmod_2 == NULL) die (2); 
  fpout = fopen(outfile, "w");
  if (fpout == NULL) die (2);

  /******************************************************************************************/
  /********   Read model files                                                        *******/
  /******************************************************************************************/                   
  cov = (double *) calloc(NDIM+1, sizeof(double));
  cov1= (double *) calloc(NDIM+1, sizeof(double));  
  cov2= (double *) calloc(NDIM+1, sizeof(double));    

  printf ("Reading model file %s\n",  modfile_1 ); 
  read_mod( fpmod_1 );
  for (i=1 ; i<= NDIM ; i++ ) cov1[i] = cov[i] ; 

  printf ("Reading model file %s\n",  modfile_2 ); 
  read_mod( fpmod_2 );
  for (i=1 ; i<= NDIM ; i++ ) cov2[i] = cov[i] ; 

  /******************************************************************************************/
  /********   Read data files                                                         *******/
  /******************************************************************************************/                   
  /* Read data file 1 */
  printf ("Reading data file %s\n",  infile_1 ); 
  read_data( fpin_1 ); 

  N1 = N;
  m1 = (int *) calloc( N, sizeof(int) ); 
  yorf_1 = (char **) calloc( N ,sizeof(char*)); 
  gene_name_1 = (char **) calloc( N ,sizeof(char*)); 
  X1 = (float **) calloc( N  ,sizeof(float*));          
  Y1 = (float **) calloc( N  ,sizeof(float*));          
  for (i=0 ; i<N ; i++ ){
    m1[i] = m[i] ; 
    X1[i] =  (float *) calloc( m1[i],sizeof(float));
    Y1[i] =  (float *) calloc( m1[i],sizeof(float));
    yorf_1[i] =  (char *) calloc( MAXLENGTH,sizeof(char));
    gene_name_1[i] =  (char *) calloc( MAXLENGTH,sizeof(char));          
    strcpy( yorf_1[i], yorf[i] );
    strcpy( gene_name_1[i], gene_name[i] );
    for (j=0 ; j<m1[i] ; j++ ){
      X1[i][j] = X[i][j] ; 
      Y1[i][j] = Y[i][j] ; 
    }
  }    
  
  /* free memory */
  for (i=1 ; i<N ; i++ ){
    free(X[i]), free(Y[i]);
    free(yorf[i]), free(gene_name[i]);
  } 
  free(m) ; free(X); free(Y) ; free(yorf) ; free(gene_name);
  
  /* Read data file 2 */
  printf ("Reading data file %s\n",  infile_2 ); 
  read_data( fpin_2 ); 

  N2 = N;
  m2 = (int *) calloc( N, sizeof(int) ); 
  yorf_2 = (char **) calloc( N ,sizeof(char*)); 
  gene_name_2 = (char **) calloc( N ,sizeof(char*)); 
  X2 = (float **) calloc( N  ,sizeof(float*));          
  Y2 = (float **) calloc( N  ,sizeof(float*));          
  for (i=0 ; i<N ; i++ ){
    m2[i] = m[i] ; 
    X2[i] =  (float *) calloc( m2[i],sizeof(float));
    Y2[i] =  (float *) calloc( m2[i],sizeof(float));
    yorf_2[i] =  (char *) calloc( MAXLENGTH,sizeof(char));
    gene_name_2[i] =  (char *) calloc( MAXLENGTH,sizeof(char));          
    strcpy( yorf_2[i], yorf[i] );
    strcpy( gene_name_2[i], gene_name[i] );
    for (j=0 ; j<m2[i] ; j++ ){
      X2[i][j] = X[i][j] ; 
      Y2[i][j] = Y[i][j] ; 
    }
  }    
  
  /* free memory */
  for (i=1 ; i<N ; i++ ){
    free(X[i]), free(Y[i]);
    free(yorf[i]), free(gene_name[i]);
  } 
  free(m) ; free(X); free(Y) ; free(yorf) ; free(gene_name);

  m = (int *) calloc (1, sizeof(int)); /* m[0] will be used in Z_gene_12 */

  /******************************************************************************************/
  /********   Locate pairs                                                            *******/
  /******************************************************************************************/                   

  printf("Locating pairs of genes found both in set 1 and set 2\n");

  partner_temp = (int **) calloc ( N1+N2, sizeof( int *) );
  for ( i=0 ; i<(N1+N2) ; i++ ){
    partner_temp[i] = (int *) calloc ( 2, sizeof(int) ) ; 
  }

  index_2 = 0 ; /* counter for current position in set 2 */
  n_partner = 0 ; /* counter for number of pairs found */
  for ( n=0 ; n < N1 ; n++ ){ /* cycle begins over genes in set 1 */
    for ( i=index_2 ; i<N2 ; i++){ /* look for matching gene name in set 2, beginning from index_2 counter */
      res = strcmp(yorf_1[n], yorf_2[i]); 
      if ( res == 0 ){ /* found a match */
	/* printf( "gene no. %d (set1, %s) is found at no. %d (set2) \n" , n+1, yorf_1[n], index_2+1 ); */
	partner_temp[n_partner][0] = n; 
	partner_temp[n_partner][1] = index_2; 
	n_partner++; 
	index_2++; /* increment before exiting */
	break; 
      }
      else if ( res < 0 ){ /* we went too far (no match found), and need to stop */
	/* printf( " No match found to gene no. %d (set1) %s (looking at %s) \n", n+1, yorf_1[n], yorf_2[i] ); */
	break; /* break without incrementing index_2 */
      }
      index_2++;
    } /* end loop over genes in set 2 */
  } /* end loop over genes in set 2 */


  partner = (int **) calloc ( n_partner, sizeof( int *) );
  for ( i=0 ; i<n_partner ; i++ ){
    partner[i] = (int *) calloc ( 2, sizeof(int) ) ;    
    partner[i][0] = partner_temp[i][0] ; 
    partner[i][1] = partner_temp[i][1] ; 
    /*    printf(" %d %d \n", partner[i][0], partner[i][1] ); */
    free( partner_temp[i] );
  }
  free( partner_temp );

  printf("Found %d pairs\n", n_partner);

  /******************************************************************************************/
  /********   Loop over pairs                                                         *******/
  /******************************************************************************************/                   

  mu_vec = (double *) calloc(4, sizeof(double));  /* declare for use with indicess 1,2 or 1,2,3 */
  grad_mu =  (double *) calloc( 4,sizeof(double)); 

  fprintf( fpout, "#GENE_NAME  DESCRIPT.     mu_X1      mu_X2       mu_Y      lambda   muRATIO  T\n"); 
  fprintf( fpout, "#---------  ---------  --------  ---------  ---------   ---------   -------  -\n"); 

  //for ( n=4998 ; n<4999 ; n++ ){
  for ( n=0 ; n<n_partner ; n++ ){

    gene_index_1 = partner[n][0]; /* index of first gene in pair */ 
    gene_index_2 = partner[n][1]; /* index of second gene in pair */

    if( disp_iter == TRUE )       
      printf("++++++gene pair %d ( %s )+++++++\n", n, yorf_1[gene_index_1] ); 

    /* Compute averges for initial guesses for mus  */

    avg_x1=0 ;
    avg_x2=0 ; 
    avg_y1=0 ;
    avg_y2=0 ; 
    for (j=0; j<m1[gene_index_1] ; j++ ){
      avg_x1 += X1[gene_index_1][j] ;
      avg_y1 += Y1[gene_index_1][j] ;
    }
    avg_x1 /= (float)(m1[gene_index_1]) ;
    avg_y1 /= (float)(m1[gene_index_1]) ;
    for (j=0; j<m2[gene_index_2] ; j++ ){
      avg_x2 += X2[gene_index_2][j] ;
      avg_y2 += Y2[gene_index_2][j] ;
    }
    avg_x2 /= (float)(m2[gene_index_2]) ;
    avg_y2 /= (float)(m2[gene_index_2]) ;

    /* Compute numerator, mu_x1 = mu_x2 constrained */
    
    if (disp_iter == TRUE )
      printf("Computing numerator, (mu_x1 = mu_x2 constrained) \n");

    mu_vec[1] = ( avg_x1 + avg_x2 )/2.; 
    mu_vec[2] = ( avg_y1 + avg_y2 )/2.; 

    if (disp_iter == TRUE )
      printf("mu_vec before minimization %f %f\n", mu_vec[1], mu_vec[2]);
    LINMIN_TOL = LINMIN_TOL_LOW;
    frprmn(mu_vec, NDIM_MU, FTOL_LOW, &iter, &numerator, Z_gene_12_xc_yc, gradZ_gene_12_xc_yc); 

#ifdef MUNONNEG
    mu_vec[1] = fabs( mu_vec[1] );
    mu_vec[2] = fabs( mu_vec[2] );

    if ( mu_vec[1] <= 1e-8 || mu_vec[2] <= 1e-8 ){
      LINMIN_TOL = LINMIN_TOL_LOW;
      frprmn(mu_vec, NDIM_MU, FTOL_LOW, &iter, &numerator, Z_gene_12_xc_yc, gradZ_gene_12_xc_yc); 
    }   

    mu_vec[1] = fabs( mu_vec[1] );
    mu_vec[2] = fabs( mu_vec[2] );

#endif

    if (disp_iter == TRUE ){
      printf("mu_vec after minimization  %f %f\n", mu_vec[1], mu_vec[2]); 
      printf("numerator %f\n", numerator); 
    }

    /* Compute denominator, mu_x1 and mu_x2 unconstrained */

    mu_vec[1] = avg_x1; 
    mu_vec[2] = avg_x2; 
    mu_vec[3] = ( avg_y1 + avg_y2 )/2.; 

    if (disp_iter == TRUE )
      printf("Computing denominator, (mu_x1, mu_x2 unconstrained) \n");

    if (disp_iter == TRUE )
      printf("mu_vec before minimization %f %f %f\n", mu_vec[1], mu_vec[2], mu_vec[3]);

    LINMIN_TOL = LINMIN_TOL_LOW;
    frprmn(mu_vec, 3, FTOL_LOW, &iter, &denominator, Z_gene_12_xu_yc, gradZ_gene_12_xu_yc); 

#ifdef MUNONNEG
    mu_vec[1] = fabs( mu_vec[1] );
    mu_vec[2] = fabs( mu_vec[2] );
    mu_vec[3] = fabs( mu_vec[3] );

    
    if ( mu_vec[1] <= 1e-8 || mu_vec[2] <= 1e-8 || mu_vec[3] <= 1e-8 ){
      LINMIN_TOL = LINMIN_TOL_LOW;
      frprmn(mu_vec, 3, FTOL_LOW, &iter, &denominator, Z_gene_12_xu_yc, gradZ_gene_12_xu_yc); 
    }   

    mu_vec[1] = fabs( mu_vec[1] );
    mu_vec[2] = fabs( mu_vec[2] );
    mu_vec[3] = fabs( mu_vec[3] );
#endif

    if (disp_iter == TRUE ){
      printf("mu_vec after minimization  %f %f %f\n", mu_vec[1], mu_vec[2], mu_vec[3]); 
      printf("denominator %f\n", denominator); 
    }

    llr = 2.*(numerator-denominator);

    if ( llr < 0 ) llr = 0 ; /* rare cases. consider setting up a tally for this */


    fprintf( fpout, "%10s %10s", yorf_1[gene_index_1], gene_name_1[gene_index_1]);
    fprintf( fpout,"%10.2f %10.2f %10.2f", mu_vec[1], mu_vec[2], mu_vec[3] );
    fprintf( fpout,"%12.6f", llr );
    
    /* tempered ratio. For first pass, simply replace old expression x->x1 y->y1 */

    /* threshold is factor * standard error of mean */
    sigma_delta_x1 = cov1[4];
    sigma_delta_x2 = cov2[4]; 
    factor = 1.;

    if ( m1[gene_index_1] >= 2 ){
      t_x1 = factor * sigma_delta_x1 / sqrt( m1[gene_index_1] - 1 ) ; 
    } else { /* there is one sample */
      t_x1 = factor * sigma_delta_x1 ;
    } 

    if ( m2[gene_index_2] >= 2 ){
      t_x2 = factor * sigma_delta_x2 / sqrt( m2[gene_index_2] - 1 ) ; 
    } else { /* there is one sample */
      t_x2 = factor * sigma_delta_x2 ;
    } 

    t_x1 = sqrt( (t_x1 * t_x1 + t_x2* t_x2)/2. ) ; 
    t_x2 = t_x1 ; 
  
    duh1 = (float)mu_vec[1]; 
    duh2 = (float)mu_vec[2]; 
    log_ratio_displayed = tempered_log_ratio ( (float)mu_vec[1],(float)mu_vec[2], t_x1, t_x2 , &return_flag ); 
    
    if ( return_flag == 'T' ) tempered_count++ ; 
    
    fprintf(fpout, "%10.4f  %c\n", log_ratio_displayed, return_flag ) ; 

  }

  printf("For display in output file, %d log ratios were tempered\n\n", tempered_count); 
 

}
