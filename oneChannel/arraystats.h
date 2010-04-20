 
/** user specified compilation parameters  **/
#undef USERHOD /* define (undef), uses (does not use) correlation of delta  */
#undef MUNONNEG /* define (undef), restricts (does not restrict) mu_x, mu_y to be non-negative on mu optimization */
/*** end user specified compilation compilation parameters **/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "util.h" 
#define NRANSI

/** global variables **/
BOOL verbose;     /* flag=1 for verbose output */
BOOL disp_iter;  /* flag=1 to display progress of iterations run */
BOOL start_mod;  /* flag=1  if starting guess to be used for optimization */
BOOL evolution;  /* flag=1 to display parameter evolution during optimization */
BOOL stop_crit; /* flag=1 to specify stop criterion for optimization */
int N;  /* number of genes, maximum no. of repeats */
char **yorf; /* unique gene name identifiers for each gene */
char **gene_name; /* common use gene name */ 
int *m; /* number of samples for each gene */
float **X, **Y, *muX, *muY; /* intensities for all genes&repeats and mu parameters for all genes*/
int *m1, *m2; /* number of samples in each case for reference comparison */
float **X1, **Y1, **X2, **Y2; /* intensities for data sets involved in reference comparison */
float LINMIN_TOL; /* tolerance for linmin.c */
double *cov;   /* parameters of covariance matrix */
double *cov1, *cov2;   /* parameters of covariance matrices for reference sample comparisons */

int gene_index; /* gene index for Z_gene.c and gradZ_gene.c */
int gene_index_1, gene_index_2; /* gene index needed for two sample comparisons */
int Nfull, *full2internal, *internal2full; 
/* total number of genes in original data set Nfull and mapping to N genes used for internal code  */
                          /* full2internal[i] = -1 if gene i was excluded */







#define NDIM_MU 2 /* number of dimensions for single gene mu_x mu_y */



#define NDIM 2

