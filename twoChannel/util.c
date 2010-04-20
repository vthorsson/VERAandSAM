/* utility functions */
#include "util.h"

/* generic error handling function */
void die (int error_code) {
  
  char usage_bmo[600], usage_sig[600];
  char usage_sig2[300]; /* core dump may result if these are too large */

  sprintf(usage_bmo, "%s%s%s%s%s", 
	  "Usage:      VERA <mergedFile> <ErrorModel>\n", 
	  "Options:    -init  <InitialErrorModel>  Use an input model to initiate optimization\n",
	  "            -crit  <number>             Optimization ceases when all changes are less than <number> (Default=0.0005)\n",      
  	  "            -evol                       Output file showing model approach to optimum\n", 
  	  "            -iter                       Display iterations (Use for debugging only)\n"  ); 

  sprintf(usage_sig, "%s%s",
	  "Usage:      SAM <mergedFile> <ErrorModel> <mergedFileSignificance>\n", 
	  "            -iter                      Display iterations (Use for debugging only)\n");

  sprintf(usage_sig2, "%s%s",
	  "Usage:      SAM2 <mergedFile_1> <ErrorModel_1> <mergedFile_2> <ErrorModel_2> <mergedFileSignificance>\n", 
	  "            -iter                      Display iterations (Use for debugging only)\n"); 

  switch (error_code) {

  case 0: 
    fprintf(stderr, "\nMemory allocation failure.\n");
    break;
  case 1: 
    fprintf(stderr, "\nCommand line error.\n\n%s\n\n", usage_bmo);
    break;
  case 2: 
    fprintf(stderr, "\nCannot read input file.\n");
    break;
  case 3:
    fprintf(stderr, "\nNumber out of range.\n");
    break;
  case 5:
    fprintf(stderr, "\nCommand line error.\n\n%s\n\n", usage_sig);
    break;
  case 6:
    fprintf(stderr, "\nCommand line error.\n\n%s\n\n", usage_sig2);
    break;
  default: 
    fprintf(stderr, "\nUnknown error.\n");  
    break;
  }

  fprintf(stderr, "Aborting program execution.\n");
  exit(1);
}

/* timing functions */

void start_time (void) {
  tk.begin_clock = tk.save_clock = clock();
  tk.begin_time = tk.save_time = time(NULL);
}

double prn_time (FILE *fp) {

  char s1[MAXSTRING], s2[MAXSTRING];
  int field_width, n1, n2;
  double clocks_per_second = (double) CLOCKS_PER_SEC, 
    user_time, real_time;

  user_time = ( clock() - tk.save_clock ) / clocks_per_second;
  real_time = difftime( time(NULL), tk.save_time );
  tk.save_clock = clock();
  tk.save_time = time(NULL);

  /* print values found, and do it neatly */
  n1 = sprintf(s1, "%.1f", user_time);
  n2 = sprintf(s2, "%.1f", real_time);
  field_width = (n1 > n2) ? n1 : n2;
  fprintf(fp, "%s%*.1f%s\n%s%*.1f%s\n\n",
	 "user_time ", field_width, user_time, " seconds",
	 "real_time ", field_width, real_time, " seconds");
  return user_time;
}

double stdout_prn_time (void) {

  char s1[MAXSTRING], s2[MAXSTRING];
  int field_width, n1, n2;
  double clocks_per_second = (double) CLOCKS_PER_SEC, 
    user_time, real_time;

  user_time = ( clock() - tk.save_clock ) / clocks_per_second;
  real_time = difftime( time(NULL), tk.save_time );
  tk.save_clock = clock();
  tk.save_time = time(NULL);

  /* print values found, and do it neatly */
  n1 = sprintf(s1, "%.1f", user_time);
  n2 = sprintf(s2, "%.1f", real_time);
  field_width = (n1 > n2) ? n1 : n2;
  printf( "%s%*.1f%s\n%s%*.1f%s\n\n",
	 "user_time ", field_width, user_time, " seconds",
	 "real_time ", field_width, real_time, " seconds");
  return user_time;
}
