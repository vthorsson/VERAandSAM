#include "arraystats.h"
#include "io.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/***************************************************************************************************/
/** read_data.c                                                               **********************/
/** read paired sample data file                                              **********************/
/***************************************************************************************************/

int Nheader; 
/* Number of lines in header */
/* declared here since it will be used by output routine write_llr */

void read_data( FILE *fpin ){

#define MAXGENES 50000  /* maximum no. genes for temporary storage */
#define MAXREPS 500 /* maxiumum no. replicates for temporary storage */
#define Mkeep 1 /* keep only genes with replicates Mkeep or greater */
#define MAXLENGTH 30 /* maximum length of gene name strings */

  int i, j, k, l, Norig, col_count;
  char line[1000000], word[1000], save_line[1000000]; 
  float **Xtemp, **Ytemp;
  char **yorf_temp, **gene_name_temp;
  int mstore[MAXGENES], internal2full_temp[MAXGENES] ;
  int Nkeep=0, mtemp=0;
  int N_col, Xcol[MAXREPS],  Ycol[MAXREPS], Fcol[MAXREPS], index_to_find;
  int file_position, line_length, nXcols, nYcols, nFcols, npairs;
  int n_header_col, n_data_col; 
  char index_to_find_string[MAXREPS], full_string[MAXREPS]; 
  char yorf_word[MAXLENGTH], gene_name_word[MAXLENGTH], *token; 
  
  /* allocate temporary memory */
  Xtemp = (float **) calloc( MAXGENES  ,sizeof(float*));
  Ytemp = (float **) calloc( MAXGENES  ,sizeof(float*));
  yorf_temp = (char **) calloc( MAXGENES  ,sizeof(char*));
  gene_name_temp = (char **) calloc( MAXGENES  ,sizeof(char*));
  for ( i=0 ; i<MAXGENES ; i++ ){
    Xtemp[i] =  (float *) calloc( MAXREPS,sizeof(float));
    Ytemp[i] =  (float *) calloc( MAXREPS,sizeof(float));
    yorf_temp[i] =  (char *) calloc( MAXLENGTH,sizeof(char));   /* is this really needed? */
    gene_name_temp[i] =  (char *) calloc( MAXLENGTH,sizeof(char));  /* ditto */
  }

  /*********************************************************************/
  /** Prescreen input file for possible file format violations        **/
  /*********************************************************************/

  /**  Check for equal number of columns (space and tab separations) ***/

  /* Process first line  */
  fgets( line, sizeof(line), fpin); 
  col_count = 0;
  token = strtok ( line, " \t\n\r|");
  if ( token == NULL ){
    printf("Error: First line has no column structure \n"); exit(1);    
  }
  col_count++; 
  for (i=0; i<10*MAXREPS ; i++){
    token = strtok ( NULL, " \t\n\r|");
    if ( token == NULL ) break ; 
    col_count++;
  }
  if ( i == 10*MAXREPS ){
    printf("Error: First line has too many columns\n"); exit(1);    
  }
  n_header_col=col_count; 

  /* begin processing data records */
  N=0; /* data record counter */
  while ( (fgets( line, sizeof(line), fpin))  != NULL ) {
    N++; 
    if  ( line[0] != '#' && (strspn(line," \t\r\n") != strlen(line)) ){   /* ignore lines beginning with '#' and blank lines */

      col_count = 0;
      token = strtok ( line, " \t\n\r");
      if ( token == NULL ){
	printf("Error:Data record no. %d has no column structure\n",N); exit(1);
      }
      col_count++; 
      for (i=0; 10*MAXREPS ; i++){
	token = strtok ( NULL, " \t\n\r");
	if ( token == NULL ) break ; 
	col_count++;
      }
      if ( i == 10*MAXREPS ){
	printf("Error:Data record no. %d has too many columns\n",N); exit(1);
      }
      n_data_col=col_count; 
      if ( n_data_col != n_header_col ){
	printf("Error:No of columns in data record no. %d differs from header\n",N); exit(1);
      }
    }
  }
  rewind(fpin);

  /****************************  Parse first header line, assumed to have column headings ***************************/

  /* find length of first line, exluding spaces at the end */
  fgets( line, sizeof(line), fpin); 
  line_length = strlen(line);  /* '\n' was included and counts in the length of the line */
  line_length-- ; /* subract out the '\n' contribution to the length */ 
  while ( line[line_length-1] == ' ' || line[line_length-1] == '\r' ) line_length-- ;  
  rewind(fpin);
  
  /***************************** find N column, if it exists ************************/
  col_count = -1; 
  strcpy( word, " ");
  while ( strcmp( word, "N") != 0 ){
    fscanf(fpin, "%s", word); /* get word string  */    
    /* printf("%s \n", word ); */
    file_position = ftell(fpin );
    if (file_position == line_length ) break;
    if ( strcmp(word,"|") != 0 ) col_count++ ; 
  }
  N_col = col_count; 
  if ( file_position == line_length ) N_col = -1 ; /* N_col=-1 indicates that there is no no. of replicates column */
  rewind(fpin);

  /**************************** find X fields ************************************/
  nXcols=0;
  for( index_to_find=0 ; index_to_find<MAXREPS ; index_to_find++ ){

    /* create string to search for, e.g. X0 */
    sprintf( index_to_find_string, "%d", index_to_find); /* copy index to find to a string */
    strcpy( full_string, "X" );
    strcat( full_string, index_to_find_string ); /* full_string is "X0" for index to find 0 etc. */
    /* printf( "full_string %s\n", full_string ); */
    
    rewind(fpin); /* spool to beginning of file */
    file_position = ftell(fpin );
    col_count = -1; 
    while ( strcmp( word, full_string) != 0 && (file_position != line_length) ){  /* walk through fields in header string */
      fscanf(fpin, "%s", word);    
      if ( strcmp(word,"|") != 0 ) col_count++ ; 
      file_position = ftell(fpin );
      /* printf("%s %d\n", word, file_position ); */
    }

    if ( strcmp( word, full_string ) == 0 ){ /* string was found on termination (possibly at end of line) */
      Xcol[index_to_find] = col_count; /* if we found a column successfully copy it into Xcol */
      nXcols++; 
      /* printf("%s is in column %d\n", full_string, Xcol[index_to_find] ); */
    } else { /* we terminated *only* because we reached the end of line -> string was not found */
      break;
    }

  }
  if ( file_position == line_length && index_to_find == 0 ){
    printf("Error: No X0 column in data file\n"); exit(1);
  }
  rewind(fpin);
  file_position = ftell(fpin );

  /*************************** find Y fields **********************************/
  nYcols=0;
  for( index_to_find=0 ; index_to_find<MAXREPS ; index_to_find++ ){

    /* create string to search for, e.g. Y0 */
    sprintf( index_to_find_string, "%d", index_to_find); /* copy index to find to a string */
    strcpy( full_string, "Y" );
    strcat( full_string, index_to_find_string ); /* full_string is "Y0" for index to find 0 etc. */
    
    rewind(fpin); /* spool to beginning of file */
    file_position = ftell(fpin );
    col_count = -1; 

    while ( strcmp( word, full_string) != 0 && (file_position != line_length) ){  /* walk through fields in header string */
      fscanf(fpin, "%s", word);    
      if ( strcmp(word,"|") != 0 ) col_count++ ; 
      file_position = ftell(fpin );
      /* printf("%s %d\n", word, file_position ); */
    }

    if ( strcmp( word, full_string ) == 0 ){ /* string was found on termination (possibly at end of line) */
      Ycol[index_to_find] = col_count; /* if we found a column successfully copy it into Ycol */
      nYcols++; 
      /*printf("%s is in column %d\n", full_string, Ycol[index_to_find] ); */
    } else { /* we terminated *only* because we reached the end of line -> string was not found */
      break;
    }

  }
  if ( file_position == line_length && index_to_find == 0 ){
    printf("Error: No Y0 column in data file\n"); exit(1);
  }
  rewind(fpin);
  file_position = ftell(fpin );
  if( nXcols != nYcols ){
    printf("Error: Numbers of X columns(%d) and Y columns (%d) don't agree\n", nXcols, nYcols );
    exit(1);
  }
  npairs = nXcols; 

  /******************************* Look for flag, F, fields *****************************************************/

  nFcols=0;
  for( index_to_find=0 ; index_to_find<MAXREPS ; index_to_find++ ){

    /* create string to search for, e.g. F0 */
    sprintf( index_to_find_string, "%d", index_to_find); /* copy index to find to a string */
    strcpy( full_string, "F" );
    strcat( full_string, index_to_find_string ); /* full_string is "F0" for index to find 0 etc. */
    
    rewind(fpin); /* spool to beginning of file */
    col_count = -1; 
    file_position = ftell(fpin );

    while ( strcmp( word, full_string) != 0 && (file_position != line_length) ){  /* walk through fields in header  */
      fscanf(fpin, "%s", word);                                                   /* terminate either when string is */
      if ( strcmp(word,"|") != 0 ) col_count++ ;                                  /* found or when line ends */
      file_position = ftell(fpin );
    }

    if ( strcmp( word, full_string ) == 0 ){ /* string was found on termination (possibly at end of line) */
      Fcol[index_to_find] = col_count; /* if we found a column successfully copy it into Fcol */
      nFcols++; 
      /* printf("%s is in column %d\n", full_string, Fcol[index_to_find] ); */
    } else { /* we terminated *only* because we reached the end of line -> string was not found */
      break;
    }
      
  }

  if ( file_position == line_length && index_to_find == 0 ){ /* did not find F0 string */
    nFcols = -1; /* denotes that there are no qualitity flags on the replicates */
  }

  rewind(fpin);
  file_position = ftell(fpin );
  if( (nXcols != nFcols) && (nFcols != -1 )){
    printf("Error: Numbers of data pairs (%d) and F columns (%d) don't agree\n", npairs, nFcols );
    exit(1);
  }

  /*** Print out summary on column detection */
  printf("Maximum %d (X,Y) pairs from column headings: ", npairs);
  for ( i=0 ; i<npairs ; i++ ) printf( "(%d,%d) ", Xcol[i], Ycol[i] ) ;
  printf("\n");
  if ( N_col != -1 ){
    printf("Taking actual number of pairs from N column (column %d)\n", N_col);
  } else { 
    printf("Excluding data pairs with '-' place holders ");
    if ( nFcols != -1 ){
      printf("or with an 'X' in F columns ( ");
      for ( i=0 ; i<npairs ; i++ ) printf( "%d ", Fcol[i] ) ;
      printf(")");
    }
  }
  printf("\n");
 
  /**  End parsing of column header (first line) search for additonal header lines beginning with hash (#) */
  rewind(fpin);
  Nheader=1;
  fgets( line, sizeof(line), fpin); /* move file pointer to end of first line */
  strcpy( line, "#" );
  while ( line[0] == '#' || (strspn(line," \t\r\n")==strlen(line))  ){
    Nheader++;
    fgets( line, sizeof(line), fpin); 
  }  
  Nheader--; /* we overcounted by including a data record */
  printf ( "Header contains %d lines\n", Nheader );

  /* Find maximum possible number of genes, Nfull */

  N=1; /* we have read in the first data line already */
  while ( (fgets( line, sizeof(line), fpin))  != NULL ) {
    if  ( line[0] != '#' && (strspn(line," \t\r\n")!=strlen(line))  ) N++ ; /* ignore lines beginning with '#' and blank lines */
  }
  /*  printf("Maximum %d genes\n", N); */

  Nfull = N;
  full2internal = (int *) calloc(Nfull, sizeof(int));  /* allocate memory for indexing to internal variables */

  rewind(fpin);

  /* Skip header  */
  for ( i=0 ; i<Nheader  ; i++ ){
    fgets( line, sizeof(line), fpin); 
  }

  /**************************** Data reading loop begins here **********************************************/

  /* Read expression levels */

  printf("Reading expression levels\n");  

  i=-1; /* data record counter initialized */
  
  while ( (fgets( line, sizeof(line), fpin))  != NULL ) {

    if  ( line[0] != '#' && (strspn(line," \t\r\n") != strlen(line)) ){   /* ignore lines beginning with '#' and blank lines */

      i++ ; /* data record found, increment data record index */
      strcpy ( save_line, line ); /* use save_line for storage, since strtok will chomp line */

      token = strtok ( line, " \t"); /* ORF field */ 
      strcpy ( yorf_word, token );
      token = strtok ( NULL, " \t"); /* NAME field */
      strcpy ( gene_name_word, token );

      if ( N_col != -1 ){ /* if there is an N no. of samples column, get no. of samples */
	for ( l=2 ; l<= N_col ; l++ ) token = strtok( NULL, " \t");
	mtemp = atoi( token ); 
      }

      if ( N_col == -1 ){ /* there is no number of replicates specified */

	mtemp=0 ;               /* counter for number of replicates initialized to 0 */
      
	for ( k=0 ; k<npairs ; k++ ){
	
	  /* Get F (flag) value. If it contains  an 'X' then skip corresponding X, Y pair */
	  strcpy( line, save_line );
	  token = strtok ( line, " \t"); 
	  for ( l=1 ; l<= Fcol[k] ; l++ ) token = strtok( NULL, " \t");
	  if ( strchr ( token, 'X' ) != NULL ) continue ; 
	  
	  /* Get X value */
	  strcpy( line, save_line );
	  token = strtok ( line, " \t"); 
	  for ( l=1 ; l<= Xcol[k] ; l++ ) token = strtok( NULL, " \t");
	  if ( strcmp( "-", token ) == 0 ){
	    continue ;
	  } else {
	    Xtemp[Nkeep][mtemp]= atof(token) ;
	  }
	
	  /* Get Y value */
	  strcpy( line, save_line );
	  token = strtok ( line, " \t"); 
	  for ( l=1 ; l<= Ycol[k] ; l++ ) token = strtok( NULL, " \t");
	  if ( strcmp( "-", token ) == 0 ){
	    continue ;
	  } else {
	    Ytemp[Nkeep][mtemp] = atof( token ) ;
	  }	  
	  mtemp++;
	}

	/* see if mtemp exceeds required amount, Mkeep */
	if ( mtemp >= Mkeep ){
	  strcpy( yorf_temp[Nkeep], yorf_word ); 
	  strcpy( gene_name_temp[Nkeep], gene_name_word ); 
	  full2internal[i] = Nkeep;
	  internal2full_temp[Nkeep] = i;
	  mstore[Nkeep] = mtemp ; 
	  Nkeep++ ; /* increment kept gene counter */
	} else { /* store excluded genes as -1 in indexing array */
	  full2internal[i] = -1;
	}

      } else if ( mtemp >= Mkeep && N_col != -1 ){ /* if N column is used for true no of samples and 
						     number of samples >= no required */

	mstore[Nkeep] = mtemp;

	for ( j=0 ; j <mtemp ; j++ ){
	  
	 /* Get F (flag) value, and make sure it does not contain 'X' */
	  /* Feb 28, 2001 : Cut out temporarily */ 
	  /*
	    strcpy( line, save_line );
	    token = strtok ( line, " \t"); 
	    for ( l=1 ; l<= Fcol[j] ; l++ ) token = strtok( NULL, " \t");
	    if ( strchr ( token, 'X' ) != NULL ){
	    printf ("Error, unexpected X Flag. Gene (%s,%s), column %d\n", yorf_word, gene_name_word, Fcol[j]);
	  }
	  */
	  
	  /* Get X value */
	  strcpy( line, save_line );
	  token = strtok ( line, " \t"); 
	  for ( l=1 ; l<= Xcol[j] ; l++ ) token = strtok( NULL, " \t");
	  if ( strcmp( "-", token ) == 0 ){
	    printf ("Error, unexpected '-' place holder. Gene (%s,%s), column %d\n", yorf_word, gene_name_word, Xcol[j]);
	  } else {
	    Xtemp[Nkeep][j] = atof(token) ;
	  }
 		
	  /* Get Y value */
	  strcpy( line, save_line );
	  token = strtok ( line, " \t"); 
	  for ( l=1 ; l<= Ycol[j] ; l++ ) token = strtok( NULL, " \t");
	  if ( strcmp( "-", token ) == 0 ){
	    printf ("Error, unexpected '-' place holder. Gene (%s,%s), column %d\n", yorf_word, gene_name_word, Ycol[j]);
	  } else {
	    Ytemp[Nkeep][j] = atof( token ) ;
	  }

	}
	strcpy( yorf_temp[Nkeep], yorf_word ); 
	strcpy( gene_name_temp[Nkeep], gene_name_word ); 
        full2internal[i] = Nkeep;
        internal2full_temp[Nkeep] = i;
        Nkeep++ ; /* increment kept gene counter */

      } else if ( mtemp < Mkeep && N_col != 1 ) { /* store excluded genes as -1 in indexing array */
        full2internal[i] = -1;
	
      }
    }
  }                                        

  printf("Found %d genes\n", Nkeep );

  /*  printf("Keep %d genes, having %d repeats or greater\n", Nkeep, Mkeep ); */

  Norig = N; /* needed for de-allocation */
  N = Nkeep;  /* keep only Nkeep genes */

  /* Memory allocation over genes */
  X = (float **) calloc(N,sizeof(float*));
  Y = (float **) calloc(N,sizeof(float*));
  m = (int *) calloc (N, sizeof(int)) ;
  internal2full = (int *) calloc (N, sizeof(int)) ;  
  yorf = (char **) calloc( N,sizeof(char*));
  gene_name = (char **) calloc( N  ,sizeof(char*));
  for ( i=0 ; i<N ; i++ ){
    yorf[i] =  (char *) calloc( MAXLENGTH,sizeof(char));
    gene_name[i] =  (char *) calloc( MAXLENGTH,sizeof(char));
  }

  /* copy unique identifiers and common gene names  */
  for ( i=0 ; i<N ; i++ ){
    strcpy( yorf[i], yorf_temp[i] );
    strcpy( gene_name[i], gene_name_temp[i] );
  }

  /* copy number of sample for genes from temporary to permanent storage */
  for (i=0 ; i<N ; i++ ) m[i] = mstore[i];
  
  /* copy internal2full indexing variable from temporary storage */
  for (i=0 ; i<N ; i++ ) internal2full[i] = internal2full_temp[i];  

  /* memory allocation, over samples */
  for (i=0 ; i<N ; i++ ){ 
    X[i] = (float *) calloc( m[i] ,sizeof(float));      
    Y[i] = (float *) calloc( m[i] ,sizeof(float));  	 
  }

  /* copy X, Y values from temporay memory */
  for (i=0 ; i<N ; i++ ){ 
    for( j=0 ; j<m[i] ; j++ ){
      X[i][j] = Xtemp[i][j];
      Y[i][j] = Ytemp[i][j];
    }
  }

  /* deallocate temporary memory */
  for (i=0; i<Norig; i++ ){
    free(Xtemp[i]);free(Ytemp[i]);
  }
  free(Xtemp);free(Ytemp);

}

/***************************************************************************************************/
/** write_llr.c                                                               **********************/
/***************************************************************************************************/

void write_llr( FILE *fpin, FILE *fpout, double *llrvals ){

  int i, internal_index;
  char line[10000], return_flag;
  float sigma_delta_x=cov[4], sigma_delta_y=cov[5], t_x, t_y ; 
  float log_ratio_displayed, factor; 
  int tempered_count=0  ; 
  float tempered_log_ratio ( float, float, float, float, char* ); 
  float max_threshold; 

  rewind( fpin );

  /* Copy header to data output file , append new column headings */

  fgets( line, sizeof(line), fpin );
  fputs( line, fpout );
  fseek( fpout, -1, SEEK_CUR); /* go back over the "\n" */
  if( line[strlen(line)-2] == '\r') fseek( fpout, -1, SEEK_CUR); /* go back over the "\r" if found */
  fputs( "      mu_X       mu_Y       lambda   muRATIO  T\n", fpout);  

  for (i=1; i<Nheader; i++ ){
    fgets( line, sizeof(line), fpin );
    fputs( line, fpout );
    if( line[strlen(line)-2] == '\r'){
      fseek( fpout, -2, SEEK_CUR); /* go back over the "\r" if found */
      fputs( "\n", fpout ); /* "replace" the line feed */
    }
  }

  /* copy lines from infile to outfile and append mu_x, mu_y, and llr score */

  i=-1 ; /* initialize full gene counter */
  while ( (fgets( line, sizeof(line), fpin))  != NULL ) {

    fprintf( fpout, "%s", line );  /* output line */
    if  ( line[0] != '#' && (strspn(line," \t\r\n") != strlen(line)) ){  
      i++ ; /* data record found, increment counter */
      fseek( fpout, -1, SEEK_CUR ); /* one character back to eliminate line feed */      
      if( line[strlen(line)-2] == '\r') fseek( fpout, -1, SEEK_CUR); /* go back over the "\r" if found */
      if ( full2internal[i] == -1 ){
	fprintf(fpout," - \n");  /* blank place holders if gene was excluded */
      }
      else {
	internal_index = full2internal[i];
	fprintf(fpout,"%10.2f %10.2f ", muX[internal_index], muY[internal_index] );
	fprintf(fpout,"%12.6f", llrvals[internal_index]);

	/* compute ratio, take conservative ratio for low intentsity signal */

	/* threshold is factor * standard error of mean */
	factor = 1.;
	if ( m[internal_index] >= 2 ){
     	  t_x = factor * sigma_delta_x / sqrt( (float)(m[internal_index] - 1) ) ; 
     	  t_y = factor * sigma_delta_y / sqrt( (float)(m[internal_index] - 1) ) ; 
	} else { /* there is one sample */
     	  t_x = factor * sigma_delta_x ;
     	  t_y = factor * sigma_delta_y ;
	} 

	t_x = sqrt( (t_x * t_x + t_y * t_y)/2. ) ; 
	t_y = t_x ; 

	log_ratio_displayed = tempered_log_ratio ( muX[internal_index], muY[internal_index], t_x, t_y , &return_flag ); 

	if ( return_flag == 'T' ) tempered_count++ ; 

	fprintf(fpout, "%10.4f  %c\n", log_ratio_displayed, return_flag ) ; 
	
      }
    }

  }    

  max_threshold = factor * ( sqrt(sigma_delta_x*sigma_delta_x + sigma_delta_y*sigma_delta_y) ) / sqrt(2.*3.)  ;  
  /* printf("4 sample threshold %f\n",  max_threshold  );  */
  printf("For display in output file, %d log ratios were tempered\n\n", tempered_count); 
 
  fclose(fpin);fclose(fpout);

}

/**************************************************************************************************/
/* tempered_log_ratio                                                       ***********************/
/**************************************************************************************************/


float tempered_log_ratio ( float x, float y, float xT, float yT , char *flag ){

  /* input variable requirement:  (x,y,xT,yT) >= 0 */ 
  float minT, xtemp, ytemp, lr_straight, lr_threshold; 

  x += 1.e-10 ; /* to avoid division by zero */ 
  y += 1.e-10 ; 

  if ( xT > 0. && yT > 0. ){
    minT = fabs( log10 ( xT / yT ) ) ; 
  } else { 
    minT = 0. ; /* in this case, the straight log ratio will always be returned */
  }

  /* determine value of threshold function at x, y  */
  
  xtemp = x ; 
  ytemp = y ; 
  if ( x < xT ) xtemp = xT ; 
  if ( y < yT ) ytemp = yT ; 

  lr_straight = log10 ( x / y )   ; 
  lr_threshold = fabs( log10( xtemp/ytemp ) ) ; 
  if ( lr_threshold < minT ) lr_threshold = minT ; 

  /* printf ( "%f %f %f %f %f %f %f \n \n", minT, xT, yT, x, y, lr_threshold, lr_straight );  */

  /* If straight log ratio is more extreme than threshold function value, replace with threshold function value */
  /* return flag T to indicate that log ratio has been tempered */

  if ( x >= y && lr_straight > lr_threshold ) {
    *flag = 'T';
    return lr_threshold ; 
  } else if ( x >= y && lr_straight <= lr_threshold ) {
    *flag = '-' ; 
    return lr_straight ; 
  } else if ( x < y && lr_straight < -lr_threshold ) {
    *flag = 'T' ; 
    return -lr_threshold ; 
  } else if ( x < y && lr_straight >= -lr_threshold ) {
    *flag = '-' ; 
    return lr_straight ;
  } else {
    printf ( "Error: No conditions were met\n"); 
    exit(1);
  }  

}


/***************************************************************************************************/
/** read_mod.c                                                                **********************/
/** read error model                                                          **********************/
/***************************************************************************************************/


void read_mod( FILE *fp){

  int i; 
  char line[10000];
  char word[100];
  float dummy;
  
  /* skip line */
  fgets( line, sizeof(line), fp ); 
  
  /* read into global covariance */
  for (i=1 ; i<=NDIM ;i++ ){  
    fscanf(fp," %f ", &dummy);
    cov[i]=(double)dummy;
  }

  /* takes care of a rare problem: when model has no variation in */
  /* x direction and single gene values have only variation in this */
  /* direction and are low we end up outside the allowable search region */
  if ( cov[4] == 0 ) cov[4] = 0.0001  ;
  if ( cov[5] == 0 ) cov[5] = 0.0001  ;

  /* skip linefeed and line  */
  fgets( line, sizeof(line), fp ); 

  fscanf(fp, "%s", word );
  if ( strcmp( word, "gene" ) == 0 ){ /* found mus in the model file */
    fgets( line, sizeof(line), fp ); /* ignore remainder of line */
    for (i=0 ; i<N ; i++ ){ /* read in mu values  */
      fscanf(fp,"%f %f %f ", &dummy, &muX[i], &muY[i]);   
/* better would be to have it verify that there are N genes */

    }
  } 

}












