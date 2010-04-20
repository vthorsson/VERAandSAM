#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXSTRING 100

/* structures */
enum BOOL {FALSE, TRUE};
typedef   enum BOOL   BOOL;

typedef struct {
  clock_t begin_clock, save_clock;
  time_t begin_time, save_time;
} time_keeper;

/* global variables */
static time_keeper tk;   /*known only to this file */

/* functions */
void die (int);
void start_time (void);
double prn_time (FILE *);
double stdout_prn_time (void);

#endif
