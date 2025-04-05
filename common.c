#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"

/*******************************************************/
double sTime()
{ static struct  timeval  this_tv;
  static struct  timezone dumbTZ;
  double t;

  gettimeofday(&this_tv, &dumbTZ);
  t = this_tv.tv_sec + 0.000001*this_tv.tv_usec;
  return t;
}


/********************************************/
#define MAX_MSG_SIZE 1024
void printTmp(char *fmt, ...)
{ char    tmp[4*MAX_MSG_SIZE];
  static  int lastSize=0;
  va_list ap;
  int     i;

  /* Grab the message that was passed: It is VERY IMPORTANT */
  /* to have the va_end(ap) after!                          */
  va_start(ap, fmt);
  vsnprintf(tmp, 4*MAX_MSG_SIZE, fmt, ap);
  va_end(ap);

  printf("%s", tmp);
  for (i=strlen(tmp); i<=lastSize; i++)
    printf(" ");
  printf("\r");
  lastSize = strlen(tmp);
  fflush(stdout);
}

