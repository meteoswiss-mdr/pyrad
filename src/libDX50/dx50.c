/*
  *************************************************************

  DX50 library main file

  *************************************************************

  Filename:         dx50.c
  Author:           Andreas Leuenberger
  Creation date:    2012-10-29
  Last update:      2013-03-08

  Copyright:        MeteoSwiss

  Project:          MALSplus
  Target:           SunOS sparc, Gnu/Linux x86_64, x86_32
  Compiler:         GCC

  *************************************************************

  Description:
  ------------
  Calls init and fini functions

  History:
  2013-02-21 V0.6  Support for reduced psr files
  2013-03-08 V0.7  Read datatype

  *************************************************************
  =============================================================
 */

#include <stdio.h>

#include "psr.h"

#ifdef __cplusplus
   extern "C" {
#endif

#define LIBNAME "libDX50"
#define VERSION_NUM  "0.9.2"
#define VERSION_DATE __DATE__
#define VERSION_TIME __TIME__

#define VERSION_STR ("Version " VERSION_NUM " (" VERSION_DATE " - " VERSION_TIME ")")

/**
 * Print version
 */
void idl_printLibVersion(void)
{
  printf("%s %s\r\n", LIBNAME, VERSION_STR);
}

/**
 * Library init function
 */
void dx50_init(void)
{
  printf("%% load %s %s\r\n", LIBNAME, VERSION_STR);
  psr_init();
}

/**
 * Library close function
 */
void dx50_close(void)
{
  psr_cleanup();
}

#ifdef __cplusplus
   } /* extern "C" */
#endif
