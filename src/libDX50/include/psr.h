/*
  *************************************************************

  Power Spectrum Recording: File header (C language code)

  *************************************************************

  Filename:         psr.h
  Author:           Andreas Leuenberger
  Creation date:    2012-10-25
  Last update:      2013-03-05

  Copyright:        MeteoSwiss

  Project:          MALSplus
  Target:           SunOS sparc, Gnu/Linux x86_64, x86_32
  Compiler:         GCC

  *************************************************************

  Description:
  ------------
  Defines for psr reading.

  *************************************************************
  =============================================================
 */

#ifndef __PSR_H
#define __PSR_H

#include <stdint.h>
#include "psr_header.h"

#ifdef __cplusplus
   extern "C" {
#endif

/* version, ...*/
#define PSR_VERSION_NUM  "2.2"
#define PSR_VERSION_DATE __DATE__
#define PSR_VERSION_TIME __TIME__

#define PSR_VERSION_STR ("PSR-tools: Version " PSR_VERSION_NUM " (" PSR_VERSION_DATE " - " PSR_VERSION_TIME ")")

#define PSR_REDUCED_STR            "PSR REDUCED V"
#define PSR_REDUCED_STR_VERSION   (PSR_REDUCED_STR PSR_VERSION_NUM)

#define FLOATSIZE 4  /* Bytes */
#define COMPLEXSIZE 8  /* Bytes */

#define MAXLINESIZE 1024

/* Defines for item header */
#define ITEM_LENGTH_OFFSET 4
#define ITEM_HEADER_OFFSET 12
#define ITEM_TYPE_ID     0x00000400
#define ITEM_TYPE_ID_RED 0x44444444

struct item_header {
  unsigned long     type;
  unsigned long long length;
};


/* Defines for CPI header */
#define CPI_HEADER_LEN  72  /* [Bytes] */
#define CPI_HEADER_STEP 4   /* [Bytes] */

#define CPI_CPI_SIZE_OFFSET           0
#define CPI_NUM_RANGEGATES_OFFSET     4
#define CPI_NUM_TXPULSES_OFFSET       8
#define CPI_TRACKINGNUM2D_OFFSET     12
#define CPI_SYNCCOUNTER_OFFSET       16
#define CPI_AZIMUTHSTART_OFFSET      20
#define CPI_AZIMUTHSTOP_OFFSET       24
#define CPI_ELEVATIONSTART_OFFSET    28
#define CPI_ELEVATIONSTOP_OFFSET     32
#define CPI_POLARIZATIONINDEX_OFFSET 36
#define CPI_AVGTXPOWER_OFFSET        40
#define CPI_TXCORRMODE_OFFSET        44
#define CPI_PRF_OFFSET               48
#define CPI_PRFMODE_OFFSET           52
#define CPI_PRFBATCHREL_OFFSET       56
#define CPI_QUALITYINDICATION_OFFSET 60
#define CPI_TESTSIGINJMASK_OFFSET    64
#define CPI_DATATYPEMASK_OFFSET      68

struct cpi_header {
  uint32_t cpi_size;
  uint32_t num_rangegates;
  uint32_t num_txpulses;       /* per range gate */
  uint32_t trackingNum2D;
  uint32_t syncCounter;
  float         azimuthStart;       /* [deg] */
  float         azimuthStop;        /* [deg] */
  float         elevationStart;     /* [deg] */
  float         elevationStop;      /* [deg] */
  uint32_t polarizationIndex;  /* 1: horizontal, 2:vertical */
  float         avgTxPower;         /* [kW] */
  uint32_t txCorrMode;         /* 1:no correction, 2:TX phase, 3:TX amplitude, 4:both */
  float         prf;                /* [Hz] */
  uint32_t prfmode;            /* 1:single: 2:dual standard, 3:dual adaptive, 3:dual fixed */
  uint32_t prfBatchRel;        /* 1:single: 2:dual standard, 3:dual adaptive, 3:dual fixed */
  uint32_t qualityIndicator;   /* 0: no error */
  uint32_t testSigInjMask;
  uint32_t dataTypeMask;
};

/* Item table from reduced psr file */
struct psr_item_table_entry {
  long offset;
  long uncompressed_len;
  float noise;
};


#define MAXFILENAME 1024

/**
 * Structure to save psr file info
 */
struct psr_file {
  char filename[MAXFILENAME];
  struct psr_file_header header;
  unsigned char header_read;
  void *filebase; /* start of file including header */
  void *base;     /* start of data (without file header) */
  long mmap_size;
  void *pitem;
  int itemnum;
  long long itemsize;
  unsigned char compressed;
  unsigned long item_type;
  struct psr_item_table_entry *item_table;
};

void psr_printVersion(void);

/**
 * C functions (called from IDL->C mapping functions)
 */
/* PSR-File header info: */
long psr_getHeaderLength(char * fname);
long psr_getNumberOfItems(char *fname);
long psr_getRangeStart(char *fname, float *value);
long psr_getRangeStop(char *fname, float *value);
long psr_getRangeStep(char *fname, float *value);
long psr_getTwoWayPathAttenuation(char *fname, float *value);
long psr_getStartTime(char *fname, char *value, int len);
long psr_getStopTime(char *fname, char *value, int len);
long psr_getRxLossH(char *fname, float *value);
long psr_getRxLossV(char *fname, float *value);
long psr_getTxLossH(char *fname, float *value);
long psr_getTxLossV(char *fname, float *value);
long psr_getRadomeLoss(char *fname, float *value);
long psr_getTxPowerSplitFactorH(char *fname, float *value);
long psr_getTxPowerSplitFactorV(char *fname, float *value);
long psr_getPWIndex(char *fname);
long psr_getNoisePowerH(char *fname, long pw, float *value);
long psr_getNoisePowerV(char *fname, long pw, float *value);
long psr_getRadarConstantH(char *fname, long pw, float *value);
long psr_getRadarConstantV(char *fname, long pw, float *value);
long psr_getdBadu2dBmOffsetH(char *fname, long pw, float *value);
long psr_getdBadu2dBmOffsetV(char *fname, long pw, float *value);
long psr_getMatchedFilterLoss(char *fname, long pw, float *value);
long psr_getMatchedFilterBandwidth(char *fname, long pw, float *value);
long psr_getTxNormPower(char *fname, long pw, float *value);

/* CPI header info of all items */
long psr_getAzimuthStartAngles(char *fname, float *data);
long psr_getAzimuthStopAngles(char *fname, float *data);
long psr_getElevationStartAngles(char *fname, float *data);
long psr_getElevationStopAngles(char *fname, float *data);
long psr_getNumberOfTxPulsesArray(char *fname, long *data);
long psr_getPrfArray(char *fname, float *data);
long psr_getNumberOfRangeGatesArray(char *fname, long *data);
long psr_getTxPowerArray(char *fname, float *data);
long psr_getNoiseArray(char *fname, float *data);
long psr_getValueArrays(char *fname, float *azStart, float *azStop,
                        float *elStart, float *elStop, long *ntxpulses, float *prf,
                        long *nrbins, float *txpower);

/* CPI header info */
long psr_printCpiHeader(char *fname, int item);
long psr_getNumberOfTxPulses(char *fname, int item);
long psr_getNumberOfRangeGates(char *fname, int item);
long psr_getAzimuthStartAngle(char *fname, int item, float *value);
long psr_getAzimuthStopAngle(char *fname, int item, float *value);
long psr_getElevationStartAngle(char *fname, int item, float *value);
long psr_getElevationStopAngle(char *fname, int item, float *value);
long psr_getPrf(char *fname, int item, float *value);
long psr_getTxPower(char *fname, int item, float *value);

/* Item / CPI */
long psr_getItemNoise(char *fname, int item, float *value);

/* PRF data */
long psr_getPowerSpectrum(char *fname, int item, int rangegate, float *complex);

/* Functions for datat decompression */
void psr_initDecompress(void);
int psr_decompress(unsigned char *buf, int cnt, int(*writefunc)(void *buf, int len) );

/* misc functions */
void psr_init(void);
int psr_checkFinished();
void psr_cleanup(void);

int psr_openFile(char *fname);
struct psr_file_header *psr_getHeaderStruct(void);
struct psr_file *psr_getFileStruct(void);

#ifdef __cplusplus
   } /* extern "C" */
#endif

#endif /* __PSR_H */
