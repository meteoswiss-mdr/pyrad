/*
  *************************************************************

  Functions that can be called from IDL

  *************************************************************

  Filename:         idl.h
  Author:           Andreas Leuenberger
  Creation date:    2012-10-25
  Last update:      2014-03-26

  Copyright:        MeteoSwiss

  Project:          MALSplus
  Target:           SunOS sparc, Gnu/Linux x86_64, x86_32
  Compiler:         GCC

  *************************************************************

  Description:
  ------------

  *************************************************************
  =============================================================
 */

#ifndef __IDL_H
#define __IDL_H

#ifdef __cplusplus
   extern "C" {
#endif

/**
 * idl->C mapping functions:
 * Cal from IDL with call_external()
 */

/*
 * =================================================================
 * -----------------                                ----------------
 * =================================================================
 */

/* Print libDX50.so version */
void idl_printLibVersion(void);

/*
 * =================================================================
 * ----------------- Rainbow file reading functions ----------------
 * =================================================================
 */

/* General file infos: (arguments: filename)  */
IDL_LONG idl_rainbow_getNumberOfSlices(int argc, void *argv[]);
IDL_LONG idl_rainbow_getRangeResolution(int argc, void *argv[]);
IDL_LONG idl_rainbow_getRangeSampling(int argc, void *argv[]);
float idl_rainbow_getAngleStep(int argc, void *argv[]);

/* General file infos: (arguments: filename, ptr to return val)  */
IDL_LONG idl_rainbow_getSliceAngles(int argc, void *argv[]);

/* Sensor info: (arguments: filename) */
float idl_rainbow_getLongitude(int argc, void *argv[]);
float idl_rainbow_getLatitude(int argc, void *argv[]);
float idl_rainbow_getAltitude(int argc, void *argv[]);
float idl_rainbow_getWavelength(int argc, void *argv[]);
float idl_rainbow_getBeamwidth(int argc, void *argv[]);

/* Sensor info: (arguments: filename, ptr to return val) */
IDL_LONG idl_rainbow_getSensorType(int argc, void *argv[]);
IDL_LONG idl_rainbow_getSensorID(int argc, void *argv[]);
IDL_LONG idl_rainbow_getSensorName(int argc, void *argv[]);

/* Slice info: (arguments: filename, slicenumber [, pointer to return value/array]) */
IDL_LONG idl_rainbow_getNumberOfAngles(int argc, void *argv[]);
IDL_LONG idl_rainbow_getScanName(int argc, void *argv[]);
IDL_LONG idl_rainbow_getNumberOfRangeBins(int argc, void *argv[]);
IDL_LONG idl_rainbow_getSliceDate(int argc, void *argv[]);
IDL_LONG idl_rainbow_getSliceTime(int argc, void *argv[]);
IDL_LONG idl_rainbow_getFixedAngle(int argc, void *argv[]);
IDL_LONG idl_rainbow_getAngles(int argc, void *argv[]);
IDL_LONG idl_rainbow_getStartAngles(int argc, void *argv[]);
IDL_LONG idl_rainbow_getStopAngles(int argc, void *argv[]);
IDL_LONG idl_rainbow_getSliceData(int argc, void *argv[]);
IDL_LONG idl_rainbow_getAntSpeed(int argc, void *argv[]);

float idl_rainbow_getNoisePowerZh(int argc, void *argv[]);
float idl_rainbow_getNoisePowerZv(int argc, void *argv[]);
IDL_LONG idl_rainbow_getRadarConstanth(int argc, void *argv[]);
IDL_LONG idl_rainbow_getRadarConstantv(int argc, void *argv[]);
IDL_LONG idl_rainbow_getPulseWidthIndex(int argc, void *argv[]);

/* compression of slice parameters to put into BLOB */
IDL_LONG idl_rainbow_compressData(int argc, void *argv[]); // fvj 15.01.2015
IDL_LONG idl_rainbow_getCompressedDataSize(int argc, void *argv[]);

/*
 * =================================================================
 * ----------------- PSR functions ---------------------------------
 * =================================================================
 */

/* PSR-File header info: */
IDL_LONG idl_psr_getNumberOfItems(int argc, void *argv[]);
IDL_LONG idl_psr_getRangeStart(int argc, void *argv[]);
IDL_LONG idl_psr_getRangeStop(int argc, void *argv[]);
IDL_LONG idl_psr_getRangeStep(int argc, void *argv[]);
IDL_LONG idl_psr_getTwoWayPathAttenuation(int argc, void *argv[]);
IDL_LONG idl_psr_getStartTime(int argc, void *argv[]);
IDL_LONG idl_psr_getStopTime(int argc, void *argv[]);
IDL_LONG idl_psr_getPWIndex(int argc, void *argv[]);
IDL_LONG idl_psr_getNoisePowerH(int argc, void *argv[]);
IDL_LONG idl_psr_getNoisePowerV(int argc, void *argv[]);
IDL_LONG idl_psr_getRadarConstantH(int argc, void *argv[]);
IDL_LONG idl_psr_getRadarConstantV(int argc, void *argv[]);
IDL_LONG idl_psr_getdBadu2dBmOffsetH(int argc, void *argv[]);
IDL_LONG idl_psr_getdBadu2dBmOffsetV(int argc, void *argv[]);
IDL_LONG idl_psr_getMatchedFilterLoss(int argc, void *argv[]);
IDL_LONG idl_psr_getRxLossH(int argc, void *argv[]);
IDL_LONG idl_psr_getRxLossV(int argc, void *argv[]);
IDL_LONG idl_psr_getTxLossH(int argc, void *argv[]);
IDL_LONG idl_psr_getTxLossV(int argc, void *argv[]);
IDL_LONG idl_psr_getRadomeLoss(int argc, void *argv[]);
IDL_LONG idl_psr_getTxPowerSplitFactorH(int argc, void *argv[]);
IDL_LONG idl_psr_getTxPowerSplitFactorV(int argc, void *argv[]);
IDL_LONG idl_psr_getMatchedFilterBandwidth(int argc, void *argv[]);
IDL_LONG idl_psr_getTxNormPower(int argc, void *argv[]);

/* CPI header info of all items */
IDL_LONG idl_psr_getAzimuthStartAngles(int argc, void *argv[]);
IDL_LONG idl_psr_getAzimuthStopAngles(int argc, void *argv[]);
IDL_LONG idl_psr_getElevationStartAngles(int argc, void *argv[]);
IDL_LONG idl_psr_getElevationStopAngles(int argc, void *argv[]);
IDL_LONG idl_psr_getNumberOfTxPulsesArray(int argc, void *argv[]);
IDL_LONG idl_psr_getPrfArray(int argc, void *argv[]);
IDL_LONG idl_psr_getNumberOfRangeGatesArray(int argc, void *argv[]);
IDL_LONG idl_psr_getTxPowerArray(int argc, void *argv[]);
IDL_LONG idl_psr_getNoiseArray(int argc, void *argv[]);
IDL_LONG idl_psr_getValueArrays(int argc, void *argv[]);

/* CPI header info */
IDL_LONG idl_psr_printCpiHeader(int argc, void *argv[]);
IDL_LONG idl_psr_getNumberOfTxPulses(int argc, void *argv[]);
IDL_LONG idl_psr_getNumberOfRangeGates(int argc, void *argv[]);
IDL_LONG idl_psr_getAzimuthStartAngle(int argc, void *argv[]);
IDL_LONG idl_psr_getAzimuthStopAngle(int argc, void *argv[]);
IDL_LONG idl_psr_getElevationStartAngle(int argc, void *argv[]);
IDL_LONG idl_psr_getElevationStopAngle(int argc, void *argv[]);
IDL_LONG idl_psr_getPrf(int argc, void *argv[]);
IDL_LONG idl_psr_getTxPowerArray(int argc, void *argv[]);

/* Item info */
IDL_LONG idl_psr_getItemNoise(int argc, void *argv[]);

/* PRF data */
IDL_LONG idl_psr_getPowerSpectrum(int argc, void *argv[]);


#ifdef __cplusplus
   } /* extern "C" */
#endif

#endif /* __IDL_H */
