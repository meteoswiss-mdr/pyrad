/*
  *************************************************************

  Power Spectrum Recording
  IDL to C interface

  *************************************************************

  Filename:         idl_psr.c
  Author:           Andreas Leuenberger
  Creation date:    2012-10-26
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

#include <stdio.h>
#include <string.h>

#include "idl_export.h"
#include "psr.h"

#ifdef __cplusplus
   extern "C" {
#endif

/*
 * ========================================================================================
 * ----------------------- PSR file header info -------------------------------------------
 * ========================================================================================
 */

/**
 * Get number of psr items of a psr file.
 * An item is a cpi interval from a start azimuth/elevation
 * angle to a stop azimuth/elevation angle.
 *
 * @param filename  Name of the psr file
 *
 * @return Number of psr items, or -1 if an error occured.
 */
IDL_LONG idl_psr_getNumberOfItems(int argc, void *argv[])
{
  IDL_STRING filename;

  if (argc != 1) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);

  return (IDL_LONG)psr_getNumberOfItems(filename.s);
}

/**
 * Get range start
 *
 * @param filename  Name of the psr file
 * @param value     Range start [km]
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getRangeStart(int argc, void *argv[])
{
  IDL_STRING filename;
  float *value;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  value = ((float *)argv[1]);
  return psr_getRangeStart(filename.s, value);
}

/**
 * Get range stop
 *
 * @param filename  Name of the psr file
 * @param value     Range stop [km]
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getRangeStop(int argc, void *argv[])
{
  IDL_STRING filename;
  float *value;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  value = ((float *)argv[1]);
  return psr_getRangeStop(filename.s, value);
}

/**
 * Get range step
 *
 * @param filename  Name of the psr file
 * @param value     Range step [km]
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getRangeStep(int argc, void *argv[])
{
  IDL_STRING filename;
  float *value;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  value = ((float *)argv[1]);
  return psr_getRangeStep(filename.s, value);
}

/**
 * Get two-way atmospheric attenuation [db/km]
 *
 * @param filename  Name of the psr file
 * @param value     attenuation
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getTwoWayPathAttenuation(int argc, void *argv[])
{
  IDL_STRING filename;
  float *value;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  value = ((float *)argv[1]);
  return psr_getTwoWayPathAttenuation(filename.s, value);
}

/**
 * Get start time
 *
 * @param filename  Name of the psr file
 * @param value     String to write the start time
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getStartTime(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_STRING value;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  value    = *((IDL_STRING *)argv[1]);
  return psr_getStartTime(filename.s, value.s, value.slen);
}

/**
 * Get stop time
 *
 * @param filename  Name of the psr file
 * @param value     String to write the stop time
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getStopTime(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_STRING value;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  value    = *((IDL_STRING *)argv[1]);
  return psr_getStopTime(filename.s, value.s, value.slen);
}

/**
 * Get receiver loss of horizontal channel [db]
 *
 * @param filename  Name of the psr file
 * @param value     attenuation
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getRxLossH(int argc, void *argv[])
{
  IDL_STRING filename;
  float *value;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  value = ((float *)argv[1]);
  return psr_getRxLossH(filename.s, value);
}

/**
 * Get receiver loss of vertical channel [db]
 *
 * @param filename  Name of the psr file
 * @param value     attenuation
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getRxLossV(int argc, void *argv[])
{
  IDL_STRING filename;
  float *value;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  value = ((float *)argv[1]);
  return psr_getRxLossV(filename.s, value);
}

/**
 * Get transmitter loss of horizontal channel [db]
 *
 * @param filename  Name of the psr file
 * @param value     attenuation
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getTxLossH(int argc, void *argv[])
{
  IDL_STRING filename;
  float *value;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  value = ((float *)argv[1]);
  return psr_getTxLossH(filename.s, value);
}

/**
 * Get transmitter loss of vertical channel [db]
 *
 * @param filename  Name of the psr file
 * @param value     attenuation
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getTxLossV(int argc, void *argv[])
{
  IDL_STRING filename;
  float *value;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  value = ((float *)argv[1]);
  return psr_getTxLossV(filename.s, value);
}

/**
 * Get radome loss [db]
 *
 * @param filename  Name of the psr file
 * @param value     attenuation
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getRadomeLoss(int argc, void *argv[])
{
  IDL_STRING filename;
  float *value;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  value = ((float *)argv[1]);
  return psr_getRadomeLoss(filename.s, value);
}

/**
 * Get transmit power split factor of horizontal channel [1]
 *
 * @param filename  Name of the psr file
 * @param value     attenuation
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getTxPowerSplitFactorH(int argc, void *argv[])
{
  IDL_STRING filename;
  float *value;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  value = ((float *)argv[1]);
  return psr_getTxPowerSplitFactorH(filename.s, value);
}

/**
 * Get transmit power split factor of vertical channel [1]
 *
 * @param filename  Name of the psr file
 * @param value     attenuation
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getTxPowerSplitFactorV(int argc, void *argv[])
{
  IDL_STRING filename;
  float *value;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  value = ((float *)argv[1]);
  return psr_getTxPowerSplitFactorV(filename.s, value);
}

/**
 * Get number pulse width index using in a psr file.
 *
 * @param filename  Name of the psr file
 *
 * @return PW index, or -1 if an error occured.
 */
IDL_LONG idl_psr_getPWIndex(int argc, void *argv[])
{
  IDL_STRING filename;

  if (argc != 1) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  return (IDL_LONG)psr_getPWIndex(filename.s);
}

/**
 * Get horizontal noise power
 *
 * @param filename  Name of the psr file
 * @param pw        Pulse width number (number from 0 to 3)
 * @param value     Variable where to write the noise power [lin adu]
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getNoisePowerH(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG pw;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  pw  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);
  return psr_getNoisePowerH(filename.s, pw, value);
}

/**
 * Get vertical noise power
 *
 * @param filename  Name of the psr file
 * @param pw        Pulse width number (number from 0 to 3)
 * @param value     Variable where to write the noise power [lin adu]
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getNoisePowerV(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG pw;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  pw  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);
  return psr_getNoisePowerV(filename.s, pw, value);
}

/**
 * Get horizontal radar constant
 *
 * @param filename  Name of the psr file
 * @param pw        Pulse width number (number from 0 to 3)
 * @param value     Variable where to write the radar constant
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getRadarConstantH(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG pw;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  pw  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);
  return psr_getRadarConstantH(filename.s, pw, value);
}

/**
 * Get vertical radar constant
 *
 * @param filename  Name of the psr file
 * @param pw        Pulse width number (number from 0 to 3)
 * @param value     Variable where to write the radar constant
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getRadarConstantV(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG pw;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  pw  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);
  return psr_getRadarConstantV(filename.s, pw, value);
}

/**
 * Get horizontal dbadu to dBm offset
 *
 * @param filename  Name of the psr file
 * @param pw        Pulse width number (number from 0 to 3)
 * @param value     Variable where to write
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getdBadu2dBmOffsetH(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG pw;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  pw  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);
  return psr_getdBadu2dBmOffsetH(filename.s, pw, value);
}

/**
 * Get vertical dbadu to dBm offset
 *
 * @param filename  Name of the psr file
 * @param pw        Pulse width number (number from 0 to 3)
 * @param value     Variable where to write
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getdBadu2dBmOffsetV(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG pw;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  pw  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);
  return psr_getdBadu2dBmOffsetV(filename.s, pw, value);
}

/**
 * Get matched filter loss [db]
 *
 * @param filename  Name of the psr file
 * @param pw        Pulse width number (number from 0 to 3)
 * @param value     Variable where to write matched filter loss
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getMatchedFilterLoss(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG pw;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  pw  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);
  return psr_getMatchedFilterLoss(filename.s, pw, value);
}

/**
 * Get matched filter bandwidth [kHz]
 *
 * @param filename  Name of the psr file
 * @param pw        Pulse width number (number from 0 to 3)
 * @param value     Variable where to write matched filter loss
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getMatchedFilterBandwidth(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG pw;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  pw  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);
  return psr_getMatchedFilterBandwidth(filename.s, pw, value);
}

/**
 * Get tx norm power [kW]
 *
 * @param filename  Name of the psr file
 * @param pw        Pulse width number (number from 0 to 3)
 * @param value     Variable where to write matched filter loss
 *
 * @return 0 on success, -1 otherwise
 */
IDL_LONG idl_psr_getTxNormPower(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG pw;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  pw  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);
  return psr_getTxNormPower(filename.s, pw, value);
}

/*
 * ========================================================================================
 * ----------------------- CPI header info of all items -----------------------------------
 * ========================================================================================
 */

/**
 * Get azimuth start angles of all slices.
 *
 * @param filename  Name of the psr file
 *
 * @return 0 on sucess, -1 otherwise.
 */
IDL_LONG idl_psr_getAzimuthStartAngles(int argc, void *argv[])
{
  IDL_STRING filename;
  float *data;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  data     =   (float *)argv[1];

  return (IDL_LONG)psr_getAzimuthStartAngles(filename.s, data);
}

/**
 * Get azimuth stop angles of all slices.
 *
 * @param filename  Name of the psr file
 *
 * @return 0 on sucess, -1 otherwise.
 */
IDL_LONG idl_psr_getAzimuthStopAngles(int argc, void *argv[])
{
  IDL_STRING filename;
  float *data;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  data     =   (float *)argv[1];

  return (IDL_LONG)psr_getAzimuthStopAngles(filename.s, data);
}


/**
 * Get elevation start angles of all slices.
 *
 * @param filename  Name of the psr file
 *
 * @return 0 on sucess, -1 otherwise.
 */
IDL_LONG idl_psr_getElevationStartAngles(int argc, void *argv[])
{
  IDL_STRING filename;
  float *data;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  data     =   (float *)argv[1];

  return (IDL_LONG)psr_getElevationStartAngles(filename.s, data);
}

/**
 * Get elevation stop angles of all slices.
 *
 * @param filename  Name of the psr file
 *
 * @return 0 on sucess, -1 otherwise.
 */
IDL_LONG idl_psr_getElevationStopAngles(int argc, void *argv[])
{
  IDL_STRING filename;
  float *data;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  data     =   (float *)argv[1];

  return (IDL_LONG)psr_getElevationStopAngles(filename.s, data);
}

/**
 * Get the number of tx pulses of all slices.
 *
 * @param filename  Name of the psr file
 *
 * @return 0 on sucess, -1 otherwise.
 */
IDL_LONG idl_psr_getNumberOfTxPulsesArray(int argc, void *argv[])
{
  IDL_STRING filename;
  long *data;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  data     =   (long *)argv[1];

  return (IDL_LONG)psr_getNumberOfTxPulsesArray(filename.s, data);
}

/**
 * Get prf of all slices.
 *
 * @param filename  Name of the psr file
 *
 * @return 0 on sucess, -1 otherwise.
 */
IDL_LONG idl_psr_getPrfArray(int argc, void *argv[])
{
  IDL_STRING filename;
  float *data;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  data     =   (float *)argv[1];

  return (IDL_LONG)psr_getPrfArray(filename.s, data);
}

/**
 * Get the number of range gates of all slices.
 *
 * @param filename  Name of the psr file
 *
 * @return 0 on sucess, -1 otherwise.
 */
IDL_LONG idl_psr_getNumberOfRangeGatesArray(int argc, void *argv[])
{
  IDL_STRING filename;
  long *data;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  data     =   (long *)argv[1];

  return (IDL_LONG)psr_getNumberOfRangeGatesArray(filename.s, data);
}

/**
 * Get averaged tx power of all slices.
 *
 * @param filename  Name of the psr file
 *
 * @return 0 on sucess, -1 otherwise.
 */
IDL_LONG idl_psr_getTxPowerArray(int argc, void *argv[])
{
  IDL_STRING filename;
  float *data;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  data     =   (float *)argv[1];

  return (IDL_LONG)psr_getTxPowerArray(filename.s, data);
}

/**
 * Get noise of all slices.
 *
 * @param filename  Name of the psr file
 *
 * @return 0 on sucess, -1 otherwise.
 */
IDL_LONG idl_psr_getNoiseArray(int argc, void *argv[])
{
  IDL_STRING filename;
  float *data;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  data     =   (float *)argv[1];

  return (IDL_LONG)psr_getNoiseArray(filename.s, data);
}

/**
 * Get prf of all slices.
 *
 * @param filename  Name of the psr file
 *
 * @return 0 on sucess, -1 otherwise.
 */
IDL_LONG idl_psr_getValueArrays(int argc, void *argv[])
{
  IDL_STRING filename;
  float *azStart;
  float *azStop;
  float *elStart;
  float *elStop;
  long *ntxpulses;
  float *prf;
  long *nrbins;
  float *txpower;

  if (argc != 9) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  azStart   = (float *)argv[1];
  azStop    = (float *)argv[2];
  elStart   = (float *)argv[3];
  elStop    = (float *)argv[4];
  ntxpulses =  (long *)argv[5];
  prf       = (float *)argv[6];
  nrbins    =  (long *)argv[7];
  txpower   = (float *)argv[8];

  return (IDL_LONG)psr_getValueArrays(filename.s, azStart, azStop,
                                      elStart, elStop, ntxpulses, prf,
                                      nrbins, txpower);
}

/*
 * ========================================================================================
 * ----------------------- CPI header info ------------ -----------------------------------
 * ========================================================================================
 */

/**
 * Print CPI header info
 *
 * @param filename  Name of the psr file
 * @param itemnum   Number of item (startign with 0)
 *
 * @return 0 on success , -1 on error.
 */
IDL_LONG idl_psr_printCpiHeader(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG itemnum;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  itemnum  = *((IDL_LONG *)argv[1]) & 0xffff;

  return (IDL_LONG)psr_printCpiHeader(filename.s, itemnum);
}

/**
 * Get number of tx pulses of item.
 *
 * @param filename  Name of the psr file
 * @param itemnum   Number of item (startign with 0)
 *
 * @return number of tx pulses, -1 on error.
 */
IDL_LONG idl_psr_getNumberOfTxPulses(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG itemnum;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  itemnum  = *((IDL_LONG *)argv[1]) & 0xffff;

  return (IDL_LONG)psr_getNumberOfTxPulses(filename.s, itemnum);
}

/**
 * Get number of range gates of item.
 *
 * @param filename  Name of the psr file
 * @param itemnum   Number of item (startign with 0)
 *
 * @return number of range gates, -1 on error.
 */
IDL_LONG idl_psr_getNumberOfRangeGates(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG itemnum;

  if (argc != 2) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename = *((IDL_STRING *)argv[0]);
  itemnum  = *((IDL_LONG *)argv[1]) & 0xffff;

  return (IDL_LONG)psr_getNumberOfRangeGates(filename.s, itemnum);
}

/**
 * Get azimuth start angle of an item.
 *
 * @param filename  Name of the psr file
 * @param itemnum   Number of item (startign with 0)
 *
 * @return 0 on success or -1 on error.
 */
IDL_LONG idl_psr_getAzimuthStartAngle(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG itemnum;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  itemnum  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);

  return psr_getAzimuthStartAngle(filename.s, itemnum, value);
}

/**
 * Get azimuth stop angle of an item.
 *
 * @param filename  Name of the psr file
 * @param itemnum   Number of item (startign with 0)
 *
 * @return 0 on success or -1 on error.
 */
IDL_LONG idl_psr_getAzimuthStopAngle(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG itemnum;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  itemnum  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);
  return psr_getAzimuthStopAngle(filename.s, itemnum, value);
}

/**
 * Get elevation start angle of an item.
 *
 * @param filename  Name of the psr file
 * @param itemnum   Number of item (startign with 0)
 *
 * @return 0 on success or -1 on error.
 */
IDL_LONG idl_psr_getElevationStartAngle(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG itemnum;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  itemnum  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);

  return psr_getElevationStartAngle(filename.s, itemnum, value);
}

/**
 * Get elevation stop angle of an item.
 *
 * @param filename  Name of the psr file
 * @param itemnum   Number of item (startign with 0)
 *
 * @return 0 on success or -1 on error.
 */
IDL_LONG idl_psr_getElevationStopAngle(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG itemnum;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  itemnum  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);

  return psr_getElevationStopAngle(filename.s, itemnum, value);
}

/**
 * Get prf an item.
 *
 * @param filename  Name of the psr file
 * @param itemnum   Number of item (startign with 0)
 *
 * @return 0 on success or -1 on error.
 */
IDL_LONG idl_psr_getPrf(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG itemnum;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  itemnum  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);

  return psr_getPrf(filename.s, itemnum, value);
}

/**
 * Get averaged tx power of an item.
 *
 * @param filename  Name of the psr file
 * @param itemnum   Number of item (startign with 0)
 *
 * @return 0 on success or -1 on error.
 */
IDL_LONG idl_psr_getTxPower(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG itemnum;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  itemnum  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);

  return psr_getTxPower(filename.s, itemnum, value);
}

/**
 * Get mean noise of an item.
 *
 * @param filename  Name of the psr file
 * @param itemnum   Number of item (startign with 0)
 *
 * @return 0 on success or -1 on error.
 */
IDL_LONG idl_psr_getItemNoise(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG itemnum;
  float *value;

  if (argc != 3) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }
  filename = *((IDL_STRING *)argv[0]);
  itemnum  = *((IDL_LONG *)argv[1]) & 0xffff;
  value = ((float *)argv[2]);

  return psr_getItemNoise(filename.s, itemnum, value);
}



/*
 * ========================================================================================
 * ----------------------- PSR data -------------------------------------------------------
 * ========================================================================================
 */

/**
 * Get power spectrum
 *
 * @param filename  Name of the psr file
 * @param itemnum   Item number (startign with 0)
 * @param rangegage Range gate number
 * @param complex   Pointer to power spectrum to save
 *
 * @return number of range gates, -1 on error.
 */
IDL_LONG idl_psr_getPowerSpectrum(int argc, void *argv[])
{
  IDL_STRING filename;
  IDL_LONG itemnum;
  IDL_LONG rangegate;
  float *complex;

  if (argc != 4) {
    printf("%s:%s: ERROR: Not enough arguments\r\n", __FILE__, __func__);
    return -1;
  }

  filename  = *((IDL_STRING *)argv[0]);
  itemnum   = *((IDL_LONG *)argv[1]) & 0xffff;
  rangegate = *((IDL_LONG *)argv[2]) & 0xffff;
  complex   = ((float *)argv[3]);

  return (IDL_LONG)psr_getPowerSpectrum(filename.s, itemnum, rangegate, complex);
}


#ifdef __cplusplus
   } /* extern "C" */
#endif
