/* File : rainbow.h
This is the python interface file for the DX50 c library */

%module librainbow
%include "cpointer.i"
%include "carrays.i"
%include "cmalloc.i"
%include "typemaps.i"
%{
#define SWIG_FILE_WITH_INIT
/* Includes the header in the wrapper code */
#include "./include/rainbow.h"
#include "./include/vars_def.h"
#include "./include/qUncompress.h"
#include "./include/qCompress.h"
%}

%include "numpy.i"

%init %{
import_array();
%}

/* get fixed angle directly as a float */
%apply float *OUTPUT { float *elevation };

/* get antenna speed directly as a float */
%apply float *OUTPUT { float *antspeed };

/* get angles, data and calibration constants as numpy array */
%apply (float* ARGOUT_ARRAY1, int DIM1) { (float *data, long nangles), (float *data, long nvals) (float *radarconstants, int len) }

/* Typemap definitions, to allow SWIG to properly handle 'char**' data types. */
%typemap(in,numinputs=0) char **ARGOUT_STRPTR (char *temp ) {	
	$1 = &temp;
}

%typemap(argout) char **ARGOUT_STRPTR {
	int i, len;
	PyObject * output_text=PyUnicode_FromString($1[0]);
	$result = output_text;
}

%apply char **ARGOUT_STRPTR { char** str };

%include "./include/rainbow.h"

%clear(float *data, long nangles);
%clear(float *data, long nvals);

%pointer_class(float, floatp);
%array_class(float,floatArray);
%allocators(char)

