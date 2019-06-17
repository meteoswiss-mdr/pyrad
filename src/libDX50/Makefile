# *****************************************************************************
#
#  Projects Makefile
#
#  Filename:         Makefile
#  Author:           Andreas Leuenberger
#  Creation date:    2012-10-26
#  Last update:      2013-02-21
#
#  Copyright:        MeteoSchweiz
#
#  Project:          MALSplus
#
#  Usage:
#  ------
#  make            : Compile libDX50.so and psr tools
#  make lib        : Compile libDX50.so
#  make psr        : psr tools: Compile psr reduce and uncompress applications
#  make clean      : Clean project directory
#
# Options:
# --------
#
# *****************************************************************************

CC = gcc
#LD = gcc

INCLUDEDIRS = -Iinclude
#QTINCLUDEDIR = -I/usr/include/qt4
QTINCLUDEDIR =
CPPINCLUDEDIRS = ${INCLUDEDIRS} ${QTINCLUDEDIR}

OBJDIR = obj

CFLAGS = -Wall -O2 -fPIC -g -m64
#LIBS = -lQtCore
LIBS = -lz -lm

LIBTARGET=libDX50.so
TARGETPATHLIB = lib

LIBTARGET_ZUEUB222=libDX50.zueub222.so
LIBTARGET_LOMWS900=libDX50.lomws900.so
LIBTARGET_DX50=libDX50.dx50.so
LIBTARGET_ZUEUX410=libDX50.zueux410.so
LIBTARGET_ZUEUB242=libDX50.zueub242.so

TARGET_RED=psrReduce
TARGET_UNCOM=psrUncompress
TARGETPATHBIN = bin

INSTALLPATH=${TARGETPATHLIB}

INITFUNC = dx50_init
FINIFUNC = dx50_close

LDFLAGS = -Wall
# Linux
LDLIBFLAGS = -shared -Wl,-soname,${LIBTARGET} -Wl,--init,${INITFUNC} -Wl,--fini,${FINIFUNC}
# SUN:
#LDLIBFLAGS = -shared -Wl,-soname,${LIBTARGET}

CSources_Common = psr.c psr_header.c

LIB_CSources  = ${CSources_Common} dx50.c idl_psr.c idl_rainbow.c rainbow_read_raw.c rainbow_compress_raw.c \
	    	qUncompress.c qCompress.c

PSR_RED_CSources   = ${CSources_Common} psr_reduce.c
PSR_UNCOM_CSources = ${CSources_Common} psr_uncompress.c

LIB_Objects := $(patsubst %.c,$(OBJDIR)/%.o,$(LIB_CSources))
PSR_RED_Objects := $(patsubst %.c,$(OBJDIR)/%.o,$(PSR_RED_CSources))
PSR_UNCOM_Objects := $(patsubst %.c,$(OBJDIR)/%.o,$(PSR_UNCOM_CSources))

.PHONY: all psr lib clean libzueub222 libzueub242 liblomws900 libdx50 force_lib_version \
	force_psr_version

all: lib psr

lib: ${TARGETPATHLIB}/${LIBTARGET}

psr: ${TARGETPATHBIN}/${TARGET_RED} ${TARGETPATHBIN}/${TARGET_UNCOM}

-include $(patsubst %.c,$(OBJDIR)/%.d,$(LIB_CSources)) \
	$(patsubst %.c,$(OBJDIR)/%.d,$(PSR_RED_CSources)) \
	$(patsubst %.c,$(OBJDIR)/%.d,$(PSR_UNCOM_CSources))

# Force compilation to update version string
${OBJDIR}/dx50.o: .FORCE
${OBJDIR}/psr.o:  .FORCE

.FORCE:

${TARGETPATHBIN}/${TARGET_RED}: ${PSR_RED_Objects}
	@echo --- linking $@
	@${CC} ${LDFLAGS} ${PSR_RED_Objects} -o $@ ${LIBS}

${TARGETPATHBIN}/${TARGET_UNCOM}: ${PSR_UNCOM_Objects}
	@echo --- linking $@
	@${CC} ${LDFLAGS} ${PSR_UNCOM_Objects} -o $@ ${LIBS}

${TARGETPATHLIB}/${LIBTARGET}: ${LIB_Objects}
	@echo --- linking $@
	@${CC} ${LDFLAGS} ${LDLIBFLAGS} ${LIB_Objects} -o $@ ${LIBS}

${OBJDIR}/%.o:%.c Makefile
	@echo --- generating $@ from $<
	@${CC} -c $< -o $@ ${CFLAGS} ${INCLUDEDIRS}

${OBJDIR}/%.d:%.c
	@echo --- generating $@ from $<
	@${CC} -M $(CFLAGS) ${INCLUDEDIRS} $< > $@.$$$$;      \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;   \
	rm -f $@.$$$$

libzueub222: ${TARGETPATHLIB}/${LIBTARGET}
	@echo --- installing ${LIBTARGET} to ${INSTALLPATH}/${LIBTARGET_ZUEUB222}
	@cp ${TARGETPATHLIB}/${LIBTARGET} ${INSTALLPATH}/${LIBTARGET_ZUEUB222}

libzueub242: ${TARGETPATHLIB}/${LIBTARGET}
	@echo --- installing ${LIBTARGET} to ${INSTALLPATH}/${LIBTARGET_ZUEUB242}
	@cp ${TARGETPATHLIB}/${LIBTARGET} ${INSTALLPATH}/${LIBTARGET_ZUEUB242}

liblomws900: ${TARGETPATHLIB}/${LIBTARGET}
	@echo --- installing ${LIBTARGET} to ${INSTALLPATH}/${LIBTARGET_LOMWS900}
	@cp ${TARGETPATHLIB}/${LIBTARGET} ${INSTALLPATH}/${LIBTARGET_LOMWS900}

libdx50: ${TARGETPATHLIB}/${LIBTARGET}
	@echo --- installing ${LIBTARGET} to ${INSTALLPATH}/${LIBTARGET_DX50}
	@cp ${TARGETPATHLIB}/${LIBTARGET} ${INSTALLPATH}/${LIBTARGET_DX50}

libzueux410: ${TARGETPATHLIB}/${LIBTARGET}
	@echo --- installing ${LIBTARGET} to ${INSTALLPATH}/${LIBTARGET_ZUEUX410}
	@cp ${TARGETPATHLIB}/${LIBTARGET} ${INSTALLPATH}/${LIBTARGET_ZUEUX410}

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d ${TARGETPATHLIB}/${LIBTARGET} \
	${TARGETPATHBIN}/${TARGET_RED} ${TARGETPATHBIN}/${TARGET_UNCOM}
	@find . -name "*~" | xargs rm -f
