#
#	Makefile for VERA and SAM
#
#

 CC     = gcc
# CC     = cc
# CC = /opt/sfw/bin/gcc # for sun01
LIBS   = -lm 
#  CFLAGS = -g
# Warning: optimization will cause conflict with errno.h on many compilers 

LOCAL_HOME =  .
UTILS = $(LOCAL_HOME)
SRC = $(LOCAL_HOME)

INCLUDES = -I$(UTILS)

OBJ_VERA   = $(SRC)/VERA.c
OBJ_SAM      = $(SRC)/SAM.c

UTILS_DOT_H = $(UTILS)/conjugate_gradient.h \
	$(UTILS)/nrutil.h \
	$(UTILS)/arraystats.h \
	$(UTILS)/objective.h \
	$(UTILS)/io.h \
	$(UTILS)/util.h


UTILS_OBJ = $(UTILS)/conjugate_gradient.o \
	$(UTILS)/nrutil.o \
	$(UTILS)/objective.o \
	$(UTILS)/io.o \
	$(UTILS)/util.o 

all: VERA SAM
.SUFFIXES : .o .c
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) $*.c -c -o $*.o

VERA: $(OBJ_VERA) $(UTILS_OBJ) $(UTILS_DOT_H)
	$(CC) $(CFLAGS) $(INCLUDES) $(UTILS_OBJ) $(OBJ_VERA) \
	-o VERA $(LIBS)

SAM: $(OBJ_SAM) $(UTILS_OBJ) $(UTILS_DOT_H)
	$(CC) $(CFLAGS) $(INCLUDES) $(UTILS_OBJ) $(OBJ_SAM) \
	-o SAM $(LIBS)

clean_all:
	rm $(UTILS_OBJ)

