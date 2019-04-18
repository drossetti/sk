#############################################################################
#
# make ARCH=XXXX [DEBUG=1 opzionale]
#	dove XXXX puo' assumere i simboli
#	
#	x86_linux
#	sparc
#	alpha_osf3
#	alpha_osf4
#	alpha_linux
#	ppc_linux
#
#############################################################################

SHELL=/bin/sh
EXE=spin
#  measure.c
SRCS=spin.c montecarlo.c init.c io.c phys.c skrand.c

SRCS2=correlation_fun.c skrand.c

#############################################################################

.PHONY: all

CPPFLAGS= -D"LOGN=11" \
	  -D"KAPPA=1000"
#-DPRINTA  -DDEBUG_LINKS 

LIBS=-lm 

ifndef ARCH
  ARCH=x86_linux
endif

ifeq "$(ARCH)" "sparc"

  CC=gcc
  CXX=g++
  CCC=g++

#  CPPFLAGS+=-DSYSV -D_POSIX_C_SOURCE -DARCH_TYPE=long 

  CPPFLAGS+=-DARCH_TYPE=long \
	-DARCH_UTYPE="unsigned long"
  LIBS += -ldl

  ifdef DEBUG
    CPPFLAGS+=-DDEBUG
    CFLAGS=-Wall -g -fno-inline -msupersparc
    CXXFLAGS=-Wall -g -fno-inline -msupersparc
    LDFLAGS+=-g
    OBJDIR=debug
  else
    CXXFLAGS=-Wall -O4 -msupersparc
    CFLAGS=-Wall -O4 -msupersparc
    LDFLAGS=
    OBJDIR=release
  endif

endif # sparc

ifeq "$(ARCH)" "ppc_linux"

  CC=gcc
  CXX=g++
  CCC=g++

  CPPFLAGS+=-DARCH_TYPE=long \
	-DARCH_UTYPE="unsigned long"

  LIBS = -lmoto -lm

  ifdef DEBUG
    CPPFLAGS+=-DDEBUG
    CFLAGS=-Wall -g
    CXXFLAGS=-Wall -g
    LDFLAGS+=-g
    OBJDIR=debug
  else
    CXXFLAGS=-Wall -O2
    CFLAGS=-Wall -O2
    LDFLAGS+=
    OBJDIR=release
  endif

endif # ppc_linux

ifeq "$(ARCH)" "alpha_linux"

  CC=gcc
  CXX=g++
  CCC=g++

  CPPFLAGS+=-DARCH_TYPE="long" \
	-DARCH_UTYPE="unsigned long"

# -O4 -tune-ev2

  ifdef DEBUG
    CPPFLAGS+=-DDEBUG
    CFLAGS=-Wall -g -fno-inline
    CXXFLAGS=-Wall -g -fno-inline
    LDFLAGS+=-g
    OBJDIR=debug
  else
    CXXFLAGS=-Wall -O2
    CFLAGS=-Wall -O2
    LDFLAGS=
    OBJDIR=release
  endif

endif # alpha

ifeq "$(ARCH)" "alpha_osf3"

  CC=cc
  CXX=test
  CCC=test

  CPPFLAGS+=-DARCH_TYPE="long" \
	-DARCH_UTYPE="unsigned long"

  ifdef DEBUG
    CPPFLAGS+=-DDEBUG
    CFLAGS= -g
    CXXFLAGS= -g
    LDFLAGS+=-g
    OBJDIR=debug
  else
    CFLAGS= -O4 -tune host
    LDFLAGS=
    OBJDIR=release
  endif

endif # alphaosf

ifeq "$(ARCH)" "alpha_osf4"

  CC=cc #-oldc
  CXX=test
  CCC=test

  CPPFLAGS+=-DARCH_TYPE="long" \
	-DARCH_UTYPE="unsigned long"

  ifdef DEBUG
    CPPFLAGS+=-DDEBUG
    CFLAGS= -g
    CXXFLAGS= -g
    LDFLAGS+=-g
    OBJDIR=debug
  else
    CFLAGS= -fast -O4 -tune host
    LDFLAGS=
    OBJDIR=release
  endif

endif # alphaosf4

ifeq "$(ARCH)" "x86_linux"

  CC=gcc
  CXX=g++
  CCC=g++

  CPPFLAGS+= -DARCH_TYPE=long -DARCH_UTYPE="unsigned long"
  CFLAGS=-Wall
  CXXFLAGS=-Wall

  ifdef DEBUG
    CPPFLAGS+=-DDEBUG
    CFLAGS  +=-g -fno-inline
    CXXFLAGS+=-g -fno-inline
    LDFLAGS  =-g
    OBJDIR   =debug
  else
    CPPFLAGS+=-DUSE_ASM_MATH
    CXXFLAGS =-O4 -march=pentiumpro -mcpu=pentiumpro
    CFLAGS   =-O4 -march=pentiumpro -mcpu=pentiumpro
#    CFLAGS   =-O4 -mpentium
    LDFLAGS  =
    OBJDIR   =release
  endif

endif # x86_linux

#############################################################################

DEPEND=$(SRCS:%.c=.%.d) $(SRCS2:%.c=.%.d)
OBJS  = ${SRCS:%.c=$(OBJDIR)/%.o}
OBJS2 = ${SRCS2:%.c=$(OBJDIR)/%.o}


.%.d: %.c
	@$(SHELL) -ec '$(CC) -M $(CPPFLAGS) $< \
		| sed '\''s/$*.o[ :]*/& $@ /g'\'' \
		| sed '\''s/$*.o/release\/$*.o debug\/$*.o/'\'' > $@'

.%.d: %.cc
	@$(SHELL) -ec '$(CCC) -M $(CPPFLAGS) $< \
		| sed '\''s/$*.o[ :]*/& $@ /g'\'' \
		| sed '\''s/$*.o/release\/$*.o debug\/$*.o/'\'' > $@'

$(OBJDIR)/%.o: %.c
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(COMPILE.c) -o $@ $<

$(OBJDIR)/%.o: %.cc
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(COMPILE.cc) -o $@ $<

#############################################################################

all: $(EXE) 

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

topo: topo.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

correlation_fun: $(OBJS2)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)


include $(DEPEND)

#############################################################################

clean:
	-rm -f .*.d $(OBJS) $(EXE)


