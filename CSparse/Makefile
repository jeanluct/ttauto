# Build the (minimal) CSparse library to use cs_dmperm.

CC = cc
CFLAGS = -O3 -ffast-math

AR = ar cr
RANLIB = ranlib

all: libcsparse.a

CS = \
	cs_compress.o \
	cs_cumsum.o \
	cs_entry.o \
	cs_fkeep.o \
	cs_permute.o \
	cs_scc.o \
	cs_dfs.o \
	cs_malloc.o \
	cs_pinv.o \
	cs_transpose.o \
	cs_dmperm.o \
	cs_maxtrans.o \
	cs_randperm.o \
	cs_util.o

$(CS): cs.h Makefile

%.o: %.c cs.h
	$(CC) $(CFLAGS) -c $<

libcsparse.a: $(CS)
	$(AR) libcsparse.a $(CS)
	$(RANLIB) libcsparse.a

clean:
	rm -f *.o

purge: distclean

distclean: clean
	rm -f *.a
