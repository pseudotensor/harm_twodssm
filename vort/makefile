#
#
#
#CFLAGS = -O3 -ipo -tpp6 -unroll -w1 -axiMKW -I/usr/local/include -L/home/bryanj/fftw/lib
CFLAGS = -O3 -ipo -tpp6 -unroll -w1 -axiMKW -I/usr/local/include 
#CFLAGS = -O0 -g
LDFLAGS = -lm
#
CC = icc

SRCS = \
main.c bounds.c damp_perturb.c dfour1.c diag.c drealft.c dump.c image.c \
init.c nrutil.c ranc.c remap_linear.c remap_flux.c remap_fourier.c restart.c \
set_arrays.c step_ie.c step_pg.c step_trans.c step_visc.c \
stepvar.c sweepx.c sweepy.c sweepym.c timestep.c util.c \
zero_arrays.c fft_grid.c tfft.c fft_mayer.c init_spec.c

OBJS = \
main.o bounds.o damp_perturb.o dfour1.o diag.o drealft.o dump.o image.o \
init.o nrutil.o ranc.o remap_linear.o remap_flux.o remap_fourier.o restart.o \
set_arrays.o step_ie.o step_pg.o step_trans.o step_visc.o \
stepvar.o sweepx.o sweepy.o sweepym.o timestep.o util.o \
zero_arrays.o fft_grid.o tfft.o fft_mayer.o init_spec.o

twodssm: $(OBJS) $(SRCS) makefile
	$(CC) $(CFLAGS) -o twodssm $(OBJS) $(LDFLAGS)

# dependencies
$(OBJS) : defs.h decs.h makefile

SRCS2 = \
test.c set_arrays.c util.c

OBJS2 = \
test.o set_arrays.o util.o

shiftcol: $(OBJS2) $(SRCS2) makefile
	$(CC) $(CFLAGS) -o shiftcol $(OBJS2) $(LDFLAGS)

# dependencies
$(OBJS2) : defs.h decs.h makefile

clean:
	rm *.o *.il 

SRC4 = r8toras.c nrutil.c
OBJ4 = r8toras.o nrutil.o
r8toras: $(OBJ4) $(SRC4) makefile
	$(CC) $(CFLAGS) -o r8toras $(OBJ4) $(LDFLAGS)

