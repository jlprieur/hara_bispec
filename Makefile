######################################################################
# Makefile for decode...
#
# JLP
# Version 14-04-2008
######################################################################
include $(JLPSRC)/jlp_make.mk
mylib=$(JLPLIB)/jlp
OBJ1=jlp_cover0.o fft_nag.o nagfft.o
OBJ3=jlp_cover_mask.o correct_bisp_auto1.o photon_model.o 
# decode_set.o
OUT= $(OBJ3) inv_bispec2.exe 
#FOURN=fourn.o
FOURN=fft_jlp.o poidev.o
JLIB=$(mylib)/jlpacc.a $(mylib)/jlputil.a \
     $(mylib)/newplot0.a $(mylib)/jlp_splot.a $(mylib)/jlpacc.a
JLIB=$(mylib)/jlpacc.a $(mylib)/jlputil.a $(mylib)/jlp_splot.a

.SUFFIXES: 
.SUFFIXES: .o .c .for .exe $(SUFFIXES)

.for.o:
	runs esoext1 -f $*.for
	$(F77) $(FFLAGS) -c $*.f
	rm $*.f

.for.exe:
	runs esoext1 -f $*.for
	$(F77) $(FFLAGS) -c $*.f
	$(F77) -g -o $(EXEC)/$*.exe $*.o cover_mask.o \
	jlp_atan2.o $(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) -lm
	rm $*.o


.c.o:
	cc -c $(CFLAGS) $*.c

.o.exe:
	cc -o $(EXEC)/$*.exe $*.o \
	$(JLIB) $(FITSLIB) $(F77LIB) $(XLIB) -lm

.c.exe:
	cc -c $(CFLAGS) $*.c
	cc -o $(EXEC)/$*.exe $*.o \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) $(XLIB) -lm

all: $(OUT)

decode_set.o: decode_set.c
	cc -c $(CFLAGS) $*.c

decode_ima.exe : decode_ima.c 
	cc -c $(CFLAGS) decode_ima.c
	cc $(CFLAGS) -o $(EXEC)/decode_ima.exe decode_ima.o $(OBJ3) \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) -lm

decode_simu2.exe : decode_simu2.c $(OBJ3) 
	cc -c $(CFLAGS) decode_simu2.c
	cc $(CFLAGS) -o $(EXEC)/decode_simu2.exe decode_simu2.o $(OBJ3) \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) -lm

decode_ttau.exe : decode_ttau.c $(OBJ1) 
	cc -c $(CFLAGS) decode_ttau.c
	cc $(CFLAGS) -o $(EXEC)/decode_ttau.exe decode_ttau.o $(OBJ1) \
	$(JLIB) $(FITSLIB) $(F77LIB) -lm

# JLP2000: with FFTW
decode_ima_1D.exe : decode_ima_1D.c $(OBJ3) 
	cc -c $(CFLAGS) decode_ima_1D.c
	cc $(CFLAGS) -o $(EXEC)/decode_ima_1D.exe decode_ima_1D.o $(OBJ3) \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) -lm

decode_1D_fits_cube.exe : decode_ima_1D.c $(OBJ3) 
	cc -c $(CFLAGS) -DFITS_CUBE decode_ima_1D.c
	cc $(CFLAGS) -o $(EXEC)/decode_1D_fits_cube.exe decode_ima_1D.o $(OBJ3)\
	$(JLIB) $(FITSLIB) $(F77LIB) -lm

decode_fits_cube.exe : decode_ima.c $(OBJ3) 
	cc -c $(CFLAGS) -DFITS_CUBE decode_ima.c
	cc $(CFLAGS) -o $(EXEC)/decode_fits_cube.exe decode_ima.o $(OBJ3) \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) -lm

decode_ima_eric.exe : decode_ima.c $(OBJ3) 
	cc -c $(CFLAGS) -DERIC_FORMAT decode_ima.c
	cc $(CFLAGS) -o $(EXEC)/decode_ima_eric.exe decode_ima.o $(OBJ3) \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) -lm

bisp_from_image.exe : bisp_from_image.c $(OBJ3) 
	cc -c $(CFLAGS) bisp_from_image.c
	cc $(CFLAGS) -o $(EXEC)/bisp_from_image.exe bisp_from_image.o $(OBJ3) \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) -lm

correct_bisp.exe : correct_bisp.c $(OBJ3) 
	cc -c $(CFLAGS) correct_bisp.c
	c++ $(CFLAGS) -o $(EXEC)/correct_bisp.exe correct_bisp.o $(OBJ3) \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) -lm

correct_bisp_auto1.o : correct_bisp_auto1.c 

bisp_to_image.exe : bisp_to_image.c $(OBJ3) 
	cc -c $(CFLAGS) bisp_to_image.c
	cc $(CFLAGS) -o $(EXEC)/bisp_to_image.exe bisp_to_image.o $(OBJ3) \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) -lm

inv_bi_98.exe : inv_bi_98.c $(OBJ3) 
	cc -c $(CFLAGS) inv_bi_98.c
	cc $(CFLAGS) -o $(EXEC)/inv_bi_98.exe inv_bi_98.o $(OBJ3) \
	$(JLIB) $(FITSLIB) $(F77LIB) -lm

inv_bispec2.exe : inv_bispec2.c 
	cc -c $(CFLAGS) inv_bispec2.c
	c++ $(CFLAGS) -o $(EXEC)/inv_bispec2.exe inv_bispec2.o $(OBJ3) \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) -lm

inv_bispec2_1D.exe : inv_bispec2_1D.c $(OBJ3) 
	cc -c $(CFLAGS) inv_bispec2_1D.c
	cc $(CFLAGS) -o $(EXEC)/inv_bispec2_1D.exe inv_bispec2_1D.o $(OBJ3) \
	$(JLIB) $(FITSLIB) $(F77LIB) -lm

jlp_cover0.o : jlp_cover0.for

jlp_cover_mask.o : jlp_cover_mask.c

jlp_bispec1.o : jlp_bispec1.for

fft_jlp.o : fft_jlp.for

photon_model.o : photon_model.c

fft2.exe: fft2.for jlp_atan2.o jlp_cover_mask.o

jlp_atan2.o: jlp_atan2.for

fft2_norm.exe: fft2_norm.for

fft_1D.exe: fft_1D.for


nagfft.o : nagfft.for

bisp_to_eric.exe : bisp_to_eric.c $(OBJ3) 
	cc -c $(CFLAGS) bisp_to_eric.c
	cc $(CFLAGS) -o $(EXEC)/bisp_to_eric.exe bisp_to_eric.o $(OBJ3) \
	$(JLIB) $(FITSLIB) $(F77LIB) -lm

kolmo.o : kolmo.c

poidev.o : poidev.c

simu2.exe : simu2.c poidev.o $(FOURN)

simu3.exe : simu3.c kolmo.o poidev.o $(FOURN)

kolmo.exe : kolmo.c $(FOURN)
