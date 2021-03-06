######################################################################
# Makefile for decode...
#
# JLP
# Version 14-03-99
######################################################################
include $(JLPSRC)/jlp_make.mk
mylib=$(JLPLIB)/jlp
OBJ2=jlp_bispec1.o jlp_cover.o fft_nag.o nagfft.o
OBJ1=jlp_cover2.o fft_nag.o nagfft.o
OBJ3_APO=jlp_cover_mask.o fft_nag.o nagfft.o lk_fmt2.o decode_set_apo.o
#OBJ3=jlp_cover_mask.o fft_nag.o nagfft.o lk_fmt2.o decode_set.o
OBJ3=jlp_cover_mask.o decode_set.o
OBJ33=jlp_cover_mask.o fft_nag.o nagfft.o lk_fmt2.o decode_set.o
#OBJ4=jlp_cover_mask.o lk_fmt2.o decode_set.o $(mylib)/fftn.a
OBJ4=jlp_cover_mask.o lk_fmt2.o decode_set.o
OUT= $(OBJ3) correct_bisp.exe
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
	$(F77) -g -o $(EXEC)/$*.exe $*.o jlp_cover_mask.o \
	jlp_atan2.o $(JLIB) $(MIDLIB) $(F77LIB) -lm
	rm $*.o


.c.o:
	cc -c $(CFLAGS) $*.c

.o.exe:
	cc -o $(EXEC)/$*.exe $*.o \
	$(JLIB) $(MIDLIB) $(F77LIB) $(XLIB) -lm

.c.exe:
	cc -c $(CFLAGS) $*.c
	cc -o $(EXEC)/$*.exe $*.o \
	$(JLIB) $(MIDLIB) $(F77LIB) $(XLIB) -lm

all: $(OUT)

decode_set.o: decode_set.c
	cc -c $(CFLAGS) $*.c

decode_simu2.exe : decode_simu2.c $(OBJ1) 
	cc -c $(CFLAGS) decode_simu2.c
	cc $(CFLAGS) -o $(EXEC)/decode_simu2.exe decode_simu2.o $(OBJ1) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

decode_ttau.exe : decode_ttau.c $(OBJ1) 
	cc -c $(CFLAGS) decode_ttau.c
	cc $(CFLAGS) -o $(EXEC)/decode_ttau.exe decode_ttau.o $(OBJ1) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

# JLP2000: with FFTW
decode_ima_1D.exe : decode_ima_1D.c $(OBJ3) 
	cc -c $(CFLAGS) decode_ima_1D.c
	cc $(CFLAGS) -o $(EXEC)/decode_ima_1D.exe decode_ima_1D.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(MATHLIB) $(F77LIB) -lm

decode_1D_fits_cube.exe : decode_ima_1D.c $(OBJ3) 
	cc -c $(CFLAGS) -DFITS_CUBE decode_ima_1D.c
	cc $(CFLAGS) -o $(EXEC)/decode_1D_fits_cube.exe decode_ima_1D.o $(OBJ3)\
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

decode_fits_cube.exe : decode_ima.c $(OBJ3) 
	cc -c $(CFLAGS) -DFITS_CUBE decode_ima.c
	cc $(CFLAGS) -o $(EXEC)/decode_fits_cube.exe decode_ima.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(MATHLIB) $(F77LIB) -lm

decode_ima_eric.exe : decode_ima.c $(OBJ3) 
	cc -c $(CFLAGS) -DERIC_FORMAT decode_ima.c
	cc $(CFLAGS) -o $(EXEC)/decode_ima_eric.exe decode_ima.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(MATHLIB) $(F77LIB) -lm

decode_ima.exe : decode_ima.c 
	cc -c $(CFLAGS) decode_ima.c
	cc $(CFLAGS) -o $(EXEC)/decode_ima.exe decode_ima.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(MATHLIB) $(F77LIB) -lm

decode_cp40.exe : decode_cp40.c $(OBJ3) lk_fmt2.o autoc_set.o 
	cc -c $(CFLAGS) decode_cp40.c
	cc $(CFLAGS) -o $(EXEC)/decode_cp40.exe decode_cp40.o $(OBJ4) \
	autoc_set.o \
        $(JLIB) $(MIDLIB) $(MATHLIB) $(F77LIB) -lm

test.exe : decode_cp40.c $(OBJ4) lk_fmt2.o autoc_set.o 
	cc -c $(CFLAGS) decode_cp40.c
	cc -bloadmap:tt $(CFLAGS) -o $(EXEC)/decode_cp40.exe decode_cp40.o $(OBJ4) \
	autoc_set.o $(JLIB) $(MIDLIB) -lm

decode_cp40_apo.exe : decode_cp40.c $(OBJ3) lk_fmt2.o autoc_set.o 
	cp decode_set.c decode_set_apo.c
	cc -c $(CFLAGS) -DAPODIZATION decode_set_apo.c
	rm -f decode_set_apo.c
	cc -c $(CFLAGS) -DAPODIZATION decode_cp40.c
	cc $(CFLAGS) -o $(EXEC)/decode_cp40_apo.exe decode_cp40.o $(OBJ3_APO) \
	autoc_set.o \
        $(JLIB) $(MIDLIB) $(F77LIB) -lm

decode_cp400.exe : decode_cp400.c $(OBJ3) lk_fmt2.o autoc_set.o 
	cc -c $(CFLAGS) decode_cp400.c
	cc $(CFLAGS) -o $(EXEC)/decode_cp400.exe decode_cp400.o $(OBJ33) \
	autoc_set.o \
        $(JLIB) $(MIDLIB) $(F77LIB) -lm

correct_bisp.exe : correct_bisp.c $(OBJ3) 
	cc -c $(CFLAGS) correct_bisp.c
	c++ $(CFLAGS) -o $(EXEC)/correct_bisp.exe correct_bisp.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

correct_bisp_hege.exe : correct_bisp_hege.c $(OBJ3) 
	cc -c $(CFLAGS) correct_bisp_hege.c
	cc $(CFLAGS) -o $(EXEC)/correct_bisp_hege.exe correct_bisp_hege.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

phot_noise_ma2.exe : phot_noise_ma2.c $(OBJ3) 
	cc -c $(CFLAGS) phot_noise_ma2.c
	cc $(CFLAGS) -o $(EXEC)/phot_noise_ma2.exe phot_noise_ma2.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

phot_noise_mask.exe : phot_noise_mask.c 
	cc -c $(CFLAGS) phot_noise_mask.c
	cc $(CFLAGS) -o $(EXEC)/phot_noise_mask.exe phot_noise_mask.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(MATHLIB) $(F77LIB) -lm

phot_noise_1D.exe : phot_noise_1D.c $(OBJ3) 
	cc -c $(CFLAGS) phot_noise_1D.c
	cc $(CFLAGS) -o $(EXEC)/phot_noise_1D.exe phot_noise_1D.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

#dec_car3.exe : dec_car3.c $(OBJ1)  ---> (now autoc_car.exe) 
autoc_car.exe : autoc_car.c $(OBJ3) autoc_set.o 
	cc -c $(CFLAGS) autoc_car.c
	cc $(CFLAGS) -o $(EXEC)/autoc_car.exe autoc_car.o $(OBJ3) \
	autoc_set.o $(JLIB) $(MIDLIB) $(F77LIB) -lm

autoc_set.o : autoc_set.c

decode_car_1D.exe : decode_car_1D.c $(OBJ3) autoc_set.o 
	cc -c $(CFLAGS) decode_car_1D.c
	cc $(CFLAGS) -o $(EXEC)/decode_car_1D.exe decode_car_1D.o $(OBJ3) \
	autoc_set.o $(JLIB) $(MIDLIB) $(F77LIB) -lm

decode_car2.exe : decode_car2.c $(OBJ3) autoc_set.o 
	cc -c $(CFLAGS) decode_car2.c
	cc $(CFLAGS) -o $(EXEC)/decode_car2.exe decode_car2.o $(OBJ3) \
	autoc_set.o $(JLIB) $(MIDLIB) $(MATHLIB) $(F77LIB) -lm

saad_car.exe : saad_car.c $(OBJ3) autoc_set.o 
	cc -c $(CFLAGS) saad_car.c
	cc $(CFLAGS) -o $(EXEC)/saad_car.exe saad_car.o $(OBJ3) \
	autoc_set.o $(JLIB) $(MIDLIB) $(F77LIB) -lm

dec_car4.exe : dec_car4.c $(OBJ1) autoc_set.o 
	cc -c $(CFLAGS) dec_car4.c
	cc $(CFLAGS) -o $(EXEC)/dec_car4.exe dec_car4.o $(OBJ1) \
	autoc_set.o $(JLIB) $(MIDLIB) $(F77LIB) -lm

decode_test.exe : decode_car2.c $(OBJ3) 
	cc -c $(CFLAGS) decode_car2.c
	cc $(CFLAGS) -o $(EXEC)/decode_test.exe decode_car2.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

bisp_to_image.exe : bisp_to_image.c $(OBJ3) 
	cc -c $(CFLAGS) bisp_to_image.c
	cc $(CFLAGS) -o $(EXEC)/bisp_to_image.exe bisp_to_image.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

inv_bi_98.exe : inv_bi_98.c $(OBJ3) 
	cc -c $(CFLAGS) inv_bi_98.c
	cc $(CFLAGS) -o $(EXEC)/inv_bi_98.exe inv_bi_98.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

inv_bispec2.exe : inv_bispec2.c 
	cc -c $(CFLAGS) inv_bispec2.c
	cc $(CFLAGS) -o $(EXEC)/inv_bispec2.exe inv_bispec2.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

inv_bispec2_1D.exe : inv_bispec2_1D.c $(OBJ3) 
	cc -c $(CFLAGS) inv_bispec2_1D.c
	cc $(CFLAGS) -o $(EXEC)/inv_bispec2_1D.exe inv_bispec2_1D.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

jlp_cover.o : jlp_cover.for

jlp_cover2.o : jlp_cover2.c

jlp_cover_mask.o : jlp_cover_mask.c

jlp_bispec1.o : jlp_bispec1.for

fft_jlp.o : fft_jlp.for

fft2.exe: fft2.for jlp_atan2.o jlp_cover_mask.o

jlp_atan2.o: jlp_atan2.for

fft2_norm.exe: fft2_norm.for

fft_1D.exe: fft_1D.for


nagfft.o : nagfft.for

bisp_to_eric.exe : bisp_to_eric.c jlp_cover_mask.o $(OBJ2) 
	cc -c $(CFLAGS) bisp_to_eric.c
	cc $(CFLAGS) -o $(EXEC)/bisp_to_eric.exe bisp_to_eric.o $(OBJ2) \
	jlp_cover_mask.o \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

kolmo.o : kolmo.c

poidev.o : poidev.c

simu2.exe : simu2.c poidev.o $(FOURN)

simu3.exe : simu3.c kolmo.o poidev.o $(FOURN)

kolmo.exe : kolmo.c $(FOURN)
