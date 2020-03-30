######################################################################
# Makefile for decode...
# Pic du Midi version
#
# JLP
# Version 21-03-96
######################################################################
OBJ2=jlp_bispec1.o jlp_cover.o fft_nag.o nagfft.o
OBJ1=jlp_cover2.o fft_nag.o nagfft.o
OBJ3=jlp_cover_mask.o fft_nag.o nagfft.o
OBJ3=jlp_cover_mask.o fft_jlp.o
OUT= $(OBJ1) decode_simu2.exe
midir=$(JLPLIB)/midas
mylib=$(JLPLIB)/jlp
JLIB=$(mylib)/jlpacc.a $(mylib)/jlputil.a \
     $(mylib)/newplot0.a $(mylib)/jlp_splot1.a $(mylib)/jlp_splot2.a \
     $(mylib)/jlpacc.a $(JLPLIB)/math/mymath.a
F77= xlf
F77= f77 
# Sun and Dec:
#F77LIB=$(F77DIR)/libF77.a $(F77DIR)/libI77.a $(F77DIR)/libU77.a 
# IBM:
#F77LIB= $(F77DIR)/libxlf.a
#Pic du Midi:
F77LIB= -lfor -lutil -lots -li 
MIDLIB=$(midir)/stlib.a $(midir)/udiolib.a $(midir)/oslib.a
ESOEXT=/midas/frozen/exec/esoext.exe
INC=-I../midincl
#FFLAGS= -O 
FFLAGS= -g 
#CCLAGS= -g -C -D$(JLPSYSTEM) 
CCLAGS= -g -D$(JLPSYSTEM) 
#CCLAGS= -O -D$(JLPSYSTEM) 

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
	cc -c $(CCLAGS) $(INC) $*.c

.c.exe:
	cc -c $(CCLAGS) $(INC) $*.c
	cc $(CCLAGS) -o $(EXEC)/$*.exe $*.o $(OBJ1) \
	jlp_atan2.o $(JLIB) $(MIDLIB) $(F77LIB) -lm

all: $(OUT)

decode_simu2.exe : decode_simu2.c $(OBJ1) 
	cc -c $(CCLAGS) $(INC) decode_simu2.c
	cc $(CCLAGS) -o $(EXEC)/decode_simu2.exe decode_simu2.o $(OBJ1) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

decode_ttau.exe : decode_ttau.c $(OBJ1) 
	cc -c $(CCLAGS) $(INC) decode_ttau.c
	cc $(CCLAGS) -o $(EXEC)/decode_ttau.exe decode_ttau.o $(OBJ1) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

decode_ima.exe : decode_ima.c $(OBJ1) 
	cc -c $(CCLAGS) $(INC) decode_ima.c
	cc $(CCLAGS) -o $(EXEC)/decode_ima.exe decode_ima.o $(OBJ1) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

decode_car_1D.exe : decode_car_1D.c $(OBJ3) 
	cc -c $(CCLAGS) $(INC) decode_car_1D.c
	cc $(CCLAGS) -o $(EXEC)/decode_car_1D.exe decode_car_1D.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

decode_car2a.exe : decode_car2.c $(OBJ3) 
	cc -c $(CCLAGS) $(INC) decode_car2.c
	cc $(CCLAGS) -o $(EXEC)/decode_car2a.exe decode_car2.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

decode_cp40.exe : decode_cp40.c $(OBJ3) lk_fmt2.c 
	cc -c $(CCLAGS) $(INC) decode_cp40.c
	cc $(CCLAGS) -o $(EXEC)/decode_cp40.exe decode_cp40.o $(OBJ3) lk_fmt2.o \
        $(JLIB) $(MIDLIB) $(F77LIB) -lm

phot_noise_mask.exe : phot_noise_mask.c $(OBJ3) 
	cc -c $(CCLAGS) $(INC) phot_noise_mask.c
	cc $(CCLAGS) -o $(EXEC)/phot_noise_mask.exe phot_noise_mask.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

phot_noise_1D.exe : phot_noise_1D.c $(OBJ3) 
	cc -c $(CCLAGS) $(INC) phot_noise_1D.c
	cc $(CCLAGS) -o $(EXEC)/phot_noise_1D.exe phot_noise_1D.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

dec_car3.exe : dec_car3.c $(OBJ1) 
	cc -c $(CCLAGS) $(INC) dec_car3.c
	cc $(CCLAGS) -o $(EXEC)/dec_car3.exe dec_car3.o $(OBJ1) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

dec_car4.exe : dec_car4.c $(OBJ1) 
	cc -c $(CCLAGS) $(INC) dec_car4.c
	cc $(CCLAGS) -o $(EXEC)/dec_car4.exe dec_car4.o $(OBJ1) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

test.exe : test.c $(OBJ1) 
	cc -c $(CCLAGS) $(INC) test.c
	cc $(CCLAGS) -o $(EXEC)/test.exe test.o $(OBJ1) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

decode_test.exe : decode_car2.c $(OBJ3) 
	cc -c $(CCLAGS) $(INC) decode_car2.c
	cc $(CCLAGS) -o $(EXEC)/decode_test.exe decode_car2.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

bisp_to_image.exe : bisp_to_image.c $(OBJ3) 
	cc -c $(CCLAGS) $(INC) bisp_to_image.c
	cc $(CCLAGS) -o $(EXEC)/bisp_to_image.exe bisp_to_image.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

inv_bispec2.exe : inv_bispec2.c $(OBJ3) 
	cc -c $(CCLAGS) $(INC) inv_bispec2.c
	cc $(CCLAGS) -o $(EXEC)/inv_bispec2.exe inv_bispec2.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

inv_bispec2_1D.exe : inv_bispec2_1D.c $(OBJ3) 
	cc -c $(CCLAGS) $(INC) inv_bispec2_1D.c
	cc $(CCLAGS) -o $(EXEC)/inv_bispec2_1D.exe inv_bispec2_1D.o $(OBJ3) \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

jlp_cover.o : jlp_cover.for

jlp_cover2.o : jlp_cover2.c

jlp_cover_mask.o : jlp_cover_mask.c

jlp_bispec1.o : jlp_bispec1.for

fft_jlp.o : fft_jlp.for

nagfft.o : nagfft.for

bisp_to_eric.exe : bisp_to_eric.c jlp_cover_mask.o $(OBJ2) 
	cc -c $(CCLAGS) $(INC) bisp_to_eric.c
	cc $(CCLAGS) -o $(EXEC)/bisp_to_eric.exe bisp_to_eric.o $(OBJ2) \
	jlp_cover_mask.o \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm
