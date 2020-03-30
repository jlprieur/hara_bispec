/********************************************************************
* compute_phase.c
*
* Conversion of real and imaginary part of a spectrum 
* to phase and modulus for further display 
*
* JLP
* Version 03-09-2008
*******************************************************************/
#include <stdio.h>
#include <math.h>
#include <jlp_ftoc.h>
#include <jlp_cover_mask.h>

#define DEBUG

/* Contained here: */
int compute_phase_and_sqmod(float *in_real, float *in_imag, float *out_phase, 
                            float *out_sqmod, int nx, int ny, float backgd); 
int clean_the_phase(float *InPhase, float *OutPhase, int nx, int ny, 
                    float backgd);
int unwrap_the_phase(float *InPhase, float *OutPhase, int nx, int ny, 
                    float backgd);

int main(int argc, char *argv[])
{
float *in_real, *in_imag, *out_phase, *out_sqmod;
float *cleaned_phase, *unwrapped_phase, backgd;
INT4 nx, ny, nx_imag, ny_imag; 
INT_PNTR pntr1;
int isize;
char real_name[61], imag_name[61], phase_name[61], sqmod_name[61];
char comments[81];

printf(" compute_phase           -- Version 03-09-2008 --\n");
printf(" Conversion of real and imag part to phase and square modulus\n");

/* Carefull: 7 parameters always, using JLP "runs" */
/* If command line with "runs" */
if(argc == 7){
 if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
 else argc = 1;
 }

if(argc != 5) {
  printf("        Fatal error: Wrong syntax, argc=%d\n",argc);
  printf(" Syntax is:  \n");
  printf(" runs compute_phase real imag phase modulus");
 }

strcpy(real_name, argv[1]);
strcpy(imag_name, argv[2]);
strcpy(phase_name, argv[3]);
strcpy(sqmod_name, argv[4]);

printf("*** INPUT: real=>%s< imag=>%s< \n*** OUTPUT: phase=>%s< square modulus=>%s<\n",
       real_name, imag_name, phase_name, sqmod_name);

/***************************************************************/
JLP_INQUIFMT();

/* Read real part */
  JLP_VM_READIMAG1(&pntr1,&nx,&ny,real_name,comments);
  in_real = (float *)pntr1;

/* Read imaginary part */
  JLP_VM_READIMAG1(&pntr1,&nx_imag,&ny_imag,imag_name,comments);
  in_imag = (float *)pntr1;

/* Check size is correct: */
  if((nx != nx_imag) || (ny != ny_imag))
    {
    fprintf(stderr, 
            "Fatal error: real and imaginary parts have inconsistent sizes!\n");
    return(-1);
    } 

/* Get memory space: */
   isize = nx * ny * sizeof(float);
   JLP_GVM(&out_phase, &isize);
   JLP_GVM(&cleaned_phase, &isize);
   JLP_GVM(&unwrapped_phase, &isize);
   JLP_GVM(&out_sqmod, &isize);

/* I set the background to backgd =  -2 PI for an easier display: */
   backgd = -2.* PI;

/* Compute phase and modulus: */
   compute_phase_and_sqmod(in_real, in_imag, out_phase, out_sqmod, nx, ny, 
                           backgd); 
/* Clean the holes in the phase, if any: */
   clean_the_phase(out_phase, cleaned_phase, nx, ny, backgd);
   unwrap_the_phase(cleaned_phase, unwrapped_phase, nx, ny, backgd);

/* Output phase and square modulus: */
/*
  sprintf(comments,"Cleaned phase (%.20s and %.20s)", real_name, imag_name);
  strcpy(phase_name,"cleaned_phase");
  JLP_WRITEIMAG(cleaned_phase,&nx,&ny,&nx,phase_name,comments);

  sprintf(comments,"Unwrapped phase (%.20s and %.20s)", real_name, imag_name);
  strcpy(phase_name,"unwrapped_phase");
  JLP_WRITEIMAG(unwrapped_phase,&nx,&ny,&nx,phase_name,comments);

  sprintf(comments,"Phase from %.20s and %.20s", real_name, imag_name);
  JLP_WRITEIMAG(out_phase,&nx,&ny,&nx,phase_name,comments);
*/
  sprintf(comments,"Unwrapped phase (%.20s and %.20s)", real_name, imag_name);
  JLP_WRITEIMAG(unwrapped_phase,&nx,&ny,&nx,phase_name,comments);

  sprintf(comments,"Square modulus from %.20s and %.20s", real_name, imag_name);
  JLP_WRITEIMAG(out_sqmod,&nx,&ny,&nx,sqmod_name,comments);

/* Free memory: */
  JLP_FVM(&out_phase);
  JLP_FVM(&cleaned_phase);
  JLP_FVM(&out_sqmod);

return(0);
}
/*****************************************************************
* Compute phase and modulus: 
*
* INPUT:
* in_real, in_imag: real and imaginary parts
* nx, ny: size of arrays
*
* OUTPUT:
* out_phase, out_sqmod: phase and square modulus
*
******************************************************************/
int compute_phase_and_sqmod(float *in_real, float *in_imag, float *out_phase, 
                            float *out_sqmod, int nx, int ny, float backgd)
{
register int i;

/* Initialize the arrays: */
   for(i = 0; i < nx * ny; i++) {
/* Power spectrum: */
     out_sqmod[i] = SQUARE(in_real[i]) + SQUARE(in_imag[i]);
/* Phase: */
     if(in_imag[i] == 0.) {
         if(in_real[i] > 0) out_phase[i] = 0.;
         else if(in_real[i] < 0) out_phase[i] = PI; 
/* I set the background to backgd =  -2 PI for an easier display: */
         else
           out_phase[i] = backgd;
     } else {
       out_phase[i] = atan2(in_imag[i], in_real[i]);
     }
     }

return(0);
}
/*****************************************************************
* Morphological cleaning 
*
* INPUT:
* InPhase: phase before processing 
* nx, ny: size of arrays
*
* OUTPUT:
* OutPhase: phase after processing 
*
******************************************************************/
int clean_the_phase(float *InPhase, float *OutPhase, int nx, int ny, 
                    float backgd)
{
float g_left, g_right, g_up, g_down, threshold = PI;
register int i, j;

for(i = 0; i < nx * ny; i++) OutPhase[i] = InPhase[i]; 

/* Fills the holes: 
*/
for(j = 1; j < ny-1; j++) { 
  for(i = 1; i < nx-1; i++) {
    if(InPhase[i + j * nx] == backgd) {
      g_left = ABS(InPhase[i + j * nx] - InPhase[i-1 + j * nx]);
      g_right = ABS(InPhase[i+1 + j * nx] - InPhase[i + j * nx]);
      g_up = ABS(InPhase[i + (j+1) * nx] - InPhase[i + j * nx]);
      g_down = ABS(InPhase[i + j * nx] - InPhase[i + (j-1) * nx]);
      if(g_left > threshold && g_right > threshold 
         && g_up > threshold && g_down > threshold) {
         OutPhase[i + j * nx] = (InPhase[i-1 + j * nx] + InPhase[i+1 + j * nx]
                               + InPhase[i + (j+1) * nx] 
                               + InPhase[i + (j-1) * nx])/4.;
      }
    }
  }
}

/* Removes isolated points: 
*/
for(j = 1; j < ny-1; j++) { 
  for(i = 1; i < nx-1; i++) {
    if((InPhase[i + j * nx] != backgd) && (InPhase[i-1 + j * nx] == backgd)
      && (InPhase[i+1 + j * nx] == backgd) 
      && (InPhase[i + (j-1) * nx] == backgd)
      && (InPhase[i + (j+1) * nx] == backgd)) {
         OutPhase[i + j * nx] = backgd; 
     
    }
  }
}

return(0);
}
/*****************************************************************
* Unwrapping of the phase 
*
* INPUT:
* InPhase: phase before processing 
* nx, ny: size of arrays
*
* OUTPUT:
* OutPhase: phase after processing 
*
******************************************************************/
int unwrap_the_phase(float *InPhase, float *OutPhase, int nx, int ny, 
                    float backgd)
{
float *tmp, g_left, g_up, threshold = 1.2*PI, c0;
float previous_good_invalue;
register int i, j;

tmp = malloc(nx * ny * sizeof(float));

for(i = 0; i < nx * ny; i++) tmp[i] = InPhase[i]; 
for(i = 0; i < nx * ny; i++) OutPhase[i] = InPhase[i]; 

/* Main loop, along the lines,
* - goes from left to right, 
* - identifies 2*PI discontinuities
* - and corrects them by adjusting c0
*/
for(j = 1; j < ny-1; j++) { 
  c0 = 0.;
  previous_good_invalue = backgd;
  for(i = 1; i < nx-1; i++) {
/* General case, add c0 to the valid phases: */
    if(tmp[i + j * nx] != backgd) {
      OutPhase[i + j * nx] += c0;
      previous_good_invalue = tmp[i + j * nx];
      }
/* Check if c0 needs to be incremented: */
    if((previous_good_invalue != backgd) && (tmp[i+1 + j * nx] != backgd)) {
    if(tmp[i + j * nx] != backgd) 
       g_left = tmp[i+1 + j * nx] - tmp[i + j * nx];
    else
       g_left = tmp[i+1 + j * nx] - previous_good_invalue;
    if(ABS(g_left) > threshold) {
      if(g_left > 0.) {
        c0 -= 2. * PI;
     } else if(g_left < 0.) {
        c0 += 2. * PI;
     }
     }
/* DEBUG:
if(j == ny/2)printf("i=%d g_left=%f c0=%f (%f %f)\n", i, g_left, c0,
                     tmp[i + j * nx], tmp[i+1 + j * nx]);
*/
   }
  }
}

/* Final checkup, along the columns, from middle to top: */
for(i = 0; i < nx * ny; i++) tmp[i] = OutPhase[i]; 

for(i = 1; i < nx-1; i++) {
  c0 = 0.;
  previous_good_invalue = backgd;
  for(j = ny/2; j < ny-1; j++) { 
/* General case, add c0 to the valid phases: */
    if(OutPhase[i + j * nx] != backgd) {
      OutPhase[i + j * nx] += c0;
      previous_good_invalue = tmp[i + j * nx];
      }
/* Check if c0 needs to be incremented: */
    if((previous_good_invalue != backgd) && (tmp[i + (j+1) * nx] != backgd)) {
    if(tmp[i + j * nx] != backgd) 
       g_up = tmp[i + (j+1) * nx] - tmp[i + j * nx];
    else
       g_up = tmp[i + (j+1) * nx] - previous_good_invalue;
    if(ABS(g_up) > threshold) {
      if(g_up > 0.) {
        c0 -= 2. * PI;
     } else if(g_up < 0.) {
        c0 += 2. * PI;
     }
     }
/* DEBUG:
if(i == 63)printf("val(i=%d,j=%d)= %f val(i,j+1)=%f g_up=%f c0=%f\n", i, j, 
                  tmp[i + j * nx], tmp[i + (j+1) * nx], g_up, c0);
*/
   }
  }
}

/* Final checkup, along the columns, from middle to bottom: */
for(i = 0; i < nx * ny; i++) tmp[i] = OutPhase[i]; 

for(i = 1; i < nx-1; i++) {
  c0 = 0.;
  previous_good_invalue = backgd;
  for(j = ny/2; j > 1; j--) { 
/* General case, add c0 to the valid phases: */
    if(OutPhase[i + j * nx] != backgd) {
      OutPhase[i + j * nx] += c0;
      previous_good_invalue = tmp[i + j * nx];
      }
/* Check if c0 needs to be incremented: */
    if((previous_good_invalue != backgd) && (tmp[i + (j-1) * nx] != backgd)) {
    if(tmp[i + j * nx] != backgd) 
       g_up = tmp[i + (j-1) * nx] - tmp[i + j * nx];
    else
       g_up = tmp[i + (j-1) * nx] - previous_good_invalue;
    if(ABS(g_up) > threshold) {
      if(g_up > 0.) {
        c0 -= 2. * PI;
     } else if(g_up < 0.) {
        c0 += 2. * PI;
     }
     }
/* DEBUG:
if(i == 63)printf("val(i=%d,j=%d)= %f val(i,j+1)=%f g_up=%f c0=%f\n", i, j, 
                  tmp[i + j * nx], tmp[i + (j+1) * nx], g_up, c0);
*/
   }
  }
}

free(tmp);

return(0);
}
