/**************************************************************** 
*  correct_bisp_hege
*  Same as correct_bisp.c, but with hege's method 
*  (fits to the power spectrum a parametric model of the photon response) 
*
* JLP
* Version 16-04-2008
*****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <jlp_ftoc.h>
#include "jlp_cover_mask.h"

/* Defined here: */
int photon_corr_mask_hege(double *bispp, double *modsq, double *phot_modsq,
                          INT4 *nx, INT4 *ny, INT4 *bisp_dim, float *xphot,
                          INT4 *nbeta, INT4 *ngamma);

int main(int argc, char *argv[])
{
/* Images : */
float *long_int, *modsq1, *bisp1, *ffield, *phot_modsq1, *mask;
double *modsq, *bispp, *phot_modsq;
float xphot, sky_level;
INT_PNTR pntr_ima;
INT4 nx_bisp, ny_bisp, nx, ny, nx1, ny1, nx_mask, ny_mask;
INT4 ir, nbeta, ngamma, bisp_dim, isize, max_nclosure;
register int i;
char *pc, in_prefix[41], out_prefix[41];
char filename[61], ff_name[61], comments[81];
char phot_modsq_name[61], mask_name[61];

printf(" Program correct_bisp_hege \n");
printf(" JLP Version 16-04-2008 \n");
printf(" For the extension, either _tt or tt_\n");
printf("(i.e. either modsq_tt, bisp1_tt, ... or tt_m, tt_b, t_l\n");


/* One or three parameters only are allowed to run the program: */
/* Carefull: 7 parameters always, using JLP "runs" */
#ifdef DEBUG
  printf(" argc=%d\n",argc);
  printf(" argv[3]=>%s<\n",argv[3]);
#endif
if((argc > 1 && argc < 6) || (argc == 7 && !strcmp(argv[3],"")))
  {
  printf("        Fatal error: Wrong syntax, argc=%d\n",argc);
  printf(" Syntax is:  \n");
  printf("runs correct_bisp in_prefix out_prefix phot_modsq flat_field ir,max_nclo nphot_per_frame\n");
  printf("or runs correct_bisp in_prefix out_prefix phot_modsq 0 ir,max_nclo nphot_per_frame\n");
  printf("or runs correct_bisp in_prefix out_prefix 0 0 ir,max_nclo nphot_per_frame \n");
  exit(-1);
  }

/* Interactive input of in_prefix and out_prefix: */
if (argc == 1)
 {
   printf(" Input file extension := "); scanf("%s",in_prefix);
/************* File extension for output files: */
   printf(" Output file extension := "); scanf("%s",out_prefix);
   printf(" Photon modsq (0 if not available) := "); scanf("%s",phot_modsq_name);
   printf(" Flat field (0 if unity flat field) := "); scanf("%s",ff_name);
   printf(" Radius of uv-coverage (IR) in pixels, max_nclosure: (12,1000 f.i.) \n");
   scanf("%d,%d",&ir,&max_nclosure);
   printf(" Number of photons per frame: \n");
   scanf("%f",&xphot);
 }
else
 {
  strcpy(in_prefix,argv[1]);
  strcpy(out_prefix,argv[2]);
  strcpy(phot_modsq_name,argv[3]);
  strcpy(ff_name,argv[4]);
  sscanf(argv[5],"%d,%d",&ir,&max_nclosure);
  sscanf(argv[6],"%f",&xphot);
 }

/* End extension strings with zero: */
  pc = in_prefix;
  while(*pc && *pc != ' ') pc++;
  *pc='\0';
  pc = out_prefix;
  while(*pc && *pc != ' ') pc++;
  *pc='\0';

#ifdef DEBUG
printf(" DEBUG Version, will read >%s< files and write >%s< \n",in_prefix);
#endif

/***************************************************************/
JLP_INQUIFMT();

/*****************************************************************/
/* Computing the spectral and bispectral lists
   corresponding to the selected uv coverage: */
/* JLP2000: I set the option to "unmasked uv-coverage" */
  strcpy(mask_name,"0");
  if(mask_name[0] == '0')
    {
/* Full u-v coverage: */
     printf(" OK: full u-v coverage.\n");
/* Allocate memory: */
     nx_mask = 2 * (ir + 2);
     ny_mask = nx_mask;
     isize = nx_mask * ny_mask * sizeof(float);
     mask = (float *)malloc(isize);
/* Create a filled mask: */
     for(i = 0; i < nx_mask * ny_mask; i++) mask[i] = 1.;
    }
  else
/* Read Fourier mask file: */
    {
     JLP_VM_READIMAG1(&pntr_ima,&nx_mask,&ny_mask,mask_name,comments);
     mask = (float *)pntr_ima;
    }

/* "Masked" u-v coverage: */
    COVERA_MASK(mask,&nx_mask,&ny_mask,&ir,&max_nclosure,&nbeta,&ngamma);

  printf(" Radius of uv-coverage (IR) in pixels: %d\n",ir);
  printf(" nbeta: %d ngamma: %d\n",nbeta,ngamma);
/* For debugging: 
output_lists_coverage(&nbeta,&ngamma);
*/

/* Free virtual memory space: */
    free(mask);

/* Mean squared modulus : */
   sprintf(filename,"%sm",in_prefix);
   printf(" Reading mean squared modulus: >%s< \n",filename);
   JLP_VM_READIMAG1(&pntr_ima,&nx,&ny,filename,comments);
   modsq1 = (float *) pntr_ima;

/* Mean bispectrum : */
   sprintf(filename,"%sb",in_prefix);
   printf(" Reading mean bispectrum: >%s< \n",filename);
   JLP_VM_READIMAG1(&pntr_ima,&nx_bisp,&ny_bisp,filename,comments);
   bisp1 = (float *) pntr_ima;
   if(nx_bisp != ngamma)
    {
     printf(" Fatal error: input bispectrum incompatible with u-v coverage!\n");
     exit(-1);
    }
   bisp_dim = ngamma;

/* Copy the three lines to double size array: */
   isize = bisp_dim * 3 * sizeof(double);
   bispp = (double *) malloc(isize);
   for(i = 0; i < 3*bisp_dim; i++) bispp[i] = (double)bisp1[i];

/******************************************************************/
/* Read flat field file: */
 if(ff_name[0] != '0')
  {
   printf(" Sorry this option (ffield) is not available yet \n"); exit(-1);
   JLP_VM_READIMAG1(&pntr_ima,&nx1,&ny1,ff_name,comments);
   if(nx1 != nx || ny1 != ny) {printf("uncompatible size for %s \n",ff_name);
                               exit(-1);}
   ffield = (float *) pntr_ima;

/* Free virtual memory space: */
   free(ffield);
  }
/********************* Unity file for ffield: *************/
  else
  {
   printf("\n OK, We assume unity ffield has been used\n");
  }

/* Compute the sky value from the long integration : */
   sprintf(filename,"%sl",in_prefix);
   JLP_VM_READIMAG1(&pntr_ima,&nx,&ny,filename,comments);
   long_int = (float *) pntr_ima;
   auto_sky(long_int,(int)nx,(int)ny,(int)nx,&sky_level);
   free(long_int);

/* Photon modsq (response to a photo-event: */
   JLP_VM_READIMAG1(&pntr_ima,&nx1,&ny1,phot_modsq_name,comments);
   if(nx1 != nx || ny1 != ny) 
         {printf("uncompatible size for %s \n",phot_modsq_name);
          exit(-1);}
   phot_modsq1 = (float *) pntr_ima;

/****************************************************************/
/* Recentre the frames: (necessary for photon_corr_mask) */
   RECENT_FFT(phot_modsq1,phot_modsq1,&nx,&ny,&nx);
   isize = nx * ny * sizeof(double);
   phot_modsq = (double *) malloc(isize);
   TO_DOUBLE(phot_modsq1,phot_modsq,&nx,&ny,&nx);
   RECENT_FFT(modsq1,modsq1,&nx,&ny,&nx);
   isize = nx * ny * sizeof(double);
   modsq = (double *) malloc(isize);
   TO_DOUBLE(modsq1,modsq,&nx,&ny,&nx);

/* Correction of central value of modsq because of a non null background */
   corr_bisp_sky(bispp,modsq,&nx,&ny,&bisp_dim,&sky_level,&nbeta,&ngamma);

/* Correction of photon bias: */
   photon_corr_mask_hege(bispp,modsq,phot_modsq,&nx,&ny,&bisp_dim,
                         &xphot,&nbeta,&ngamma);

/****************************************************************/
/* Now output the results : */
   sprintf(comments,"correct_bisp: sky_level= %.3f xphot=%.1f",sky_level,xphot);

/* Mean squared modulus : */
/* Recentre the frames: */
   TO_SINGLE(modsq,modsq1,&nx,&ny,&nx);
   RECENT_FFT(modsq1,modsq1,&nx,&ny,&nx);
   sprintf(filename,"%sm",out_prefix);
   JLP_WRITEIMAG(modsq1,&nx,&ny,&nx,filename,comments);

/* Bispectrum : */
   sprintf(filename,"%sb",out_prefix);
   JLP_D_WRITEIMAG(bispp,&bisp_dim,&ny_bisp,&bisp_dim,filename,comments);

return(0);
}
/**************************************************************** 
* Photon noise correction with Hege's method 
*
*  Input:
* modsq[]: mean normalized squared modulus
* phot_modsq[]: power spectrum of photon response  
* xphot: mean photon flux per frame 
* bispp[]: bispectrum (sorted in a "standard" way, and thus different
*         from what is needed by photon_corr...)
*
* Photon noise correction (cf JOSA 2, 14, Wirnitzer):
* (When spectrum normalized to one in the center)
* <i(u)>**2 = E(D...)/N**2 - 1/N
*
* <i(3)(u,v)> = E(D(3)(u,v)/N**3 - E(D(2)(u)/N**2)/N - E(D(2)(v)... +2/N**3)
*****************************************************************/
int photon_corr_mask_hege(double *bispp, double *modsq, double *phot_modsq,
                          INT4 *nx, INT4 *ny, INT4 *bisp_dim, float *xphot,
                          INT4 *nbeta, INT4 *ngamma)
{
float *xmodsq;
float xphot2, *phot_re;
double w2, work, work1, epsilon=1.e-8;
INT4 k1, k2, k3, iix, iiy, nb;
register INT4 i, ng, iklm;

xmodsq = (float *)malloc((*nbeta + 1) * sizeof(float));

/*******************************************************/
/* first division by the power spectrum phot_modsq: */
  for( i = 0; i < (*ny) * (*nx); i++) 
   if(phot_modsq[i] > 0.)
         modsq[i] = modsq[i] / (phot_modsq[i] + epsilon);
   else
         modsq[i] = 0.;

if((phot_re = (float*)malloc((*nbeta + 1) * sizeof(float))) == NULL)
 {
 printf("photon_corr_mask/fatal error, allocating memory space \n");
 free(xmodsq);
 exit(-1);
 };

/*******************************************************/
/* First the bispectrum, then the power spectrum */

  xphot2 = *xphot * *xphot;
/* IXY(1,NB) and IXY(2,NB) are the X,Y coordinates of the
* element NB of the spectral list.
* As they (i.e. IXY) can be negative, and that zero frequency is at (0,0),
* we do the following transformation:
*/
  for(nb = 0; nb <= *nbeta; nb++)
     {
/* Get coordinates (i,j) from nb index: */
       COVER_IXY(&iix, &iiy, &nb);
       iix = ((iix + *nx) % *nx);
       iiy = ((iiy + *ny) % *ny);
/* xmodsq is used to store the mean (not normalized) square modulus: */
       xmodsq[nb] = modsq[iix + iiy * *nx];
/* The photon response is assumed to be symmetric and real, hence its FT is real: */
       phot_re[nb] = sqrt(phot_modsq[iix + iiy * *nx]);
#ifdef DEBUG
       if(nb < 4) printf(" nb = %d iix, iiy: %d %d \n",nb,iix,iiy);
       if(nb < 4) printf(" xmodsq = %f \n",xmodsq[nb]);
#endif
     }

/* JLP96: */
/* Number of photons:  xphot**2 + xphot - modsq_0 = 0
* Delta = b**2 - 4 a*c = 1 + 4 * modsq_0
* Solutions: (- b +/- sqrt(Delta))/2a
*  i.e.,     (-1 + sqrt(1 + 4 * modsq_0) )/2
*/
  work = (-1. + sqrt((double)(1. + 4. * xmodsq[0])))/2.;
  printf(" My estimate of xphot (from modsq[0]) is: %f, whereas yours is: %f \n",
           work,*xphot);
/* End of JLP96. */

/****************************************************************/
/* Phase factor of the bispectrum (with bispectral list):
* and correction from photon noise effects:
*/
 w2 = 2. *  *xphot;
 iklm = 0;
for(ng = 0; ng < *ngamma; ng++)
 {
     cover_klm0(&k1,1,ng);
     cover_klm0(&k2,2,ng);
     cover_klm0(&k3,3,ng);

/* JLP2000: division by the bispectrum of the photon response: */
     work1 = phot_re[k1] * phot_re[k2] * phot_re[k3] + epsilon;
     bispp[ng] /= work1;
     bispp[ng + *bisp_dim] /= work1;

#ifdef DEBUG
     if(ng < 4) printf(" ng = %d k, l, m: %d %d %d \n",ng,k1,k2,k3);
     if(ng < 4) printf(" xmodsq1, xmodsq2, xmodsq3 %f %f %f \n",xmodsq[k1],xmodsq[k2],xmodsq[k3]);
#endif
/* Photon noise correction 
*/
     work1 = - (xmodsq[k1] + xmodsq[k2] + xmodsq[k3]) + w2;
#ifdef DEBUG
     if(ng < 4) printf(" ng= %d corr= %f \n",ng,work1);
#endif
     bispp[ng] += work1;
/* Phase factor: */
/* JLP96: Do not normalize it (to be able to visualize it with bisp_to_image)
     w1 = bispp[ng]*bispp[ng] + bispp[ng + *bisp_dim] * bispp[ng + *bisp_dim];
     if(w1 > 0) 
       {
       w1 = sqrt(w1);
       bispp[ng] /= w1;
       bispp[ng + *bisp_dim] /= w1;
       }
*/

  }
/*******************************************************/
/* Then correcting the squared modulus: */
  for( i = 0; i < *ny * *nx; i++)
     {
/* Photon noise correction: (biased_sq = xphot + xphot_sq * unbiased_sq)*/
      modsq[i] -= *xphot;
    }

free(phot_re);
free(xmodsq);
return(0);
}
