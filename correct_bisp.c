/**************************************************************** 
*  correct_bisp
*  From phot_noise_mask.c
*  Correction for sky background using long exposure.
*
* JLP
* Version 13-05-2009
*****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <jlp_ftoc.h>
#include "jlp_cover_mask.h" 

#define DEBUG

int main(int argc, char *argv[])
{
/* Images : */
float *long_int, *modsq, *bispp, *ffield, *phot_modsq, *mask;
float xphot, fraction, sky_level;
INT_PNTR pntr_ima;
INT4 nx_bisp, ny_bisp, nx, ny, nx1, ny1, nx_mask, ny_mask;
INT4 ir, nbeta, ngamma, bisp_dim, isize, max_nclosure;
int status;
register int i;
char *pc, in_ext[41], out_ext[41];
char filename[61], ff_name[61], comments[81];
char phot_modsq_name[61], mask_name[61];

printf(" Program correct_bisp \n");
printf(" JLP Version 13-05-2009 \n");
printf(" For the extension, either _tt or tt_\n");
printf("(i.e. either modsq_tt, bisp1_tt, ... or tt_m, tt_b, t_l\n");


/* One or three parameters only are allowed to run the program: */
/* Careful: 7 parameters always, using JLP "runs" */
#ifdef DEBUG
  printf(" argc=%d\n",argc);
  printf(" argv[3]=>%s<\n",argv[3]);
#endif
if((argc > 1 && argc < 5) || (argc == 7 && !strcmp(argv[4],"")))
  {
  printf("        Fatal error: Wrong syntax, argc=%d\n",argc);
  printf(" Syntax is:  \n");
  printf("runs correct_bisp in_exten out_exten phot_modsq flat_field ir,max_nclo,fraction\n");
  printf("or runs correct_bisp in_exten out_exten phot_modsq 0 ir,max_nclo,fraction\n");
  printf("or runs correct_bisp in_exten out_exten 0 0 ir,max_nclo,fraction \n");
  exit(-1);
  }

/* Interactive input of in_exten and out_exten: */
if (argc == 1)
 {
   printf(" Input file extension := "); scanf("%s",in_ext);
/************* File extension for output files: */
   printf(" Output file extension := "); scanf("%s",out_ext);
   printf(" Photon modsq (0 if not available) := "); scanf("%s",phot_modsq_name);
   printf(" Flat field (0 if unity flat field) := "); scanf("%s",ff_name);
   printf(" Radius of uv-coverage (IR) in pixels, max_nclosure, fraction: (12,1000,1. f.i.) \n");
   scanf("%d,%d,%f",&ir,&max_nclosure,&fraction);
 }
else
 {
  strcpy(in_ext,argv[1]);
  strcpy(out_ext,argv[2]);
  strcpy(phot_modsq_name,argv[3]);
  strcpy(ff_name,argv[4]);
  sscanf(argv[5],"%d,%d,%f",&ir,&max_nclosure,&fraction);
 }

/* End extension strings with zero: */
  pc = in_ext;
  while(*pc && *pc != ' ') pc++;
  *pc='\0';
  pc = out_ext;
  while(*pc && *pc != ' ') pc++;
  *pc='\0';

#ifdef DEBUG
printf(" DEBUG Version, will read >%s< files and write >%s< \n",in_ext,out_ext);
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
/* Output spectral and bispectral lists 
*/
#ifdef DEBUG
  output_lists_coverage(&nbeta,&ngamma);
#endif

/* Free virtual memory space: */
    free(mask);

/* JLP99: I add the possibility of using tt_l, tt_m, tt_b
instead of long_tt modsq_tt bisp1_tt */

/* Mean squared modulus : */
   if(in_ext[0] == '_')sprintf(filename,"modsq%s",in_ext);
   else sprintf(filename,"%sm",in_ext);
   printf(" Reading mean squared modulus: >%s< \n",filename);
   JLP_VM_READIMAG1(&pntr_ima,&nx,&ny,filename,comments);
   modsq = (float *) pntr_ima;

/* Mean bispectrum : */
   if(in_ext[0] == '_')sprintf(filename,"bisp1%s",in_ext);
   else sprintf(filename,"%sb",in_ext);
   printf(" Reading mean bispectrum: >%s< \n",filename);
   JLP_VM_READIMAG1(&pntr_ima,&nx_bisp,&ny_bisp,filename,comments);
   bispp = (float *) pntr_ima;
   if(nx_bisp != ngamma)
    {
     printf(" Fatal error: input bispectrum incompatible with u-v coverage!\n");
     exit(-1);
    }
   bisp_dim = ngamma;

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

#if 0
/* Old version: */
/* Compute the sky value from the long integration : */
   if(in_ext[0] == '_')sprintf(filename,"long%s",in_ext);
   else sprintf(filename,"%sl",in_ext);
   JLP_VM_READIMAG1(&pntr_ima,&nx,&ny,filename,comments);
   long_int = (float *) pntr_ima;
   auto_sky(long_int,nx,ny,nx,&sky_level);
   free(long_int);
#endif

/* Photon modsq (response to a photo-event: */
   JLP_VM_READIMAG1(&pntr_ima,&nx1,&ny1,phot_modsq_name,comments);
   if(nx1 != nx || ny1 != ny) 
         {printf("uncompatible size for %s \n",phot_modsq_name);
          exit(-1);}
   phot_modsq = (float *) pntr_ima;

/* Correction of central value of modsq because of a non null background */
printf("JLPPP/DDEBUG: NO CORRECTION HERE, FOR DEBUG.... \n");
/* WARNING: should recenter the frames for this old version:
* and double precision arrays!
   corr_bisp_sky(bispp,modsq,&nx,&ny,&bisp_dim,&sky_level,&nbeta,&ngamma);
*/

   printf(" Gain correction to be applied to phot_modsq = %f\n",fraction);

/* Correction of photon bias: (OLD METHOD) */
#if 0
/* Recentre the frames: (necessary for photon_corr_mask) */
   RECENT_FFT_DOUBLE(modsq,modsq,&nx,&ny,&nx);
/* Before 2008: */
   photon_corr_mask(bispp,modsq,phot_modsq,&nx,&ny,&bisp_dim,
                    &xphot,&nbeta,&ngamma,&fraction);
/* Recentre the frames: */
   RECENT_FFT_DOUBLE(modsq,modsq,&nx,&ny,&nx);
   sprintf(comments,"correct_bisp: sky_level= %.3f xphot=%.1f",sky_level,xphot);
#else
/* After 2008: new JLP's automatic method */
   status = photon_corr_auto1(bispp,modsq,phot_modsq,&nx,&ny,&bisp_dim,
                              &xphot,&nbeta,&ngamma);
   sprintf(comments,"correct_bisp: xphot=%.1f status=%d", xphot, status);
#endif
/****************************************************************/
/* Now output the results : */

   if(!status) {
/* Mean squared modulus : */
   if(out_ext[0] == '_')sprintf(filename,"modsq%s",out_ext);
   else sprintf(filename,"%sm",out_ext);
   JLP_WRITEIMAG(modsq,&nx,&ny,&nx,filename,comments);

/* Bispectrum : */
   if(out_ext[0] == '_')sprintf(filename,"bisp1%s",out_ext);
   else sprintf(filename,"%sb",out_ext);
   JLP_WRITEIMAG(bispp,&bisp_dim,&ny_bisp,&bisp_dim,filename,comments);

/* Estimated photon response : */
   if(out_ext[0] == '_')sprintf(filename,"photon_resp%s",out_ext);
   else sprintf(filename,"%sphoton_resp",out_ext);
   JLP_WRITEIMAG(phot_modsq,&nx,&ny,&nx,filename,comments);
   }

return(status);
}
