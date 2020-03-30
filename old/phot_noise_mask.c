/**************************************************************** 
  phot_noise_mask.c
 Correction for photon noise, using long exposure and flat field.

 When no mask is used, it is better to use photon_corr.for (...)

 JLP
 Version 04-06-99
*****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <jlp_ftoc.h>

main(argc, argv)
int argc;
char **argv;
{
/* Images : */
float *long_int, *modsq1, *bisp1, *ffield, *mask;
double *modsq, *bispp;
float xphot, fraction;
INT_PNTR pntr_ima;
INT4 nx_bisp, ny_bisp, nx, ny, nx_ff, ny_ff, nx_mask, ny_mask;
INT4 ir, nbeta, ngamma, bisp_dim, isize, max_nclosure;
register int i;
char *pc, in_ext[41], out_ext[41], buffer[81];
char filename[61], ff_name[61], ff_comments[81], comments[81];
char mask_name[61], mask_comments[81];

printf(" Program phot_noise_mask \n");
printf(" JLP Version 12-10-99 \n");
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
  printf(" runs phot_noise in_exten out_exten mask flat_field ir,max_nclo \n");
  printf(" or runs phot_noise in_exten out_exten mask 0 ir,max_nclo nphot_per_frame\n");
  printf(" or runs phot_noise in_exten out_exten 0 0 ir,max_nclo nphot_per_frame \n");
  exit(-1);
  }

/* Interactive input of in_exten and out_exten: */
if (argc == 1)
 {
   printf(" Input file extension := "); scanf("%s",in_ext);
/************* File extension for output files: */
   printf(" Output file extension := "); scanf("%s",out_ext);
   printf(" Fourier mask (0 if no mask) := "); scanf("%s",mask_name);
   printf(" Flat field (0 if unity flat field) := "); scanf("%s",ff_name);
   printf(" Radius of uv-coverage (IR) in pixels and max_nclosure: (12,1000 f.i.) \n");
   scanf("%d,%d",&ir,&max_nclosure);
 }
else
 {
  strcpy(in_ext,argv[1]);
  strcpy(out_ext,argv[2]);
  strcpy(mask_name,argv[3]);
  strcpy(ff_name,argv[4]);
  sscanf(argv[5],"%d,%d",&ir,&max_nclosure);
 }

/* End extension strings with zero: */
  pc = in_ext;
  while(*pc && *pc != ' ') pc++;
  *pc='\0';
  pc = out_ext;
  while(*pc && *pc != ' ') pc++;
  *pc='\0';

#ifdef DEBUG
printf(" DEBUG Version, will read >%s< files and write >%s< \n",in_ext);
#endif

/***************************************************************/
JLP_BEGIN();
JLP_INQUIFMT();


/*****************************************************************/
/* Computing the spectral and bispectral lists
   corresponding to the selected uv coverage: */
  if(mask_name[0] == '0')
    {
/* Full u-v coverage: */
    printf(" OK: full u-v coverage.\n");
/*
    COVERA(&ir,&nbeta,&ngamma);
*/
/* Free virtual memory space: */
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
    JLP_VM_READIMAG1(&pntr_ima,&nx_mask,&ny_mask,mask_name,mask_comments);
    mask = (float *)pntr_ima;
    }

/* "Masked" u-v coverage: */
    COVERA_MASK(mask,&nx_mask,&ny_mask,&ir,&max_nclosure,&nbeta,&ngamma);

  printf(" Radius of uv-coverage (IR) in pixels: %d\n",ir);
  printf(" nbeta: %d ngamma: %d\n",nbeta,ngamma);

/* Free virtual memory space: */
    free(mask);

/* JLP99: I add the possibility of using tt_l, tt_m, tt_b
instead of long_tt modsq_tt bisp1_tt */

/* Mean squared modulus : */
   if(in_ext[0] == '_')sprintf(filename,"modsq%s",in_ext);
   else sprintf(filename,"%sm",in_ext);
   printf(" Reading mean squared modulus: >%s< \n",filename);
   JLP_VM_READIMAG1(&pntr_ima,&nx,&ny,filename,comments);
   modsq1 = (float *) pntr_ima;

/* Mean bispectrum : */
   if(in_ext[0] == '_')sprintf(filename,"bisp1%s",in_ext);
   else sprintf(filename,"%sb",in_ext);
   printf(" Reading mean bispectrum: >%s< \n",filename);
   JLP_VM_READIMAG1(&pntr_ima,&nx_bisp,&ny_bisp,filename,comments);
   bisp1 = (float *) pntr_ima;
   if(nx_bisp != ngamma)
     {
      printf(" Fatal error: input bispectrum incompatible with u-v coverage!\n");
      exit(-1);
     }
   bisp_dim = ngamma;

/* Copy only the first two lines to double size array: */
   isize = bisp_dim * 2 * sizeof(double);
   bispp = (double *) malloc(isize);
   for(i = 0; i < 2*bisp_dim; i++) bispp[i] = (double)bisp1[i];

/******************************************************************/
/* Read flat field file: */
 if(ff_name[0] != '0')
  {
   printf(" Sorry this option (ffield) is not available yet \n"); exit(-1);
   JLP_VM_READIMAG1(&pntr_ima,&nx_ff,&ny_ff,ff_name,ff_comments);
   ffield = (float *) pntr_ima;

/* Long integration : */
   if(in_ext[0] == '_')sprintf(filename,"long%s",in_ext);
   else sprintf(filename,"%sl",in_ext);
   JLP_VM_READIMAG1(&pntr_ima,&nx,&ny,filename,comments);
   long_int = (float *) pntr_ima;

/* Free virtual memory space: */
   free(long_int);
   free(ffield);
  }
/********************* Unity file for ffield: *************/
  else
  {
   printf("\n OK, We assume unity ffield has been used\n");
   printf(" Enter number of photons/frame and fraction of modsq to be used:\n");
   fraction = 0;
   if(argc == 1)
     scanf("%f,%f",&xphot,&fraction);
   else
     sscanf(argv[6],"%f,%f",&xphot,&fraction);
  }

/****************************************************************/
/* Recentre the frames: (necessary for photon_corr_mask) */
   RECENT_FFT(modsq1,modsq1,&nx,&ny,&nx);
   isize = nx * ny * sizeof(double);
   modsq = (double *) malloc(isize);
   to_double(modsq1,modsq,&nx,&ny,&nx);

   printf(" Correction with xphot = %f and fraction = %f\n",xphot,fraction);
/* JLP2000: program to be updated ... */
   photon_corr_mask(bispp,modsq,phot_modsq,&nx,&ny,&bisp_dim,
                    &xphot,&nbeta,&ngamma,&fraction);

/****************************************************************/
/* Now output the results : */
   sprintf(comments," Photon correction: xphot= %f",xphot);

/* Mean squared modulus : */
/* Recentre the frames: */
   to_single(modsq,modsq1,&nx,&ny,&nx);
   RECENT_FFT(modsq1,modsq1,&nx,&ny,&nx);
   if(out_ext[0] == '_')sprintf(filename,"modsq%s",out_ext);
   else sprintf(filename,"%sm",out_ext);
   JLP_WRITEIMAG(modsq1,&nx,&ny,&nx,filename,comments);

/* Bispectrum : */
   if(out_ext[0] == '_')sprintf(filename,"bisp1%s",out_ext);
   else sprintf(filename,"%sb",out_ext);
   for(i = 0; i < 2*bisp_dim; i++) bisp1[i] = (float)bispp[i];
   JLP_WRITEIMAG(bisp1,&bisp_dim,&ny_bisp,&bisp_dim,filename,comments);

/* End : */
  JLP_END();
}
