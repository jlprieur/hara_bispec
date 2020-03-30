/**************************************************************** 
* phot_noise_1D
*
* From phot_noise_mask.c
* Correction for photon noise, using long exposure and flat field.
*
* When no mask is used, it is better to use photon_corr.for (...)
*
* JLP
* Version 21-01-96
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
float *modsq, *long_int;
float *bisp1, *ffield, *mask;
long int pntr_ima, nx_bisp, ny_bisp, nx, ny, nx_ff, ny_ff, nx_mask;
int ir, nbeta, ngamma, bisp_dim, isize, max_nclosure;
register int i;
char *pc, in_ext[41], out_ext[41], buffer[81];
char filename[61], ff_name[61], ff_comments[81], comments[81];
char mask_name[61], mask_comments[81];

printf(" Program phot_noise_1D \n");
printf(" JLP Version 06-02-96 \n");

/* One or three parameters only are allowed to run the program: */
/* Carefull: always 7 parameters, using JLP "runs" */
#ifdef DEBUG
  printf(" argc=%d\n",argc);
  printf(" argv[3]=>%s<\n",argv[3]);
#endif
if((argc > 1 && argc < 6) || (argc == 7 && !strcmp(argv[3],"")))
  {
  printf("        Fatal error: Wrong syntax, argc=%d\n",argc);
  printf(" Syntax is:  \n");
  printf(" runs phot_noise_1D in_exten out_exten mask flat_field ir,max_nclo \n");
  printf(" or runs phot_noise in_exten out_exten mask 0 ir,max_nclo\n");
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
      isize = nx_mask * sizeof(float);
      JLP_GVM(&mask,&isize);

/* Create a filled mask: */
      for(i = 0; i < nx_mask; i++) mask[i] = 1.;

    }
  else
/* Read Fourier mask file: */
    {
    printf(" Sorry this option (Fourier mask) is not allowed here! \n");
    exit(-1);
/*
    JLP_VM_READIMAG(&pntr_ima,&nx_mask,&ny_mask,mask_name,mask_comments);
    JLP_FROM_MADRID(&pntr_ima,&mask);
*/
    }

/* "Masked" u-v coverage: */
    COVERA_MASK_1D(mask,&nx_mask,&ir,&max_nclosure,&nbeta,&ngamma);

  printf(" Radius of uv-coverage (IR) in pixels: %d\n",ir);
  printf(" nbeta: %d ngamma: %d\n",nbeta,ngamma);

/* Free virtual memory space: */
    JLP_FVM(&mask);

/* Mean squared modulus : */
   sprintf(filename,"modsq%s",in_ext);
   printf(" Reading mean squared modulus: >%s< \n",filename);
   JLP_VM_READIMAG(&pntr_ima,&nx,&ny,filename,comments);
   JLP_FROM_MADRID(&pntr_ima,&modsq);

/* Mean bispectrum : */
   sprintf(filename,"bisp1%s",in_ext);
   printf(" Reading mean bispectrum: >%s< \n",filename);
   JLP_VM_READIMAG(&pntr_ima,&nx_bisp,&ny_bisp,filename,comments);
   JLP_FROM_MADRID(&pntr_ima,&bisp1);
   if(ny_bisp != ny)
     {
      printf(" ny_bisp = %d whereas ny_modsq = %d\n",ny_bisp,ny);
      printf(" Fatal error: input bispectrum incompatible square modulus!\n");
      exit(-1);
     }
   if(nx_bisp != 3*ngamma)
     {
      printf(" nx_bisp = %d whereas 3 * ngamma = %d\n",nx_bisp,3*ngamma);
      printf(" Fatal error: input bispectrum incompatible with u-v coverage!\n");
      exit(-1);
     }
   bisp_dim = 3*ngamma;

/******************************************************************/
/* Read flat field file: */
 if(ff_name[0] != '0')
   {printf(" Sorry this option (ffield) is not available yet \n"); exit(-1);}

 if(ff_name[0] != '0')
  {
   JLP_VM_READIMAG(&pntr_ima,&nx_ff,&ny_ff,ff_name,ff_comments);
   JLP_FROM_MADRID(&pntr_ima,&ffield);

/* Long integration : */
   sprintf(filename,"long%s",in_ext);
   JLP_VM_READIMAG(&pntr_ima,&nx,&ny,filename,comments);
   JLP_FROM_MADRID(&pntr_ima,&long_int);

/* Free virtual memory space: */
   JLP_FVM(&long_int);
   JLP_FVM(&ffield);
  }
/********************* Unity file for ffield: *************/
  else
  {
   printf("\n OK, We assume unity ffield has been used\n");
  }

/****************************************************************/
/* Recentre the frames: (necessary for photon_corr_1D) */
   RECENT_FFT_1D(modsq,modsq,&nx,&ny,&nx);

   photon_corr_1D(bisp1,modsq,&nx,&ny,&bisp_dim,&nbeta,&ngamma);

/****************************************************************/
/* Now output the results : */

/* Mean squared modulus : */
/* Recentre the frames: */
   RECENT_FFT_1D(modsq,modsq,&nx,&ny,&nx);
   sprintf(filename,"modsq%s",out_ext);
   JLP_WRITEIMAG(modsq,&nx,&ny,&nx,filename,comments);

/* Bispectrum : */
   sprintf(filename,"bisp1%s",out_ext);
   JLP_WRITEIMAG(bisp1,&nx_bisp,&ny_bisp,&nx_bisp,filename,comments);

/* End : */
  JLP_END();
}
