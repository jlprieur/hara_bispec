/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* bisp_from_image (from decode_ima.c, version 13/10/1999)
* To compute power spectrum and bispectrum of *.PHIN 
*
* JLP
* Version 13-04-2008
---------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <jlp_ftoc.h>

#define DEBUG 

/* Defined in "decode_set.c" */
int compute_uv_coverage(FILE *fp1, char *fmask_name,
                        INT4 ir, INT4 *nbeta, INT4 *ngamma, INT4 max_nclosure);
int prepare_output_bisp(double *bispp, double *modsq, double *snrm,
                        INT4 nx, INT4 ny, float xframes, float xphot,
                        INT4 nbeta, INT4 ngamma, INT4 photon_correction,
                        FILE *fp1);

/* Defined in "jlp_cover_mask.c" */ 
int bispec3(double *re, double *im, double *modsq, double *snrm,
            INT4 *nx, INT4 *ny, double *bispp, INT4 *ir,
            INT4 *nbeta, INT4 *ngamma);

int main(int argc, char *argv[])
{
/* Images nx*ny, irmax=25, ngmax=187566 : */

double *image2, *modsq, *snrm, *bispp, *mim;
float *image1, xframes, xphot;
INT_PNTR pntr_ima;
INT4 max_nclosure, descr_length, istatus;
INT4  isize, ir, nbeta, ngamma, photon_correction;
INT4  nx, ny, nxy, ny_bisp;
register int i;
char infile[61], comments[81], buffer[81];
char descr_name[61], descr_value[81], logfile[81], *pc;
char output_prefix[61], outfile[41], outcomments[81], fmask_name[61];
FILE *fp1;

  printf(" Program bisp_from_image (to process *.PHIN).  Version of 13-04-2008\n");
/**************************************************************/
/* Opening logfile: */
sprintf(logfile,"decode_im%s.log",output_prefix);
if((fp1 = fopen(logfile,"w")) == NULL)
   {
   printf("decode_ima/Fatal error opening logfile >%s< \n",logfile);
   exit(-1);
   }
  fprintf(fp1," Program bisp_from_image  Version 13-04-2008 \n");
  fprintf(fp1," To compute the bispectrum of one image only\n");


/* Input parameters: */
if (argc == 7 && argv[3]) argc = 4;
if (argc != 4)
  {
  printf("argc = %d \n",argc);
  printf("\nUSAGE:\n");
  printf("decode_ima input_image Radius_uv_coverage,nclosure output_file_prefix \n");
  printf(" Example1: decode_ima ttv2.PHIS 12,100 ads10\n ");
  printf(" Example2: decode_ima test.PHIS 15,200,mask_uv st\n ");
  printf("\n(Since possibility of adding uv-mask file name) \n");
  exit(-1);
  }

/* Generic name for the data image files: */
   strcpy(infile,argv[1]);

/* Radius of uv-coverage: */
/* Syntax: ir,max_nclosure
 or:       ir,max_nclosure,fmask_name 
*/ 
  *fmask_name = '\0';
  strcpy(buffer,argv[2]);
  sscanf(buffer,"%d,%d",&ir,&max_nclosure);
/* Look for possible presence of mask name: */
  pc = buffer;
  while(*pc && *pc != ',') pc++;
  if(*pc == ',') pc++;
  while(*pc && *pc != ',') pc++;
  if(*pc == ',') {pc++; strcpy(fmask_name, pc);}
  if(*fmask_name) printf(" Mask of uv-coverage : %s\n",fmask_name);
  else printf(" Full pupil was used: no mask\n");

/* Output file prefix: */
  strcpy(output_prefix,argv[3]);
  pc = output_prefix;
  while(*pc && *pc != ' ') pc++;
  *pc='\0';

printf("OK: input image: >%s<   output prefix: >%s<\n", infile, output_prefix);
printf("OK: ir=%d max_nclosure=%d fmask_name >%s<\n", 
       ir, max_nclosure, fmask_name);
fprintf(fp1,"Input image: >%s< output prefix: >%s<\n", infile, output_prefix);
fprintf(fp1,"Parameters: ir=%d max_nclosure=%d fmask_name >%s<\n", 
        ir, max_nclosure, fmask_name);

/*****************************************************************/
/* Inquire format of the images */
JLP_INQUIFMT();

/***************************************************************/
/* Reading the data file */
  istatus = JLP_VM_READIMAG1(&pntr_ima,&nx,&ny,infile,comments);
  if(istatus)
     {
     printf(" Fatal error reading image %s\n",infile);
     exit(-1);
     }
  image1 = (float *)pntr_ima;

/* Allocation of memory: */
 nxy = nx * ny;
 isize = nxy * sizeof(double);
 modsq = (double *) malloc(isize);
 snrm = (double *) malloc(isize);
 image2 = (double *) malloc(isize);
 mim = (double *) malloc(isize);
 if(modsq == NULL || snrm == NULL || image2 == NULL
    || mim == NULL)
 {
 printf(" Fatal error allocating memory space (modsq, ...) \n");
 exit(-1);
 }

/* Transfer from float image1 to double image2: */
   for(i = 0; i < nxy; i++) image2[i] = image1[i];

/* Erasing the arrays : */
    for(i = 0; i < nxy; i++) 
      {
      modsq[i] = 0.;
      snrm[i] = 0.;
      mim[i] = 0.;
      }

/* Computing the uv coverage: */
  printf(" Radius of uv-coverage (IR) in pixels: %d\n",ir);
  compute_uv_coverage(fp1,fmask_name,ir,&nbeta,&ngamma,max_nclosure);

/* Allocating memory space for the bispectrum: */
  isize = 4 * ngamma * sizeof(double);
  bispp = (double*) malloc(isize);
  for(i=0; i < 4 * ngamma; i++) { bispp[i]=0.;}

/* Resetting the imaginary part for the FFT: */
  for(i = 0; i < nxy; i++) mim[i] = 0.;

/* Fourrier Transform: */
/*
  FFT_2D_FAST(image2,mim,&nx,&ny);
*/
  fftw_double(image2, mim, (int)nx, (int)ny,1);

/* Processing this image now:  bispec1 is with photon noise correction
bispec3 is without*/
  /*
   bispec1(image2,mim,modsq,snrm,&nx,&ny,&nx,bispp,&ir,&nbeta,&ngamma);
   */
/* bispec3 in "jlp_cover_mask.c" */
   bispec3(image2,mim,modsq,snrm,&nx,&ny,bispp,&ir,&nbeta,&ngamma);

/* Compute snrm, SNR of bispectrum, etc:  (xphot is not used) */
   xframes = 1.; xphot = 1.; photon_correction = 0;

#ifdef DEBUG
  printf(" Now sequence is: re[i],im[i],sumsq_re[i],sumsq_im[i],re[i+1],im[i+1]...\n");
  for(i = 0; i < 5; i++)
      printf("bispp[%d]: (%e;%e)\n",i,bispp[4*i],bispp[4*i+1]);
#endif

   prepare_output_bisp(bispp,modsq,snrm,nx,ny,xframes,xphot,
                         nbeta,ngamma,photon_correction,fp1);

#ifdef DEBUG  
  printf(" New sequence is: re[i],re[i+1]...im[i],im[i+1]...sumsq_re[i],sumsq_re[i+1]...sumsq_im[i],sumsq_im[i+1]\n");
  for(i = 0; i < 5; i++)
      printf("bispp[%d]: (%e;%e)\n",i,bispp[i],bispp[i+ngamma]);
#endif


/*********************************************************/
/* Now output of the results : */

/* Comments: */
   sprintf(outcomments,"From %s", infile);

/* Leave blanks at the end, for Fortran interface: */
   sprintf(descr_name,"IR  ");
   sprintf(descr_value,"%d  ",ir);
   descr_length = 4;
   JLP_WDESCR(descr_name,descr_value,&descr_length,&istatus);
/* Leave blanks at the end, for Fortran interface: */
   sprintf(descr_name,"MAX_NCLO  ");
   sprintf(descr_value,"%d  ",max_nclosure);
   descr_length = 8;
   JLP_WDESCR(descr_name,descr_value,&descr_length,&istatus);

/* Long integration : */
   sprintf(outfile,"%s_l",output_prefix);
   JLP_WRITEIMAG(image1,&nx,&ny,&nx,outfile,outcomments);

/* Mean squared modulus : */
   RECENT_FFT_DOUBLE(modsq,modsq,&nx,&ny,&nx);
   sprintf(outfile,"%s_m",output_prefix);
   JLP_D_WRITEIMAG(modsq,&nx,&ny,&nx,outfile,outcomments);

/* SNR of squared modulus (actually 1/sigma if no photon correction): */
   RECENT_FFT_DOUBLE(snrm,snrm,&nx,&ny,&nx);
   sprintf(outfile,"%s_snrm",output_prefix);
   JLP_D_WRITEIMAG(snrm,&nx,&ny,&nx,outfile,outcomments);

/* Bispectrum : */
   ny_bisp = 3;
   printf(" ngamma = %d, ny = %d \n", ngamma, ny);
   sprintf(outfile,"%s_b",output_prefix);
   JLP_D_WRITEIMAG(bispp,&ngamma,&ny_bisp,&ngamma,outfile,outcomments);

  printf(" Successful end \n");

/* End : */
  free(mim);

fclose(fp1);

return(0);
}
