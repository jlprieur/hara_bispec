/**************************************************************** 
* decode_set.c
* Set of routines used
* to decode data from the CAR, CP40 or MAMA photon counting cameras 
* and compute bispectrum.
*
* Contains:
* void swap_int(i)
* void swap_lint(i)
* void inv_ulong_int(i)
* int create_ffield(ffield,in_ffield,apodi,nx_ff,ny_ff,mini_value)
* int jlp_hamming(apodi,nx,ny,idim)
* int jlp_hamming1(apodi,nx,ny,idim)
* int jlp_hamming2(apodi,nx,ny,idim)
* int jlp_blackman(apodi,nx,ny,idim)
* int compute_uv_coverage(fp1,fmask_name,ir,nbeta,ngamma,max_nclosure)
* int compute_uv_coverage_1D(fp1,fmask_name,ir,nbeta,ngamma,max_nclosure)
* int compute_ffield(ffield,ff_name, nx_ff, ny_ff,  
*                    ixstart, iystart, nx, ny_ima, file_ext, fp1)
* int prepare_output_bisp(bispp, modsq, snrm, nx, ny, xframes, xphot, 
*                        nbeta, ngamma, photon_correction,fp1)
* int prepare_output_bisp_1D(bispp, modsq, snrm, nx, ny, xframes, xphot, 
*                        nbeta, ngamma, photon_correction,fp1)
* int output_bisp(bispp, bisp1, ngamma)
* int output_bisp_1D(bispp, bisp1, ngamma, ny)
* int julian(aa,mm,idd,time,djul)
* int inv_julian(date_in_years,aa,mm,idd,time,djul)
*
* JLP
* Version 25-02-97
*****************************************************************/
/*
#define DEBUG
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <jlp_ftoc.h>

#ifndef linux
#define PI 3.14159
#endif

/*
#define APODIZATION
*/
#ifdef APODIZATION
/* BLACKMAN or HAMMING */
/* When HAMMING is defined, uses Hamming instead of Blackman window */
#define HAMMING
#endif

/* Defined here: */
void swap_int(unsigned short *i);
void swap_lint(unsigned long *i);
void inv_ulong_int(unsigned long *i);
int create_ffield(float *ffield, float *in_ffield, float *apodi,
                  INT4 nx_ff, INT4 ny_ff, INT4 ixstart, INT4 iystart,
                  INT4 ixend, INT4 iyend, float mini_value);
int jlp_hamming(float *apodi, INT4 nx, INT4 ny, INT4 idim, 
                INT4 ixcent, INT4 iycent, INT4 xwidth, INT4 yheight);
int jlp_hamming1(float *apodi, INT4 nx, INT4 ny, INT4 idim);
int jlp_hamming2(float *apodi, INT4 nx, INT4 ny, INT4 idim);
int jlp_blackman(float *apodi, INT4 nx, INT4 ny, INT4 idim);
int compute_uv_coverage(FILE *fp1, char *fmask_name,
                        INT4 ir, INT4 *nbeta, INT4 *ngamma, INT4 max_nclosure);
int compute_uv_coverage_1D(FILE *fp1, char *fmask_name,
                        INT4 ir, INT4 *nbeta, INT4 *ngamma, INT4 max_nclosure);
int compute_ffield(float **ffield, char *ff_name, INT4 *nx_ff, INT4 *ny_ff, 
                   INT4 ixcent, INT4 iycent, INT4 nx, INT4 ny, 
                   char *file_ext, FILE *fp1, float mini_value);
int prepare_output_bisp(double *bispp, double *modsq, double *snrm, 
                        INT4 nx, INT4 ny, float xframes, float xphot, 
                        INT4 nbeta, INT4 ngamma, INT4 photon_correction, 
                        FILE *fp1);
int prepare_output_bisp_1D(double *bispp, double *modsq, double *snrm, 
                           INT4 nx, INT4 ny, float xframes, float xphot, 
                           INT4 nbeta, INT4 ngamma, INT4 photon_correction,
                           FILE *fp1);
int test_machine();
int bin_long_disp(long a);
int julian(double aa, INT4 mm, INT4 idd, double time, double *djul);
int inv_julian(double *date_in_years, double *aa, int *mm, int *idd,
               double time, double djul);
int output_bisp(double *bispp, float *bisp1, int ngamma);
int output_bisp_1D(double *bispp, float *bisp1, int ngamma, int nx);

/**************************************************************
* Swap two bytes of an unsigned short integer
***************************************************************/
void swap_int(unsigned short *i)
{
union {
unsigned short ii;
char         ch[2];
      } tmp;
char ch0; 

tmp.ii = *i;

#ifdef DEBUG
printf(" swap_int/before: *i= %d ch[O,1] %d %d \n",*i,tmp.ch[0],tmp.ch[1]);
printf(" swap_int/before: tmp.ii= %d  \n",tmp.ii);
#endif

ch0 = tmp.ch[0]; 
tmp.ch[0] = tmp.ch[1]; 
tmp.ch[1] = ch0; 
*i = tmp.ii;

#ifdef DEBUG
printf(" swap_int/after: *i= %d ch[O,1] %d %d \n",*i,tmp.ch[0],tmp.ch[1]);
#endif
}
/**************************************************************
* Swap two halves of a lont integer
* 1234 => 3412
***************************************************************/
void swap_lint(unsigned long *i)
{
union {
unsigned long ii;
short         ch[2];
      } tmp;
short ch0; 

tmp.ii = *i;

ch0 = tmp.ch[0]; 
tmp.ch[0] = tmp.ch[1]; 
tmp.ch[1] = ch0; 
*i = tmp.ii;

}
/**************************************************************
* Inversion of an unsigned long integer
* From [0 1 2 3] to [3 2 1 0]
*
***************************************************************/
void inv_ulong_int(unsigned long *i)
{
union {
unsigned long ii;
char ch[4];
      } tmp;
char ch0, ch1; 

tmp.ii = *i;
ch0 = tmp.ch[0]; 
ch1 = tmp.ch[1]; 
tmp.ch[0] = tmp.ch[3]; 
tmp.ch[1] = tmp.ch[2]; 
tmp.ch[2] = ch1; 
tmp.ch[3] = ch0; 
*i = tmp.ii;
}
/****************************************************************
* create_ffield
*
* From a mean flat field we create an array with the
* correct size and starting pixels which fit to the data. 
*
* Called by "compute_ffield"
*
* INPUT:
* in_ffield[nx_ffield,ny_ffield]
* ixstart,iystart,ixend,iyend: coordinates of working window 
*
* OUTPUT:
* ffield[nx_ffield,ny_ffield]
*
****************************************************************/
int create_ffield(float *ffield, float *in_ffield, float *apodi,
                  INT4 nx_ff, INT4 ny_ff, INT4 ixstart, INT4 iystart,
                  INT4 ixend, INT4 iyend, float mini_value)
{
register int i, j;
INT4 nvalues, nxy;
double sum;
float ff_mean;

/* Transfer: */
nxy = nx_ff * ny_ff;
for(j = 0; j < nxy; j++) ffield[j] = in_ffield[j];

/* First normalization: (take into account only non zero values) */
sum = 0.; nvalues=0;
for(j = iystart; j < iyend; j++)
  {
  for(i = ixstart; i < ixend; i++)
    {if(ffield[i + j * nx_ff] > 0)
       {
       sum += ffield[i + j * nx_ff];
       nvalues++;
       }
    }
  }
if(sum == 0)
  {
  printf("create_ffield/Fatal error, total sum of input flat field frame is zero!\n");
  exit(-1);
  }
ff_mean = sum / (float)(nvalues);
printf("create_ffield/Flat field raw average value is %f\n",ff_mean);

for(j = 0; j < nxy; j++) ffield[j] /= ff_mean;
 
/* Inversion of flat field: */
/* For the CAR, when the relative value is under mini_value=0.6, pb even with
   apodization. As this means that we are anyway outside of the 
   sensitive zone, we set the maximum to that value (or 1/value here): */ 
for(j = 0; j < nxy; j++) 
     {if(ffield[j] < mini_value) 
         ffield[j] = 1./mini_value;
     else 
         ffield[j] = 1./ffield[j];
     }

/* Apodization: */
#ifdef APODIZATION
for(j = 0; j < nxy; j++) ffield[j] *= apodi[j];
#endif

/* Second normalization: (take into account only values close to mean)*/
sum = 0.; nvalues=0;
for(j = iystart; j < iyend; j++)
  {
  for(i = ixstart; i < ixend; i++)
    {if(ffield[i + j * nx_ff] < 2. && ffield[i + j * nx_ff] > 0.5)
       {
       sum += ffield[i + j * nx_ff];
       nvalues++;
       }
    }
  }
ff_mean = sum / (float)(nvalues);
for(j = 0; j < nxy; j++) ffield[j] /= ff_mean;

#ifdef APODIZATION
printf("create_ffield/Inverse flat field average value (after apodization) is %f\n",
       ff_mean);
#else
printf("create_ffield/Inverse flat field average value (no apodization) is %f\n",
       ff_mean);
#endif

return(0);
}
/***************************************************************
* Creates a pseudo-Hamming apodization filter
* (Filled with ones in the middle): 
*
* xwidth, ywidth: size of working window
* xcent, ycent: coordinates of the center of the working window
**************************************************************/
int jlp_hamming(float *apodi, INT4 nx, INT4 ny, INT4 idim, 
                INT4 ixcent, INT4 iycent, INT4 xwidth, INT4 yheight)
{
double argx, w1, radius, radmin, radmax, edge_width, mean_width;
INT4 jrad;
register int i, j;

/* We take a subwindow to fit better the shape of the
   sensitive part of the camera: */
for(j = 0; j < ny; j++)
   for(i = 0; i < nx; i++)
    apodi[i + j * idim] = 1.;

/* Previously:
width = (double)nx / 6.;
radmin = (double)(nx/2) - width;
radmax = (double)(nx/2);
w1 = PI / (double)width;
*/
mean_width = (xwidth + yheight ) / 2.;
edge_width = mean_width / 6.;
radmin = mean_width / 2. - edge_width;
radmax = mean_width / 2.;
w1 = PI / (double)edge_width;

/* Works with circular filter: */
for(j = 0; j < ny; j++)
 {
   jrad = (j - iycent) * (j - iycent);
   for(i = 0; i < nx; i++) 
    {
     radius = jrad + (i - ixcent) * (i - ixcent);
     radius = sqrt(radius);
      if(radius > radmin && radius < radmax)
      {
       argx = w1 * (radius - radmin);
       apodi[i + j * idim] *= (0.54 + 0.46 * cos(argx));
      }
      else if(radius >= radmax)
      {
       apodi[i + j * idim] = 0.;
      }
    }
 }

return(0);
}
/***************************************************************
* Creates Hamming apodization filter on a subwindow 
* (when using full size (1024x1024)): 
**************************************************************/
int jlp_hamming1(float *apodi, INT4 nx, INT4 ny, INT4 idim)
{
double argx, argy, wx1, wy1;
float apodi_y;
INT4 istart, jstart, nx1, ny1;
register int i, j;

/* We take a subwindow to fit better the shape of the
   sensitive part of the camera: */
for(j = 0; j < ny; j++)
   for(i = 0; i < nx; i++)
    apodi[i + j * idim] = 0.;

ny1 = (ny * 5) / 6;
nx1 = (nx * 5) / 6;
jstart = (ny - ny1) / 2;
istart = (nx - nx1) / 2;
wy1 = 2. * PI / (double)ny1;
wx1 = 2. * PI / (double)nx1;
for(j = jstart; j < jstart + ny1; j++)
 {
   argy = wy1 * (double)(j - ny/2);
   apodi_y = 0.54 + 0.46 * cos(argy);
   for(i = istart; i < istart + nx1; i++)
    {
    argx = wx1 * (double)(i - nx/2);
    apodi[i + j * idim] = apodi_y * (0.54 + 0.46 * cos(argx));
    }
 }

return(0);
}
/***************************************************************
* Creates Hamming apodization filter on the whole frame:
**************************************************************/
int jlp_hamming2(float *apodi, INT4 nx, INT4 ny, INT4 idim)
{
double argx, argy, wx1, wy1;
float apodi_y;
register int i, j;

wy1 = 2. * PI / (double)ny;
wx1 = 2. * PI / (double)nx;
for(j = 0; j < ny; j++)
 {
   argy = wy1 * (double)(j - ny/2);
   apodi_y = 0.54 + 0.46 * cos(argy);
   for(i = 0; i < nx; i++)
    {
    argx = wx1 * (double)(i - nx/2);
    apodi[i + j * idim] = apodi_y * (0.54 + 0.46 * cos(argx));
    }
 }

return(0);
}
/***************************************************************
* Creates Blackman apodization file: 
****************************************************************/
int jlp_blackman(float *apodi, INT4 nx, INT4 ny, INT4 idim)
{
double argx, argy;
double wx1, wy1;
float apodi_y;
register int i, j;
INT4 istart, jstart, nx1, ny1;

printf(" nx=%d, ny=%d \n",nx,ny);

/* Normally, should be: */
/*
for(j = 0; j < ny; j++)
 {
   argy = 2. * PI * (double)(j - ny/2)/(double)ny;
   apodi_y = 0.42 + 0.5 * cos(argy) + 0.08 * cos(2. * argy);
   for(i = 0; i < nx; i++)
    {
    argx = 2. * PI * (double)(i - nx/2)/(double)nx;
    apodi[i + j * idim] = apodi_y * 
                      (0.42 + 0.5 * cos(argx) + 0.08 * cos(2. * argx));
    }
 }
*/

/* But we take a subwindow to fit better the shape of the 
   sensitive part of the camera: */
for(j = 0; j < ny; j++)
   for(i = 0; i < nx; i++)
    apodi[i + j * idim] = 0.; 

ny1 = (ny * 5) / 6;
nx1 = (nx * 5) / 6;
jstart = (ny - ny1) / 2;
istart = (nx - nx1) / 2;
wy1 = 2. * PI / (double)ny1;
wx1 = 2. * PI / (double)nx1;
for(j = jstart; j < jstart + ny1; j++)
 {
   argy = wy1 * (double)(j - ny/2);
   apodi_y = 0.42 + 0.5 * cos(argy) + 0.08 * cos(2. * argy);
   for(i = istart; i < istart + nx1; i++)
    {
    argx = wx1 * (double)(i - nx/2);
    apodi[i + j * idim] = apodi_y * 
                      (0.42 + 0.5 * cos(argx) + 0.08 * cos(2. * argx));
    }
 }

return(0);
}
/****************************************************************
*  compute_uv_coverage
*****************************************************************/
int compute_uv_coverage(FILE *fp1, char *fmask_name,
                        INT4 ir, INT4 *nbeta, INT4 *ngamma, INT4 max_nclosure)
{
 float *fmask;
 INT_PNTR pntr_ima;
 char fmask_comments[81];
 INT4 isize, nx_fmask, ny_fmask;
 register int i;

  if(fmask_name[0] == '\0')
    {
/* Old call: (either in jlp_cover1.for or in jlp_cover2.c)
    COVERA(&ir,nbeta,ngamma);
*/
/* Full u-v coverage: */
    printf(" Full u-v coverage.\n");
    fprintf(fp1," Full u-v coverage.\n");
    nx_fmask = 2 * ir + 1; ny_fmask = nx_fmask;
    isize = nx_fmask * ny_fmask * sizeof(float);
    fmask = (float*) malloc(isize);
    for(i = 0; i < nx_fmask * ny_fmask; i++) fmask[i] = 1.; 
    }
  else
/* Read Fourier mask file: */
    {
    JLP_VM_READIMAG1(&pntr_ima,&nx_fmask,&ny_fmask,fmask_name,fmask_comments);
    fmask = (float *)pntr_ima;
    if( nx_fmask < (2 * ir + 1) || ny_fmask < (2 * ir + 1))
      {
      fprintf(fp1," Fatal error/ Input mask is too small relative to ir \n");
      fprintf(fp1,"     (ir = %d, nx_fmask = %d, ny_fmask = %d)\n",
              ir, nx_fmask, ny_fmask);
      printf(" Fatal error/ Input mask is too small relative to ir \n");
      printf("     (ir = %d, nx_fmask = %d, ny_fmask = %d)\n",
              ir, nx_fmask, ny_fmask);
      exit(-1);
      }
    fprintf(fp1," Masked u-v coverage. Mask= %s\n",fmask_name);
    }

/* "Masked" u-v coverage: */
    COVERA_MASK(fmask,&nx_fmask,&ny_fmask,&ir,&max_nclosure,nbeta,ngamma);
/* Free virtual memory space: */
    free(fmask);

  fprintf(fp1," Radius of uv-coverage (IR) in pixels: %d\n",ir);
  fprintf(fp1," nbeta: %d ngamma: %d\n",*nbeta,*ngamma);
#ifdef DEBUG
  printf(" nbeta: %d ngamma: %d\n",*nbeta,*ngamma);
#endif

return(0);
}
/****************************************************************
*  compute_uv_coverage_1D
*****************************************************************/
int compute_uv_coverage_1D(FILE *fp1, char *fmask_name,
                        INT4 ir, INT4 *nbeta, INT4 *ngamma, INT4 max_nclosure)
{
 float *fmask;
 INT_PNTR pntr_ima;
 char fmask_comments[81];
 INT4 isize, nx_fmask, ny_fmask;
 register int i;

  if(fmask_name[0] == '\0')
    {
/* Full u-v coverage: */
    printf(" Full u-v coverage.\n");
    fprintf(fp1," Full u-v coverage.\n");
    nx_fmask = 2 * ir + 1; 
    isize = nx_fmask * sizeof(float);
    fmask = (float*) malloc(isize);
    for(i = 0; i < nx_fmask; i++) fmask[i] = 1.; 
    }
  else
/* Read Fourier mask file: */
    {
    JLP_VM_READIMAG1(&pntr_ima,&nx_fmask,&ny_fmask,fmask_name,fmask_comments);
    fmask = (float *)pntr_ima;
    if( nx_fmask < (2 * ir + 1) || ny_fmask != 1)
      {
      fprintf(fp1," Fatal error/ Input mask is too small relative to ir (or ny not equal to one)\n");
      fprintf(fp1,"     (ir = %d, nx_fmask = %d, ny_fmask = %d)\n",
              ir, nx_fmask, ny_fmask);
      printf(" Fatal error/ Input mask is too small relative to ir (or ny not equal to one)\n");
      printf("     (ir = %d, nx_fmask = %d, ny_fmask = %d)\n",
              ir, nx_fmask, ny_fmask);
      exit(-1);
      }
    fprintf(fp1," Masked u-v coverage. Mask= %s\n",fmask_name);
    }

/* "Masked" u-v coverage: */
    COVERA_MASK_1D(fmask,&nx_fmask,&ir,&max_nclosure,nbeta,ngamma);
/* Free virtual memory space: */
    free(fmask);

  fprintf(fp1," Radius of uv-coverage (IR) in pixels: %d\n",ir);
  fprintf(fp1," nbeta: %d ngamma: %d\n",*nbeta,*ngamma);
#ifdef DEBUG
  printf(" nbeta: %d ngamma: %d\n",*nbeta,*ngamma);
#endif

return(0);
}
/****************************************************************
* Flat field 
*
* JLP96: I use the routine (in decode_cp40) as follows: 
*     nx_ff= x size of ffield (full detector format)
*     ny_ff= y size of ffield (full detector format)
*     nx, ny= size of working window
****************************************************************/
int compute_ffield(float **ffield, char *ff_name, INT4 *nx_ff, INT4 *ny_ff, 
                   INT4 ixcent, INT4 iycent, INT4 nx, INT4 ny, 
                   char *file_ext, FILE *fp1, float mini_value)
{
char ff_comments[81], buffer[81]; 
float *in_ffield, *apodi;
INT_PNTR pntr_ima;
INT4 isize, ixstart, iystart, ixend, iyend;
register int i;

/* Boundaries of working window: */
ixstart = ixcent - nx / 2;
iystart = iycent - ny / 2;
ixend = ixstart + nx;
iyend = iystart + ny;

if(ixstart < 0 || ixend > *nx_ff) 
 {printf("create_ffield/Fatal error: ixstart=%d, ixend=%d xwidth=%d nx_ff=%d\n",
         ixstart,ixend,nx,*nx_ff);
  exit(-1);
 }
if(iystart < 0 || iyend > *ny_ff) 
 {printf("create_ffield/Fatal error: iystart=%d, iyend=%d ywidth=%d ny_ff=%d\n",
          iystart,iyend,ny,*ny_ff);
  exit(-1);
 }

/****************************************************************/
if(ff_name[0] != '0')
  {
/* Displays message if apodization: */
#ifdef APODIZATION
     buffer[0] = '\0';
#else
     strcpy(buffer,"no");
#endif

   printf(" OK: %s flat field apodization\n",buffer);
   fprintf(fp1," OK: %s flat field apodization\n",buffer);

#ifdef APODIZATION
/* Creates Hamming apodization file: */
  isize = *nx_ff * *ny_ff * sizeof(float);
  apodi = (float *)isize);
  for(i=0; i < *nx_ff * *ny_ff; i++) apodi[i] = 0.;
#ifdef HAMMING
  jlp_hamming(apodi,*nx_ff,*ny_ff,*nx_ff,ixcent,iycent,nx,ny);
#else
  jlp_blackman(apodi,*nx_ff,*ny_ff,*nx_ff);
#endif

#ifdef DEBUG
#ifdef HAMMING
   strcpy(outfile,"hamming");
   strcpy(outcomments,"Hamming filter: 0.54 + 0.46 cos(pi t / T)");
#else
   strcpy(outfile,"blackman");
   strcpy(outcomments,"Blackman filter: 0.42 + 0.5 cos(pi t / T) + 0.08 cos(2pit/T)");
#endif
   JLP_WRITEIMAG(apodi,nx_ff,ny_ff,nx_ff,outfile,outcomments);
#endif

/* End of case (#ifdef APODIZATION) */
#endif
/****************************************************************/

/* Read flat field file: */
  JLP_VM_READIMAG1(&pntr_ima,nx_ff,ny_ff,ff_name,ff_comments);
  in_ffield = (float *)pntr_ima;

/* Allocate memory space for flat field: */
 isize = *nx_ff * *ny_ff * sizeof(float);
 *ffield = (float *) malloc(isize);

/* Calling "create_ffield" */
  create_ffield(*ffield,in_ffield,apodi,*nx_ff,*ny_ff,ixstart,iystart,
                ixend,iyend,mini_value);

  printf(" Ffield successfully built \n");

/* For debug only:
*/
#ifdef DEBUG
  sprintf(outfile,"ff%s",file_ext);
  strcpy(outcomments,"inverse ffield reduced");
  JLP_WRITEIMAG(*ffield,nx_ff,ny_ff,nx_ff,outfile,outcomments);
#endif

/* Free virtual memory space: */
  free(in_ffield);
#ifdef APODIZATION
  free(apodi);
#endif

}
/********************* Unity file for ffield: *************/
else
{
  fprintf(fp1,"\n OK, Take unity file for ffield (no apodization)\n");
  printf("\n OK, Take unity file for ffield (no apodization)\n");
 isize = *nx_ff * *ny_ff *sizeof(float);
 *ffield = (float *) malloc(isize);
 for(i = 0; i < *nx_ff * *ny_ff; i++) (*ffield)[i] = 1.;
}

return(0);
}
/*******************************************************************
* prepare_output_bisp
* double precision version (from decode_cp40)
********************************************************************/
int prepare_output_bisp(double *bispp, double *modsq, double *snrm, 
                        INT4 nx, INT4 ny, float xframes, float xphot, 
                        INT4 nbeta, INT4 ngamma, INT4 photon_correction, 
                        FILE *fp1)
{
register int i;
INT4 bisp_dim;
double w1, w2, w3;

/* Compute the mean: */
  for(i = 0; i < nx * ny; i++) {
     modsq[i] /= xframes;
     snrm[i] /= xframes;
        }
  for(i=0; i< 4 * ngamma; i++) bispp[i] /= xframes;
/****************************************************************/
/* SNR of modsq. First step: sigma**2 computation 
* (Keep sigma**2, for consistency with normalisation in photon_corr)
*  AFTER MEAN COMPUTATION AND BEFORE PHOTON NOISE CORRECTION!!! */
   for(i = 0; i < nx * ny; i++) {
        snrm[i] = snrm[i] - modsq[i]*modsq[i];
        }

/****************************************************************/
   if(photon_correction) 
       {
       printf(" photon correction with xphot = %f \n",xphot);
       fprintf(fp1," photon correction with xphot = %f \n",xphot);
       photon_corr(bispp,modsq,snrm,&nx,&ny,&xphot,&nbeta,&ngamma);
       }

   bisp_dim = ngamma;
   rearrange_mask(bispp,&ngamma,&bisp_dim);

/****************************************************************/
/* SNR of modsq. Second step: modsq[i]/sqrt(variance) 
*  AFTER PHOTON NOISE CORRECTION!!! */
   for(i = 0; i < nx * ny; i++) 
      {
      if(snrm[i] <= 1.e-4) snrm[i] = 1.e-4;
      snrm[i] = 1. / sqrt(snrm[i]);
      }
/* If no photon correction: (1./sigma), otherwise (modsq/sigma) */
   if(photon_correction)
   {
     for(i = 0; i < nx * ny; i++) snrm[i] *= modsq[i];
   }

/***************** Bispectrum: **********************************/
   for(i=0; i<ngamma; i++) {
/* First computing the variance: */
        w1 =  bispp[2*bisp_dim + i] - bispp[i]*bispp[i];
        w2 =  bispp[3*bisp_dim + i] - bispp[bisp_dim + i]*bispp[bisp_dim + i];
/* Then the sigma (real and imag together): */
        w1 = w1 + w2;
        if(w1 < 1.e-10) w1 = 1.e-10;
        w1 = sqrt(w1);
/* Phase factor of the bispectrum */
        w3 = bispp[i]*bispp[i] + bispp[bisp_dim + i]*bispp[bisp_dim + i];
        w3 = sqrt(w3);
/* Normalization: */
/*
        bispp[i]=bispp[i]/w3;
        bispp[bisp_dim + i] = bispp[bisp_dim + i]/w3;
*/
/* SNR of bispectrum in 3rd line (modulus divided by sigma): */
        bispp[2*bisp_dim + i] = w3/w1;
        }

return(0);
}
/*******************************************************************
* prepare_output_bisp_1D
* double precision version 
********************************************************************/
int prepare_output_bisp_1D(double *bispp, double *modsq, double *snrm, 
                           INT4 nx, INT4 ny, float xframes, float xphot, 
                           INT4 nbeta, INT4 ngamma, INT4 photon_correction,
                           FILE *fp1)
{
register int i, ix, ix_b;
INT4 bisp_dim;
double w1, w2, w3;

if(photon_correction)
{printf("prepare_output_bisp_1D/ Fatal warning: photon correction not allowed here!");
 exit(-1);
}

/* Compute the mean: */
  for(i = 0; i < nx * ny; i++) {
     modsq[i] /= xframes;
     snrm[i] /= xframes;
        }
  for(i = 0; i < 4 * ngamma * nx; i++) bispp[i] /= xframes;
/****************************************************************/
/* SNR of modsq. */ 
   for(i = 0; i < nx * ny; i++) {
        snrm[i] = snrm[i] - modsq[i]*modsq[i];
        }

/****************************************************************/
   bisp_dim = ngamma;

/****************************************************************/
/* SNR of modsq. Second step: modsq[i]/sqrt(variance)  */
   for(i = 0; i < nx * ny; i++) 
      {
      if(snrm[i] <= 1.e-4) snrm[i] = 1.e-4;
      snrm[i] = 1. / sqrt(snrm[i]);
      }

/***************** Bispectrum: **********************************/
   for(ix = 0; ix < nx; ix++) 
   {
   ix_b = ix * bisp_dim * 4;
   for(i=0; i<ngamma; i++) {
/* First computing the variance: */
        w1 =  bispp[2*bisp_dim + i + ix_b] - bispp[i + ix_b]*bispp[i + ix_b];
        w2 =  bispp[3*bisp_dim + i + ix_b] - bispp[bisp_dim + i + ix_b]*bispp[bisp_dim + i + ix_b];
/* Then the sigma (real and imag together): */
        w1 = w1 + w2;
        if(w1 < 1.e-10) w1 = 1.e-10;
        w1 = sqrt(w1);
/* Phase factor of the bispectrum */
        w3 = bispp[i + ix_b]*bispp[i + ix_b] + bispp[bisp_dim + i + ix_b]*bispp[bisp_dim + i + ix_b];
        w3 = sqrt(w3);
/* Normalization: */
/*
        bispp[i + ix_b]=bispp[i + ix_b]/w3;
        bispp[bisp_dim + i + ix_b] = bispp[bisp_dim + i + ix_b]/w3;
*/
/* SNR of bispectrum in 3rd line (modulus divided by sigma): */
        bispp[2*bisp_dim + i + ix_b] = w3/w1;
        }
    }

return(0);
}
/*************************************************************
* Tests to see if the current machine is working properly
* (i.e., compatible with fast autocorrelation used by autoc_car.c)
**************************************************************/
int test_machine()
{
typedef union {
long lg;
short sh[2];
} coord_photon;

coord_photon a, b, c;

/* Load values in a: */
  a.sh[0] = 121;
  a.sh[1] = 112;
  printf(" a = (121) (112) ");
  bin_long_disp(a.lg);
  printf(" a.sh1 = %d  a.sh2 = %d \n",a.sh[1],a.sh[2]);

/* Load values in b: */
  b.sh[0] = 221;
  b.sh[1] = 202;
  printf(" b = (221) (202) ");
  bin_long_disp(b.lg);
  printf(" b.sh1 = %d  b.sh2 = %d \n",b.sh[1],b.sh[2]);

/* Try a + b : */
  c.lg = a.lg + b.lg;
  printf(" Now trying full operation: \n ");
  printf(" c = a + b = (342) (314) ? ");
  bin_long_disp(c.lg);
  printf(" c.sh0 = %d  c.sh1 = %d \n",c.sh[0],c.sh[1]);

  c.sh[0] = a.sh[0] + b.sh[0];
  c.sh[1] = a.sh[1] + b.sh[1];
  printf(" Now separate operation: \n ");
  printf(" (c0 = a0+b0 | c1 = a1+b1) = (342) (314) ? ");
  bin_long_disp(c.lg);
  printf(" c.sh0 = %d  c.sh1 = %d \n",c.sh[0],c.sh[1]);
 
/* Try a - b :
    and adding 10000000 00000000 10000000 00000000
  to avoid negative values: */
  c.lg = a.lg - b.lg + 0x8000 + 0x80000000;
  printf(" Now trying full operation: \n ");
  printf(" c = a - b = (32768-100) (32768-90) ? ");
  bin_long_disp(c.lg);
  printf(" c.sh0 = %d  c.sh1 = %d \n",c.sh[0],c.sh[1]);

  c.sh[0] = a.sh[0] - b.sh[0] + 0x8000;
  c.sh[1] = a.sh[1] - b.sh[1] + 0x8000;
  printf(" Now separate operation: \n ");
  printf(" (c0 = a0-b0 | c1 = a1-b1) = (32768-100) (32768-90) ? ");
  bin_long_disp(c.lg);
  printf(" c.sh0 = %d  c.sh1 = %d \n",c.sh[0],c.sh[1]);
return(0);
}
/*********************************************
* Displays a long integer with binary code
**********************************************/
int bin_long_disp(long a)
{
int i, j;

/* Long integer is four byte long: */
 for (j = 0 ; j < 4 ; j++)
  {
    printf ("  ");    
    for (i = 0 ; i < 8 ; i++)
     {
/* Mask with 8 * 16**7 = 10000000 00000000 00000000 00000000 */
      printf ("%1ld", (a & 0x80000000) >> 31);
      a = a << 1;
     }
  }
 printf ("\n");    
return(0);
}
/*********************************************************************
* Subroutine JULIAN to compute the Julian day of an observation:
*
* The Julian day begins at Greenwich mean noon (at 12 U.T.)
*
* Here also the Gregorian calendar reform is taken into account.
* Thus the day following 1582 October 4 is 1582 October 15.
*
* The B.C. years are counted astronomically. Thus the year
* before the year +1 is the year 0.
*
* INPUT:
* AA, MM, IDD, TIME : year,month, day, time of the observation
* DJUL : Julian day
**********************************************************************/
int julian(double aa, INT4 mm, INT4 idd, double time, double *djul)
{
double day1, year1, date_obs, date_reform;
INT4 month1, ia1, ib1;

  day1 = time/24. + (double)idd;
/* First the year after the 1st March ... */
  if(mm > 2)
    {
     year1 = aa;
     month1 = mm;
    }
   else
    {
     year1 = aa - 1;
     month1 = mm + 12;
    }
 
/* Then check if after the Gregorian reform: */
    date_obs = aa + ((int)(275 * mm / 9) 
               - 2. * (int) ((mm + 9) / 12) + idd - 30 ) / 365.; 
    date_reform = 1582. + 289./365.;
    if(date_obs >= date_reform)
      {
       ia1 = (int) (year1 / 100.);
       ib1 = 2 - ia1 + (int) (((float)ia1)/4.);
      }
    else
       ib1 = 0;
 
/* Now final formula: */
      *djul = (int)(365.25 * year1) + (int)(30.6001 * (month1 + 1))
              + day1 + 1720994.5 + ib1;
 
return(0);
}
/******************************************************************
* Inverse Julian date to standard calendar (and also provide
* date in years, for binary studies)
*
* INPUT:
* time, djul: time and Julian date
*
* OUTPUT:
* aa: year
* mm: month
* idd: day number
******************************************************************/
int inv_julian(double *date_in_years, double *aa, int *mm, int *idd,
               double time, double djul)
{
double aa_0, time_0, djul_0;
int mm_0, idd_0, ifound, max_days[13];
register int i;

/* Loop from 1990 to 2000, to compute the Julian day of Jan 1st at 0H: */
time_0 = 0.;
mm_0 = 1;
idd_0 = 1;
ifound = 0;
for(i = 1990; i <= 2001; i++)
   {
   aa_0 = (double)i;
   julian(aa_0,mm_0,idd_0,time_0,&djul_0);
   if(djul_0 >= djul) {ifound = 1; break;}
   }
if(!ifound) return(-1);
aa_0--;
*aa = aa_0;

/* Date in years: */
   julian(aa_0,mm_0,idd_0,time_0,&djul_0);
  *date_in_years = *aa + (djul + time / 24. - djul_0)/365.25;

/* Look for month now: */
for(i = 1; i <= 12; i++)
   {
   mm_0 = i;
   julian(aa_0,mm_0,idd_0,time_0,&djul_0);
   if(djul_0 >= djul) break;
   }
mm_0--;
*mm = mm_0;

/* Look for day now: */
/* Load number of days: */
for(i = 1; i < 5; i++) {max_days[2 * i - 1] = 31; max_days[2 * i] = 30;}
for(i = 4; i < 6; i++) {max_days[2 * i] = 31; max_days[2 * i + 1] = 30;}
max_days[12] = 31;
/* Case of April: */
if(*aa == 1992 || *aa == 1996) max_days[4] = 29; else max_days[4] = 28;

/* Debug: */
/*
for(i = 1; i < 13; i++) printf(" max_days %d month %d \n",max_days[i],i);
*/

/* Compute the date of the first day of the month:
   idd_0 = 1;
   julian(aa_0,mm_0,idd_0,time_0,&djul_0);
   *idd = (int)(djul - djul_0) + 1;
* Wrong sometimes...
*/

/* Other way: */
for(i = 1; i <= max_days[*mm]; i++)
   {
   idd_0 = i;
   julian(aa_0,mm_0,idd_0,time_0,&djul_0);
   if(djul_0 >= djul) break;
   }
*idd = idd_0 - 1;

/* Debug: 
  printf(" Full date is: %02d-%02d-%04d \n",*idd,*mm,(int)(*aa));
*/

return(0);
}
/***********************************************************
* output_bisp
* Prepare output of bispectrum
*
* Go from double to single precision
* and contraction from 4 to 3 parameters
***********************************************************/
int output_bisp(double *bispp, float *bisp1, int ngamma)
{
register int ng, k;

   for(ng = 0; ng < ngamma; ng++)
      {
      for(k = 0; k < 3; k++)
        {
        bisp1[ng + k*ngamma] = (float)bispp[ng + k*ngamma];
        }
      }
return(0);
}
/***********************************************************
* output_bisp_1D
* Prepare output of bispectrum
*
* Go from double to single precision
* and contraction from 4 to 3 parameters
***********************************************************/
int output_bisp_1D(double *bispp, float *bisp1, int ngamma, int nx)
{
register int ix, ng, k;
int ix_b, ix_bb;

   for(ix = 0; ix < nx; ix++)
   {
   ix_b = ngamma * 3 * ix;
   ix_bb = ngamma * 4 * ix;
   for(ng = 0; ng < ngamma; ng++)
      {
      for(k = 0; k < 3; k++)
        {
        bisp1[ng + k*ngamma + ix_b] = (float)bispp[ng + k*ngamma + ix_bb];
        }
      }
   }
return(0);
}
