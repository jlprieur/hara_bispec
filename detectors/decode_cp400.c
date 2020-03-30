/**************************************************************** 
*  decode_cp40.c
* To decode data from the CP40 photon counting camera 
* and compute bispectrum.
*
* No photon noise correction. 
* Flat field correction (before compression)
*
* JLP
* Version 14-09-96
* Same as decode_cp40, but changed to allow for i^2(u)i^*(u) (cf. Letter)
*****************************************************************/
/************************** REMEMBER::: ***********************
* To dump a file on the screen and see Ascii and Hexa codes:
*      od -cx file_name
* each line with 32 bytes
* first block of 512 bytes from first up to line starting with 1000 (Hexa)
* then from 2000 (2nd block), 3000 (3rd block), etc.
* (WARNING: all the lines are not displayed if filled with zeroes...)
*/
/*
#define DEBUG
*/
#define NBLOCKMAX 10000 
/* Useful window for CP40 (January 1995): (40,10) to (360,280) 
 i.e., xwidth=320  ywidth=260 , xcent= 200 ycent= 140*/
/* be carefull coordinates are indeed between [0,512[ ....
*/
/* Good set with center at: 200,140, NX=320, NY=260 */
/* Good "square" set with center at: 168,140, NX=256, NY=256 */

/* Origin of date is 01-01-1904 at 0H TU
 which corresponds to Julian day 2416480.500 */
#define DATE_ORIGIN 2416480.50

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <jlp_ftoc.h>

/* LK format: */
#include "lk_fmt2.h"
camera *quelle_CP40;

#undef DOUBLE_PRECISION

main(argc, argv)
int argc;
char **argv;
{
/* Images : */
double *image1, *re1, *re2, *im1, *im2;
double *modsq, *snrm, *long_int, *bispp;
float *modsq_sngl, *snrm_sngl, *long_sngl, *bispp_sngl;
float *ffield, *autoco, *interco;
float xframes, w1, w2, integ, xphot;
double aa, time;
int mm, idd, ixstart, iystart;
int photon_correction, long_integration_only, distor_corr;
int nphot_wanted, nphot_header, nphot_found, nframes;
int status, status1, status2, ireduc, ireduc_ff, isize, ixcent, iycent;
/* Buffer ixy_1 allows to store 100 frames: */
long int nx, ny, nx_ff, ny_ff, *ixy_1, descr_length;
long int ir, max_nclosure, nbeta, ngamma, nx_auto, ny_auto, max_nphot;
register long int i;
char *pc, file_ext[41], buffer[81], logfile[41];
char fmask_name[61], descr_name[61], descr_value[81];
char in_name[61], ff_name[61], outfile[61], outcomments[81];
FILE *fp1;

printf(" Program decode_cp40 \n");
printf(" JLP Version 24-05-96 \n");
printf(" Enter: Xcenter=0 and Ycenter=0 to start at pixel(0,0)\n");

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
  printf(" runs decode_cp40 in_photon_file output_file_exten");
  printf(" reduc_fact,ixcent,iycent,nx,ny,ir,max_nclosure,total_nphot,distortion_correction_flag ");
  printf(" flat_field Fourier_mask\n");
  printf(" Example: (96mar22: channel 2: 488,512) \n"); 
  printf(" (96feb27: channel 4: 190,143) \n"); 
  printf("  runs decode_cp40 test tt 16,488,512,20,20,12,900 ff_oc45 mask1 \n\n"); 
  printf(" or without flat field and without mask: \n");
  printf("  runs decode_cp40 test tt 16,488,512,20,20,12,1000 0 0 \n\n"); 
  exit(-1);
  }

/* Interactive input of in_name and file_ext: */
if (argc == 1)
 {
   printf(" Input file := "); scanf("%s",in_name);
/************* File extension for output files: */
   printf(" Output file extension := "); scanf("%s",file_ext);
/* For some unknown reason, should read nothing here... */
  gets(buffer);
 }
else
 {
  strcpy(in_name,argv[1]);
  strcpy(file_ext,argv[2]);
 }
  pc = file_ext;
  while(*pc && *pc != ' ') pc++;
  *pc='\0';

#ifdef DEBUG
printf(" DEBUG Version, will read >%s< \n",in_name);
#endif

/**************************************************************/
/* Opening logfile: */
sprintf(logfile,"decode_cp40%s.log",file_ext);
if((fp1 = fopen(logfile,"w")) == NULL)
   {
   printf("decode_cp40/Fatal error opening logfile >%s< \n",logfile);
   exit(-1);
   }

fprintf(fp1," Program decode_cp40 \n");
fprintf(fp1," JLP Version 24-05-96 \n");

/***************************************************************/
/* Read header first (for interactive input of the parameters...) */
  status = rcp40_header(in_name,&aa,&mm,&idd,&time,&integ,&nphot_header,fp1);
  w1 = nphot_header/(50. * integ);
  w2 = nphot_header/integ;
  printf(" Average flux in 20 msec was: %f photons (i.e. %f photons/s)\n",
          w1,w2);
  fprintf(fp1," Average flux in 20 msec was: %f photons (i.e. %f photons/s)\n",
          w1,w2);
/* Set maximum number of photons per frame:
 (Maximum number of photons per 20msec frame): */
  max_nphot = w2 * 1.5;

/***************************************************************/
/* In case of interactive input of parameters, it is better to ask
 the following when the header has been read */ 
if (argc == 1)
  {
  printf(" Reduction factor of 4 implies frame size of 256x256 \n\n");
  printf(" If ir, radius of uv-coverage <= 0, long integration only \n\n");
  printf(" Select reduction factor of the output frames, (4, 8, 16 or 32),");
  printf(" IXcenter, IYcenter, NX, NY, radius of uv coverage, max_nclosure, \n");
  printf(" total number of photons to be processed");
  printf(" and photon noise correction (1 for yes, 0 for no):\n");
  printf("     (Example: 16,488,512,12,2000,250,1) \n");
  gets(buffer);sscanf(buffer,"%d,%d,%d,%d,%d,%d,%d,%d,%d",
               &ireduc,&ixcent,&iycent,&nx,&ny,&ir,&max_nclosure,
               &nphot_wanted,&distor_corr);
  photon_correction = 0;
  }
else
  {
  sscanf(argv[3],"%d,%d,%d,%d,%d,%d,%d,%d,%d",
               &ireduc,&ixcent,&iycent,&nx,&ny,&ir,&max_nclosure,
               &nphot_wanted,&distor_corr);
  photon_correction = 0;
  }

  long_integration_only = (ir <= 3) ? 1 : 0;

/* Logfile: */
  fprintf(fp1," Reduction factor of the output frames: %d \n",ireduc);
  fprintf(fp1," Radius of uv coverage=%d ",ir);
  fprintf(fp1," Max number of closure relations: %d \n",max_nclosure);
  fprintf(fp1," Total number of photons to be processed: %d\n",nphot_wanted);

/* Just to check (for the user): */ 
  if(photon_correction)
     {
     printf(" OK photon noise correction \n");
     fprintf(fp1," OK photon noise correction \n");
     }
  else
     {
     printf(" OK no correction for photon noise\n");
     fprintf(fp1," OK no correction for photon noise\n");
     }

  if(distor_corr) buffer[0] = '\0'; else strcpy(buffer,"no"); 
  printf(" OK: %s correction of geometrical distortion \n",buffer);
  fprintf(fp1," OK: %s correction of geometrical distortion \n",buffer);

/* Check size of output images */
  if(nx <=0 || ny <=0)
  {printf(" Fatal error for the output size: nx = %d, ny = %d\n",nx,ny);
  fprintf(fp1," Fatal error for the output size: nx = %d, ny = %d\n",nx,ny);
   exit(-1);}
/* Take even numbers (otherwise pb in "recentre"..) */
  nx = (nx / 2) * 2; ny = (ny / 2) * 2;
  printf(" Output image size: nx = %d ny = %d \n",nx,ny);
  fprintf(fp1," Output image size: nx = %d ny = %d \n",nx,ny);

/* Photon counting: */
  printf(" Will process %d photons \n",nphot_wanted);
  fprintf(fp1," Will process %d photons \n",nphot_wanted);

  if(nphot_wanted > nphot_header) 
    { printf(" decode_rcar/Warning, only %d photons in file!\n",nphot_header);
     fprintf(fp1," decode_rcar/Warning, only %d photons in file!\n",nphot_header);
     nphot_wanted = nphot_header;
    }

/* Flat field: */
if (argc == 1)
 {
   printf(" Flat field file := ");scanf("%s",ff_name);
   printf(" Fourier mask := ");scanf("%s",fmask_name);
 }
else
 {
  strcpy(ff_name,argv[4]);
  strcpy(fmask_name,argv[5]);
 }

/**********************************************************/
JLP_BEGIN();
JLP_INQUIFMT();

/* Boundaries: */
if(ixcent == 0 && iycent == 0)
{
ixstart = 0;
iystart = 0;
ixcent = (ireduc * nx) / 2;
iycent = (ireduc * ny) / 2;
}
else
{
ixstart = ixcent - (ireduc * nx) / 2; 
iystart = iycent - (ireduc * ny) / 2; 
}
/* Logfile: */
  fprintf(fp1," IXcenter=%d, IYcenter=%d, NX=%d NY=%d\n",
           ixcent,iycent,nx,ny);

/* Set camera type: */
quelle_CP40 = malloc(30*sizeof(int));
set_CP4T(quelle_CP40);

/*****************************************************************/
/* Computing the spectral and bispectral lists
   corresponding to the selected uv coverage: */

  printf(" Radius of uv-coverage (IR) in pixels: %d\n",ir);
  printf(" Max number of closure relations: %d \n",max_nclosure);

/* Case of bispectrum computation: */
if(ir > 0)
  {
  compute_uv_coverage(fp1,fmask_name,ir,&nbeta,&ngamma,max_nclosure);

  isize = nx * ny * sizeof(double);
    JLP_GVM(&re1,&isize);
    JLP_GVM(&im1,&isize);
    JLP_GVM(&re2,&isize);
    JLP_GVM(&im2,&isize);
    JLP_GVM(&modsq,&isize);
    JLP_GVM(&snrm,&isize);
  isize = 4 * ngamma * sizeof(double);
    JLP_GVM(&bispp,&isize);
/* Erasing the arrays : */
    for(i = 0; i < nx * ny; i++) {
      re1[i] = 0.; im1[i] = 0.;
      re2[i] = 0.; im2[i] = 0.;
      modsq[i] = 0.;
      snrm[i] = 0.;
      }
    for(i=0; i < 4 * ngamma; i++) { bispp[i]=0.;}
  }
else
  { nbeta = ir; ngamma = ir;
/* Compute autocorrelation and intercorrelation if ir < 0: */
  if(ir < 0)
    {
/* Autocorrelation has double size in x and y: */
    nx_auto = 2 * nx; ny_auto = 2 * ny; 
    isize = nx_auto * ny_auto * sizeof(double);
/* "modsq" will contain autocorrelation, "snrm" the intercorrelation: */
    JLP_GVM(&modsq,&isize);
    JLP_GVM(&snrm,&isize);
/* Erasing the arrays : */
    for(i = 0; i < nx_auto * ny_auto; i++) { modsq[i] = 0.; snrm[i] = 0.;}
    }
  }

isize = nx * ny * sizeof(double);
  JLP_GVM(&image1,&isize);
  JLP_GVM(&long_int,&isize);
isize = 100 * max_nphot * sizeof(long);
  JLP_GVM(&ixy_1,&isize);

/* Erasing the arrays : */
    for(i = 0; i < nx * ny; i++) { image1[i] = 0.; long_int[i] = 0.;}
    for(i = 0; i < 100 * max_nphot; i++) ixy_1[i] = -2; 

/******************************************************************/
/*************** Flat field *************************************/
/* Calling compute_ffield with compression factor 
*  ("ireduc_ffield") set to unity: */
  nx_ff = 400; ny_ff = 300; ireduc_ff = 1;
  compute_ffield(&ffield,ff_name,&nx_ff,&ny_ff,&ireduc_ff,nx_ff,ny_ff,
                  file_ext,fp1,0.3);

/* Long integration: 
* debug only:
  status = rcp40_long(long_int,ixy_1,max_nphot,nx,ny,ixstart,iystart,ireduc,
                      in_name,outcomments,nphot_wanted,&nphot_found);
*/

/* Photon processing: */
  status = rcp40_phot(image1,re1,im1,re2,im2,modsq,snrm,bispp,long_int,
               ffield,ixy_1,max_nphot,nx,ny,nx,ixstart,iystart,ir,ireduc,nbeta,
               ngamma,in_name,nphot_wanted,&nphot_found,&nframes,nx_ff,ny_ff,
               distor_corr);
  printf(" Output from rcp40: %d photons and %d frames processed\n",
          nphot_found, nframes);
  fprintf(fp1," Output from rcp40: %d photons and %d frames processed\n",
          nphot_found, nframes);

/* Check that some frames have been processed */
   xframes=(float)nframes;

/* Avoid division by zero ... */
/* Long integration only: */
   if(long_integration_only)
       {
       xphot = nphot_found;
       fprintf(fp1," xframes=%f \n",xframes);
       printf(" xframes=%f \n",xframes);
       if(xframes <= 0) {xframes = 1.; xphot = 1.;}
       }
   else
     {
      if(xframes > 0. && modsq[0] > 0.)
        {
         printf(" Normalization: xframes=%f modsq[0]=%f \n",
            xframes,modsq[0]);
         xphot = sqrt((double)(modsq[0]/xframes));
         printf(" xphot (1)= %f\n",xphot);
         xphot *= sqrt((double)(nx * ny));
         printf(" xphot (2)= %f\n",xphot);
        }
      else
        {
        fprintf(fp1," Error (for normalization): xframes=%f modsq[0]=%f \n",
            xframes,modsq[0]);
        printf(" Error (for normalization): xframes=%f modsq[0]=%f \n",
            xframes,modsq[0]);
        xframes = 1.;
        xphot = 1.; 
        }
     }
  sprintf(outcomments,"%s %02d-%02d-%04d %.2fH %.1fs %.1fph %dfr %dph",
          in_name,idd,mm,(int)aa,time,integ,xphot,nframes,nphot_found);

/* Mean of the frames: */
   for(i = 0; i < nx * ny; i++) long_int[i] /= xframes;
/*
   if(!long_integration_only)
     {
     prepare_output_bisp(bispp,modsq,snrm,nx,ny,xframes,xphot,
                         nbeta,ngamma,photon_correction,fp1);
     }
*/
/****************************************************************/
/* Now output the results : */

isize = nx * ny * sizeof(float);

/* Long integration : */
  JLP_GVM(&long_sngl,&isize);
  to_single(long_int,long_sngl,&nx,&ny,&nx);
  JLP_FVM(&long_int);
  sprintf(outfile,"long%s",file_ext);
  fprintf(fp1," %s",outfile);
/* To write descriptors, prepare the descriptor string BEFORE
 writing the image file */
/* Leave blanks at the end, for Fortran interface: */
  sprintf(descr_name,"REDUC  ");
  sprintf(descr_value,"%d  ",ireduc);
  descr_length = 2;
  JLP_WDESCR(descr_name,descr_value,&descr_length,&status);
  JLP_WRITEIMAG(long_sngl,&nx,&ny,&nx,outfile,outcomments);

/* modsq, snrm, bispp: */
  if(!long_integration_only)
  {
/* Leave blanks at the end, for Fortran interface: */
   sprintf(descr_name,"IR  ");
   sprintf(descr_value,"%d  ",ir);
   descr_length = 4;
   JLP_WDESCR(descr_name,descr_value,&descr_length,&status);
/* Leave blanks at the end, for Fortran interface: */
   sprintf(descr_name,"MAX_NCLO  ");
   sprintf(descr_value,"%d  ",max_nclosure);
   descr_length = 8;
   JLP_WDESCR(descr_name,descr_value,&descr_length,&status);
/* */
   JLP_GVM(&modsq_sngl,&isize);
   to_single(modsq,modsq_sngl,&nx,&ny,&nx);
   JLP_FVM(&modsq);
   JLP_GVM(&snrm_sngl,&isize);
   to_single(snrm,snrm_sngl,&nx,&ny,&nx);
   JLP_FVM(&snrm);
   isize = 3 * ngamma * sizeof(float);
   JLP_GVM(&bispp_sngl,&isize);
   output_bisp(bispp, bispp_sngl, ngamma);
   JLP_FVM(&bispp);

/* Recentre the frames: */
   RECENT_FFT(modsq_sngl,modsq_sngl,&nx,&ny,&nx);
   RECENT_FFT(snrm_sngl,snrm_sngl,&nx,&ny,&nx);

/* Mean squared modulus : */
   sprintf(outfile,"modsq%s",file_ext);
   fprintf(fp1," Output of %s",outfile);
   JLP_WRITEIMAG(modsq_sngl,&nx,&ny,&nx,outfile,outcomments);

/* SNR of squared modulus : */
   sprintf(outfile,"snrm%s",file_ext);
   fprintf(fp1," %s",outfile);
   JLP_WRITEIMAG(snrm_sngl,&nx,&ny,&nx,outfile,outcomments);

/* Bispectrum : */
   sprintf(outfile,"bisp1%s",file_ext);
   ny=3;
   fprintf(fp1," %s\n",outfile);
   JLP_WRITEIMAG(bispp_sngl,&ngamma,&ny,&ngamma,outfile,outcomments);
   }
/* Autocorrelation: */
   else if(ir < 0)
   {
   isize = nx_auto * ny_auto * sizeof(float);
   JLP_GVM(&autoco,&isize);
   to_single(modsq,autoco,&nx_auto,&ny_auto,&nx_auto);
   JLP_FVM(&modsq);
   JLP_GVM(&interco,&isize);
   to_single(snrm,interco,&nx_auto,&ny_auto,&nx_auto);
   JLP_FVM(&snrm);
/* Symmetry of autocorrelation and intercorrelation */
   autoco_sym(autoco,nx_auto,ny_auto,nx_auto);
   autoco_sym(interco,nx_auto,ny_auto,nx_auto);

/* Normalization: */
   status1 = normalize_float(autoco,nx_auto,ny_auto,nx_auto);
   status2 = normalize_float(interco,nx_auto,ny_auto,nx_auto);
/****************************************************************/
/* Now output the results : */

/* JLP96:
   if(!status1 && !status2)
*/
   if(status1 != 234)
    {
/* Autocorrelation : */
/* JLP96: don't output this file, since it is not used 
   sprintf(outfile,"autoc%s",file_ext);
   fprintf(fp1," Output of %s",outfile);
   JLP_WRITEIMAG(autoco,&nx_auto,&ny_auto,&nx_auto,outfile,outcomments);
*/

/* Intercorrelation : */
/* JLP96: don't output this file, since it is not used 
   sprintf(outfile,"interc%s",file_ext);
   fprintf(fp1," Output of %s",outfile);
   JLP_WRITEIMAG(interco,&nx_auto,&ny_auto,&nx_auto,outfile,outcomments);
*/

/* Autocorrelation corrected by subtracting the intercorrelation: */
   for(i = 0; i < nx_auto * ny_auto; i++) {
      autoco[i] -= interco[i];
      }
   sprintf(outfile,"autocc%s",file_ext);
   fprintf(fp1," Output of %s",outfile);
   JLP_WRITEIMAG(autoco,&nx_auto,&ny_auto,&nx_auto,outfile,outcomments);
   }
   else
   {
    printf(" Sorry, since pb in 'normalize', do not output autocor and interc\n"
);
   }

   }
   

   fprintf(fp1," Successful end \n");
   printf(" Successful end \n");

/* End : */
end1:
  JLP_END();
  printf(" Output logfile: %s\n",logfile);
fclose(fp1);
}
/***********************************************************
* output_bisp
* Prepare output of bispectrum
*
* Go from double to single precision
* and contraction from 4 to 3 parameters
***********************************************************/
int output_bisp(bispp, bisp1, ngamma)
double *bispp;
float *bisp1;
int ngamma;
{
register int ng, k;

   for(ng = 0; ng < ngamma; ng++)
      {
      for(k = 0; k < 3; k++)
        {
        bisp1[ng + k*ngamma] = bispp[ng + k*ngamma];
        }
      }
return(0);
}
/********************************************************
*
* Processing the data frame by frame: 
*********************************************************/
int process_frame(image, re1, im1, re2, im2,  modsq, snrm, bispp, long_int, 
                  nx, ny, idim, ir, nbeta, ngamma, ixy_1, max_nphot, 
                  nphot, iframe)
double image[], re1[], im1[], re2[], im2[];
double modsq[], snrm[], bispp[], long_int[];
long nx, ny, idim, ir, nbeta, ngamma, *ixy_1, nphot, iframe, max_nphot;
{
float w1;
register int i;

/* JLP96: Do not copy the long integration here, too slow otherwise... */
/* JLP98
  for(i = 0; i < nx * ny; i++) long_int[i] += image[i];
*/

/* JLP96: I add this test to allow for long integration only
*         when ir_max=0: */
  if(nbeta > 0)
   {
/* Computing the image to the square (stored in re2): */
    for(i = 0; i < nx * ny; i++) {re2[i] = image[i]*image[i];}

/* Resetting the imaginary part for the FFT: */
    for(i = 0; i < nx * ny; i++) {im1[i] = 0.; im2[i] = 0.;}

/* Fourier Transform: */
    FFT_2D_FAST(image,im1,&nx,&ny);
/* Fourier Transform: */
    FFT_2D_FAST(re2,im2,&nx,&ny);

/* Processing this image now: bispec3 is without photon noise correction */
/* JLP97
    bispec3(image,im,modsq,snrm,&nx,&ny,bispp,&ir,&nbeta,&ngamma);
*/
/* Store result in modsq (real part) and snrm (imag. part): */
    for(i = 0; i < nx * ny; i++) {
         modsq[i] += re1[i] * re2[i] + im1[i] * im2[i]; 
         snrm[i] += im1[i] * re2[i] - re1[i] * im2[i]; 
         }

/* Erasing work image: */
  for(i = 0; i < nx * ny; i++) {image[i] = 0.;}
   }
/* When nbeta is negative, compute mean autocorrelation (in modsq) 
* and intercorrelation (in snrm): */
  else if(nbeta < 0)
   {
    autocor_cp40(modsq,snrm,ixy_1,max_nphot,nphot,iframe,nx,ny);
   }
   
return(0);
}
/****************************************************************/
/*
* Read format of CP40 data 
* and create a long integration
*/
/****************************************************************/
/* Subroutine rdcp40 to read CP40 format */
int rcp40_header(in_name,aa,mm,idd,time,integ,nphot_header,fp1)
double *aa, *time;
float *integ;
int *mm, *idd, *nphot_header;
char in_name[];
FILE *fp1;
{
FILE *fd;
float date;
double date_in_years; 
int nbytes_to_read, nbytes, nvalues, nblock, mode_test; 
long ix, iy, not_found;
char *buffer, keyword[9];
unsigned long s_date, integ_time, nphot;
register int i, j, k;
long itest, istart;
char ctest[5];
phot_buf_rec *photon_buffer;

/* Header of 32 bytes */
typedef struct {
unsigned long integ_time;          /* duree en msec */
unsigned long date;                /* date of observation 
                             (in seconds starting from 01-01-1904 at 0H TU)
                             (i.e. Julian day 2416480.500) */
unsigned long nber_of_photons;     /* number of photons */
long keyword1;            /* "FORM" in ASCII */
long keyword2;            /* "CP4T" in ASCII */
long nber_of_images;      /* number of images */
short refNum;             /* always 0 */
long read_already;        /* not used */ 
short everything_read;    /* not used */
} LK_HEADER;

union{
long lg[PHOT_BUF_SIZE/sizeof(long)];
char ch[PHOT_BUF_SIZE];
} buff;

LK_HEADER *chead;
void swap_lint(), inv_lint();

strcpy(keyword,"FORMCP4T");

photon_buffer = malloc(30*sizeof(float));

#ifdef DEBUG
printf(" \n rdcar/reading file : %s \n",in_name);
#endif

/* Opens the input file */
if((fd = fopen(in_name,"r")) == NULL)
  {
  printf("rdcar_rdheader/error opening input file: >%s< \n",in_name);
  return(-1);
  }

/***************************************************************/
/* Read header (file starts with "FORMCP4T", generally in the second
   block of 512 bytes): */
  nbytes_to_read = PHOT_BUF_SIZE;
  k=0;
  not_found = 1;
  while(not_found && k < 5)
  {
   nbytes = fread(buff.ch,sizeof(char),nbytes_to_read,fd);
#ifdef DEBUG
   printf("        %d bytes read \n", nbytes);
#endif
/* I remove this since I increased the size of the buffer */
/*
   if(nbytes != nbytes_to_read)
    {
     printf("rdcar/error reading header: \n");
     printf("       Only %d bytes read (nbytes_to_read=%d) \n", nbytes,
            nbytes_to_read);
     return(-2);
    }
*/

/* Decode header: first look for 'F', 
   then compares the following letters to the keyword FORMCP4T */
  for(i = 0; i < nbytes; i++) 
    {
    if(buff.ch[i] == 'F')
      {
      buffer = &buff.ch[i];
      for(j = 0; j < 8; j++) 
         if(buffer[j] != keyword[j]) break; 
#ifdef DEBUG
         printf(" i=%d j=%d buffer=%s \n",i,j,buffer);
#endif
/* OK: Successfull research: */
      if(j == 8) 
        {
         not_found = 0; istart = i - 12; 
         chead = (LK_HEADER *)&(buff.ch[istart]);
#ifdef DEBUG
         printf(" i=%d, istart=%d block #%d (of %d bytes) \n",i,istart,k,PHOT_BUF_SIZE);
#endif
/* Just to exit from the loop: */
         i = nbytes + 1; break;
        }
      } 
    }

/* Goes to next block */
 k++;
 }

/*************************************************************/
/* Check if it is OK: */
 if(not_found) 
  {
  printf(" Sorry, keyword has not been found \n");
  return(-1);
  }
 else
   {
   printf(" OK, header found \n");

/*************** date ****************************************/
   s_date = chead->date;
/* When working with DEC computer, reverses long integers... */
#ifdef dec
   inv_lint(&s_date);
#endif
/* s_date was in seconds, I convert it to days: */
   date = (float)s_date / 86400.;
/* Time in hours: */
   *time = 24. * (date - (int)date);
/* JLP96: Look for actual date */
   date += DATE_ORIGIN;
   inv_julian(&date_in_years,aa,mm,idd,*time,(double)date);
   printf(" Date is: %f or %02d-%02d-%04d \n",date_in_years,*idd,*mm,(int)(*aa));
   fprintf(fp1," Date is: %f or %02d-%02d-%04d \n",date_in_years,*idd,*mm,(int)(*aa));

/*************** integration time *******************************/
   integ_time = chead->integ_time;
#ifdef dec
   inv_lint(&integ_time);
#endif
/* integration time was in milliseconds, I convert it to seconds: */
   *integ = (float)integ_time / 1000.;
   printf(" integration time =%f (sec) \n",*integ);

  printf(" Input file:%s date=%02d-%02d-%04d time=%.2fH \n",
          in_name,*idd,*mm,(int)(*aa),time);
  fprintf(fp1," Input file:%s date=%02d-%02d-%04d time=%.2fH \n",
          in_name,*idd,*mm,(int)(*aa),time);

/*************** number of photons *******************************/
   nphot = chead->nber_of_photons;
#ifdef dec
   inv_lint(&nphot);
#endif
   *nphot_header = nphot;
   printf(" Integ. time=%.2f nphotons=%d \n", integ,*nphot_header);
   fprintf(fp1," Integ. time=%.2f nphotons=%d \n", integ,*nphot_header);
   }

/***************************************************************/
/* Shift values of array buff.ch to the origin: */ 
   nvalues = PHOT_BUF_SIZE - (istart + 32);
   for(i = 0; i < nvalues; i++) buff.ch[i] = buff.ch[i + istart + 32];  
/* Going back to long values: */
   nvalues = nvalues/sizeof(long);

/* Closes the input file */
  fclose(fd);
  return(0);
}
/****************************************************************/
/*
* Read format of CP40 data 
* and create a long integration
* (used for debugging)
*/
/****************************************************************/
/* Subroutine rdcp40 to read CP40 format */
int rcp40_long(real_array,ixy_1,max_nphot,nx,ny,ixstart,iystart,ireduc,
               in_name,comments,nphot_wanted,nphot_found)
double real_array[];
long nx, ny, ixstart, iystart, ireduc, nphot_wanted, *nphot_found; 
long *ixy_1, max_nphot;
char in_name[], comments[];
{
FILE *fd;
int nbytes_to_read, nbytes, nvalues, nblock, nexpected_blocks, mode_test; 
long ix, iy, not_found, iphot, nb_images;
char *buffer, keyword[9];
register int i, j, k;
long itest, istart;
char ctest[5];
phot_buf_rec *photon_buffer;

/* Header of 32 bytes */
typedef struct {
unsigned long integ_time;          /* duree en msec */
unsigned long date;                /* date of observation 
                             (in seconds starting from 01-01-1904 at 0H TU)
                             (i.e. Julian day 2416480.500) */
unsigned long nber_of_photons;     /* number of photons */
long keyword1;            /* "FORM" in ASCII */
long keyword2;            /* "CP4T" in ASCII */
long nber_of_images;      /* number of images */
short refNum;             /* always 0 */
long read_already;        /* not used */ 
short everything_read;    /* not used */
} LK_HEADER;

union{
long lg[PHOT_BUF_SIZE/sizeof(long)];
char ch[PHOT_BUF_SIZE];
} buff;

LK_HEADER *chead;
void swap_lint(), inv_lint();

strcpy(keyword,"FORMCP4T");

photon_buffer = malloc(30*sizeof(float));

#ifdef DEBUG
printf(" \n rcp40_long/reading file : %s \n",in_name);
#endif

/* Opens the input file */
if((fd = fopen(in_name,"r")) == NULL)
  {
  printf("rcp40_long/error opening input file: >%s< \n",in_name);
  return(-1);
  }

/***************************************************************/
/* Read header (file starts with "FORMCP4T", generally in the second
   block of 512 bytes): */
  nbytes_to_read = PHOT_BUF_SIZE;
  k=0;
  not_found = 1;
  while(not_found && k < 5)
  {
   nbytes = fread(buff.ch,sizeof(char),nbytes_to_read,fd);
#ifdef DEBUG
   printf("        %d bytes read \n", nbytes);
#endif
/* I remove this since I increased the size of the buffer */
/*
   if(nbytes != nbytes_to_read)
    {
     printf("rdcar/error reading header: \n");
     printf("       Only %d bytes read (nbytes_to_read=%d) \n", nbytes,
            nbytes_to_read);
     return(-2);
    }
*/

/* Decode header: first look for 'F', 
   then compares the following letters to the keyword FORMCP4T */
  for(i = 0; i < nbytes; i++) 
    {
    if(buff.ch[i] == 'F')
      {
      buffer = &buff.ch[i];
      for(j = 0; j < 8; j++) 
         if(buffer[j] != keyword[j]) break; 
#ifdef DEBUG
         printf(" i=%d j=%d buffer=%s \n",i,j,buffer);
#endif
/* OK: Successfull research: */
      if(j == 8) 
        {
         not_found = 0; istart = i - 12; 
         chead = (LK_HEADER *)&(buff.ch[istart]);
#ifdef DEBUG
         printf(" i=%d, istart=%d block #%d (of %d bytes) \n",i,istart,k,PHOT_BUF_SIZE);
#endif
/* Just to exit from the loop: */
         i = nbytes + 1; break;
        }
      } 
    }

/* Goes to next block */
 k++;
 }

/*************************************************************/
/* Check if it is OK: */
 if(not_found) 
  {
  printf(" Sorry, keyword has not been found \n");
  return(-1);
  }
 else
  {
   printf(" OK, header found \n");
  }

/***************************************************************/
/* Shift values of array buff.ch to the origin: */ 
   nvalues = PHOT_BUF_SIZE - (istart + 32);
   for(i = 0; i < nvalues; i++) buff.ch[i] = buff.ch[i + istart + 32];  
/* Going back to long values: */
   nvalues = nvalues/sizeof(long);

/***************************************************************/
/* Read the data: */
iphot = 0;
nb_images = 0;
mode_test = 0;
/* Assume that less than MAX_IMAGES images within a block: */
photon_buffer->debuts_image = malloc(MAX_IMAGES*sizeof(float));

/* Assuming 4 bytes per photon: */
  nexpected_blocks = nphot_wanted * 4 / PHOT_BUF_SIZE;

/* Set nvalues to nphot_wanted if few photons are wanted: */
  if(nexpected_blocks <= 0) 
     {
     nexpected_blocks = 1;
     nvalues = nphot_wanted * 4;
     }

for(nblock = 0; nblock < NBLOCKMAX; nblock++)
{
/* Set photon_buffer data pointer and number of values: */
photon_buffer->adr = buff.lg;
photon_buffer->photons = nvalues;

/* Displays current value of nblock: */
#ifdef DEBUG
  if((nblock % 10) == 1)
     printf(" Block #%d  of %d bytes /%d blocks to be used\n",
              nblock,PHOT_BUF_SIZE,nexpected_blocks);
#else
  if((nblock % 20) == 1)
     printf(" Block #%d  of %d bytes /%d blocks to be used\n",
              nblock,PHOT_BUF_SIZE,nexpected_blocks);
#endif

CP40_vers_YX9M(photon_buffer, quelle_CP40, real_array, nx, ny, 
               ixstart, iystart, ireduc, ixy_1, max_nphot, &iphot, mode_test, 
               nvalues);
/* Always null...
lk_nphot += photon_buffer->photons;
*/
nb_images += photon_buffer->nb_images;

/* Exit from loop if number of photons has been reached: */
   if(iphot >= nphot_wanted) break;

/* End of current block: read next values */
  nvalues = fread(buff.lg,sizeof(long),PHOT_BUF_SIZE/sizeof(long),fd);
  if(nvalues != PHOT_BUF_SIZE/sizeof(long))
          printf(" rdcp40/Warning, only %d values read in block #%d\n",
                   nvalues,nblock);
  if(nvalues <= 0 )
  {
  printf(" rdcp40/end of file: >%s< \n",in_name);
  nblock = NBLOCKMAX;
  }
/* Going to process these data */
}

/* Exit before reaching nphot: */
if(iphot < nphot_wanted)
printf(" Sorry only %d photons read, (should have been %d photons)\n",
         iphot,nphot_wanted);

/* Closes the input file */
  fclose(fd);
  printf("nphot=%d nb_images=%d corresponding exposure time=%.2fs \n",
          iphot,nb_images,0.02*(float)nb_images);
  *nphot_found = iphot;
  return(0);
}
/****************************************************************/
/*
* Read format of CP40 data 
* and process data 
*/
/****************************************************************/
int rcp40_phot(image1,re1,im1,re2,im2,modsq,snrm,bispp,long_int,ffield,ixy_1,
               max_nphot,nx,ny,idim,ixstart,iystart,ir,ireduc,nbeta,ngamma,
               in_name,nphot_wanted,nphot_found,nframes,nx_ff,ny_ff,
               distor_corr)
double image1[], re1[], im1[], re2[], im2[];
double modsq[], snrm[], bispp[], long_int[];
float ffield[];
long nx, ny, idim, ir, nbeta, ngamma;
int ixstart, iystart, ireduc, nx_ff, ny_ff;
int *nphot_found, nphot_wanted, *nframes;
long distor_corr, *ixy_1, max_nphot;
char in_name[];
{
FILE *fd;
int nbytes_to_read, nbytes, nvalues, nblock, nexpected_blocks, mode_test; 
long ix, iy, not_found, iphot, nb_images;
char *buffer, keyword[9];
register int i, j, k;
long itest, istart;
char ctest[5];
phot_buf_rec *photon_buffer;

/* Header of 32 bytes */
typedef struct {
unsigned long integ_time;          /* duree en msec */
unsigned long date;                /* date of observation 
                             (in seconds starting from 01-01-1904 at 0H TU)
                             (i.e. Julian day 2416480.500) */
unsigned long nber_of_photons;     /* number of photons */
long keyword1;            /* "FORM" in ASCII */
long keyword2;            /* "CP4T" in ASCII */
long nber_of_images;      /* number of images */
short refNum;             /* always 0 */
long read_already;        /* not used */ 
short everything_read;    /* not used */
} LK_HEADER;

union{
long lg[PHOT_BUF_SIZE/sizeof(long)];
char ch[PHOT_BUF_SIZE];
} buff;

LK_HEADER *chead;
void swap_lint(), inv_lint();

strcpy(keyword,"FORMCP4T");

photon_buffer = malloc(30*sizeof(float));

#ifdef DEBUG
printf(" \n rcp40_phot/reading file : %s \n",in_name);
#endif

/* Opens the input file */
if((fd = fopen(in_name,"r")) == NULL)
  {
  printf("Fatal error (rcp40_phot)/error opening input file: >%s< \n",in_name);
  exit(-1);
  }

/***************************************************************/
/* Read header (file starts with "FORMCP4T", generally in the second
   block of 512 bytes): */
  nbytes_to_read = PHOT_BUF_SIZE;
  k=0;
  not_found = 1;
  while(not_found && k < 5)
  {
   nbytes = fread(buff.ch,sizeof(char),nbytes_to_read,fd);
#ifdef DEBUG
   printf("        %d bytes read \n", nbytes);
#endif

/* I remove this since I increased the size of the buffer */
/*
   if(nbytes != nbytes_to_read)
    {
     printf("rdcar/error reading header: \n");
     printf("       Only %d bytes read (nbytes_to_read=%d) \n", nbytes,
            nbytes_to_read);
     return(-2);
    }
*/

/* Decode header: first look for 'F', 
   then compares the following letters to the keyword FORMCP4T */
  for(i = 0; i < nbytes; i++) 
    {
    if(buff.ch[i] == 'F')
      {
      buffer = &buff.ch[i];
      for(j = 0; j < 8; j++) 
         if(buffer[j] != keyword[j]) break; 
#ifdef DEBUG
         printf(" i=%d j=%d buffer=%s \n",i,j,buffer);
#endif
/* OK: Successfull research: */
      if(j == 8) 
        {
         not_found = 0; istart = i - 12; 
         chead = (LK_HEADER *)&(buff.ch[istart]);
#ifdef DEBUG
         printf(" i=%d, istart=%d block #%d (of %d bytes) \n",i,istart,k,PHOT_BUF_SIZE);
#endif
/* Just to exit from the loop: */
         i = nbytes + 1; break;
        }
      } 
    }

/* Goes to next block */
 k++;
 }

/*************************************************************/
/* Check if it is OK: */
 if(not_found) 
  {
  printf(" Sorry, keyword has not been found \n");
  return(-1);
  }
 else
  {
   printf(" OK, header found \n");
  }

/***************************************************************/
/* Shift values of array buff.ch to the origin: */ 
   nvalues = PHOT_BUF_SIZE - (istart + 32);
   for(i = 0; i < nvalues; i++) buff.ch[i] = buff.ch[i + istart + 32];  
/* Going back to long values: */
   nvalues = nvalues/sizeof(long);

/***************************************************************/

/* Assuming 4 bytes per photon: */
  nexpected_blocks = nphot_wanted * 4 / PHOT_BUF_SIZE;

/* Set nvalues to nphot_wanted if few photons are wanted: */
  if(nexpected_blocks <= 0) 
     {
     nexpected_blocks = 1;
     nvalues = nphot_wanted * 4;
     }

/***************************************************************/
/* Read the data: */
iphot = 0;
nb_images = 0;
mode_test = 0;
/* Assume that less than MAX_IMAGES images within a block: */
photon_buffer->debuts_image = malloc(MAX_IMAGES*sizeof(float));

for(nblock = 0; nblock < NBLOCKMAX; nblock++)
{
/* Set photon_buffer data pointer and number of values: */
photon_buffer->adr = buff.lg;
photon_buffer->photons = nvalues;

/* Displays current value of nblock: */
#ifdef DEBUG
  if((nblock % 10) == 1)
     printf(" Block #%d  of %d bytes /%d blocks to be used\n",
              nblock,PHOT_BUF_SIZE,nexpected_blocks);
#else
  if((nblock % 20) == 1)
     printf(" Block #%d  of %d bytes /%d blocks to be used\n",
              nblock,PHOT_BUF_SIZE,nexpected_blocks);
#endif

CP40_vers_YX9M(photon_buffer, quelle_CP40, image1, nx, ny, 
               ixstart,iystart,ireduc,re1,im1,re2,im2,modsq,snrm,bispp,
               long_int,ffield,nx_ff,ny_ff,
               ir,nbeta,ngamma,
               &iphot,ixy_1,max_nphot,distor_corr,mode_test);
/* Always null ...
lk_nphot += photon_buffer->photons;
*/
nb_images += photon_buffer->nb_images;

/* Exit from loop if number of photons has been reached: */
   if(iphot >= nphot_wanted) break;

/* End of current block: read next values */
  nvalues = fread(buff.lg,sizeof(long),PHOT_BUF_SIZE/sizeof(long),fd);
  if(nvalues != PHOT_BUF_SIZE/sizeof(long))
          printf(" rdcp40/Warning, only %d values read in block #%d\n",
                   nvalues,nblock);
  if(nvalues <= 0 )
  {
  printf(" rdcp40/end of file: >%s< \n",in_name);
  nblock = NBLOCKMAX;
  }
/* Going to process these data */
}

/* Exit before reaching nphot: */
if(iphot < nphot_wanted)
printf(" Sorry only %d photons read, (should have been %d photons)\n",
         iphot,nphot_wanted);

/* Closes the input file */
  fclose(fd);
  printf("nphot=%d nb_images=%d corresponding exposure time=%.2fs \n",
          iphot,nb_images,0.02*(float)nb_images);
  *nphot_found = iphot;
  *nframes = nb_images;
  return(0);
}
