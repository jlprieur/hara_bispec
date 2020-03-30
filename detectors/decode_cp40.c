/**************************************************************** 
*  decode_cp40.c
* To decode data from the CP40 photon counting camera 
* and compute bispectrum.
*
* No photon noise correction. 
* Flat field correction (before compression)
*
* JLP
* Version 01-03-01
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
/* CP4T format: */
#define CP4T
#ifdef CP4T
/* For February/March 1996: */
#define CP40KWD "FORMCP4T" 
#else
/* For January 1994: (Old format) */
#define CP40KWD "FORMCQ40" 
#endif

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
#include <jlp_fftw.h>

/* PC/LINUX 2001: integer inversion like DEC */
#define swap_like_dec 1

void swap_lint(long *i), inv_ulong_int(unsigned long *i);
int rcp40_header(char *in_name, double *aa, int *mm, int *idd, double *time,
                 float *integ, int *nphot_header, FILE *fd, FILE *fp1);
int rcp40_phot(double *image1, double *re1, double *im1, 
               double *re2, double *im2, double *modsq, double *snrm, 
               double *bispp, double *long_int, float *ffield, long *ixy_1,
               long max_nphot, long nx, long ny, long idim, int ixstart, 
               int iystart, long ir, int ireduc, long nbeta, long ngamma,
               char *in_name, int nphot_wanted, int *nphot_found, int *nframes,
               int nx_ff, int ny_ff, long distor_corr, FILE *fp1);
/* NB: another version of process_frame is in "decode_cp400.c" */
int process_frame(double *image, double *re1, double *im1, 
                  double *re2, double *im2, double *modsq, double *snrm, 
                  double *bispp, double *long_int, 
                  long nx, long ny, long idim, long ir, long nbeta, long ngamma,
                  long  *ixy_1, long max_nphot, long nphot, long iframe);

/* LK format: */
#include "lk_fmt2.h"
camera *quelle_CP40;

#undef DOUBLE_PRECISION

void main(argc, argv)
int argc;
char **argv;
{
/* Images : */
double *image1, *re1, *re2, *im1, *im2;
double *modsq, *snrm, *long_int, *bispp;
float *ffield, *autoco, *interco;
float xframes, w1, w2, integ, xphot;
double aa, time;
int mm, idd, ixstart, iystart, i1, i2;
int photon_correction, long_integr_only, autocc_computation, distor_corr;
int modsq_computation; 
int nphot_wanted, nphot_header, nphot_found, nframes;
int isize, status, status1, status2, ireduc, ixcent, iycent;
/* Buffer ixy_1 allows to store 100 frames: */
long nx_ff, ny_ff, *ixy_1;
INT4 nx, ny, descr_length, nbeta, ngamma;
INT4 ir, max_nclosure, nx_auto, ny_auto, max_nphot;
register long int i;
char *pc, file_ext[41], buffer[81], logfile[41];
char fmask_name[61], descr_name[61], descr_value[81];
char in_name[61], ff_name[61], outfile[61], outcomments[81];
FILE *fd, *fp1;

printf(" Program decode_cp40 /JLP Version 01-03-01 \n");
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
  printf(" ireduc,ixcent,iycent,nx,ny,ir,max_nclosure,total_nphot,distortion_correction_flag ");
  printf(" flat_field Fourier_mask\n");
  printf(" Example: (96mar22: channel 2: 488,512) \n"); 
  printf(" (96feb27: channel 4: 190,143) \n"); 
  printf("  runs decode_cp40 test tt 16,488,512,20,20,12,900,10000,0 ff_oc45 mask1 \n\n"); 
  printf(" or without flat field and without mask: \n");
  printf("  runs decode_cp40 test tt 16,488,512,20,20,12,900,10000,0 0 0 \n\n"); 
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
sprintf(logfile,"dec_cp40%s.log",file_ext);
if((fp1 = fopen(logfile,"w")) == NULL)
   {
   printf("decode_cp40/Fatal error opening logfile >%s< \n",logfile);
   exit(-1);
   }

fprintf(fp1," Program decode_cp40 \n");
fprintf(fp1," JLP Version 01-03-01 \n");

/***************************************************************/
/* Open the input file */
if((fd = fopen(in_name,"r")) == NULL)
  {
  printf("Fatal error opening input file: >%s< \n",in_name);
  exit(-1);
  }

/* Read header first (for interactive input of the parameters...) */
  status = rcp40_header(in_name,&aa,&mm,&idd,&time,&integ,&nphot_header,fd,fp1);
/* Closes the input file */
  fclose(fd);

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
  printf(" If ir, radius of uv-coverage == 0, long integration only \n\n");
  printf(" If ir, radius of uv-coverage < 0, autoccorelation computation \n\n");
  printf(" Select reduction factor of the output frames, (4, 8, 16 or 32),");
  printf(" IXcenter, IYcenter, NX, NY, radius of uv coverage, max_nclosure, \n");
  printf(" total number of photons to be processed");
  printf(" and photon noise correction (1 for yes, 0 for no):\n");
  printf("     (Example: 16,488,512,12,2000,250,1) \n");
  gets(buffer);sscanf(buffer,"%d,%d,%d,%d,%d,%d,%d,%d,%d",
               &ireduc,&ixcent,&iycent,&i1,&i2,&ir,&max_nclosure,
               &nphot_wanted,&distor_corr);
  nx = i1; ny = i2;
  photon_correction = 0;
  }
else
  {
  sscanf(argv[3],"%d,%d,%d,%d,%d,%d,%d,%d,%d",
               &ireduc,&ixcent,&iycent,&i1,&i2,&ir,&max_nclosure,
               &nphot_wanted,&distor_corr);
  nx = i1; ny = i2;
  photon_correction = 0;
  }
  printf(" OK: ireduc=%d ixcent=%d iycent=%d nx=%d ny=%d \n",
          ireduc, ixcent, iycent, nx, ny);
  printf(" OK: nphot=%d distor_corr=%d\n", nphot_wanted, distor_corr);

  long_integr_only = (ir == 0) ? 1 : 0;
  autocc_computation = (ir < 0) ? 1 : 0;
  modsq_computation = (ir > 0) ? 1 : 0;

/* Logfile: */
  fprintf(fp1," Reduction factor of the output frames: %d \n",ireduc);
  fprintf(fp1," Radius of uv coverage=%d ",(int)ir);
  fprintf(fp1," Max number of closure relations: %d \n",(int)max_nclosure);
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
  {printf(" Fatal error for the output size: nx = %d, ny = %d\n",
          (int)nx,(int)ny);
  fprintf(fp1," Fatal error for the output size: nx = %d, ny = %d\n",
          (int)nx,(int)ny);
   exit(-1);}
/* Take even numbers (otherwise pb in "recentre"..) */
  nx = (nx / 2) * 2; ny = (ny / 2) * 2;
  printf(" Output image size: nx = %d ny = %d \n",(int)nx,(int)ny);
  fprintf(fp1," Output image size: nx = %d ny = %d \n",(int)nx,(int)ny);

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
   if(fmask_name[0] == '0') *fmask_name = '\0';
 }
else
 {
  strcpy(ff_name,argv[4]);
  strcpy(fmask_name,argv[5]);
  if(fmask_name[0] == '0') *fmask_name = '\0';
 }

  if(*ff_name == '0')
     fprintf(fp1," No flat field correction\n");
  else
     fprintf(fp1," OK flat field correction with: %s\n",ff_name);

  if(*fmask_name == '\0')
     fprintf(fp1," No Fourier mask\n");
  else
     fprintf(fp1," OK Fourier mask is: %s\n",fmask_name);

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
if(ixstart < 0) {ixstart = 0; ixcent = (ireduc * nx) / 2;}
if(iystart < 0) {iystart = 0; iycent = (ireduc * ny) / 2;}
}
/* Logfile: */
  fprintf(fp1," IXcenter=%d, IYcenter=%d, NX=%d NY=%d\n",
           ixcent,iycent,(int)nx,(int)ny);

/* Set camera type: */
quelle_CP40 = malloc(30*sizeof(int));
#ifdef CP4T
set_CP4T(quelle_CP40);
#else
set_CQ40(quelle_CP40);
#endif

/* Prepare future FFT's */
fftw_setup(nx,ny);

/*****************************************************************/
/* Computing the spectral and bispectral lists
   corresponding to the selected uv coverage: */

  printf(" Radius of uv-coverage (IR) in pixels: %d\n",(int)ir);
  printf(" Max number of closure relations: %d \n",(int)max_nclosure);

/* Case of bispectrum computation: */
if(ir > 0)
  {
  compute_uv_coverage(fp1,fmask_name,ir,&nbeta,&ngamma,max_nclosure);

  isize = nx * ny * sizeof(double);
    im1 = (double *)malloc(isize);
    re1 = (double *)malloc(isize);
    im2 = (double *)malloc(isize);
    re2 = (double *)malloc(isize);
    modsq = (double *)malloc(isize);
    snrm = (double *)malloc(isize);
  isize = 4 * ngamma * sizeof(double);
    bispp = (double *)malloc(isize);
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
    modsq = (double *)malloc(isize);
    snrm = (double *)malloc(isize);
/* Erasing the arrays : */
    for(i = 0; i < nx_auto * ny_auto; i++) { modsq[i] = 0.; snrm[i] = 0.;}
    }
  }

isize = nx * ny * sizeof(double);
    image1 = (double *)malloc(isize);
    long_int = (double *)malloc(isize);
isize = 100 * max_nphot * sizeof(long);
    ixy_1 = (long *)malloc(isize);

/* Erasing the arrays : */
    for(i = 0; i < nx * ny; i++) { image1[i] = 0.; long_int[i] = 0.;}
    for(i = 0; i < 100 * max_nphot; i++) ixy_1[i] = -2; 

/******************************************************************/
/*************** Flat field *************************************/
/* Computing flat-field */ 
  nx_ff = 400; ny_ff = 300; 
  compute_ffield(&ffield,ff_name,&nx_ff,&ny_ff,
                  ixcent,iycent,(int)(nx*ireduc),(int)(ny*ireduc),
                  file_ext,fp1,0.3);

/* Photon processing: */
  status = rcp40_phot(image1,re1,im1,re2,im2,modsq,snrm,bispp,long_int,
               ffield,ixy_1,max_nphot,nx,ny,nx,ixstart,iystart,ir,ireduc,nbeta,
               ngamma,in_name,nphot_wanted,&nphot_found,&nframes,nx_ff,ny_ff,
               distor_corr, fp1);
  printf(" Output from rcp40: %d photons and %d frames processed\n",
          nphot_found, nframes);
  fprintf(fp1," Output from rcp40: %d photons and %d frames processed\n",
          nphot_found, nframes);

/* Check that some frames have been processed */
   xframes=(float)nframes;

/* Avoid division by zero ... */
/* Long integration only or autoccorelation: */
   if(!modsq_computation)
       {
       xphot = nphot_found;
       fprintf(fp1," xframes=%f \n",xframes);
       printf(" xframes=%f \n",xframes);
       if(xframes <= 0) {xframes = 1.; xphot = 1.;}
       }
/* modsq computation: */
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
   if(modsq_computation)
     {
/* Normalizes FFT (for compatibility with FFT_2D instead of FFT_2D_FAST...*/
     jlp_normalize_fft(bispp,modsq,snrm,&nx,&ny,&nbeta,&ngamma);

     prepare_output_bisp(bispp,modsq,snrm,nx,ny,xframes,xphot,
                         nbeta,ngamma,photon_correction,fp1);
     }
/****************************************************************/
/* Now output the results : */

isize = nx * ny * sizeof(float);

/* Long integration : */
  sprintf(outfile,"long%s",file_ext);
  fprintf(fp1," %s",outfile);
/* To write descriptors, prepare the descriptor string BEFORE
 writing the image file */
/* Leave blanks at the end, for Fortran interface: */
  sprintf(descr_name,"REDUC  ");
  sprintf(descr_value,"%d  ",ireduc);
  descr_length = 2;
/*
  JLP_WDESCR(descr_name,descr_value,&descr_length,&status);
*/
  JLP_D_WRITEIMAG(long_int,&nx,&ny,&nx,outfile,outcomments);

/* modsq, snrm, bispp: */
  if(modsq_computation)
  {
/* Leave blanks at the end, for Fortran interface: */
   sprintf(descr_name,"IR  ");
   sprintf(descr_value,"%d  ",(int)ir);
   descr_length = 4;
/*
   JLP_WDESCR(descr_name,descr_value,&descr_length,&status);
*/
/* Leave blanks at the end, for Fortran interface: */
   sprintf(descr_name,"MAX_NCLO  ");
   sprintf(descr_value,"%d  ",(int)max_nclosure);
   descr_length = 8;
/*
   JLP_WDESCR(descr_name,descr_value,&descr_length,&status);
*/

/* Recentre the frames: */
   RECENT_FFT_DOUBLE(modsq,modsq,&nx,&ny,&nx);
   RECENT_FFT_DOUBLE(snrm,snrm,&nx,&ny,&nx);

/* Mean squared modulus : */
   sprintf(outfile,"modsq%s",file_ext);
   fprintf(fp1," Output of %s",outfile);
   JLP_D_WRITEIMAG(modsq,&nx,&ny,&nx,outfile,outcomments);

/* SNR of squared modulus (actually 1/sigma if no photon correction): */
   sprintf(outfile,"snrm%s",file_ext);
   fprintf(fp1," %s",outfile);
   JLP_D_WRITEIMAG(snrm,&nx,&ny,&nx,outfile,outcomments);

/* Bispectrum : */
   sprintf(outfile,"bisp1%s",file_ext);
   ny=3;
   fprintf(fp1," %s\n",outfile);
   JLP_D_WRITEIMAG(bispp,&ngamma,&ny,&ngamma,outfile,outcomments);
   }
/* Autocorrelation: */
   else if(autocc_computation)
   {
   isize = nx_auto * ny_auto * sizeof(float);
   autoco = (float *)malloc(isize);
   to_single(modsq,autoco,&nx_auto,&ny_auto,&nx_auto);
   free(modsq);
   interco = (float *)malloc(isize);
   to_single(snrm,interco,&nx_auto,&ny_auto,&nx_auto);
   free(snrm);
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
  JLP_END();
  printf(" Output logfile: %s\n",logfile);
fclose(fp1);
}
/********************************************************
*
* Processing the data frame by frame: 
*
* Called by "lk_fmt2.c"
* (Another version is in "decode_cp400.c")
*********************************************************/
int process_frame(double *image, double *re1, double *im1, 
                  double *re2, double *im2, double *modsq, double *snrm, 
                  double *bispp, double *long_int, 
                  long nx, long ny, long idim, long ir, long nbeta, long ngamma,                  long  *ixy_1, long max_nphot, long nphot, long iframe)
{
register int i;

/* JLP96: Do not copy the long integration here, too slow otherwise... */
/* JLP98
  for(i = 0; i < nx * ny; i++) long_int[i] += image[i];
*/

/* JLP96: I add this test to allow for long integration only
*         when ir_max=0: */
  if(nbeta > 0)
   {
/* Resetting the imaginary part for the FFT: */
    for(i = 0; i < nx * ny; i++) {im1[i] = 0.;}

/* Fourier Transform: */
#ifdef NAG_FFT
    FFT_2D_FAST(image,im1,&nx,&ny);
/* Please note that output arrays should be multiplied to sqrt(nx * ny) to
* obtain same results as FFT_2D ...*/
#else
/* GNU Singleton FFT
    fftn_2D_double(image,im1,nx,ny,1);
*/
    fftw_double(image,im1,nx,ny,1);
#endif

/* Processing this image now: bispec3 is without photon noise correction */
    bispec3(image,im1,modsq,snrm,&nx,&ny,bispp,&ir,&nbeta,&ngamma);

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
int rcp40_header(char *in_name, double *aa, int *mm, int *idd, double *time,
                 float *integ, int *nphot_header, FILE *fd, FILE *fp1)
{
float date;
double date_in_years; 
int nbytes_to_read, nbytes, nvalues; 
long not_found;
char *buffer, keyword[9];
unsigned long s_date, integ_time, nphot;
register int i, j, k;
long istart;
phot_buf_rec *photon_buffer;

/* Header of 32 bytes */
typedef struct {
unsigned long integ_time;          /* duree en msec */
unsigned long date;                /* date of observation 
                             (in seconds starting from 01-01-1904 at 0H TU)
                             (i.e. Julian day 2416480.500) */
unsigned long nber_of_photons;     /* number of photons */
unsigned long keyword1;            /* "FORM" in ASCII */
unsigned long keyword2;            /* "CP4T" in ASCII */
unsigned long nber_of_images;      /* number of images */
short refNum;             /* always 0 */
unsigned long read_already;        /* not used */ 
short everything_read;    /* not used */
} LK_HEADER;

union{
long lg[PHOT_BUF_SIZE];
char ch[PHOT_BUF_SIZE*4];
} buff;

LK_HEADER *chead;

/* Load FORMCQ40 or FORMCP4T : */
strcpy(keyword,CP40KWD);

photon_buffer = malloc(30*sizeof(float));

#ifdef DEBUG
printf(" \n rcp40_header/reading header from file : %s \n",in_name);
#endif

/***************************************************************/
/* Read header (file starts with "FORMCP4T", or "FORMCQ40" 
* generally in the second block of 512 bytes): */
  nbytes_to_read = PHOT_BUF_SIZE*4;
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
     printf("rcp40_header/error reading header: \n");
     printf("       Only %d bytes read (nbytes_to_read=%d) \n", nbytes,
            nbytes_to_read);
     return(-2);
    }
*/

/* Decode header: first look for 'F', 
   then compares the following letters to the keyword FORMCP4T or FORMCQ40 */
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
         printf(" i=%d, istart=%d block #%d (of %d bytes) \n",
                  i,(int)istart,k,(int)PHOT_BUF_SIZE*4);
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
#ifdef DEBUG
   printf(" OK, header found \n");
#endif

/*************** date ****************************************/
   s_date = chead->date;
/* When working with DEC computer, reverses long integers... */
#ifdef swap_like_dec
   inv_ulong_int(&s_date);
#endif
/* s_date was in seconds, I convert it to days: */
   date = (float)s_date / 86400.;
/* Time in hours: */
   *time = 24. * (date - (int)date);
/* JLP96: Look for actual date */
   date += DATE_ORIGIN;
   inv_julian(&date_in_years,aa,mm,idd,*time,(double)date);
   printf(" Date is: %f or %02d-%02d-%04d \n",
            date_in_years,*idd,*mm,(int)(*aa));
   fprintf(fp1," Date is: %f or %02d-%02d-%04d \n",
            date_in_years,*idd,*mm,(int)(*aa));

/*************** integration time *******************************/
   integ_time = chead->integ_time;
#ifdef swap_like_dec
   inv_ulong_int(&integ_time);
#endif
/* integration time was in milliseconds, I convert it to seconds: */
   *integ = (float)integ_time / 1000.;
   printf(" integration time =%f (sec) \n",*integ);

  printf(" Input file:%s date=%02d-%02d-%04d time=%.2fH \n",
          in_name,*idd,*mm,(int)(*aa),*time);
  fprintf(fp1," Input file:%s date=%02d-%02d-%04d time=%.2fH \n",
          in_name,*idd,*mm,(int)(*aa),*time);

/*************** number of photons *******************************/
   nphot = chead->nber_of_photons;
#ifdef swap_like_dec
   inv_ulong_int(&nphot);
#endif
   *nphot_header = nphot;
/* JLP96: these data are not correct:
   printf(" Integ. time=%.2f nphotons=%d \n", integ,*nphot_header);
   fprintf(fp1," Integ. time=%.2f nphotons=%d \n", integ,*nphot_header);
*/
   }

/***************************************************************/
/* Shift values of array buff.ch to the origin: */ 
   nvalues = PHOT_BUF_SIZE*4 - (istart + 32);
   for(i = 0; i < nvalues; i++) buff.ch[i] = buff.ch[i + istart + 32];  

  return(0);
}
/****************************************************************/
/*
* Read format of CP40 data 
* and process data 
*/
/****************************************************************/
int rcp40_phot(double *image1, double *re1, double *im1, 
               double *re2, double *im2, double *modsq, double *snrm, 
               double *bispp, double *long_int, float *ffield, long *ixy_1,
               long max_nphot, long nx, long ny, long idim, int ixstart, 
               int iystart, long ir, int ireduc, long nbeta, long ngamma,
               char *in_name, int nphot_wanted, int *nphot_found, int *nframes,
               int nx_ff, int ny_ff, long distor_corr, FILE *fp1)
{
FILE *fd;
int nvalues, nblock, nexpected_blocks, mode_test; 
long iphot, nb_images, istart;
char keyword[9];
register int i;
phot_buf_rec *photon_buffer;
/* JLP 2001: new declarations to call rcp40_header: */
double aa, time; 
int mm, idd, nphot_header, status; 
float integ; 

union{
long lg[PHOT_BUF_SIZE];
char ch[PHOT_BUF_SIZE*4];
} buff;

/* Load FORMCQ40 or FORMCP4T : */
strcpy(keyword,CP40KWD);

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

/* Read header first */
  status = rcp40_header(in_name,&aa,&mm,&idd,&time,&integ,&nphot_header,fd,fp1);
  if(status) {fclose(fd); return(-1);}

/***************************************************************/
/* Shift values of array buff.ch to the origin: */ 
   nvalues = PHOT_BUF_SIZE*4 - (istart + 32);
   for(i = 0; i < nvalues; i++) buff.ch[i] = buff.ch[i + istart + 32];  

/***************************************************************/

/* 4 bytes per photon: */
  nexpected_blocks = nphot_wanted / PHOT_BUF_SIZE;

/* Set nvalues to nphot_wanted if few photons are wanted: */
  if(nexpected_blocks <= 0) nexpected_blocks = 1;

/***************************************************************/
/* Read the data: */
iphot = 0;
nb_images = 0;
mode_test = 0;
/* Assume that less than MAX_IMAGES images within a block: */
photon_buffer->debuts_image = malloc(MAX_IMAGES*sizeof(float));

for(nblock = 0; nblock < NBLOCKMAX; nblock++)
{
/* Set photon_buffer data pointer and number of photons (4-bytes) to be read: */
photon_buffer->adr = buff.lg;
photon_buffer->photons = nphot_wanted;

/* Displays current value of nblock: */
#ifdef DEBUG
  if((nblock % 10) == 1)
     printf(" Block #%d  of %d bytes /%d blocks to be used\n",
              nblock,(int)PHOT_BUF_SIZE*4,nexpected_blocks);
#else
  if((nblock % 20) == 1)
     printf(" Block #%d  of %d bytes /%d blocks to be used\n",
              nblock,(int)PHOT_BUF_SIZE*4,nexpected_blocks);
#endif

CP40_vers_YX9M(photon_buffer, quelle_CP40, image1, nx, ny, 
               ixstart,iystart,ireduc,
               re1,im1,re2,im2,
               modsq,snrm,bispp,
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
  nvalues = fread(buff.lg,4,PHOT_BUF_SIZE,fd);
  if(nvalues != PHOT_BUF_SIZE)
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
         (int)iphot,nphot_wanted);

/* Closes the input file */
  fclose(fd);
  printf("nphot=%d nb_images=%d corresponding exposure time=%.2fs \n",
          (int)iphot,(int)nb_images,0.02*(float)nb_images);
  *nphot_found = iphot;
  *nframes = nb_images;
  return(0);
}
