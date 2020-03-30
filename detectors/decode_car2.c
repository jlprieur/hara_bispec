/**************************************************************** 
  decode_car2.c
 To decode data from the CAR photon counting camera 
 and compute bispectrum.

 Photon noise correction. 
 Uses twice the data (interleaved by half a frame).
 Flat field correction (before compression)

 JLP
 Version 16-07-01
*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <jlp_ftoc.h>
#include <jlp_fftw.h>

#define DEBUG
/*
*/
#define swap_like_dec

/* In "decode_set.c": */
int prepare_output_bisp(double *bispp, double *modsq, double *snrm,
                        INT4 nx, INT4 ny, float xframes, float xphot,
                        INT4 nbeta, INT4 ngamma, INT4 photon_correction,
                        FILE *fp1);
void inv_ulong_int(unsigned long *i);

/* Contained here: */
int rdcar_long(double *long_int, INT4 nx, INT4 ny, char *in_name,
               INT4 nphot_wanted, INT4 *nphot_found);
int rdcar_phot(double *image1, double *image2, double *mim, 
               double *modsq, double *snrm, double *bispp, double *long_int,
               float *ffield, INT4 nx, INT4 ny, 
               INT4 ixcent, INT4 iycent, INT4 ir, INT4 nbeta, INT4 ngamma, 
               char *in_name, INT4 npack, INT4 nphot_wanted, INT4 *nphot_found,
               INT4 *nframes, INT4 nframes_expected, INT4 idim_ff, INT4 ireduc);
int process_frame(double *image, double *mim, double *modsq, 
                 double *snrm, double *bispp, double *long_int, 
                 INT4 nx, INT4 ny, INT4 ir, INT4 nbeta, INT4 ngamma);

#undef DOUBLE_PRECISION

#define NBLOCKMAX 100000

/* Maximum size for the images */
#define IDIM 256
/* CAR detector has 1024x1024 pixels, but useful size is 864x864 : */
/* #define CAR_WIDTH 864 */ 
#define CAR_WIDTH 1024

int main(int argc, char **argv)
{
/* Images : */
double *image1, *image2, *mim;
double *modsq, *snrm, *long_int, *bispp;
float *float_array;
INT_PNTR pntr_ima;
float *ffield, *fmask;
float xframes, w1, w2, date, integ, xphot;
double aa, time, date_in_years;
long int mm, idd;
INT4 nframes_expected, photon_correction;
INT4 npack, nphot_wanted, nphot_header, nphot_found, nframes;
INT4 status, ireduc, isize, ixcent, iycent;
INT4 nx, ny, nx_ff, ny_ff, nx_fmask, ny_fmask;
INT4 ir, max_nclosure, nbeta, ngamma, descr_length;
register long int i;
char *pc, file_ext[41], buffer[81], logfile[41];
char fmask_name[61], fmask_comments[81], descr_name[61], descr_value[81];
char in_name[61], ff_name[61], outfile[61], outcomments[81];
FILE *fp1;

printf(" Program decode_car2  (%d x %d maxi) \n",IDIM,IDIM);
printf(" (Uses twice the data (interleaved by half a frame)) \n");
printf(" JLP Version 15-07-01 \n");

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
  printf(" runs decode_car2 in_photon_file output_file_exten");
  printf(" reduc_fact,ixcent,iycent,nx,ny,ir,max_nclosure,total_nphot,frame_nphot ");
  printf(" flat_field Fourier_mask\n");
  printf(" Example: \n"); 
  printf("  runs decode_car test tt 16,488,512,256,256,12,900,100 ff_oc45 mask1 \n\n"); 
  printf(" or without flat field and without mask: \n");
  printf("  runs decode_car test tt 16,488,512,256,256,12,1000,100 0 0 \n\n"); 
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
sprintf(logfile,"decode_car2%s.log",file_ext);
if((fp1 = fopen(logfile,"w")) == NULL)
   {
   printf("decode_car2/Fatal error opening logfile >%s< \n",logfile);
   exit(-1);
   }

fprintf(fp1," Program decode_car2  (%d x %d maxi) \n",IDIM,IDIM);
fprintf(fp1," (Uses twice the data (interleaved by half a frame)) \n");
fprintf(fp1," JLP Version 15-07-01 \n");

/***************************************************************/
/* Read header first (for interactive input of the parameters...) */
  status = rdcar_header(in_name,&date_in_years,&time,&idd,&mm,&aa,
                        &integ,&nphot_header);
  printf(" Input file:%s date=%.4f time=%.2fH itime=%.2f nphot=%d \n",
          in_name,date_in_years,time,integ,nphot_header);
  fprintf(fp1," Input file:%s date=%.4f (%02d-%02d-%4d) \
time=%.2fH itime=%.2f nphot=%d \n",
          in_name,date_in_years,(int)idd,(int)mm,(int)aa,time,integ,nphot_header);

  w1 = nphot_header/(50. * integ);
  w2 = nphot_header/integ;
  printf(" Average flux in 20 msec was: %f photons (i.e. %f photons/s)\n",
          w1,w2);
  fprintf(fp1," Average flux in 20 msec was: %f photons (i.e. %f photons/s)\n",
          w1,w2);

/***************************************************************/
/* In case of interactive input of parameters, it is better to ask
 the following when the header has been read */ 
if (argc == 1)
  {
  printf(" Reduction factor of 4 implies frame size of 256x256 \n\n");
  printf(" Select reduction factor of the output frames, (4, 8, 16 or 32),");
  printf(" IXcenter, IYcenter, nx, ny, radius of uv coverage, max_nclosure, \n");
  printf(" total number of photons to be processed,");
  printf(" number of photons per frame ");
  printf(" and photon noise correction (1 for yes, 0 for no):\n");
  printf("     (Example: 16,488,512,32,32,12,2000,250,1) \n");
  gets(buffer);sscanf(buffer,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
               &ireduc,&ixcent,&iycent,&nx,&ny,&ir,&max_nclosure,
               &nphot_wanted,&npack,&photon_correction);
  }
else
  {
  sscanf(argv[3],"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
               &ireduc,&ixcent,&iycent,&nx,&ny,&ir,&max_nclosure,
               &nphot_wanted,&npack,&photon_correction);
  }

/* Logfile: */
  fprintf(fp1," Reduction factor of the output frames: %d \n",ireduc);
  fprintf(fp1," IXcenter=%d, IYcenter=%d, radius of uv coverage=%d\n",
           ixcent,iycent,ir);
  fprintf(fp1," Max number of closure relations: %d \n",max_nclosure);
  fprintf(fp1," Total number of photons to be processed: %d\n",nphot_wanted);
  fprintf(fp1," Number of photons per frame: %d \n",npack);

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
/* Set size of output images */
/* Take even numbers (otherwise pb in "recentre"..) */
  nx = (nx / 2) * 2; ny = (ny / 2) * 2;
  printf(" Output image size: nx = %d ny = %d \n",nx,ny);
  fprintf(fp1," Output image size: nx = %d ny = %d \n",nx,ny);

/* Photon counting: */
  printf(" Will process %d photons in packs of %d photons \n",
          nphot_wanted,npack);
  fprintf(fp1," Will process %d photons in packs of %d photons \n",
          nphot_wanted,npack);

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

/*****************************************************************/
/* Computing the spectral and bispectral lists
   corresponding to the selected uv coverage: */

  printf(" Radius of uv-coverage (IR) in pixels: %d\n",ir);
  printf(" Max number of closure relations: %d \n",max_nclosure);

  if(fmask_name[0] == '0')
    {
/* Old call: (either in jlp_cover1.for or in jlp_cover2.c)
    COVERA(&ir,&nbeta,&ngamma);
*/
/* Full u-v coverage: */
    printf(" Full u-v coverage.\n");
    fprintf(fp1," Full u-v coverage.\n");
    nx_fmask = 2 * ir + 1; ny_fmask = nx_fmask;
    isize = nx_fmask * ny_fmask * sizeof(float);
    fmask = (float *)malloc(isize);
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
    COVERA_MASK(fmask,&nx_fmask,&ny_fmask,&ir,&max_nclosure,&nbeta,&ngamma);
/* Free virtual memory space: */
    free(fmask);

  fprintf(fp1," Radius of uv-coverage (IR) in pixels: %d\n",ir);
  fprintf(fp1," nbeta: %d ngamma: %d\n",nbeta,ngamma);
#ifdef DEBUG
  printf(" nbeta: %d ngamma: %d\n",nbeta,ngamma);
#endif

isize = nx * ny * sizeof(double);
  image1 = (double *)malloc(isize);
  image2 = (double *)malloc(isize);
  mim = (double *)malloc(isize);
  long_int = (double *)malloc(isize);
  modsq = (double *)malloc(isize);
  snrm = (double *)malloc(isize);
isize = 4 * ngamma * sizeof(double);
  bispp = (double *)malloc(isize);

/* Erasing the arrays : */
    for(i = 0; i < nx * ny; i++) {
      image1[i] = 0.;
      image2[i] = 0.;
      mim[i] = 0.;
      long_int[i] = 0.;
      modsq[i] = 0.;
      snrm[i] = 0.;
      }
    for(i=0; i < 4 * ngamma; i++) { bispp[i]=0.;}

/******************************************************************/
/*************** Flat field *************************************/
/* Computing flat-field */
#ifdef FFIELD_CORRECTION
  nx_ff = 1024; ny_ff = 1024;
  compute_ffield(&ffield,ff_name,&nx_ff,&ny_ff,
                  ixcent,iycent,(int)(nx*ireduc),(int)(ny*ireduc),
                  file_ext,fp1,0.3);
#endif

/* Since we process twice the data (with a shift of half a frame: */
nframes_expected = 2 * nphot_wanted / npack;

/* Long integration: 
  status = rdcar_long(long_int,nx,ny,in_name,nphot_wanted,&nphot_found);
*/
/* Photon processing: */
  status = rdcar_phot(image1,image2,mim,modsq,snrm,bispp,long_int,
               ffield,nx,ny,ixcent,iycent,ir,nbeta,ngamma,in_name,npack,
               nphot_wanted,&nphot_found,&nframes,nframes_expected,
               nx_ff,ireduc);
  printf(" Output from rdcar: %d photons and %d frames processed\n",
          nphot_found, nframes);
  fprintf(fp1," Output from rdcar: %d photons and %d frames processed\n",
          nphot_found, nframes);

/* Check that some frames have been processed */
   xframes=(float)nframes;

/* Avoid division by zero ... */
   if(xframes > 0. && modsq[0] > 0.)
     {
      printf(" Normalization: xframes=%.1f modsq[0]=%f \n",
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
     if(xframes <= 0) xframes = 1.;
      xphot = 1.; 
     }
  sprintf(outcomments,"%s %d-08-93 %.2fH %.1fs %.1fph %dfr %dph",
          in_name,(int)date,time,integ,xphot,nframes,nphot_found);

/****************************************************************/
   printf(" xphot = %f \n",xphot);
   fprintf(fp1," xphot = %f \n",xphot);

/* Normalizes FFT (for compatibility with FFT_2D instead of FFT_2D_FAST...*/
   jlp_normalize_fft(bispp,modsq,snrm,&nx,&ny,&nbeta,&ngamma);

/* Compute snrm, SNR of bispectrum, etc:  (xphot is not used) */
   prepare_output_bisp(bispp,modsq,snrm,nx,ny,xframes,xphot,
                         nbeta,ngamma,photon_correction,fp1);

/* Mean of the frames: */
   for(i = 0; i < nx*ny; i++) {
        long_int[i] /= xframes;
        }

/****************************************************************/
/* Now output the results : */

/* Long integration : */
  sprintf(outfile,"long%s",file_ext);
  fprintf(fp1," %s",outfile);
  JLP_D_WRITEIMAG(long_int,&nx,&ny,&nx,outfile,outcomments);

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

/* Mean squared modulus : */
   RECENT_FFT_DOUBLE(modsq,modsq,&nx,&ny,&nx);
   sprintf(outfile,"modsq%s",file_ext);
   JLP_D_WRITEIMAG(modsq,&nx,&ny,&nx,outfile,outcomments);

/* SNR of squared modulus (actually 1/sigma if no photon correction): */
   RECENT_FFT_DOUBLE(snrm,snrm,&nx,&ny,&nx);
   sprintf(outfile,"snrm%s",file_ext);
   JLP_D_WRITEIMAG(snrm,&nx,&ny,&nx,outfile,outcomments);

/* Bispectrum : */
   isize = 3 * ngamma * sizeof(float);
   float_array = (float *)malloc(isize);
   output_bisp(bispp,float_array,ngamma);
   free(bispp);
   sprintf(outfile,"bisp1%s",file_ext);
   ny=3;
   JLP_WRITEIMAG(float_array,&ngamma,&ny,&ngamma,outfile,outcomments);
   free(float_array);

   fprintf(fp1," Successful end \n");
   printf(" Successful end \n");

/* End : */
  JLP_END();
  printf(" Output logfile: %s\n",logfile);
fclose(fp1);
return(0);
}
/****************************************************************
*
Format CAR data 
Long integration only: 
(With debugging mode)
*
****************************************************************/
/* Subroutine rdcar_long to read CAR format */
int rdcar_long(double *long_int, INT4 nx, INT4 ny, char *in_name,
               INT4 nphot_wanted, INT4 *nphot_found)
{
FILE *fd;
int nbytes_to_read, nbytes, nvalues, nblock, nblock_max; 
int ix, iy, itime;
register int i, j;

union{
unsigned long ulg[256];
char ch[1024];
} buff;


/* Opens the input file */
if((fd = fopen(in_name,"r")) == NULL)
  {
  printf("rdcar/error opening input file: >%s< \n",in_name);
  return(-1);
  }

/***************************************************************/
/* Skip header (32 bytes containing "FORMYCAR"): */
  nbytes_to_read = 32;
  nbytes = fread(buff.ch,sizeof(char),nbytes_to_read,fd);
  if(nbytes != nbytes_to_read)
    {
     printf("rdcar/error skipping header: \n");
     printf("       Only %d bytes read \n", nbytes);
     return(-2);
    }

/***************************************************************/
/* Read the data: */
*nphot_found = 0;
nblock_max = 1 + nphot_wanted / 256;
if(nblock_max > NBLOCKMAX) nblock_max = NBLOCKMAX;
for(nblock = 0; nblock < nblock_max; nblock++)
{
/* Read next values */
  nvalues = fread(buff.ulg,sizeof(unsigned long),256,fd);
  if(nvalues != 256)
          printf(" rdcar/Warning, only %d values read in block #%d\n",
                   nvalues,nblock);
  if(nvalues <= 0 )
  {
  printf("rdcar/end of file: >%s< \n",in_name);
  nblock = NBLOCKMAX;
  }

/* Displays current value of nblock: */
#ifdef DEBUG
  if((nblock % 500) == 1)
     printf(" Block #%d  of 1024 bytes\n",nblock);
#else
  if((nblock % 1000) == 1)
     printf(" Block #%d  of 1024 bytes\n",nblock);
#endif

/* Reads photon coordinates and fills image array: */
/* Values are tyx, tyx, ... 
   with t on the first 12 bits, y on the next 10 bits and the x on 10 bits 
   msb first then lsb*/
for(i = 0; i < nvalues; i++) 
  {

#ifdef swap_like_dec
    inv_ulong_int(&buff.ulg[i]);
#endif
/* Use a mask: (since sometimes does not fill with zeroes) */
    itime = (buff.ulg[i] >> 20) & 4095;
    iy = (buff.ulg[i] >> 10) & 1023; 
    ix = buff.ulg[i] & 1023; 

/* Store photon at the coordinates location: */ 
   if(ix > 0 && ix < nx && iy > 0 && iy < ny)
           {
            j=ix + iy * nx;
/* Long integration: */
            long_int[j]++;
            (*nphot_found)++;
           }
/* End of loop on nvalues */ 
  }
/* End of current block: read next values */
}

/* Closes the input file */
  fclose(fd);
  return(0);
}
/****************************************************************
*
Format CAR data 
Photon processing
*
****************************************************************/
/* Subroutine rdcar_phot to read CAR format */
int rdcar_phot(double *image1, double *image2, double *mim, 
               double *modsq, double *snrm, double *bispp, double *long_int,
               float *ffield, INT4 nx, INT4 ny, 
               INT4 ixcent, INT4 iycent, INT4 ir, INT4 nbeta, INT4 ngamma, 
               char *in_name, INT4 npack, INT4 nphot_wanted, INT4 *nphot_found,
               INT4 *nframes, INT4 nframes_expected, INT4 idim_ff, INT4 ireduc)
{
FILE *fd;
int nbytes_to_read, nbytes, nvalues, nblock, nblock_max; 
int ix, iy, ix_start, iy_start, ix_end, iy_end, itime, iph1, iph2, iframe;
register int ii, j, j_ff;

union{
unsigned long ulg[256];
char ch[1024];
} buff;

/* Computes reduction factor: */
ix_start = ixcent - (nx * ireduc)/2; 
if(ix_start < 0) ix_start = 0;
iy_start = iycent - (ny * ireduc)/2; 
if(iy_start < 0) iy_start = 0;
ix_end = ix_start + nx * ireduc;
iy_end = iy_start + ny * ireduc;
#ifdef DEBUG
 printf("rdcar_phot/ ix_start = %d, iy_start=%d ",ix_start,iy_start);
 printf("ix_end = %d, iy_end=%d\n",ix_end,iy_end);
#endif

/* Opens the input file */
if((fd = fopen(in_name,"r")) == NULL)
  {
  printf("rdcar/error opening input file: >%s< \n",in_name);
  return(-1);
  }

/***************************************************************/
/* Skip header (32 bytes containing "FORMYCAR"): */
  nbytes_to_read = 32;
  nbytes = fread(buff.ch,sizeof(char),nbytes_to_read,fd);
  if(nbytes != nbytes_to_read)
    {
     printf("rdcar/error skipping header: \n");
     printf("       Only %d bytes read \n", nbytes);
     return(-2);
    }

/***************************************************************/
/* Read the data: */
*nphot_found = 0; iframe = 0; iph1 = 0; iph2 = 0;
nblock_max = 1 + nphot_wanted / 256;
if(nblock_max > NBLOCKMAX) nblock_max = NBLOCKMAX;
for(nblock = 0; nblock < nblock_max; nblock++)
{
/* Read next values */
  nvalues = fread(buff.ulg,sizeof(unsigned long),256,fd);
  if(nvalues != 256)
          printf(" rdcar/Warning, only %d values read in block #%d\n",
                   nvalues,nblock);
  if(nvalues <= 0 )
  {
  printf("rdcar/end of file: >%s< \n",in_name);
  nblock = NBLOCKMAX;
  }

/* Displays current value of nblock: */
#ifdef DEBUG
  if((nblock % 500) == 1)
     printf(" Block #%d  of 1024 bytes\n",nblock);
#else
  if((nblock % 1000) == 1)
     printf(" Block #%d  of 1024 bytes\n",nblock);
#endif

/* Reads photon coordinates and fills image array: */
/* Values are tyx, tyx, ... 
   with t on the first 12 bits, y on the next 10 bits and the x on 10 bits 
   msb first then lsb*/
for(ii = 0; ii < nvalues; ii++) 
  {

#ifdef swap_like_dec
    inv_ulong_int(&buff.ulg[ii]);
#endif
/* Use a mask: (since sometimes does not fill with zeroes) */
    itime = (buff.ulg[ii] >> 20) & 4095;
    iy = (buff.ulg[ii] >> 10) & 1023; 
    ix = buff.ulg[ii] & 1023; 

/* Compute coordinate in flat field file */ 
#ifdef FFIELD_CORRECTION
    j_ff = ix + iy * idim_ff;
#ifdef DEBUG
    if(j_ff < 0 || j_ff > idim_ff*idim_ff)
      printf(" DEBUG/error with ffield ..., ix=%d iy=%d \n",ix,iy);
#endif
#endif

/* Store photon at the coordinates location: */ 
   if(ix >= ix_start && ix < ix_end && iy >= iy_start && iy < iy_end)
           {
/* Division by "ireduc" to reduce size of output images */
/* JLP2001/ Warning, do not try to save time, and use all the following,
* otherwise round-off problems... */
            ix=(ix - ix_start)/ireduc;
            iy=(iy - iy_start)/ireduc;
            j = ix + iy*nx;
/* Building the image: */
#ifdef FFIELD_CORRECTION
            image1[j] += ffield[j_ff];
            image2[j] += ffield[j_ff];
#else
            (image1[j])++; (image2[j])++;
#endif
            (*nphot_found)++; iph1++; iph2++;
           }
/*****************************************************************/
/* Test to see if npack has been reached for the current pack1 
   (pack full, i.e. image1 with npack photons)
 */ 
    if(iph1 == npack)
    {
     if((iframe % 200) == 0) 
          printf(" rdcar_phot/Processing frame #%d/%d (pack1) \n",
                   iframe,nframes_expected);
     process_frame(image1, mim, modsq, snrm, bispp, long_int, 
                  nx, ny, ir, nbeta, ngamma);
    iph1 = 0; iframe++;
/* Reset image1 to zero: */
     for (j = 0; j < nx * ny; j++) image1[j] = 0.;
    }

/* Shift image2 by half a frame at the beginning : */
    if( iframe == 0 && iph1 == npack/2)
      {
       iph2 = 0;
/* Reset image2 to zero: */
       for (j = 0; j < nx * ny; j++) image2[j] = 0.;
      }

/* Test to see if npack has been reached for the current pack2 
   (pack full, i.e. image2 with npack photons)
 */ 
    if(iph2 == npack)
    {
     if((iframe % 200) == 101) 
          printf(" rdcar_phot/Processing frame #%d/%d (pack2) \n",
                   iframe,nframes_expected);
     process_frame(image2, mim, modsq, snrm, bispp, long_int, 
                  nx, ny, ir, nbeta, ngamma);
    iph2 = 0; iframe++;
/* Reset image to zero: */
     for (j = 0; j < nx * ny; j++) image2[j] = 0.;
    }

/* End of loop on ii (nvalues) */ 
  }
/* End of current block: read next values */
}

/* In case no frame has been processed,
   we want to recover the long integration at least... */
if(!iframe)
   {
     printf(" rdcar_phot/Before exiting from the program, save long integration\n");
     process_frame(image1, mim, modsq, snrm, bispp, long_int, 
                  nx, ny, ir, nbeta, ngamma);
     iframe = 1;
   }

/* Closes the input file */
  *nframes = iframe;
  fclose(fd);
  return(0);
}
/********************************************************
*
* Processing the data frame by frame: 
*********************************************************/
int process_frame(double *image, double *mim, double *modsq, 
                  double *snrm, double *bispp, double *long_int, 
                  INT4 nx, INT4 ny, INT4 ir, INT4 nbeta, INT4 ngamma)
{
register long int i;

/* Long integration : */
  for(i = 0; i < nx * ny; i++) long_int[i] += image[i];

/* Resetting the imaginary part for the FFT: */
  for(i = 0; i < nx * ny; i++) {mim[i] = 0.;}

/* Fourier Transform: (Cf ../fft/fftw_set.c)*/
  fftw_double(image,mim,(int)nx,(int)ny,1);
/* Please note that output arrays should be multiplied to sqrt(nx * ny) to
* obtain same results as FFT_2D ...*/

/* Processing this image now: bispec3 is without photon noise correction */
/* bispec3 in "jlp_cover_mask.c" */
   bispec3(image,mim,modsq,snrm,&nx,&ny,bispp,&ir,&nbeta,&ngamma);

return(0);
}
