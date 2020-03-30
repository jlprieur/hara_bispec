/**************************************************************** 
  decode_car1.c
 To decode data from the CAR photon counting camera 
 and compute bispectrum. 1-D only: spectrum computed along the lines.
 Same as decode_car2 but for 1-D spectra and 2-D bispectra.
 Allow for rotations of the initial spectrum.

 Photon noise correction. 
 Uses twice the data (interleaved by half a frame).
 Flat field correction (before compression)

 SYNTAX: 
 runs decode_car1 in_photon_file output_file_exten
   reduc_fact,ixcent,iycent,ir,max_nclosure,
   total_nphot,frame_nphot,photon_corr,rotation_angle
   flat_field Fourier_mask

 EXAMPLE:
 runs decode_car test tt 16,488,512,12,900,100,0,36.0 ff_oc45 mask1
 or without flat field and without mask:
  runs decode_car test tt 16,488,512,12,1000,100,0,36.0 0 0

 COMMENTS:
 Rotation angle is in degrees: it corresponds to the angle of the
 spectrum relative to the X axis (increasing pixels) and Y axis (increasing rows)
 The goal is to obtain a spectrum with emission lines on one row.

 JLP
 Version 09-08-95
*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <jlp_ftoc.h>

#undef DOUBLE_PRECISION

/*
#define DEBUG
*/

/*
#define APODIZATION
*/
#ifdef APODIZATION
/* Default is BLACKMAN */
/* When HAMMING is defined, uses Hamming instead of Blackman window */
#define HAMMING
#endif

/* Number of frames for intermediate output (in case of crash) */ 
#define NFRAMES_TO_SAVE 10000 

#define NBLOCKMAX 100000
/* Origin of date is 01-01-1904 at 0H TU
 which corresponds to Julian day 2416480.500 */
#define DATE_ORIGIN 2416480.50
/* Maximum size for the images */
#define IDIM 256
/* CAR detector has 1024x1024 pixels, but useful size is 864x864 : */
#define CAR_WIDTH 864 

/* Rotation matrix: */
static float c1, c2, c3, c4, c5, c6;
#define PI  3.14159

/* Work array for FFT: */
static double *nag_work;

main(argc, argv)
int argc;
char **argv;
{
/* Images : */
double *image1, *image2, *im, *modsq, *snrm, *long_int, *bispp;
float *modsq1, *snrm1, *long_int1;
float *bisp1, *in_ffield, *ffield, *apodi, *fmask;
float xframes, w1, w2, w3, date, time, integ, xphot; 
double rot_angle;
int nframes_expected, photon_correction;
int npack, nphot_wanted, nphot_header, nphot_found, nframes;
int status, ix_reduc, iy_reduc, ireduc_ff, isize, ixcent, iycent;
long int pntr_ima, nx, ny, nx_ff, ny_ff, nx_fmask, ny_fmask;
long int ir, max_nclosure, nbeta, ngamma, bisp_dim, iy_b;
register long int i, iy, ng;
char *pc, file_ext[41], buffer[81], logfile[41];
char fmask_name[61], fmask_comments[81];
char in_name[61], ff_name[61], ff_comments[81], outfile[61], outcomments[81];
FILE *fp1;

printf(" Program decode_car1  (%d x %d maxi) \n",IDIM,IDIM);
printf(" (Uses twice the data (interleaved by half a frame)) \n");
printf(" JLP Version 09-08-95 \n");


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
  printf(" runs decode_car1 in_photon_file output_file_exten");
  printf(" reduc_fact,ixcent,iycent,ir,max_nclosure,total_nphot,nphot_per_frame,photon_correction,rotation_angle ");
  printf(" flat_field Fourier_mask\n");
  printf(" (Rotation angle is in degrees): \n"); 
  printf(" Example: \n"); 
  printf("  runs decode_car test tt 16,488,512,12,900,10000,130,1,,30.5 ff_oc45 mask1 \n\n"); 
  printf(" or without flat field and without mask: \n");
  printf("  runs decode_car test tt 16,488,512,12,1000,13000,100,0,-12.4 0 0 \n\n"); 
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
sprintf(logfile,"decode_car1%s.log",file_ext);
if((fp1 = fopen(logfile,"w")) == NULL)
   {
   printf("decode_car1/Fatal error opening logfile >%s< \n",logfile);
   exit(-1);
   }

fprintf(fp1," Program decode_car1  (%d x %d maxi) \n",IDIM,IDIM);
fprintf(fp1," (Uses twice the data (interleaved by half a frame)) \n");
fprintf(fp1," JLP Version 09-08-95 \n");

/***************************************************************/
/* Read header first (for interactive input of the parameters...) */
  status = rdcar_header(in_name,&date,&integ,&nphot_header);
#ifdef DEBUG
  printf(" date = %f \n",date);
  printf(" (int)date = %d \n",(int)date);
#endif
  time = 24.*(date - (float)((int)date));
  printf(" Input file:%s date=%d-08-93 time=%.2fH itime=%.2f nphot=%d \n",
          in_name,(int)date,time,integ,nphot_header);
  fprintf(fp1," Input file:%s date=%d-08-93 time=%.2fH itime=%.2f nphot=%d \n",
          in_name,(int)date,time,integ,nphot_header);
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
  printf(" Select reduction factor of the output frames in X and Y, (4, 8, 16 or 32),");
  printf(" IXcenter, IYcenter, radius of uv coverage, max_nclosure, \n");
  printf(" total number of photons to be processed,");
  printf(" number of photons per frame, ");
  printf(" photon noise correction (1 for yes, 0 for no),\n");
  printf(" and rotation angle (degrees):\n");
  printf("     (Example: 16,488,512,12,1000,200000,250,1,12.5) \n");
  gets(buffer);sscanf(buffer,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf",
               &ix_reduc,&iy_reduc,&ixcent,&iycent,&ir,&max_nclosure,
               &nphot_wanted,&npack,&photon_correction,&rot_angle);
  }
else
  {
  sscanf(argv[3],"%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf",
               &ix_reduc,&iy_reduc,&ixcent,&iycent,&ir,&max_nclosure,
               &nphot_wanted,&npack,&photon_correction,&rot_angle);
  }

/* Logfile: */
  fprintf(fp1," Reduction factor of the output frames: ix_reduc=%d iy_reduc=%d\n",
               ix_reduc,iy_reduc);
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
/* Note that we have nx = 108, ny = 108, when ix_reduc = iy_reduc = 8  */
  nx = CAR_WIDTH / ix_reduc; ny = CAR_WIDTH / iy_reduc;
  if(nx <=0 || nx > IDIM || ny <=0 || ny > IDIM)
  {printf(" Fatal error for the output size: nx = %d, ny = %d\n",nx,ny);
  fprintf(fp1," Fatal error for the output size: nx = %d, ny = %d\n",nx,ny);
   exit(-1);}
/* Take even numbers (otherwise pb in "RECENT_FFT"..) */
  nx = (nx / 2) * 2; ny = (ny / 2) * 2;
  printf(" Output image size: nx = %d ny = %d \n",nx,ny);
  fprintf(fp1," Output image size: nx = %d ny = %d \n",nx,ny);

/* Rotation if needed: */
  fprintf(fp1," Rotation angle: %f \n",rot_angle);
/* Conversion to radians: */
  rot_angle *= (PI/180.);
  c2 = cos(rot_angle);  c3 = sin(rot_angle); 
  c5 = -sin(rot_angle); c6 = cos(rot_angle); 
/* Compute c1 and c4 to maintain center in the middle: */
  c1 = (float)ixcent * (1. - c2) - (float)iycent * c3; 
  c4 = (float)iycent * (1. - c6) - (float)ixcent * c5; 
  fprintf(fp1," Rotation matrix: c1=%f c2=%f c3=%f  \n",c1,c2,c3);
  fprintf(fp1,"                  c4=%f c5=%f c6=%f  \n",c4,c5,c6);

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
/* JLP prolog: */
JLP_BEGIN();
JLP_INQUIFMT();

/*****************************************************************/
/* Computing the spectral and bispectral lists
   corresponding to the selected uv coverage: */

  printf(" Radius of uv-coverage (IR) in pixels: %d\n",ir);
  printf(" Max number of closure relations: %d \n",max_nclosure);

  if(fmask_name[0] == '0')
    {
/* Full u-v coverage: */
    printf(" Full u-v coverage.\n");
    fprintf(fp1," Full u-v coverage.\n");
    nx_fmask = 2 * ir + 1; ny_fmask = 1;
    isize = nx_fmask * sizeof(float);
    JLP_GVM(&fmask, &isize);
    for(i = 0; i < nx_fmask; i++) fmask[i] = 1.; 
    }
  else
/* Read Fourier mask file: */
    {
    JLP_VM_READIMAG(&pntr_ima,&nx_fmask,&ny_fmask,fmask_name,fmask_comments);
    JLP_FROM_MADRID(&pntr_ima,&fmask);
    if( nx_fmask < (2 * ir + 1) || ny_fmask < 2)
      {
      fprintf(fp1,"Fatal error/ Wrong size for Fourier mask\n");
      fprintf(fp1," nx should be larger than to ir  and ny = 1\n");
      fprintf(fp1,"     (ir = %d, nx_fmask = %d, ny_fmask = %d)\n",
              ir, nx_fmask, ny_fmask);
      printf("Fatal error/ Wrong size for Fourier mask\n");
      printf(" nx should be larger than to ir  and ny = 1\n");
      printf("     (ir = %d, nx_fmask = %d, ny_fmask = %d)\n",
              ir, nx_fmask, ny_fmask);
      exit(-1);
      }
    fprintf(fp1," Masked u-v coverage. Mask= %s\n",fmask_name);
    }

/* "Masked" u-v coverage: */
    COVERA_MASK_1D(fmask,&nx_fmask,&ir,&max_nclosure,&nbeta,&ngamma);
/* Free virtual memory space: */
    JLP_FVM(&fmask);

  fprintf(fp1," Radius of uv-coverage (IR) in pixels: %d\n",ir);
  fprintf(fp1," nbeta: %d ngamma: %d\n",nbeta,ngamma);
#ifdef DEBUG
  printf(" nbeta: %d ngamma: %d\n",nbeta,ngamma);
#endif

isize = nx * sizeof(double);
  JLP_GVM(&nag_work,&isize);
isize = nx * ny * sizeof(double);
  JLP_GVM(&image1,&isize);
  JLP_GVM(&image2,&isize);
  JLP_GVM(&im,&isize);
  JLP_GVM(&modsq,&isize);
  JLP_GVM(&snrm,&isize);
  JLP_GVM(&long_int,&isize);
isize = 4 * ngamma * ny * sizeof(double);
  JLP_GVM(&bispp,&isize);

/* Erasing the arrays : */
    for(i = 0; i < nx * ny; i++) {
      image1[i] = 0.;
      image2[i] = 0.;
      im[i] = 0.;
      long_int[i] = 0.;
      modsq[i] = 0.;
      snrm[i] = 0.;
      }
    for(i=0; i < 4 * ngamma * ny; i++) { bispp[i]=0.;}

/******************************************************************/
/*************** Flat field *************************************/
if(ff_name[0] != '0')
 {

/* Read flat field file: */
  JLP_VM_READIMAG(&pntr_ima,&nx_ff,&ny_ff,ff_name,ff_comments);
  JLP_FROM_MADRID(&pntr_ima,&in_ffield);

/* Check input size: */
  if(nx_ff != ny_ff)
  {printf(" Fatal error/Ffield has wrong size: it should be square!\n");
  fprintf(fp1," Fatal error/Ffield has wrong size: it should be square!\n");
   exit(-1);
  }
/* Compute reduction of the flat field: */
  ireduc_ff = 1024 / nx_ff;
  if(nx_ff * ireduc_ff != 1024)
  {printf(" Fatal error/Ffield has wrong size: it should be a submultiple of the 1024x1024 detector format)!\n");
  fprintf(fp1," Fatal error/Ffield has wrong size: it should be a submultiple of the 1024x1024 detector format)!\n");
   exit(-1);
  }

/* Allocate memory space for flat field: */
 isize = nx_ff * ny_ff *sizeof(float);
 JLP_GVM(&ffield,&isize);

/****************************************************************/
/* Displays message if apodization: */
#ifdef APODIZATION
     printf(" Flat field apodization\n");
     fprintf(fp1," Flat field apodization\n");
#else
     printf(" No flat field apodization\n");
     fprintf(fp1," No flat field apodization\n");
#endif

#ifdef APODIZATION
/* Creates Hamming apodization file: */
  JLP_GVM(&apodi,&isize);
#ifdef HAMMING
  jlp_hamming(apodi,nx_ff,ny_ff,nx_ff);
#else
  jlp_blackman(apodi,nx_ff,ny_ff,nx_ff);
#endif

#ifdef DEBUG
#ifdef HAMMING
   strcpy(outfile,"hamming");
   strcpy(outcomments,"Hamming filter: 0.54 + 0.46 cos(pi t / T)");
#else
   strcpy(outfile,"blackman");
   strcpy(outcomments,"Blackman filter: 0.42 + 0.5 cos(pi t / T) + 0.08 cos(2pit/T)");
#endif
   JLP_WRITEIMAG(apodi,&nx_ff,&ny_ff,&nx_ff,outfile,outcomments);
#endif

/* End of case (#ifdef APODIZATION) */
#endif
/****************************************************************/

  create_ffield(ffield,in_ffield,apodi,nx_ff,ny_ff,ireduc_ff,0.6);

  printf(" Ffield successfully built \n");
  sprintf(outfile,"ff%s",file_ext);
  strcpy(outcomments,"inverse ffield reduced");
  JLP_WRITEIMAG(ffield,&nx_ff,&ny_ff,&nx_ff,outfile,outcomments);

/* Free virtual memory space: */
  JLP_FVM(&in_ffield);
#ifdef APODIZATION
  JLP_FVM(&apodi);
#endif

}
/********************* Unity file for ffield: *************/
else
{
/* Should be a power of two, to avoid overflow problems in rdcar_long... */ 
/* JLP94: I set it to 4 since I had problems with 1024/nx... */
 ireduc_ff = 4;
 nx_ff = 1024 / ireduc_ff; ny_ff = nx_ff;
 isize = nx_ff * ny_ff *sizeof(float);
 JLP_GVM(&ffield,&isize);
 for(i = 0; i < nx_ff * ny_ff; i++) ffield[i] = 1.;
 fprintf(fp1,"\n OK, Take unity file for ffield \n");
 printf("\n OK, Take unity file for ffield \n");
}

/* Since we process twice the data (with a shift of half a frame: */
nframes_expected = 2 * nphot_wanted / npack;

/* Long integration: 
  status = rdcar_long(long_int,nx,ny,nx,in_name,nphot_wanted,&nphot_found);
*/
/* Photon processing: */
  status = rdcar_phot(image1,image2,im,modsq,snrm,bispp,long_int,
               ffield,nx,ny,nx,ixcent,iycent,ir,nbeta,ngamma,in_name,npack,
               nphot_wanted,&nphot_found,&nframes,nframes_expected,
               nx_ff,ireduc_ff);
  printf(" Output from rdcar: %d photons and %d frames processed\n",
          nphot_found, nframes);
  fprintf(fp1," Output from rdcar: %d photons and %d frames processed\n",
          nphot_found, nframes);

/* Check that some frames have been processed */
   xframes=(float)nframes;

/* Avoid division by zero ... */
   xphot = (double)npack; 
   if(xframes < 1.) xframes = 1.;
   printf(" Normalization: xframes=%f xphot=%f \n",xframes,xphot);
   fprintf(fp1," Normalization: xframes=%f xphot=%f \n",xframes,xphot);

  sprintf(outcomments,"%s %d-08-93 %.2fH %.1fs %.1fph %dfr %dph",
          in_name,(int)date,time,integ,xphot,nframes,nphot_found);

/* Mean of the frames: */
   for(i = 0; i < nx * ny; i++) {
        long_int[i] /= xframes;
        modsq[i] /= xframes;
        snrm[i] /= xframes;
        }
   for(i=0; i< 4 * ngamma * ny; i++) bispp[i] /= xframes;
/****************************************************************/
/* SNR of modsq. First step: sigma**2 computation 
* (Keep sigma**2, for consistency with normalisation in photon_corr)
*  AFTER MEAN COMPUTATION AND BEFORE PHOTON NOISE CORRECTION!!! */
   for(i = 0; i < nx * ny; i++) {
        snrm[i] = snrm[i] - modsq[i]*modsq[i];
        }

/****************************************************************/
/* Normalizes FFT (for compatibility with FFT_2D instead of FFT_2D_FAST...*/
   jlp_normalize_fft_1D(bispp,modsq,snrm,&nx,&ny,&nbeta,&ngamma);
/* Option removed: */
/*
   if(photon_correction) 
          photon_corr(bispp,modsq,snrm,&nx,&ny,&xphot,&nbeta,&ngamma);
*/

/****************************************************************/
/* SNR of modsq. Second step: modsq[i]/sqrt(variance) 
*  AFTER PHOTON NOISE CORRECTION!!! */
   for(i = 0; i < nx * ny; i++) 
      {
      if(snrm[i] <= 1.e-4) snrm[i]=1.e-4;
      snrm[i] = 1. / sqrt((double)snrm[i]);
      }
/* JLP96: */
     for(i = 0; i < nx * ny; i++) snrm[i] *= modsq[i];

/***************** Bispectrum: **********************************/
   bisp_dim = 4 * ngamma;
   for(iy=0; iy<ny; iy++) 
     {
     iy_b = iy * bisp_dim;
     for(ng=0; ng<ngamma; ng++) 
       {
/* First computing the variance: */
        w1 =  bispp[ng + 2*ngamma + iy_b] - bispp[ng + iy_b]*bispp[ng + iy_b];
        w2 =  bispp[ng + 3*ngamma + iy_b] 
                    - bispp[ng + ngamma + iy_b]*bispp[ng + ngamma + iy_b];
/* Then the sigma (real and imag together): */
        w1 = w1 + w2;
        if(w1 < 1.e-10) w1 = 1.e-10;
        w1 = sqrt((double)w1);
/* Phase factor of the bispectrum */
        w3 = bispp[ng + iy_b]*bispp[ng + iy_b] 
              + bispp[ng + ngamma + iy_b]*bispp[ng + ngamma + iy_b];
        if(w3 < 1.e-10) w3 = 1.e-10;
        w3 = sqrt((double)w3);
/* Normalization: 
* Commented out...
        bispp[ng + iy_b] /= w3;
        bispp[ng + ngamma + iy_b] /= w3;
*/
/* SNR of bispectrum in 3rd line: */
        bispp[ng + 2*ngamma + iy_b] = w3/w1;
        }
    }

/****************************************************************/
/* Now output the results : */

/* Go from double to single precision: */
/* Contraction from 4 to 3 parameters for the bispectrum: */
   isize = 3 * ngamma * ny * sizeof(double);
   JLP_GVM(&bisp1,&isize);
   output_bisp_1D(bispp,bisp1,ngamma,ny,bisp_dim);
   JLP_FVM(&bispp);
/* nx * ny arrays: */
   isize = nx * ny * sizeof(float);
   JLP_GVM(&modsq1,&isize);
   to_single(modsq,modsq1,&nx,&ny,&nx);
   JLP_FVM(&modsq);
   JLP_GVM(&snrm1,&isize);
   to_single(snrm,snrm1,&nx,&ny,&nx);
   JLP_FVM(&snrm);
   JLP_GVM(&long_int1,&isize);
   to_single(long_int,long_int1,&nx,&ny,&nx);
   JLP_FVM(&long_int);

/* Recentre the frames: */
   RECENT_FFT_1D(modsq1,modsq1,&nx,&ny,&nx);
   RECENT_FFT_1D(snrm1,snrm1,&nx,&ny,&nx);

/* Mean squared modulus : */
   sprintf(outfile,"modsq%s",file_ext);
   fprintf(fp1," Output of %s",outfile);
   JLP_WRITEIMAG(modsq1,&nx,&ny,&nx,outfile,outcomments);

/* SNR of squared modulus : */
   sprintf(outfile,"snrm%s",file_ext);
   fprintf(fp1," %s",outfile);
   JLP_WRITEIMAG(snrm1,&nx,&ny,&nx,outfile,outcomments);

/* Long integration : */
   sprintf(outfile,"long%s",file_ext);
   fprintf(fp1," %s",outfile);
   JLP_WRITEIMAG(long_int1,&nx,&ny,&nx,outfile,outcomments);


/* Bispectrum : */
   sprintf(outfile,"bisp1%s",file_ext);
   isize = 3*ngamma;
   fprintf(fp1," %s\n",outfile);
   JLP_WRITEIMAG(bisp1,&isize,&ny,&isize,outfile,outcomments);

   fprintf(fp1," Successful end \n");
   printf(" Successful end \n");

/* End : */
end1:
  JLP_END();
  printf(" Output logfile: %s\n",logfile);
fclose(fp1);
}
/***********************************************************************/
/****************************************************************/
/*
Format CAR : read header 
*/
/****************************************************************/
/* Subroutine rdcar_header to read CAR format */
int rdcar_header(in_name,date,integ,nphot)
char in_name[];
float *date, *integ;
int *nphot;
{
FILE *fd;
int nbytes_to_read, nbytes, nvalues; 
char cbuff[9];
unsigned long s_date, integ_time, nphot1;

/* Header of 32 bytes */
typedef struct {
unsigned long integ_time;          /* duree en msec */
unsigned long date;                /* date of observation 
                             (in seconds starting from 01-01-1904 at 0H TU)
                             (i.e. Julian day 2416480.500) */
unsigned long nber_of_photons;     /* number of photons */
long keyword1;         /* "FORM" in ASCII */
long keyword2;         /* "YCAR" in ASCII */
long nber_of_images;      /* used only for the CP40, here always 0 */
short refNum;             /* always 0 */
long read_already;        /* not used */ 
short everything_read;    /* not used */
} CAR_HEADER;

union{
long lg[256];
char ch[1024];
} buff;

CAR_HEADER *chead;
void swap_lint(), inv_lint();

#ifdef DEBUG
printf(" \n rdcar/reading file : %s \n",in_name);
#endif

/* Opens the input file */
if((fd = fopen(in_name,"r")) == NULL)
  {
  printf("rdcar/error opening input file: >%s< \n",in_name);
  return(-1);
  }

/***************************************************************/
/* Read header (32 bytes containing "FORMYCAR"): */
  nbytes_to_read = 32;
  nbytes = fread(buff.ch,sizeof(char),nbytes_to_read,fd);
#ifdef DEBUG
   printf("        %d bytes read \n", nbytes);
#endif
   if(nbytes != nbytes_to_read)
    {
     printf("rdcar/error reading header: \n");
     printf("       Only %d bytes read \n", nbytes);
     return(-2);
    }

/* Decode header: check that 'FORMYCAR' is at the right spot: */
  strncpy(cbuff,&buff.ch[12],8);
  cbuff[8]='\0';
  if(strcmp(cbuff,"FORMYCAR") != 0)
     {
      printf("rdcar/Fatal error, wrong header: >%s< \n",cbuff);
      return(-1);
     }
  else
/* OK: Successfull research: */
     {
      chead = (CAR_HEADER *)&(buff.ch[0]);
     }

/*************** date ****************************************/
   s_date = chead->date;
/* When working with DEC computer, reverses long integers... */
#ifdef dec
   inv_lint(&s_date);
#endif
/* s_date was in seconds, I convert it to days: */
   *date = (float)s_date / 86400.;
/* Then I subtract 2449199.500 which corresponds to the 31-07-93 */
    *date = *date - (2449199.500 - DATE_ORIGIN);
   printf(" date =%.4f-08-93\n",*date);

/*************** integration time *******************************/
   integ_time = chead->integ_time;
#ifdef dec
   inv_lint(&integ_time);
#endif
/* integration time was in milliseconds, I convert it to seconds: */
   *integ = (float)integ_time / 1000.;
   printf(" integration time =%f (sec) \n",*integ);

/*************** number of photons *******************************/
   nphot1 =chead->nber_of_photons;
#ifdef dec
   inv_lint(&nphot1);
#endif
   *nphot = nphot1;
   printf(" nphot =%u (photons)\n",*nphot);

/* Closes the input file */
  fclose(fd);
  return(0);
}
/**************************************************************
* Swap two bytes of an unsigned short integer
***************************************************************/
void swap_int(i)
unsigned short *i;
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
***************************************************************/
void swap_lint(i)
unsigned long *i;
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
***************************************************************/
void inv_lint(i)
unsigned long *i;
{
union {
unsigned long ii;
char         ch[4];
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
*
Format CAR data 
Long integration only: 
(With debugging mode)
*
****************************************************************/
/* Subroutine rdcar_long to read CAR format */
int rdcar_long(long_int,nx,ny,idim,in_name,nphot_wanted,nphot_found)
double long_int[];
long int nx, ny, idim;
int *nphot_found, nphot_wanted;
char in_name[];
{
FILE *fd;
int nbytes_to_read, nbytes, nvalues, nblock, nblock_max; 
int ix, iy, itime, ix_rot, iy_rot;
register int i, j, k;

union{
long lg[256];
char ch[1024];
} buff;

void swap_lint(), inv_lint();

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
  nvalues = fread(buff.lg,sizeof(long),256,fd);
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

#ifdef DEBUG
     if(i < 5 && nblock < 1)
       {
       printf(" i=%d array=%x (Hexa) ",i,buff.lg[i]);
       }
#endif

#ifdef dec
    inv_lint(&buff.lg[i]);
#endif
/* Use a mask: (since sometimes does not fill with zeroes) */
    itime = (buff.lg[i] >> 20) & 4095;
    iy = (buff.lg[i] >> 10) & 1023; 
    ix = buff.lg[i] & 1023; 

/* WARNING: first rotation, then compression... (otherwise artefacts...)*/
/* If rotation, compute new position: */
/* Add 0.5 to avoid steps... */
/* WARNING: should use new variables, since ix is modified in first line... */
    ix_rot = (int)(0.5 + c1 + c2 * (float)ix + c3 * (float)iy);
    iy_rot = (int)(0.5 + c4 + c5 * (float)ix + c6 * (float)iy);

/* Division by 4 to reduce size of output images
*/
    ix = (ix_rot >> 2) & 255; iy = (iy_rot >> 2) & 255;

#ifdef DEBUG
     if(i < 5 && nblock < 1)
       {
       printf(" i=%d array=%x (Hexa) ",i,buff.lg[i]);
       printf("   itime=%x ix=%x iy=%x (Hexa) \n",itime,ix,iy);
       printf("   time=%d x = %d, y=%d (dec) \n",itime,ix,iy);
       }
#endif

/* Store photon at the coordinates location: */ 
   if(ix > 0 && ix < nx && iy > 0 && iy < ny)
           {
            j=ix + iy * idim;
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
int rdcar_phot(image1,image2,im,modsq,snrm,bispp,long_int,
               ffield,nx,ny,idim,ixcent,iycent,ir,nbeta,ngamma,in_name,npack,
               nphot_wanted,nphot_found,nframes,nframes_expected,
               idim_ff,ireduc_ff)
double image1[], image2[], im[];
double modsq[], snrm[], long_int[], bispp[];
float ffield[]; 
long nx, ny, idim, ir, nbeta, ngamma;
int ixcent, iycent, idim_ff, ireduc_ff;
int npack, *nphot_found, nphot_wanted, *nframes, nframes_expected;
char in_name[];
{
FILE *fd;
int nbytes_to_read, nbytes, nvalues, nblock, nblock_max; 
int ix, iy, ix_rot, iy_rot, ixstart, iystart, ixend, iyend; 
int itime, iph1, iph2, iframe, ix_reduc, iy_reduc, x_roundoff, y_roundoff;
register int ii, i, j, j_ff, k;

union{
long lg[256];
char ch[1024];
} buff;

void swap_lint(), inv_lint();

/* Computes reduction factors: */
ix_reduc = CAR_WIDTH / nx;
iy_reduc = CAR_WIDTH / ny;
x_roundoff = ix_reduc / 2;
y_roundoff = iy_reduc / 2;
ixstart = (ixcent / ix_reduc) - nx/2; 
iystart = (iycent / iy_reduc) - ny/2; 
ixend = ixstart + nx;
iyend = iystart + ny;
#ifdef DEBUG
 printf("rdcar_phot/ ixstart = %d, iystart=%d ",ixstart,iystart);
 printf("ixend = %d, iyend=%d\n",ixend,iyend);
#endif
if(ixstart < 0 || iystart < 0 ) 
 {printf("rdcar_phot/Fatal error, ixstart = %d, iystart=%d ",ixstart,iystart);
  exit(-1);
 }

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
  nvalues = fread(buff.lg,sizeof(long),256,fd);
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

#ifdef dec
    inv_lint(&buff.lg[ii]);
#endif
/* Use a mask: (since sometimes does not fill with zeroes) */
    itime = (buff.lg[ii] >> 20) & 4095;
    iy = (buff.lg[ii] >> 10) & 1023; 
    ix = buff.lg[ii] & 1023; 

/* Compute coordinate in flat field file */ 
/* JLP95: ffield removed: 
    j_ff = (ix / ireduc_ff) + (iy / ireduc_ff) * idim_ff;
#ifdef DEBUG
    if(j_ff < 0 || j_ff > idim_ff* idim_ff)
      printf(" DEBUG/error with ffield ..., ix=%d iy=%d \n",ix,iy);
#endif
*/

/* WARNING: first rotation, then compression... (otherwise artefacts...)*/
/* If rotation, compute new position: */
/* Add 0.5 to avoid steps... */
/* WARNING: should use new variables, since ix is modified in first line... */
    ix_rot = (int)(0.5 + c1 + c2 * (float)ix + c3 * (float)iy);
    iy_rot = (int)(0.5 + c4 + c5 * (float)ix + c6 * (float)iy);

/* Division by "ireduc" to reduce size of output images
*/
    ix = (ix_rot + x_roundoff) / ix_reduc; iy = (iy_rot + y_roundoff) / iy_reduc;

/* Store photon at the coordinates location: */ 
   if(ix >= ixstart && ix < ixend && iy >= iystart && iy < iyend)
           {
            j=(ix - ixstart) + (iy - iystart) * idim;
/* Building the image: */
/* JLP95: FField removed:
            image1[j] += ffield[j_ff];
            image2[j] += ffield[j_ff];
*/
            image1[j] ++; image2[j] ++;
            (*nphot_found)++; iph1++; iph2++;
           }
/*****************************************************************/
/* Test to see if npack has been reached for the current pack1 
   (pack full, i.e. image1 with npack photons)
 */ 
    if(iph1 == npack)
    {
#ifndef DEBUG
     if((iframe % 200) == 0) 
#endif
     printf(" rdcar_phot/Processing frame #%d/%d (pack1) \n",
                 iframe,nframes_expected);
     process_frame(image1, im, modsq, snrm, bispp, long_int, 
                  nx, ny, idim, ir, nbeta, ngamma);
    iph1 = 0; iframe++;
/* Reset image1 to zero: */
     for (j = 0; j < nx * ny; j++) image1[j] = 0.;
    }

/* Shift image2 by half a frame: */
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
#ifndef DEBUG
     if((iframe % 200) == 101) 
#endif
          printf(" rdcar_phot/Processing frame #%d/%d (pack2) \n",
                   iframe,nframes_expected);
     process_frame(image2, im, modsq, snrm, bispp, long_int, 
                  nx, ny, idim, ir, nbeta, ngamma);
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
     process_frame(image1, im, modsq, snrm, bispp, long_int, 
                  nx, ny, idim, ir, nbeta, ngamma);
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
int process_frame(image, im, modsq, snrm, bispp, long_int, 
                  nx, ny, idim, ir, nbeta, ngamma)
double image[], im[], modsq[], snrm[], long_int[], bispp[];
long nx, ny, idim, ir, nbeta, ngamma;
{
float w1;
register int i, iy;

/* Long integration : */
  for(i = 0; i < nx * ny; i++) long_int[i] += image[i];

/* Resetting the imaginary part for the FFT: */
  for(i = 0; i < nx * ny; i++) {im[i] = 0.;}

/* Fourier Transform: */
  for(iy = 0; iy < ny; iy++) 
     FFT_1D_FAST(&image[iy*nx],&im[iy*nx],&nx,nag_work);

/* Please note that output arrays should be multiplied to sqrt(nx) to
* obtain same results as FFT_2D ...*/

/* Processing this image now: bispec_1D is without photon noise correction */
   bispec_1D(image,im,modsq,snrm,&nx,&ny,bispp,&ir,&nbeta,&ngamma);

return(0);
}
/****************************************************************/
int create_ffield(ffield,in_ffield,apodi,nx_ff,ny_ff,ireduc_ff,mini_value)
float ffield[], in_ffield[], apodi[],mini_value;
int nx_ff, ny_ff, ireduc_ff;
{
register int j;
int nvalues, nxy;
double sum;
float ff_mean;

/* Transfer to new array: */
nxy = nx_ff * ny_ff;
  for(j = 0; j < nxy; j++) ffield[j] = in_ffield[j];

/* First normalization: (take into account only non zero values) */
sum = 0.; nvalues=0;
for(j = 0; j < nxy; j++) 
    {if(ffield[j] > 0)
       {
       sum += ffield[j];
       nvalues++;
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
   apodization. As this means that we are anyway outside the 
   sensitive zone, we set the maximum to that value (or 1/value here): */ 
for(j = 0; j < nxy; j++) 
     if(ffield[j] < mini_value) 
         ffield[j] = 1./mini_value;
     else 
         ffield[j] = 1./ffield[j];

/* Apodization: */
#ifdef APODIZATION
for(j = 0; j < nxy; j++) ffield[j] *= apodi[j];
#endif

/* Second normalization: (take into account only values close to mean)*/
sum = 0.; nvalues=0;
for(j = 0; j < nxy; j++) 
    {if(ffield[j] < 2. && ffield[j] > 0.5)
       {
       sum += ffield[j];
       nvalues++;
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
**************************************************************/
int jlp_hamming(apodi,nx,ny,idim)
int nx, ny,idim;
float *apodi;
{
double argx, w1, radius, radmin, radmax, width;
int nx1, ny1, jrad;
register int i, j;

/* We take a subwindow to fit better the shape of the
   sensitive part of the camera: */
for(j = 0; j < ny; j++)
   for(i = 0; i < nx; i++)
    apodi[i + j * idim] = 1.;

width = (double)nx / 6.;
radmin = (double)(nx/2) - width;
radmax = (double)(nx/2);
w1 = PI / (double)width;
/* Works with circular filter: */
for(j = 0; j < ny; j++)
 {
   jrad = (j - ny/2)*(j - ny/2);
   for(i = 0; i < nx; i++) 
    {
     radius = jrad + (i - nx/2) * (i - nx/2);
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
int jlp_hamming1(apodi,nx,ny,idim)
int nx, ny,idim;
float *apodi;
{
double argx, argy, wx1, wy1;
float apodi_y;
int istart, jstart, nx1, ny1;
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
int jlp_hamming2(apodi,nx,ny,idim)
int nx, ny,idim;
float *apodi;
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
/***************************************************************/
/* Creates Blackman apodization file: */
int jlp_blackman(apodi,nx,ny,idim)
int nx, ny,idim;
float *apodi;
{
double argx, argy;
double wx1, wy1;
float apodi_y;
register int i, j;
int istart, jstart, nx1, ny1;

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
/***********************************************************
* output_bisp_1D
* Prepare output of bispectrum
*
* Go from double to single precision
* and contraction from 4 to 3 parameters
***********************************************************/
int output_bisp_1D(bispp, bisp1, ngamma, ny, bisp_dim)
double *bispp;
float *bisp1;
int ngamma, ny, bisp_dim;
{
register int ng, iy, k;

for(iy = 0; iy < ny; iy++)
   {
   for(ng = 0; ng < ngamma; ng++)
      {
      for(k = 0; k < 3; k++)
        {
        bisp1[ng + k*ngamma + iy*ngamma*3] = 
                              bispp[ng + k*ngamma + iy*bisp_dim];
        }
      }
   }
return(0);
}
