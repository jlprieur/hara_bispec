/**************************************************************** 
  dec_car4.c
 To decode data from the CAR photon counting camera 
 and compute long integration and autocorrelation
 Fast version of decode_car2 which computes FFT's and bispectrum. 

 Smaller files for autocorrelation, intercorrelation.

 Flat field correction

 dec_car3.c (now called autoc_car) is cpu optimized, compared to this version.

 JLP
 Version 02-09-93
*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <jlp_ftoc.h>

/*
#define DEBUG
#define FFIELD_CORRECTION 
*/

#define APODIZATION
#ifdef APODIZATION
/* When HAMMING is defined, uses Hamming instead of Blackman window */
#define HAMMING
#endif

/* Block size (unit is long integer length) */
#define BSIZE 8192 
#define NBLOCKMAX 100000
/* Origin of date is 01-01-1904 at 0H TU
 which corresponds to Julian day 2416480.500 */
#define DATE_ORIGIN 2416480.50
/* Maximum size for the images */
#define IDIM 256
/* Size (in x) for the autocorrelation and intercorrelation */
#define DIM_AUTO 64 
/* CAR detector has 1024x1024 pixels : */
#define CAR_WIDTH 1024 

#ifdef ibm
#define COVERA covera
#define FFT_2D fft_2d
#define BISPEC1 bispec1
#define BISPEC2 bispec2
#else
#define COVERA covera_
#define FFT_2D fft_2d_
#define BISPEC1 bispec1_
#define BISPEC2 bispec2_
#endif

main(argc, argv)
int argc;
char **argv;
{
/* Images : */
float *autoco, *interco, *long_int;
float *val_1;
long int *ix_1, *iy_1;
float *in_ffield, *ffield, *apodi;
float w1, w2, date, time, integ;
int npack, nphot_wanted, nphot_header, nphot_found;
int status, ireduc, ireduc_ff, isize, ixcent, iycent;
long int pntr_ffield, nx, ny, nx_auto, ny_auto, nx_ff, ny_ff;
register int i;
char *pc, file_ext[41], buffer[81], logfile[41];
char in_name[61], ff_name[61], ff_comments[81], outfile[61], outcomments[81];
FILE *fp1;

printf(" Program dec_car4  (%d x %d maxi) \n",IDIM,IDIM);
printf(" JLP Version 02-10-93 \n");

/*
test_machine();
*/

/* One or three parameters only are allowed to run the program: */
/* Carefull: 7 parameters always, using JLP "runs" */
#ifdef DEBUG
  printf(" argc=%d\n",argc);
  printf(" argv[3]=>%s<\n",argv[3]);
#endif
if((argc > 1 && argc < 5) || (argc == 7 && !strcmp(argv[3],"")))
  {
  printf("        Fatal error: Wrong syntax, argc=%d\n",argc);
  printf(" Syntax is:  \n");
  printf(" runs dec_car3 in_photon_file output_file_exten");
  printf(" reduc_fact,total_nphot flat_field\n");
  printf(" Example: runs dec_car3 test tt 16,1000,100 ff_oc45 \n\n"); 
  printf(" or without flat field: runs dec_car3 test tt 16,1000,100 0 \n\n"); 
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
sprintf(logfile,"dec_car4%s.log",file_ext);
if((fp1 = fopen(logfile,"w")) == NULL)
   {
   printf("dec_car4/Fatal error opening logfile >%s< \n",logfile);
   exit(-1);
   }

fprintf(fp1," Program dec_car4  (%d x %d maxi) \n",IDIM,IDIM);
fprintf(fp1," Sliding window method\n");
fprintf(fp1," JLP Version 09-10-93 \n");

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
  printf(" Select reduction factor of the output frames, (4, 8, 16 or 32),");
  printf(" total number of photons to be processed,");
  printf(" number of photons per frame \n");
  printf("     (Example: 16,2000,250) \n");
  gets(buffer);sscanf(buffer,"%d,%d,%d",
               &ireduc,&nphot_wanted,&npack);
  }
else
  {
  sscanf(argv[3],"%d,%d,%d",
               &ireduc,&nphot_wanted,&npack);
  }

/* Logfile: */
  fprintf(fp1," Reduction factor of the output frames: %d \n",ireduc);
  fprintf(fp1," IXcenter=%d, IYcenter=%d\n",
           ixcent,iycent);
  fprintf(fp1," Total number of photons to be processed: %d\n",nphot_wanted);
  fprintf(fp1," Number of photons per frame: %d \n",npack);

/* Just to check (for the user): */ 
#ifdef APODIZATION
     printf(" Flat field apodization\n");
     fprintf(fp1," Flat field apodization\n");
#else
     printf(" No flat field apodization\n");
     fprintf(fp1," No flat field apodization\n");
#endif

/* Set size of output images */
/* Note that we have nx = 128, ny = 128, when ireduc = 8  */
  nx = CAR_WIDTH / ireduc; ny = CAR_WIDTH / ireduc;
  if(nx <=0 || nx > IDIM || ny <=0 || ny > IDIM)
  {printf(" Fatal error for the output size: nx = %d, ny = %d\n",nx,ny);
  fprintf(fp1," Fatal error for the output size: nx = %d, ny = %d\n",nx,ny);
   exit(-1);}
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
 }
else
 {
  strcpy(ff_name,argv[4]);
 }

/**********************************************************/
JLP_BEGIN();
JLP_INQUIFMT();


/*****************************************************************/

isize = nx * ny * sizeof(float);
  JLP_GVM(&long_int,&isize);

/* Autocorrelation is reduced to the upper half (y >=0) */ 
nx_auto = DIM_AUTO;
ny_auto = DIM_AUTO/2;
isize = nx_auto * ny_auto * sizeof(float);
  JLP_GVM(&autoco,&isize);
  JLP_GVM(&interco,&isize);

isize = BSIZE * sizeof(long);
  JLP_GVM(&ix_1,&isize);
  JLP_GVM(&iy_1,&isize);

#ifdef FFIELD_CORRECTION
isize = BSIZE * sizeof(float);
  JLP_GVM(&val_1,&isize);
    for(i = 0; i < BSIZE; i++) val_1[i] = 0.;
#endif

/* Erasing the arrays : */
    for(i = 0; i < BSIZE; i++) {
      ix_1[i] = 0;
      iy_1[i] = 0;
      }
    for(i = 0; i < nx_auto * ny_auto; i++) {
      autoco[i] = 0.;
      interco[i] = 0.;
      }
    for(i = 0; i < nx * ny; i++) {
      long_int[i] = 0.;
      }

/******************************************************************/
#ifdef FFIELD_CORRECTION
/*************** Flat field *************************************/
if(ff_name[0] != '0')
 {
#ifdef APODIZATION
/* Creates Hamming apodization file: */
  isize = nx * ny * sizeof(float);
  JLP_GVM(&apodi,&isize);
#ifdef HAMMING
  jlp_hamming(apodi,nx,ny,nx);
#else
  jlp_blackman(apodi,nx,ny,nx);
#endif

#ifdef DEBUG
#ifdef HAMMING
   strcpy(outfile,"hamming");
   strcpy(outcomments,"Hamming filter: 0.54 + 0.46 cos(pi t / T)");
#else
   strcpy(outfile,"blackman");
   strcpy(outcomments,"Blackman filter: 0.42 + 0.5 cos(pi t / T) + 0.08 cos(2pit/T)");
#endif
   JLP_WRITEIMAG(apodi,&nx,&ny,&nx,outfile,outcomments);
#endif

/* End of case (#ifdef APODIZATION) */
#endif

/****************************************************************/
/* Read flat field file: */
  JLP_VM_READIMAG(&pntr_ffield,&nx_ff,&ny_ff,ff_name,ff_comments);
  JLP_FROM_MADRID(&pntr_ffield,&in_ffield);

  if(nx_ff != ny_ff || nx_ff < nx)
  {printf(" Fatal error/Ffield has wrong size: it should be square and larger than output image size)!\n");
  fprintf(fp1," Fatal error/Ffield has wrong size: it should be square and larger than output image size)!\n");
   exit(-1);
  }
/* Reduction of the flat field to fit the image size: */
  ireduc_ff = nx_ff / nx;
  if(nx_ff != ireduc_ff * nx)
  {printf(" Fatal error/Ffield has wrong size: it should be a multiple of the output image size)!\n");
  fprintf(fp1," Fatal error/Ffield has wrong size: it should be a multiple of the output image size)!\n");
   exit(-1);
  }

/* Allocate memory space for flat field: */
 isize = nx * ny *sizeof(float);
 JLP_GVM(&ffield,&isize);

  create_ffield(ffield,in_ffield,apodi,nx_ff,ny_ff,nx,ny,ireduc_ff,0.6);

  printf(" Ffield successfully built \n");
  strcpy(outfile,"ffield_car2_debug");
  strcpy(outcomments,"inverse ffield reduced");
  JLP_WRITEIMAG(ffield,&nx,&ny,&nx,outfile,outcomments);

/* Free virtual memory space: */
  JLP_FVM(&in_ffield);
#ifdef APODIZATION
  JLP_FVM(&apodi);
#endif

}
/********************* Unity file for ffield: *************/
else
{
  fprintf(fp1,"\n OK, Take unity file for ffield \n");
  printf("\n OK, Take unity file for ffield \n");
 isize = nx * ny *sizeof(float);
 JLP_GVM(&ffield,&isize);
 for(i = 0; i < nx * ny; i++) ffield[i] = 1.;
}
#endif

/* Photon processing: */
  status = rdcar_phot(ix_1,iy_1,val_1,autoco,interco,
               long_int,ffield,nx,ny,nx,nx_auto,ny_auto,nx_auto,
               ixcent,iycent,in_name,npack,
               nphot_wanted,&nphot_found);
  printf(" Output from rdcar: %d photons processed\n", nphot_found);
  fprintf(fp1," Output from rdcar: %d photons processed\n", nphot_found);

  sprintf(outcomments,"%s %d-08-93 %.2fH %.1fs %dph %dph",
          in_name,(int)date,time,integ,npack,nphot_found);


/* Mean of the frames: */
   w1 = (float)nphot_found/(float)npack;
   for(i = 0; i < nx * ny; i++) {
        long_int[i] /= w1;
        }
/* "Symmetry" of autocorrelation and intercorrelation */
   autoco_sym1(autoco,nx_auto,ny_auto,nx_auto);
   autoco_sym1(interco,nx_auto,ny_auto,nx_auto);

/* Normalization: */
   normalize(autoco,nx_auto,ny_auto,nx_auto);
   normalize(interco,nx_auto,ny_auto,nx_auto);

/****************************************************************/
/* Now output the results : */

/* Autocorrelation : */
   sprintf(outfile,"autoc%s",file_ext);
   fprintf(fp1," Output of %s",outfile);
   JLP_WRITEIMAG(autoco,&nx_auto,&ny_auto,&nx_auto,outfile,outcomments);

/* Intercorrelation : */
   sprintf(outfile,"interc%s",file_ext);
   fprintf(fp1," Output of %s",outfile);
   JLP_WRITEIMAG(interco,&nx_auto,&ny_auto,&nx_auto,outfile,outcomments);

/* Autocorrelation corrected by subtracting the intercorrelation: */
   for(i = 0; i < nx_auto * ny_auto; i++) {
      autoco[i] -= interco[i];
      }
   sprintf(outfile,"autocc%s",file_ext);
   fprintf(fp1," Output of %s",outfile);
   JLP_WRITEIMAG(autoco,&nx_auto,&ny_auto,&nx_auto,outfile,outcomments);

/* Long integration : */
   sprintf(outfile,"long%s",file_ext);
   fprintf(fp1," %s",outfile);
   JLP_WRITEIMAG(long_int,&nx,&ny,&nx,outfile,outcomments);

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
register int i, j, k;

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
Photon processing
*
****************************************************************/
/* Subroutine rdcar_phot to read CAR format */
int rdcar_phot(ix_1,iy_1,val_1,autoco,interco,long_int,
               ffield,nx,ny,idim,nx_auto,ny_auto,idim_auto,
               ixcent,iycent,in_name,npack,
               nphot_wanted,nphot_found)
long ix_1[], iy_1[];
float val_1[], ffield[];
float autoco[], interco[], long_int[];
long nx_auto, ny_auto, idim_auto;
long nx, ny, idim;
int ixcent, iycent;
int npack, *nphot_found, nphot_wanted;
char in_name[];
{
FILE *fd;
int nbytes_to_read, nbytes, nvalues, nblock, nblock_max; 
int ix, iy, itime, ireduc;
register int ii, i, j, k;

union{
long lg[BSIZE];
char ch[4*BSIZE];
} buff;

void swap_lint(), inv_lint();

/* Computes reduction factor: */
ireduc = CAR_WIDTH / nx;

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
nblock_max = 1 + nphot_wanted / BSIZE ;
if(nblock_max > NBLOCKMAX) nblock_max = NBLOCKMAX;
for(nblock = 0; nblock < nblock_max; nblock++)
{
/* Read next values */
  nvalues = fread(buff.lg,sizeof(long),BSIZE,fd);
  if(nvalues != BSIZE)
          printf(" rdcar/Warning, only %d values read in block #%d\n",
                   nvalues,nblock);
  if(nvalues <= 0 )
  {
  printf("rdcar/end of file: >%s< \n",in_name);
  nblock = NBLOCKMAX;
  }

/* Displays current value of nblock: */
#ifdef DEBUG
  if((nblock % 50) == 1)
     printf(" Block #%d  of %d bytes\n",nblock,4*BSIZE);
#else
  if((nblock % 100) == 1)
     printf(" Block #%d  of %d bytes\n",nblock,4*BSIZE);
#endif

/* Reads photon coordinates and fills buffer array: */
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

/* Division by "ireduc" to reduce size of output images
*/
    ix_1[ii] = ix / ireduc; iy_1[ii] = iy / ireduc;

/* Building the long integration: */
#ifdef FFIELD_CORRECTION
/* Store photon at the coordinates location: */ 
         j = ix_1[ii] + iy_1[ii] * idim;
         long_int[j] += ffield[j];
         val_1[ii] = ffield[j];
#else
         long_int[ix_1[ii] + iy_1[ii] * idim]++;
#endif
         (*nphot_found)++; 

/* End of loop on ii (nvalues) */ 
  }

/* Now process buffer arrays: */
   autocor_photon(ffield,ix_1,iy_1,val_1,autoco,interco,nx,ny,idim,
                  nx_auto,ny_auto,idim_auto,npack,nvalues);
/* End of current block: read next values */
}

/* Closes the input file */
  fclose(fd);
  return(0);
}
/****************************************************************/
int create_ffield(ffield,in_ffield,apodi,nx_ff,ny_ff,nx,ny,ireduc_ff,mini_value)
float ffield[], in_ffield[], apodi[], mini_value;
int nx, ny, nx_ff, ny_ff, ireduc_ff;
{
register int i, j, ii, jj;
int nvalues;
double sum;
float ff_mean;

/* Transfer now to smaller array (and mean for each pixel): */
for(j = 0; j < ny_ff; j++)
  {
  jj = j / ireduc_ff;
  for(i = 0; i < nx_ff; i++)
    {
    ii = i / ireduc_ff;
    ffield[ii + jj * nx] = in_ffield[i + j * nx_ff];
    }
  }

/* First normalization: (take into account only non zero values) */
sum = 0.; nvalues=0;
for(j = 0; j < nx * ny; j++) 
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

for(j = 0; j < nx * ny; j++) ffield[j] /= ff_mean;
 
/* Inversion of flat field: */
/* For the CAR, when the relative value is under mini_value=0.6, pb even with
   apodization. As this means that we are anyway outside the 
   sensitive zone, we set the maximum to that value (or 1/value here): */ 
for(j = 0; j < nx * ny; j++) 
     if(ffield[j] < mini_value) 
         ffield[j] = 1./mini_value;
     else 
         ffield[j] = 1./ffield[j];

/* Apodization: */
#ifdef APODIZATION
for(j = 0; j < nx * ny; j++) ffield[j] *= apodi[j];
#endif

/* Second normalization: (take into account only values close to mean)*/
sum = 0.; nvalues=0;
for(j = 0; j < nx * ny; j++) 
    {if(ffield[j] < 2. && ffield[j] > 0.5)
       {
       sum += ffield[j];
       nvalues++;
       }
    }
ff_mean = sum / (float)(nvalues);
for(j = 0; j < nx * ny; j++) ffield[j] /= ff_mean;
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
double PI = 3.14159;
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
double PI = 3.14159;
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
double PI = 3.14159;
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
double PI = 3.14159, wx1, wy1;
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
/*************************************************************
* autocor_photon
*
*************************************************************/
int autocor_photon(ffield,ix_1,iy_1,val_1,autoco,interco,nx,ny,idim,
                   nx_auto,ny_auto,idim_auto,npack,nvalues)
long ix_1[], iy_1[];
float ffield[], val_1[];
float autoco[], interco[];
long int nx, ny, idim, nx_auto, ny_auto, idim_auto;
int npack, nvalues;
{
float val;
int xmax, ymax; 
register long ixc, iyc, iwx, iwy;
register int i, j, ii;

/* Center coordinates for autocorrelation: */
ixc = nx_auto/2;
/* Please note that */
iyc = 0;

xmax = nx_auto/2; 
ymax = ny_auto; 

for(ii = 0; ii < nvalues-npack; ii++) 
  {

#ifdef FFIELD_CORRECTION
   val = val_1[ii];
#endif

/* Building the autocorrelation: */
   for(i = ii+1; i < ii + npack - 1; i++)
      {
       iwx = ix_1[i] - ix_1[ii];
       iwy = iy_1[i] - iy_1[ii];
/* Take only upper half autocorrelation (symmetry relative to center) */
       if(iwy < 0) {iwy = -iwy; iwx = -iwx;}
       if(iwx > -xmax && iwx < xmax && iwy < ymax)
          { 
          j = iwx + ixc + iwy * idim_auto;
#ifdef FFIELD_CORRECTION
          autoco[j] += val * val_1[i];
#else
          autoco[j]++;
#endif
          }
      }

/* End of loop on ii (nvalues) */ 
  }

/* Building the intercorrelation: (correlation with [4*npack,5*npack]) */
for(ii = 0; ii < nvalues - 5*npack; ii++) 
  {
#ifdef FFIELD_CORRECTION
   val = val_1[ii];
#endif
   for(i = ii + 4 * npack; i < ii + 5 * npack - 2; i++)
      {
       iwx = ix_1[i] - ix_1[ii];
       iwy = iy_1[i] - iy_1[ii];
/* Take only upper half intercorrelation (symmetry relative to center) */
       if(iwy < 0) {iwy = -iwy; iwx = -iwx;}
       if(iwx > -xmax && iwx < xmax && iwy < ymax)
           {
           j = iwx + ixc + iwy * idim_auto;
#ifdef FFIELD_CORRECTION
           interco[j] += val_1[i] * val;
#else
           interco[j]++;
#endif
           }
      }
  }

 return(0);
}
/***************************************************
* Symmetry of autocorrelation 
* Special version for half frame
*
****************************************************/
int autoco_sym1(autoco,nx_auto,ny_auto,idim_auto)
int nx_auto, ny_auto, idim_auto;
float autoco[];
{
register int i, j;

/* X axis: */
  for(i = 1; i <= nx_auto/2; i++)
    {
    autoco[i] += autoco[nx_auto - i];
    autoco[nx_auto - i] =  autoco[i];
    }

return(0);
}
/*
    autoco[i + j * idim_auto] += 
           autoco[(nx_auto - i) + (ny_auto - j) * idim_auto];
    autoco[(nx_auto - i) + (ny_auto - j) * idim_auto] =  
           autoco[i + j * idim_auto];
*/
/*******************************************************************/
int normalize(array,nx,ny,idim)
float array[];
int nx, ny, idim;
{
register int i, j, jj;
double sum;

/* First computing the total sum: */
sum = 0.;
 for(j = 0; j < ny; j++) 
  {
   jj = j * idim;
   for(i = 0; i < nx; i++) 
    {
     sum += array[i + jj];
    }
  }

/* Then dividing by this sum: */
 if(sum == 0) 
    {
    printf("normalize/error: total sum is null!\n");
    return(-1);
    }

 for(j = 0; j < ny; j++) 
  {
   jj = j * idim;
   for(i = 0; i < nx; i++) 
    {
     array[i + jj] /= sum;
    }
  }

 return(0);
}
/*********************************************
* Displays a long integer with binary code 
**********************************************/
int bin_long_disp(a)
long a;
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
}
/*************************************************************
* Tests to see if the current machine is working properly 
**************************************************************/
int test_machine()
{
long i1, i2;

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
  
}
