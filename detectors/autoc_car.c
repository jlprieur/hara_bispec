/**************************************************************** 
 autoc_car.c
 (formely called: dec_car3.c)
 To decode data from the CAR photon counting camera 
 and compute long integration and autocorrelation
 Fast version of decode_car2 which computes FFT's and bispectrum. 

 Flat field correction (if defined)
 cpu time optimized compared to dec_car4.c

 JLP
 Version 01-03-01
*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <jlp_ftoc.h>

/*
#define DEBUG
*/
#define swap_like_dec 
/*
#define FFIELD_CORRECTION 
*/

/* Block size (unit is sizeof(long)) */
#define BSIZE 8192 
#define NBLOCKMAX 100000
/* Maximum size for the images */
#define IDIM 512 
/* CAR detector has 1024x1024 pixels : */
#define CAR_WIDTH 1024 

/* Contained here: */
int rdcar_phot(INT4 *ixy_1, float *val_1, float *autoco, float *interco, 
               float *long_int, float *ffield, INT4 nx, INT4 ny, INT4 idim, 
               INT4 ireduc, INT4 ixcent, INT4 iycent, char *in_name, INT4 npack,
               INT4 nphot_wanted, INT4 *nphot_found);

main(argc, argv)
int argc;
char **argv;
{
/* Images : */
float *autoco, *interco, *long_int;
float *val_1;
INT4 *ixy_1, idd, mm;
float *ffield;
float w1, w2, integ;
double aa, time, date_in_years;
INT4 npack, nphot_wanted, nphot_header, nphot_found;
INT4  status, status1, status2, ireduc, isize, ixcent, iycent;
INT4 nx, ny, nx_auto, ny_auto;
register int i;
char *pc, file_ext[41], logfile[41];
char in_name[61], ff_name[61], outfile[61], outcomments[81];
FILE *fp1;

printf(" Program autoc_car  (%d x %d maxi) \n",IDIM,IDIM);
printf(" JLP Version 01-03-01 \n");

/*
* Tests to see if the current machine is working properly 
* (i.e., compatible with fast autocorrelation used by autoc_car.c)
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
  printf(" runs autoc_car in_photon_file output_file_exten");
  printf(" reduc_fact,total_nphot,phot_per_frame flat_field\n");
  printf(" Example: runs autoc_car test tt 16,1000,100 ff_oc45 \n\n"); 
  printf(" or without flat field: runs autoc_car test tt 16,1000,100 0 \n\n"); 
  exit(-1);
  }

/* Interactive input of in_name and file_ext: */
if (argc == 1)
 {
   printf(" Input file := "); scanf("%s",in_name);
/************* File extension for output files: */
   printf(" Output file extension := "); scanf("%s",file_ext);
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
sprintf(logfile,"autoc_car%s.log",file_ext);
if((fp1 = fopen(logfile,"w")) == NULL)
   {
   printf("autoc_car/Fatal error opening logfile >%s< \n",logfile);
   exit(-1);
   }

fprintf(fp1," Program autoc_car  (%d x %d maxi) \n",IDIM,IDIM);
fprintf(fp1," Sliding window method\n");
fprintf(fp1," JLP Version 14-05-96 \n");

/***************************************************************/
/* Read header first (for interactive input of the parameters...) */
  status = rdcar_header(in_name,&date_in_years,&time,&idd,&mm,&aa,
                        &integ,&nphot_header);
  if(status) exit(-1);
  printf(" Input file:%s date=%.4f time=%.2fh expo=%.1fsec nphot=%d\n",
          in_name,date_in_years,time,integ,nphot_header);
  fprintf(fp1," Input file:%s date=%.4f (%02d-%02d-%4d) \
time=%.2fH itime=%.2f nphot=%d \n",
          in_name,date_in_years,idd,mm,(int)aa,time,integ,nphot_header);
  w1 = nphot_header/(50. * integ);
  w2 = nphot_header/integ;
  printf(" Average flux in 20 msec was: %.3f photons (i.e. %.3f photons/s)\n",
          w1,w2);
  fprintf(fp1," Average flux in 20 msec was: %.3f photons (i.e. %.3f photons/s)\n",
          w1,w2);

/***************************************************************/
/* In case of interactive input of parameters, it is better to ask
 the following when the header has been read */ 
if (argc == 1)
  {
  printf(" Select reduction factor of the output frames, (4, 8, 16 or 32),");
  printf(" ixcent, iycent, nx,(power of two), ny (power of two),\
 total number of photons to be processed,");
  printf(" number of photons per frame \n");
  printf("     (Example: 16,512,512,64,64,2000,250) \n");
  scanf("%d,%d,%d,%d,%d,%d,%d",
        &ireduc,&ixcent,&iycent,&nx,&ny,&nphot_wanted,&npack);
  }
else
  {
  sscanf(argv[3],"%d,%d,%d,%d,%d,%d,%d",
               &ireduc,&ixcent,&iycent,&nx,&ny,&nphot_wanted,&npack);
  }

/* Logfile: */
  fprintf(fp1," Reduction factor of the output frames: %d \n",ireduc);
  fprintf(fp1," IXcenter=%d, IYcenter=%d\n",
           ixcent,iycent);
  fprintf(fp1," Total number of photons to be processed: %d\n",nphot_wanted);
  fprintf(fp1," Number of photons per frame: %d \n",npack);

/* Set size of output images */
/* Take powers of two */
  for(i = IDIM; i > 0; i /= 2)
   { 
    if(nx >= i) {nx = i; break;} 
   }
  for(i = IDIM; i > 0; i /= 2)
   { 
    if(ny >= i) {ny = i; break;} 
   }
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

isize = nx * ny * sizeof(float);
  long_int = (float *) malloc(isize);

/* To avoid tests in most frequently used loop, take double size: */
nx_auto = 2 * nx;
ny_auto = 2 * ny;
isize = nx_auto * ny_auto * sizeof(float);
  autoco = (float *) malloc(isize);
  interco = (float *) malloc(isize);

/* Buffer to store coordinates of a whole block */
isize = BSIZE * sizeof(INT4);
  ixy_1 = (INT4 *) malloc(isize);

#ifdef FFIELD_CORRECTION
isize = BSIZE * sizeof(float);
  val_1 = (float *) malloc(isize);
    for(i = 0; i < BSIZE; i++) val_1[i] = 0.;
#endif

/* Erasing the arrays : */
    for(i = 0; i < BSIZE; i++) {
      ixy_1[i] = 0;
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
  compute_ffield(ffield,ff_name,&nx_ff,&ny_ff,ixcent,iycent,nx,ny,
                 file_ext,fp1,0.6);
#endif

/* Photon processing: */
  status = rdcar_phot(ixy_1,val_1,autoco,interco,
               long_int,ffield,nx,ny,nx,ireduc,ixcent,iycent,in_name,npack,
               nphot_wanted,&nphot_found);
  printf(" Output from rdcar: %d photons processed\n", nphot_found);
  fprintf(fp1," Output from rdcar: %d photons processed\n", nphot_found);

  sprintf(outcomments,"%s %02d-%02d-%04d %.2fH %.1fs %dph %dph",
          in_name,idd,mm,(int)aa,time,integ,npack,nphot_found);

/* Mean of the frames: */
     w1 = (float)nphot_found/(float)npack;
/* Test just in case of long integrations, when npack is not reached: */
   if(w1 > 1.0)
     {
     for(i = 0; i < nx * ny; i++) {
        long_int[i] /= w1;
        }
     }

/* Symmetry of autocorrelation and intercorrelation */
   autoco_sym(autoco,nx_auto,ny_auto,nx_auto);
   autoco_sym(interco,nx_auto,ny_auto,nx_auto);

/* Normalization: */
   status1=normalize_float(autoco,nx_auto,ny_auto,nx_auto);
   status2=normalize_float(interco,nx_auto,ny_auto,nx_auto);

/****************************************************************/
/* Now output the results : */

   if(!status1 && !status2)
    { 
#ifdef FULL_OUTPUT
/* Autocorrelation : */
   sprintf(outfile,"autoc%s",file_ext);
   fprintf(fp1," Output of %s",outfile);
   JLP_WRITEIMAG(autoco,&nx_auto,&ny_auto,&nx_auto,outfile,outcomments);

/* Intercorrelation : */
   sprintf(outfile,"interc%s",file_ext);
   fprintf(fp1," Output of %s",outfile);
   JLP_WRITEIMAG(interco,&nx_auto,&ny_auto,&nx_auto,outfile,outcomments);
#endif

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
    printf(" Sorry, since pb in 'normalize', do not output autocor and interc\n");
   }

/* Long integration : */
   sprintf(outfile,"long%s",file_ext);
   fprintf(fp1," %s",outfile);
   JLP_WRITEIMAG(long_int,&nx,&ny,&nx,outfile,outcomments);

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
Photon processing
*
****************************************************************/
/* Subroutine rdcar_phot to read CAR format */
int rdcar_phot(INT4 *ixy_1, float *val_1, float *autoco, float *interco, 
               float *long_int, float *ffield, INT4 nx, INT4 ny, INT4 idim, 
               INT4 ireduc, INT4 ixcent, INT4 iycent, char *in_name, INT4 npack,
               INT4 nphot_wanted, INT4 *nphot_found)
{
FILE *fd;
INT4 nbytes_to_read, nbytes, nvalues, nblock, nblock_max; 
INT4 good_values, ix, iy, itime, idim2, ix_start, iy_start, ix_end, iy_end;
register int ii;

union{
unsigned long ulg[BSIZE];
char ch[4*BSIZE];
} buff;

void inv_ulong_int(unsigned long *i);

/* X size of autocorrelations: */ 
idim2 = idim * 2;
ix_start = ixcent - (nx * ireduc)/2;
if(ix_start < 0) ix_start = 0;
ix_end = ix_start + nx * ireduc;
iy_start = iycent - (ny * ireduc)/2;
if(iy_start < 0) iy_start = 0;
iy_end = iy_start + ny * ireduc;
 printf(" ix_start,end: %d,%d iy_start, end: %d,%d \n",
         ix_start, ix_end, iy_start, iy_end);

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
  nvalues = fread(buff.ulg,sizeof(unsigned long),BSIZE,fd);
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
good_values = 0;
for(ii = 0; ii < nvalues; ii++) 
  {

#ifdef swap_like_dec
    inv_ulong_int(&buff.ulg[ii]);
#endif
/* Use a mask: (since sometimes does not fill with zeroes) */
    itime = (buff.ulg[ii] >> 20) & 4095;
    iy = (buff.ulg[ii] >> 10) & 1023; 
    ix = buff.ulg[ii] & 1023; 

/* Select coordinates within the working frame: */
    if(ix >= ix_start && ix < ix_end && iy >= iy_start && iy < iy_end) 
      {
        ix -= ix_start; iy -= iy_start;
/* Division by "ireduc" to reduce size of output images
*/
         ix /= ireduc; iy /= ireduc;
/* Storing ix and iy in a buffer such that ix and iy can be separated
   for further additions and subtractions (since idim is a power of two): */
         ixy_1[ii] = ix + iy * idim2;
/* Building the long integration: */
#ifdef FFIELD_CORRECTION
/* Store photon at the coordinates location: */ 
         j = ix + iy * idim;
         long_int[j] += ffield[j];
         val_1[ii] = ffield[j];
#else
         long_int[ix + iy * idim] += 1.;
#endif
         good_values++;
     }

/* End of loop on ii (nvalues) */ 
  }
   *nphot_found += good_values; 

/* Now process buffer arrays: */
   car_autocor_photon(ffield,ixy_1,val_1,autoco,interco,idim,ny,npack,good_values);
/* End of current block: read next values */
}

/* Closes the input file */
  fclose(fd);
  return(0);
}
/*************************************************************
* car_autocor_photon
* Sliding window method
*
*************************************************************/
int car_autocor_photon(ffield,ixy_1,val_1,autoco,interco,idim,ny,npack,nvalues)
INT4 ixy_1[];
float ffield[], val_1[];
float autoco[], interco[];
INT4 npack, idim, ny, nvalues;
{
register INT4 ixyc, iw1, idim2;
register int i, ii;

/* Autocorrelation has double size (2*idim) */
idim2 = 2 * idim;
ixyc = idim + ny * idim2;
#ifdef DEBUG
  printf("DEBUG/ ixyc = %d (idim2 = %d) \n",idim2,ixyc);
  bin_int4_disp(ixyc);
#endif

for(ii = 0; ii < nvalues-npack; ii++) 
  {
/* The separation (powers of two) allows coordinate difference 
  in one operation:*/
   iw1 = ixyc - ixy_1[ii];

#ifdef FFIELD_CORRECTION
   val = val_1[ii];
#endif

/* Building the autocorrelation: 
* Value of autoc in (i,j shifted by ixc, iyc)
* equals Sum_{ixy_1[i]-iyx_1[ii] = (i,j)} of val((ixy_1[i] - ixy_1[ii]) + ixyc)*/
   for(i = ii+1; i < ii + npack - 1; i++)
      {
#ifdef FFIELD_CORRECTION
       autoco[ixy_1[i] + iw1] += val * val_1[i];
#else
       autoco[ixy_1[i] + iw1]++;
#endif
      }

/*
       j = (ix_1[i] - ix_1[ii]) + ixc + ((iy_1[i] - iy_1[ii]) + iyc) * idim2;
       j = ix_1[i] + iw + iy_1[i] * idim2;
       autoco[j] += val * val_1[i];
       autoco[ix_1[i] + iw + iy_1[i] * idim2] += val * val_1[i];
*/

/* End of loop on ii (nvalues) */ 
  }

/* Building the intercorrelation: (correlation with [4*npack,5*npack]) */
for(ii = 0; ii < nvalues - 5*npack; ii++) 
  {
   iw1 = ixyc - ixy_1[ii];
#ifdef FFIELD_CORRECTION
   val = val_1[ii];
#endif
   for(i = ii + 4 * npack; i < ii + 5 * npack - 2; i++)
      {
#ifdef FFIELD_CORRECTION
       interco[ixy_1[i] + iw1] += val_1[i] * val;
#else
       interco[ixy_1[i] + iw1]++;
#endif
      }
  }

 return(0);
}
/*********************************************
* Displays a long integer with binary code 
**********************************************/
int bin_int4_disp(a)
INT4 a;
{
INT4 i, j;

/* Long integer is four byte long: */
 for (j = 0 ; j < 4 ; j++)
  {
    printf ("  ");     
    for (i = 0 ; i < 8 ; i++)
     {
/* Mask with 8 * 16**7 = 10000000 00000000 00000000 00000000 */
      printf ("%1d", (a & 0x80000000) >> 31);
      a = a << 1;
     }
  }
 printf ("\n");     
return(0);
}
