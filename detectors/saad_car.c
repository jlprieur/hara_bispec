/**************************************************************** 
 saad_car.c
 To decode data from the CAR photon counting camera 
 and compute an image with the Shift And Add method 

 Flat field correction (before compression)

 JLP
 Version 06-03-97
*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <jlp_ftoc.h>

#define DEBUG 1

#define NBLOCKMAX 100000

/* Maximum size for the images */
#define IDIM 256
/* CAR detector has 1024x1024 pixels, but useful size is 864x864 : */
#define CAR_WIDTH 864 

main(argc, argv)
int argc;
char **argv;
{
/* Images : */
float *image1, *long_int;
float *float_array, *pntr_image;
float *in_ffield, *ffield;
float xframes, w1, w2, w3, date, integ, xphot;
double aa, time, date_in_years;
long int mm, idd;
int nframes_expected, photon_correction;
int npack, nphot_wanted, nphot_header, nphot_found, nframes;
int status, ireduc, isize, ixcent, iycent, ixmax, iymax;
long int pntr_ima, nx, ny, nx_ff, ny_ff;
long int descr_length;
register long int i;
char *pc, file_ext[41], buffer[81], logfile[41];
char descr_name[61], descr_value[81];
char in_name[61], ff_name[61], ff_comments[81], outfile[61], outcomments[81];
FILE *fp1;

printf(" Program saad_car  (%d x %d maxi) \n",IDIM,IDIM);
printf(" JLP Version 06-03-97 \n");
printf(" You must know the position of the maximum (ixmax,iymax) in \
the coordinates of the final image\n"); 


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
  printf(" runs saad_car in_photon_file output_file_exten");
  printf(" reduc_fact,ixcent,iycent,ixmax,iymax,total_nphot,frame_nphot flat_field\n");
  printf(" Example: \n"); 
  printf("  runs saad_car test tt 16,488,512,234,228,900,100 ff_oc45 \n\n"); 
  printf(" or without flat field: \n");
  printf("  runs saad_car test tt 16,488,512,234,228,1000,100 0 \n\n"); 
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
sprintf(logfile,"saad_car%s.log",file_ext);
if((fp1 = fopen(logfile,"w")) == NULL)
   {
   printf("saad_car/Fatal error opening logfile >%s< \n",logfile);
   exit(-1);
   }

fprintf(fp1," Program saad_car  (%d x %d maxi) \n",IDIM,IDIM);
fprintf(fp1," JLP Version 06-03-97 \n");

/***************************************************************/
/* Read header first (for interactive input of the parameters...) */
  status = rdcar_header(in_name,&date_in_years,&time,&idd,&mm,&aa,
                        &integ,&nphot_header);
  printf(" Input file:%s date=%.4f time=%.2fH itime=%.2f nphot=%d \n",
          in_name,date_in_years,time,integ,nphot_header);
  fprintf(fp1," Input file:%s date=%.4f (%02d-%02d-%4d) \
time=%.2fH itime=%.2f nphot=%d \n",
          in_name,date_in_years,idd,mm,(int)aa,time,integ,nphot_header);

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
  printf(" IXcenter, IYcenter, IXmax, IYmax, nx, ny, \n");
  printf(" total number of photons to be processed,");
  printf(" number of photons per frame \n");
  printf("     (Example: 16,488,512,234,288,32,32,2000,250) \n");
  gets(buffer);sscanf(buffer,"%d,%d,%d,%d,%d,%d,%d,%d,%d",
               &ireduc,&ixcent,&iycent,&ixmax,&iymax,
               &nx,&ny,&nphot_wanted,&npack);
  }
else
  {
  sscanf(argv[3],"%d,%d,%d,%d,%d,%d,%d,%d,%d",
               &ireduc,&ixcent,&iycent,&ixmax,&iymax,
               &nx,&ny,&nphot_wanted,&npack);
  }

/* Logfile: */
  fprintf(fp1," Reduction factor of the output frames: %d \n",ireduc);
  fprintf(fp1," IXcenter=%d, IYcenter=%d\n",ixcent,iycent);
  fprintf(fp1," IXmax=%d, IYmax=%d\n",ixmax,iymax);
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
    { printf(" saad_rcar/Warning, only %d photons in file!\n",nphot_header);
     fprintf(fp1," saad_rcar/Warning, only %d photons in file!\n",nphot_header);
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
  JLP_GVM(&image1,&isize);
  JLP_GVM(&long_int,&isize);

/* Erasing the arrays : */
    for(i = 0; i < nx * ny; i++) {
      image1[i] = 0.;
      long_int[i] = 0.;
      }
/******************************************************************/
/*************** Flat field *************************************/
/* Computing flat-field */
  nx_ff = 1024; ny_ff = 1024;
  compute_ffield(&ffield,ff_name,&nx_ff,&ny_ff,
                  ixcent,iycent,(int)(nx*ireduc),(int)(ny*ireduc),
                  file_ext,fp1,0.3);

nframes_expected = nphot_wanted / npack;

/* Photon processing: */
  status = rdcar_saad(image1,long_int,ffield,nx,ny,nx,
                      ixcent,iycent,ixmax,iymax,in_name,npack,nphot_wanted,
                      &nphot_found,&nframes,nframes_expected,nx_ff);
  printf(" Output from rdcar_saad: %d photons and %d frames processed\n",
          nphot_found, nframes);
  fprintf(fp1," Output from rdcar_saad: %d photons and %d frames processed\n",
          nphot_found, nframes);

/* Check that some frames have been processed */
   xframes=(float)nframes;
   xphot = nphot_found; 
  sprintf(outcomments,"%s %d-08-93 %.2fH %.1fs %.1fph %dfr %dph",
          in_name,(int)date,time,integ,xphot,nframes,nphot_found);

/****************************************************************/
   printf(" xphot = %f \n",xphot);
   fprintf(fp1," xphot = %f \n",xphot);

/* Mean of the frames: */
   for(i = 0; i < nx*ny; i++) {
        long_int[i] /= xframes;
        }

/****************************************************************/
/* Now output the results : */

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
/****************************************************************
*
Format CAR data 
Photon processing
*
****************************************************************/
/* Subroutine rdcar_saad to read CAR format */
int rdcar_saad(image1,long_int,ffield,nx,ny,idim,
               ixcent,iycent,ixmax,iymax,in_name,npack,
               nphot_wanted,nphot_found,nframes,nframes_expected,idim_ff)
float image1[], long_int[], ffield[];
long nx, ny, idim;
int ixcent, iycent, ixmax, iymax, idim_ff;
int npack, *nphot_found, nphot_wanted, *nframes, nframes_expected;
char in_name[];
{
FILE *fd;
int nbytes_to_read, nbytes, nvalues, nblock, nblock_max; 
int ix, iy, ixstart, iystart, ixend, iyend, itime, iph1, iframe, ireduc;
register int ii, i, j, j_ff, k;

union{
long lg[256];
char ch[1024];
} buff;

/* Computes reduction factor: */
ireduc = CAR_WIDTH / nx;
ixstart = (ixcent / ireduc) - nx/2; 
iystart = (iycent / ireduc) - ny/2; 
ixend = ixstart + nx;
iyend = iystart + ny;
#ifdef DEBUG
 printf("rdcar_saad/ ixstart = %d, iystart=%d ",ixstart,iystart);
 printf("ixend = %d, iyend=%d\n",ixend,iyend);
#endif
if(ixstart < 0 || iystart < 0 ) 
 {printf("rdcar_saad/Fatal error, ixstart = %d, iystart=%d ",ixstart,iystart);
  exit(-1);
 }

/* Opens the input file */
if((fd = fopen(in_name,"r")) == NULL)
  {
  printf("rdcar_saad/error opening input file: >%s< \n",in_name);
  return(-1);
  }

/***************************************************************/
/* Skip header (32 bytes containing "FORMYCAR"): */
  nbytes_to_read = 32;
  nbytes = fread(buff.ch,sizeof(char),nbytes_to_read,fd);
  if(nbytes != nbytes_to_read)
    {
     printf("rdcar_saad/error skipping header: \n");
     printf("       Only %d bytes read \n", nbytes);
     return(-2);
    }

/***************************************************************/
/* Read the data: */
*nphot_found = 0; iframe = 0; iph1 = 0; 
nblock_max = 1 + nphot_wanted / 256;
if(nblock_max > NBLOCKMAX) nblock_max = NBLOCKMAX;
for(nblock = 0; nblock < nblock_max; nblock++)
{
/* Read next values */
  nvalues = fread(buff.lg,sizeof(long),256,fd);
  if(nvalues != 256)
          printf(" rdcar_saad/Warning, only %d values read in block #%d\n",
                   nvalues,nblock);
  if(nvalues <= 0 )
  {
  printf("rdcar_saad/end of file: >%s< \n",in_name);
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
    j_ff = ix + iy * idim_ff;
#ifdef DEBUG
    if(j_ff < 0 || j_ff > idim_ff*idim_ff)
      printf(" DEBUG/error with ffield ..., ix=%d iy=%d \n",ix,iy);
#endif

/* Division by "ireduc" to reduce size of output images
*/
    ix /= ireduc; iy /= ireduc;

/* Store photon at the coordinates location: */ 
   if(ix >= ixstart && ix < ixend && iy >= iystart && iy < iyend)
           {
            j=(ix - ixstart) + (iy - iystart) * idim;
/* Building the image: */
            image1[j] += ffield[j_ff];
/*
            image1[j] ++; 
*/
            (*nphot_found)++; iph1++; 
           }
/*****************************************************************/
/* Test to see if npack has been reached
   (pack full, i.e. image1 with npack photons)
 */ 
    if(iph1 == npack)
    {
#ifndef DEBUG
     if((iframe % 200) == 0) 
#endif
          printf(" rdcar_saad/Processing frame #%d/%d \n",
                   iframe,nframes_expected);
     process_frame(image1, long_int, nx, ny, idim, iframe, ixmax, iymax);
    iph1 = 0; iframe++;
/* Reset image1 to zero: */
     for (j = 0; j < nx * ny; j++) image1[j] = 0.;
    }

/* End of loop on ii (nvalues) */ 
  }
/* End of current block: read next values */
}

/* In case no frame has been processed,
   we want to recover the long integration at least... */
if(!iframe)
  {
  printf("rdcar_saad/Before exiting from the program, save long integration\n");
  process_frame(image1, long_int, nx, ny, idim, iframe, ixmax, iymax);
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
int process_frame(image, long_int, nx, ny, idim, iframe, ixmax, iymax)
float image[], long_int[];
long nx, ny, idim, iframe, ixmax, iymax;
{
static float offset_max;
float max0, offset;
int i1, j1, ioffset, joffset, iloc, ixmax1, iymax1;
register long int i, j;

/* Looking for the maximum : */
  max0 = image[0];
  iloc = 0;
  for(i = 0; i < nx * ny; i++) 
   { if(image[i] > max0)
        { max0 = image[i]; iloc = i;}
   } 
  iymax1 = (int)(iloc / nx);
  ixmax1 = iloc  - iymax1 * nx;

  ioffset = ixmax - ixmax1;
  joffset = iymax - iymax1;
  offset = ioffset * ioffset + joffset * joffset;
  offset = sqrt((double)offset);

#ifdef DEBUG
  printf("#%d iloc=%d ixmax1=%d iymax1=%d ioffset=%d joffset=%d offset=%.1f\n",
           iframe, iloc, ixmax1, iymax1, ioffset, joffset, offset); 
#endif

/* Long integration : */
  offset_max = ((float)nx)/10.;
  if(offset < offset_max)
  {
  for(j = 0; j < ny; j++)
    {
    j1 = j + joffset;
    for(i = 0; i < nx; i++) 
      {
      i1 = i + ioffset;
      long_int[i1 + j1 * nx] += image[i + j * nx];
      }
    }
  }

return(0);
}
