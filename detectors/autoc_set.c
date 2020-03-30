/**************************************************************** 
* autoc_set.c
*
* Set of routines used by dec_car3_1D, dec_car4 and autoc_car (formely dec_car3)
* To compute long integration and autocorrelation
*
* Contains:
* int autocor_photon(ffield,ixy_1,max_nphot,val_1,autoco,interco,idim,npack,nvalues)
* int autoco_sym(autoco,nx,ny,idim)
* int normalize_float(array,nx,ny,idim)
* int rdcar_header(in_name,date_in_years,time,idd, mm,aa,integ,nphot)
*
* JLP
* Version 01-03-01
*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <jlp_ftoc.h>

/* Origin of date is 01-01-1904 at 0H TU
 which corresponds to Julian day 2416480.500 */
#define DATE_ORIGIN 2416480.50

/* When working with DEC or LINUX/PC computer, reverses long integers... */
#define swap_like_dec

/*************************************************************
* autocor_photon
* Sliding window method of width "npack" on a buffer 
* containing "nvalues" coordinates ixy_1[] (with values val_1[] if ffield corr.)
*
* And building the intercorrelation (correlation of current coordinates
* with window located "4*npacks" further ([4*npack,5*npack])
*
* (called by autoc_car.c)
*************************************************************/
int autocor_photon(ffield,ixy_1,max_nphot,val_1,autoco,interco,idim,npack,nvalues)
INT4 ixy_1[], max_nphot;
float ffield[], val_1[];
float autoco[], interco[];
INT4 npack, idim, nvalues;
{
float val;
register INT4 ixyc, iw1, idim2;
register int i, j, ii;

/* Autocorrelation has double size (2*idim) */
idim2 = 2 * idim;
ixyc = idim + idim * idim2;
#ifdef DEBUG
  printf(" ixyc = (idim2 = %d) \n",idim2,ixyc);
  bin_long_disp(ixyc);
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
/***************************************************************
* autocor_cp40_old
*
* Compute autocorrelation 
* and intercorrelation with frame shifted of 100 msec
* Old version with a shift every frame, but using all the frames
****************************************************************/
int autocor_cp40_old(autoco, interco, ixy_1, max_nphot, nphot, nx, ny)
double *autoco, *interco;
INT4 *ixy_1, max_nphot, nphot, nx, ny;
{
register INT4 ixyc, iw1, idim2;
register int i, j, k, ii;

/* Autocorrelation has double size (2*idim) */
idim2 = 2 * nx;
ixyc = nx + ny * idim2;
#ifdef DEBUG
  printf(" ixyc = (idim2 = %d)  ",idim2,ixyc);
  bin_long_disp(ixyc);
#endif

for(ii = 0; ii < nphot; ii++)
  {
/* The separation (powers of two) allows coordinate difference
  in one operation:*/
/* ix and iy has been stored in a buffer ixy_1 
   such that ix and iy can be separated
   for further additions and subtractions (since idim is a power of two): */
   iw1 = ixyc - ixy_1[ii];

/* Building the autocorrelation:
* Value of autoc in (i,j shifted by ixc, iyc)
* equals Sum_{ixy_1[i]-iyx_1[ii] = (i,j)} of val((ixy_1[i] - ixy_1[ii]) + ixyc)
*/
   for(i = 0; i < nphot; i++)
      {
       autoco[ixy_1[i] + iw1]++;
      }
/* Remove autocorrelation of photon by itself: */
    autoco[ixyc]--;

/* End of loop on ii (nphot values) */
  }

/* Building the intercorrelation: (correlation with [4*max_nphot,5*max_nphot]) 
 Note that the next loops will be effective only when at least 5 frames
 have been entered: (since array initialized with -2) */

for(ii = 0; ii < nphot; ii++)
  {
   iw1 = ixyc - ixy_1[ii];
   for(i = 4*max_nphot; i < 5*max_nphot; i++)
      {
      if(ixy_1[i] < 0) break;
       interco[ixy_1[i] + iw1]++;
      }
  }

/* Set the next value after the last one to -1: */
ixy_1[nphot] = -1;

/* Shift of "max_nphot" values in ixy_1 to compute next intercorrelation: */
for(i = 0; i < max_nphot; i++) 
   {
   for(k = 3; k >= 0; k--)
       ixy_1[i + (k+1)*max_nphot] = ixy_1[i + k*max_nphot];
   }

return(0);
}
/***************************************************************
* autocor_cp40
*
* Compute autocorrelation 
* and intercorrelation with frame shifted of 100 msec
* New version without shifting every frame, but loosing 5% of the frames
* for the intercorrelation
* iframe : (from 1 to 100)
****************************************************************/
int autocor_cp40(autoco, interco, ixy_1, max_nphot, nphot, iframe, nx, ny)
double *autoco, *interco;
INT4 *ixy_1, max_nphot, nphot, nx, ny, iframe;
{
INT4 istart;
register INT4 ixyc, iw1, idim2;
register int i, j, k, ii;

/* Autocorrelation has double size (2*idim) */
idim2 = 2 * nx;
ixyc = nx + ny * idim2;
#ifdef DEBUG
  printf(" ixyc = (idim2 = %d)  ",idim2,ixyc);
  bin_long_disp(ixyc);
#endif

istart = iframe * max_nphot;

for(ii = istart; ii < istart + nphot; ii++)
  {
/* The separation (powers of two) allows coordinate difference
  in one operation:*/
/* ix and iy has been stored in a buffer ixy_1 
   such that ix and iy can be separated
   for further additions and subtractions (since idim is a power of two): */
   iw1 = ixyc - ixy_1[ii];

/* Building the autocorrelation:
* Value of autoc in (i,j shifted by ixc, iyc)
* equals Sum_{ixy_1[i]-iyx_1[ii] = (i,j)} of val((ixy_1[i] - ixy_1[ii]) + ixyc)
*/
   for(i = istart; i < istart + nphot; i++)
      {
       autoco[ixy_1[i] + iw1]++;
      }
/* Remove autocorrelation of photon by itself: */
    autoco[ixyc]--;

/* End of loop on ii (nphot values) */
  }

/* Building the intercorrelation: (correlation with [4*max_nphot,5*max_nphot]) 
*/

if(iframe > 4)
{
for(ii = istart; ii < istart + nphot; ii++)
  {
   iw1 = ixyc - ixy_1[ii];
   for(i = istart - 5*max_nphot; i < istart - 4*max_nphot; i++)
      {
/* Negative value is the signal to say the frame is finished */
      if(ixy_1[i] < 0) break;
       interco[ixy_1[i] + iw1]++;
      }
  }
}

/* Set the next value after the last one to -1: */
ixy_1[istart + nphot] = -1;

return(0);
}

/***************************************************
* Symmetry of autocorrelation 
*
****************************************************/
int autoco_sym(float *autoco, INT4 nx, INT4 ny, INT4 idim)
{
INT4 nx2, ny2;
register int i, j;

nx2 = nx / 2;
ny2 = ny / 2;
/* Cannot be made for the lower left edges */
for(j = 0; j < ny; j++) autoco[j * idim] = 0.;
for(i = 0; i < nx; i++) autoco[i] = 0.;

/* Scan all the lines from j=1 to j=ny/2-1: */
for(j = 1; j < ny2; j++)
  for(i = 1; i < nx; i++)
    {
    autoco[i + j * idim] += autoco[(nx - i) + (ny - j) * idim];
    autoco[(nx - i) + (ny - j) * idim] =  autoco[i + j * idim];
    }

/* Process special case of X axis (Y=0 frequencies, i.e., central line): */
j = ny2;
  for(i = 1; i < nx2; i++)
    {
    autoco[i + ny2 * idim] += autoco[(nx - i) + ny2 * idim];
    autoco[(nx - i) + ny2 * idim] =  autoco[i + ny2 * idim];
    }

/* Also central value: */
    autoco[nx2 + ny2 * idim] *= 2.;

return(0);
}
/*******************************************************************
* To normalize a float array
*
*******************************************************************/
int normalize_float(float *array, INT4 nx, INT4 ny, INT4 idim)
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
    printf("normalize_float/error: total sum is null!\n");
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
/****************************************************************
*
* Subroutine rdcar_header to read header in CAR format
*
****************************************************************/
int rdcar_header(in_name,date_in_years,time,idd, mm,aa,integ,nphot)
char in_name[];
float *integ;
double *time, *date_in_years, *aa;
INT4 *idd, *mm, *nphot;
{
FILE *fd;
INT4 nbytes_to_read, nbytes, nvalues; 
char cbuff[9];
float date;
unsigned long s_date, integ_time, nphot1;
register int i, j, k;

/* Header of 32 bytes */
typedef struct {
/* March 2001: do not replace unsigned long by INT4!!!! 
*/
unsigned long integ_time;          /* duree en msec */
unsigned long date;                /* date of observation 
                             (in seconds starting from 01-01-1904 at 0H TU)
                             (i.e. Julian day 2416480.500) */
unsigned long nber_of_photons;     /* number of photons */
unsigned long keyword1;            /* "FORM" in ASCII */
unsigned long keyword2;            /* "YCAR" in ASCII */
unsigned long nber_of_images;      /* used only for the CP40, here always 0 */
short refNum;             /* always 0 */
unsigned long read_already;        /* not used */ 
short everything_read;    /* not used */
} CAR_HEADER;

union{
INT4 lg[256];
char ch[1024];
} buff;

CAR_HEADER *chead;
void inv_ulong_int(unsigned long *i), swap_lint(unsigned short *i);

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
#ifdef swap_like_dec
   inv_ulong_int(&s_date);
#endif
/* s_date was in seconds, I convert it to days: */
   date = (float)s_date / 86400.;
/* Time in hours: */
   *time = 24. * (date - (int)date);
/* JLP96: Look for actual date */
   date += DATE_ORIGIN;
   inv_julian(date_in_years,aa,mm,idd,*time,(double)date);
   printf(" Date is: %f or %02d-%02d-%04d \n",
          *date_in_years,(int)*idd,(int)*mm,(int)*aa);
/*************** integration time *******************************/
   integ_time = chead->integ_time;
#ifdef swap_like_dec
   inv_ulong_int(&integ_time);
#endif
/* integration time was in milliseconds, I convert it to seconds: */
   *integ = (float)integ_time / 1000.;

/*************** number of photons *******************************/
   nphot1 =chead->nber_of_photons;
#ifdef swap_like_dec
   inv_ulong_int(&nphot1);
#endif
   *nphot = nphot1;

/* Closes the input file */
  fclose(fd);
  return(0);
}
