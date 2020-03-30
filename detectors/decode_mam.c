/* decode_mam.c
 To decode data from mama photon counting camera 
 and compute autocorrelation and bispectrum.

 Photon noise correction 

 JLP
 Version 17-10-91
 */

#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <jlp_ftoc.h>

/*
#define DEBUG
*/

#define MAXPHOTONS 120000 
#define IDIM 256
/* Note that DIM = IDIM * IDIM     */
#define DIM 65536 
/* Note that DIM2 = (2*IDIM + 1) * (2*IDIM + 1)     */
#define DIM2 263169

main()
{
/* Images IDIM*IDIM: */

float out_auto[DIM], image[DIM], long_int[DIM];
int xwork[MAXPHOTONS], ywork[MAXPHOTONS], autocor[DIM2];
float sum, w1, w2;
FILE *fd1;
int ivalues, nvalues, nph, ixy, idim=IDIM;
unsigned char pacframe[MAXPHOTONS*2];
int i1, ix, iy, ixc, iyc, block_size, imax, kmax;
int in_descr, npack, nphotons;
int nskip, nframes, nx, ny, iframe, iph, istat;
register int i, j, k;
unsigned length1;
char infile[41], nphfile[41], nphcomments[81], buffer[81];
char outfile[41], outcomments[81], cdescr[1025];

  printf(" Program decode_mam  Version 07-11-91 (%d x %d) \n",idim,idim);

/* Erasing the arrays : */
    for(i=0; i<DIM; i++) {
      long_int[i]=0.;
      out_auto[i]=0.;
      }
    for(i=0; i<DIM2; i++) autocor[i]=0;

/* Reading the file with the number of photons */
  JLP_INQUIFMT();
  printf(" Total number of photons to be processed?\n");
  gets(buffer);sscanf(buffer,"%d",&nphotons);
  printf(" Number of photons per pack ?\n");
  gets(buffer);sscanf(buffer,"%d",&npack);
  npack = (npack > MAXPHOTONS - 10) ? MAXPHOTONS - 10 : npack;

  /*
  printf(" Size of the image: nx,ny ?\n");
  gets(buffer);sscanf(buffer,"%d,%d",&nx,&ny);
  */
  nx = 256; ny = 256;

  printf(" Will process %d photons in packs of %d photons \n",nphotons,npack);
  printf(" Output size is nx = %d, ny=%d \n",nx,ny);
  if(nx <=0 || ny <=0) 
  {printf(" Fatal error for the input size\n");
   exit(-1);}

/* Opening the data file */
  printf(" Input data file ?\n");
  gets(infile);
  printf(" Reading %s \n",infile);
  fd1=fopen(infile,"r");

/* Reading the first four-byte of the data file (integration time)*/
  nvalues = fread(pacframe,1,4,fd1);
  i = pacframe[0];
  j = pacframe[1];
  i = (i << 8) | j;
  j = pacframe[2];
  i = (i << 8) | j;
  j = pacframe[3];
  i = (i << 8) | j;
  printf(" Integration time is %d \n",i);

  printf(" Current pointer is %d \n",ftell(fd1));
  printf(" Number of bytes to skip?\n");
  gets(buffer);sscanf(buffer,"%d",&nskip);
  if(nskip >= 1) 
    {
    i = nskip;
    while(i > 2*MAXPHOTONS) 
     {
      nvalues = fread(pacframe,1,2*MAXPHOTONS,fd1);
      i = i - 2*MAXPHOTONS; 
     }
      nvalues = fread(pacframe,1,i,fd1);
    }
  printf(" Current pointer is %d \n",ftell(fd1));

/* Reading the first 2*(npack-1) data values in the data file */
  nvalues = fread(pacframe,1,2*(npack-1),fd1);

#ifdef DEBUG
  printf(" pacframe[0] %d \n",pacframe[0]);
  printf(" Current pointer is %d \n",ftell(fd1));
#endif

  if(nvalues != 2*(npack-1)) 
   {printf(" Fatal error reading data file: nvalues=%d \n",nvalues);
    goto end1;}

/* Transfer to work array: */
   for(i=0; i<=npack-2; i++)
     {ywork[i] = pacframe[i*2+1];
     xwork[i] = pacframe[i*2];
     ixy = IDIM*ywork[i] + xwork[i];
#ifdef DEBUG
     printf(" pacframe[%d]= %d pacframe[%d]= %d\n",2*i,pacframe[2*i],2*i+1,
	      pacframe[2*i+1]);
     printf(" x, y %d %d \n",xwork[i],ywork[i]);
#endif
/* Long integration : */
#ifdef DEBUG
     if(ixy < 0 || ixy >= DIM)printf(" error: ixy = %d \n",ixy);
#endif
     long_int[ixy]++;}
  
/****************************************************************/
/* Start main loop: */
  iph = npack-1;
  while (iph < nphotons) {

/* Read by sets of (MAXPHOTONS-npack)*2 */
    block_size = (nphotons - iph > (MAXPHOTONS - npack)) ? 
		 (MAXPHOTONS - npack) : nphotons - iph;

    nvalues = fread(pacframe,1,2*block_size,fd1);
    if(nvalues == 0) 
     {printf(" End of data file: nphotons=%d \n",iph);
      break;}

    printf("iph=%d, block_size=%d nvalues=%d \n",iph,block_size,nvalues);

/* Main loop on all the values of packframe */
    kmax = nvalues/2;
    for(k=0; k < kmax; k++)
      {
      iph++;
/* Fills in last value of work arrays: */
      i1 = npack-1+k;
      xwork[i1] = pacframe[2*k];
      ywork[i1] = pacframe[2*k+1];
      ixy = IDIM*ywork[i1] + xwork[i1];

#ifdef DEBUG
    printf(" x, y %d %d \n",xwork[i1],ywork[i1]);
#endif

/* Long integration : */
      long_int[ixy]++;

/* Auto-correlation: */
      auto_cor(xwork,ywork,autocor,npack,k,nx,ny);

      }

/* Shift the array down by kmax pixels: */
      imax = npack+kmax-1;
      for(i = kmax; i < imax; i++) {
      i1=i-kmax;
      xwork[i1] = xwork[i]; 
      ywork[i1] = ywork[i]; 
      }

   }

   nphotons = iph;

/****************************************************************/
/* Now output the results : */

/* Autocorrelation : */
   symmetry_autocor(out_auto,autocor,nx,ny);
   strcpy(outfile,"autocor.bdf                               ");
   sprintf(outcomments," nph: %d; npack: %d ",nphotons,npack);
   JLP_WRITEIMAG(out_auto,&nx,&ny,&idim,outfile,outcomments);

/* Long integration : */

#ifdef DEBUG
   sum =0;
   for(i=0; i<DIM; i++) sum=sum+long_int[i];
   printf(" sum = %f \n",sum);
#endif

   strcpy(outfile,"long.bdf                                ");
   sprintf(outcomments," nph: %d; npack: %d ",nphotons,npack);
   JLP_WRITEIMAG(long_int,&nx,&ny,&idim,outfile,outcomments);

  printf(" Successful end \n");

/* End : */
  end1: fclose(fd1);
  JLP_END();
}
/*******************************************************
*  Auto-correlation: 
* [0:nx, 0:ny]
*******************************************************/
int auto_cor(xwork,ywork,autocor,npack,k,nx,ny)
int xwork[], ywork[], npack;
int autocor[];
int nx, ny, k;
{
register int i;
int idim2, imax, ix, iy, ixy;

      imax = npack+k;
      idim2 = 2*IDIM + 1; 
      for(i = k+1; i < imax; i++) {
/* Compute the coordinates of the contribution of this photon to
the autocorrelation: */
        ix = nx + xwork[i] - xwork[k]; 
        iy = ny + ywork[i] - ywork[k]; 
        ixy = idim2 * iy + ix; 
        autocor[ixy]++ ;
       }

return(0);
}
/*******************************************************
*  Output format of auto-correlation: 
*  Sum (x,y) and (-x,-y), and then symmetry (x,y) (-x,-y) 
*******************************************************/
int symmetry_autocor(out_auto,autocor,nx,ny)
float out_auto[];
int autocor[];
int nx, ny;
{
int idim2, idim22, ixc, iyc, ioff, joff;
register int i, j;
/* Symmetry of autocorrelation 
 Pixel (i,j) is linked with (idim2-i,idim2-j)
 Goes from [0,0] to [512,512]
 Central row at y=256, central column at x=256
*/
   idim2 = 2*IDIM + 1;
   idim22 = 2*IDIM;
   ixc = nx;
   iyc = ny;
   for(j=0; j<iyc; j++)
   {
    for(i=0; i<idim2; i++)
    {
    autocor[(idim22 - i) + (idim22 - j)*idim2] = autocor[i + j*idim2] 
	      + autocor[(idim22 - i) + (idim22 - j)*idim2];
    autocor[i + j*idim2] = autocor[(idim22 - i) + (idim22 - j)*idim2];
    }
   }
/* Correction of central row: */
    j=iyc;
    for(i=0; i<ixc; i++)
    {
    autocor[(idim22 - i) + (idim22 - j)*idim2] = autocor[i + j*idim2] 
	      + autocor[(idim22 - i) + (idim22 - j)*idim2];
    autocor[i + j*idim2] = autocor[(idim22 - i) + (idim22 - j)*idim2];
    }

/* Troncation now: */
/* As for autocor the central row is at y=256, central column at x=256
  and for out_auto there are at y=128, and x=128
*/
   ioff = nx/2;
   joff = ny/2;
   for(j=0; j<ny; j++)
   {
    for(i=0; i<nx; i++)
    {
     out_auto[i+j*IDIM]=autocor[(i+ioff)+(j+joff)*idim2];
    }
   }
return(0);
}
