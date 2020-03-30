/***********************************************************************
* decode_simu2.c
* To decode simulated data from "simu2.c" and compute modsq and bispectrum
* Photon noise correction 
*
* JLP
* Version 17-07-2008
*************************************************************************/
#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <jlp_ftoc.h>
#include <jlp_cover_mask.h>

/*
#define DEBUG
*/
#define MAXPHOTONS 24000 
/* For IR=30:
#define NGMAX 388400 
For IR=25:
*/
#define NGMAX 187566 
#define IDIM 128
#define DIM (IDIM * IDIM) 

int main(int argc, char *argv[])
{
/* Images IDIM*IDIM, irmax=25, ngmax=187566 : */

double image[DIM], modsq[DIM], snrm[DIM], mre[DIM], mim[DIM], long_int[DIM];
double yce1[4 * NGMAX], im[DIM];
float xframes, *nphotons, xphotons, w1, w2, w3;
int fd1, nvalues, nph, ix, iy;
INT4 nx1, ny1;
INT_PNTR pntr;
/*
short int pacframe[MAXPHOTONS*2+1];
*/
char pacframe[MAXPHOTONS*2+1];
int ir, max_nclosure, nbeta, ngamma, maxframes;
int nframes, nx, ny, kod, iframe, istat, nval;
char infile[41], nphfile[41], nphcomments[81];
char outfile[41], outcomments[81];
/*
int open(), read(), close();
*/
register int i;

nx = IDIM; ny = IDIM;
  printf(" Program decode_simu2  Version 25-05-2008 (%d x %d maxi) \n",
         IDIM, IDIM);

/* Input of parameters: */
if (argc != 5 )
  {
printf("\nUSAGE:\n");
printf("decode_simu2 data_file photons_file ir,max_nclosure Max_Number_of_frames\n ");
  exit(-1);
  }

/* Data file name: */
 strcpy(infile,argv[1]);

/* Photon file name: */
 strcpy(nphfile,argv[2]);

/* Radius of uv-coverage: */
 nval = sscanf(argv[3], "%d,%d", &ir, &max_nclosure);
 printf("nval=%d\n", nval);
 if(nval != 2) {
 fprintf(stderr,"Fatal error: ir=%d max_nclosure=%d\n", ir, max_nclosure);
 return(-1);
 }

/* Maximum number of frames: */
 sscanf(argv[4],"%d",&maxframes);

/* Erasing the arrays : */
    for(i = 0;i < nx * ny;i++) {
      long_int[i]=0.;
      modsq[i]=0.;
      snrm[i]=0.;
      mre[i]=0.;
      mim[i]=0.;
      }
    for(i=0; i < 4 * NGMAX; i++) { yce1[i]=0.;}

/* Computing the uv coverage: */
  printf(" Radius of uv-coverage (IR) in pixels: %d\n",ir);
  COVERA(&ir, &max_nclosure, &nbeta, &ngamma);
  if(ngamma > NGMAX) {
   fprintf(stderr,"Fatal error: ngamma=%d > NGMAX=%d\n", ngamma, NGMAX);
   exit(-1);
  }

/* Reading the file with the number of photons */
  JLP_INQUIFMT();
  printf(" Reading %s \n",nphfile);
  istat = JLP_VM_READIMAG1(&pntr,&nx1,&ny1,nphfile,nphcomments);
  if(istat) { 
     printf(" Fatal error accessing frame >%s< \n",nphfile);
     return(-1);
     }
  nphotons = (float *)pntr;
  nframes = nx1; 

  printf(" Max number of frames: %d \n",maxframes);
  nframes = (nframes > maxframes) ? maxframes : nframes;
  printf(" Will process %d frames \n",nframes);

/* Opening the data file */
  printf(" Reading %s \n",infile);
  fd1=open(infile,O_RDONLY);

/* Main loop processing all the frames included in the data file */
  xphotons = 0.;
  for(iframe = 1; iframe <= nframes; iframe++) {

/* Resetting the image for next step: */
    for(i = 0; i < nx * ny;i++) {image[i] = 0.; im[i] = 0.;}

/* Reading the data for image #iframe */
    nph=(int) *(nphotons+iframe-1);
    xphotons += (float)nph;
/* modulus 50 ...*/
     if((iframe % 50) == 1)
       printf(" Frame #%d photons: %d\n",iframe,nph);

/*
    nvalues=read(fd1,pacframe,(nph * 2 + 1) * sizeof(short));
* Now characters:
*/
    nvalues=read(fd1, pacframe, nph * 2 + 1);

#ifdef DEBUG
 printf(" Pacframe[0]: %d , nvalues %d \n",(int)pacframe[0],nvalues);
#endif

    if(nvalues == (nph * 2 + 1)) {
       for(i = 1; i <= nph; i++) {
/* Read the coordinates from pacfrme array: */
	  ix = (int)pacframe[2 * i - 1];
          iy = (int)pacframe[2 * i];
	  image[ix + iy * nx] += 1.;
	  }
/* Output the first image, to check : */
	  if(iframe == 1){
	  strcpy(outfile,"decode2_frame1");
	  sprintf(outcomments,"First image, nph:%d decode_simu2",nph);
          JLP_D_WRITEIMAG(image,&nx,&ny,&nx,outfile,outcomments); 
          }

/* Long integration : */
	  for(i = 0; i < nx * ny; i++) long_int[i] += image[i];

/* Fourrier Transform: */
          kod = 1;
          FFT_2D_DOUBLE(image, im, &nx, &ny, &nx, &kod);

/* Output the first image, to check : */
	  if(iframe == 1){
          RECENT_FFT_DOUBLE(image,image,&nx,&ny,&nx);
	  strcpy(outfile,"decode2_fft");
	  sprintf(outcomments,"First image, nph:%d decode_simu2",nph);
          JLP_D_WRITEIMAG(image,&nx,&ny,&nx,outfile,outcomments); 
          }
	  /*
	  printf(" image[0] %f \n",image[0]);
	  printf(" im[0] %f \n",im[0]);
          */

/* As this FFT (fourn1) divides by nx*ny we correct this back: 
	  w2 = (float)(nx * ny);
	  for(i=0; i < nx * ny; i++) 
          {image[i] = image[i] * w2;
           im[i] = im[i] * w2;
          }
	  */

/* Processing this image now:
*/
/*
int bispec3(double *re, double *im, double *modsq, double *snrm,
            INT4 *nx, INT4 *ny, double *bispp, INT4 *ir,
            INT4 *nbeta, INT4 *ngamma);
*/
          bispec3(image,im,modsq,snrm,&nx,&ny,yce1,&ir,
		   &nbeta,&ngamma);
       }
    else
       {printf(" Fatal error reading image #%d \n",iframe);
        printf(" Pacframe[0]: %d , nvalues %d, nph %d \n",
                 (int)pacframe[0],nvalues,nph);
	goto end1;}
   }

/* Recentre the frames: */
   RECENT_FFT_DOUBLE(modsq,modsq,&nx,&ny,&nx);
   RECENT_FFT_DOUBLE(snrm,snrm,&nx,&ny,&nx);
   RECENT_FFT_DOUBLE(mre,mre,&nx,&ny,&nx);
   RECENT_FFT_DOUBLE(mim,mim,&nx,&ny,&nx);

/* Computing the mean number of photons per frame: */
   xframes=(float)nframes;
   xphotons = xphotons/xframes;

/* Mean of the frames: */
   for(i = 0; i < nx * ny; i++) {

/* Photon noise correction in "jlp_bispec" (17-07-90) */
	modsq[i] /= xframes; 
	if(modsq[i] < 0.) modsq[i]=1.e-18;

/* SNR of modsq: */
	snrm[i] = snrm[i]/xframes - SQUARE(modsq[i]);
	if(snrm[i] <= 1.e-4) snrm[i]=1.e-4;
        snrm[i]= modsq[i]/sqrt((double)snrm[i]);

	mre[i]=mre[i]/xframes;
	mim[i]=mim[i]/xframes;
	long_int[i]=long_int[i]/xframes;
	}

   for(i = 0; i < ngamma; i++) { 
/* First computing the mean: */
        yce1[4*i] /= xframes;
        yce1[4*i + 1] /= xframes;
/* Sum of squares (real, imag): */
        yce1[4*i + 2] /= xframes;
        yce1[4*i + 3] /= xframes;
/* Then the variance:*/
        w1 =  yce1[4*i + 2] - SQUARE(yce1[4*i]);
        w2 =  yce1[4*i + 3] - SQUARE(yce1[4*i + 1]);
/* Then the sigma (real and imag together): */
        w1 = w1 + w2;
        if(w1 < 1.e-10) w1 = 1.e-10;
        w1 = sqrt((double)w1);
/* Phase factor of the bispectrum */
	w3 = SQUARE(yce1[4*i]) + SQUARE(yce1[4*i + 1]);
        w3 = sqrt((double)w3);
        if(w3 < 1.e-10) w3 = 1.e-10;
        yce1[4*i] /= w3;
        yce1[4*i + 1] /= w3;
/* SNR of bispectrum in 3rd line: */
        yce1[4*i + 2] = w3/w1;
        }

/* Now output of the results : */

/* Mean squared modulus : */
   strcpy(outfile,"modsq                               ");
   sprintf(outcomments,"simu2/ nfr: %d; nph: %f ",nframes,xphotons);
   JLP_D_WRITEIMAG(modsq,&nx,&ny,&nx,outfile,outcomments);

/* SNR of squared modulus : */
   strcpy(outfile,"snrm                               ");
   JLP_D_WRITEIMAG(snrm,&nx,&ny,&nx,outfile,outcomments);

/* Long integration : */
   strcpy(outfile,"long                                ");
   JLP_D_WRITEIMAG(long_int,&nx,&ny,&nx,outfile,outcomments);

/* Bispectrum : */
   strcpy(outfile,"bisp1                               ");
   ny=3;
   JLP_D_WRITEIMAG(yce1,&ngamma,&ny,&ngamma,outfile,outcomments);

  printf(" Successful end \n");

/* End : */
  end1:close(fd1);

return(0);
}
