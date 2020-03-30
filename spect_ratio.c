/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 spect_ratio.c
 To compute the spectral ratio from long integration and mean
 power spectrum 

 JLP
 Version 16-04-93
---------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <jlp_ftoc.h>

#define DEBUG 1

main(argc, argv)
int argc;
char **argv;
{
float *modsq, *long_int, *specrat, *im_work;
long int  nx, ny, pntr, isize, nxy;
register int i, j;
char modsq_name[61], long_name[61], comments[81], buffer[81];
char specrat_name[41], outcomments[81];

  printf(" Program spect_ratio  Version 16-04-93\n");

/* Input parameters:  should be 4*/
if (argc != 4)
  {
  printf("argc = %d \n",argc);
  printf("\nUSAGE:\n");
  printf("spect_ratio modsq long_int spectral_ratio ");
  exit(-1);
  }

/* File names: */
  strcpy(modsq_name,argv[1]);
  strcpy(long_name,argv[2]);
  strcpy(specrat_name,argv[3]);

/*****************************************************************/
   JLP_INQUIFMT();

/* Reading the input files */
  JLP_VM_READIMAG(&pntr,&nx,&ny,modsq_name,comments,&istatus);
  JLP_FROM_MADRID(&pntr,&modsq);
  JLP_VM_READIMAG(&pntr,&nx1,&ny1,long_name,comments,&istatus);
  JLP_FROM_MADRID(&pntr,&long_int);

/* Check size: */
  if(nx != nx1 || ny != ny1)
     {
     printf(" Error: incompatible image size (modsq and long_int! \n");
     istatus = 1;
     }

/* Allocation of memory: */
 isize = nx * ny * sizeof(float);
 JLP_GVM(&specrat,&isize);
 JLP_GVM(&im_work,&isize);

/* Erasing the arrays : */
    nxy = nx * ny;
    for(i = 0; i < nxy; i++) 
      {
      long_int[i]=0.;
      modsq[i]=0.;
      im_work[i]=0.;
      }

    strcpy(outcomments,"   ");
    JLP_WRITEIMAG(image,&nx,&ny,&nx,outfile,outcomments);
   }
#endif

/* Output index loop every 100 frames: */
    if((nframes % 100) == 1)
       printf(" Processing frame # %d \n",nframes);

/* Long integration : */
  for(i = 0; i < nxy; i++) long_int[i]=long_int[i]+image[i];

/* Resetting the imaginary part for the FFT: */
  for(i = 0; i < nxy;i++) {im_work[i] = 0.;}

/* Fourrier Transform: */
  kod=1;
  FFT_2D(image,im_work,&nx,&ny,&nx,&kod);

/* Processing this image now:  bispec1 is with photon noise correction
bispec2 is without*/
  /*
   BISPEC1(image,im_work,modsq,snrm,&nx,&ny,&nx,yce1,&ir,
	   &nbeta,&ngamma);
   */
   BISPEC2(image,im_work,modsq,snrm,&nx,&ny,&nx,yce1,&ir,
	   &nbeta,&ngamma);

/************************ End of loop with ifil (all input files processed) */
   JLP_FVM(&pntr_image);
   }

/* Check that some frames have been selected */
   xframes=(float)nframes;
   if(!nframes) 
     {printf("    Sorry no frame has been selected \n");
     goto end1;
     }
   else
     printf("    %d frames have been selected \n",nframes);

/* Recentre the frames: */
   RECENT_FFT(modsq,modsq,&nx,&ny,&nx);
   RECENT_FFT(snrm,snrm,&nx,&ny,&nx);
   RECENT_FFT(mre,mre,&nx,&ny,&nx);
   RECENT_FFT(mim,mim,&nx,&ny,&nx);


/* Mean of the frames: */
   for(i = 0; i < nxy; i++) {

/* SNR of modsq: */
	modsq[i]=modsq[i]/xframes;
	snrm[i]=snrm[i]/xframes - modsq[i]*modsq[i];
	if(snrm[i] <= 1.e-4) snrm[i]=1.e-4;
        snrm[i]= modsq[i]/sqrt((double)snrm[i]);

	mre[i]=mre[i]/xframes;
	mim[i]=mim[i]/xframes;
	long_int[i]=long_int[i]/xframes;
	}

   for(i=0; i<ngamma; i++) { 
/* First computing the mean: */
        yce1[i]=yce1[i]/xframes;
        yce1[NGMAX+i]=yce1[NGMAX+i]/xframes;
/* Sum of squares (real, imag): */
        yce1[2*NGMAX+i]=yce1[2*NGMAX+i]/xframes;
        yce1[3*NGMAX+i]=yce1[3*NGMAX+i]/xframes;
/* Then the variance:*/
        w1 =  yce1[2*NGMAX+i] - yce1[i]*yce1[i];
        w2 =  yce1[3*NGMAX+i] - yce1[NGMAX+i]*yce1[NGMAX+i];
/* Then the sigma (real and imag together): */
        w1 = w1 + w2;
        if(w1 < 1.e-10) w1 = 1.e-10;
        w1 = sqrt((double)w1);
/* Phase factor of the bispectrum */
	w3 = yce1[i]*yce1[i] + yce1[NGMAX+i]*yce1[NGMAX+i];
        w3 = sqrt((double)w3);
/* Normalization: */
/*
        yce1[i]=yce1[i]/w3;
        yce1[NGMAX+i]=yce1[NGMAX+i]/w3;
*/
/* SNR of bispectrum in 3rd line: */
        yce1[2*NGMAX+i]=w3/w1;
        }

/* Now output of the results : */

/* Mean squared modulus : */
   sprintf(outfile,"modsq%s",file_ext);
   if(argc == 5)
/* Case of a list: */
      sprintf(outcomments,"%s nfr: %d, nfiles= %d",
           infile[0],nframes,nfiles);
   else
      sprintf(outcomments,"%s nfr: %d, from %d to %d, nfiles= %d",
           infile[0],nframes,iz_start,iz_end,nfiles);

   JLP_WRITEIMAG(modsq,&nx,&ny,&nx,outfile,outcomments);

/* SNR of squared modulus : */
   sprintf(outfile,"snrm%s",file_ext);
   JLP_WRITEIMAG(snrm,&nx,&ny,&nx,outfile,outcomments);

/* Long integration : */
   sprintf(outfile,"long%s",file_ext);
   JLP_WRITEIMAG(long_int,&nx,&ny,&nx,outfile,outcomments);

/* Bispectrum : */
   sprintf(outfile,"bisp1%s",file_ext);
   ny=3;
   JLP_WRITEIMAG(yce1,&ngamma,&ny,&ngmax,outfile,outcomments);

  printf(" Successful end \n");

/* End : */
  end1:close(fd1);
  JLP_FVM(&long_int);
  JLP_FVM(&modsq);
  JLP_FVM(&snrm);
  JLP_FVM(&mre);
  JLP_FVM(&mim);
  JLP_FVM(&im_work);
  JLP_FVM(&apodi);
  JLP_END();
}
/* Creates Hamming apodization file: */
int jlp_hamming(apodi,nx,ny,idim)
int nx, ny,idim;
float *apodi;
{
double argx, argy;
double PI = 3.14159;
float apodi_y;
register int i, j;

printf(" nx=%d, ny=%d \n",nx,ny);

for(j = 0; j < ny; j++)
 {
   argy = PI * (j - ny/2)/ny; 
   apodi_y = 0.54 + 0.46 * cos(argy);
   for(i = 0; i < nx; i++)
    {
    argx = PI * (i - nx/2)/nx; 
    apodi[i + j * idim] = apodi_y * (0.54 + 0.46 * cos(argx));
    }
 }
  
return(0);
}
