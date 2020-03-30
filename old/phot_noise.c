/**************************************************************** 
 phot_noise.c
 Correction for photon noise, using long exposure and flat field.

 JLP
 Version 22-09-93
*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <jlp_ftoc.h>

main(argc, argv)
int argc;
char **argv;
{
/* Images : */
double *image1, *image2, *im;
float *modsq, *snrm, *long_int;
float *yce1, *in_ffield, *ffield, *apodi;
float xframes, w1, w2, w3, date, time, integ, xphot;
int npack, nphot_wanted, nphot_header, nphot_found, nframes;
int status, ireduc, ireduc_ff, isize, ixcent, iycent;
long int pntr_ffield, nx, ny, nx_ff, ny_ff;
long int ir, nbeta, ngamma, yce_dim;
register int i;
char *pc, file_ext[41], buffer[81], logfile[41];
char in_name[61], ff_name[61], ff_comments[81], outfile[61], outcomments[81];

printf(" Program phot_noise \n");
printf(" JLP Version 22-09-93 \n");


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
  printf(" runs phot_noise in_file_exten flat_field");
  exit(-1);
  }

/* Interactive input of in_name and file_ext: */
if (argc == 1)
 {
   printf(" Input file extension := "); scanf("%s",in_name);
/************* File extension for output files: */
   printf(" flat field := "); scanf("%s",file_ext);
 }
else
 {
  strcpy(file_ext,argv[2]);
  strcpy(in_name,argv[1]);
 }
  pc = file_ext;
  while(*pc && *pc != ' ') pc++;
  *pc='\0';

#ifdef DEBUG
printf(" DEBUG Version, will read >%s< \n",in_name);
#endif

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
  w1 = nphot_header/(50. * integ);
  w2 = nphot_header/integ;
  printf(" Average flux in 20 msec was: %f photons (i.e. %f photons/s)\n",
          w1,w2);

/***************************************************************/
JLP_BEGIN();
JLP_INQUIFMT();

/*****************************************************************/
/* Computing the spectral and bispectral lists
   corresponding to the selected uv coverage: */
  printf(" Radius of uv-coverage (IR) in pixels: %d\n",ir);
  COVERA(&ir,&nbeta,&ngamma);
#ifdef DEBUG
  printf(" nbeta: %d ngamma: %d\n",nbeta,ngamma);
#endif

isize = nx * ny * sizeof(double);
  JLP_GVM(&image1,&isize);
  JLP_GVM(&image2,&isize);
  JLP_GVM(&im,&isize);
isize = nx * ny * sizeof(float);
  JLP_GVM(&long_int,&isize);
  JLP_GVM(&modsq,&isize);
  JLP_GVM(&snrm,&isize);
isize = 4 * ngamma * sizeof(float);
  JLP_GVM(&yce1,&isize);

/* Erasing the arrays : */
    for(i = 0; i < nx * ny; i++) {
      image1[i] = 0.;
      image2[i] = 0.;
      im[i] = 0.;
      long_int[i] = 0.;
      modsq[i] = 0.;
      snrm[i] = 0.;
      }
    for(i=0; i < 4 * ngamma; i++) { yce1[i]=0.;}

/******************************************************************/
/* Read flat field file: */
  JLP_VM_READIMAG(&pntr_ffield,&nx_ff,&ny_ff,ff_name,ff_comments);
  JLP_FROM_MADRID(&pntr_ffield,&in_ffield);

  if(nx_ff != ny_ff || nx_ff < nx)
  {printf(" Fatal error/Ffield has wrong size: it should be square and larger than output image size)!\n");
   exit(-1);
  }
/* Reduction of the flat field to fit the image size: */
  ireduc_ff = nx_ff / nx;
  if(nx_ff != ireduc_ff * nx)
  {printf(" Fatal error/Ffield has wrong size: it should be a multiple of the output image size)!\n");
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

}
/********************* Unity file for ffield: *************/
else
{
  printf("\n OK, Take unity file for ffield \n");
 isize = nx * ny *sizeof(float);
 JLP_GVM(&ffield,&isize);
 for(i = 0; i < nx * ny; i++) ffield[i] = 1.;
}

/* Since we process twice the data (with a shift of half a frame: */
nframes_expected = 2 * nphot_wanted / npack;

/* Long integration: 
  status = rdcar_long(long_int,nx,ny,nx,in_name,nphot_wanted,&nphot_found);
*/
/* Photon processing: */
  status = rdcar_phot(image1,image2,im,modsq,snrm,yce1,long_int,
               ffield,nx,ny,nx,ixcent,iycent,ir,nbeta,ngamma,in_name,npack,
               nphot_wanted,&nphot_found,&nframes,nframes_expected);
  printf(" Output from rdcar: %d photons and %d frames processed\n",
          nphot_found, nframes);

/* Check that some frames have been processed */
   xframes=(float)nframes;

   xphot = sqrt((double)(modsq[0]/xframes));
  sprintf(outcomments,"%s %d-08-93 %.2fH %.1fs %.1fph %dfr %dph",
          in_name,(int)date,time,integ,xphot,nframes,nphot_found);


/* Mean of the frames: */
   for(i = 0; i < nx * ny; i++) {
        long_int[i]=long_int[i]/xframes;
        modsq[i]=modsq[i]/xframes;
        snrm[i]=snrm[i]/xframes;
        }
/****************************************************************/
   printf(" xphot = %f \n",xphot);
/* Normalizes FFT (for compatibility with FFT_2D instead of FFT_2D_FAST...*/
   jlp_normalize_fft(yce1,modsq,snrm,&nx,&ny,&nbeta,&ngamma);
   if(photon_correction) photon_corr(yce1,modsq,snrm,&nx,&ny,&xphot,&nbeta,&ngamma);
   yce_dim = ngamma;
   rearrange(yce1,&ngamma,&yce_dim);

/****************************************************************/
/* SNR of modsq: */
   for(i = 0; i < nx * ny; i++) {
        snrm[i] = modsq[i] / sqrt((double)snrm[i]);
        }

/****************************************************************/
/* Now output the results : */

/* Recentre the frames: */
   RECENT_FFT(modsq,modsq,&nx,&ny,&nx);
   RECENT_FFT(snrm,snrm,&nx,&ny,&nx);

/* Mean squared modulus : */
   sprintf(outfile,"modsq%s",file_ext);
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
   JLP_WRITEIMAG(yce1,&ngamma,&ny,&yce_dim,outfile,outcomments);

   printf(" Successful end \n");

/* End : */
end1:
  JLP_END();
  printf(" Output logfile: %s\n",logfile);
}
