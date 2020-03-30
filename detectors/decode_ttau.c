/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 decode_ttau.c
 Decode IR data from CFHT CIRCUS
 3-D Fits array

 JLP
 Version 06-01-93
---------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <jlp_ftoc.h>

#define DEBUG 1
#define MAXFILES 50 
/* For IR=30:
#define NGMAX 388400 
For IR=25:
*/
#define NGMAX 187566 
#define IDIM 32 
/* Circus camera: */
#define  NPLANES 600 
/* Note that DIM = IDIM * IDIM     */
#define DIM 1024 

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
/* Images IDIM*IDIM, irmax=25, ngmax=187566 : */

float image[DIM], modsq[DIM], snrm[DIM], mre[DIM], mim[DIM], long_int[DIM];
float select_array[NPLANES], yce1[4 * NGMAX], im[DIM], *pntr_image;
float xframes, w1, w2, w3;
long int  fd1, ivalues, nvalues, pntr, nph, idim=IDIM;
long int  ir, nbeta, ngamma, ngmax=NGMAX;
long int  nframes, iz_start[MAXFILES], iz_end[MAXFILES];
long int  iz1, nx, ny, kod, i, j, ifil, nfiles;
long int  istatus, dflag, vm_flag, nplanes;
unsigned  length1;
char infile[MAXFILES][61], comments[81], buffer[81], jlp_descr[200];
char file_list[41];
char file_ext[61], select_file[MAXFILES][61], outfile[41], outcomments[81], *pc;
FILE *fp;

  printf(" Program decode_ttau  Version 08-01-93 (%d x %d) \n",idim,idim);

/* Input parameters:  can be 4,5 or 6*/
if (argc != 4 && argc != 5 && argc != 6 )
  {
  printf("argc = %d \n",argc);
  printf("\nUSAGE:\n");
  printf("decode_ttau data_file Radius_uv_coverage iz_start,iz_end file_extension [selection_file] \n ");
  printf(" Example1: decode_ttau ttau1.fits 12 1,300 _t01 \n ");
  printf(" Example2: decode_ttau ttau1.fits 15 301,600 _s01 ttau1_sel_fits\n ");
  printf(" Example3 (with list of files in test.dat): decode_ttau test.dat 15 _r02\n ");
  exit(-1);
  }

/* Case of list of files: ****************************************/
 if(argc == 4)
  {
  strcpy(file_list,argv[1]);
  if((fp = fopen(file_list,"r")) == NULL)
    {printf(" Fatal error, cannot open file >%s< \n",file_list);
     exit(-1);
    }

/* Radius of uv-coverage: */
  sscanf(argv[2],"%d",&ir);

/* File extension: */
  strcpy(file_ext,argv[3]);
  pc = file_ext;
  while(*pc && *pc != ' ') pc++;
  *pc='\0';

/* Loop on the lines: */
  nfiles = -1;
  while(fgets(buffer,81,fp) != NULL)
  {
  nfiles++;
  printf("buffer= >%s< \n",buffer);
  sscanf(buffer,"%s %d,%d %s",infile[nfiles],&iz_start[nfiles],
                            &iz_end[nfiles],select_file[nfiles]);
  if(*select_file[nfiles])
    printf(" OK, selection file entered: >%s<\n",select_file[nfiles]);
  else
    {
    printf(" No selection file, will take all data within boundaries\n");
    }
  }

  fclose(fp);
  }
 else
  {
  nfiles = 0;
/* Data file name: */
  strcpy(infile[nfiles],argv[1]);

/* Radius of uv-coverage: */
  sscanf(argv[2],"%d",&ir);

/* iz_start, iz_end (working boundaries): */
  sscanf(argv[3],"%d,%d",&iz_start[nfiles],&iz_end[nfiles]);

/* File extension: */
  strcpy(file_ext,argv[4]);
  pc = file_ext;
  while(*pc && *pc != ' ') pc++;
  *pc='\0';

/* Selected frames: */
  if(argc > 5)
    {
    strcpy(select_file[nfiles],argv[5]);
    printf(" OK, selection file entered: >%s<\n",select_file[nfiles]);
    }
   else
    {
    *select_file[nfiles]='\0';
    printf(" No selection file, will take all data within boundaries\n");
    }
  } /* end of (argc != 3) */
  nfiles++;

/*****************************************************************/
   JLP_INQUIFMT();

/* Erasing the arrays : */
    for(i=0;i<DIM;i++) {
      long_int[i]=0.;
      modsq[i]=0.;
      snrm[i]=0.;
      mre[i]=0.;
      mim[i]=0.;
      }
    for(i=0; i < 4 * NGMAX; i++) { yce1[i]=0.;}

/* Computing the uv coverage: */
  printf(" Radius of uv-coverage (IR) in pixels: %d\n",ir);
  COVERA(&ir,&nbeta,&ngamma);

/* Maximum number of frames: */
  nframes = 0;
/***************************************************************/
/* Main loop on the frames: */
  for(ifil = 0; ifil < nfiles; ifil++)
  {
  printf(" Will process %s from %d to %d (inclusive)\n",
           infile[ifil],iz_start[ifil],iz_end[ifil]);

/* Computes the selection array (not selected if 0, selected otherwise) */
    if(*select_file[ifil])
/* If selection file reads directly the selection values : */
    { idim = NPLANES; 
     JLP_READIMAG(select_array,&nplanes,&i,&idim,select_file[ifil],comments); 
     if(nplanes != 600 || i != 1)
        {
        printf("Fatal error, nplanes= %d i= %d ... not Circus camera? \n",
                nplanes,i); 
        exit(-1);
        }
    }
/* Otherwise, assume that all values are correct: */ 
    else
    { nplanes = NPLANES;
      for(i = 0; i < nplanes; i++) select_array[i] = 1.;
    }

/* Reading the data file */
  for(iz1 = iz_start[ifil]-1 ; iz1 < iz_end[ifil]; iz1++) {
/* Internal index is from 0 to nplanes-1 for jlp0_rdfits... (so iz-1)*/

/* Works with selected frames only: */
    if(select_array[iz1] != 0.)
    {

    nframes++; 
    idim = IDIM; dflag = 0; vm_flag = 0;
/* dflag set to -1 to prevent error warning messages: */
    if(nframes > 1) dflag = -1;
    jlp0_rdfits(&pntr_image,image,&nx,&ny,&iz1,&idim,
                 infile[ifil],comments,jlp_descr,&dflag,&istatus,&vm_flag);
    if(istatus)
     {
     printf(" Fatal error reading image #%d \n",nframes);
     printf(" istatus = %d \n",istatus);
     nframes--;
     goto end1;
     }

/* Output index loop every 100 frames: */
    if((nframes % 100) == 1)
       printf(" Processing frame # %d (plane #%d)\n",nframes,iz1);

/* Output of the first image, to check : */
  if(nframes == 1)
  {
  strcpy(outfile,"decode_ttau_file1                         ");
  sprintf(outcomments,"First image (plane #%d) ",iz1);
  JLP_WRITEIMAG(image,&nx,&ny,&idim,outfile,outcomments); 
  }

/* Long integration : */
  for(i=0; i<DIM; i++) long_int[i]=long_int[i]+image[i];

/* Resetting the imaginary part for the FFT: */
  for(i=0;i<DIM;i++) {im[i] = 0.;}

/* Fourrier Transform: */
  kod=1;
  FFT_2D(image,im,&nx,&ny,&idim,&kod);

/* Processing this image now:  bispec1 is with photon noise correction
bispec2 is without*/
  /*
   BISPEC1(image,im,modsq,snrm,&nx,&ny,&idim,yce1,&ir,
	   &nbeta,&ngamma);
   */
   BISPEC2(image,im,modsq,snrm,&nx,&ny,&idim,yce1,&ir,
	   &nbeta,&ngamma);
   }
/* End of selected plane */
   }
/************************ End of loop with iz1 (all input frames processed) */
  }
/************************ End of loop with ifil (all input files processed) */

/* Check that some frames have been selected */
   xframes=(float)nframes;
   if(!nframes) 
     {printf("    Sorry no frame has been selected \n");
     goto end1;
     }
   else
     printf("    %d frames have been selected \n",nframes);

/* Recentre the frames: */
   RECENT_FFT(modsq,modsq,&nx,&ny,&idim);
   RECENT_FFT(snrm,snrm,&nx,&ny,&idim);
   RECENT_FFT(mre,mre,&nx,&ny,&idim);
   RECENT_FFT(mim,mim,&nx,&ny,&idim);


/* Mean of the frames: */
   for(i=0; i<DIM; i++) {

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
   sprintf(outcomments,"%s nfr: %d, from %d to %d, nfiles= %d",
           infile[0],nframes,iz_start[0],iz_end[0],nfiles);
   JLP_WRITEIMAG(modsq,&nx,&ny,&idim,outfile,outcomments);

/* SNR of squared modulus : */
   sprintf(outfile,"snrm%s",file_ext);
   JLP_WRITEIMAG(snrm,&nx,&ny,&idim,outfile,outcomments);

/* Long integration : */
   sprintf(outfile,"long%s",file_ext);
   JLP_WRITEIMAG(long_int,&nx,&ny,&idim,outfile,outcomments);

/* Bispectrum : */
   sprintf(outfile,"bisp1%s",file_ext);
   ny=3;
   JLP_WRITEIMAG(yce1,&ngamma,&ny,&ngmax,outfile,outcomments);

  printf(" Successful end \n");

/* End : */
  end1:close(fd1);
  JLP_END();
}
