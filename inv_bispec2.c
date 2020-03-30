/****************************************************************
* inv_bispec2.c
* Same as inv_bispec_30
* But use C routines for uv coverage (jlp_cover_mask.c)
* Warning: This prog does not use BMASK (---> check if this is right later!)
*
* PROGRAM INV_BISPEC.FOR : 2DPA (2D PHASE A)
* To compute the phase of the spectrum from the bispectrum,
* assuming a full pupil.
*
* Contains CREYCE_SIMU, CREYCE_OBSERV, CREYCE2, TRANS, RECURSIVE, EPLUS,
*          MODIF_WEIGHTS_SIGM1, MODIF_WEIGHTS_SIGM2
*          NOYAU, PROJ_EPLUS, ATRANS, ATRANS_A, CGRADIENT, ERROR, SORTIE
*          BISP_WEIGHT,BISP_WEIGHT2
*
* Syntax:
*
* Example1:
* RUNS INV_BISPEC 0 12,20,0.08,220,0.004,0.8 out_name real_fft ima_fft 
*
*   with 0=simulation option
*   12=radius of uv cover., 20 = weights or phase error, 
*   220=nbeta, 0.004=low_modulus, 0.8= max sigma
*   real_fft and ima_fft= FFT of image to be reconstructed
*
* Example2:
* RUNS INV_BISPEC 1 12,20,0.08,220,0.004,0.8 out_name modsq bisp1
*
*   with 1=real data, 
*   12=radius of uv cover., 20 = weights or phase error, 
*   220=nbeta, 0.004=low_modulus, 0.8 = max sigma
*   modsq=mean squared modulus of FFT
*   bisp1=mean bispectrum
*
* Example3:
* RUNS INV_BISPEC 2 12,20,0.08,220,0.004,0.8 out_name modsq real_fft,ima_fft bisp1
*
*   with 2=simulated data, 
*   12=radius of uv cover., 20 = weights or phase error, 
*   220=nbeta, 0.004=low_modulus, 0.8= max sigma
*   modsq=mean squared modulus of FFT (simulated) 
*   real_fft and ima_fft= FFT of image to be reconstructed
*   bisp1=mean bispectrum (simulated)
*
* INPUT files (in the case of real data, i.e., observations):
*  modsq: square modulus
*  bisp: bispectrum
*
* INPUT files (in the case of simulations):
*  modsq: square modulus
*  bisp: bispectrum
*  real_fft: real part of the spectrum of the object
*  ima_fft: imaginary part of the spectrum of the object
*
* OUPUT files:
*  out_name.SNR
*  out_name.SIG
*
* JLP
* Version: 25-01-00
*          (with comments added in April 2008)
* From A. Lannes P.FOR (January 1988)
* From inv_bispec1.for (February 1994)
*****************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jlp_ftoc.h>
#include "jlp_cover_mask.h" 

#ifndef PI
#define PI 3.141592654
#endif
/* Don't forget the parenthesis... */
#define DEGRAD (PI/180.)
/* Minimum value of sigma (for the central frequencies): */
/* This implies that SNR max = 1 / MIN_SIGMA */
#define MIN_SIGMA 1.E-2

/*
#define DEBUG
*/

/* For ir = 25, nbmax = 980  and ngmax = 187566
*  For ir = 30, nbmax = 1410 and ngmax = 388400
*  For ir = 35, nbmax = 1926 and ngmax = 724838
*  For ir = 40, nbmax = 2512 and ngmax = 1233164
*/
/*
#define IRMAX 25
#define NBMAX 980

#define IRMAX 30
#define NBMAX 1410

#define IRMAX 35
#define NBMAX 1926

#define IRMAX 40
#define NBMAX 2512
*/

/* For ir=30 :
#define IRMAX 30
#define MIRMAX -IRMAX
#define NBMAX 1410
#define NGMAX 388400
#define IDIM 256
*/
 
/* For ir=50 :  (to be checked) (successfully tested with 50,100 on 21/01/00)*/
#define IRMAX 50
#define MIRMAX -IRMAX
#define NBMAX 3922 
#define NGMAX 379368 
#define IDIM 256
 
int CREYCE_SIMU(float **yce_re, float **yce_im, float **xc_re, float **xc_im, 
                float **xcr_re, float **xcr_im, int **negative_modulus, 
                INT4 ir, float cte, INT4 nbeta, INT4 ngamma, INT4 ngamma_max, 
                float lower_modulus, char *real_name, char *imag_name, 
                float **weight, float **bmask);
int CREYCE_OBSERV(float **yce_re, float **yce_im, float **xc_re,
            float **xc_im, float **xcr_re, float **xcr_im, 
            int **negative_modulus, INT4 ir, float cte,
            INT4 nbeta, INT4 ngamma, INT4 ngamma_max,
            float lower_modulus, char *modsq_name, char *bisp_name,
            float **weight, float **bmask);
int CREYCE2(float **yce_re, float **yce_im, float **xc_re, float **xc_im, 
            float **xcr_re, float **xcr_im, int **negative_modulus, 
            INT4 ir, float cte, INT4 nbeta, INT4 ngamma, INT4 ngamma_max,
            float lower_modulus, char *modsq_name,
            char *real_name, char *imag_name, char *bisp_name,
            float **weight, float **bmask);
int CREYCE_LOAD_BISP(float **yce_re, float **yce_im, float *bisp1,
                        INT4 ngamma, INT4 nx1);
int CREYCE_BISP_WEIGHTS(INT4 nbeta, INT4 ngamma, INT4 ngamma_max, INT4 ir, 
                        float cte, float *bmask, float *ro, float *bisp1, 
                        float *weight);
int RECURSIVE(float *sigm, INT4 nbeta, INT4 ngamma, float sigma_null,
              INT4 ifermax, float *bmask, float *yce_re, float *yce_im,
              float *xc_re, float *xc_im, float *weight);
int BISP_WEIGHT11(INT4 nbeta, INT4 ngamma, INT4 ir, float cte,
                  float *bmask, float *ro, float *weight);
int BISP_WEIGHT0(INT4 ngamma, INT4 ir, float cte, float *weight);
int BISP_WEIGHT2(float *bisp1, INT4 nbeta, INT4 ngamma, INT4 ngamma_max, 
                 INT4 ir, float cte, float *weight);
int BISP_WEIGHT22(float *bisp1, INT4 nbeta, INT4 ngamma, INT4 ir,
                  float cte, float *bmask, float *weight);
int MODIF_WEIGHTS_SIGM2(float *sigm, INT4 nbeta, INT4 ngamma, float sig_max,
                        float lower_modulus, INT4 ifermax, float *bmask,
                        float *weight);
int MODIF_WEIGHTS_SIGM1(float *sigm, INT4 nbeta, INT4 ngamma, float sig_max,
                        float lower_modulus, INT4 ifermax,float *bmask,
                        float *weight);
int TRANS(INT4 ir, INT4 nbeta, float *bmask, float *xc_re, float *xc_im);
int COMPLEX_PROD2(float c1_re, float c1_im, float c2_re, float c2_im,
                  int sign1, int sign2, double *c12_re, double *c12_im);
int COMPLEX_PROD3(float c1_re, float c1_im, float c2_re, float c2_im,
                  float c3_re, float c3_im, int sign1, int sign2, int sign3,
                  double *c123_re, double *c123_im);
int OUTPUT_SNR1(float *sigm, INT4 nbeta, char *fname);
int ZERO_PHASE(float *xc_re, float *xc_im, INT4 nbeta);
int ERROR_SIMU(float *xc_re, float *xc_im, float *xcr_re, float *xcr_im,
               INT4 nbeta, float *err_spec, float *err_phas,
               char *errname, float *bmask);
int NORMALIZE_L1(float *array, INT4 npts, double *norm);
int NORMALIZE_L2(float *array, INT4 npts, double *norm);
int ATRANS(INT4 nbeta, INT4 ngamma, float *qmoy, INT4 ifermax,
           float *bmask, float *yce_re, float *yce_im,
           float *weight, float *xc_re, float *xc_im,
           float *xx, INT4 idim_xx);
int CGRADIENT(INT4 nbeta, INT4 ngamma, INT4 *it, float *qmoy,
              INT4 ifermax, float *bmask, float *yce_re, float *yce_im,
              float *weitht, float *xc_re, float *xc_im,
              float *xx, INT4 idim_xx);
int LSQUARES1(float exit_tolerance, INT4 nbeta, INT4 ngamma, INT4 ifermax,
              float *bmask, float *yce_re, float *yce_im,
              float *weight, float *xc_re, float *xc_im,
              float *xx, INT4 idim_xx, INT4 ittmax);
int SORTIE(float *xc_re, float *xc_im,
           INT4 nbeta, INT4 isortie, char *fname, float *bmask);
int PROJ_EPLUS(INT4 nbeta, INT4 n1, INT4 n2, float *xx,
               INT4 idim_xx, float *vecxy);
int NOYAU(INT4 nbeta, float *bmask, float *vecxy, INT4 idim_xx);
int EPLUS(INT4 nbeta, float *bmask, float *xc_re, float *xc_im, 
          float *xx, INT4 idim_xx, float *vecxy);
int FINAL_ERROR(INT4 ir, INT4 nbeta, float *final_err_spec,
                float *final_err_phas, char *errname, float *bmask,
                float *xcr_re, float *xcr_im, float *xc_re, float *xc_im,
                float *xx, INT4 idim_xx, float *vecxy, float *sigm);
int ATRANS_A(INT4 nbeta, INT4 ngamma, INT4 n1, INT4 n2, INT4 ifermax,
                float *bmask, float *weight, float *xx, INT4 idim_xx);
int jlp_atan2c(float *xxc, float cc_im, float cc_re);
int EIGEN_VALUE1(INT4 nbeta, INT4 ngamma, float *xlambda1, INT4 ifermax,
                 float *bmask, float *weight, float *xx, INT4 idim_xx, 
                 float *vecxy);
int EIGEN_VALUE2(INT4 nbeta, INT4 ngamma, float xlambda1, float *xlambda2, 
                 INT4 ifermax, float *bmask, float *weight, float *xx, 
                 INT4 idim_xx, float *vecxy);

/* Logfile: */
 static FILE *fp1;

/* Extension for output file names: */
 static char out_prefix[41];

/* xc: Spectrum (list): phase factor
* ro: Spectrum (list): modulus
* yce: Bispectrum (list)
* xcr: Reference spectrum (list), phase factor  (not known for real observ.)
* vecxy: Kernel vectors
*/
   static float *mask, *xx, *vecxy;
   static float *weight, *re, *im, *ro;
   static INT4 nx, ny;
   static INT4 ifermax = 10000;
/* ifermax: simply to give a limit to the number of closure relations to
* be taken into account (when this number would be really too big...)
* (Max. number of closure phase relations allowed for each pixel of
* the uv-coverage).
* IFERMAX is an internal upper limit, whereas NCLOSURE_MAX is set 
* when computing the data.
*/
int main(int argc, char *argv[])
{
 INT4 idim_xx, ir, nbeta, nbeta_max, istat;
 INT4 nx_mask, ny_mask, nclosure_max, ngamma, ngamma_max, ittmax;
 float *xc_re, *xc_im, *yce_re, *yce_im, *xcr_re, *xcr_im;
 float *bmask, *sigm, sig_max;
 float lower_modulus, cte, exit_tolerance, sigma_null;
 float ini_err_spec, ini_err_phas, final_err_spec, final_err_phas;
 int *negative_modulus;
 char fname[60], date[40], errname[40], logfile_name[40]; 
 char modsq_name[60], real_name[60], imag_name[60], bisp_name[60], *pc;
 int iopt, isize;
 register INT4 i;

printf(" inv_bispec2 \n");

/* Opening logfile: */
   strcpy(logfile_name,"inv_bispec2.log");
   if((fp1 = fopen(logfile_name,"w")) == NULL)
      {
       printf("Fatal error opening logfile %s \n",logfile_name);
       exit(-1);
      }

  JLP_CTIME(date,&istat);
  fprintf(fp1," Program inv_bispec2, version 10/03/2008, date=%s \n",date);
  printf(" Program inv_bispec2, version 10/03/2008, date=%s \n",date);

 
/* One or three parameters only are allowed to run the program: */
/* Carefull: 7 parameters always, using JLP "runs" */
#ifdef DEBUG
  printf(" argc=%d\n",argc);
  printf(" argv[4]=>%s<\n",argv[4]);
#endif


/* With "runs" argc == 7, otherwise argc should be 1 or 6: */
if (argc == 7 && argv[5]) argc = 6;
if(argc != 6 && argc != 1) 
  {
  printf("     Fatal error: Wrong syntax, argc=%d\n", argc);
  printf(" Syntax is:  \n");
  printf(" runs inv_bispec2 ioption parameters  out_prefix modsq bisp\n");
  printf(" Example (simulation): RUNS INV_BISPEC 0 12,20,0.08,220,0.004,0.8 tt_ real imag\n");
  printf(" Example (observations): RUNS INV_BISPEC 1 12,20,0.08,220,0.004,0.8 tt_ modsq bisp\n");
  printf(" Example (simulation): RUNS INV_BISPEC 2 12,20,0.08,220,0.004,0.8 tt_ modsq bisp real,imag\n");
  printf("(parameters: ir,cte,nbeta,lower_modulus,sig_max,sigma_null,nclosure_max)\n");
  exit(-1);
  }

*modsq_name = '\0';
*bisp_name = '\0';
*real_name = '\0';
*imag_name = '\0';

/* Interactive input of parameters: */
if (argc == 1)
 {
  printf(" Menu : \n 0=simulation (re, im) \n");
  printf(" 1= observations (mosq, bisp) \n");
  printf(" 2=simulation (modsq, bisp, re, im)\n");
  printf(" Enter the number of the option : ");
  scanf("%d",&iopt);

/* cte = angular reference for bispectral perturbation simulations (degrees) 
* Used for weight computation and weight selection...
* This reference is linked to the square norm of IR frequency
* Formerly: cte = 20.
* 
*/ 
  printf(" Radius (IR) of uv-coverage, phase error");
  printf("(0 if all weights=1, <0 if modulus weighted)");
  printf(" nbeta, lower_modulus, max sigma (for weights), sigma_value_when_modulus_is_null,");
  printf(" nclosure_max (used for input data computation)");
  printf("  (Ex: 12,20,0.1,220,0.,0.8,0.7,1000)");
  scanf("%d,%f,%d,%f,%f,%f,%d",&ir,&cte,
          &nbeta,&lower_modulus,&sig_max,&sigma_null,&nclosure_max);
/* Input file names: */
  printf(" Output extension for filenames:= "); scanf("%s",out_prefix);

  switch(iopt)
   {
/* iopt=0 Simulations with spectrum only: */
   case 0:
     printf(" Input real part:= "); scanf("%s",real_name);
     printf(" Input imaginary part:= "); scanf("%s",imag_name);
     break;
/* iopt=2 Simulations with all files: */
   case 2:
     printf(" Input real part:= "); scanf("%s",real_name);
     printf(" Input imaginary part:= "); scanf("%s",imag_name);
     printf(" Input square modulus file:= "); scanf("%s",modsq_name);
     printf(" Input bispectrum file:= "); scanf("%s",bisp_name);
     break;
/* iopt=1 or default: Real observations: */
   case 1:
   default:
     printf(" Input square modulus file:= "); scanf("%s",modsq_name);
     printf(" Input bispectrum file:= "); scanf("%s",bisp_name);
     break;
   }
 }
else
 {
  sscanf(argv[1],"%d",&iopt);
  sscanf(argv[2],"%d,%f,%d,%f,%f,%f,%d",&ir,&cte,
          &nbeta,&lower_modulus,&sig_max,&sigma_null,&nclosure_max);
/* Input file names: */
  sscanf(argv[3],"%s",out_prefix);
  printf(" argv[6]=%s argv[7]=%s \n", argv[6], argv[7]);
  switch(iopt)
   {
/* iopt=0 : simulations: */
/* Syntax: inv_bispec2 out_prefix iopt parameters real imag
*/
   case 0:
     strcpy(real_name,argv[4]);
     strcpy(imag_name,argv[5]);
     break;
/* iopt=2 : simulations with all files: */
/* Syntax: inv_bispec2 out_prefix iopt parameters modsq bisp real,imag 
*/
   case 2:
     strcpy(modsq_name,argv[4]);
     strcpy(bisp_name,argv[7]);
/* Only 7 parameters are allowed in the command line, hence real,imag files
* are loaded in the 7th parameter: */
     pc = argv[6];
     while(*pc != ',' && *pc) pc++;
     if(*pc == ',')
         {*pc = '\0'; 
          strcpy(real_name,argv[5]);
          pc++; 
          strcpy(imag_name,pc);}
     else{printf("Fatal error reading file names (imag_name)\n"); exit(-1);}
     printf("OK: real=%s imag=%s\n",real_name,imag_name);
     break;
/* iopt=1 or default: Real observations: */
/* Syntax: inv_bispec2 out_prefix iopt parameters modsq bisp 
*/
   case 1:
   default:
     strcpy(modsq_name,argv[4]);
     strcpy(bisp_name,argv[5]);
     break;
   }

 }
/* Removing blanks at the end: */
  pc = bisp_name; while(*pc && *pc != ' ') pc++; *pc='\0';
  pc = real_name; while(*pc && *pc != ' ') pc++; *pc='\0';
  pc = imag_name; while(*pc && *pc != ' ') pc++; *pc='\0';
  pc = modsq_name; while(*pc && *pc != ' ') pc++; *pc='\0';
  pc = out_prefix; while(*pc && *pc != ' ') pc++; *pc='\0';

/* Output to logfile: */
  fprintf(fp1," Option = %d \n",iopt);
  fprintf(fp1," Maximum number of closure relations allowed by the program: ifermax=%d \n",
          ifermax);
  printf(" Option = %d \n",iopt);
  printf(" Maximum number of closure relations allowed by the program: ifermax=%d \n",
          ifermax);
  fprintf(fp1,"ir: %d, cte: %6.2f, lower_modulus: %12.5e \n",
          ir, cte, lower_modulus);
  fprintf(fp1," Sig_max (for weights): %10.3e  Sig_null: %10.3e \n",
          sig_max,sigma_null);

#ifdef DEBUG
printf(" DEBUG Version, will read >%s< \n",modsq_name);
#endif

/*
#define SIMU
*/
#ifdef SIMU
   JLP_RANDOM_INIT(1);
#endif


/* Max. number of closure phase relations allowed for each pixel of
* the uv-coverage:
*/
   ifermax = NBMAX;
 
/* File size */
   printf(" WARNING: all the files should not have odd");
   printf(" numbers for NX and NY (size in X and Y) \n");
   JLP_INQUIFMT();

/****************************************************
* UV coverage computation: */
/* mask is the mask of accessible frequencies (here all frequencies): */
  nx_mask = 2 * ir + 4;
  ny_mask = nx_mask;
  isize = nx_mask * ny_mask * sizeof(float);
  mask = (float *)malloc(isize);
  for(i = 0; i < nx_mask * ny_mask; i++) mask[i] = 1.;
/* IFERMAX is an internal upper limit, whereas NCLOSURE_MAX is set 
* when computing the data:
*/
  if(nclosure_max < 0)
    {
    printf(" Fatal error: nclosure_max must be strictly positive! \n");
    fprintf(fp1," Fatal error: nclosure_max must be strictly positive! \n");
    fclose(fp1); exit(-1);
    }
  COVERA_MASK(mask,&nx_mask,&ny_mask,&ir,&nclosure_max,&nbeta_max,&ngamma_max);
/* For debug if necessary:
  output_lists_coverage(&nbeta_max, &ngamma_max);
*/
  free(mask);
 
  if(nbeta > nbeta_max)
   {
    fprintf(fp1,"inv_bispec2/Error: NBETA (wanted) = %d whereas NBETA_MAX = %d\n",
            nbeta,nbeta_max);
    fprintf(fp1," I correct it to BETA_MAX\n");
    printf("inv_bispec2/Error: NBETA (wanted) = %d whereas NBETA_MAX = %d\n",
            nbeta,nbeta_max);
    printf(" I correct it to BETA_MAX\n");
    nbeta = nbeta_max;
   }

/* Please note that ngamma = ngt(nbeta) */
     cover_ngt1(&ngamma,nbeta);
   printf(" NBETA_MAX (input data): %d NGAMMA_MAX (input data): %d \n",
            nbeta_max,ngamma_max);
   printf(" Maximum number of closure relations (input data): %d \n",
            nclosure_max); 
   printf(" NBETA (used here): %d, corresponding NGAMMA: %d \n",nbeta,ngamma);

   fprintf(fp1," NBETA_MAX (input data): %d NGAMMA_MAX (input data): %d \n",
            nbeta_max,ngamma_max);
   fprintf(fp1," Maximum number of closure relations (input data): %d \n",
            nclosure_max); 
   fprintf(fp1," NBETA (used here): %d, corresponding NGAMMA: %d \n",
            nbeta,ngamma);

/* Max number of iterations in the main loop: */
  ittmax = 4;
 
  fprintf(fp1,
   " Angular constant for bispectral noise estimation: %8.3e (degrees)\n",cte);
 
/* REAL VECXY(NBMAX+1,2) */
   isize = (nbeta + 1) * 2 * sizeof(float);
   vecxy = (float *)malloc(isize);
/* REAL X(NBMAX+1,4) */
   isize = (nbeta + 1) * 4 * sizeof(float);
   xx = (float *)malloc(isize);
 
/* First dimension of array xx  */
   idim_xx = nbeta + 1;

/* Creates YCE : bispectrum phasor,
* and initial set of weights */
/* iopt=0 : simulations */
  if(iopt == 0)
     {
     CREYCE_SIMU(&yce_re,&yce_im,&xc_re,&xc_im,&xcr_re,&xcr_im,
                 &negative_modulus,ir,cte,nbeta,ngamma,ngamma_max,
                 lower_modulus, real_name,imag_name,&weight,&bmask);
     strcpy(fname,real_name);
     }
/* iopt=1 or default: Real observations: */
  else if(iopt == 1)
     {
     CREYCE_OBSERV(&yce_re,&yce_im,&xc_re,&xc_im,&xcr_re,&xcr_im, 
             &negative_modulus,ir,cte,nbeta,ngamma,ngamma_max,
             lower_modulus,modsq_name,bisp_name,&weight,&bmask);
     strcpy(fname,modsq_name);
     }
/* iopt=2 : simulations with all files */
  else
     {
     CREYCE2(&yce_re,&yce_im,&xc_re,&xc_im,&xcr_re,&xcr_im,
             &negative_modulus,ir,cte,nbeta,ngamma,ngamma_max,
             lower_modulus,modsq_name,real_name,imag_name,bisp_name,
             &weight,&bmask);
     strcpy(fname,modsq_name);
     }

/******************************************************
* To initialize the solution,
* we solve the problem with the recursive method (Weigelt,...) 
* The initial spectral phase is stored in XC(.)
*/
     isize = (nbeta + 1) * sizeof(float);
     sigm = (float *)malloc(isize);
/* Compute the recursive solution and derives the standard
* deviation sigm of the spectral terms */
     RECURSIVE(sigm,nbeta,ngamma,sigma_null,ifermax,
               bmask,yce_re,yce_im,xc_re,xc_im,weight);
 
/* Modification of the weights, according to the recursive solution:
* JLP93:
*/
     if(cte == 0.)
         MODIF_WEIGHTS_SIGM1(sigm,nbeta,ngamma,sig_max,
              lower_modulus,ifermax,bmask,weight);
     else
         MODIF_WEIGHTS_SIGM2(sigm,nbeta,ngamma,sig_max,
              lower_modulus,ifermax,bmask,weight);

/* Output of SNR map: */
     OUTPUT_SNR1(sigm,nbeta,fname);

/* Output the errors of initial bispectrum: */
#ifdef DEBUG1
     strcpy(errname,"bisp_error1.dat");
     nbx = nbeta;
     nby = nbeta/2;
     memsize = nbx * nby * sizeof(float);
     work_array = (float *)malloc(isize);
     ERROR_BISPECT(nbeta,ngamma,errname,work_array,nbx,nby,
                   ifermax,yce_re,yce_im,xc_re,xc_im);
     strcpy(errname,"bisperr1");
     strcpy(comments,"Quadratic errors");
     JLP_WRITEIMAG(work_array,&nbx,&nby,&nbx,errname,comments);
     free(work_array);
#endif

/* Initial error estimation: */
 if(iopt != 1)
   {
   strcpy(errname,"error1.dat");
   ERROR_SIMU(xc_re,xc_im,xcr_re,xcr_im,nbeta,&ini_err_spec,
                 &ini_err_phas,errname,bmask);
   printf("Comparison with the model (since simulation) \n");
   printf("Initial rms error of the spectrum: %12.5e \n",ini_err_spec);
   printf("Initial rms error of the phase factor of the spect.: %12.5e \n",
           ini_err_phas);
   fprintf(fp1,"Comparison with the model (since simulation) \n");
   fprintf(fp1,"Initial rms error of the spectrum: %12.5e \n",ini_err_spec);
   fprintf(fp1,"Initial rms error of the phase factor of the spect.: %12.5e \n",
           ini_err_phas);
   }
 
/* Projection of the initial solution onto E^+ */
    EPLUS(nbeta,bmask,xc_re,xc_im,xx,idim_xx,vecxy);
 
/* Other alternative (Null phases)
    ZERO_PHASE(xc_re,xc_im,nbeta);
*/
 
/* SORTIE DES PARTIES REELLES ET IMAGINAIRES DE XC
* In rei and imi (i for initial)
*/
   SORTIE(xc_re,xc_im,nbeta,1,fname,bmask);

/*
* exit_tolerance = exit test for the main iteration (in degrees)
* This means that when the largest angular correction
* is smaller than "exit_tolerance", we exit from the main loop:
* Formerly:   exit_tolerance = 0.1
* From Aug 2008, I set this parameter internally: */
  exit_tolerance = 1.e-8;
  printf(" Exit tolerance in degrees (to exit from loop) : %12.3E \n",
          exit_tolerance);
  fprintf(fp1," Exit tolerance in degrees (to exit from loop) : %12.3E \n",
          exit_tolerance);

/* Main loop: least-squares non-linear fit */
printf(" JLP AUG08 ****************** \n");
   LSQUARES1(exit_tolerance,nbeta,ngamma,ifermax,
             bmask,yce_re,yce_im,weight,xc_re,xc_im,xx,idim_xx,ittmax);
 
/* Files ref et imf (f for final) */
   SORTIE(xc_re,xc_im,nbeta,2,fname,bmask);

/* Calage en translation */
      TRANS(ir,nbeta,bmask,xc_re,xc_im);

/* Final global errors and conclusions */
  if(iopt != 1)
    {
     strcpy(errname,"error2.dat");
     FINAL_ERROR(ir,nbeta,&final_err_spec,&final_err_phas,errname,bmask,
                 xcr_re,xcr_im,xc_re,xc_im,xx,idim_xx,vecxy,sigm);
 
     printf(" Comparison with the model (since simulation) \n");
     printf(" Initial & final rms error of the spectrum: %8.3f %8.3f\n",
             ini_err_spec,final_err_spec);
     printf(" Initial & final rms error of the phase factor (spect): %8.3f %8.3f\n",
             ini_err_phas,final_err_phas);
     fprintf(fp1," Comparison with the model (since simulation) \n");
     fprintf(fp1," Initial & final rms error of the spectrum: %8.3f %8.3f\n",
             ini_err_spec,final_err_spec);
     fprintf(fp1," Initial & final rms error of the phase factor (spect): %8.3f %8.3f\n",
             ini_err_phas,final_err_phas);
    }
 

/* Output the errors of final bispectrum: */
#ifdef DEBUG1
     strcpy(errname,"bisp_error2.dat");
     nbx = nbeta;
     nby = nbeta/2;
     memsize = nbx * nby * sizeof(float);
     work_array = (float *)malloc(isize);
     ERROR_BISPECT(nbeta,ngamma,errname,work_array,nbx,nby,
                   ifermax,yce_re,yce_im,xc_re,xc_im);
     strcpy(errname,"bisperr2");
     strcpy(comments,"Quadratic errors");
     JLP_WRITEIMAG(work_array,&nbx,&nby,&nbx,errname,comments);
     free(work_array);
#endif

/* Eigenvalues: */
#ifdef EIGEN_VALUES
   EIGEN_VALUE1(nbeta,ngamma,&xlambda1,ifermax,
                bmask,weight,xx,idim_xx,vecxy);
   EIGEN_VALUE2(nbeta,ngamma,xlambda1,&xlambda2,ifermax,
                bmask,weight,xx,idim_xx,vecxy);

   printf(" Largest eigen value: %12.5e \n",xlambda1);
   printf(" Smallest eigen value: %12.5e \n",xlambda2);
   fprintf(fp1," Largest eigen value: %12.5e \n",xlambda1);
   fprintf(fp1," Smallest eigen value: %12.5e \n",xlambda2);
/* Condition number: */
   if(xlambda2 < 1.e-12) xlambda2 = 1.e-12;
   printf(" Condition number: %12.5e \n",xlambda1/xlambda2);
   fprintf(fp1," Condition number: %12.5e \n",xlambda1/xlambda2);
#endif

 fclose(fp1);
 printf(" Log File in \"inv_bispec2.log\" \n");
 return(0);
}
/*******************************************************************
* CREYCE2
* Reads the modulus (squared) and bispectrum derived from simulations,
* re and im (available since it is a simulation)
*
* OUTPUT:
* weight: weights used for spectral list
*******************************************************************/
int CREYCE2(float **yce_re, float **yce_im, float **xc_re, float **xc_im, 
            float **xcr_re, float **xcr_im, int **negative_modulus, 
            INT4 ir, float cte, INT4 nbeta, INT4 ngamma, INT4 ngamma_max,
            float lower_modulus, char *modsq_name,
            char *real_name, char *imag_name, char *bisp_name,
            float **weight, float **bmask)
{
FILE *fp2;
double errmod, errbisp, errbispw, yy1_re, yy1_im;
float xw0, xw1, xw2, w1, xr, xi, xm, work, small_sqmodulus;
INT4 nx1, ny1, isize, inull, iix, iiy, ixc, iyc, irem;
INT_PNTR pntr;
int ir2, irad2;
/* nb and ng cannot be declared as "register" 
   since they are used as &nb or &ng in function calls... */
INT4 nb, ng, k, l, m;
float *modsq, *bisp1;
char comments[80];
register int i, j;
/*  COMPLEX YCE(*),YY1
*/
 
/* Input of the modulus: */
  printf(" Modulus (squared) of the FFT of the image (centered in the frame)\n");
  JLP_VM_READIMAG1(&pntr,&nx,&ny,modsq_name,comments);
  modsq = (float *)pntr;
  fprintf(fp1," Sq. modulus: %.14s  Comments: %.30s\n",modsq_name,comments);
 
/* Reading RE and IM : */
  printf(" Real part of the FFT of the image (centered in the frame)\n");
  JLP_VM_READIMAG1(&pntr,&nx,&ny,real_name,comments);
  re = (float *)pntr;
  fprintf(fp1," Real part of the FFT: %.14s  Comments: %.30s \n",
              real_name,comments);

  printf(" Imaginary part of the FFT of the image (centered in the frame)\n");
  JLP_VM_READIMAG1(&pntr,&nx,&ny,imag_name,comments);
  im = (float *)pntr;
  fprintf(fp1," Imag. part of the FFT: %.14s  Comments: %.30s \n",
              imag_name,comments);
 
/* Allocation of memory
*/
   isize = (nbeta + 1) * sizeof(float);
   *xc_re = (float *)malloc(isize);
   *xc_im = (float *)malloc(isize);
   *xcr_re = (float *)malloc(isize);
   *xcr_im = (float *)malloc(isize);
   *bmask = (float *)malloc(isize);
   *negative_modulus = (int *)malloc(isize);
/* ro is declared as a static array... */
   ro = (float *)malloc(isize);

/* Assume that zero frequency is at IXC,IYC:  */
   ixc = nx / 2;
   iyc = ny / 2;
   xw0 = SQUARE(re[ixc + iyc * nx]);
   xw1 = SQUARE(im[ixc + iyc * nx]); 
   xw1 = sqrt((double)(xw0 + xw1));
   xw2 = sqrt( (double)modsq[ixc + iyc * nx]);
   printf(" Central value of input and computed modulus: %12.5e %12.5e\n",
            xw2,xw1);
   fprintf(fp1," Central value of input and computed modulus: %12.5e %12.5e\n",
            xw2,xw1);

/* Look for the minimum value of the modulus: */
   small_sqmodulus = ABS(modsq[ixc + iyc * nx]);
   ir2 = ir * ir;
   for(j = 0; j < ny; j++) 
     for(i = 0; i < nx; i++) {
     irad2 = SQUARE(i - ixc) + SQUARE(j - iyc);
     if(irad2 < ir2) 
         small_sqmodulus = MINI(small_sqmodulus,ABS(modsq[i + j * nx]));
     }
   printf(" Smallest positive value of square modulus: %12.5e\n", 
             small_sqmodulus);
 
/* Opening error file for modulus: */
   if((fp2 = fopen("err_mod.dat","w")) == NULL)
      {
       printf("CREYCE2/Fatal error opening error file for modulus \"err_mod.dat\"\n");
       fprintf(fp1,"CREYCE2/Fatal error opening error file for modulus \"err_mod.dat\"\n");
       fclose(fp1); exit(-1);
      }
   fprintf(fp2,
" 0 0. 1. 1. Spectral list, modulus error, input mod., theoretical mod.\n");

/* Loading modulus and phasor to ro and xc */
     inull = 0;
     errmod = 0.;
     for(nb = 0; nb <= nbeta; nb++)
        {
/* Get coordinates (iix,iiy) from nb index: */
         COVER_IXY(&iix,&iiy,&nb);
/* Add (ixc,iyc) vector, since centered FFT... */
         iix += ixc;
         iiy += iyc;
         xr = re[iix + iiy * nx];
         xi = im[iix + iiy * nx];
         xm = sqrt((double)(xr * xr + xi * xi));
/* Phase factor: */
          if(xm == 0.)
            {
            (*xc_re)[nb] = 1.;
            (*xc_im)[nb] = 0.;
            }
          else
            {
            (*xc_re)[nb] = xr / xm;
            (*xc_im)[nb] = xi / xm;
            }
/* Modulus: */
          w1 = modsq[iix + iiy * nx];
/* JLP2008: DO NOT DISCARD NEGATIVE MODULII! */
/*
          if(w1 < 0.) w1 = 0.;
*/
          if(w1 <= 0.) {
              w1 = small_sqmodulus;
              inull++;
              (*negative_modulus)[nb] = 1;
           } else {
              (negative_modulus)[nb] = 0;
           }
          ro[nb] = sqrt((double)w1);
          work = xm/xw1 - ro[nb]/xw2;
          fprintf(fp2,"%d %f %f %f\n",nb,work,ro[nb]/xw2,xm/xw1);
          work *= work;
          errmod += work;

/* End of "for ... nb" loop */
    }

/* Close error file: */
    fclose(fp2);

    errmod = sqrt(errmod/(float)nbeta);
    printf(" sqerror of the modulus : %17.10e \n",errmod);
    if(inull > 0)
       {
       printf(" CREYCE2/Modulus negative or null for %d values\n",inull);
       fprintf(fp1," CREYCE2/Modulus negative or null for %d values\n",inull);
       }

    printf(" rms error of the modulus : %11.4e \n",errmod);
    fprintf(fp1," rms error of the modulus : %11.4e \n",errmod);
 
/* Output of some errors of the modulus: */
        for(nb = 0; nb < 5; nb++)
           {
/* Get coordinates (iix,iiy) from nb index: */
           COVER_IXY(&iix,&iiy,&nb);
/* Add (ixc,iyc) vector, since centered FFT... */
           iix += ixc;
           iiy += iyc;
           xr = re[iix + iiy * nx];
           xi = im[iix + iiy * nx];
           xm = sqrt((double)(xr * xr + xi * xi))/xw1;
           printf(" Index=%d, normalized RO and computed modulus: %11.4e %11.4e \n",
           nb, ro[nb]/xw2, xm);
           fprintf(fp1," Index=%d, normalized RO and computed modulus: %11.4e %11.4e \n",
           nb, ro[nb]/xw2, xm);
           }
 
/* Truncation: */
        irem = 0;
/* Mask to discard some values of the spectral list: */
        for(nb = 0; nb <= nbeta; nb++)
        {
/* Check that RO[NB] greater than LOWER_MODULUS:
* Remember that LOWER_MODULUS can be negative... */
          if(ro[nb] < lower_modulus || ro[nb] <= 0.)
            {
             (*bmask)[nb] = 0.;
             irem++;
            } else {
            (*bmask)[nb] = 1.;
            }
       }
 
        printf(" Number of terms of the spectral list after correction: %d \n",
                 nbeta-irem);
        fprintf(fp1," Number of terms of the spectral list after correction: %d \n",
                 nbeta-irem);

/* Spectrum phasor blocked in translation: xc[.] */
        TRANS(ir,nbeta,*bmask,*xc_re,*xc_im);

/* Reference (since simulation) */
        for(nb = 0; nb <= nbeta; nb++)
        {
          (*xcr_re)[nb] = (*xc_re)[nb];
          (*xcr_im)[nb] = (*xc_im)[nb];
/* Set initial solution XC to zero:
            (*xc_re)[nb] = 1.;
            (*xc_im)[nb] = 0.;
*/
        }

/***********************************************************************
* Bispectrum phasor: YCE
*/
        printf(" Bispectrum of the image (phase term) (ordered with bispectral list) \n");
        JLP_VM_READIMAG1(&pntr,&nx1,&ny1,bisp_name,comments);
        bisp1 = (float *) pntr;
        fprintf(fp1," Bispectrum: %.14s  Comments: %.30s\n",bisp_name,comments);
        if( (ny1 < 2) || (nx1 != ngamma_max))
          {
          printf(" CREYCE2/FATAL ERROR: Size of bispectrum inconsistent with IRMAX\n");
          fprintf(fp1," CREYCE2/FATAL ERROR: Size of bispectrum inconsistent with IRMAX\n");
          fclose(fp1); exit(-1);
          }

/* Load BISP array to YCE array: */
        CREYCE_LOAD_BISP(yce_re,yce_im,bisp1,ngamma,ngamma_max);

/* REAL WEIGHT(NGMAX) */
        isize = (ngamma + 1) * sizeof(float);
        *weight = (float *) malloc(isize);

/* Load BISP array to YCE array: */
        CREYCE_BISP_WEIGHTS(nbeta,ngamma,ngamma_max,ir,cte,*bmask,ro,bisp1,
                            *weight);
 
/* Compute the errors (since XCR is available): */
        errbisp = 0.;
        errbispw = 0.; 
        for(ng = 1; ng <= ngamma; ng++)
          {
           cover_klm1(&k,1,ng);
           cover_klm1(&l,2,ng);
           cover_klm1(&m,3,ng);
/* Theoretical bispectrum: yy1 = xcr[k]*xcr[l]*conjg(xcr[m]);
*/
/*
          yy1_re =   (*xcr_re)[k]*(*xcr_re)[l]*(*xcr_re)[m]
                   - (*xcr_im)[k]*(*xcr_im)[l]*(*xcr_re)[m]
                   + (*xcr_re)[k]*(*xcr_im)[l]*(*xcr_im)[m]
                   + (*xcr_im)[k]*(*xcr_re)[l]*(*xcr_im)[m];
          yy1_im =   (*xcr_re)[k]*(*xcr_im)[l]*(*xcr_re)[m]
                   + (*xcr_im)[k]*(*xcr_re)[l]*(*xcr_re)[m]
                   - (*xcr_re)[k]*(*xcr_re)[l]*(*xcr_im)[m]
                   + (*xcr_im)[k]*(*xcr_im)[l]*(*xcr_im)[m];
*/
         COMPLEX_PROD3((*xcr_re)[k],(*xcr_im)[k],(*xcr_re)[l],(*xcr_im)[l],
                       (*xcr_re)[m],(*xcr_im)[m],1,1,-1,&yy1_re,&yy1_im);
/* True error: */
          work = (yy1_re - (*yce_re)[ng])*(yy1_re - (*yce_re)[ng])
                 + (yy1_im - (*yce_im)[ng])*(yy1_im - (*yce_im)[ng]);
          errbisp += work;
          errbispw += (work * (*weight)[ng]);
          }

        errbisp = sqrt(errbisp/(double)ngamma);
        errbispw = sqrt(errbispw);
        printf(" rms error of the bispectrum : %11.4e \n",errbisp);
        printf(" weighted rms error of the bispectrum : %11.4e \n",errbispw);
        fprintf(fp1," rms error of the bispectrum : %11.4e \n",errbisp);
        fprintf(fp1," weighted rms error of the bispectrum : %11.4e \n",errbispw);
 
/* Just debug mode: */
    for(ng = 1; ng <= 5; ng++) {
/*  yy1=xcr(k)*xcr(l)*conjg(xcr(m)) */
       cover_klm1(&k,1,ng);
       cover_klm1(&l,2,ng);
       cover_klm1(&m,3,ng);
/*
       yy1_re =   (*xcr_re)[k]*(*xcr_re)[l]*(*xcr_re)[m]
                - (*xcr_im)[k]*(*xcr_im)[l]*(*xcr_re)[m]
                + (*xcr_re)[k]*(*xcr_im)[l]*(*xcr_im)[m]
                + (*xcr_im)[k]*(*xcr_re)[l]*(*xcr_im)[m];
       yy1_im =   (*xcr_re)[k]*(*xcr_im)[l]*(*xcr_re)[m]
                + (*xcr_im)[k]*(*xcr_re)[l]*(*xcr_re)[m]
                - (*xcr_re)[k]*(*xcr_re)[l]*(*xcr_im)[m]
                + (*xcr_im)[k]*(*xcr_im)[l]*(*xcr_im)[m];
*/
       COMPLEX_PROD3((*xcr_re)[k],(*xcr_im)[k],(*xcr_re)[l],(*xcr_im)[l],
                     (*xcr_re)[m],(*xcr_im)[m],1,1,-1,&yy1_re,&yy1_im);
#ifdef DEBUG
       printf(" ng=%d k=%d l=%d m=%d \n xcr_re[k]=%f ",ng,k,l,m,(*xcr_re)[k]);
       printf(" xcr_im[k]=%f xcr_re[l]=%f ",(*xcr_im)[k],(*xcr_re)[l]);
       printf(" xcr_im[l]=%f xcr_re[m]=%f xcr_im[m]=%f \n",
              (*xcr_im)[l],(*xcr_re)[m], (*xcr_im)[m]);
#endif
       printf(" Input, computed bisp & snr: %d (%8.5f,%8.5f) (%8.5f,%8.5f) %6.3f\n",
              ng, (*yce_re)[ng], (*yce_im)[ng], yy1_re, yy1_im,
              bisp1[ng - 1 + 2 * ngamma_max]);
       fprintf(fp1," Input, computed bisp & snr: %d (%8.5f,%8.5f) (%8.5f,%8.5f) %6.3f\n",
               ng, (*yce_re)[ng], (*yce_im)[ng], yy1_re, yy1_im,
               bisp1[ng - 1 + 2 * ngamma_max]);
      }
 
/* Free memory: */
free(bisp1);

return(0);
}
/*******************************************************************
* TRANS
* fait le calage en translation du facteur de phase spectral XC(.)
* En sortie, la phase de nb=1 est zero degre.
*
* A translation in the direct plane of (Dx, Dy) induces an increase of
* the phase of u_k by (2 pi (u_kx Dx + u_ky Dy)):
* TF[f(x - Dx)](u_kx, u_ky) = 
*                  exp(-2 i pi (u_kx Dx + u_ky Dy)) * TF[f(x)](u_kx,u_ky)
*
*  XC(1+1) = exp i Beta(1)
*
* Array TX is defined as:
*  TX(1) = exp -i Beta(1)
*  TX(-1)= exp i Beta(1)
*  TX(2) = exp -i 2 Beta(1)
*   ....
*  TX(N) = exp -i N Beta(1)
*
* OUTPUT:
*  If coordinates of frequency #NB are: (IXY1,IXY2)
*  XC[NB] = XC[NB] * TX(IXY1) * TY(IXY2)
*  i.e.:
*  exp i Beta_new = exp i Beta_old * exp -i IXY1 Beta(1) * exp -i IXY2 Beta(2)
*******************************************************************/
int TRANS(INT4 ir, INT4 nbeta, float *bmask, float *xc_re, float *xc_im)
{ 
float tx_re[IRMAX+1], tx_im[IRMAX+1], ty_re[IRMAX+1], ty_im[IRMAX+1]; 
float cx_re, cx_im, cy_re, cy_im;
double ccx_re, ccx_im, ccy_re, ccy_im, work_re, work_im;
INT4 iix, iiy, nb;
register int i;
 
    if(ir > IRMAX)
      {
        printf(" TRANS/Fatal error, maximum IR= %d \n",IRMAX);
        fprintf(fp1," TRANS/Fatal error, maximum IR= %d \n",IRMAX);
        fclose(fp1); exit(-1);
      }

/* Take the conjugate values of xc[1], xc[2] as a starting point, 
* to obtain a null phase  for xc[1] and xc[2]i at the end of this routine
* i.e., xc[1]=(1,0) and xc[2]=(1,0) 
*/
/* cx=conjg(xc(1+1)) */
  cx_re = xc_re[1];
  cx_im = -xc_im[1];
/* cy=conjg(xc(2+1)) */
  cy_re = xc_re[2];
  cy_im = -xc_im[2];
 
/* Tableau de la forme lineaire concerne du noyau */
  ccx_re = 1.; ccx_im = 0.;
  ccy_re = 1.; ccy_im = 0.;
  tx_re[0] = 1.; tx_im[0] = 0.;
  ty_re[0] = 1.; ty_im[0] = 0.;
 
   for(i = 1; i <= ir; i++)
     {
/* CCX=CCX*CX  Complex multiplication. */
       COMPLEX_PROD2(ccx_re,ccx_im,cx_re,cx_im,1,1,&ccx_re,&ccx_im);
       tx_re[i] = ccx_re;
       tx_im[i] = ccx_im;

/* CCY=CCY*CY  Complex multiplication. */
       COMPLEX_PROD2(ccy_re,ccy_im,cy_re,cy_im,1,1,&ccy_re,&ccy_im);
       ty_re[i] = ccy_re;
       ty_im[i] = ccy_im;
     }
 
/* Calage en translation */
     for(nb = 0; nb <= nbeta; nb++)
       {
/* XC(NB) = XC(NB) * TX(IXY(1,NB)) * TY(IXY(2,NB)) */
       COVER_IXY(&iix,&iiy,&nb);
              if(iix >= 0)
                {
                  cx_re = tx_re[iix];
                  cx_im = tx_im[iix];
                }
              else
                {
                  cx_re = tx_re[-iix];
                  cx_im = -tx_im[-iix];
                }
            cy_re = ty_re[iiy];
            cy_im = ty_im[iiy];

/* XC[NB] = XC[NB] * CX * CY; Complex multiplication... */
/*
            work_re =   xc_re[nb] * cx_re * cy_re
                      - xc_im[nb] * cx_im * cy_re
                      - xc_re[nb] * cx_im * cy_im
                      - xc_im[nb] * cx_re * cy_im;
            work_im =   xc_re[nb] * cx_im * cy_re
                      + xc_im[nb] * cx_re * cy_re
                      + xc_re[nb] * cx_re * cy_im
                      - xc_im[nb] * cx_im * cy_im;
*/
            COMPLEX_PROD3(xc_re[nb],xc_im[nb],cx_re,cx_im,cy_re,cy_im,
                          1,1,1, &work_re, &work_im);
            xc_re[nb] = work_re * bmask[nb]; 
            xc_im[nb] = work_im * bmask[nb]; 
      }
 
 return(0);
}
/*******************************************************************
* CREYCE_LOAD_BISP
* Load input data from file to YCE array
* yce index will be like fortran arrays: between 1 and ngamma
*******************************************************************/
int CREYCE_LOAD_BISP(float **yce_re, float **yce_im, float *bisp1,
                     INT4 ngamma, INT4 ngamma_max)
{
float xr, xi, xm;
INT4 ng, inull;
int isize;

/* COMPLEX YCE(NGMAX) */
  isize = (ngamma + 1) * 2 * sizeof(float);
  *yce_re = (float *)malloc(isize);
  *yce_im = (float *)malloc(isize);
  inull = 0.;
  for(ng = 0; ng < ngamma; ng ++)
  {
      xr = bisp1[ng];
      xi = bisp1[ng + ngamma_max];
/* Modulus */
      xm = sqrt((double)(xr*xr + xi*xi));
/* Phasor: */
      if(xm == 0.)
         {
         (*yce_re)[ng+1] = 1.;
         (*yce_im)[ng+1] = 0.;
         inull++;
/*         
         fprintf(fp1," CREYCE_LOAD_BISP/Warning: BISPEC(%d) is null! \n",ng); 
         printf(" CREYCE_LOAD_BISP/Warning: BISPEC(%d) is null! \n",ng); 
*/
         }
      else
         {
         (*yce_re)[ng+1] = xr / xm;
         (*yce_im)[ng+1] = xi / xm;
         }
    }
 
  if(inull > 0)
    {
    printf(" CREYCE_LOAD_BISP/Warning: null bispectrum for %d values \n",
           inull);
    fprintf(fp1," CREYCE_LOAD_BISP/Warning: null bispectrum for %d values \n",
           inull);
    }

/* Debug: */
  for(ng = 1; ng <= 5; ng ++)
    { 
     printf(" YCE(%d) = %12.5f %12.5f \n",ng,(*yce_re)[ng],(*yce_im)[ng]);
     fprintf(fp1," YCE(%d) = %12.5f %12.5f \n",ng,(*yce_re)[ng],(*yce_im)[ng]);
    }
 
return(0);
}
/**********************************************************************
* NORMALIZE_L1
* Norm 1 is the initial value of the sum of the absolute values of array.
* WARNING: it works with arrays from 1 to ngamma
*         Index from 1 to npts inclusive!!!
***********************************************************************/
int NORMALIZE_L1(float *array, INT4 npts, double *norm)
{
double ssum;
register int i;

 ssum=0.;
 for(i = 1; i <= npts; i++)
  {
  if(array[i] > 0)
    ssum += array[i];
  else
    ssum -= array[i];
  }

/* Normalisation of the array:
*/
    if(ssum == 0.)
      {
        printf(" NORMALIZE_L1/Fatal error: norm is null\n");
        fprintf(fp1," NORMALIZE_L1/Fatal error: norm is null\n");
        fclose(fp1); exit(-1);
      }

    for(i = 1; i <= npts; i++) array[i] /= ssum;

    *norm = ssum;

return(0);
}
/**********************************************************
* To normalize a vector with L2 norm
* Called to normalize the two vectors of the Kernel
* Norm 2 is the square root of the sum of the squared values of array.
* WARNING: index from 1 to npts inclusive!!!
**********************************************************/
int NORMALIZE_L2(float *array, INT4 npts, double *norm)
{
double ssum;
register int i;

 ssum=0.;
 for(i = 1; i <= npts; i++)
    ssum += (array[i] * array[i]);


/* Normalisation of the array:
*/
  if(ssum == 0.)
    {
      printf(" NORMALIZE_L2/Fatal error: norm is null\n");
      fprintf(fp1," NORMALIZE_L2/Fatal error: norm is null\n");
      fclose(fp1); exit(-1);
    }


 *norm = sqrt(ssum);
 for(i = 1; i <= npts; i++) array[i] /= *norm;

return(0);
}
/*********************************************************************
* Load the initial weights for computing the recursive solution
* Called by CREYCE2, CREYCE_SIMU and CREYCE_OBSERV
*
* OUTPUT:
* weight: weights to be used for computing the recursive solution
*********************************************************************/
int CREYCE_BISP_WEIGHTS(INT4 nbeta, INT4 ngamma, INT4 ngamma_max, INT4 ir, 
                        float cte, float *bmask, float *ro, float *bisp1, 
                        float *weight)
{
double xnorm;
INT4 k, kk, ng;

/* Computing the weights: */
        if(cte > 0)
          {
           BISP_WEIGHT11(nbeta,ngamma,ir,cte,bmask,ro,weight);
           NORMALIZE_L1(weight,ngamma,&xnorm);
          }
        else if (cte < 0)
          {
/* With SNR stored in 3rd line of bispectrum: */
          BISP_WEIGHT2(bisp1,nbeta,ngamma,ngamma_max,ir,cte,weight);
/* With modulus only: 
          BISP_WEIGHT22(bisp1,nbeta,ngamma,ir,cte,bmask,weight);
*/ 
/*        BISP_WEIGHT1(nbeta,ngamma,ir,cte,*bmask,ro,*weight);
*/
          NORMALIZE_L1(weight,ngamma,&xnorm);
          }
        else 
          {
           printf(" Weights set to unity, and then normalized \n");
           fprintf(fp1," Weights set to unity, and then normalized \n");
           for(ng = 1; ng <= ngamma; ng++)
              {
               weight[ng]=1.;
               for(kk = 1; kk <= 3; kk++) 
                  {
                  cover_klm1(&k,kk,ng);
                  if(bmask[k] == 0.) weight[ng]=0.;
                  }
              }
           NORMALIZE_L1(weight,ngamma,&xnorm);
          }
return(0);
}
/*****************************************************************
* BISP_WEIGHT0
* Computing the weights:
* JLP Version, not very good...
*
* OUTPUT:
* weight: weights used for the bispectral list
*****************************************************************/
int BISP_WEIGHT0(INT4 ngamma, INT4 ir, float cte, float *weight)
{
double work;
float cte1, sigb, sig2, wmin, wmax;
INT4 ng, k, kk, iix, iiy;

 cte1 = cte * DEGRAD / (float)(ir*ir);
 
 for(ng = 1; ng <= ngamma; ng++)
   {
        wmin = ir * ir;
        wmax = 0.;
        for(kk = 1; kk <= 3; kk++) 
         {
/*  K = KLM(KK,NG)*/
          cover_klm1(&k,kk,ng);
          COVER_IXY(&iix,&iiy,&k);
          work = iix*iix + iiy*iiy;
          if(wmin > work) wmin = work;
          if(wmax < work) wmax = work;
        }

/* To make the axes sensitive to the other dimension: 
*       wmin=max(wmin,0.1)
*        wmax=wmin*wmax
*/

/* sigb = ecart type (ou standard deviation) de gamma
*/
          sigb = cte1 * wmax;
 
/* Valeur moyenne du carre de [2 sin( (gamma - gamma barre)/2 ) ]
* cette valeur moyenne majore la precedente
*/
          work = -0.5*(sigb*sigb);
          sig2 = 1. - exp(work);
  
/* WEIGHT[NG] : (do not limit the weight values ...)
*/
         if(sig2 > 1.E-9)
              weight[ng] = 1./sig2;
         else
             {
             printf(" BISP_WEIGHT0/Warning: Null sigma for NG = %d \n",ng);
             fprintf(fp1," BISP_WEIGHT0/Warning: Null sigma for NG = %d \n",ng);
             weight[ng] = 1.E+9;
             }
    }
 
 printf(" WEIGHT0/New version. Maxi...\n");
 fprintf(fp1," WEIGHT0/New version. Maxi...\n");

return(0);
}
/*****************************************************************
* BISP_WEIGHT11
* Computing the weights:
* Old version (before 25-07-91) 
*
*****************************************************************/
int BISP_WEIGHT11(INT4 nbeta, INT4 ngamma, INT4 ir, float cte,
                  float *bmask, float *ro, float *weight)
{
double work;
float cte1, sigb, sig2;
INT4 ng, is2, k, l, m, kk, iix, iiy;

 cte1 = cte * DEGRAD / (float)(ir*ir);
 printf(" BISP_WEIGHT11/ cte1 = %f\n",cte1);
 
    for(ng = 1; ng <= ngamma; ng++)
      {
       is2 = 0;
       cover_klm1(&k,1,ng);
       cover_klm1(&l,2,ng);
       cover_klm1(&m,3,ng);
       work = ro[k] * ro[l] * ro[m];
       work *= (bmask[k] * bmask[l] * bmask[m]);
       if(work == 0)
           weight[ng]=0.;
       else
           {
            for(kk = 1; kk <= 3; kk++) 
             {
              cover_klm1(&k,kk,ng);
              COVER_IXY(&iix,&iiy,&k);
              is2 += iix*iix + iiy*iiy;
             }
/* sigb = ecart type (ou standard deviation) de gamma
*/
          sigb = cte1 * (float)is2;
 
/* Valeur moyenne du carre de [2 sin( (gamma - gamma barre)/2 ) ]
* cette valeur moyenne majore la precedente
*/
          work = -0.5*(sigb*sigb);
          sig2 = 1. - exp(work);
  
/* WEIGHT[NG] : (do not limit the weight values ...)
*/
         if(sig2 > 1.E-9)
              weight[ng] = 1./sig2;
         else
             {
             printf(" BISP_WEIGHT11/Warning: Null sigma for NG = %d (sigb = %f)\n",ng,sigb);
             fprintf(fp1," BISP_WEIGHT11/Warning: Null sigma for NG = %d (sigb = %f)\n",ng,sigb);
             weight[ng] = 1.E+9;
             }
        }
    }
 
 printf(" WEIGHT11/Old version.\n");
 fprintf(fp1," WEIGHT11/Old version.\n");

return(0);
}
/*****************************************************************
* BISP_WEIGHT2
* Computing the weights with SNR stored in the bispectrum file (line 3)
*
******************************************************************/
int BISP_WEIGHT2(float *bisp1, INT4 nbeta, INT4 ngamma, INT4 ngamma_max, 
                 INT4 ir, float cte, float *weight)
{
double sum;
float mean_snr; 
float lower_itt, upper_itt;
INT4 ng;

 auto_scale1(&bisp1[2*ngamma_max],ngamma_max,1,ngamma_max,
             &lower_itt,&upper_itt);

 sum = 0.;
 for(ng = 0; ng < ngamma; ng++) {
    weight[ng+1] = ABS(bisp1[ng + 2 * ngamma_max]);
/* OLD version:
    weight[ng+1] = bisp1[ng + 2 * ngamma_max];
    if(weight[ng+1] < 0) {
       weight[ng+1] = 0.;
       printf("JLPAUG08/ negative weight for ng=%d (bisp=%e)\n", 
               ng, bisp1[ng + 2 * ngamma]);
      }
*/
/* JLP93 : keep BISP_SNR, do not put any **0.5 **1.5 or **2, or anything else
*/
    sum += weight[ng+1];
    }
 
/* Diagnostic: */
  mean_snr = sum/(float)ngamma;

  printf(" BISP_WEIGHT2/SNR weight Initial sum: %12.5f Mean bisp. SNR: %12.5f \n",
          sum,mean_snr);
  fprintf(fp1," BISP_WEIGHT2/SNR weight Initial sum: %12.5f Mean bisp. SNR: %12.5f \n",
           sum,mean_snr);
return(0);
}
/*****************************************************************
* BISP_WEIGHT22
* Computing the weights taking the modulus into account
*
******************************************************************/
int BISP_WEIGHT22(float *bisp1, INT4 nbeta, INT4 ngamma, INT4 ir,
                  float cte, float *bmask, float *weight)
{
double sum;
float wrange, wmax, wmin, mean_snr; 
INT4 nb, ng, ng1, ng2;

/* JLP93: add modulus 
* First get the relative weights according to the modulus: 
*/
   wmax = 0.;
   wmin = ro[0];
   sum = 0.;
   ng1 = 1;
     for(nb = 3; nb <= nbeta; nb++) {
/*  IG2=NGT(NB)
*/
       cover_ngt1(&ng2,nb);
          for(ng = ng1; ng <= ng2; ng++)
            {
            weight[ng] = bmask[nb] * ro[nb];
            sum += weight[ng];
            if(weight[ng] != 0.) wmin = MINI(wmin, weight[ng]);
            wmax = MAXI(wmax, weight[ng]);
            }
         ng1 = ng2+1;
       }

    wrange=wmax/wmin;
    printf("BISP_WEIGHT22/Modulus: Initial range %f MIN, MAX: %f %f \n",
            wrange, wmin, wmax);
    fprintf(fp1,"BISP_WEIGHT22/Modulus: Initial range %f MIN, MAX: %f %f \n",
             wrange, wmin, wmax);

/* Truncation (max range to 100):
*/
    wmin = wmax / 100.;
    for(ng = 1; ng <= ngamma; ng++) weight[ng] = MAXI(wmin, weight[ng]);

/* Diagnostic: */
   mean_snr = sum/(float)ngamma;
   printf(" BISP_WEIGHT22/SNR weight Initial sum: %f, Mean bisp. SNR: %f \n",
            sum,mean_snr);
   fprintf(fp1," BISP_WEIGHT22/SNR weight Initial sum: %f, Mean bisp. SNR: %f \n",
            sum,mean_snr);

return(0);
}
/*******************************************************************
* CREYCE_SIMU : create the phasors of the spectrum and of the bispectrum
* Output in real/imaginary (from artificial bispectrum)
* Centers the Fourier transform
********************************************************************/
int CREYCE_SIMU(float **yce_re, float **yce_im, float **xc_re, float **xc_im, 
                float **xcr_re, float **xcr_im, int **negative_modulus, 
                INT4 ir, float cte, INT4 nbeta, INT4 ngamma, INT4 ngamma_max, 
                float lower_modulus, char *real_name, char *imag_name, 
                float **weight, float **bmask)
{
char comments[81];
INT4 ixc, iyc, nb, ng, inull, iix, iiy, is2;
INT4 irem, k, l, m, kk; 
float xr, xi, xm, cte1, work, work_re, work_im;
float cos1, sin1, dgamma;
double yre, yim;
INT_PNTR pntr;
int isize;

/* Reading RE and IM : */
  printf(" Real part of the FFT of the image (centered in the frame)\n");
  JLP_VM_READIMAG1(&pntr,&nx,&ny,real_name,comments);
  re = (float *)pntr;
  fprintf(fp1," Real part of the FFT: %.14s  Comments: %.30s \n",
              real_name,comments);

  printf(" Imaginary part of the FFT of the image (centered in the frame)\n");
  JLP_VM_READIMAG1(&pntr,&nx,&ny,imag_name,comments);
  im = (float *)pntr;
  fprintf(fp1," Imag. part of the FFT: %.14s  Comments: %.30s \n",
              imag_name,comments);
 
/* Allocation of memory: */
  isize = (nbeta + 1) * sizeof(float);
  *xc_re = (float *)malloc(isize);
  *xc_im = (float *)malloc(isize);
  *xcr_re = (float *)malloc(isize);
  *xcr_im = (float *)malloc(isize);
  *bmask = (float *)malloc(isize);
  *negative_modulus = (int *)malloc(isize);
/* ro is declared as a static array... */
  ro = (float *)malloc(isize);

/* Assume that zero frequency is at IXC,IYC:  */
   ixc = nx / 2;
   iyc = ny / 2;
 
/* Loading modulus and phasor to ro and xc */
     inull = 0;
     for(nb = 0; nb <= nbeta; nb++)
        {
/* Get coordinates (iix,iiy) from nb index: */
         COVER_IXY(&iix,&iiy,&nb);
/* Add (ixc,iyc) vector, since centered FFT... */
         iix += ixc;
         iiy += iyc;
         xr = re[iix + iiy * nx];
         xi = im[iix + iiy * nx];
/* Modulus: */
         xm = sqrt(xr * xr + xi * xi);
         ro[nb] = xm;
/* Phase factor: */
          if(xm == 0.)
            {
            (*xc_re)[nb] = 1.;
            (*xc_im)[nb] = 0.;
            (*negative_modulus)[nb] = 1;
            inull++;
            }
          else
            {
            (*xc_re)[nb] = xr / xm;
            (*xc_im)[nb] = xi / xm;
            (*negative_modulus)[nb] = 0;
            }

/* End of "for ... nb" loop */
    }

    if(inull > 0)
       {
       printf(" CREYCE_SIMU/Modulus null for %d values\n",inull);
       fprintf(fp1," CREYCE_SIMU/Modulus null for %d values\n",inull);
       }
/* Output of some values of the modulus: */
    for(nb = 0; nb < 4; nb++)
       {
       printf(" RO(%d) = %f \n",nb,ro[nb]);
       fprintf(fp1," RO(%d) = %f \n",nb,ro[nb]);
       printf(" XC(%d) = (RE, IM): %f %f\n",nb,(*xc_re)[nb],(*xc_im)[nb]);
       fprintf(fp1," XC(%d) = (RE, IM): %f %f\n",nb,(*xc_re)[nb],(*xc_im)[nb]);
       }
 
/* Truncation: */
        irem = 0;
/* Mask to discard some values of the spectral list: */
        for(nb = 0; nb <= nbeta; nb++)
        {
/* Check that RO[NB] greater than LOWER_MODULUS:
* Remember that LOWER_MODULUS can be negative... */
          if(ro[nb] < lower_modulus || ro[nb] <= 0.)
            {
            (*bmask)[nb] = 0.;
            irem++;
            }
          else
            (*bmask)[nb] = 1.;
       }
 
        printf(" Number of terms of the spectral list after correction: %d \n",
                 nbeta-irem);
        fprintf(fp1," Number of terms of the spectral list after correction: %d \n",
                 nbeta-irem);
 
/* Spectrum phasor blocked in translation: xc[.] */
        TRANS(ir,nbeta,*bmask,*xc_re,*xc_im);

/* Reference (since simulation) */
        for(nb = 0; nb <= nbeta; nb++)
        {
          (*xcr_re)[nb] = (*xc_re)[nb];
          (*xcr_im)[nb] = (*xc_im)[nb];
/* Set initial solution XC to zero:
            (*xc_re)[nb] = 1.;
            (*xc_im)[nb] = 0.;
*/
        }

/***********************************************************************
* Bispectrum phasor: YCE
*/
/* COMPLEX YCE(NGMAX) */
   isize = (ngamma + 1) * sizeof(float);
   *yce_re = (float *)malloc(isize);
   *yce_im = (float *)malloc(isize);
   *weight = (float *)malloc(isize);

/* Load BISP array to YCE array: */
    for(ng = 1; ng <= ngamma; ng++)
      {
/*  yce[ng]=xc(k)*xc(l)*conjg(xc(m)) */
       cover_klm1(&k,1,ng);
       cover_klm1(&l,2,ng);
       cover_klm1(&m,3,ng);
/*
       (*yce_re)[ng] =   (*xcr_re)[k]*(*xcr_re)[l]*(*xcr_re)[m]
                       - (*xcr_im)[k]*(*xcr_im)[l]*(*xcr_re)[m]
                       + (*xcr_re)[k]*(*xcr_im)[l]*(*xcr_im)[m]
                       + (*xcr_im)[k]*(*xcr_re)[l]*(*xcr_im)[m];
       (*yce_im)[ng] =   (*xcr_re)[k]*(*xcr_im)[l]*(*xcr_re)[m]
                       + (*xcr_im)[k]*(*xcr_re)[l]*(*xcr_re)[m]
                       - (*xcr_re)[k]*(*xcr_re)[l]*(*xcr_im)[m]
                       + (*xcr_im)[k]*(*xcr_im)[l]*(*xcr_im)[m];
*/
      COMPLEX_PROD3((*xcr_re)[k],(*xcr_im)[k],(*xcr_re)[l],(*xcr_im)[l],
                    (*xcr_re)[m],(*xcr_im)[m],1,1,-1, &yre, &yim);
      (*yce_re)[ng] = yre;
      (*yce_im)[ng] = yim;
       }
/* Just debug mode: */
     for(ng = 1; ng <= 5; ng++)
       {
        printf(" yce(%d) = (RE, IM): %f %f \n",ng,(*yce_re)[ng],(*yce_im)[ng]);
        fprintf(fp1," yce(%d) = (RE, IM): %f %f \n",ng,(*yce_re)[ng],(*yce_im)[ng]);
       }
 
/* Load the weights: */ 
     if(cte >= 0)
        CREYCE_BISP_WEIGHTS(nbeta,ngamma,ngamma_max,ir,cte,*bmask,ro,
                            NULL,*weight);
     else if (cte < 0)
       {
/* Version with SNR stored in 3rd line of bispectrum: */
        printf(" CREYCE_SIMU/Fatal error, this option (cte<0) is not possible here\n");
        fprintf(fp1," CREYCE_SIMU/Fatal error, this option (cte<0) is not possible here\n");
        fclose(fp1); exit(-1);
       }
 
/* Perturbation of YCE */
 cte1 = cte * DEGRAD / (float)(ir*ir);
 
 for(ng = 1; ng <= ngamma; ng++)
   {
     is2 = 0;
     for(kk = 1; kk <= 3; kk++) 
      {
/*  K = KLM(KK,NG)*/
       cover_klm1(&k,kk,ng);
/* IS2 = IS2 + IXY(1,K)**2 + IXY(2,K)**2
*/
       COVER_IXY(&iix,&iiy,&k);
       is2 += (iix*iix + iiy*iiy);
     }

/* Random generation (Gaussian law, (1.,0.)) of DGAMMA = DELTA GAMMA */
     JLP_RANDOM_GAUSS(&work);
     dgamma = work * cte1 * (float)is2;
     cos1=cos((double)dgamma);
     sin1=sin((double)dgamma);
/* yce[ng] = yce[ng]*cmplx(cos(dgamma),sin(dgamma)) */
     work_re = (*yce_re)[ng] * cos1 - (*yce_im)[ng] * sin1;
     work_im = (*yce_re)[ng] * sin1 + (*yce_im)[ng] * cos1;
     (*yce_re)[ng] = work_re;
     (*yce_im)[ng] = work_im;
  }
 
return(0);
}
/*******************************************************************
* CREYCE_OBSERV
* Reads a real spectrum (square modulus only) and bispectrum
* (Case of observations)
*
* OUTPUT:
*  bmask: mask applied to the spectral list to discard bad data
*******************************************************************/
int CREYCE_OBSERV(float **yce_re, float **yce_im, float **xc_re,
            float **xc_im, float **xcr_re, float **xcr_im, 
            int **negative_modulus, INT4 ir, float cte,
            INT4 nbeta, INT4 ngamma, INT4 ngamma_max,
            float lower_modulus, char *modsq_name, char *bisp_name,
            float **weight, float **bmask)
{
char comments[81];
INT4 ixc, iyc, nb, inull, iix, iiy, nx1, ny1;
INT4 irem, isize; 
int ir2, irad2;
float *modsq, *bisp1; 
double work, small_sqmodulus;
INT_PNTR pntr;
register int i, j;

/* Input of the modulus: */
  printf(" Modulus (squared) of the FFT of the image (centered in the frame)\n");
  JLP_VM_READIMAG1(&pntr,&nx,&ny,modsq_name,comments);
  modsq = (float *)pntr;
  fprintf(fp1," Sq. modulus: %.14s  Comments: %.30s\n",modsq_name,comments);
 
/* Allocation of memory: */
   isize = nx * ny * sizeof(float);
   re = (float *)malloc(isize);
   im = (float *)malloc(isize);
   isize = (nbeta + 1) * sizeof(float);
   *xc_re = (float *)malloc(isize);
   *xc_im = (float *)malloc(isize);
   *xcr_re = (float *)malloc(isize);
   *xcr_im = (float *)malloc(isize);
   *bmask = (float *)malloc(isize);
   *negative_modulus = (int *)malloc(isize);
/* ro is declared as a static array... */
   ro = (float *)malloc(isize);

/* Look for the minimum value of the modulus: */
   small_sqmodulus = ABS(modsq[0]);
   ir2 = ir * ir;
   for(j = 0; j < ny; j++) 
     for(i = 0; i < nx; i++) {
     irad2 = SQUARE(i - ixc) + SQUARE(j - iyc);
     if(irad2 < ir2) 
        small_sqmodulus = MINI(small_sqmodulus,ABS(modsq[i + j * nx]));
     }
   printf(" Smallest positive value of square modulus: %12.5e\n", 
            small_sqmodulus);

/* Initializing the real and imaginary part of the spectrum: */
  inull = 0;
  for(j = 0; j < ny; j++)
    {
    for(i = 0; i < nx; i++)
        {
           work = modsq[i + j * nx];
/* JLP2008: DO NOT DISCARD NEGATIVE MODULII ...
*/
           if(work < 0.) {
              work = small_sqmodulus;
              (*negative_modulus)[nb] = 1;
           } else {
              (*negative_modulus)[nb] = 0;
           }
           re[i + j * nx] = sqrt(work);
           im[i + j * nx] = 0.;
        }
    }
 
/* Assume that zero frequency is at IXC,IYC:  */
   ixc = nx / 2;
   iyc = ny / 2;
 
/* Loading modulus and phasor to ro and xc */
     for(nb = 0; nb <= nbeta; nb++)
        {
/* Get coordinates (iix,iiy) from nb index: */
         COVER_IXY(&iix,&iiy,&nb);
/* Add (ixc,iyc) vector, since centered FFT... */
         iix += ixc;
         iiy += iyc;
/* Modulus (here it is the real part) : */
         ro[nb] = re[iix + iiy * nx];
/* Phase factor: */
         (*xc_re)[nb] = 1.;
         (*xc_im)[nb] = 0.;
        }
 
    if(inull > 0)
       {
       printf(" CREYCE_OBSERV/Modulus null for %d values\n",inull);
       fprintf(fp1," CREYCE_OBSERV/Modulus null for %d values\n",inull);
       }

/* Output of some values of the modulus: */
#ifdef DEBUG
    for(nb = 0; nb < 5; nb++)
       {
       printf(" RO(%d) = %f \n",nb,ro[nb]);
       fprintf(fp1," RO(%d) = %f \n",nb,ro[nb]);
       }
#endif
 
/* Truncation: */
    irem = 0;
/* Mask to discard some values of the spectral list: */
    for(nb = 0; nb <= nbeta; nb++)
    {
/* Check that RO[NB] greater than LOWER_MODULUS:
* Remember that LOWER_MODULUS can be negative... */
      if(ro[nb] < lower_modulus || ro[nb] < 0.)
        {
         printf("CREYCE_OBSERV/Discarding value ro[%d]=%e (bmask)\n", 
                nb, ro[nb]);
         (*bmask)[nb] = 0.;
         irem++;
        }
      else
        (*bmask)[nb] = 1.;
   }
 
   printf(" Number of items of the spectral list after bmask selection: %d (%d removed)\n",
           nbeta-irem, irem);
   fprintf(fp1," Number of terms of the spectral list after bmask selection: %d(%d removed) \n",
            nbeta-irem, irem);
 
/* Spectrum phasor blocked in translation: xc[.] */
   TRANS(ir,nbeta,*bmask,*xc_re,*xc_im);

/***********************************************************************
* Bispectrum phasor: YCE
*/
    printf(" Bispectrum of the image (phase term) (NOT centered in the frame) \n");
/* Three lines (real, imag, snr) for the input bispectrum: */
    JLP_VM_READIMAG1(&pntr,&nx1,&ny1,bisp_name,comments);
    bisp1 = (float *)pntr;
    fprintf(fp1," Bispectrum: %.14s  Comments: %.30s\n",bisp_name,comments);
    if( (ny1 < 2) || (nx1 != ngamma_max))
      {
      printf("CREYCE_OBSERV/FATAL ERROR: Size of bispectrum inconsistent with IRMAX\n");
      fprintf(fp1,"CREYCE_OBSERV/FATAL ERROR: Size of bispectrum inconsistent with IRMAX\n");
      fclose(fp1); exit(-1);
      }

/* Load BISP array to YCE array: */
    CREYCE_LOAD_BISP(yce_re,yce_im,bisp1,ngamma,ngamma_max);

/* REAL WEIGHT(NGMAX) */
    isize = (ngamma + 1) * sizeof(float);
    *weight = (float *) malloc(isize);

/* Load BISP array to YCE array: */
    CREYCE_BISP_WEIGHTS(nbeta,ngamma,ngamma_max,ir,cte,*bmask,ro,bisp1,*weight);
 
/* Free memory: */
   free(bisp1);

return(0);
}
/*******************************************************************
* RECURSIVE: Initial solution of the phasor.
*                      YCE(.)  to  XC(.)
*
* exp(i*YCE[NG]) = exp (i * phase_K) * exp (i * phase_L) * exp (i * phase_M)
* with L=NB or M=NB
*
* INPUT:
*  yce_re, yce_im: bispectrum
*  weight: weights used for computing the recursive solution
*
* OUTPUT:
*  xc_re[nb], xc_im[nb]: spectrum corresponding to the recursive solution
*                        (distributed along the spectral list)
*  sigm[nb]: standard deviation of the spectral terms obtained with
*            the recursive solution (distributed along the spectral list)
*******************************************************************/
int RECURSIVE(float *sigm, INT4 nbeta, INT4 ngamma, float sigma_null,
              INT4 ifermax, float *bmask, float *yce_re, float *yce_im,
              float *xc_re, float *xc_im, float *weight)
{
double cc1_re, cc1_im;
/* Double precision is necessary for large sums !!!
*/
double sumxr1, sumxi1, sumsqr, sumsqi, sigr, sigi, sumweights, rnorm;
float xw;
INT4 k, l, m, ng1, ng2, nb, ng, indx, nval, inull;
/* JLP96: */
FILE *fp2;
fp2 = fopen("recurs.dat","w");
 
xc_re[0] = 1.;
xc_im[0] = 0.;
xc_re[1] = 1.;
xc_im[1] = 0.;
xc_re[2] = 1.;
xc_im[2] = 0.;
/* In the central part of the spectrum, good SNR...
*/
sigm[0] = MIN_SIGMA;
sigm[1] = MIN_SIGMA;
sigm[2] = MIN_SIGMA;

ng1 = 1;
inull = 0;
 
for(nb = 3; nb <= nbeta; nb++)
 {
 sumxi1 = 0.;
 sumxr1 = 0.;
 sumsqr = 0.;
 sumsqi = 0.;
 sumweights = 0.;
/* JLP94: sigma has a meaning only if more than 2 values are involved in the
* mean... 
*/
 nval = 0;
/* NG2=NGT(NB) */
 cover_ngt1(&ng2,nb);
 indx = 0;
/* For each nb index, the bispectral list is contained in the range ng1, ng2:
*/
 for(ng = ng1; ng <= ng2; ng++)
   {
   indx++;
   if(indx <= ifermax)
     {
/*  K=KLM(1,NG) L=KLM(2,NG) M=KLM(3,NG)
*/
      cover_klm1(&k,1,ng);
      cover_klm1(&l,2,ng);
      cover_klm1(&m,3,ng);
/* JLPAUG08: I replace sqrt(weight) with weight here
* (and put SQUARE(weight) in ATRANS, ... */
      xw = weight[ng];
      if(xw > 0.)
         {
          if(m == nb)
             {
/* Case 1 (m = nb): cc1=xc(k)*xc(l)*conjg(yce[ng]) */
             COMPLEX_PROD3(xc_re[k],xc_im[k],xc_re[l],xc_im[l],
                  yce_re[ng],yce_im[ng],1,1,-1,&cc1_re,&cc1_im);
             }
          else
             {
/* Case 2 (l = nb): cc1=conjg(xc(k))*xc(m)*yce[ng] */
             COMPLEX_PROD3(xc_re[k],xc_im[k],xc_re[m],xc_im[m],
                    yce_re[ng],yce_im[ng],-1,1,1,&cc1_re,&cc1_im);
             }
          nval++;
          sumsqr += (cc1_re * cc1_re * xw);
          sumsqi += (cc1_im * cc1_im * xw);
          sumxr1 += (cc1_re * xw);
          sumxi1 += (cc1_im * xw);
          sumweights += xw;
/* End of test on xw (weight) */
         }
/* End of test on ifermax  (max number of closure relations) */
       }
/* End of ng loop */
    }
 
/* Now since the phase of xc[nb] has been found, we normalize the phasor: */
    rnorm = sqrt(sumxr1 * sumxr1 + sumxi1 * sumxi1);
    if(rnorm > 0.)
       {
        xc_re[nb] = sumxr1/rnorm;
        xc_im[nb] = sumxi1/rnorm;

/* JLP94: sigma has a meaning only if more than 3 values are involved in the
* mean... */
        if(nval >= 3)
          {
/* SIGR**2 = Sum of (weight_i *(X_i - mean)**2) / Sum of (weight_i) 
* or Sum of (weight_i X_i**2)/Sum of weight_i - mean**2
* Here the weighted mean of real(CC) is SUMXR1/SUMWEIGHTS:
* (Checked again that's OK in 1996)
*/
           sigr = sumsqr / sumweights - SQUARE(sumxr1 / sumweights);
           sigi = sumsqi / sumweights - SQUARE(sumxi1 / sumweights);
/* Note that double precision variables are needed for large sums !!!
*/
           if(sigr < 0. || sigi < 0.) {
/* Just in case of a (round-off??) problem:
*/
             printf("\n\n");
             printf(" Round-off problem in RECURSIVE... for NB= %d \n",nb);
             fprintf(fp1," Round-off problem in RECURSIVE... for NB= %d \n",nb);
             printf(" sigi: %f sigr: %f \n",sigi,sigr);
             fprintf(fp1," sigi: %f sigr: %f \n",sigi,sigr);
             printf(" rnorm: %f sumweights: %f \n",rnorm,sumweights);
             printf(" sumsqr: %f sumsqi: %f \n",sumsqr,sumsqi);
             printf(" sumxr1: %f sumxi1: %f \n",sumxr1,sumxi1);
             printf(" ng1: %d ng2: %d \n",ng1,ng2);
             printf(" weight(ng1): %f  weight(ng2): %f\n",
                      weight[ng1],weight[ng2]);
             printf("\n\n");
/* JLP97: I comment this out:
             fclose(fp1); exit(-1);
 and replace with:
*/
              sigm[nb] = 0.;
           } else {
              sigm[nb] = sqrt(sigi+sigr);
/* JLP96: */
              fprintf(fp2,"%d phase:|(%.2f,%.2f)|=%f sig:|%.3f,%.3f|=%f nval=%d\n",
                  nb,xc_re[nb],xc_im[nb],rnorm,sqrt(sigr),sqrt(sigi),sigm[nb],nval);
           }
/* End of sigr < 0 or sigi < 0 */
/* If less than 3 values were used for the mean computation:
* Assume that it is the same as for the previous spectral term
*/
       } else {
          printf("JLPPP: less than 3 values for nb=%d (nval=%d)\n", nb, nval);
          sigm[nb] = sigm[nb-1];
       }
/* rnorm = 0, i.e., the recursive process has been interrupted: */
     } else {
       inull++;
/* We arbitrarily set the phase to zero: */
       xc_re[nb] = 1.;
       xc_im[nb] = 0.;
/* Maximum sigma is one, but I set it to SIGMA_NULL since we assume that
* phase is not important if modulus is too small
       printf("JLPPP: recursive process interrupt for nb=%d (nval=%d)\n", 
               nb, nval);
       printf(" rnorm: %f sumweights: %f \n",rnorm,sumweights);
       printf(" sumsqr: %f sumsqi: %f \n",sumsqr,sumsqi);
       printf(" sumxr1: %f sumxi1: %f \n",sumxr1,sumxi1);
       printf(" ng1: %d ng2: %d  (ngamma=%d)\n", ng1, ng2, ngamma);
       printf(" weight(ng1): %f  weight(ng2): %f\n", weight[ng1],weight[ng2]);
       printf("\n\n");
*/
       sigm[nb] = sigma_null;
       sigm[nb] = 1.;
/* Allows breaks only if in the list of null modulus: */
       if(bmask[nb] != 0.) {
/* Then add this frequency to BMASK, since it won't be able to recover the
* phase of this frequency:
*/
           bmask[nb] = 0.;
           printf(" RECURSIVE/Warning: too many null modulii: the phase of xc[%d] will remain undetermined.\n",nb);  
           fprintf(fp1," RECURSIVE/Warning: too many null modulii: the phase of xc[%d] will remain undetermined.\n",nb);  
#ifdef DEBUG
/* Print the values of the modulii that may be null: */
           printf(" ng1=%d, ng2=%d sumxr1=%e sumxi1=%e rnorm=%e\n",
                       ng1,ng2,sumxr1,sumxi1,rnorm);
           for(ng = ng1; ng <= ng2; ng++) {
                cover_klm1(&k,1,ng);
                cover_klm1(&l,2,ng);
                cover_klm1(&m,3,ng);
                printf(" ng=%d, k=%d, l=%d, m=%d \n",ng,k,l,m);
                printf(" ro[k]=%.4e, ro[l]=%.4e, ro[m]=%.4e weight[ng]=%.4e\n",
                   ro[k],ro[l],ro[m],weight[ng]);
           }
#endif
        }
/* End of bmask != 0 */
     }
/* End of interruption of recursive process */
   ng1 = ng2+1;
/* End of nb loop */
   }

/* JLP96: */
  fclose(fp2);
 
/* Set the first 4 values to zero, to avoid problems... */
sigm[0] = MIN_SIGMA;
sigm[1] = MIN_SIGMA;
sigm[2] = MIN_SIGMA;
sigm[3] = MIN_SIGMA;
sigm[4] = MIN_SIGMA;

/* Applies a threshold on sigm, to avoid SIGM = 0, 
* and problems when computing 1/SIGM in SNR map:
*/
for(nb = 5; nb <= nbeta; nb++) sigm[nb] = MAXI(MIN_SIGMA, sigm[nb]);

/* Diagnostic: */
 if(inull > 0)
   {
   printf(" RECURSIVE/Warning: Null modulus (i.e., the spectral list has been reduced) for %d values \n",inull);
   fprintf(fp1," RECURSIVE/Warning: Null modulus (i.e., the spectral list has been reduced) for %d values \n",inull);
   printf(" WARNING: In that case SIGMA was set to %f (and SNR to %f) \n",
           sigma_null, 1./sigma_null); 
   fprintf(fp1," WARNING: In that case SIGMA was set to %f (and SNR to %f) \n",
           sigma_null, 1./sigma_null); 
   }
 
return(0);
}
/**************************************************************
* COMPLEX_PROD2
* Complex product of two complex values
* INPUT:
* c1_re, c1_im, c2_re, c2_im
* sign1, sign2 ("sign" of imaginary part: 1 if simple product, -1 if conjugate)
* OUTPUT:
* (c1 * c2)_re
* (c1 * c2)_im
**************************************************************/
int COMPLEX_PROD2(float c1_re, float c1_im, float c2_re, float c2_im,
                  int sign1, int sign2, double *c12_re, double *c12_im)
{
/* The values are modified only in this routine (not outside) */
  c1_im *= (float)sign1;
  c2_im *= (float)sign2;

/* c12 = c1 * c2: Complex multiplication... */
 *c12_re =   c1_re * c2_re - c1_im * c2_im; 
 *c12_im =   c1_re * c2_im + c1_im * c2_re;
return(0);
}
/**************************************************************
* COMPLEX_PROD3
* Complex product of three complex values
* INPUT:
* c1_re, c1_im, c2_re, c2_im, c3_re, c3_im
* sign1, sign2, sign3 ("sign" of imaginary part: 
*                             1 if simple product, -1 if conjugate)
* OUTPUT:
* (c1 * c2 * c3)_re
* (c1 * c2 * c3)_im
**************************************************************/
int COMPLEX_PROD3(float c1_re, float c1_im, float c2_re, float c2_im,
                  float c3_re, float c3_im, int sign1, int sign2, int sign3,
                  double *c123_re, double *c123_im)
{
/* The values are modified only in this routine (not outside) */
  c1_im *= (float)sign1;
  c2_im *= (float)sign2;
  c3_im *= (float)sign3;

/* c123 = c1 * c2 * c3: Complex multiplication... */
 *c123_re =   c1_re * c2_re * c3_re
            - c1_im * c2_im * c3_re
            - c1_re * c2_im * c3_im
            - c1_im * c2_re * c3_im;
 *c123_im =   c1_re * c2_im * c3_re
            + c1_im * c2_re * c3_re
            + c1_re * c2_re * c3_im
            - c1_im * c2_im * c3_im;
return(0);
}
/*******************************************************************
* MODIF_WEIGHTS_SIGM2: to modify the weights according to the sigma found
* when computing the recursive solution
* SIGM2 is specially designed with weights according to snr_bispectrum
*
* INPUT:
* sig_max: Upper threshold of sigma 
*          (truncation used to reduce the range of weights)
*
* OUTPUT
* weight: weights used for the bispectral list
*******************************************************************/
int MODIF_WEIGHTS_SIGM2(float *sigm, INT4 nbeta, INT4 ngamma, float sig_max,
                        float lower_modulus, INT4 ifermax, float *bmask,
                        float *weight)
{
FILE *fp2;
INT4 ng1, ng2, ng, nb, irem;
double xnorm;
 
ng1 = 1;
irem = 0;
 
/* Opening sigma file: */
if((fp2 = fopen("sigma.dat","w")) == NULL)
   {
   printf("MODIF_WEIGHTS_SIGM2/Fatal error opening sigma file \"sigma.dat\"\n");
   fprintf(fp1,"MODIF_WEIGHTS_SIGM2/Fatal error opening sigma file \"sigma.dat\"\n");
   fclose(fp1); exit(-1);
   }

for(nb = 3; nb <= nbeta; nb++)
   {
   cover_ngt1(&ng2,nb);
   fprintf(fp2,"%d %f \n",nb,sigm[nb]);

/* Truncation to reduce the range of weights: */
   if(sigm[nb] >= sig_max || bmask[nb] == 0.)
      {
      irem++;
      bmask[nb] = 0.;
      for(ng = ng1; ng <= ng2; ng++) weight[ng] = 0.;
      }
   else
      {
        for(ng = ng1; ng <= ng2; ng++)
          {
/* Previous weight divided by this sigma: 
* I have tried without dividing: gives worse results 
* I have tried with **1, **0.8, **0.2 instead: gives worse results 
* JLP93
*/
           weight[ng] /= sqrt((double)sigm[nb]);
/* Normalisation to the number of bispectral terms:
* I have tried without dividing, and it gives better results 
* JLP93               weight[ng] /= (float)(ng2-ng1+1);
*/
          }
/* End of ng loop */
     }
/* End of (sigm[nb] >= sig_max || bmask[nb] == 0) */
  ng1 = ng2 +1;
  }
/* End of nb loop */
 
/* Close sigma file: */
fclose(fp2);

/* Normalisation of the weights: */
NORMALIZE_L1(weight,ngamma,&xnorm);

/* Diagnostic: */
  if(irem > 0)
  {
   printf("MODIF_WEIGHTS_SIGM2/Warning: %d discarded BETA values because of modulus < %f or sigma > %f \n",
          irem,lower_modulus,sig_max); 
   fprintf(fp1,"MODIF_WEIGHTS_SIGM2/Warning: %d discarded BETA values because of modulus < %f or sigma > %f \n",
          irem,lower_modulus,sig_max); 
  }

return(0);
}
/*******************************************************************
* MODIF_WEIGHTS_SIGM1: to modify the weights according to the sigma found
* when computing the recursive solution
* SIGM1 is specially designed for unity weights
*
* INPUT:
* sig_max: Upper threshold of sigma 
*          (truncation used to reduce the range of weights)
*
* OUTPUT
* weight: weights used for the bispectral list
*******************************************************************/
int MODIF_WEIGHTS_SIGM1(float *sigm, INT4 nbeta, INT4 ngamma, float sig_max,
                        float lower_modulus, INT4 ifermax,float *bmask,
                        float *weight)
{
FILE *fp2;
INT4 ng1, ng2, ng, nb, irem;
double xnorm;
 
ng1 = 1;
irem = 0;
 
/* Opening sigma file: */
if((fp2 = fopen("sigma.dat","w")) == NULL)
   {
   printf("MODIF_WEIGHTS_SIGM1/Fatal error opening sigma file \"sigma.dat\"\n");
   fprintf(fp1,"MODIF_WEIGHTS_SIGM1/Fatal error opening sigma file \"sigma.dat\"\n");
   fclose(fp1); exit(-1);
   }
  
for(nb = 3; nb <= nbeta; nb++)
   {
   cover_ngt1(&ng2,nb);
   fprintf(fp2,"%d %f \n",nb,sigm[nb]);

/* Truncation to reduce the range of weights: */
   if(sigm[nb] >= sig_max || bmask[nb] == 0.)
      {
      irem++;
      bmask[nb] = 0.;
      for(ng = ng1; ng <= ng2; ng++) weight[ng] = 0.;
      }
   else
      {
        for(ng = ng1; ng <= ng2; ng++)
          {
/* JLP93
* (Unity weights: gives worse results when not dividing by SIGM)
*/
           weight[ng] /= sigm[nb];
/*
* Normalisation to the number of bispectral terms:
* (Unity weights: gives worse results when dividing by nber_terms**2 only)
*/
           weight[ng] /= (float)(ng2-ng1+1)*(ng2-ng1+1);
          }
/* End of ng loop */
     }
/* End of (sigm[nb] >= sig_max || bmask[nb] == 0) */
  ng1 = ng2 +1;
  }
/* End of nb loop */
 
/* Close sigma file: */
fclose(fp2);


/* Normalisation of the weights: */
NORMALIZE_L1(weight,ngamma,&xnorm);

/* Diagnostic: */
  if(irem > 0)
  {
   printf("MODIF_WEIGHTS_SIGM1/Warning: %d discarded BETA values because of modulus < %f or sigma > %f \n",
          irem,lower_modulus,sig_max); 
   fprintf(fp1,"MODIF_WEIGHTS_SIGM1/Warning: %d discarded BETA values because of modulus < %f or sigma > %f \n",
          irem,lower_modulus,sig_max); 
  }

return(0);
}
/**********************************************************************
* Output SNR map according to the adopted weights
* Generating SNR and sigma maps (needed by DIANE)
**********************************************************************/
int OUTPUT_SNR1(float *sigm, INT4 nbeta, char *fname)
{
char name[60],comments[81];
INT4 nb, ix, iy, ixc, iyc, iix, iiy;
register int i, j;
float sigg;
 
/* Assume that zero frequency is at IXC,IYC:  */
   ixc = nx / 2;
   iyc = ny / 2;
 
/* Initializing RE and IM:
* 1/SIGM will be in RE and SIGM in IM:
*/
  for(j = 0; j < ny; j++)
    {
    for(i = 0; i < nx; i++)
        {
        re[i + j * nx] = 0.;
        im[i + j * nx] = 0.;
        }
    }
 
/* Generating SNR and sigma maps (needed for deconvolution by DIANE):
*/
   for(nb = 0; nb <= nbeta; nb++)
      {
/* Get coordinates (iix,iiy) from nb index: */
      COVER_IXY(&ix,&iy,&nb);
/* Add (ixc,iyc) vector, since centered FFT... */
      iix = ix + ixc;
      iiy = iy + iyc;
      sigg = sigm[nb];
/* SNR: */ 
      re[iix + iiy * nx] = 1. / sigg;
/* SIGMA: */
      im[iix + iiy * nx] = sigg * ro[nb];
/* Same values for the symmetrical frequency relative to the zero frequency: */
      iix = ixc -ix;
      iiy = iyc - iy;
      re[iix + iiy * nx] = 1. / sigg;
      im[iix + iiy * nx] = sigg * ro[nb];
      }
 
/* Output of SNR image (i.e., 1/sigm[nb]): 
*/
    sprintf(name,"%s.SNR",out_prefix);
    sprintf(comments," 1/SIGMA_PHASE of : %.20s",fname);
    printf(" Output of SNR=1/SIGMA_PHASE in image file: %s \n",name);
    JLP_WRITEIMAG(re,&nx,&ny,&nx,name,comments);

/* Output of SIGMA (i.e., sigm[nb] * ro[nb])
*/
    sprintf(name,"%s.SIG",out_prefix);
    sprintf(comments," MODULUS * SIGMA_PHASE of : %.20s",fname);
    printf(" Output of MODULUS * SIGMA_PHASE in image file: %s \n",name);
    JLP_WRITEIMAG(im,&nx,&ny,&nx,name,comments);

return(0);
}
/*******************************************************************
* ERROR_SIMU
* Compute the absolute errors for simulations
* (since model is known and stored in XCR)
*
* E is the mean of the errors of the phase terms 
*       weighted by the value of MODSQ**2:
*   or the mean of the errors of the full term (amplitude and phase)
* E = sqrt( Sum of (RO**2 * |XCR-XC|**2) / Sum of RO**2)
* whereas EP is raw:
* EP = sqrt( Sum of |XCR-XC|**2 / NBETA-2)
*******************************************************************/
int ERROR_SIMU(float *xc_re, float *xc_im, float *xcr_re, float *xcr_im,
               INT4 nbeta, float *err_spec, float *err_phas,
               char *errname, float *bmask)
{
FILE *fp2;
float cc_re, cc_im;
double r2, rr2, sr2;
INT4 indx, nb;
 
if((fp2 = fopen(errname,"w")) == NULL)
   {
   printf("ERROR_SIMU/Fatal error opening error file %s\n",errname);
   fprintf(fp1,"ERROR_SIMU/Fatal error opening error file %s\n",errname);
   fclose(fp1); exit(-1);
   }
  
fprintf(fp2," 0 0. 0. 0. Index, Phase error, Full quad. error, Cumul. mean phase error\n");

/* Global errors:
* err_spec : rms error of the spectrum
* err_phas: rms error of the phasor of the spectrum
*/
 *err_spec = 0.;
 *err_phas = 0.;
 sr2 = 0.;
 indx = 0;
 for(nb = 3; nb <= nbeta; nb++)
  {
   if(bmask[nb] != 0.)
     {
      indx++;
      r2 = ro[nb]*ro[nb];
      sr2 += r2;
      cc_re = xcr_re[nb] - xc_re[nb];
      cc_im = xcr_im[nb] - xc_im[nb];
      rr2 = cc_re * cc_re + cc_im * cc_im;
      *err_spec += r2 * rr2;
      *err_phas += rr2;
      fprintf(fp2,"%d %12.5e %12.5e %12.5e\n",nb,sqrt(rr2),sqrt(r2*rr2),
              sqrt((double)(*err_phas/(double)(nb-2))));
     }
  }
 
  if(sr2 == 0.)
    {
     *err_spec = 1000.;
     printf(" ERROR_SIMU/Warning: sr2 = 0. !! \n");
     fprintf(fp1," ERROR_SIMU/Warning: sr2 = 0. !! \n");
    }
/* err_spec = sqrt( Sum of (ro**2 * |xcr-xc|**2) / Sum of ro**2) */
   else
     *err_spec = sqrt((double)(*err_spec/sr2));

/* err_phas = sqrt( Sum of |xcr-xc|**2 / nbeta-2) */
     *err_phas = sqrt((double)(*err_phas/(double)indx));
 
fclose(fp2);
return(0);
}
/*******************************************************************
* EPLUS projects ALPHA_0 onto E^+
* xx(idim_x,4)
* vecxy(idim_xx,2) 
*******************************************************************/
int EPLUS(INT4 nbeta, float *bmask, float *xc_re, float *xc_im, 
          float *xx, INT4 idim_xx, float *vecxy)
{
float cc_re, cc_im, xxc, alpha;
INT4 nb;
 
/* Vectors of the base of the kernel */
  NOYAU(nbeta,bmask,vecxy,idim_xx);
 
/* Alpha_0 en xx(.,1) */
 for(nb = 0; nb <= nbeta; nb++)
   {
    cc_re = xc_re[nb];
    cc_im = xc_im[nb];
/* X(NB+1,1) = BMASK[NB]*ATAN2( IMAG(CC), REAL(CC) ) */
    jlp_atan2c(&xxc,cc_im,cc_re);
    xx[nb] = bmask[nb] * xxc;
   }
/* Alpha_0^+  en xx(.,1) */
    PROJ_EPLUS(nbeta,1,1,xx,idim_xx,vecxy);
 
/* Corresponding initial phasor: */
 for(nb = 0; nb <= nbeta; nb++)
   {
/* ALPHA=X(NB+1,1) */
    alpha = xx[nb];
    xc_re[nb] = bmask[nb] * cos((double)alpha); 
    xc_im[nb] = bmask[nb] * sin((double)alpha); 
   }
 
return(0);
}
/*******************************************************************
* NOYAU:
* Define the vectors of the base of the kernel
*  vecxy(idim_xx,2)
*******************************************************************/
int NOYAU(INT4 nbeta, float *bmask, float *vecxy, INT4 idim_xx)
{
INT4 nb, iix, iiy;
double xnorm, sum;
 
/* The base vectors are formed with the coordinates of the uv vectors. 
*/
for(nb = 0; nb <= nbeta; nb++)
  { 
/*  VECXY(NB,1)=BMASK(NB)*IXY(1,NB)
*  VECXY(NB,2)=BMASK(NB)*IXY(2,NB)
*/
    COVER_IXY(&iix,&iiy,&nb);
    vecxy[nb] = bmask[nb] * (float)iix;
    vecxy[nb + idim_xx] = bmask[nb] * (float)iiy;
  }

/* Normalization of the first vector: */
  NORMALIZE_L2(&vecxy[0],nbeta,&xnorm);

/* Orthogonalization (Schmidt process): */
/* First compute scalar product: */
sum = 0.;
for(nb = 0; nb <= nbeta; nb++)
   sum += (vecxy[nb] * vecxy[nb + idim_xx]);

/* Then removes to the second vector its projection unto the first one: */
/* VECXY(NB+1,2)=VECXY(NB+1,2)-SUM*VECXY(NB+1,1) */
for(nb = 0; nb <= nbeta; nb++)
   vecxy[nb + idim_xx] = vecxy[nb + idim_xx] - sum * vecxy[nb];

/* Normalization of the second vector: */
  NORMALIZE_L2(&vecxy[idim_xx],nbeta,&xnorm);
        
return(0);
}
/*******************************************************************
* PROJ_EPLUS
* Projects X(.,N1) onto E^+ and stores result in X(.,N2)
* X(IDIM_X,4),VECXY(IDIM_X,2)
*******************************************************************/
int PROJ_EPLUS(INT4 nbeta, INT4 n1, INT4 n2, float *xx,
               INT4 idim_xx, float *vecxy)
{
double xx1, yy1, delta;
INT4 nb, nn1, nn2;
 
/* Scalar product of X(.,N1) with the base vectors of the kernel */
xx1 = 0.;
yy1 = 0.;
nn1 = n1 - 1;
nn2 = n2 - 1;
 
for(nb = 1; nb <= nbeta; nb++)
  {
   xx1 += xx[nb + nn1 * idim_xx] * vecxy[nb];
   yy1 += xx[nb + nn1 * idim_xx] * vecxy[nb + idim_xx];
  }

 delta = sqrt(xx1 * xx1 + yy1 * yy1);
/* 
 printf("xx1 = %17.10e, yy1= %17.10e \n",xx1,yy1);
 printf("Distance to E+: DELTA = %17.10e \n",delta);
*/
 
/* Projection onto E^+ and stored in X(.,N2) */
for(nb = 1; nb <= nbeta; nb++)
   xx[nb + nn2 * idim_xx] = xx[nb + nn1 * idim_xx]
           - xx1 * vecxy[nb] - yy1 * vecxy[nb + idim_xx];

return(0);
}
/*******************************************************************
* SORTIE: write output files
* 
* INPUT:
* isortie: flag (if isortie=1 rei,imi, if isortie=2 RFT, IFT) 
*******************************************************************/
int SORTIE(float *xc_re, float *xc_im,
           INT4 nbeta, INT4 isortie, char *fname, float *bmask)
{
float cc_re, cc_im, xr, xi, xm;
char name[60], comments[81];
INT4 nb, ixc, iyc, iix, iiy, ix, iy; 
register int i, j;

/* Assume that zero frequency is at IXC,IYC:  */
   ixc = nx / 2;
   iyc = ny / 2;
 
/* Initializing the real and imaginary parts: */
  for(j = 0; j < ny; j++)
    {
    for(i = 0; i < nx; i++)
        {
           re[i + j * nx] = 0.;
           im[i + j * nx] = 0.;
        }
    }
 
/* Computing real and imaginary parts: */
for(nb = 0; nb <= nbeta; nb++)
   {
/* JLP94, try...
*  XM = BMASK[NB]*RO[NB]
*/
     xm = ro[nb];
     cc_re = xc_re[nb];
     cc_im = xc_im[nb];
     xr = xm * cc_re;
     xi = xm * cc_im;
     COVER_IXY(&ix,&iy,&nb);
     iix = ixc + ix;
     iiy = iyc + iy;
     re[iix + iiy * nx] = xr;
     im[iix + iiy * nx] = xi;
     iix = ixc - ix;
     iiy = iyc - iy;
     re[iix + iiy * nx] = xr;
     im[iix + iiy * nx] = -xi;
   }
 
/* Output to image files: 
* In rei and imi if isortie=1
* In ref and imf if isortie=1
*/ 
  if(isortie == 1)
    {
/* Initial output (recursive solution) */
    strcpy(name,"rei");
    sprintf(comments,"FFT (re) of: %.20s",fname);
    printf(" Output of %s \n",comments);
    JLP_WRITEIMAG(re,&nx,&ny,&nx,name,comments);

    strcpy(name,"imi");
    sprintf(comments,"FFT (im) of: %.20s",fname);
    printf(" Output of %s \n",comments);
    JLP_WRITEIMAG(im,&nx,&ny,&nx,name,comments);
    }
  else
    {
/* Final output (least squares solution) */
    sprintf(name,"%s.RFT",out_prefix);
    sprintf(comments,"FFT (re) of: %.20s",fname);
    printf(" Output of %s \n",comments);
    JLP_WRITEIMAG(re,&nx,&ny,&nx,name,comments);

    sprintf(name,"%s.IFT",out_prefix);
    sprintf(comments,"FFT (im) of: %.20s",fname);
    printf(" Output of %s \n",comments);
    JLP_WRITEIMAG(im,&nx,&ny,&nx,name,comments);
  }

return(0);
}
/**********************************************************************
* Main loop
* Least square minimization, non-linear least square fit
*
* INPUT:
* exit_tolerance = exit test for the main iteration (in degrees)
*                  This means that when the largest angular correction
*                  is smaller than "exit_tolerance", we exit from the main loop
*
* xx(.,4)
**********************************************************************/
int LSQUARES1(float exit_tolerance, INT4 nbeta, INT4 ngamma, INT4 ifermax,
              float *bmask, float *yce_re, float *yce_im,
              float *weight, float *xc_re, float *xc_im,
              float *xx, INT4 idim_xx, INT4 ittmax)
{
double xre, xim;
float w1, sup_phi, qmoy, qmoy0, delq, cos1, sin1, exit_tolerance_radians;
INT4 nb, itt, it;

exit_tolerance_radians = exit_tolerance * DEGRAD;
 
 qmoy0 = 0.;
 
/* Main loop: */
for(itt = 1; itt <= ittmax; itt++)
  {
/* Computes phase error term X(.,1) for PHI^+
* Iterative solution by conjugate gradients.
*/ 
   CGRADIENT(nbeta,ngamma,&it,&qmoy,ifermax,bmask,
             yce_re,yce_im,weight,xc_re,xc_im,xx,idim_xx);

   printf(" rms bisp error : %12.5e, itt = %d Internal it = %d done\n",
         qmoy,itt,it);
   fprintf(fp1," rms bisp error : %12.5e, itt = %d Internal it = %d done\n",
         qmoy,itt,it);
 
/* sup_phi : Maximum of the solution PHI^+ X(.,1) (L1 norm).
*/ 
   sup_phi = 0.;
   for(nb = 1; nb <= nbeta; nb++)
      {
      w1 = xx[nb] * bmask[nb];
      if(w1 < 0) w1 = -w1;
      if(w1 > sup_phi) sup_phi = w1;
      }
 
/* Convergence test:  */ 
   if (sup_phi < exit_tolerance_radians)
      {
       printf(" Normal exit: sup_phi < exit_tolerance \n"); 
       fprintf(fp1," Normal exit: sup_phi < exit_tolerance \n"); 
/* Exit from main loop: */
       break;
      }
 
/* Generating the new value of the spectral term XC(.) */ 
   for(nb = 1; nb <= nbeta; nb++)
      {
      w1 = xx[nb] * bmask[nb];
/* xc[nb]=xc[nb]*cmplx(cos(w1),sin(w1)) */
      cos1 = cos((double)w1);
      sin1 = sin((double)w1);
      COMPLEX_PROD2(xc_re[nb],xc_im[nb],cos1,sin1,
                    1, 1, &xre, &xim); 
      xc_re[nb] = xre;
      xc_im[nb] = xim;
      }
 
/* Test on the relative variation of qmoy 
*  between two successive iterations: */ 
   delq = (qmoy0 - qmoy) / qmoy;
   if(delq < 0) delq = -delq;
   if(delq < 2.E-4 && itt > 10)
     {
      printf(" Warning: Exit because the relative DELQ is too small\n");
      fprintf(fp1," Warning: Exit because the relative DELQ is too small\n");
/* Exit from main loop: */
     break;
    }
/* New value of qmoy0 for the next iteration: */
   else
    qmoy0 = qmoy;
/* End of main loop (itt index) */
  }

return(0);
}
/*******************************************************************
* CGRADIENT 
* Resolution of the linear system: 
*                  [ATwA] X(.,1)  = X(.,3)
* with conjugate gradients
* X(.,1) is the unknown (and then solution) PHI 
* X(.,2) is the direction D_n
* X(.,3) is the second member of the equation, and then the residual
* X(.,4) is ATwA D_n   (called Z_n)
*******************************************************************/
int CGRADIENT(INT4 nbeta, INT4 ngamma, INT4 *it, float *qmoy,
              INT4 ifermax, float *bmask, float *yce_re, float *yce_im,
              float *weitht, float *xc_re, float *xc_im,
              float *xx, INT4 idim_xx)
{
double ss, r2_rn, r2_rnplus1, r2_phiplus1;
float omega_n, work, gamma_n;
INT4 itmax, nb;

/* A good problem, with tolerance of 1.E-08, implies 4 iterations,
* Badly conditionned problem with cond_number=25 
* implies around 40 iterations, so:
*/
 itmax = 40;
 
/* STEP 0
* The starting solution PHI_0 is null : X(.,1)=0.
*/
/* Compute the initial residual R_0 and store it in X(.,3)
*/
  ATRANS(nbeta,ngamma,qmoy,ifermax,bmask,yce_re,yce_im,
         weight,xc_re,xc_im,xx,idim_xx);
 r2_rn = 0.;
 for(nb = 1; nb <= nbeta; nb++)
   {
/* The initial solution is set to 0. */
    xx[nb] = 0.;
/* The first direction D_0 = R_0 is copied to X(.,2)
* Note that X(NB+1,3) is null when BMASK[NB] is null.
*/
    xx[nb + idim_xx] = xx[nb + 2 * idim_xx];
/* Compute R2_RN : square of the norm of R_0
*  R2_RN is the square norm of R_N
*/
    r2_rn += xx[nb + 2 * idim_xx] * xx[nb + 2 * idim_xx];
  }
/* End of loop on nb */

  if(r2_rn < 1.e-15)
   {
   printf(" CGRADIENT/Error: Square norm of Residual_{N}=0 !");
   return(-1);
   }
 
/* Main loop: */
 for(*it = 1; *it <= itmax; (*it)++)
  { 
/* STEP 1
* Compute Z_N = [AT A] D_N and store it in X(.,4)
*/
   ATRANS_A(nbeta,ngamma,2,4,ifermax,bmask,weight,xx,idim_xx);

/* Compute OMEGA_N = R2_RN / (D_N scalar Z_N) */
/* Warning: SS is very small */
   ss = 0.;
   for(nb = 1; nb <= nbeta; nb++)
     {
     ss += xx[nb + idim_xx] * xx[nb + 3 * idim_xx] * bmask[nb];
     }
   omega_n = r2_rn / ss;

/* Compute the residual and the next value of PHI
*  R_{N+1} put to X(.,3) and  PHI_{N+1} put to X(.,1) ;
*  R2_RNPLUS1 is the square norm of R_{N+1} */
   r2_rnplus1 = 0.;
   r2_phiplus1 = 0.;
   for(nb = 1; nb <= nbeta; nb++)
     {
/* R_[N+1] = R_N - OMEGA_N * Z_N */
      xx[nb + 2 * idim_xx] -= (omega_n * xx[nb + 3 * idim_xx]);
/* PHI_[N+1] = PHI_N + OMEGA_N * D_N */
      xx[nb] += (omega_n * xx[nb + idim_xx]);
      r2_rnplus1 += (xx[nb + 2 * idim_xx] * xx[nb + 2 * idim_xx] * bmask[nb]);
      r2_phiplus1 += (xx[nb] * xx[nb] * bmask[nb]);
     }
 
/* STEP2 : Exit test */
   if(r2_phiplus1 < 1.e-15)
     {
     printf(" CGRADIENT/Error: Square norm of Delta_{PHI+1}=0 !\n");
     fprintf(fp1," CGRADIENT/Error: Square norm of Delta_{PHI+1}=0 !\n");
     return(-2);
     }

/* Normal exit: */
    work = r2_rnplus1 / r2_phiplus1;
    if(work < 1.e-5) return(0);
 
/* STEP 3
* GAMMA_N = || R_[N+1] ||**2 / || R_N ||**2
*/
    gamma_n = r2_rnplus1 / r2_rn;
/* For next step || R_[N+1] ||**2 becomes || R_N ||**2
*/
    r2_rn = r2_rnplus1;
/* Compute next conjugate direction D_{N+1} which is stored in X(.,1)
* JLP91    sum = 0;
*/
  for(nb = 1; nb <= nbeta; nb++)
    {
        xx[nb + idim_xx] = xx[nb + 2 * idim_xx] + gamma_n * xx[nb + idim_xx];
/* JLP91 Test if D_N and D_N+1 are orthogonal directions:
* JLP91    sum += xx[nb] * xx[nb + 2 * idim_xx];
*/
    }
/* JLP91 Test if ATwA D_N and D_N+1 are orthogonal directions:
* JLP91   printf(" D_n+1 scalar D_n is : %12.5e\n",sum);
*/

/* End of main loop
*/
   }
 
/* Abnormal exit: */
 printf(" Internal loop: not accurate enough, IT = %d, Test = %12.5e\n",
 itmax,work);
 fprintf(fp1," Internal loop: not accurate enough, IT = %d, Test = %12.5e\n",
 itmax,work);

return(-3);
}
/*******************************************************************
* ATRANS compute the second member of the system to be solved.
* This member is stored in X(.,3) which is the initial
* residual R_0
*
*  R_0 = AT ( PSI - A PHI_0)
* Here the residual is IMAG(exper.bispectrum * CONJ(estimated bispectrum))) :
*******************************************************************/
int ATRANS(INT4 nbeta, INT4 ngamma, float *qmoy, INT4 ifermax,
           float *bmask, float *yce_re, float *yce_im,
           float *weight, float *xc_re, float *xc_im,
           float *xx, INT4 idim_xx)
{
double c1_re, c1_im, c2_re, c2_im, c3_re, c3_im, yy1_re, yy1_im;
double q2;
INT4 ng1, ng2, nb, ng, indx, k, l, m;
 
/* Initialization to zero: */
for(nb = 1; nb <= nbeta; nb++)
   xx[nb + 2 * idim_xx] = 0.;
 
/* Weighted quadratic measure: */
   q2 = 0.;
 
 ng1 = 1;
 for(nb = 3; nb <= nbeta; nb++)
   {
/* NG2=NGT(NB) */
    cover_ngt1(&ng2,nb);
    indx = 0;
    for(ng = ng1; ng <= ng2; ng++)
      {
      indx++;
      if(indx <= ifermax && bmask[nb] != 0.)
        {
/*  K=KLM(1,NG) L=KLM(2,NG) M=KLM(3,NG)
*/
         cover_klm1(&k,1,ng);
         cover_klm1(&l,2,ng);
         cover_klm1(&m,3,ng);

/* Computed bispectrum (from the spectrum of the previous iteration):
*  c1=xc(k+1)*xc(l+1)*conjg(xc(m+1))
*/
         COMPLEX_PROD3(xc_re[k],xc_im[k],xc_re[l],xc_im[l],
                    xc_re[m],xc_im[m],1,1,-1,&c1_re,&c1_im);
/* Experimental bispectrum:
*/
         c2_re = yce_re[ng];
         c2_im = yce_im[ng];
/* Corresponding error (between experimental bispectrum and computed
* bispectrum from previous iteration):
*/
         c3_re = c2_re - c1_re;
         c3_im = c2_im - c1_im;
/* JLP93: Put WEIGHT=g**2 for AT A and for A (Use also g**2 for the
* estimation of convergence...)
*/
          q2 += (c3_re * c3_re + c3_im * c3_im) * SQUARE(weight[ng]);
/* Residual is IMAG(YCE[NG]*CONJ(C1)) :
* JLP93: Put WEIGHT=g**2 for AT A and for A: 
*  yy1=imag(c2*conjg(c1))*SQUARE(weight[ng])
*/
          COMPLEX_PROD2(c2_re,c2_im,c1_re,c1_im,1,-1,&yy1_re,&yy1_im);
          yy1_im *= SQUARE(weight[ng]);
/* Computing AT [YY1] by adding the contribution to X(.,3):
*/
          xx[k + 2 * idim_xx] += yy1_im;
          xx[l + 2 * idim_xx] += yy1_im;
          xx[m + 2 * idim_xx] -= yy1_im;
/* End of if(indx < ifermax): */
         }
/* End of loop on ng: */
      }
    ng1 = ng2 + 1;
/* End of loop on nb: */
   }
 
/* Mean quadratic error:
* It is a weighted error, and the sum of the weights is 1, so:
*/
   *qmoy = sqrt(q2);
/* 
   printf(" Mean weighted error of the solution on the bispectrum : %f\n",
           *qmoy);
   fprintf(fp1," Mean weighted error of the solution on the bispectrum : %f\n",
           *qmoy);
*/
 
return(0);
}
/*******************************************************************
* ATRANS_A compute X(.,N2) = [AT A] X(.,N1)
*******************************************************************/
int ATRANS_A(INT4 nbeta, INT4 ngamma, INT4 n1, INT4 n2, INT4 ifermax,
             float *bmask, float *weight, float *xx, INT4 idim_xx)
{
INT4 nn1, nn2, ng1, ng2, nb, ng, indx, k, l, m;
float yy1;
 
nn1 = n1 - 1;
nn2 = n2 - 1;

/* Initialization to zero: */
for(nb = 1; nb <= nbeta; nb++)
   xx[nb + nn2 * idim_xx] = 0.;
 
 ng1 = 1;
 for(nb = 3; nb <= nbeta; nb++)
   {
/* NG2=NGT(NB) */
    cover_ngt1(&ng2,nb);
    indx = 0;
    for(ng = ng1; ng <= ng2; ng++)
      {
      indx++;
      if(indx <= ifermax && bmask[nb] != 0.)
        {
/*  K=KLM(1,NG) L=KLM(2,NG) M=KLM(3,NG)
*/
         cover_klm1(&k,1,ng);
         cover_klm1(&l,2,ng);
         cover_klm1(&m,3,ng);

/* JLP93: Put WEIGHT=g**2 for AT A */
         yy1 = (xx[k + nn1 * idim_xx] + xx[l + nn1 * idim_xx]
                - xx[m + nn1 * idim_xx]) * SQUARE(weight[ng]);
/* Computing AT [YY1] by adding the contribution to X(.,N2): */
         xx[k + nn2 * idim_xx] += yy1;
         xx[l + nn2 * idim_xx] += yy1;
         xx[m + nn2 * idim_xx] -= yy1;
/* End of if(indx < ifermax): */
        }
/* End of loop on ng: */
      }
    ng1 = ng2 + 1;
/* End of loop on nb: */
  }
 
return(0);
}
/**********************************************************************
* FINAL_ERROR 
* Compute errors when object is available (simulations)
*
**********************************************************************/
int FINAL_ERROR(INT4 ir, INT4 nbeta, float *final_err_spec,
                float *final_err_phas, char *errname, float *bmask,
                float *xcr_re, float *xcr_im, float *xc_re, float *xc_im,
                float *xx, INT4 idim_xx, float *vecxy, float *sigm)
{
float cc_re, cc_im, w1, w2;
INT4 nb, irs, i, nbs;

ERROR_SIMU(xc_re,xc_im,xcr_re,xcr_im,nbeta,final_err_spec,
           final_err_phas,errname,bmask);
 
/* Ecart de phase angulaire point par point exprime en radian
*/ 
  for(nb = 0; nb <= nbeta; nb++)
    {
     cc_re = xcr_re[nb];
     cc_im = xcr_im[nb];
     jlp_atan2c(&w1,cc_im,cc_re);
     cc_re = xc_re[nb];
     cc_im = xc_im[nb];
     jlp_atan2c(&w2,cc_im,cc_re);
     xx[nb] = bmask[nb] * ( w1 - w2);
    }
 
/* Projection of this error X(.,1) onto E^+ */ 
  PROJ_EPLUS(nbeta,1,1,xx,idim_xx,vecxy);
 
/* Ecart angulaire correspondant aux u croissants en norme le long
* d'un rayon; finalement exprime en degre */ 
  for(irs = 1; irs <= ir; irs++)
  {
/* NBS = NBCOUV(IRS,0) */ 
   i = 0;
   COVER_NBCOUV(&nbs,&irs,&i,&ir);
   w1 = bmask[nbs] * xx[nbs] / DEGRAD;
   printf(" #%d Spect. index %d Phase err (deg) %8.3f Sigma %8.3f\n",
          irs,nbs,w1,sigm[nbs]);
   fprintf(fp1," #%d Spect. index %d Phase err (deg) %8.3f Sigma %8.3f\n",
          irs,nbs,w1,sigm[nbs]);
  }
 
return(0);
}
/**************************************************************/
#ifdef TTT
/**********************************************************************
* Output bispectrum errors:
* Warning: This prog does not use BMASK (---> check if this is right later!)
**********************************************************************/
        SUBROUTINE ERROR_BISPECT(NBETA,NGAMMA,ERRNAME,
     1                           ARRAY,NBX,NBY,IFERMAX,YCE,XC)
        INTEGER*4 IFERMAX,NBX,NBY
        REAL ARRAY(NBX,*)
        REAL*8 ERROR1
        CHARACTER ERRNAME*(*) 
        COMPLEX XC(*),YCE(*),C1,C2
        REAL MOD1,MOD2,RE1,IM1,ERROR0,ERROR2
 
/* Erase input array: */
        DO J=1,NBY
          DO I=1,NBX
            ARRAY(I,J)=0.
          END DO
        END DO

        WRITE(6,42) ERRNAME
        WRITE(2,42) ERRNAME
42      FORMAT(' Output of bispectrum errors (unweighted) in: ',A)
 
/* Opening output ASCII file:
*        OPEN(3,FILE=ERRNAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
*/

/*        WRITE(3,52) NG,ERROR1,ERROR2,RE1,IM1
*/
52      FORMAT('0 0. 0. 0. 0. Bisp.list, error,',
     1         ' cumul.error, error_re, error_im') 
 
        ERROR1=0
        NG1=1
        DO 3 NB=3,NBETA
/*  NG2=NGT(NB)
*/
          CALL COVER_NGT(NG2,NB)
          INDEX=0
            DO 2 NG=NG1,NG2
              INDEX=INDEX+1
                IF(INDEX.LE.IFERMAX)THEN
/*  K=KLM(1,NG)
*  L=KLM(2,NG)
*  M=KLM(3,NG)
*/
                  CALL COVER_KLM(K,1,NG)
                  CALL COVER_KLM(L,2,NG)
                  CALL COVER_KLM(M,3,NG)
/* Computed bispectrum (from the spectrum of the previous iteration):
* minus experimental bispectrum:
*/
                  C1=XC(K+1)*XC(L+1)*CONJG(XC(M+1))
                  C2=YCE[NG]
                  MOD1=REAL(C1)**2+IMAG(C1)**2
                  MOD2=REAL(C2)**2+IMAG(C2)**2
                  RE1=REAL(C1)/MOD1-REAL(C2)/MOD2
                  IM1=IMAG(C1)/MOD1-IMAG(C2)/MOD2
                  ERROR0=RE1*RE1+RE2*RE2
                  ARRAY(NB,NG-NG1+1)=SQRT(ERROR0)
                  ERROR1=ERROR1+ERROR0
                  IF(NG.GT.2)THEN
                    ERROR2 = SQRT(ERROR1/FLOAT(NG-2))
                  ELSE
                    ERROR2 = 0.
                  ENDIF
/*                  WRITE(3,53) NG,SQRT(ERROR0),ERROR2,RE1,IM1
*/
53                FORMAT(I5,1X,4(E12.4,1X))
                ENDIF
2            CONTINUE
           NG1=NG2+1
3       CONTINUE

        ERROR2 = SQRT(ERROR1/FLOAT(NG-2))
        WRITE(6,68) ERROR2
        WRITE(2,68) ERROR2
68      FORMAT(/,' Final bispectrum error (rms unweighted):',G10.4,/)

/*        CLOSE(3)
*/

        RETURN
        END
/*****************************************************************
* BISP_WEIGHT1
* Computing the weights:
* Old: Takes the modulus into account 
* New: Takes a simulated modulus into acount
*
* OUTPUT
* weight: weights used for the bispectral list
******************************************************************/
        SUBROUTINE BISP_WEIGHT1(NBETA,NGAMMA,IR,CTE,BMASK,RO,WEIGHT)
        REAL RO(*),BMASK(*)
        REAL WEIGHT(*),CTE,W1,WMIN,WMAX,RANGE_MAX1,RANGE_MAX2
        INTEGER NGAMMA,IR

/* 20 is too small, even 40 is a bit small...
*/
        RANGE_MAX1=800.
        RANGE_MAX2=100000.
        WRITE(2,23) RANGE_MAX1,RANGE_MAX2
        WRITE(6,23) RANGE_MAX1,RANGE_MAX2
23      FORMAT(' Weight3/Modulus, R1*R2*R3/sum*radius**2',/,
     1  ' Maximum range per column, and global: ',2(1X,F8.2))

/* First get the relative weights according to the modulus: 
*/
        WMAX=0.
        WMIN=RO(0+1)*RO(0+1)*RO(0+1)
        DO 203 NG=1,NGAMMA
/*  WORK=BMASK(KLM(1,NG))*BMASK(KLM(2,NG))*BMASK(KLM(3,NG))
*/
          CALL COVER_KLM(K,1,NG)
          CALL COVER_KLM(L,2,NG)
          CALL COVER_KLM(M,3,NG)
          WORK=BMASK(K+1)*BMASK(L+1)*BMASK(M+1)
          WEIGHT[NG]=WORK*RO(K+1)*RO(L+1)*RO(M+1)
          IF(WEIGHT[NG].NE.0.)WMIN=MIN(WMIN,WEIGHT[NG])
          WMAX=MAX(WMAX,WEIGHT[NG])
/*          WEIGHT[NG]=SQRT(WEIGHT[NG])
*         WEIGHT[NG]=MIN(RO(K+1),RO(L+1))
*         WEIGHT[NG]=MIN(WEIGHT[NG],RO(M+1))
*/
203     CONTINUE 

        WRANGE=WMAX/WMIN
        WRITE(2,24) WRANGE
        WRITE(6,24) WRANGE
24      FORMAT(' Initial range of the weights:',G12.5)

/* Now modulation according to the column (ie NB)
*/
        WMAX=0.
        WMIN=RO(0+1)*RO(0+1)*RO(0+1)
        IG1=1
        DO NB=3,NBETA
/*  IG2=NGT(NB)
*/
          CALL COVER_NGT(IG2,NB)

/* Truncation in the column to regularize the problem:
*/
          WWMIN=RO(0+1)*RO(0+1)*RO(0+1)
          WWMAX=0.
            DO NG=IG1,IG2
              IF(WEIGHT[NG].NE.0.)WWMIN=MIN(WWMIN,WEIGHT[NG])
              WWMAX=MAX(WWMAX,WEIGHT[NG])
            END DO

          IF(WWMAX.NE.0.)THEN

          WWMAX=WWMIN*RANGE_MAX1
          SUM0=0.
            DO NG=IG1,IG2
/* To check (JLP91)
*/
              IF(WEIGHT[NG].GT.WWMAX)THEN
               WRITE(2,69) NB,WEIGHT[NG],WWMAX
               WRITE(6,69) NB,WEIGHT[NG],WWMAX
69             FORMAT(' Column troncation, NB, wold,new',I5,2(1X,1PG12.5))
              ENDIF
              WEIGHT[NG]=MIN(WEIGHT[NG],WWMAX) 
              SUM0=SUM0+WEIGHT[NG]
            END DO


/* Normalizes the weights to obtain a 1/||frequency||**2 law 
* for the spectral list: 
* **4 is good, but weight range is too large for conditionning...
*  XW=4 + IXY(1,NB)**2+IXY(2,NB)**2
*/
          CALL COVER_IXY(IXY1,IXY2,NB)
          XW=4 + IXY1**2 + IXY2**2
          XW=SUM0*XW*XW
            DO NG=IG1,IG2
              WEIGHT[NG]=WEIGHT[NG]/XW
              IF(WEIGHT[NG].NE.0.)WMIN=MIN(WMIN,WEIGHT[NG])
              WMAX=MAX(WMAX,WEIGHT[NG])
            ENDDO

          ENDIF

          IG1=IG2+1
        ENDDO

        WRANGE=WMAX/WMIN
        WRITE(2,25) WRANGE
        WRITE(6,25) WRANGE
25      FORMAT(' Range of the weights before troncation: ',G12.5)
        PRINT *,' WMAX,WMIN',WMAX,WMIN

/* Reajustment of the range:
*/
        WMAX=RANGE_MAX2*WMIN
        PRINT *,' WMAX',WMAX
        DO NG=1,NGAMMA
         IF(WEIGHT[NG].GT.WMAX)PRINT *,' NG,WEIGHT[NG]',WEIGHT[NG],NG
         WEIGHT[NG]=MIN(WEIGHT[NG],WMAX)
        ENDDO

        RETURN
        END
#endif
/**********************************************************
* Largest eigen value: method of the power.
**********************************************************/
int EIGEN_VALUE1(INT4 nbeta, INT4 ngamma, float *xlambda1, INT4 ifermax,
                 float *bmask, float *weight, float *xx, INT4 idim_xx, 
                 float *vecxy)
{
double xnorm;
INT4 nb, i;

/* Initialization: */
 for(nb = 1; nb <= nbeta; nb++) xx[nb + 3 * idim_xx] = bmask[nb];

/* Method of the power, with 20 iterations: */
 for(i = 1; i <= 20; i++)
  {
/* Projection of this error onto E^+ */
   PROJ_EPLUS(nbeta,4,2,xx,idim_xx,vecxy);

/* Normalization:  */
   NORMALIZE_L2(&xx[1 + idim_xx],nbeta,&xnorm);

/*
    printf(" Eigen_value1/ Iteration %d Norm: %12.5e \n",i-1,xnorm);
*/

/* Computing ATwA of X(.,2). 
* Output in X(.,4)
*/
    ATRANS_A(nbeta,ngamma,2,4,ifermax,bmask,weight,xx,idim_xx);

/* End of loop on "i": */
  }

  *xlambda1 = xnorm;

return(0);
}
/**********************************************************
* Smallest eigen value:
* Method of the power, applied to I - AtwA/lambda1
**********************************************************/
int EIGEN_VALUE2(INT4 nbeta, INT4 ngamma, float xlambda1, float *xlambda2, 
                 INT4 ifermax, float *bmask, float *weight, float *xx, 
                 INT4 idim_xx, float *vecxy)
{
double xnorm;
INT4 nb, i;

/* Initialization: */
 for(nb = 1; nb <= nbeta; nb++) xx[nb + idim_xx] = bmask[nb];

/* Method of the power, with 20 iterations: */
 for(i = 1; i <= 20; i++)
  {
/* Projection of this error onto E^+
*        xx(1 + idim_xx] = 0.;
*        xx(2 + idim_xx] = 0.;
*/
   PROJ_EPLUS(nbeta,2,2,xx,idim_xx,vecxy);

/* Normalization:  */
   NORMALIZE_L2(&xx[1 + idim_xx],nbeta,&xnorm);

/*
    printf(" Eigen_value2/ Iteration %d Norm: %12.5e \n",i-1,xnorm);
*/

/* Computing ATwA of X(.,2). 
* Output in X(.,4)
*/
    ATRANS_A(nbeta,ngamma,2,4,ifermax,bmask,weight,xx,idim_xx);

/* Now computing I - ATwA/lambda1 :
*/
   for(nb = 1; nb <= nbeta; nb++)
        xx[nb + idim_xx] -= xx[nb + 3 * idim_xx]/xlambda1;

/* End of loop on "i": */
  }

/* mu = 1 - lambda2/lambda1
* Thus lambda2=lambda1*(1-mu)
*/
  *xlambda2 = xlambda1 * (1. - xnorm);

return(0);
}
/************************************************************
* Set complex array XC to (1.,0.)
*
************************************************************/
int ZERO_PHASE(float *xc_re, float *xc_im, INT4 nbeta)
{
INT4 nb;

for(nb = 0; nb <= nbeta; nb++)
  {
   xc_re[nb] = 1.;
   xc_im[nb] = 0.;
  }

return(0);
}
/***********************************************************
* jlp_atan2cc
* Arc Tangent (radians)
************************************************************/
int jlp_atan2cc(float *xxc, float cc_im, float cc_re)
{
/* Call jlp_atan2 in "jlp_atan2.for" 
* which calls FORTRAN routine ATAN2: */
/* jlp_atan2(xxc,&cc_im,&cc_re);
*/
return(0);
}
/***********************************************************
* jlp_atan2c
* C Version of Fortran ATAN2
* Arc Tangent (radians)
* Slightly different (leads to differences of 1/10000 in final images)
************************************************************/
int jlp_atan2c(float *xxc, float cc_im, float cc_re)
{
*xxc = atan((double)(cc_im / cc_re));
/* To be compatible with fortran ATAN2(cc_im,cc_re): */
 if(cc_re < 0) 
    {
    *xxc += PI;
    if(*xxc > PI) *xxc -= 2. * PI;
    }
return(0);
}
/**********************************************************************/
