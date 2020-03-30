/****************************************************************
* inv_bispec2_1D.c
* Same as inv_bispec2 but for 1D only.
* To compute the phase of the spectrum from that of the bispectrum,
*
* Internal format: FFT along the lines (for faster processing), 
* but in input/output: FFT along the columns. 
*
* Use C routines for uv coverage (jlp_cover_mask.c)
* Warning: This prog does not use BMASK (---> check if this is right later!)
*
* Contains CREYCE_SIMU_1D, CREYCE1_1D, CREYCE2_1D, TRANS_1D, RECURSIVE_1D, 
*          EPLUS_1D, WEIGHTS_SIGM1, WEIGHTS_SIGM2,
*          NOYAU_1D, PROJ_EPLUS_1D, ATRANS_1D, ATRANS_A_1D, CGRADIENT_1D, 
*          ERROR_SIMU_1D, SORTIE_1D, BISP_WEIGHT_1D, BISP_WEIGHT2_1D
*
* Syntax:
*
* Example1:
* RUNS INV_BISPEC2_1D 0 12,20,0.08,220,0.004,0.8 real_fft ima_fft 
*
*   with 0=simulation option
*   12=radius of uv cover., 20 = weights or phase error, 0.08=exit_tolerance,
*   220=nbeta, 0.004=low_modulus, 0.8= max sigma
*   real_fft and ima_fft= FFT of image to be reconstructed
*
* Example2:
* RUNS INV_BISPEC2_1D 1 12,20,0.08,220,0.004,0.8 modsq bisp1
*
*   with 1=real data, 
*   12=radius of uv cover., 20 = weights or phase error, 0.08=exit_tolerance,
*   220=nbeta, 0.004=low_modulus, 0.8 = max sigma
*   modsq=mean squared modulus of FFT
*   bisp1=mean bispectrum
*
* Example3:
* RUNS INV_BISPEC2_1D 2 12,20,0.08,220,0.004,0.8 modsq real_fft ima_fft bisp1
*
*   with 2=simulated data, 
*   12=radius of uv cover., 20 = weights or phase error, 0.08=exit_tolerance,
*   220=nbeta, 0.004=low_modulus, 0.8= max sigma
*   modsq=mean squared modulus of FFT (simulated) 
*   real_fft and ima_fft= FFT of image to be reconstructed
*   bisp1=mean bispectrum (simulated)
*
* JLP
* Version: 20-01-00
* From inv_bispec2.c (of september 1995)
*****************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jlp_ftoc.h>

/*
#define PI 3.141592654
*/
/* Don't forget the parenthesis... */
#define DEGRAD (PI/180.)
/* Minimum value of sigma (for the central frequencies): */
#define MIN_SIGMA 5.E-2

#define DEBUG
/* To debug this program: (reference phase is available for simulations only)*/
#define USE_REFERENCE_PHASE 

/* For ir = 25, nbmax = 980  and ngmax = 187566   FOR 2D!!!!!
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

/* For ir=30 :  (wrong since for 2D!!!)
#define IRMAX 30
#define MIRMAX -IRMAX
#define NBMAX 1410
#define NGMAX 388400
#define IDIM 256
*/
 
/* For ir=50 : (to be checked): */
#define IRMAX 50
#define MIRMAX -IRMAX
#define NBMAX 1410
#define NGMAX 388400
#define IDIM 256
 
int CREYCE_SIMU_1D(float **yce_re, float **yce_im, float **xc_re,
                   float **xc_im, float **xcref_re, 
                   float **xcref_im, INT4 ir, float cte, 
                   INT4 nbeta, INT4 ngamma, 
                   float lower_modulus, char *real_name, char *imag_name,
                   float **weight, float **bmask);
int CREYCE1_1D(float **yce_re, float **yce_im, float **xc_re, 
               float **xc_im, float **xcref_re,
               float **xcref_im, INT4 ir, float cte, 
               INT4 nbeta, INT4 ngamma, INT4 ngamma_max, 
               float lower_modulus, char *modsq_name, char *bisp_name, 
               float **weight, float **bmask);
int CREYCE2_1D(float **yce_re, float **yce_im, float **xc_re,
            float **xc_im, float **xcref_re,
            float **xcref_im, INT4 ir, float cte, 
            INT4 nbeta, INT4 ngamma, INT4 ngamma_max,
            float lower_modulus, char *modsq_name, 
            char *real_name, char *imag_name, char *bisp_name, 
            float **weight, float **bmask);
int RECURSIVE_1D(float *sigm, INT4 nbeta, INT4 ngamma, float sigma_null, 
                 INT4 ifermax, float *bmask, float *yce_re, float *yce_im, 
                 float *xc_re, float *xc_im, float *xcref_re, float *xcref_im, 
                 float *weight);
int BISP_WEIGHT11(INT4 nbeta, INT4 ngamma, INT4 ir, float cte,
                  float *bmask, float *ro, float *weight);
int BISP_WEIGHT0(INT4 ngamma, INT4 ir, float cte, float *weight);
int BISP_WEIGHT2(float *bisp1, INT4 nbeta, INT4 ngamma, INT4 ir,
                 float cte, float *weight);
int BISP_WEIGHT22(float *bisp1, INT4 nbeta, INT4 ngamma, INT4 ir,
                  float cte, float *bmask, float *weight);
int TRANS_1D(INT4 ir, INT4 nbeta, float *bmask, float *xc_re, float *xc_im,
             int zero_at_one);
int CREYCE_LOAD_BISP_1D(float **yce_re, float **yce_im, float *bisp1, 
                        INT4 ngamma, INT4 dim_bispec_input);
int COMPLEX_PROD2(float c1_re, float c1_im, float c2_re, float c2_im, 
                  int sign1, int sign2, float *c12_re, float *c12_im);
int COMPLEX_PROD3(float c1_re, float c1_im, float c2_re, float c2_im,
                  float c3_re, float c3_im, int sign1, int sign2, int sign3, 
                  float *c123_re, float *c123_im);
int OUTPUT_SNR1_1D(float *sigm, INT4 nbeta, char *fname);
int ZERO_PHASE_1D(float *xc_re, float *xc_im, INT4 nbeta);
int ERROR_SIMU_1D(float *xc_re, float *xc_im, float *xcref_re, float *xcref_im, 
                  INT4 nbeta, float *err_spec, float *err_phas,
                  char *errname, float *bmask);
int NORMALIZE_L1(float *array, INT4 npts, float *norm, int iy);
int NORMALIZE_L2(float *array, INT4 npts, float *norm, int iy);
int ATRANS_1D(INT4 nbeta, INT4 ngamma, float *qmoy, INT4 ifermax, 
              float *bmask, float *yce_re, float *yce_im,
              float *xc_re, float *xc_im,
              float *weight, float *xx, INT4 idim_xx);
int CGRADIENT_1D(INT4 nbeta, INT4 ngamma, INT4 *it, float *qmoy, 
                 INT4 ifermax, float *bmask, float *yce_re, float *yce_im, 
                 float *xc_re, float *xc_im,
                 float *weight, float *xx, INT4 idim_xx);
int LSQUARES1_1D(float exit_tolerance, INT4 nbeta, INT4 ngamma, INT4 ifermax,
                 float *bmask, float *yce_re, float *yce_im,
                 float *xc_re, float *xc_im,
                 float *weight, float *xx, INT4 idim_xx, INT4 ittmax);
int SORTIE_1D(float *xc_re, float *xc_im,
              INT4 nbeta, INT4 isortie, char *fname, float *bmask);
int PROJ_EPLUS_1D(INT4 nbeta, INT4 n1, INT4 n2, float *xx, 
                  INT4 idim_xx, float *vecxy);
int NOYAU_1D(INT4 nbeta, float *bmask, float *vecxy, INT4 idim_xx, int iy);
int EPLUS_1D(float *xc_re, float *xc_im, INT4 nbeta, float *bmask, 
             float *xx, INT4 idim_xx, float *vecxy);
int FINAL_ERROR_1D(INT4 ir, INT4 nbeta, float *final_err_spec, 
                   float *final_err_phas, char *errname, float *bmask,
                   float *xc_re, float *xc_im, float *xcref_re, float *xcref_im,
                   float *xx, INT4 idim_xx, float *vecxy, float *sigm);
int ATRANS_A_1D(INT4 nbeta, INT4 ngamma, INT4 n1, INT4 n2, INT4 ifermax,
                float *bmask, float *weight, float *xx, INT4 idim_xx);
int WEIGHTS_SIGM2(float *sigm, INT4 nbeta, INT4 ngamma, float sig_max,
                 float lower_modulus, INT4 ifermax, float *bmask,
                 float *weight);
int WEIGHTS_SIGM1(float *sigm, INT4 nbeta, INT4 ngamma, float sig_max,
                  float lower_modulus, INT4 ifermax,float *bmask,
                  float *weight);
int jlp_atan2c(float *xxc, float cc_im, float cc_re);

/* Logfile: */
 static FILE *fp1;

/* NBCOUV( X DE -IRMAX A IRMAX,  Y DE 0 A IRMAX )
* IXY( 1 POUR X ; 2 POUR Y,  NB DE 0 A  NBMAX)
*/
   static INT4 nbcouv[2*IRMAX+1][IRMAX+1], ixy[2][NBMAX+1];
 
/* klm(1,.) == k ; klm(2,.) == l ; klm(3,.) == m ;
*/
   static INT4 klm[3*NGMAX], ngt[NBMAX];
/* xc: Spectrum (list): phase factor
* ro: Spectrum (list): modulus
* yce: Bispectrum (list)
* xcref: Reference spectrum (list), phase factor  (not known for real observ.)
* vecxy: Kernel vector(s)
*/
   static float *ro, *xx, *vecxy, *bad_line;
   static float *re, *im;
   static INT4 nx, ny, ifermax = 10000;
/* ifermax: simply to give a limit to the number of closure relations to
* be taken into account (when this number would be really too big...)
* (Max. number of closure phase relations allowed for each pixel of
* the uv-coverage).
* IFERMAX is an internal upper limit, whereas NCLOSURE_MAX is set 
* when computing the data.
*/
main(argc, argv)
int argc;
char **argv;
{
 INT4 idim_xx, ir, nbeta, nbeta_max, istat, dim_spec;
 INT4 nx_mask, nclosure_max, ngamma, ngamma_max, ittmax;
 float *yce_re, *yce_im, *mask, *sigm, sig_max, *weight, *bmask;
 float *xcref_re, *xcref_im, *xc_re, *xc_im;
 float lower_modulus, cte, exit_tolerance, sigma_null;
 float ini_err_spec, ini_err_phas, final_err_spec, final_err_phas;
 float xlambda1, xlambda2;
 char fname[60], comments[80], date[40], errname[40], logfile_name[40]; 
 char modsq_name[60], real_name[60], imag_name[60], bisp_name[60], *pc;
 INT4 iopt, isize;
 register int i;
/* COMPLEX replaced by [2] array: */
 float cc[2];

  JLP_BEGIN();

/* Opening logfile: */
   strcpy(logfile_name,"inv_bispec2_1D.log");
   if((fp1 = fopen(logfile_name,"w")) == NULL)
      {
       printf("Fatal error opening logfile %s \n",logfile_name);
       exit(-1);
      }

  JLP_CTIME(date,&istat);
  fprintf(fp1," Program inv_bispec2_1D, version 23-01-00, date=%s \n",date);
  printf(" Program inv_bispec2_1D, version 23-01-00, date=%s \n",date);

 
/* One or three parameters only are allowed to run the program: */
/* Carefull: 7 parameters always, using JLP "runs" */
#ifdef DEBUG
  printf(" argc=%d\n",argc);
  printf(" argv[3]=>%s<\n",argv[3]);
#endif

/* debug: 
  printf("     argv[4]=%s argv[5]=%s \n", argv[4], argv[5]);
*/

/* With "runs" argc == 7, otherwise argc should be 1 or 5: */
if((argc < 5 && argc != 1) || (argc > 5 && argc !=7) 
/* If equal to 7, argv[4] must be filled: */
   || (argc == 7 && !strcmp(argv[4],"")))
  {
  printf("     Fatal error: Wrong syntax, argc=%d argv[4] = %s\n",argc,argv[4]);
  printf(" Syntax is:  \n");
  printf(" runs inv_bispec2_1D ioption parameters modsq bisp1 \n");
  printf(" Example: runs inv_bispec2_1D 1 12,20,0.08,220,0.004,0.8 modsq bisp1\n");
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
* exit_tolerance = exit test for the main iteration (in degrees)
* This means that when the largest angular correction
* is larger than "exit_tolerance", we exit from the main loop:
* Formerly:   exit_tolerance = 0.1
*/ 
  printf(" Radius (IR) of uv-coverage, phase error");
  printf("(0 if all weights=1, <0 if modulus weighted), exit_tolerance (degrees),");
  printf(" nbeta, lower_modulus, max sigma, sigma_value_when_modulus_is_null,");
  printf(" nclosure_max (used for input data computation)");
  printf("  (Ex: 12,20,0.1,220,0.,0.8,0.7,1000)");
  scanf("%d,%f,%f,%d,%f,%f,%f,%d",&ir,&cte,&exit_tolerance,
          &nbeta,&lower_modulus,&sig_max,&sigma_null,&nclosure_max);
/* Input file names: */
  switch(iopt)
   {
/* Simulations: */
   case 0:
     printf(" Input real part:= "); scanf("%s",real_name);
     printf(" Input imaginary part:= "); scanf("%s",imag_name);
     break;
/* Simulations with all files: */
   case 2:
     printf(" Input square modulus file:= "); scanf("%s",modsq_name);
     printf(" Input real part:= "); scanf("%s",real_name);
     printf(" Input imaginary part:= "); scanf("%s",imag_name);
     printf(" Input bispectrum file:= "); scanf("%s",bisp_name);
     break;
/* Default is observations: */
   default:
     printf(" Input square modulus file:= "); scanf("%s",modsq_name);
     printf(" Input bispectrum file:= "); scanf("%s",bisp_name);
     break;
   }
 }
else
 {
  sscanf(argv[1],"%d",&iopt);
  sscanf(argv[2],"%d,%f,%f,%d,%f,%f,%f,%d",&ir,&cte,&exit_tolerance,
          &nbeta,&lower_modulus,&sig_max,&sigma_null,&nclosure_max);
/* Input file names: */
  switch(iopt)
   {
/* Simulations: */
   case 0:
     strcpy(real_name,argv[3]);
     strcpy(imag_name,argv[4]);
     break;
/* Simulations with all files: */
   case 2:
     strcpy(modsq_name,argv[3]);
     strcpy(real_name,argv[4]);
     strcpy(imag_name,argv[5]);
     strcpy(bisp_name,argv[6]);
     break;
   default:
     strcpy(modsq_name,argv[3]);
     strcpy(bisp_name,argv[4]);
     break;
   }

 }
/* Removing blanks at the end: */
  pc = bisp_name; while(*pc && *pc != ' ') pc++; *pc='\0';
  pc = real_name; while(*pc && *pc != ' ') pc++; *pc='\0';
  pc = imag_name; while(*pc && *pc != ' ') pc++; *pc='\0';
  pc = modsq_name; while(*pc && *pc != ' ') pc++; *pc='\0';

/* Output to logfile: */
  fprintf(fp1," Option = %d \n",iopt);
  fprintf(fp1," Maximum number of closure relations allowed by the program: ifermax=%d \n",
          ifermax);
  printf(" Option = %d \n",iopt);
  printf(" Maximum number of closure relations allowed by the program: ifermax=%d \n",
          ifermax);
  fprintf(fp1,"ir: %d, cte: %6.2f, exit_tolerance: %10.3e, lower_modulus: %12.5e \n",
          ir,cte,exit_tolerance,lower_modulus);
  fprintf(fp1," Sig_max: %10.3e  Sig_null: %10.3e \n",
          sig_max,sigma_null);

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
  nx_mask = 2 * ir + 1;
  isize = nx_mask * sizeof(float);
  mask = (float *)malloc(isize);
  for(i = 0; i < nx_mask; i++) mask[i] = 1.;
/* IFERMAX is an internal upper limit, whereas NCLOSURE_MAX is set 
* when computing the data:
*/
  if(nclosure_max < 0)
    {
    printf(" Fatal error: nclosure_max must be strictly positive! \n");
    fprintf(fp1," Fatal error: nclosure_max must be strictly positive! \n");
    fclose(fp1); exit(-1);
    }
  COVERA_MASK_1D(mask,&nx_mask,&ir,&nclosure_max,&nbeta_max,&ngamma_max);
  free(mask);
 
  if(nbeta > nbeta_max)
   {
    fprintf(fp1,"inv_bispec2_1D/Error: NBETA (wanted) = %d whereas NBETA_MAX = %d\n",
            nbeta,nbeta_max);
    fprintf(fp1," I correct it to BETA_MAX\n");
    printf("inv_bispec2_1D/Error: NBETA (wanted) = %d whereas NBETA_MAX = %d\n",
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
  ittmax = 15;
 
  fprintf(fp1,
   " Angular constant for bispectral noise estimation: %8.3e (degrees)\n",cte);
  fprintf(fp1," EXIT_TOLERANCE (for exit check) : %12.3E \n",exit_tolerance);
 
/* Allocation of memory
*/
   dim_spec = nbeta + 1;
   isize = dim_spec * sizeof(float);
/* REAL VECXY(NBMAX+1) */
   vecxy = (float *)malloc(isize);
/* REAL X(NBMAX+1,4) */
   isize = dim_spec * 4 * sizeof(float);
   xx = (float *)malloc(isize);
 
/* First dimension of array xx  */
   idim_xx = nbeta + 1;

/* Create YCE : bispectrum phasor, and initial set of weights */
/* Simulations (re,im): */
  if(iopt == 0)
     {
     CREYCE_SIMU_1D(&yce_re,&yce_im,&xc_re,&xc_im,&xcref_re,&xcref_im,
                    ir,cte,nbeta,ngamma,lower_modulus,
                    real_name,imag_name,&weight,&bmask);
     strcpy(fname,real_name);
     }
/* Real observations (modsq,bisp): */
  else if(iopt == 1)
     {
     CREYCE1_1D(&yce_re,&yce_im,&xc_re,&xc_im,&xcref_re,&xcref_im,
                ir,cte,nbeta,ngamma,ngamma_max,
                lower_modulus,modsq_name,bisp_name,&weight,&bmask);
     strcpy(fname,modsq_name);
     }
/* Simulations (re,im) and (modsq,bisp): */
  else
     {
     CREYCE2_1D(&yce_re,&yce_im,&xc_re,&xc_im,&xcref_re,&xcref_im,
                ir,cte,nbeta,ngamma,ngamma_max,
                lower_modulus,modsq_name,real_name,imag_name,bisp_name,
                &weight,&bmask);
     strcpy(fname,modsq_name);
     }

/******************************************************
* To initialize the solution,
* we solve the problem with the recursive method (Weigelt,...) 
* The initial spectral phase is stored in XC(.)
*/
  dim_spec = nbeta + 1;
  isize = dim_spec * ny * sizeof(float);
  sigm = (float *)malloc(isize);
  RECURSIVE_1D(sigm,nbeta,ngamma,sigma_null,ifermax,
               bmask,yce_re,yce_im,xc_re,xc_im,xcref_re,xcref_im,weight);
 
/* Other alternative (Null phases)
  ZERO_PHASE_1D(xc_re,xc_im,nbeta);
*/
 
/* New version of the weights:
* JLP93:
*/
  if(cte == 0.)
      WEIGHTS_SIGM1(sigm,nbeta,ngamma,sig_max,
                    lower_modulus,ifermax,bmask,weight);
  else
      WEIGHTS_SIGM2(sigm,nbeta,ngamma,sig_max,
                    lower_modulus,ifermax,bmask,weight);

/* Output of SNR map: */
  OUTPUT_SNR1_1D(sigm,nbeta,fname);

/* Output the errors of initial bispectrum: */
#ifdef JLP94 
  strcpy(errname,"bisp_error1.dat");
  nbx = nbeta;
  nby = nbeta/2;
  isize = nbx * nby * sizeof(float);
  work_array = (float *)malloc(isize);
  ERROR_BISPECT_1D(nbeta,ngamma,errname,work_array,nbx,nby,
                   ifermax,yce_re,yce_im);
  strcpy(errname,"bisperr1");
  strcpy(comments,"Quadratic errors");
  JLP_WRITEIMAG(work_array,&nbx,&nby,&nbx,errname,comments);
  free(work_array);
#endif

/* Initial error estimation: */
 if(iopt != 1)
   {
   strcpy(errname,"error1.dat");
   ERROR_SIMU_1D(xc_re,xc_im,xcref_re,xcref_im,nbeta,&ini_err_spec,
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
 
/* Translation of nx/2 to put the spectrum in the middle line
* Not necessary if zero phase input
*/
 if(iopt != 1)
   TRANS_1D(ir,nbeta,bmask,xc_re,xc_im,0);
 
/* SORTIE DES PARTIES REELLES ET IMAGINAIRES DE XC
* In rei and imi (i for initial)
*/
    SORTIE_1D(xc_re,xc_im,nbeta,1,fname,bmask);

/* JLP99 */
 if(ir > 0) {JLP_END(); exit(0);}

/* Projection of the initial solution onto E^+ */
    EPLUS_1D(xc_re,xc_im,nbeta,bmask,xx,idim_xx,vecxy);

/* Main loop: least square non linear fit. */
   exit_tolerance *= DEGRAD;
   LSQUARES1_1D(exit_tolerance,nbeta,ngamma,ifermax,
             bmask,yce_re,yce_im,xc_re,xc_im,weight,xx,idim_xx,ittmax);
 
/* Files ref et imf (f for final) */
   SORTIE_1D(xc_re,xc_im,nbeta,2,fname,bmask);

/* Calage en translation */
   TRANS_1D(ir,nbeta,bmask,xcref_re,xcref_im,1);

/* Final global errors and conclusions */
  if(iopt != 1)
    {
     strcpy(errname,"error2.dat");
     FINAL_ERROR_1D(ir,nbeta,&final_err_spec,&final_err_phas,errname,bmask,
                 xc_re,xc_im,xcref_re,xcref_im,xx,idim_xx,vecxy,sigm);
 
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
#ifdef JLP94 
     strcpy(errname,"bisp_error2.dat");
     nbx = nbeta;
     nby = nbeta/2;
     isize = nbx * nby * sizeof(float);
     work_array = (float *)malloc(isize);
     ERROR_BISPECT_1D(nbeta,ngamma,errname,work_array,nbx,nby,
                   ifermax,yce_re,yce_im);
     strcpy(errname,"bisperr2");
     strcpy(comments,"Quadratic errors");
     JLP_WRITEIMAG(work_array,&nbx,&nby,&nbx,errname,comments);
     free(work_array);
#endif

 fclose(fp1);
 printf(" Log File in \"inv_bispec2_1D.log\" \n");
 JLP_END();
 exit(0);
}
/*******************************************************************
* CREYCE2_1D
* Reads the modulus (squared) and bispectrum derived from simulations,
* re and im (available since it is a simulation)
*******************************************************************/
int CREYCE2_1D(float **yce_re, float **yce_im, float **xc_re,
            float **xc_im, float **xcref_re,
            float **xcref_im, INT4 ir, float cte, 
            INT4 nbeta, INT4 ngamma, INT4 ngamma_max,
            float lower_modulus, char *modsq_name, 
            char *real_name, char *imag_name, char *bisp_name, 
            float **weight, float **bmask)
{
FILE *fp2;
double errmod, errbisp, errbispw;
float xw0, xw1, xw2, w1, xr, xi, xm, yy1_re, yy1_im, work;
INT_PNTR pntr_ima;
INT4 nx1, ny1; 
INT4 isize, inull, iix, iy, iy_s, iy_b, ixc, irem;
/* nb and ng cannot be declared as "register" 
   since they are used as &nb or &ng in function calls... */
INT4 nb, ng, k, l, m, kk, dim_spec, dim_bispec, line_to_output;
float *re_rev, *im_rev, *modsq_rev;
float *modsq, *bisp1, xnorm, ro_max;
char comments[80];
register int i, j;
/*  COMPLEX YCE(*),YY1
*/
 
/* Input of the modulus: */
  printf(" Modulus (squared) of the FFT of the image along the columns\n");
/* Warning, input spectroscopic mode: line/columns are reversed: */
  JLP_VM_READIMAG1(&pntr_ima,&ny,&nx,modsq_name,comments);
  modsq_rev = (float *)pntr_ima;
  fprintf(fp1," Sq. modulus: %.14s  Comments: %.30s\n",modsq_name,comments);
   isize = nx * ny * sizeof(float);
   modsq = (float *)malloc(isize);
/* Warning, input spectroscopic mode: line/columns are reversed: */
  for(j = 0; j < ny; j++)
    for(i = 0; i < nx; i++)
     modsq[i + j * nx] = modsq_rev[j + i * ny];

/* Reading RE and IM : */
  printf(" Real part of the FFT of the image along the columns\n");
/* Warning, input spectroscopic mode: line/columns are reversed: */
  JLP_VM_READIMAG1(&pntr_ima,&ny,&nx,real_name,comments);
  re_rev = (float *)pntr_ima;
  fprintf(fp1," Real part of the FFT: %.14s  Comments: %.30s \n",
              real_name,comments);
  isize = nx * ny * sizeof(float);
  re = (float *) malloc(isize);
  for(j = 0; j < ny; j++)
    for(i = 0; i < nx; i++)
     re[i + j * nx] = re_rev[j + i * ny];

  printf(" Imaginary part of the FFT of the image along the columns\n");
  JLP_VM_READIMAG1(&pntr_ima,&ny,&nx,imag_name,comments);
/* Warning, input spectroscopic mode: line/columns are reversed: */
  im_rev = (float *)pntr_ima;
  fprintf(fp1," Imag. part of the FFT: %.14s  Comments: %.30s \n",
              imag_name,comments);
  isize = nx * ny * sizeof(float);
  im = (float *) malloc(isize);
  for(j = 0; j < ny; j++)
    for(i = 0; i < nx; i++)
     im[i + j * nx] = im_rev[j + i * ny];

/* Allocation of memory: */
   dim_spec = nbeta + 1;
   isize = dim_spec * ny * sizeof(float);
   *xc_re = (float *)malloc(isize);
   *xc_im = (float *)malloc(isize);
   *xcref_re = (float *)malloc(isize);
   *xcref_im = (float *)malloc(isize);
   ro = (float *)malloc(isize);
   *bmask = (float *)malloc(isize);
   bad_line = (float *)malloc(isize);

/* Opening error file for modulus: */
   if((fp2 = fopen("err_mod.dat","w")) == NULL)
      {
       printf("CREYCE2_1D/Fatal error opening error file for modulus \"err_mod.dat\"\n");
       fprintf(fp1,"CREYCE2_1D/Fatal error opening error file for modulus \"err_mod.dat\"\n");
       fclose(fp1); exit(-1);
      }
   fprintf(fp2,
" 0 0. 1. 1. (line#%d), modulus error, input mod., theoretical mod.\n",ny/2);

/* Assume that zero frequency is at IXC,IY:  */
   ixc = nx / 2;
/********************************************************************/
/****************** Loop on iy: *************************************/
 for(iy = 0; iy < ny; iy++)
    {
    iy_s = iy * dim_spec;
     if(iy == ny/2)
     {
     xw0 = re[ixc + iy * nx]; xw0 *= xw0;
     xw1 = im[ixc + iy * nx]; xw1 *= xw1;
     xw1 = sqrt((double)(xw0 + xw1));
     xw2 = sqrt( (double)modsq[ixc + iy * nx]);
     printf(" Central value of input and computed modulus: %12.5e %12.5e\n",
            xw2,xw1);
     fprintf(fp1," Central value of input and computed modulus: %12.5e %12.5e\n",
            xw2,xw1);
     }
 
/* Loading modulus and phasor to ro and xc */
     inull = 0;
     errmod = 0.;
     for(nb = 0; nb <= nbeta; nb++)
        {
/* Get coordinates (iix,0) from nb index: */
         COVER_IXY_1D(&iix,&nb);
/* Add (ixc,0) vector, since centered FFT... */
         iix += ixc;
         xr = re[iix + iy * nx];
         xi = im[iix + iy * nx];
         xm = sqrt((double)(xr * xr + xi * xi));
/* Phase factor: */
          if(xm == 0.)
            {
            (*xc_re)[nb + iy_s] = 1.;
            (*xc_im)[nb + iy_s] = 0.;
            inull++;
            }
          else
            {
            (*xc_re)[nb + iy_s] = xr / xm;
            (*xc_im)[nb + iy_s] = xi / xm;
            }

/* Modulus: */
          w1 = modsq[iix + iy * nx];
          if(w1 < 0.) w1 = 0.;
          ro[nb + iy_s] = sqrt((double)w1);
          if(w1 == 0.) inull++;
/* Write error of central line to file: */
          if(iy == ny/2)
            {
            work = (xm/xw1 - ro[nb + iy_s] / xw2);
            fprintf(fp2,"%d %12.5e %12.5e %12.5e\n",nb,work,ro[nb]/xw2,xm/xw1);
            work *= work;
            errmod += work;
            }

/* End of "for ... nb" loop */
    }

#ifdef DEBUG
    if(iy == ny/2)
    {
      errmod = sqrt(errmod/(float)nbeta);
      printf(" sqrror of the modulus : %17.10e \n",errmod);
      if(inull > 0)
        {
        printf(" CREYCE2_1D/Modulus null for %d values (iy=%d)\n",inull,iy);
        fprintf(fp1," CREYCE2_1D/Modulus null for %d values (iy=%d)\n",inull,iy);
        }

      printf(" rms error of the modulus : %11.4e \n",errmod);
      fprintf(fp1," rms error of the modulus : %11.4e \n",errmod);
 
/* Output of some errors of the modulus: */
        for(nb = 0; nb < 5; nb++)
           {
/* Get coordinates (iix,0) from nb index: */
           COVER_IXY_1D(&iix,&nb);
/* Add (ixc,0) vector, since centered FFT... */
           iix += ixc;
           xr = re[iix + iy * nx];
           xi = im[iix + iy * nx];
           xm = sqrt((double)(xr * xr + xi * xi))/xw1;
           printf(" Index=%d, normalized RO and computed modulus: %11.4e %11.4e \n",
           nb, ro[nb + iy_s]/xw2, xm);
           fprintf(fp1," Index=%d, normalized RO and computed modulus: %11.4e %11.4e \n",
           nb, ro[nb + iy_s]/xw2, xm);
           }
/* End of case iy==ny/2 */
      }
#endif
 
/* Truncation: */
        irem = 0;
/* Mask to discard some values of the spectral list: */
        ro_max = ro[0 + iy_s];
        if(ro_max <= 0.) ro_max = 1.;
        for(nb = 0; nb <= nbeta; nb++)
        {
          work = ro[nb + iy_s] / ro_max;
/* Check that RO[NB] greater than LOWER_MODULUS:
* Remember that LOWER_MODULUS can be negative... */
          if(work < lower_modulus || work <= 0.)
            {
             (*bmask)[nb + iy_s] = 0.;
             irem++;
            }
          else
            (*bmask)[nb + iy_s] = 1.;
        }
 
 if(irem)
   {
   printf("Line#%d/Number of terms of the spectral list after correction: %d\n",
            iy,nbeta-irem);
   fprintf(fp1,"Line#%d/Number of terms of the spectral list after correction:%d\n",
            iy,nbeta-irem);
   }

/* End of loop on iy: */
  }
/********************************************************************/

/* Spectrum phasor blocked in translation: xc[.] */
  TRANS_1D(ir,nbeta,*bmask,*xc_re,*xc_im,1);

/* Reference (since simulation) */
  for(iy = 0; iy < ny; iy++)
  {
  iy_s = iy * dim_spec;
  for(nb = 0; nb <= nbeta; nb++)
    {
    (*xcref_re)[nb + iy_s] = (*xc_re)[nb + iy_s];
    (*xcref_im)[nb + iy_s] = (*xc_im)[nb + iy_s];
    }
  }

/* Close error file: */
  fclose(fp2);

/***********************************************************************
* Bispectrum phasor: YCE
*/
  printf(" Bispectrum of the image (phase term)\n");
  JLP_VM_READIMAG1(&pntr_ima,&nx1,&ny1,bisp_name,comments);
  bisp1 = (float *)pntr_ima;
  fprintf(fp1," Bispectrum: %.14s  Comments: %.30s\n",bisp_name,comments);
  if( (ny1 < 1) || (nx1 != 3 * ngamma_max))
    {
    printf(" CREYCE2_1D/FATAL ERROR: inconsistent size of bispectrum: nx1=%d ngamma_max= %d\n",nx1,ngamma_max);
    fprintf(fp1, " CREYCE2_1D/FATAL ERROR: inconsistent size of bispectrum: nx1=%d ngamma_max= %d\n",nx1,ngamma_max);
    fclose(fp1); exit(-1);
    }

/* Load BISP array to YCE array: */
  CREYCE_LOAD_BISP_1D(yce_re,yce_im,bisp1,ngamma,3*ngamma);

/* REAL WEIGHT(NGMAX) */
  dim_bispec = ngamma + 1;
  isize = dim_bispec * ny * sizeof(float);
  *weight = (float *) malloc(isize);

/* Computing the weights: */
  if(cte > 0)
    {
    for(iy = 0; iy < ny; iy++)
      {
      iy_s = iy * dim_spec;
      iy_b = iy * dim_bispec;
      BISP_WEIGHT11(nbeta,ngamma,ir,cte,
                   &(*bmask)[iy_s],&ro[iy_s],&(*weight)[iy_b]);
      NORMALIZE_L1(&(*weight)[iy_b],ngamma,&xnorm,iy);
      bad_line[iy] = (xnorm > 0) ? 0 : 1;
      }
    }
  else if (cte < 0)
    {
    for(iy = 0; iy < ny; iy++)
      {
      iy_s = iy * dim_spec;
      iy_b = iy * dim_bispec;
/* Version with SNR stored in 3rd line of bispectrum: */
      BISP_WEIGHT2(&bisp1[iy_b],nbeta,ngamma,ir,cte,&(*weight)[iy_b]);
/*    BISP_WEIGHT1(nbeta,ngamma,ir,cte,(*bmask),weight);
*/
      NORMALIZE_L1(&(*weight)[iy_b],ngamma,&xnorm,iy);
      bad_line[iy] = (xnorm > 0) ? 0 : 1;
      }
    }
  else 
    {
     printf(" Weights set to unity, and then normalized \n");
     fprintf(fp1," Weights set to unity, and then normalized \n");
     for(iy = 0; iy < ny; iy++)
       {
       iy_b = iy * dim_bispec;
       iy_s = iy * dim_spec;
       for(ng = 0; ng < ngamma; ng++)
         {
         (*weight)[ng + iy_b]=1.;
         for(kk = 1; kk <= 3; kk++) 
            {
/*  K = KLM(KK,NG)*/
            cover_klm0(&k,kk,ng);
            if((*bmask)[k + iy_s] == 0.) 
                     (*weight)[ng + iy_b]=0.;
            }
         }
       NORMALIZE_L1(&(*weight)[iy_b],ngamma,&xnorm,iy);
       bad_line[iy] = (xnorm > 0) ? 0 : 1;
       }
    }
 
/* Compute the errors for line #0(since xcref is available): */
#ifdef DEBUG
  line_to_output = 0;
  iy = line_to_output;
  iy_s = iy * dim_spec; 
  iy_b = iy * dim_bispec; 
  errbisp = 0.;
  errbispw = 0.; 
  for(ng = 0; ng < ngamma; ng++)
    {
     cover_klm0(&k,1,ng);
     cover_klm0(&l,2,ng);
     cover_klm0(&m,3,ng);
/* Theoretical bispectrum: yy1 = xcref[k]*xcref[l]*conjg(xcref[m]);
*/
     COMPLEX_PROD3((*xcref_re)[k + iy_s],(*xcref_im)[k + iy_s],
                   (*xcref_re)[l + iy_s],(*xcref_im)[l + iy_s],
                   (*xcref_re)[m + iy_s],(*xcref_im)[m + iy_s],
                   1,1,-1,&yy1_re,&yy1_im);
#ifdef DEBUG_2
     if(ng < 10) printf("ng=%d comp. bisp=(%e,%e) input=(%e,%e)\n",
                        ng,yy1_re,yy1_im,
                        (*yce_re)[ng + iy_b],(*yce_im)[ng + iy_b]);
#endif

/* True error: */
     work = (yy1_re - (*yce_re)[ng + iy_b])*(yy1_re - (*yce_re)[ng + iy_b])
           + (yy1_im - (*yce_im)[ng + iy_b])*(yy1_im - (*yce_im)[ng + iy_b]);
     errbisp += work;
     errbispw += (work * (*weight)[ng + iy_b]);
     }

   errbisp = sqrt(errbisp/(double)ngamma);
   errbispw = sqrt(errbispw);
   printf("\n************ Line #%d  *************\n",iy);
   printf(" rms error of the bispectrum : %11.4e (iy=%d)\n",errbisp,iy);
   printf(" weighted rms error of the bispectrum : %11.4e \n",errbispw);
   fprintf(fp1," rms error of the bispectrum : %11.4e (iy=%d)\n",errbisp,iy);
   fprintf(fp1," weighted rms error of the bispectrum : %11.4e \n",errbispw);
 
/* Just debug mode: */
   for(ng = 0; ng < 5; ng++)
    {
/*  yy1=xcref(k)*xcref(l)*conjg(xcref(m)) */
     cover_klm0(&k,1,ng);
     cover_klm0(&l,2,ng);
     cover_klm0(&m,3,ng);
     COMPLEX_PROD3((*xcref_re)[k + iy_s],(*xcref_im)[k + iy_s],
                   (*xcref_re)[l + iy_s],(*xcref_im)[l + iy_s],
                   (*xcref_re)[m + iy_s],(*xcref_im)[m + iy_s],
                   1,1,-1,&yy1_re,&yy1_im);
     printf(" ng=%d k=%d l=%d m=%d \n xcref_re[k]=%f ",
             ng,k,l,m,(*xcref_re)[k + iy_s]);
     printf(" xcref_im[k]=%f xcref_re[l]=%f ",
             (*xcref_im)[k + iy_s],(*xcref_re)[l + iy_s]);
     printf(" xcref_im[l]=%f xcref_re[m]=%f xcref_im[m]=%f \n",
            (*xcref_im)[l + iy_s],(*xcref_re)[m + iy_s],(*xcref_im)[m + iy_s]);
     printf(" Input, computed bisp: %d (%8.5f,%8.5f) (%8.5f,%8.5f)\n",
            ng, (*yce_re)[ng + iy_b],(*yce_im)[ng + iy_b],yy1_re,yy1_im);
    }
   printf("******************************\n");
#endif
 
return(0);
}
/*******************************************************************
* TRANS_1D
* fait le calage en translation du facteur de phase spectral XC(.)
* En sortie, la phase de nb=1 est zero degre.
*
* A translation in the direct plane of Dx induces an increase of
* the phase of u_k by (2 pi u_k Dx):
* TF[f(x - Dx)](u_k) = exp(-2 i pi u_k Dx) * TF[f(x)](u_k)  
*
*  XC(1) = exp i Beta(1)
*
* Array TX is defined as:
*  TX(1) = exp -i Beta(1)
*  TX(-1)= exp i Beta(1)
*  TX(2) = exp -i 2 Beta(1)
*   ....
*  TX(N) = exp -i N Beta(1)
*
* OUTPUT:
*  If coordinates of frequency #NB is: IXY1
*  XC[NB] = XC[NB] * TX(IXY1)
*  i.e.:
*  exp i Beta_new = exp i Beta_old * exp -i IXY1 Beta(1) 
*******************************************************************/
int TRANS_1D(INT4 ir, INT4 nbeta, float *bmask, float *xc_re, float *xc_im,
             int zero_at_one)
{ 
float tx_re[IRMAX+1], tx_im[IRMAX+1];
float cx_re, cx_im, ccx_re, ccx_im;
float work_re, work_im;
double aang;
INT4 iix, nb;
int iy, dim_spec, iy_s;
register int i;
 
/* If a reference phase has been used for each line, return
* without doing anything. Otherwise it would un-do what has been done... */
#ifdef USE_REFERENCE_PHASE 
 if(zero_at_one)
 return(0);
#endif

dim_spec = nbeta + 1;

 if(ir > IRMAX)
   {
     printf(" TRANS_1D/Fatal error, maximum IR= %d \n",IRMAX);
     fprintf(fp1," TRANS_1D/Fatal error, maximum IR= %d \n",IRMAX);
     fclose(fp1); exit(-1);
   }

/* Main loop on all the lines: */
 for(iy = 0; iy < ny; iy++)
   {
   iy_s = iy * dim_spec;
/* Take the conjugate value of xc[1] as a starting point,
* to obtain a null phase  for xc[1] at the end of this routine
* (i.e., xc[1]=(1,0) ) 
*/ 
/* cx=conjg(xc(1)) */
if(zero_at_one)
   {
   cx_re = xc_re[1 + iy_s];
   cx_im = -xc_im[1 + iy_s];
   }
else
   {
   aang = -2.*PI*(float)(nx/2);
   aang = PI;
   cx_re = cos(aang);
   cx_im = sin(aang);
   } 
/* Tableau de la forme lineaire concerne du noyau */
   ccx_re = 1.; ccx_im = 0.;
   tx_re[0] = 1.; tx_im[0] = 0.;
 
   for(i = 1; i <= ir; i++)
    {
/* CCX=CCX*CX  Complex multiplication. */
     COMPLEX_PROD2(ccx_re,ccx_im,cx_re,cx_im,1,1,&ccx_re,&ccx_im);
     tx_re[i] = ccx_re;
     tx_im[i] = ccx_im;
    }
 
/* Calage en translation */
     for(nb = 0; nb <= nbeta; nb++)
       {
/* XC(NB) = XC(NB) * TX(IXY(1,NB)) */
       COVER_IXY_1D(&iix,&nb);
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
/* XC[NB] = XC[NB] * CX  Complex multiplication... */
       work_re = xc_re[iy_s + nb] * cx_re - xc_im[iy_s + nb] * cx_im;
       work_im = xc_re[iy_s + nb] * cx_im + xc_im[iy_s + nb] * cx_re;
       xc_re[iy_s + nb] = work_re * bmask[iy_s + nb]; 
       xc_im[iy_s + nb] = work_im * bmask[iy_s + nb]; 
      }
/* End of loop on the lines (iy) */
   }
 return(0);
}
/*******************************************************************
* CREYCE_LOAD_BISP_1D
* Load input data from file to YCE array
* yce index will be like fortran arrays: between 1 and ngamma
*******************************************************************/
int CREYCE_LOAD_BISP_1D(float **yce_re, float **yce_im, float *bisp1, 
                        INT4 ngamma, INT4 dim_bispec_input)
{
float xr, xi, xm;
INT4 ng, iy, iy_b, iy_b1, dim_bispec, inull, isize;

dim_bispec = ngamma + 1;

/* COMPLEX YCE(NGMAX) */
   isize = dim_bispec * ny * sizeof(float);
   *yce_re = (float *)malloc(isize);
   *yce_im = (float *)malloc(isize);

for(iy = 0; iy < ny; iy++)
  {
  inull = 0.;
  iy_b = iy * dim_bispec;
  iy_b1 = iy * dim_bispec_input;
  for(ng = 0; ng < ngamma; ng ++)
  {
/* Load real and imaginary parts: */
      xr = bisp1[ng + iy_b1];
      xi = bisp1[ng + ngamma + iy_b1];
#ifdef DEBUG
      if(ng < 5 && iy == 0) printf(" xr, xi, snr: %f %f %f (iy=%d ng=%d)\n",
                                    xr,xi,bisp1[ng+2*ngamma+iy_b1],iy,ng);
#endif
/* Modulus */
      xm = sqrt((double)(xr*xr + xi*xi));
/* Phasor: */
      if(xm == 0.)
         {
         (*yce_re)[ng + iy_b] = 1.;
         (*yce_im)[ng + iy_b] = 0.;
         inull++;
#ifdef DEBUG
         fprintf(fp1," CREYCE_LOAD_BISP_1D/Warning: BISPEC(%d) is null! \n",ng);
         printf(" CREYCE_LOAD_BISP_1D/Warning: BISPEC(%d) is null (line#%d)!\n",ng,iy); 
#endif
         }
      else
         {
         (*yce_re)[ng + iy_b] = xr / xm;
         (*yce_im)[ng + iy_b] = xi / xm;
         }
    }
 
#ifndef DEBUG
  if(inull > 0)
    {
    printf(" CREYCE_LOAD_BISP_1D/Warning: null bispectrum for %d values (line #%d)\n",
           inull,iy);
    fprintf(fp1," CREYCE_LOAD_BISP_1D/Warning: null bispectrum for %d values (line #%d)\n",
           inull,iy);
    }
#endif

/* Debug: */
#ifdef DEBUG
  if(iy == 0)
  {
  for(ng = 0; ng < 5; ng ++)
    { 
     printf(" YCE(%d,iy=%d) = %12.5f %12.5f \n",
                  ng,iy,(*yce_re)[ng + iy_b],(*yce_im)[ng + iy_b]);
     fprintf(fp1," YCE(%d,iy=%d) = %12.5f %12.5f \n",
                  ng,iy,(*yce_re)[ng + iy_b],(*yce_im)[ng + iy_b]);
    }
  }
#endif

/* End of loop on iy: */ 
  }

return(0);
}
/**********************************************************************
* NORMALIZE_L1
* Norm 1 is the initial value of the sum of the absolute values of array.
***********************************************************************/
int NORMALIZE_L1(float *array, INT4 npts, float *norm, int iy)
{
double ssum;
INT4 status;
register int i;

 ssum=0.;
 for(i = 0; i < npts; i++)
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
        printf(" NORMALIZE_L1/Warning: norm is null for line #%d\n",iy);
        fprintf(fp1," NORMALIZE_L1/Warning: norm is null for line #%d\n",iy);
        *norm = 0.;
        status = -1;
      }
    else
      {
       for(i = 0; i < npts; i++) array[i] /= ssum;
       *norm = ssum;
       status = 0;
      }

return(status);
}
/**********************************************************
* To normalize a vector with L2 norm
* Called to normalize the two vectors of the Kernel
* Norm 2 is the square root of the sum of the squared values of array.
* WARNING: index from 1 to npts inclusive!!!
**********************************************************/
int NORMALIZE_L2(float *array, INT4 npts, float *norm, int iy)
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
      printf(" NORMALIZE_L2/Error: norm is null for line #%d\n",iy);
      fprintf(fp1," NORMALIZE_L2/Error: norm is null for line #%d\n",iy);
    }


 *norm = sqrt(ssum);
 for(i = 1; i <= npts; i++) array[i] /= *norm;

return(0);
}
/*****************************************************************
* BISP_WEIGHT0
* Computing the weights:
* JLP Version, not very good...
*****************************************************************/
int BISP_WEIGHT0(INT4 ngamma, INT4 ir, float cte, float *weight)
{
double work;
float cte1, sigb, sig2, wmin, wmax;
INT4 ng, k, kk, iix, iiy;

 cte1 = cte * DEGRAD / (float)(ir*ir);
 
 for(ng = 0; ng < ngamma; ng++)
   {
        wmin = ir * ir;
        wmax = 0.;
        for(kk = 1; kk <= 3; kk++) 
         {
/*  K = KLM(KK,NG)*/
          cover_klm0(&k,kk,ng);
          COVER_IXY_1D(&iix,&k);
          iiy = 0;
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
INT4 ng, is2, k, l, m, kk, iix;

 cte1 = cte * DEGRAD / (float)(ir*ir);
 
    for(ng = 0; ng < ngamma; ng++)
      {
       is2 = 0;
       cover_klm0(&k,1,ng);
       cover_klm0(&l,2,ng);
       cover_klm0(&m,3,ng);
       work = ro[k] * ro[l] * ro[m];
       work *= (bmask[k] * bmask[l] * bmask[m]);
       if(work == 0)
           weight[ng]=0.;
       else
         {
            for(kk = 1; kk <= 3; kk++) 
             {
              cover_klm0(&k,kk,ng);
              COVER_IXY_1D(&iix,&k);
              is2 += iix*iix;
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
             printf(" BISP_WEIGHT11/Warning: Null sigma for NG = %d \n",ng);
             fprintf(fp1," BISP_WEIGHT11/Warning: Null sigma for NG = %d \n",ng);
             weight[ng] = 1.E+9;
             }
         }
    }
 
return(0);
}
/*****************************************************************
* BISP_WEIGHT2
* Computing the weights with SNR  stored in bispectrum file (line 3)
*
******************************************************************/
int BISP_WEIGHT2(float *bisp1, INT4 nbeta, INT4 ngamma, INT4 ir,
                 float cte, float *weight)
{
double sum;
float mean_snr; 
INT4 ng;

 sum = 0.;
 for(ng = 0; ng < ngamma; ng++)
    {
     weight[ng+1] = bisp1[ng + 2 * ngamma];
     if(weight[ng+1] < 0) weight[ng+1] = 0.;
/* JLP93 : keep BISP_SNR, do not put any **0.5 **1.5 or **2, or anything else
*/
     sum += weight[ng+1];
    }
 
/* Diagnostic: */
   mean_snr = sum/(float)ngamma;

#ifdef DEBUG
   printf(" BISP_WEIGHT2/SNR weight Initial sum: %12.5f Mean bisp. SNR: %12.5f \n",
             sum,mean_snr);
   fprintf(fp1," BISP_WEIGHT2/SNR weight Initial sum: %12.5f Mean bisp. SNR: %12.5f \n",
             sum,mean_snr);
#endif
return(0);
}
/*****************************************************************
* BISP_WEIGHT22
* Computing the weights with SNR  stored in bispectrum file (line 3)
* Same as WEIGHT2 but takes modulus into account
* (Not very good... but best results are obtained without calling
*   SIGM2 afterwards, i.e. without dividing by SIGMA or Nber_of_terms)
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
   ng1 = 0;
       for(nb = 2; nb <= nbeta; nb++)
       {
/*  IG2=NGT(NB)
*/
       cover_ngt1(&ng2,nb);
          for(ng = ng1; ng < ng2; ng++)
            {
            weight[ng] = bmask[nb] * ro[nb];
            if(weight[ng] != 0. && wmin > weight[ng]) wmin = weight[ng];
            if(wmax < weight[ng]) wmax = weight[ng];
            }
         ng1 = ng2;
       }

       wrange=wmax/wmin;
#ifdef DEBUG
       printf(" WEIGHT22/Modulus: Initial range %f MIN, MAX: %f %f \n",
              wrange, wmin, wmax);
       fprintf(fp1," WEIGHT22/Modulus: Initial range %f MIN, MAX: %f %f \n",
              wrange, wmin, wmax);
#endif

/* Truncation (max range to 100):
*/
        wmin = wmax / 100.;
        for(ng = 0; ng < ngamma; ng++)
           if(weight[ng] < wmin) weight[ng] = wmin;

/* Diagnostic: */
   mean_snr = sum/(float)ngamma;
#ifdef DEBUG
   printf(" BISP_WEIGHT22/SNR weight Initial sum: %f, Mean bisp. SNR: %f \n",
            sum,mean_snr);
   fprintf(fp1," BISP_WEIGHT22/SNR weight Initial sum: %f, Mean bisp. SNR: %f \n",
            sum,mean_snr);
#endif

return(0);
}
/*******************************************************************
* CREYCE_SIMU_1D : create the phasors of the spectrum and of the bispectrum
* Output in real/imaginary (from artificial bispectrum)
* Centers the Fourier transform
********************************************************************/
int CREYCE_SIMU_1D(float **yce_re, float **yce_im, float **xc_re,
                   float **xc_im, float **xcref_re, 
                   float **xcref_im, INT4 ir, float cte, 
                   INT4 nbeta, INT4 ngamma, 
                   float lower_modulus, char *real_name, char *imag_name,
                   float **weight, float **bmask)
{
char comments[81];
INT4 ixc, nb, ng, inull, iix, is2, isize, line_to_output;
INT4 iy_s, iy_b, irem, k, l, m, kk, dim_spec, dim_bispec; 
INT_PNTR pntr_ima;
register int i, j, iy;
float *re_rev, *im_rev;
float xr, xi, xm, cte1, work, work_re, work_im, xnorm, ro_max;
float cos1, sin1, dgamma;

/* Reading RE and IM : */
  printf(" Real part of the FFT of the image along the columns\n");
/* Warning, input spectroscopic mode: line/columns are reversed: */
  JLP_VM_READIMAG1(&pntr_ima,&ny,&nx,real_name,comments);
  re_rev = (float *)pntr_ima;
  fprintf(fp1," Real part of the FFT: %.14s  Comments: %.30s \n",
              real_name,comments);
  isize = nx * ny * sizeof(float);
  re = (float *) malloc(isize);
  for(j = 0; j < ny; j++)
    for(i = 0; i < nx; i++)
     re[i + j * nx] = re_rev[j + i * ny];

  printf(" Imaginary part of the FFT of the image along the columns\n");
  JLP_VM_READIMAG1(&pntr_ima,&nx,&ny,imag_name,comments);
/* Warning, input spectroscopic mode: line/columns are reversed: */
  im_rev = (float *)pntr_ima;
  fprintf(fp1," Imag. part of the FFT: %.14s  Comments: %.30s \n",
              imag_name,comments);
  isize = nx * ny * sizeof(float);
  im = (float *) malloc(isize);
  for(j = 0; j < ny; j++)
    for(i = 0; i < nx; i++)
     im[i + j * nx] = im_rev[j + i * ny];
 
/* Allocation of memory: */
   dim_spec = nbeta + 1;
   isize = dim_spec * ny * sizeof(float);
   *xc_re = (float *)malloc(isize);
   *xc_im = (float *)malloc(isize);
   *xcref_re = (float *)malloc(isize);
   *xcref_im = (float *)malloc(isize);
   ro = (float *)malloc(isize);
   *bmask = (float *)malloc(isize);
   bad_line = (float *)malloc(isize);

/* Assume that zero frequency is at IXC in each line  */
   ixc = nx / 2;
 
/* Loading modulus and phasor to ro and xc */
  inull = 0;
  line_to_output = 0;
/********************************************************************/
/**************** Loop on all the lines *****************************/
  for(iy = 0; iy < ny; iy++)
    {
     iy_s = iy * dim_spec;
     for(nb = 0; nb <= nbeta; nb++)
        {
/* Get coordinates (iix,0) from nb index: */
         COVER_IXY_1D(&iix,&nb);
/* Add (ixc,0) vector, since centered FFT... */
         iix += ixc;
         xr = re[iix + iy * nx];
         xi = im[iix + iy * nx];
/* Modulus: */
         xm = sqrt(xr * xr + xi * xi);
         ro[nb + iy_s] = xm;
/* Phase factor: */
          if(xm == 0.)
            {
            (*xc_re)[nb + iy_s] = 1.;
            (*xc_im)[nb + iy_s] = 0.;
            inull++;
            }
          else
            {
            (*xc_re)[nb + iy_s] = xr / xm;
            (*xc_im)[nb + iy_s] = xi / xm;
            }

/* End of "for ... nb" loop */
       }

#ifdef DEBUG
    if(inull > 0)
       {
       printf(" CREYCE_SIMU/Modulus null for %d values (iy=%d)\n",inull,iy);
       fprintf(fp1," CREYCE_SIMU/Modulus null for %d values (iy=%d)\n",inull,iy);
       }
#endif
#ifdef DEBUG_2
/* Output of some values of the modulus: */
    if(iy == line_to_output)
      {
      for(nb = 0; nb < 4; nb++)
        {
        printf(" line #%d: RO(%d) = %f \n",iy,nb,ro[nb]);
        fprintf(fp1," line #%d: RO(%d) = %f \n",iy,nb,ro[nb]);
        printf(" line #%d: XC(%d) = (RE, IM): %f %f\n",
                 iy,nb,(*xc_re)[nb],(*xc_im)[nb]);
        fprintf(fp1," line #%d: XC(%d) = (RE, IM): %f %f\n",
                 iy,nb,(*xc_re)[nb],(*xc_im)[nb]);
        }
      }
#endif
 
/* Truncation: */
        irem = 0;
/* Mask to discard some values of the spectral list: */
        ro_max = ro[0 + iy_s];
        if(ro_max <= 0.) ro_max = 1.;
        for(nb = 0; nb <= nbeta; nb++)
        {
          work = ro[nb + iy_s] / ro_max;
/* Check that RO[NB] greater than LOWER_MODULUS:
* Remember that LOWER_MODULUS can be negative... */
          if(work < lower_modulus || work <= 0.)
            {
             (*bmask)[nb + iy_s] = 0.;
             irem++;
            }
          else
            (*bmask)[nb + iy_s] = 1.;
        }
 
#ifdef DEBUG
       if(iy == line_to_output)
        {
        printf(" Number of terms of the spectral list after correction: %d (line#%d)\n",
                 nbeta-irem,iy);
        fprintf(fp1," Number of terms of the spectral list after correction: %d (line #%d)\n",
                 nbeta-irem,iy);
        }
#endif
 
/* End of loop on iy: */
  }
/********************************************************************/

/* Spectrum phasor blocked in translation: xc[.] */
  TRANS_1D(ir,nbeta,*bmask,*xc_re,*xc_im,1);

/* Reference (since simulation) */
  for(iy = 0; iy < ny; iy++)
  {
  iy_s = iy * dim_spec;
  for(nb = 0; nb <= nbeta; nb++)
    {
    (*xcref_re)[nb + iy_s] = (*xc_re)[nb + iy_s];
    (*xcref_im)[nb + iy_s] = (*xc_im)[nb + iy_s];
    }
  }

/***********************************************************************
* Bispectrum phasor: YCE
*/
/* COMPLEX YCE(NGMAX) */
   dim_bispec = ngamma + 1;
   isize = dim_bispec * ny * sizeof(float);
   *yce_re = (float *)malloc(isize);
   *yce_im = (float *)malloc(isize);
   *weight = (float *)malloc(isize);

/* Main loop on all the lines: */
 for(iy = 0; iy < ny; iy++)
   {
    iy_s = iy * dim_spec;
    iy_b = iy * dim_bispec;
    for(ng = 0; ng < ngamma; ng++)
      {
/*  yce[ng]=xc(k)*xc(l)*conjg(xc(m)) */
       cover_klm0(&k,1,ng);
       cover_klm0(&l,2,ng);
       cover_klm0(&m,3,ng);
       COMPLEX_PROD3((*xcref_re)[k + iy_s],(*xcref_im)[k + iy_s],
                     (*xcref_re)[l + iy_s],(*xcref_im)[l + iy_s],
                     (*xcref_re)[m + iy_s],(*xcref_im)[m + iy_s],
                     1,1,-1,&(*yce_re)[ng + iy_b],&(*yce_im)[ng + iy_b]);
       }

/* Just debug mode: 
*/
#ifdef DEBUG_2
     if(iy == line_to_output)
     {
     for(ng = 0; ng < 5; ng++)
       {
        cover_klm0(&k,1,ng);
        cover_klm0(&l,2,ng);
        cover_klm0(&m,3,ng);
       COMPLEX_PROD3((*xcref_re)[k + iy_s],(*xcref_im)[k + iy_s],
                     (*xcref_re)[l + iy_s],(*xcref_im)[l + iy_s],
                     (*xcref_re)[m + iy_s],(*xcref_im)[m + iy_s],
                     1,1,-1,&work_re,&work_im);
        printf(" line #%d: k=%d l=%d m=%d yce(%d) = (RE, IM): %f %f \n",
                         iy,k,l,m,ng,(*yce_re)[ng+iy_b],(*yce_im)[ng+iy_b]);
       printf(" work_re=%e work_im=%e \n",work_re,work_im);
       }
    }
#endif
 
/* End of loop on iy */
  }

/* Computing the weights: */
     if(cte > 0)
       {
        for(iy = 0; iy < ny; iy++)
          {
           iy_s = iy * dim_spec;
           iy_b = iy * dim_bispec;
           BISP_WEIGHT11(nbeta,ngamma,ir,cte,
                         &(*bmask)[iy_s],&ro[iy_s],&(*weight)[iy_b]);
           NORMALIZE_L1(&(*weight)[iy_b],ngamma,&xnorm,iy);
           bad_line[iy] = (xnorm > 0) ? 0 : 1;
          }
       }
     else if (cte < 0)
       {
/* Version with SNR stored in 3rd line of bispectrum: */
        printf(" CREYCE_SIMU/Fatal error, this option (cte<0) is not possible here\n");
        fprintf(fp1," CREYCE_SIMU/Fatal error, this option (cte<0) is not possible here\n");
        fclose(fp1); exit(-1);
       }
     else 
       {
        printf(" Weights set to unity, and then normalized \n");
        fprintf(fp1," Weights set to unity, and then normalized \n");
        for(iy = 0; iy < ny; iy++)
          {
          iy_s = iy * dim_spec;
          iy_b = iy * dim_bispec;
          for(ng = 0; ng < ngamma; ng++)
             {
             (*weight)[ng + iy_b]=1.;
             for(kk = 1; kk <= 3; kk++) 
               {
/*  K = KLM(KK,NG)*/
               cover_klm0(&k,kk,ng);
               if((*bmask)[k + iy_s] == 0.) (*weight)[ng + iy_b]=0.;
               }
             }
         NORMALIZE_L1(&(*weight)[iy_b],ngamma,&xnorm,iy);
         bad_line[iy] = (xnorm > 0) ? 0 : 1;
         }
       }
 

if(cte != 0)
  {
/* Perturbation of YCE */
  cte1 = cte * DEGRAD / (float)(ir*ir);
 
  for(iy = 0; iy < ny; iy++)
   {
   iy_b = iy * dim_bispec;
   for(ng = 0; ng < ngamma; ng++)
     {
     is2 = 0;
     for(kk = 1; kk <= 3; kk++) 
       {
/*  K = KLM(KK,NG)*/
       cover_klm0(&k,kk,ng);
/* IS2 = IS2 + IXY(1,K)**2 + IXY(2,K)**2
*/
       COVER_IXY_1D(&iix,&k);
       is2 += (iix*iix);
       }

/* Random generation (Gaussian law, (1.,0.)) of DGAMMA = DELTA GAMMA */
     JLP_RANDOM_GAUSS(&work);
     dgamma = work * cte1 * (float)is2;
     cos1=cos((double)dgamma);
     sin1=sin((double)dgamma);
/* yce[ng] = yce[ng]*cmplx(cos(dgamma),sin(dgamma)) */
     work_re = (*yce_re)[ng + iy_b] * cos1 - (*yce_im)[ng + iy_b] * sin1;
     work_im = (*yce_re)[ng + iy_b] * sin1 + (*yce_im)[ng + iy_b] * cos1;
     (*yce_re)[ng + iy_b] = work_re;
     (*yce_im)[ng + iy_b] = work_im;
     }
   }
/* End of (cte != 0) */
  }
 
return(0);
}
/*******************************************************************
* CREYCE1
* Reads a real spectrum (square modulus only) and bispectrum
* (Case of observations)
*******************************************************************/
int CREYCE1_1D(float **yce_re, float **yce_im, float **xc_re,
               float **xc_im, float **xcref_re,
               float **xcref_im, INT4 ir, float cte, INT4 nbeta, 
               INT4 ngamma, INT4 ngamma_max, 
               float lower_modulus, char *modsq_name, char *bisp_name, 
               float **weight, float **bmask)
{
char comments[81];
INT4 ixc, nb, ng, k, kk, inull, iix, iy, iy_s, iy_b, nx1, ny1;
INT4 irem, isize, dim_spec, dim_bispec; 
INT4 line_to_output;
INT_PNTR pntr_ima;
float *modsq, *bisp1, *modsq_rev, xnorm, ro_max;
double work;
register int i, j;

/* Input of the modulus: */
  printf(" Modulus (squared) of the FFT of the image along the columns\n");
/* Warning, input spectroscopic mode: line/columns are reversed: */
  JLP_VM_READIMAG1(&pntr_ima,&ny,&nx,modsq_name,comments);
  modsq_rev = (float *)pntr_ima;
  fprintf(fp1," Sq. modulus: %.14s  Comments: %.30s\n",modsq_name,comments);
   isize = nx * ny * sizeof(float);
   modsq = (float *)malloc(isize);
/* Warning, input spectroscopic mode: line/columns are reversed: */
  for(j = 0; j < ny; j++)
    for(i = 0; i < nx; i++)
     modsq[i + j * nx] = modsq_rev[j + i * ny];
 
/* Allocation of memory: */
   isize = nx * ny * sizeof(float);
   re = (float *)malloc(isize);
   im = (float *)malloc(isize);
   dim_spec = nbeta + 1;
   isize = dim_spec * ny * sizeof(float);
   *xc_re = (float *)malloc(isize);
   *xc_im = (float *)malloc(isize);
   *xcref_re = (float *)malloc(isize);
   *xcref_im = (float *)malloc(isize);
   ro = (float *)malloc(isize);
   *bmask = (float *)malloc(isize);
   bad_line = (float *)malloc(isize);

/* Initializing the real and imaginary part of the spectrum: */
  inull = 0;
  for(j = 0; j < ny; j++)
    {
    for(i = 0; i < nx; i++)
        {
           work = modsq[i + j * nx];
           if(work < 0.)
             {
              work = 0.; 
              inull++;
              }
           re[i + j * nx] = sqrt(work);
           im[i + j * nx] = 0.;
        }
    }
  if(inull > 0)
    {
     printf(" CREYCE1_1D/Modulus negative or null for %d values\n",inull);
     fprintf(fp1," CREYCE1_1D/Modulus negative or null for %d values\n",inull);
    }

 
/* Assume that zero frequency is at IXC,IY:  */
   ixc = nx / 2;
 
/* Loading modulus and phasor to ro and xc */
/* Loop on iy: */
  for(iy = 0; iy < ny; iy++)
   {
     iy_s = iy * dim_spec;
     for(nb = 0; nb <= nbeta; nb++)
        {
/* Get coordinates (iix,0) from nb index: */
         COVER_IXY_1D(&iix,&nb);
/* Add (ixc,0) vector, since centered FFT... */
         iix += ixc;
/* Modulus (here it is the real part) : */
         ro[nb + iy_s] = re[iix + iy * nx];
/* Initialization of the phase factor: */
         (*xc_re)[nb + iy_s] = 1.;
         (*xc_im)[nb + iy_s] = 0.;
/* As the reference is not available, I set it to the null phase: 
* (used in RECURSIVE as the initial value of the first terms) */
         (*xcref_re)[nb + iy_s] = 1.;
         (*xcref_im)[nb + iy_s] = 0.;
        }
 
/* Output some values of the modulus: */
    line_to_output = 0;
#ifdef DEBUG_2
    if(iy == line_to_output)
    {
    for(nb = 0; nb < 5; nb++)
       {
       printf(" iy=%d RO(%d) = %f (line #%d)\n",iy,nb,ro[nb],iy);
       fprintf(fp1," iy=%d RO(%d) = %f (line #%d)\n",iy,nb,ro[nb],iy);
       }
    }
#endif
 
 
/* Truncation: */
    irem = 0;
/* Mask to discard some values of the spectral list: */
    ro_max = ro[0 + iy_s];
    if(ro_max <= 0.) ro_max = 1.;
    for(nb = 0; nb <= nbeta; nb++)
    {
      work = ro[nb + iy_s] / ro_max;
/* Check that RO[NB] greater than LOWER_MODULUS:
* Remember that LOWER_MODULUS can be negative... */
      if(work < lower_modulus || work <= 0.)
        {
         (*bmask)[nb + iy_s] = 0.;
         irem++;
        }
      else
        (*bmask)[nb + iy_s] = 1.;
   }
 
#ifdef DEBUG
    if(iy == line_to_output)
    {
    printf(" lower_modulus = %f \n",lower_modulus);
    printf(" Number of terms of the spectral list after correction: %d (line #%d)\n",
             nbeta-irem,iy);
    fprintf(fp1," Number of terms of the spectral list after correction: %d (line #%d)\n",
             nbeta-irem,iy);
    }
#endif
 
/* End of loop on iy: */
  }

/* Spectrum phasor blocked in translation: xc[.] */
    TRANS_1D(ir,nbeta,*bmask,*xcref_re,*xcref_im,1);

/***********************************************************************
* Bispectrum phasor: YCE
*/
    printf(" Bispectrum of the image (phase term) (NOT centered in the frame) \n");
/* Three lines (real, imag, snr) for the input bispectrum: */
    JLP_VM_READIMAG1(&pntr_ima,&nx1,&ny1,bisp_name,comments);
    bisp1 = (float *)pntr_ima;
    fprintf(fp1," Bispectrum: %.14s  Comments: %.30s\n",bisp_name,comments);
    if( (ny1 != ny) || (nx1 != 3 * ngamma_max))
      {
      printf("CREYCE1_1D/FATAL ERROR: Size of bispectrum inconsistent with IRMAX\n");
      fprintf(fp1,"CREYCE1_1D/FATAL ERROR: Size of bispectrum inconsistent with IRMAX\n");
      fclose(fp1); exit(-1);
      }

/* Load BISP array to YCE array: */
    CREYCE_LOAD_BISP_1D(yce_re,yce_im,bisp1,ngamma,3*ngamma);

/* REAL WEIGHT(NGMAX) */
    dim_bispec = ngamma + 1;
    isize = dim_bispec * ny * sizeof(float);
    *weight = (float *)malloc(isize);
    if((*weight) == NULL)
     {
      printf("CREYCE1_1D/Fatal error allocating memory for \"weight\" \n");
      exit(-1);
     }

/* Computing the weights: */
    if(cte > 0)
      {
        printf(" cte=%f, ir=%d, nbeta=%d ngamma=%d\n",cte,ir,nbeta,ngamma);
        for(iy = 0; iy < ny; iy++)
          {
           iy_s = iy * dim_spec;
           iy_b = iy * dim_bispec;
           BISP_WEIGHT11(nbeta,ngamma,ir,cte,
                         &(*bmask)[iy_s],&ro[iy_s],&(*weight)[iy_b]);
           NORMALIZE_L1(&(*weight)[iy_b],ngamma,&xnorm,iy);
           bad_line[iy] = (xnorm > 0) ? 0 : 1;
          }
      }
    else if (cte < 0)
      {
        for(iy = 0; iy < ny; iy++)
          {
           iy_b = iy * dim_bispec;
/* Version with SNR stored in 3rd line of bispectrum: */
           BISP_WEIGHT2(&bisp1[iy_b],nbeta,ngamma,ir,cte,&(*weight)[iy_b]);
/*         BISP_WEIGHT1(nbeta,ngamma,ir,cte,(*bmask),&(*weight)[iy_b]);
*/
           NORMALIZE_L1(&(*weight)[iy_b],ngamma,&xnorm,iy);
           bad_line[iy] = (xnorm > 0) ? 0 : 1;
           if(xnorm <= 0) printf("CREYCE1_1D/bad_line #%d\n",iy);
          }
      }
    else 
      {
       printf(" Weights set to unity, and then normalized \n");
       fprintf(fp1," Weights set to unity, and then normalized \n");
       dim_bispec = ngamma + 1;
       for(iy = 0; iy < ny; iy++)
         {
         iy_s = iy * dim_spec;
         iy_b = iy * dim_bispec;
         for(ng = 0; ng < ngamma; ng++)
           {
           (*weight)[ng + iy_b]=1.;
           for(kk = 1; kk <= 3; kk++) 
              {
/*  K = KLM(KK,NG)*/
              cover_klm0(&k,kk,ng);
              if((*bmask)[k + iy_s] == 0.) (*weight)[ng + iy_b]=0.;
              }
          }
        NORMALIZE_L1(&(*weight)[iy_b],ngamma,&xnorm,iy);
        bad_line[iy] = (xnorm > 0) ? 0 : 1;
        }
     }
 
return(0);
}
/*******************************************************************
* RECURSIVE_1D: Initial solution of the phasor.
*                      YCE(.)  to  XC(.)
*
* exp(i*YCE[NG]) = exp (i * phase_K) * exp (i * phase_L) * exp (-i * phase_M)
* with L=NB or M=NB
*******************************************************************/
int RECURSIVE_1D(float *sigm, INT4 nbeta, INT4 ngamma, float sigma_null, 
                 INT4 ifermax, float *bmask, float *yce_re, float *yce_im, 
                 float *xc_re, float *xc_im, float *xcref_re, float *xcref_im, 
                 float *weight)
{
float cc1_re, cc1_im;
/* Double precision is necessary for large sums !!!
*/
long double sumxr1, sumxi1, sumsqr, sumsqi, sigr, sigi, sumweights;
long double rnorm;
float xw;
int nb_maxi;
INT4 k, l, m, ng1, ng2, nb, ng, indx, nval, inull;
INT4 iy, iy_s, iy_b, dim_spec, dim_bispec, nbad_lines;
 
/* Main loop on iy: */
nbad_lines = 0;
dim_spec = nbeta + 1;
dim_bispec = ngamma + 1;
for(iy = 0; iy < ny; iy++)
  {
  iy_s = iy * dim_spec;
  iy_b = iy * dim_bispec;
  if(bad_line[iy]) 
    {
     nbad_lines++;
/* In case of a bad line, will set all phase values to zero: */ 
     nb_maxi = nbeta;
    }
  else
/* First 2 values remain undermined (dim(Ker A) = 1 when 1D case): */ 
    nb_maxi = 1;

/* Initialize the first nb_maxi values to null phase: */
   for(nb = 0; nb <= nb_maxi; nb++)
    {
#ifdef USE_REFERENCE_PHASE
     xc_re[nb + iy_s] = xcref_re[nb + iy_s];
     xc_im[nb + iy_s] = xcref_im[nb + iy_s];
#else
     xc_re[nb + iy_s] = 1.;
     xc_im[nb + iy_s] = 0.;
#endif
    }

/* Process the bispectrum if it is not a bad line: */

  if(!bad_line[iy]) 
  {
/* Initialize the first three values of sigma, to avoid problems... */
   for(nb = 0; nb <= 2; nb++) sigm[nb + iy_s] = MIN_SIGMA;


  ng1 = 0;
  inull = 0;
 
  for(nb = 2; nb <= nbeta; nb++)
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
   for(ng = ng1; ng < ng2; ng++)
     {
     indx++;
     if(indx <= ifermax)
       {
/*  K=KLM(1,NG) L=KLM(2,NG) M=KLM(3,NG)
*/
        cover_klm0(&k,1,ng);
        cover_klm0(&l,2,ng);
        cover_klm0(&m,3,ng);
/* JLP2000
        xw = sqrt((double)weight[ng + iy_b]);
*/
        xw = weight[ng + iy_b];
        if(xw > 0.)
         {
          nval++;
          sumweights += xw;
          if(m == nb)
             {
/* Case 1 (m = nb): cc1=xc(k)*xc(l)*conjg(yce[ng]) */
             COMPLEX_PROD3(xc_re[k + iy_s],xc_im[k + iy_s],
                           xc_re[l + iy_s],xc_im[l + iy_s],
                           yce_re[ng + iy_b],yce_im[ng + iy_b],
                           1,1,-1,&cc1_re,&cc1_im);
             }
          else
             {
/* Case 2 (l = nb): cc1=yce[ng]*conjg(xc(k))*xc(m) */
             COMPLEX_PROD3(xc_re[k + iy_s],xc_im[k + iy_s],
                           xc_re[m + iy_s],xc_im[m + iy_s],
                           yce_re[ng + iy_b],yce_im[ng + iy_b],
                           -1,1,1,&cc1_re,&cc1_im);
             }
#ifdef DEBUG_2
/* Diagnostic to check the validity of the process: */
         if(iy == 0 && nb == 2)
           {
           printf("-----------------------------------\n");
           printf(" k=%d l=%d m=%d ng1=%d ng2=%d iy_s=%d iy_b=%d\n",
                  k,l,m,ng1,ng2,iy_s,iy_b); 
           printf(" xc_re[k]=%e xc_im[k]=%e yce_re[0]=%e yce_im[0]=%e\n",
                    xc_re[k],xc_im[k],yce_re[ng],yce_re[ng]);
           printf("xcref_re[2]=%e xcref_im[2]=%e cc1_re=%e cc1_im=%e xw=%e\n",
                   xcref_re[2],xcref_im[2],cc1_re,cc1_im,xw);
           printf("-----------------------------------\n");
           }
#endif
          sumsqr += (cc1_re * cc1_re * xw);
          sumsqi += (cc1_im * cc1_im * xw);
          sumxr1 += (cc1_re * xw);
          sumxi1 += (cc1_im * xw);
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
        xc_re[nb + iy_s] = sumxr1/rnorm;
        xc_im[nb + iy_s] = sumxi1/rnorm;

/* JLP94: sigma has a meaning only if more than 2 values are involved in the
* mean... */
        if(nval > 1)
          {
/* SIGR**2 = Sum of (weight_i *(X_i - mean)**2) / Sum of (weight_i) 
* or Sum of (weight_i X_i**2)/Sum of weight_i - mean**2
* Here the weighted mean of real(CC) is SUMXR1/SUMWEIGHTS:
*/
           sigr = sumsqr / sumweights
                - (sumxr1 / sumweights) * (sumxr1 / sumweights);
           sigi = sumsqi / sumweights
                 - (sumxi1 / sumweights) * (sumxi1 / sumweights);
/* Note that double precision variables are needed for large sums !!!
*/
           if(sigr < 0. || sigi < 0.)
/* Just in case of a (round-off??) problem:
*/
              {
/* Error message if absolute value is large: */
               if(sigr < -1e-10 || sigi < -1e-10)
               {
               printf(" Round-off problem in RECURSIVE... for NB= %d (iy=%d)\n",
                      nb,iy);
               fprintf(fp1," Round-off problem in RECURSIVE... for NB= %d (iy=%d)\n",
                      nb,iy);
               printf(" sigi: %g sigr: %g \n",(float)sigi,(float)sigr);
               fprintf(fp1," sigi: %g sigr: %g \n",(float)sigi,(float)sigr);
               printf(" sumweights: %g sumsqr: %g sumsqi: %g \n",
                      (float)sumweights,(float)sumsqr,(float)sumsqi);
               printf(" ng1=%d, ng2=%d sumxr1=%e sumxi1=%e rnorm=%e\n",
                       ng1,ng2,(float)sumxr1,(float)sumxi1,(float)rnorm);
               printf(" weight(ng1): %f  weight(ng2): %f\n",
                        weight[ng1 + iy_s],weight[ng2 + iy_s]);
               sigm[nb + iy_s] = sigma_null;
               }
/* If very small, no error message */
               else
                sigm[nb + iy_s] = 0.;
               }
           else
               {
               sigm[nb + iy_s] = sqrt(sigi+sigr);
               }
/* End of sigr < 0 or sigi < 0 */

          }
/* If only one value has been taken into account for the mean computation:
*/
       else
          {
          sigm[nb + iy_s] = sigma_null;
          }
       }
/* rnorm = 0, i.e., the recursive process has been interrupted: */
     else
       {
/* Diagnostic if it is not in the list of null modulus: */
           if(bmask[nb + iy_s] != 0.)
             {
/* Then add this frequency to BMASK, since it won't be able to recover the
* phase of this frequency:
*/
              bmask[nb + iy_s] = 0.;
/* JLP99/Debug: */
if(iy < 5)
{
              printf(" RECURSIVE/Warning: too many null modulii: the phase of xc[%d] will remain undetermined.\n",nb);  
              fprintf(fp1," RECURSIVE/Warning: too many null modulii: the phase of xc[%d] will remain undetermined.\n",nb);  
/* Print the values of the modulii which may be null: */
              printf(" ng1=%d, ng2=%d sumxr1=%e sumxi1=%e rnorm=%e\n",
                       ng1,ng2,(float)sumxr1,(float)sumxi1,(float)rnorm);
              for(ng = ng1; ng < ng2; ng++)
                {
                cover_klm0(&k,1,ng);
                cover_klm0(&l,2,ng);
                cover_klm0(&m,3,ng);
                xw = weight[ng + iy_b];
                printf(" ng=%d, k=%d, l=%d, m=%d xw=%e\n",ng,k,l,m,xw);
                printf(" xc_re[k]=%e, xc_re[l]=%e, xc_re[m]=%e yce_re[ng]=%e\n",
                         xc_re[k+iy_s],xc_re[l+iy_s],xc_re[m+iy_s],
                         yce_re[ng + iy_b]);
                printf(" xc_im[k]=%e, xc_im[l]=%e, xc_im[m]=%e yce_im[ng]=%e\n",
                         xc_im[k+iy_s],xc_im[l+iy_s],xc_im[m+iy_s],
                         yce_re[ng + iy_b]);
                printf(" ro[k]=%e, ro[l]=%e, ro[m]=%e\n",
                         ro[k+iy_s],ro[l+iy_s],ro[m+iy_s]);
                if(m == nb)
                    COMPLEX_PROD3(xc_re[k + iy_s],xc_im[k + iy_s],
                                  xc_re[l + iy_s],xc_im[l + iy_s],
                                  yce_re[ng + iy_b],yce_im[ng + iy_b],
                                  1,1,-1,&cc1_re,&cc1_im);
                else
                    COMPLEX_PROD3(xc_re[k + iy_s],xc_im[k + iy_s],
                                  xc_re[m + iy_s],xc_im[m + iy_s],
                                  yce_re[ng + iy_b],yce_im[ng + iy_b],
                                  -1,1,1,&cc1_re,&cc1_im);
                printf(" cc1_re=%e, cc1_im=%e\n",cc1_re,cc1_im);
                }
/* JLP99/debug */
}
           }
/* End of bmask != 0 */
         inull++;
/* We arbitrarily set the phase to zero: */
         xc_re[nb + iy_s] = 1.;
         xc_im[nb + iy_s] = 0.;
/* Maximum sigma is one, but I set it to SIGMA_NULL since we assume that
* phase is not important if modulus is too small
*/
         sigm[nb + iy_s] = sigma_null;
         }
/* End of interruption of recursive process */
   ng1 = ng2;
/* End of nb loop */
   }
 
/* Set the first three values to zero, to avoid problems... */
   for(nb = 0; nb <= 2; nb++) sigm[nb + iy_s] = MIN_SIGMA;

/* Diagnostic: */
 if(inull > 0 && iy < 10)
   {
   printf(" RECURSIVE/Warning: Null modulus (i.e., the spectral list has been reduced) for %d values (iy=%d)\n",inull,iy);
   fprintf(fp1," RECURSIVE/Warning: Null modulus (i.e., the spectral list has been reduced) for %d values (iy=%d)\n",inull,iy);
   printf(" WARNING: In that case SIGMA was set to %f (and SNR to %f) \n",
           sigma_null, 1./sigma_null); 
   fprintf(fp1," WARNING: In that case SIGMA was set to %f (and SNR to %f) \n",
           sigma_null, 1./sigma_null); 
   }
 
/* End of case !bad_line[iy]: */
   }
/* End of loop on iy: */
 }

#ifdef DEBUG
 printf(" Bad lines: %d over %d lines\n",nbad_lines,ny);
 fprintf(fp1," Bad lines: %d over %d lines\n",nbad_lines,ny);
#endif

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
                  int sign1, int sign2, float *c12_re, float *c12_im)
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
                  float *c123_re, float *c123_im)
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
* WEIGHTS_SIGM2: Weights according to Sigma of recursive solution
* Specially designed with weights according to snr_bispectrum
*******************************************************************/
int WEIGHTS_SIGM2(float *sigm, INT4 nbeta, INT4 ngamma, float sig_max,
                 float lower_modulus, INT4 ifermax, float *bmask,
                 float *weight)
{
INT4 ng1, ng2, ng, nb, irem, iy, iy_s, iy_b, dim_spec, dim_bispec;
float xnorm;
 
/* Loop on all the lines: */
dim_spec = nbeta + 1;
dim_bispec = ngamma + 1;
for(iy = 0; iy < ny; iy++)
{
iy_s = iy * dim_spec;
iy_b = iy * dim_bispec;
ng1 = 0;
irem = 0;
sigm[3+iy_s] = MIN_SIGMA;
sigm[4+iy_s] = MIN_SIGMA;
 
for(nb = 2; nb <= nbeta; nb++)
   {
   cover_ngt1(&ng2,nb);
/* Caution to avoid SIGM = 0, and problems when computing 1/SIGM in SNR map:
*/
   if(sigm[nb+iy_s] < MIN_SIGMA) sigm[nb+iy_s] = MIN_SIGMA;

/* Truncation to reduce the range of weights: */
   if(sigm[nb+iy_s] >= sig_max || bmask[nb + iy_s] == 0.)
      {
      irem++;
      bmask[nb + iy_s] = 0.;
      for(ng = ng1; ng < ng2; ng++) weight[ng + iy_b] = 0.;
      }
   else
      {
        for(ng = ng1; ng < ng2; ng++)
          {
/* Previous weight divided by this sigma: 
* I have tried without dividing: gives worse results 
* I have tried with **1, **0.8, **0.2 instead: gives worse results 
* JLP93
*/
           weight[ng+iy_b] /= sqrt((double)sigm[nb+iy_s]);
/* Normalisation to the number of bispectral terms:
* I have tried without dividing, and it gives better results 
* JLP93               weight[ng+iy_b] /= (float)(ng2-ng1);
*/
/* End of ng loop */
          }
/* End of (sigm[nb+iy_b] >= sig_max || bmask[nb+iy_s] == 0) */
     }
  ng1 = ng2;
/* End of nb loop */
  }
 
/* Normalisation of the weights: */
NORMALIZE_L1(&weight[iy_b],ngamma,&xnorm,iy);
bad_line[iy] = (xnorm > 0) ? 0 : 1;

/* Diagnostic: */
  if(irem > 0 && iy < 10)
  {
   printf("WEIGHTS_SIGM2/Warning: %d discarded BETA values because of modulus < %f or sigma > %f (line #%d)\n",
          irem,lower_modulus,sig_max,iy); 
   fprintf(fp1,"WEIGHTS_SIGM2/Warning: %d discarded BETA values because of modulus < %f or sigma > %f (line #%d)\n",
          irem,lower_modulus,sig_max,iy); 
  }

/* End of loop on the lines (iy): */
  }

return(0);
}
/*******************************************************************
* WEIGHTS_SIGM1: Weights according to Sigma of recursive solution
* specially designed for unity weights
*******************************************************************/
int WEIGHTS_SIGM1(float *sigm, INT4 nbeta, INT4 ngamma, float sig_max,
                  float lower_modulus, INT4 ifermax, float *bmask,
                  float *weight)
{
FILE *fp2;
INT4 ng1, ng2, ng, nb, irem, iy, iy_s, iy_b, dim_spec, dim_bispec;
INT4 line_to_output;
float xnorm;
 
/* Opening sigma file: */
if((fp2 = fopen("sigma.dat","w")) == NULL)
   {
   printf("WEIGHTS_SIGM1/Fatal error opening sigma file \"sigma.dat\"\n");
   fprintf(fp1,"WEIGHTS_SIGM1/Fatal error opening sigma file \"sigma.dat\"\n");
   fclose(fp1); exit(-1);
   }
line_to_output = 0;
  
/* Loop on all the lines: */
dim_spec = nbeta + 1;
dim_bispec = ngamma + 1;
for(iy = 0; iy < ny; iy++)
{
iy_s = iy * dim_spec;
iy_b = iy * dim_bispec;

ng1 = 0;
irem = 0;
sigm[3 + iy_s] = MIN_SIGMA;
sigm[4 + iy_s] = MIN_SIGMA;
 
for(nb = 2; nb <= nbeta; nb++)
   {
   cover_ngt1(&ng2,nb);
/* Caution to avoid SIGM = 0, and problems when computing 1/SIGM in SNR map:
*/
   if(sigm[nb + iy_s] < MIN_SIGMA) sigm[nb + iy_s] = MIN_SIGMA;
   if(iy == line_to_output) fprintf(fp2,"%d %f \n",nb,sigm[nb + iy_s]);

/* Truncation to reduce the range of weights: */
   if(sigm[nb + iy_s] >= sig_max || bmask[nb + iy_s] == 0.)
      {
      irem++;
      bmask[nb + iy_s] = 0.;
      for(ng = ng1; ng < ng2; ng++) weight[ng + iy_b] = 0.;
      }
   else
      {
        for(ng = ng1; ng < ng2; ng++)
          {
/* JLP93
* (Unity weights: gives worse results when not dividing by SIGM)
*/
           weight[ng + iy_b] /= sigm[nb + iy_s];
/*
* Normalisation to the number of bispectral terms:
* (Unity weights: gives worse results when dividing by nber_terms**2 only)
*/
           weight[ng + iy_b] /= (float)(ng2-ng1)*(ng2-ng1);
/* End of ng loop */
          }
/* End of if */
     }
/* End of (sigm[nb + iy_s] >= sig_max || bmask[nb + iy_s] == 0) */
  ng1 = ng2;
  }
/* End of nb loop */
 
/* Close sigma file: */
 if(iy == line_to_output) fclose(fp2);

/* Normalisation of the weights: */
NORMALIZE_L1(&weight[iy_b],ngamma,&xnorm, iy);
bad_line[iy] = (xnorm > 0) ? 0 : 1;

/* Diagnostic: */
  if(irem > 0 && iy < 10)
  {
   printf("WEIGHTS_SIGM1/Warning: %d discarded BETA values because of modulus < %f or sigma > %f  at line #%d\n",
          irem,lower_modulus,sig_max,iy); 
   fprintf(fp1,"WEIGHTS_SIGM1/Warning: %d discarded BETA values because of modulus < %f or sigma > %f  at line #%d\n",
          irem,lower_modulus,sig_max,iy); 
  }
/* End of loop on iy: */
}

return(0);
}
/**********************************************************************
* Output SNR map according to the adopted weights
* Generating SNR and sigma maps (needed by DIANE)
**********************************************************************/
int OUTPUT_SNR1_1D(float *sigm, INT4 nbeta, char *fname)
{
char name[60],comments[81];
INT4 nb, ng, ix, iy, ixc, iix, iy_s, dim_spec;
register int i, j;
float sigg, *work;
 
/* Assume that zero frequency is at IXC,IY:  */
   ixc = nx / 2;
 
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
 
/* Generating SNR and sigma maps (needed for the deconvolution by DIANE):
*/
 dim_spec = nbeta + 1;
 for(iy = 0; iy < ny; iy++)
   {
   iy_s = iy * dim_spec;
   for(nb = 0; nb <= nbeta; nb++)
      {
/* Get coordinates (iix,0) from nb index: */
      COVER_IXY_1D(&ix,&nb);
/* Add (ixc,0) vector, since centered FFT... */
      iix = ix + ixc;
      sigg = sigm[nb + iy_s];
      if(sigg < MIN_SIGMA) sigg = MIN_SIGMA;
      re[iix + iy * nx] = 1. / sigg;
      im[iix + iy * nx] = sigg * ro[nb + iy_s];
/* Same values for the symmetrical frequency relative to the zero frequency: */
      iix = ixc -ix;
      re[iix + iy * nx] = 1. / sigg;
      im[iix + iy * nx] = sigg * ro[nb + iy_s];
      }
   }
 
/* Output of SNR image: 
*/

/* Switch X/Y for output: */
    work = (float *)malloc(nx * ny * sizeof(float));
    for(iy = 0; iy < ny; iy++)
     for(ix = 0; ix < nx; ix++)
       work[iy + ix * ny] = re[ix + iy * nx];

    strcpy(name,"snr1");
    sprintf(comments," 1/SIGMA_PHASE of : %.20s",fname);
    printf(" Output of SNR=1/SIGMA_PHASE in image file: %s \n",name);
    JLP_WRITEIMAG(work,&ny,&nx,&ny,name,comments);

/* Output of SIGMA
*/

/* Switch X/Y for output: */
    work = (float *)malloc(nx * ny * sizeof(float));
    for(iy = 0; iy < ny; iy++)
     for(ix = 0; ix < nx; ix++)
       work[iy + ix * ny] = im[ix + iy * nx];

    strcpy(name,"sigma");
    sprintf(comments," MODULUS * SIGMA_PHASE of : %.20s",fname);
    printf(" Output of MODULUS * SIGMA_PHASE in image file: %s \n",name);
    JLP_WRITEIMAG(work,&ny,&nx,&ny,name,comments);

free(work);

return(0);
}
/*******************************************************************
* ERROR_SIMU_1D
* Compute the absolute errors for simulations
* (since model is known and stored in xcref)
*
* E is the mean of the errors of the phase terms 
*       weighted by the value of MODSQ**2:
*   or the mean of the errors of the full term (amplitude and phase)
* E = sqrt( Sum of (RO**2 * |xcref-XC|**2) / Sum of RO**2)
* whereas EP is raw:
* EP = sqrt( Sum of |xcref-XC|**2 / NBETA-2)
*******************************************************************/
int ERROR_SIMU_1D(float *xc_re, float *xc_im, float *xcref_re, float *xcref_im, 
                  INT4 nbeta, float *err_spec, float *err_phas,
                  char *errname, float *bmask)
{
FILE *fp2;
float cc_re, cc_im;
double r2, rr2, sr2;
INT4 indx, nb, iy, iy_s, dim_spec, line_to_output;
 
/* Output errors for central line only: */
line_to_output = 0;
iy = line_to_output;
dim_spec = nbeta + 1;
iy_s = iy * dim_spec;

if((fp2 = fopen(errname,"w")) == NULL)
   {
   printf("ERROR_SIMU_1D/Fatal error opening error file %s\n",errname);
   fprintf(fp1,"ERROR_SIMU_1D/Fatal error opening error file %s\n",errname);
   fclose(fp1); exit(-1);
   }
  
fprintf(fp2," 0 0. 0. 0. Index, Phase error, Full quad. error, Cumul. mean phase error (line #%d)\n",
        iy);

/* Global errors:
* err_spec : rms error of the spectrum
* err_phas: rms error of the phasor of the spectrum
*/
 *err_spec = 0.;
 *err_phas = 0.;
 sr2 = 0.;
 indx = 0;
 for(nb = 2; nb <= nbeta; nb++)
  {
   if(bmask[nb + iy_s] != 0.)
     {
      indx++;
      r2 = ro[nb + iy_s]*ro[nb +iy_s];
      sr2 += r2;
      cc_re = xcref_re[nb + iy_s] - xc_re[nb + iy_s];
      cc_im = xcref_im[nb + iy_s] - xc_im[nb + iy_s];
      rr2 = cc_re * cc_re + cc_im * cc_im;
      *err_spec += r2 * rr2;
      *err_phas += rr2;
      if(nb > 2) fprintf(fp2,"%d %12.5e %12.5e %12.5e\n",nb,sqrt(rr2),sqrt(r2*rr2),
              sqrt((double)(*err_phas/(double)(nb-2))));
     }
  }
 
  if(sr2 == 0.)
    {
     *err_spec = 1000.;
     printf(" ERROR_SIMU_1D/Warning: sr2 = 0. !! \n");
     fprintf(fp1," ERROR_SIMU_1D/Warning: sr2 = 0. !! \n");
    }
/* err_spec = sqrt( Sum of (ro**2 * |xcref-xc|**2) / Sum of ro**2) */
   else
     *err_spec = sqrt((double)(*err_spec/sr2));

/* err_phas = sqrt( Sum of |xcref-xc|**2 / nbeta-2) */
     *err_phas = sqrt((double)(*err_phas/(double)indx));
 
fclose(fp2);

return(0);
}
/*******************************************************************
* EPLUS_1D projects ALPHA_0 onto E^+
* xx(idim_xx,4)
* vecxy(idim_xx) 
*******************************************************************/
int EPLUS_1D(float *xc_re, float *xc_im, INT4 nbeta, float *bmask, 
             float *xx, INT4 idim_xx, float *vecxy)
{
float cc_re, cc_im, xxc, cos1, sin1, alpha;
INT4 nb, iy, iy_s, dim_spec;
 
/* Vectors of the base of the kernel */
iy = 0;
  NOYAU_1D(nbeta,bmask,vecxy,idim_xx,iy);
 
dim_spec = nbeta + 1;

/* Alpha_0 en xx(.,1) */
for(iy = 0; iy < ny; iy++)
  {
  iy_s = iy * dim_spec;
  for(nb = 0; nb <= nbeta; nb++)
    {
    cc_re = xc_re[nb + iy_s];
    cc_im = xc_im[nb + iy_s];
/* X(NB+1,1) = BMASK[NB]*ATAN2( IMAG(CC), REAL(CC) ) */
    jlp_atan2c(&xxc,cc_im,cc_re);
    xx[nb] = bmask[nb + iy_s] * xxc;
   }
/* Alpha_0^+  en xx(.,1) */
    PROJ_EPLUS_1D(nbeta,1,1,xx,idim_xx,vecxy);
 
/* Corresponding initial phasor: */
  for(nb = 0; nb <= nbeta; nb++)
    {
/* ALPHA=X(NB+1,1) */
    alpha = xx[nb];
    xc_re[nb + iy_s] = bmask[nb + iy_s] * cos((double)alpha); 
    xc_im[nb + iy_s] = bmask[nb + iy_s] * sin((double)alpha); 
    }
/* End of loop on iy: */
  }
 
return(0);
}
/*******************************************************************
* NOYAU_1D
* All solutions which derive from a translation are possible
* hence the kernel is formed with (1,2,3,....,nbeta)
* (See TRANS_1D for justification)
*
* Define the vector of the base of the kernel
*  vecxy(idim_xx)
*******************************************************************/
int NOYAU_1D(INT4 nbeta, float *bmask, float *vecxy, INT4 idim_xx, int iy)
{
float xnorm;
INT4 nb, iix;
double sum;
 
/* The base vector is formed with the coordinates of the uv vectors. 
*/
for(nb = 0; nb <= nbeta; nb++)
  { 
/*  VECXY(NB,1)=BMASK(NB)*IXY(1,NB)
*  VECXY(NB,2)=BMASK(NB)*IXY(2,NB)
*/
    COVER_IXY_1D(&iix,&nb);
    vecxy[nb] = bmask[nb] * (float)iix;
  }

/* Normalization of the first vector: */
  NORMALIZE_L2(vecxy,nbeta,&xnorm,iy);

return(0);
}
/*******************************************************************
* PROJ_EPLUS_1D
* Projects X(.,N1) onto E^+ and stores result in X(.,N2)
* X(IDIM_X,4),VECXY(IDIM_X,2)
*******************************************************************/
int PROJ_EPLUS_1D(INT4 nbeta, INT4 n1, INT4 n2, float *xx, 
                  INT4 idim_xx, float *vecxy)
{
double xx1, delta;
INT4 nb, nn1, nn2;
 
/* Scalar product of X(.,N1) with the base vectors of the kernel */
xx1 = 0.;
nn1 = n1 - 1;
nn2 = n2 - 1;
 
for(nb = 1; nb <= nbeta; nb++)
  {
   xx1 += xx[nb + nn1 * idim_xx] * vecxy[nb];
  }

 delta = sqrt(xx1 * xx1);
/*
 printf("Distance to E+: DELTA = %17.10e \n",delta);
*/
 
/* Projection onto E^+ and stored in X(.,N2) */
for(nb = 1; nb <= nbeta; nb++)
   xx[nb + nn2 * idim_xx] = xx[nb + nn1 * idim_xx] - xx1 * vecxy[nb];

return(0);
}
/*******************************************************************
* SORTIE_1D : CREE LES FICHIERS DE SORTIE RES ET IMS
* Output: FFT along the columns
*******************************************************************/
int SORTIE_1D(float *xc_re, float *xc_im, INT4 nbeta, 
              INT4 isortie, char *fname, float *bmask)
{
float cc_re, cc_im, xr, xi, xm, *reversed_output;
char name[60], comments[81];
INT4 dim_spec, iy_s, nb, ixc, iix, ix, iy; 
register int i, j;

/* Assume that zero frequency is at IXC,IY:  */
   ixc = nx / 2;
 
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
dim_spec = nbeta + 1;
for(iy = 0; iy < ny; iy++)
{
  iy_s = iy * dim_spec;
  for(nb = 0; nb <= nbeta; nb++)
    {
/* JLP94, try...
*  XM = BMASK[NB]*RO[NB]
*/
     xm = ro[nb + iy_s];
     cc_re = xc_re[nb + iy_s];
     cc_im = xc_im[nb + iy_s];
     xr = xm * cc_re;
     xi = xm * cc_im;
     COVER_IXY_1D(&ix,&nb);
     iix = ixc + ix;
     re[iix + iy * nx] = xr;
     im[iix + iy * nx] = xi;
     iix = ixc - ix;
     re[iix + iy * nx] = xr;
     im[iix + iy * nx] = -xi;
   }
/* End of loop on iy: */ 
}

/* Output to image files: 
* In rei and imi if isortie=1
* In ref and imf if isortie=1
*/ 
  reversed_output = (float *)malloc(nx*ny*sizeof(float));
  if(isortie == 1)
    {
/* Initial output (recursive solution) */
    strcpy(name,"rei");
    sprintf(comments," FFT (rei) of : %.20s",fname);
    printf(" Output of %s \n",comments);
    for(j = 0; j < ny; j++)
      for(i = 0; i < nx; i++)
        reversed_output[j + i * ny] = re[i + j * nx];
    JLP_WRITEIMAG(reversed_output,&ny,&nx,&ny,name,comments);
/*
    JLP_WRITEIMAG(re,&nx,&ny,&nx,name,comments);
*/
    strcpy(name,"imi");
    sprintf(comments," FFT (imi) of : %.20s",fname);
    printf(" Output of %s \n",comments);
    for(j = 0; j < ny; j++)
      for(i = 0; i < nx; i++)
        reversed_output[j + i * ny] = im[i + j * nx];
    JLP_WRITEIMAG(reversed_output,&ny,&nx,&ny,name,comments);
/*
    JLP_WRITEIMAG(im,&nx,&ny,&nx,name,comments);
*/
    }
  else
    {
/* Final output (least squares solution) */
    strcpy(name,"ref");
    sprintf(comments," FFT (ref) of : %.20s",fname);
    printf(" Output of %s \n",comments);
    for(j = 0; j < ny; j++)
      for(i = 0; i < nx; i++)
        reversed_output[j + i * ny] = re[i + j * nx];
    JLP_WRITEIMAG(reversed_output,&ny,&nx,&ny,name,comments);
/*
    JLP_WRITEIMAG(re,&nx,&ny,&nx,name,comments);
*/

    strcpy(name,"imf");
    sprintf(comments," FFT (imf) of : %.20s",fname);
    printf(" Output of %s \n",comments);
    for(j = 0; j < ny; j++)
      for(i = 0; i < nx; i++)
        reversed_output[j + i * ny] = im[i + j * nx];
    JLP_WRITEIMAG(reversed_output,&ny,&nx,&ny,name,comments);
/*
    JLP_WRITEIMAG(im,&nx,&ny,&nx,name,comments);
*/
    }
  free(reversed_output);

return(0);
}
/**********************************************************************
* Main loop
* Least square minimization, non-linear least square fit
*
* xx(.,4)
**********************************************************************/
int LSQUARES1_1D(float exit_tolerance, INT4 nbeta, INT4 ngamma, INT4 ifermax,
                 float *bmask, float *yce_re, float *yce_im,
                 float *xc_re, float *xc_im,
                 float *weight, float *xx, INT4 idim_xx, INT4 ittmax)
{
float w1, sup_phi, qmoy, qmoy0, delq, cos1, sin1;
INT4 nb, itt, it;
 
 qmoy0 = 0.;
 
/* Main loop: */
for(itt = 1; itt <= ittmax; itt++)
  {
/* Computes phase error term X(.,1) for PHI^+
* Iterative solution by conjugate gradients.
*/ 
   CGRADIENT_1D(nbeta,ngamma,&it,&qmoy,ifermax,bmask,
                yce_re,yce_im,xc_re,xc_im,weight,xx,idim_xx);

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
   if (sup_phi < exit_tolerance)
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
                    1,1,&xc_re[nb],&xc_im[nb]); 
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
* CGRADIENT_1D 
* Resolution of the linear system: 
*                  [ATwA] X(.,1)  = X(.,3)
* with conjugate gradients
* X(.,1) is the unknown (and then solution) PHI 
* X(.,2) is the direction D_n
* X(.,3) is the second member of the equation, and then the residual
* X(.,4) is ATwA D_n   (called Z_n)
*******************************************************************/
int CGRADIENT_1D(INT4 nbeta, INT4 ngamma, INT4 *it, float *qmoy, 
                 INT4 ifermax, float *bmask, float *yce_re, float *yce_im, 
                 float *xc_re, float *xc_im,
                 float *weight, float *xx, INT4 idim_xx)
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
  ATRANS_1D(nbeta,ngamma,qmoy,ifermax,bmask,yce_re,yce_im,
            xc_re,xc_im,weight,xx,idim_xx);
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
   printf(" CGRADIENT_1D/Error: Square norm of Residual_{N}=0 !");
   return(-1);
   }
 
/* Main loop: */
 for(*it = 1; *it <= itmax; (*it)++)
  { 
/* STEP 1
* Compute Z_N = [AT A] D_N and store it in X(.,4)
*/
   ATRANS_A_1D(nbeta,ngamma,2,4,ifermax,bmask,weight,xx,idim_xx);

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
     printf(" CGRADIENT_1D/Error: Square norm of Delta_{PHI+1}=0 !\n");
     fprintf(fp1," CGRADIENT_1D/Error: Square norm of Delta_{PHI+1}=0 !\n");
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
* ATRANS_1D compute the second member of the system to be solved.
* This member is stored in X(.,3) which is the initial
* residual R_0
*
*  R_0 = AT ( PSI - A PHI_0)
* Here the residual is IMAG(exper.bispectrum * CONJ(estimated bispectrum))) :
*******************************************************************/
int ATRANS_1D(INT4 nbeta, INT4 ngamma, float *qmoy, INT4 ifermax, 
              float *bmask, float *yce_re, float *yce_im,
              float *xc_re, float *xc_im,
              float *weight, float *xx, INT4 idim_xx)
{
float c1_re, c1_im, c2_re, c2_im, c3_re, c3_im, yy1_re, yy1_im;
double q2;
INT4 ng1, ng2, nb, ng, indx, k, l, m;
 
/* Initialization to zero: */
for(nb = 1; nb <= nbeta; nb++)
   xx[nb + 2 * idim_xx] = 0.;
 
/* Weighted quadratic measure: */
   q2 = 0.;
 
 ng1 = 0;
 for(nb = 2; nb <= nbeta; nb++)
   {
/* NG2=NGT(NB) */
    cover_ngt1(&ng2,nb);
    indx = 0;
    for(ng = ng1; ng < ng2; ng++)
      {
      indx++;
      if(indx <= ifermax && bmask[nb] != 0.)
        {
/*  K=KLM(1,NG) L=KLM(2,NG) M=KLM(3,NG)
*/
         cover_klm0(&k,1,ng);
         cover_klm0(&l,2,ng);
         cover_klm0(&m,3,ng);

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
          q2 += (c3_re * c3_re + c3_im * c3_im) * weight[ng];
/* Residual is IMAG(YCE[NG]*CONJ(C1)) :
* JLP93: Put WEIGHT=g**2 for AT A and for A: 
*  yy1=imag(c2*conjg(c1))*weight[ng]
*/
          COMPLEX_PROD2(c2_re,c2_im,c1_re,c1_im,1,-1,&yy1_re,&yy1_im);
          yy1_im *= weight[ng];
/* Computing AT [YY1] by adding the contribution to X(.,3):
*/
          xx[k + 2 * idim_xx] += yy1_im;
          xx[l + 2 * idim_xx] += yy1_im;
          xx[m + 2 * idim_xx] -= yy1_im;
/* End of if(indx < ifermax): */
         }
/* End of loop on ng: */
      }
    ng1 = ng2;
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
* ATRANS_A_1D compute X(.,N2) = [AT A] X(.,N1)
*******************************************************************/
int ATRANS_A_1D(INT4 nbeta, INT4 ngamma, INT4 n1, INT4 n2, INT4 ifermax,
                float *bmask, float *weight, float *xx, INT4 idim_xx)
{
INT4 nn1, nn2, ng1, ng2, nb, ng, indx, k, l, m;
float yy1;
 
nn1 = n1 - 1;
nn2 = n2 - 1;

/* Initialization to zero: */
for(nb = 1; nb <= nbeta; nb++)
   xx[nb + nn2 * idim_xx] = 0.;
 
 ng1 = 0;
 for(nb = 2; nb <= nbeta; nb++)
   {
/* NG2=NGT(NB) */
    cover_ngt1(&ng2,nb);
    indx = 0;
    for(ng = ng1; ng < ng2; ng++)
      {
      indx++;
      if(indx <= ifermax && bmask[nb] != 0.)
        {
/*  K=KLM(1,NG) L=KLM(2,NG) M=KLM(3,NG)
*/
         cover_klm0(&k,1,ng);
         cover_klm0(&l,2,ng);
         cover_klm0(&m,3,ng);

/* JLP93: Put WEIGHT=g**2 for AT A */
         yy1 = (xx[k + nn1 * idim_xx] + xx[l + nn1 * idim_xx]
                - xx[m + nn1 * idim_xx]) * weight[ng];
/* Computing AT [YY1] by adding the contribution to X(.,N2): */
         xx[k + nn2 * idim_xx] += yy1;
         xx[l + nn2 * idim_xx] += yy1;
         xx[m + nn2 * idim_xx] -= yy1;
/* End of if(indx < ifermax): */
        }
/* End of loop on ng: */
      }
    ng1 = ng2;
/* End of loop on nb: */
  }
 
return(0);
}
/**********************************************************************
* FINAL_ERROR_1D 
* Compute errors when object is available (simulations)
*
**********************************************************************/
int FINAL_ERROR_1D(INT4 ir, INT4 nbeta, float *final_err_spec, 
                   float *final_err_phas, char *errname, float *bmask,
                   float *xc_re, float *xc_im, float *xcref_re, float *xcref_im,
                   float *xx, INT4 idim_xx, float *vecxy, float *sigm)
{
float cc_re, cc_im, w1, w2;
INT4 nb, irs, i, nbs;

 ERROR_SIMU_1D(xc_re,xc_im,xcref_re,xcref_im,nbeta,final_err_spec,
               final_err_phas,errname,bmask);
 
/* Ecart de phase angulaire point par point exprime en radian
*/ 
  for(nb = 0; nb <= nbeta; nb++)
    {
     cc_re = xcref_re[nb];
     cc_im = xcref_im[nb];
     jlp_atan2c(&w1,cc_im,cc_re);
     cc_re = xc_re[nb];
     cc_im = xc_im[nb];
     jlp_atan2c(&w2,cc_im,cc_re);
     xx[nb] = bmask[nb] * ( w1 - w2);
    }
 
/* Projection of this error X(.,1) onto E^+ */ 
  PROJ_EPLUS_1D(nbeta,1,1,xx,idim_xx,vecxy);
 
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
/************************************************************
* Set complex array XC to (1.,0.)
*
************************************************************/
int ZERO_PHASE_1D(float *xc_re, float *xc_im, INT4 nbeta)
{
register int iy;
int nb, iy_s, dim_spec;

dim_spec = nbeta + 1;

for(iy = 0; iy < ny; iy++)
{
iy_s = iy * dim_spec;
for(nb = 0; nb <= nbeta; nb++)
  {
   xc_re[iy_s + nb] = 1.;
   xc_im[iy_s + nb] = 0.;
  }
}

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
/*********************************************************************
* Sum all lines  of re and im and put result in these lines (debugging mode):
*
**********************************************************************/
int SUM_LINES_REIM()
{
register int iy, ix;
double sum_re, sum_im;

sum_re = 0.;
sum_im = 0.;
for(iy = 0; iy < ny; iy++)
   {
   for(ix = 0; ix < nx; ix++)
      {
      sum_re += re[ix + iy * nx];
      sum_im += im[ix + iy * nx];
      }
   }

/* Now fill the arrays with these sums: */
for(iy = 0; iy < ny; iy++)
   {
   for(ix = 0; ix < nx; ix++)
      {
      re[ix + iy * nx] = sum_re;
      im[ix + iy * nx] = sum_im;
      }
   }

return(0);
}
/*********************************************************************
* Write content of 0th line to all lines (for re and im)
*
**********************************************************************/
int EQUAL_LINES_REIM()
{
register int iy, ix;

for(iy = 0; iy < ny; iy++)
   {
   for(ix = 0; ix < nx; ix++)
      {
      re[ix + iy * nx] = re[ix];
      im[ix + iy * nx] = im[ix];
      }
   }

return(0);
}
