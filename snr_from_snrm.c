/**************************************************************** 
*  snr_from_snrm 
*  Infer the SNR and SIG maps from the *_snrm.fits file created 
*  by vcrb with the "bispectrum option"
*
* JLP
* Version 09-04-2008
*****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <jlp_ftoc.h>

#define DEBUG

/* Define minimum value of sigma, in order to avoid
* problems with zeroes and small values of sigma:
*/
#define MIN_SIGMA 1.e-2

int main(int argc, char *argv[])
{
/* Images : */
float *modsq1, *snrm1, *snr1, *sig1, ss;
INT_PNTR pntr_ima;
INT4 nx, ny;
register int i;
char *pc, in_ext[41], out_ext[41];
char modsq_filename[61], snrm_filename[61];
char filename[61], comments[81];

printf(" Program snr_from_snrm \n");
printf(" JLP Version 10-04-08 \n");
printf(" For the extension, either _tt or tt_\n");
printf("(i.e. either modsq_tt, bisp1_tt, ... or tt_m, tt_b, t_l\n");


/* One or three parameters only are allowed to run the program: */
/* Carefull: 7 parameters always, using JLP "runs" */
#ifdef DEBUG
  printf(" argc=%d\n",argc);
  printf(" argv[2]=>%s<\n",argv[2]);
#endif
if((argc != 1) && (argc != 3) && (argc == 7 && !strcmp(argv[2],"")))
  {
  printf("        Fatal error: Wrong syntax, argc=%d\n",argc);
  printf(" Syntax is:  \n");
  printf("runs snr_from_snrm in_exten out_exten\n");
  exit(-1);
  }

/* Interactive input of in_exten and out_exten: */
if (argc == 1)
 {
   printf(" Input file extension := "); scanf("%s",in_ext);
/************* File extension for output files: */
   printf(" Output file extension := "); scanf("%s",out_ext);
 }
else
 {
  strcpy(in_ext,argv[1]);
  strcpy(out_ext,argv[2]);
 }

/* End extension strings with zero: */
  pc = in_ext;
  while(*pc && *pc != ' ') pc++;
  *pc='\0';
  pc = out_ext;
  while(*pc && *pc != ' ') pc++;
  *pc='\0';

#ifdef DEBUG
printf(" DEBUG Version, will read >%s< files and write >%s< \n",in_ext,out_ext);
#endif

/***************************************************************/
JLP_INQUIFMT();

/* JLP99: I add the possibility of using tt_l, tt_m, tt_b
instead of long_tt modsq_tt bisp1_tt */

/* Mean modsq file : */
   if(in_ext[0] == '_')sprintf(modsq_filename,"m%s",in_ext);
   else sprintf(modsq_filename,"%sm",in_ext);
   printf(" Reading modsq file: >%s< \n",modsq_filename);
   JLP_VM_READIMAG1(&pntr_ima,&nx,&ny,modsq_filename,comments);
   modsq1 = (float *) pntr_ima;

/* Mean snrm file : */
   if(in_ext[0] == '_')sprintf(snrm_filename,"snrm%s",in_ext);
   else sprintf(snrm_filename,"%ssnrm",in_ext);
   printf(" Reading snrm file: >%s< \n",snrm_filename);
   JLP_VM_READIMAG1(&pntr_ima,&nx,&ny,snrm_filename,comments);
   snrm1 = (float *) pntr_ima;

/* Allocation of memory: */
  sig1 = (float *)malloc(nx * ny * sizeof(float));
  snr1 = (float *)malloc(nx * ny * sizeof(float));

/* Computation of sig1 and snr1: */
for(i = 0; i < nx * ny; i++) sig1[i] = snrm1[i];
for(i = 0; i < nx * ny; i++) {
   ss = MAXI(sig1[i], MIN_SIGMA);
   snr1[i] = snrm1[i] * modsq1[i]/ss;
   }

/* SIG file */ 
   sprintf(filename,"%s.SIG",out_ext);
   sprintf(comments,"From %s", modsq_filename);
   JLP_WRITEIMAG(sig1,&nx,&ny,&nx,filename,comments);

/* SNR file */ 
   sprintf(filename,"%s.SNR",out_ext);
   sprintf(comments,"From %s and %s", modsq_filename, snrm_filename);
   JLP_WRITEIMAG(snr1,&nx,&ny,&nx,filename,comments);

free(sig1);
free(snr1);

return(0);
}
