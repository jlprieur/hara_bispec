/********************************************************************
* Conversion of bispectrum and modulus from JLP to Eric Anterrieu's format 
*
* Calls: 
* cover_mask_to_eric, covera_mask modsq_to_eric, bisp_to_eric 
* (in jlp_cover_mask.c)
*  
* JLP
* Version 22-11-93
*******************************************************************/
#include <stdio.h>
#include <math.h>
#include <jlp_ftoc.h>

#define DEBUG

main()
{
float *mask, *modsq, *bisp;
int nx_mask, ny_mask, ir, max_nclosure, nbeta, ngamma;
int nx_modsq, ny_modsq, nx_bisp, ny_bisp; 
int pntr1, isize;
register int i;
char mask_name[61], mask_comments[81], buffer[81];
char modsq_name[61], modsq_comments[81];
char bisp_name[61], bisp_comments[81];
char bisp_list[61], spec_list[61], dimension_file[61];
char bisp_ascii_file[61], mod_ascii_file[61], out_ext[41], *pc;

printf(" bisp_to_eric           -- Version 22-11-93 --\n");
printf(" Conversion of square modulus and bispectrum to Eric's format \n");
printf(" To combine many sets of masked data, run \"merge_to_eric\" afterwards\n");

  printf(" Output extension (Creturn if only one set of data):\n"); 
  gets(out_ext);
/* Replaces the first blank character encountered by zero (end of string): */
  pc = out_ext;
  while(*pc != ' ' && *pc) pc++;
  *pc = '\0'; 

/***************************************************************/
JLP_INQUIFMT();

  printf(" Radius of uv-coverage (IR) in pixels and max_nclosure:\n"); 
  gets(buffer); sscanf(buffer,"%d,%d",&ir,&max_nclosure);
#ifdef DEBUG
  printf(" OK, ir=%d max_nclosure=%d \n",ir,max_nclosure);
#endif

/* Read mask file: */
  printf(" Input mask (u-v coverage) (0 if no mask) =< "); gets(mask_name);
  if(mask_name[0] == '0')
    {
      nx_mask = 2 * (ir + 1);
      ny_mask = nx_mask;
      isize = nx_mask * ny_mask * sizeof(float);
      JLP_GVM(&mask,&isize);

/* Create a filled mask: */
      for(i = 0; i < nx_mask * ny_mask; i++) mask[i] = 1.;
    }
  else
    {
     JLP_VM_READIMAG(&pntr1,&nx_mask,&ny_mask,mask_name,mask_comments);
     JLP_FROM_MADRID(&pntr1,&mask);
    }

/* Compute u-v coverage with the mask and within a disk of radius ir: */
  COVERA_MASK(mask,&nx_mask,&ny_mask,&ir,&max_nclosure,&nbeta,&ngamma);

/* Free memory: */
  JLP_FVM(&mask);

/* Output spectral, bispectral list in Eric's format */
sprintf(bisp_list,"list_bispectre%s.txt",out_ext);
sprintf(spec_list,"list_spectre%s.txt",out_ext);
sprintf(dimension_file,"dimensions%s.txt",out_ext);
cover_mask_to_eric(bisp_list,spec_list,dimension_file,&nbeta,&ngamma);

/* Read square modulus file */
  printf(" Input square modulus =< "); gets(modsq_name);
  JLP_VM_READIMAG(&pntr1,&nx_modsq,&ny_modsq,modsq_name,modsq_comments);
  JLP_FROM_MADRID(&pntr1,&modsq);

/* Check size is correct: */
  if(nx_modsq < 2 * ir || ny_modsq < 2 * ir)
    {
    printf("Fatal error: modulus file is too small, does not fit u-v coverage\n");
    exit(-1);
    } 

/* Now output to Eric's format ASCII file: */
   sprintf(mod_ascii_file,"module_spectre_exp%s.txt",out_ext);
   modsq_to_eric(modsq,&nx_modsq,&ny_modsq,mod_ascii_file,&nbeta);

/* Free memory: */
  JLP_FVM(&modsq);

/* Read bispectrum file */
  printf(" Input bispectrum =< "); gets(bisp_name);
  JLP_VM_READIMAG(&pntr1,&nx_bisp,&ny_bisp,bisp_name,bisp_comments);
  JLP_FROM_MADRID(&pntr1,&bisp);

/* Check size is correct: */
  if(nx_bisp != ngamma)
    {
    printf("Fatal error: bispectrum size inconsistent with the current u-v coverage\n");
    exit(-1);
    } 

/* Now output to Eric's format ASCII file: */
   sprintf(bisp_ascii_file,"phaseur_bispectre_exp%s.txt",out_ext);
   bisp_to_eric(bisp,&nx_bisp,bisp_ascii_file,&ngamma);

/* Free memory: */
  JLP_FVM(&bisp);

exit(0);
}
