/********************************************************************
* bisp_to_image.c
*
* Conversion of a part of the bispectrum to an image file 
* for further display 
*
* bisp(u,v) = spec(u_1,0) * spec(0,v_2) * spec(-u_1,-v_2) 
* Hence we display the image bisp(u_1,0,0,v_2), 
* which is a function of the two coordinates u_1 and v_2
*
* Calls: 
* covera_mask, bisp_to_2D_image 
* (in jlp_cover_mask.c)
*  
* JLP
* Version 09-06-2008
*******************************************************************/
#include <stdio.h>
#include <math.h>
#include <jlp_ftoc.h>
#include <jlp_cover_mask.h>

#define DEBUG

int main(int argc, char *argv[])
{
float *mask, *modsq, *bisp, *out_image;
INT4 nx_mask, ny_mask, ir, max_nclosure, nbeta, ngamma;
INT4 nx_modsq, ny_modsq, nx_bisp, ny_bisp, nx_out, ny_out; 
INT_PNTR pntr1;
INT4 ShouldNormalize, isize, iopt, bispectrum_is_1D, line_to_output;
register int i;
char mask_name[61], mask_comments[81], buffer[81];
char modsq_name[61], modsq_comments[81];
char bisp_name[61], bisp_comments[81];
char out_name[61], out_comments[140], line_string[20];

printf(" bisp_to_image           -- Version 03-09-2008 --\n");
printf(" Conversion of square modulus and bispectrum to a 2-D image file \n");

/* Seven parameters are needed to run the program: */
/* Carefull: 7 parameters always, using JLP "runs" */
if(argc != 1 && !strcmp(argv[6],""))
  {
  printf("        Fatal error: Wrong syntax, argc=%d\n",argc);
  printf(" Syntax is:  \n");
  printf(" runs bisp_to_image n 0 12,200 modsq_11b bisp1_11b test");
  printf(" Is 1-Dimensionnal? n (or y)\n");
  printf(" 0=real part 1=imag part 2=snr 3=phase 4=modulus \n");
  printf(" ir,nclosure  in_modulus in_bispectrum output\n");
 }

/* 1D or 2D ? */
bispectrum_is_1D = 0;
printf(" Is the bispectrum 1-Dimensionnal (slit spectrum) ? (y)\n");
if(argc == 1) gets(buffer); else strcpy(buffer,argv[1]);
if(buffer[0] != 'N' && buffer[0] !='n') 
  {
  bispectrum_is_1D = 1;
  printf(" OK, 1D bispectrum\n");
  printf(" Let (u_1,0) and (v_1,0) the two vectors u and v\n");
  printf(" We create here the image bisp(u_1,0,v_1,0) for a given line\n");
  }
else
  {
  printf(" OK, 2D bispectrum\n");
  printf("    Let (u_1,u_2) and (v_1,v_2) the two vectors u and v\n");
  printf("    We create here the image bisp(u_1,0,0,v_2)\n");
  printf("    which is a function of the two coordinates u_1 and v_2 \n");
  }

/* Main menu: */
printf(" Menu: \n");
printf(" 0 = Real part of bispectrum\n");
printf(" 1 = Imaginary part of bispectrum\n");
printf(" 2 = SNR of bispectrum\n");
printf(" 3 = Phase of bispectrum\n");
printf(" 4 = Modulus of bispectrum\n");
printf(" Enter your choice =< ");
if(argc == 1) gets(buffer); else strcpy(buffer,argv[2]);
sscanf(buffer,"%d",&iopt);
if(iopt < 0 || iopt > 4)
  {
   printf(" Fatal error: option should be between 0 and 4 (here you have selected %d)\n",iopt);
   return(-1);
  }
else
  printf(" OK: option=%d\n",iopt);

/***************************************************************/
JLP_INQUIFMT();

  printf(" Radius of uv-coverage (IR) in pixels and max_nclosure:\n"); 
  if(argc == 1) gets(buffer); else strcpy(buffer,argv[3]);
  sscanf(buffer,"%d,%d",&ir,&max_nclosure);
#ifdef DEBUG
  printf(" OK, ir=%d max_nclosure=%d \n",ir,max_nclosure);
#endif

/* Read mask file: */
/* JLP96: commented out:
  printf(" Input mask (u-v coverage) (0 if no mask) =< "); gets(mask_name);
*/
  strcpy(mask_name,"0");
  if(mask_name[0] == '0')
    {
      nx_mask = 2 * (ir + 1);
      if(bispectrum_is_1D) ny_mask = 1;
            else ny_mask = nx_mask;
      isize = nx_mask * ny_mask * sizeof(float);
      JLP_GVM(&mask,&isize);

/* Create a filled mask: */
      for(i = 0; i < nx_mask * ny_mask; i++) mask[i] = 1.;
    }
  else
    {
     JLP_VM_READIMAG1(&pntr1,&nx_mask,&ny_mask,mask_name,mask_comments);
     mask = (float*)pntr1;
     if(bispectrum_is_1D && ny_mask != 1)
        {printf(" Fatal error: mask is not 1D-image!\n"); return(-1);}
    }

/* "Masked" u-v coverage: */
  if(bispectrum_is_1D)
/* Compute u-v coverage with the mask and within a segment of half width ir: */
    COVERA_MASK_1D(mask,&nx_mask,&ir,&max_nclosure,&nbeta,&ngamma);
  else
/* Compute u-v coverage with the mask and within a disk of radius ir: */
    COVERA_MASK(mask,&nx_mask,&ny_mask,&ir,&max_nclosure,&nbeta,&ngamma);

/* Free memory: */
  JLP_FVM(&mask);

/* Read square modulus file */
  if(iopt == 2)
    printf(" Input SNR of square modulus (snrm) =< "); 
  else
    printf(" Input square modulus (modsq) =< "); 
  if(argc == 1) gets(modsq_name); else strcpy(modsq_name,argv[4]);
  JLP_VM_READIMAG1(&pntr1,&nx_modsq,&ny_modsq,modsq_name,modsq_comments);
  modsq = (float *)pntr1;

/* Check size is correct: */
  if(nx_modsq < 2 * ir || ny_modsq < 2 * ir)
    {
    printf("Fatal error: modulus file is too small, does not fit u-v coverage\n");
    return(-1);
    } 

/* Read bispectrum file */
  printf(" Input bispectrum =< "); 
  if(argc == 1) gets(bisp_name); else strcpy(bisp_name,argv[5]);
  JLP_VM_READIMAG1(&pntr1,&nx_bisp,&ny_bisp,bisp_name,bisp_comments);
  bisp = (float *)pntr1;

/* Check size is correct: (JLP96: I add the possibility of ny=3 and nx=ngamma */
  if(nx_bisp !=  3*ngamma && nx_bisp != ngamma)
    {
    printf("Fatal error: bispectrum size inconsistent with the current u-v coverage\n");
    printf(" nx_bisp = %d, whereas 3*ngamma = %d \n",nx_bisp,3*ngamma);
    return(-1);
    } 

/* Out image name: */
  printf(" Output image name =< "); 
  if(argc == 1) gets(out_name); else strcpy(out_name,argv[6]);

  if(bispectrum_is_1D)
    {
    printf(" Number of the line you want to output: ");
    if(argc == 1) gets(buffer); else strcpy(buffer,argv[7]);
    sscanf(buffer,"%d",&line_to_output);
    if(line_to_output < 0 || line_to_output >= ny_modsq)
      {
      printf("Error: line_to_output=%d !\n",line_to_output); 
      line_to_output = ny_modsq/2;
      printf("I take default value: line_to_output=%d\n",line_to_output);
      }
    sprintf(line_string,"line#%d ",line_to_output);
    nx_out = 2 * (nbeta + 1);
    ny_out = nx_out;
    }
  else
    {
    nx_out = nx_modsq;
    ny_out = ny_modsq;
    strcpy(line_string," ");
    }

/* Write comments and check if iopt is correct, otherwise set it to "real": */
switch (iopt)
   {
   case 0:
     sprintf(out_comments,"%sReal from %s and %s",
                 line_string,modsq_name,bisp_name);
     break;
   case 1:
     sprintf(out_comments,"%sImag. from %s and %s",
                 line_string,modsq_name,bisp_name);
     break;
   case 2:
     sprintf(out_comments,"%sSNR from %s and %s",
                 line_string,modsq_name,bisp_name);
     break;
   case 3:
     sprintf(out_comments,"%sPhase from %s",
                 line_string,bisp_name);
     break;
   case 4:
     sprintf(out_comments,"%sModulus from %s",
                 line_string,bisp_name);
     break;
   default:
     sprintf(out_comments,"%sReal from %s and %s",
                 line_string,modsq_name,bisp_name);
     printf("Incorrect option, I set it to default (i.e., real part)\n");
     iopt = 0;
     break;
   }
/* Get memory space: */
   isize = nx_out * ny_out * sizeof(float);
   JLP_GVM(&out_image,&isize);

/* Clear initial image: */
   for(i = 0; i < nx_out * ny_out; i++) out_image[i] = 0.;

/* Now compute bispectrum slice: 
* (if iopt=0, real part, if iopt=1 imag part, if iopt=2, SNR, if iopt=3, Phase,
* if iopt=4 Modulus)
*/
ShouldNormalize = 0;
   if(bispectrum_is_1D)
      {
      bisp1D_to_2D_image(out_image,&nx_out,&ny_out,modsq,&nx_modsq,
                        bisp,&nx_bisp,&nbeta,&ngamma,&line_to_output,
                        &ShouldNormalize,&iopt);
      }
   else
      {
      bisp2D_to_2D_image(out_image,modsq,&nx_modsq,&ny_modsq,
                    bisp,&nx_bisp,&nbeta,&ngamma,&ShouldNormalize,&iopt);
      }

/* Output image: */
  JLP_WRITEIMAG(out_image,&nx_out,&ny_out,&nx_out,out_name,out_comments);

/* Free memory: */
  JLP_FVM(&modsq);
  JLP_FVM(&bisp);

return(0);
}
