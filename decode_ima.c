/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 decode_ima.c
 To compute power spectrum and bispectrum of elementary frames

 Possibility of apodization and zero-mean frame processing
(when #define APO_AND_ZEROMEAN)

Accoding to external definition of FITS_CUBE, ERIC_FORMAT or nothing,
3 versions: decode_fits_cube.exe, decode_ima_eric.exe, decode_ima.exe
(see decode.make)

 JLP
 Version 13-10-99
---------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <jlp_ftoc.h>

/*
#define APO_AND_ZEROMEAN
*/
/*
#define DEBUG 1
*/
/* refp: 17850 files! */
#define MAXFILES 18000 

/* To allow for FITS cube format: 
#define FITS_CUBE 1 
(defined in decode.make file to generate decode_fits_cube.exe)
*/
/* To allow for Eric Aristidi's format: 
#define ERIC_FORMAT 1 
(defined in decode.make file to generate decode_aris.exe)
*/
/* To localize the photons: */
/*
#define PHOTOCOUNTING 1 
*/

int bispec3(double *re, double *im, double *modsq, double *snrm,
            INT4 *nx, INT4 *ny, double *bispp, INT4 *ir, 
            INT4 *nbeta, INT4 *ngamma);
int ERIC_READIMAG1(float *tab, INT4 *nx, INT4 *ny, char *filename);

main(argc, argv)
int argc;
char **argv;
{
/* Images nx*ny, irmax=25, ngmax=187566 : */

double sum1, mean1;
double *image, *modsq, *snrm, *long_int, *bispp, *mim;
float *float_array, *pntr_image, *work;
float xframes, xphot;
INT_PNTR pntr_ima, pntr_cube;
INT4  fd1, ivalues, nvalues, pntr, nph;
INT4 max_nclosure, descr_length;
INT4  isize, ir, nbeta, ngamma, photon_correction;
INT4  nframes, iz_start, iz_end, istatus, istart, iend, jstart, jend;
INT4  nx1, ny1, nx1_good, ny1_good, nx, ny, nz, nxy, ifil, nfiles, iw;
INT4 idim, vm_flag, dflag, ifil_start, ifil_end, kod; 
register int i, j ,i1, j1;
char infile[MAXFILES][61], generic_name[61], comments[81], buffer[81];
char file_list[41], descr_name[61], descr_value[81], logfile[81], *pc;
char file_ext[61], outfile[41], outcomments[81], fmask_name[61];
FILE *fp, *fp1;
#ifdef PHOTOCOUNTING
float threshold;
int iside;
#endif

#ifdef APO_AND_ZEROMEAN
float *apodi;
#endif

#ifdef FITS_CUBE 
  printf(" Program decode_fits_cube  Version 13-10-99 \n");
  printf(" To decode elementary frames\n");
#elif ERIC_FORMAT
  printf(" Program decode_ima_eric  Version 03-02-98 \n");
  printf(" To decode elementary frames\n");
  printf(" Eric Aristidi's format\n");
#else
  printf(" Program decode_ima  Version 03-02-98 \n");
  printf(" To decode elementary frames  (JLP_FMT image formats)\n");
#endif
/**************************************************************/
/* Opening logfile: */
sprintf(logfile,"decode_im%s.log",file_ext);
if((fp1 = fopen(logfile,"w")) == NULL)
   {
   printf("decode_ima/Fatal error opening logfile >%s< \n",logfile);
   exit(-1);
   }
  fprintf(fp1," Program decode_ima  Version 22-09-98 \n");
  fprintf(fp1," To decode elementary frames \n");


/***************************************************************
* Hamming apodization file: routine defined in "decode_set.c"
**************************************************************/
#ifdef APO_AND_ZEROMEAN
  printf(" Apodization with Hamming filter and zero-mean processing\n");
  printf("  (Previously, error in Hamming... nx -> nx/2))\n");
#endif
  printf("  (Should be like:  t0001.fits, t0002.fits, etc)\n");
  printf("  Maximum number of frames: %d\n",MAXFILES);

/* Input parameters:  can be 4,5 or 6*/
if (argc != 5 && argc != 6 )
  {
  printf("argc = %d \n",argc);
  printf("\nUSAGE:\n");
  printf("decode_ima generic_name[or selection file] Radius_uv_coverage,nclosure");
  printf(" output_file_extension ixstart,ixend,iystart,iyend [iz_start,iz_end] \n");
  printf(" Example1: decode_ima pdm 12,100 _s1 1,32,1,32 1,50\n ");
  printf(" Example2  (with list of files in test.dat): ");
  printf("decode_ima test.dat 15,200,mask_uv _s2 32,63,1,32 \n ");
  printf("\n(Possibility of adding uv-mask file name) \n");
  exit(-1);
  }

/* Radius of uv-coverage: */
/* Syntax: ir,max_nclosure
 or:       ir,max_nclosure,fmask_name 
*/ 
  *fmask_name = '\0';
  strcpy(buffer,argv[2]);
  sscanf(buffer,"%d,%d",&ir,&max_nclosure);
/* Look for possible presence of mask name: */
  pc = buffer;
  while(*pc && *pc != ',') pc++;
  if(*pc == ',') pc++;
  while(*pc && *pc != ',') pc++;
  if(*pc == ',') {pc++; strcpy(fmask_name, pc);}
  if(*fmask_name) printf(" Mask of uv-coverage : %s\n",fmask_name);
  else printf(" Full pupil was used: no mask\n");

/* Working window: */
/*
  sscanf(argv[4],"%d,%d,%d,%d,%f",&istart,&iend,&jstart,&jend,&threshold);
*/
  sscanf(argv[4],"%d,%d,%d,%d",&istart,&iend,&jstart,&jend);
  nx = iend - istart + 1;
  ny = jend - jstart + 1;
  nxy = nx * ny;
  if(nx <= 1 || ny <= 1)
   {printf(" Fatal error, sorry, wrong window nx=%d ny=%d (!)\n", nx, ny);
    exit(-1);
   }

/* Output file extension: */
  strcpy(file_ext,argv[3]);
  pc = file_ext;
  while(*pc && *pc != ' ') pc++;
  *pc='\0';

#ifdef FITS_CUBE
/* Generic name for the data image files: */
   strcpy(infile[0],argv[1]);

/* iz_start, iz_end (first and last frames): */
   sscanf(argv[5],"%d,%d",&iz_start,&iz_end);
   nfiles = iz_end-iz_start+1;
   if(nfiles > MAXFILES) 
     {
     printf("WARNING: %d files are wanted whereas maximum number is %d files\n",
             nfiles,MAXFILES);
     nfiles = MAXFILES; iz_end = iz_start + nfiles - 1;
     }
   printf("iz_start=%d, iz_end=%d \n",iz_start,iz_end);
#else
/* Case of list of files (Example 2): ****************************************/
 if(argc == 5)
  {
  strcpy(file_list,argv[1]);
  if((fp = fopen(file_list,"r")) == NULL)
    {printf(" Fatal error, cannot open file >%s< \n",file_list);
     exit(-1);
    }

/* Reading the list, with a loop on the lines: */
   nfiles = -1;
   while(fgets(buffer,81,fp) != NULL)
   {
   nfiles++;
   if(nfiles > MAXFILES) 
     {
     printf("WARNING: %d files are wanted whereas maximum number is %d files\n",
            nfiles,MAXFILES);
     nfiles = MAXFILES; break;
     }
   sscanf(buffer,"%s",infile[nfiles]);
   }

   nfiles++;
   fclose(fp);
  }
/* Generic name: (Example 1) *************************************************/
 else
  {
/* Generic name for the data image files: */
   strcpy(generic_name,argv[1]);

/* iz_start, iz_end (first and last frames): */
   sscanf(argv[5],"%d,%d",&iz_start,&iz_end);
   nfiles = iz_end-iz_start+1;
   if(nfiles > MAXFILES) 
     {
     printf("WARNING: %d files are wanted whereas maximum number is %d files\n",
             nfiles,MAXFILES);
     nfiles = MAXFILES; iz_end = iz_start + nfiles - 1;
     }

/* Generates the file names: */
   for(i=iz_start; i<=iz_end; i++)
    sprintf(infile[i-iz_start],"%s%02d",generic_name,i);

  } /* end of (argc != 5) */
/* End of .not. FITS_CUBE: */
#endif

/* Good values (checked experimentally, JLP feb98) */
#ifdef PHOTOCOUNTING
  threshold = 20;
  iside = 2;
  printf("PHOTOCOUNTING/Threshold = %f, iside = %d (without smoothing) \n",threshold,iside);
#endif

/*****************************************************************/
/* Begining of the main work: */
   JLP_INQUIFMT();

#ifdef FITS_CUBE
    dflag = 0; vm_flag = 1; nz = 0;
    JLP_VM_RDFITS_3D(&pntr_cube,&nx1,&ny1,&nz,infile[0],comments,descr_value,
                     &dflag,&istatus);
    pntr_image = (float *)pntr_cube;
    if(istatus)
     {
     printf(" (FITS_CUBE)Fatal error reading image %s \n",infile[0]);
     printf(" istatus = %d \n",istatus);
     exit(-1);
     }
     idim = nx1;
     printf(" nx1=%d ny1=%d nz=%d\n",nx1, ny1, nz);
     if(iz_end > nz) iz_end = nz;
#ifdef DEBUG
     printf(" OK: nx=%d ny=%d\n",nx, ny);
#endif
#elif ERIC_FORMAT
 nx1 = 128;
 ny1 = 128;
 isize = nx1 * ny1 * sizeof(float);
 pntr_image = (float *) malloc(isize);
#endif

/* Allocation of memory: */
 isize = nx * ny * sizeof(double);
 long_int = (double *) malloc(isize);
 modsq = (double *) malloc(isize);
 snrm = (double *) malloc(isize);
 image = (double *) malloc(isize);
 mim = (double *) malloc(isize);
 if(long_int == NULL || modsq == NULL || snrm == NULL || image == NULL
    || mim == NULL)
 {
 printf(" Fatal error allocating memory space (long_int, ...) \n");
 exit(-1);
 }

#ifdef PHOTOCOUNTING
 isize = nx * ny * sizeof(float);
 work = (float *) malloc(isize);
#endif

/* Erasing the arrays : */
    for(i = 0; i < nxy; i++) 
      {
      long_int[i]=0.;
      modsq[i]=0.;
      snrm[i]=0.;
      mim[i]=0.;
      }

/* Computing the uv coverage: */
  printf(" Radius of uv-coverage (IR) in pixels: %d\n",ir);
  compute_uv_coverage(fp1,fmask_name,ir,&nbeta,&ngamma,max_nclosure);

/* Allocating memory space for the bispectrum: */
  isize = 4 * ngamma * sizeof(double);
  bispp = (double*) malloc(isize);
  for(i=0; i < 4 * ngamma; i++) { bispp[i]=0.;}

#ifdef APO_AND_ZEROMEAN
/* Creates Hamming apodization file: */
 isize = nx * ny * sizeof(float);
 apodi = (float *) malloc(isize);
  jlp_hamming(apodi,nx,ny,nx);

#ifdef DEBUG
   strcpy(outfile,"hamming");
   strcpy(outcomments,"Hamming filter: 0.54 + 0.46 cos(pi t / T)");
   JLP_WRITEIMAG(apodi,&nx,&ny,&nx,outfile,outcomments);
#endif
#endif

/***************************************************************/
/* Main loop on the frames: */
  nframes = 0;
#ifdef FITS_CUBE
  ifil_start = iz_start;
  ifil_end = iz_end+1;
#else
  ifil_start = 0;
  ifil_end = nfiles;
#endif

  for(ifil = ifil_start; ifil < ifil_end; ifil++)
  {
#ifdef FITS_CUBE
#ifdef DEBUG
  printf(" Processing image #%d \n",ifil);
#endif
  fprintf(fp1," Processing image #%d \n",ifil);
#else
#ifdef DEBUG
  printf(" Processing %s \n",infile[ifil]);
#endif
  fprintf(fp1," Processing %s \n",infile[ifil]);
#endif

/* Reading the data file */
  nframes++; 
#ifdef FITS_CUBE
/* JLP99 */
  pntr_image = (float *)pntr_cube;
  pntr_image = &pntr_image[nx1*ny1*(ifil-1)];
#elif ERIC_FORMAT
  istatus = ERIC_READIMAG1(pntr_image,&nx1,&ny1,infile[ifil]);
#else
  istatus = JLP_VM_READIMAG1(&pntr_ima,&nx1,&ny1,infile[ifil],comments);
  if(istatus)
     {
     printf(" Fatal error reading image %s\n",infile[ifil]);
     exit(-1);
     }
  pntr_image = (float *)pntr_ima;
#ifdef DEBUG
  printf(" Successful reading of image %s\n",infile[ifil]);
#endif
#endif

/* To get photon detection: */
#ifdef PHOTOCOUNTING
  seuillage_old(pntr_image,nx1,threshold,iside);
/*
  seuillage(pntr_image,work,nx1,threshold,iside);
*/
#endif

/* At first iteration, sets nx1_good and ny1_good: */
  if(ifil == ifil_start) {nx1_good = nx1; ny1_good = ny1;}

/* Check size: */
  if(nx1 != nx1_good || ny1 != ny1_good)
     {
     printf(" Error: image size incompatible with previous images! \n");
     fprintf(fp1," Error: image size incompatible with previous images! \n");
     istatus = 1;
     }

/* Error handling: */
  if(istatus)
     {
     printf(" Fatal error reading image #%d \n",nframes);
     fprintf(fp1," Fatal error reading image #%d \n",nframes);
     printf(" istatus = %d \n",istatus);
     nframes--;
     goto end1;
     }

/* Transfer to image array, using working window: */
  for(j=0; j<ny; j++) 
    {
    j1 = j + jstart;
    for(i=0; i<nx; i++) 
      {
      i1 = i + istart;
      image[i + j * nx] = (double)pntr_image[i1 + j1 * nx1];
      }
    }

#ifdef APO_AND_ZEROMEAN
/* Compute mean and sum: */
  sum1 = 0.;
  for(i=0; i < nxy; i++) sum1 += image[i]; 
  mean1 = sum1 / (double)nxy;

#ifdef DEBUG
 printf(" sum = %f mean = %f \n",sum1,mean1);
#endif

/* Subtraction of mean to get zero-mean frame: */
  for(i=0; i < nxy; i++) image[i] -= mean1; 

/* Normalization to unity: */
  if(sum1 != 0) for(i=0; i < nxy; i++) image[i] /= sum1; 

/* Multiply with Hamming filter: */
  for(i=0; i < nxy; i++) image[i] *= apodi[i]; 
#endif

/* Output of the first image as a diagnostic: */
#ifdef DEBUG
  if(nframes == 1)
   {
    strcpy(outfile,"first_image");
    strcpy(outcomments,"   ");
    printf(" Ouput of %s (nx=%d, ny=%d, istart=%d, jstart=%d)\n",
            outfile,nx,ny,istart,jstart);
    JLP_D_WRITEIMAG(image,&nx,&ny,&nx,outfile,outcomments);
   }
#endif

/* Output index loop every 100 frames: */
    if((nframes % 100) == 1)
       printf(" Processing frame # %d / %d\n",nframes,nfiles);

/* Long integration : */
  for(i = 0; i < nxy; i++) long_int[i] += image[i];

/* Resetting the imaginary part for the FFT: */
  for(i = 0; i < nxy;i++) {mim[i] = 0.;}

/* Fourrier Transform: */
  kod = 1;
  FFT_2D_DOUBLE(image, mim, &nx, &ny, &nx, &kod);

/* Processing this image now:  bispec1 is with photon noise correction
bispec3 is without*/
  /*
   bispec1(image,mim,modsq,snrm,&nx,&ny,&nx,bispp,&ir,&nbeta,&ngamma);
   */
/* bispec3 in "jlp_cover_mask.c" */
   bispec3(image,mim,modsq,snrm,&nx,&ny,bispp,&ir,&nbeta,&ngamma);

/****************** End of loop with ifil (all input files processed) ****/
#ifndef FITS_CUBE
#ifndef ERIC_FORMAT
   free(pntr_image);
#endif
#endif
   }

/* Check that some frames have been selected */
   xframes=(float)nframes;
   if(!nframes) 
     {printf("    Sorry no frame has been processed \n");
     goto end1;
     }
   else
     printf("    %d frames have been processed \n",nframes);

/* Compute snrm, SNR of bispectrum, etc:  (xphot is not used) */
   xphot = 1.; photon_correction = 0;

#ifdef DEBUG_2
  printf(" Now sequence is: re[i],im[i],sumsq_re[i],sumsq_im[i],re[i+1],im[i+1]...\n");
  for(i = 0; i < 5; i++)
      printf("bispp[%d]: (%e;%e)\n",i,bispp[4*i],bispp[4*i+1]);
#endif

   prepare_output_bisp(bispp,modsq,snrm,nx,ny,xframes,xphot,
                         nbeta,ngamma,photon_correction,fp1);

#ifdef DEBUG_2  
  printf(" New sequence is: re[i],re[i+1]...im[i],im[i+1]...sumsq_re[i],sumsq_re[i+1]...sumsq_im[i],sumsq_im[i+1]\n");
  for(i = 0; i < 5; i++)
      printf("bispp[%d]: (%e;%e)\n",i,bispp[i],bispp[i+ngamma]);
#endif


/* Mean of the frames: */
   for(i = 0; i < nxy; i++) {
	long_int[i] /= xframes;
	}
/*********************************************************/
/* Now output of the results : */

/* Comments: */
#ifdef FITS_CUBE
   sprintf(outcomments,"%s nfr=%d from %d to %d, nfiles=%d",
           infile[0],nframes,iz_start,iz_end,nfiles);
#else
   if(argc == 5)
/* Case of a list: */
      sprintf(outcomments,"%s nfiles=%d",infile[0],nfiles);
   else
      sprintf(outcomments,"%s nfr=%d from %d to %d, nfiles=%d",
           infile[0],nframes,iz_start,iz_end,nfiles);
#endif

/* Long integration : */
  sprintf(outfile,"long%s",file_ext);
  fprintf(fp1," %s",outfile);
  JLP_D_WRITEIMAG(long_int,&nx,&ny,&nx,outfile,outcomments);

/* Leave blanks at the end, for Fortran interface: */
   sprintf(descr_name,"IR  ");
   sprintf(descr_value,"%d  ",ir);
   descr_length = 4;
   JLP_WDESCR(descr_name,descr_value,&descr_length,&istatus);
/* Leave blanks at the end, for Fortran interface: */
   sprintf(descr_name,"MAX_NCLO  ");
   sprintf(descr_value,"%d  ",max_nclosure);
   descr_length = 8;
   JLP_WDESCR(descr_name,descr_value,&descr_length,&istatus);

/* Mean squared modulus : */
   RECENT_FFT_DOUBLE(modsq,modsq,&nx,&ny,&nx);
   sprintf(outfile,"modsq%s",file_ext);
   JLP_D_WRITEIMAG(modsq,&nx,&ny,&nx,outfile,outcomments);

/* SNR of squared modulus (actually 1/sigma if no photon correction): */
   RECENT_FFT_DOUBLE(snrm,snrm,&nx,&ny,&nx);
   sprintf(outfile,"snrm%s",file_ext);
   JLP_D_WRITEIMAG(snrm,&nx,&ny,&nx,outfile,outcomments);

/* Bispectrum : */
/*
   isize = 3 * ngamma * sizeof(float);
   float_array = (float *) malloc(isize);
   if(float_array == NULL)
     {
     printf("Fatal error allocating memory space for float_array, isize=%d\n",
             isize);
     exit(-1);
     }
   output_bisp(bispp,float_array,ngamma);
#ifdef DEBUG_2
  printf(" New sequence is: re[i],re[i+1]...im[i],im[i+1]...snr[i],snr[i+1]\n");
  for(i = 0; i < 5; i++)
      printf("bispp[%d]: (%e;%e)\n",i,bispp[i],bispp[i+ngamma]);
#endif

   printf("bispp[0]=%f, float_array[0]=%f\n",bispp[0],float_array[0]);
   free(bispp);
   sprintf(outfile,"bisp1%s",file_ext);
   ny=3;
   printf(" OK ?\n");
   JLP_WRITEIMAG(float_array,&ngamma,&ny,&ngamma,outfile,outcomments);
   printf(" OK2 ?\n");
   free(float_array);
*/
   iw=3;
   printf(" ngamma = %d, ny = %d \n", ngamma, ny);
   sprintf(outfile,"bisp1%s",file_ext);
   JLP_D_WRITEIMAG(bispp,&ngamma,&iw,&ngamma,outfile,outcomments);

  printf(" Successful end \n");

/* End : */
  end1:close(fd1);
  free(mim);
#ifdef PHOTOCOUNTING
  free(work);
#endif

#ifdef APO_AND_ZEROMEAN
  free(apodi);
#endif

  JLP_END();
  fclose(fp1);
}
/* -----------------------------------
* Lecture des fichiers non formattes
*
* Eric Aristidi format
* 128x128 pixels
* ----------------------------------*/
int ERIC_READIMAG1(float *tab, INT4 *nx, INT4 *ny, char *filename)
{
 int dim;
 register int i;
 FILE *fp,*fopen();

 *nx = 128; *ny = 128;
 dim = *nx * *ny;

 if ((fp = fopen(filename,"r")) != NULL)
 {
  for (i=0; i<dim; i++) tab[i] = (float) getc(fp);
  fclose(fp);
 }
 else 
 {
  printf("ERIC_READIMAG1/Fatal error opening file >%s< \n",filename);
  exit(-1);
 }

/* tab[i + j * (*nx)] = 
Attention, inversion en Y par rapport aux fichiers TIFF
*/

/* Threshold of 10 for Mira's data (JLP: September 2000)*/
  for (i=0; i<dim; i++) 
     {
     if(tab[i] <= 10.) tab[i] = 0.;
       else tab[i] -= 10.;
     }

return(0);
}
/****************************************************************
Salut,

Bruno Lopez est en train de rediger un article general sur Mira  
regroupant les observations de poussiere par optique adaptative  
qu'il a faites a l'ESO, nos observations du compagnon et une  
modelisation de l'enveloppe. Tu es co-auteur.
A ce titre je te forwarde le mail qu'il m'a passe :
> Salut Eric,
>
> On a repris avec Pierre C. la redaction
> de l'article sur o Ceti.
>
> Peux tu te charger de pousser Jean-Louis Prieur a nous fabriquer
> une image a partir des donnees du Pic

l'idee est la suivante : deux articles en preparation. (en plus de  
celui sur les binaires avec Marco). Un est celui dont je parle au  
debut de ce mail, porte sur l'astrophysique de l'etoile et tres  
succint quant aux conditions d'observations et aux techniques de  
depouillement. Ce qu'il faudrait c'est juste une image pas trop  
pourrie. L'autre papier regroupperait l'ensemble des observations de  
miras doubles (parmi celles du catalogue Hipparcos, 4 ont ete  
detectees doubles) et les variations photometriques du compagnon  
(marcel est en train de travailler a adapter sa technique  
probabiliste au comptage de photons).

Je viens donc te remettre la pression sur Mira (ca tombe bien,  
celle sur l'article des doubles est en ce moment sur Marco). J'avais  
fait un bout de code C qui "dirac-ise" les photons pour simplifier  
les problemes de soustraction de biais. Peux-etre peux tu appliquer  
cette diracisation et ensuite tes programmes de traitement sur la  
premiere serie d'images (celle du debut de la sequence ou le  
compagnon est le mieux visible).

Je te joins le bout de code que j'avais ecrit :
ca seuille d'abord l'image a dark+3 sigma (ca fait 20 env.)
puis ca detecte les photons, ca calcule un barycentre et ca pose
un Dirac dont l'amplitude vaut le max de la "tache" du photon dans
l'image.

avec :

dim : dimension de l'image
ima : float*, image avec la convention
      pixel i,j -> *(ima+i*dim+j)

retour de ima dans le meme tableau.
hope it helps...

Eric

-------

*/
seuillage_old(ima,dim,threshold,iside)
int	dim, iside;
float	*ima, threshold;
{
/* iside=2, threshold=20, valeurs qui vont bien */
	int i,j,imax,jmax;
	int ii,jj;
	float pmax;
	
	for (i=iside-1; i<dim-iside; i++) 
         for (j=iside-1; j<dim-iside; j++) 
          if(*(ima+i*dim+j) < threshold) *(ima+i*dim+j)=0;
/* Centrage */
          else
	  {
/* Look for the maximum: */
		pmax=-1;
		for (ii=0; ii<=iside; ii++)
		for (jj=0; jj<=iside; jj++)
		if ( pmax < *(ima+(i+ii)*dim+j+jj) )
		{
			pmax = *(ima+(i+ii)*dim+j+jj);
			imax=i+ii;
			jmax=j+jj;
		};

		for (ii=-1; ii<iside; ii++)	/* on fait le vide */
		for (jj=-1; jj<iside; jj++)
		*(ima+(imax+ii)*dim+jmax+jj)=0;

/* We put a photon at this location:
JLP: to allow for correct photon correction, I tried to set it to 1.0
--> worse!
*/
		*(ima+imax*dim+jmax)=pmax;
		
	  }

/* gestion des bords */
	for (i=dim-iside; i<dim; i++) for (j=0; j<dim; j++) *(ima+i*dim+j)=0;
	for (i=0; i<dim; i++) for (j=dim-iside; j<dim; j++) *(ima+i*dim+j)=0;
	for (i=0; i<iside-1; i++) for (j=0; j<dim; j++)  *(ima+i*dim+j)=0; 
	for (i=0; i<dim; i++) for (j=0; j<iside-1; j++)  *(ima+i*dim+j)=0; 
	
return(0);
}
/********************************************************
* New version with data transfer (avoid setting to zero around maximum
* which introduces artefacts...)
* But actually worse!
* Tried also with smoothing: not better.
********************************************************/
seuillage(ima,work,dim,threshold,iside)
int	dim, iside;
float	*ima, *work, threshold;
{
/* iside=2, threshold=20, valeurs qui vont bien */
	int i,j,imax,jmax;
	int ii,jj;
	float pmax;

/* Smoothing the input image: */
	for (i=0; i<dim-1; i++) 
           {
           for (j=0; j<dim-1; j++) 
                {
                *(work + i*dim + j ) = (*(ima+i*dim+j) + *(ima+(i+1)*dim+j)
                             + *(ima+(i+1)*dim + j+1) + *(ima+i*dim+j+1))/4.; 
                }
           }
/* Filling edges: */
	for (i=0; i<dim; i++) *(work+i*dim+dim-1)=0.; 
	for (j=0; j<dim; j++) *(work+(dim-1)*dim+j)=0.; 

/* Erasing image: */
	for (i=0; i<dim; i++) for (j=0; j<dim; j++) *(ima+i*dim+j)=0.; 

	for (i=iside-1; i<dim-iside; i++) 
         for (j=iside-1; j<dim-iside; j++) 
          if(*(work+i*dim+j) > threshold)
/* Look for the maximum: */
	  {
		pmax=-1;
		for (ii=0; ii<=iside; ii++)
		for (jj=0; jj<=iside; jj++)
		if ( pmax < *(work+(i+ii)*dim+j+jj) )
		{
			pmax = *(work+(i+ii)*dim+j+jj);
			imax=i+ii;
			jmax=j+jj;
		};

/* We set to zero pixels around location: */
		for (ii=-1; ii<iside; ii++)
		for (jj=-1; jj<iside; jj++)
		*(ima+(imax+ii)*dim+jmax+jj)=0;

/* We put a photon at this location:
JLP: to allow for correct photon correction, I tried to set it to 1.0
--> worse!
*/
		*(ima+imax*dim+jmax)=pmax;
		
	  }

return(0);
}

