/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 decode_ima_1D.c
 To compute power spectrum and bispectrum of elementary frames
 1-D only: spectrum computed along the columns.
 Same as decode_ima but for 1-D spectra and 2-D bispectra.
 Allow for rotations of the initial spectrum.

 Possibility of apodization and zero-mean frame processing
(when #define APO_AND_ZEROMEAN)

According to external definition of FITS_CUBE, ERIC_FORMAT or nothing,
3 versions: decode_fits_1D_cube.exe, decode_ima_1D_eric.exe, decode_ima_1D.exe
(see decode.make)

 JLP
 Version 23-01-00
---------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <jlp_ftoc.h>
#include <jlp_fftw.h>

/*
#define APO_AND_ZEROMEAN
*/
/*
#define DEBUG 1
*/
#define MAXFILES 8200 

/* To allow for FITS cube format: 
#define FITS_CUBE 1 
(defined in decode.make file to generate decode_1D_fits_cube.exe)
*/
/* To allow for Eric Aristidi's format: 
#define ERIC_FORMAT 1 
(defined in decode.make file to generate decode_1D_aris.exe)
*/
/* To localize the photons: */
/*
#define PHOTOCOUNTING 1 
*/
#define PI  3.14159

main(argc, argv)
int argc;
char **argv;
{
/* Images nx*ny, irmax=25, ngmax=187566 : */

double sum1, mean1;
double *image, *modsq, *snrm, *long_int, *bispp, *mim;
float *float_array, *image1, *work;
float xframes, xphot;
INT_PNTR pntr;
/* Rotation matrix: */
float c1, c2, c3, c4, c5, c6;
int  fd1, ivalues, nvalues, nph;
int max_nclosure, ixcent, iycent;
int  isize, ir, nbeta, nxy, photon_correction;
int  nframes, iz_start, iz_end, istart, iend, jstart, jend;
INT4  nx1, ny1, nx1_good, ny1_good, nx, ny, nz, kod, ngamma, nx_bisp; 
INT4 istatus, descr_length;
int ifil, nfiles, idim, vm_flag, dflag; 
double rot_angle;
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
  printf(" Program decode_fits_1D_cube  Version 03-10-98 \n");
  printf(" To decode elementary frames\n");
#elif ERIC_FORMAT
  printf(" Program decode_ima_1D_eric  Version 03-10-98 \n");
  printf(" To decode elementary frames\n");
  printf(" Eric Aristidi's format\n");
#else
  printf(" Program decode_ima_1D  Version 24-01-2000 \n");
  printf(" To decode elementary frames  (JLP_FMT image formats)\n");
#endif
/**************************************************************/
/* Opening logfile: */
sprintf(logfile,"decode_im%s.log",file_ext);
if((fp1 = fopen(logfile,"w")) == NULL)
   {
   printf("decode_ima_1D/Fatal error opening logfile >%s< \n",logfile);
   exit(-1);
   }
  fprintf(fp1," Program decode_ima_1D  Version 23-01-2000 \n");
  fprintf(fp1," To decode elementary frames \n");

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
  printf("decode_ima_1D generic_name[or selection file] Radius_uv_coverage,nclosure,rotation_angle");
  printf(" output_file_extension ixstart,ixend,iystart,iyend [iz_start,iz_end] \n");
  printf(" Example1: decode_ima_1D pdm 12,100,30.5 _s1 1,32,1,32 1,50\n ");
  printf(" Example2  (with list of files in test.dat): ");
  printf("decode_ima_1D test.dat 15,200,19.5,mask_uv _s2 32,63,1,32 \n ");
  printf("\n(Possibility of adding uv-mask file name) \n");
  exit(-1);
  }

/* Radius of uv-coverage: */
/* Syntax: ir,max_nclosure,rotation_angle
 or:       ir,max_nclosure,rotation_angle,fmask_name 
*/ 
  *fmask_name = '\0';
  strcpy(buffer,argv[2]);
  sscanf(buffer,"%d,%d,%lf",&ir,&max_nclosure,&rot_angle);
/* Look for possible presence of mask name: */
  pc = buffer;
  while(*pc && *pc != ',') pc++;
  if(*pc == ',') pc++;
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
#ifdef DEBUG
   printf("buffer= >%.20s< \n",buffer);
#endif
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
    jlp0_rdfits(image1,float_array,&nx1,&ny1,&nz,&idim,
                 infile[0],comments,descr_value,&dflag,&istatus,&vm_flag);
    if(istatus)
     {
     printf(" Fatal error reading image %s \n",infile[0]);
     printf(" istatus = %d \n",istatus);
     exit(-1);
     }
    free(image1);
#ifdef DEBUG
     printf(" OK: nx1=%d ny1=%d nz=%d idim=%d\n",nx1, ny1, nz, idim);
     printf(" OK: nx=%d ny=%d\n",nx, ny);
#endif
#elif ERIC_FORMAT
 nx1 = 128;
 ny1 = 128;
 isize = nx1 * ny1 * sizeof(float);
 image1 = (float *)malloc(isize);
#endif

/* Allocation of memory: */
 isize = nx * ny * sizeof(double);
 long_int = (double *)malloc(isize);
 modsq = (double *)malloc(isize);
 snrm = (double *)malloc(isize);
 mim = (double *)malloc(isize);
 image = (double *)malloc(isize);

#ifdef PHOTOCOUNTING
 isize = nx * ny * sizeof(float);
 work = (float *)malloc(isize);
#endif

/* Erasing the arrays : */
    for(i = 0; i < nxy; i++) 
      {
      long_int[i]=0.;
      modsq[i]=0.;
      snrm[i]=0.;
      mim[i]=0.;
      }

/* Rotation if needed: */
  fprintf(fp1," Rotation angle: %f \n",rot_angle);
/* Conversion to radians: */
  rot_angle *= (PI/180.);
  c2 = cos(rot_angle);  c3 = -sin(rot_angle);
  c5 = sin(rot_angle); c6 = cos(rot_angle);
/* Compute c1 and c4 to maintain center in the middle: */
  ixcent = nx/2; iycent = ny/2;
  c1 = (float)ixcent * (1. - c2) - (float)iycent * c3;
  c4 = (float)iycent * (1. - c6) - (float)ixcent * c5;
  fprintf(fp1," Rotation matrix: c1=%f c2=%f c3=%f  \n",c1,c2,c3);
  fprintf(fp1,"                  c4=%f c5=%f c6=%f  \n",c4,c5,c6);

/* Computing the uv coverage: */
  printf(" Radius of uv-coverage (IR) in pixels: %d\n",ir);
  compute_uv_coverage_1D(fp1,fmask_name,ir,&nbeta,&ngamma,max_nclosure);
/* For debugging purpose: */
  output_lists_coverage_1D(&nbeta,&ngamma);

/* Allocating memory space for the bispectrum: */
  isize = 4 * ngamma * nx * sizeof(double);
  bispp= (double *)malloc(isize);
  for(i=0; i < 4 * ngamma * nx; i++) { bispp[i]=0.;}

#ifdef APO_AND_ZEROMEAN
/* Creates Hamming apodization file: */
 isize = nx * ny * sizeof(float);
 apodi = (float *)malloc(isize);
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
  for(ifil = 0; ifil < nfiles; ifil++)
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
    dflag = 0; vm_flag = 1;
/* dflag set to -1 to prevent error warning messages: */
    if(nframes > 1) dflag = -1;
    jlp0_rdfits(image1,float_array,&nx1,&ny1,&ifil,&idim,
                 infile[0],comments,descr_value,&dflag,&istatus,&vm_flag);
    if(istatus)
     {
     printf(" Fatal error reading image #%d \n",nframes);
     printf(" istatus = %d \n",istatus);
     nframes--;
     goto end1;
     }
#elif ERIC_FORMAT
  ERIC_READIMAG1(image1,&nx1,&ny1,infile[ifil]);
#else
  JLP_VM_READIMAG1(&pntr,&nx1,&ny1,infile[ifil],comments);
  image1 = (float *)pntr;
#endif


/* To get photon detection: */
#ifdef PHOTOCOUNTING
  seuillage_old(image1,nx1,threshold,iside);
/*
  seuillage(image1,work,nx1,threshold,iside);
*/
#endif

/* At first iteration, sets nx1_good and ny1_good: */
  if(!ifil) {nx1_good = nx1; ny1_good = ny1;}

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


/* Rotation */
  rotate_ima(image1, nx1, ny1, image, nx, ny, 
             istart, jstart, c1, c2, c3, c4, c5, c6, rot_angle);

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
    sprintf(" Ouput of %s \n",outfile);
    JLP_D_WRITEIMAG(image,&nx,&ny,&nx,outfile,outcomments);
   }
#endif

/* Output index loop every 100 frames: */
    if((nframes % 100) == 1)
       printf(" Processing frame # %d / %d\n",nframes,nfiles);

/* Long integration : */
  for(i = 0; i < nxy; i++) long_int[i] += image[i];

/* JLP99 */
/* Resetting the imaginary part for the FFT: */
  for(i = 0; i < nxy;i++) {mim[i] = 0.;}

/* Fourrier Transform: */
  kod = 1;
  FFT_1D_Y(image,mim,&nx,&ny,&nx,&kod);

/* bispec_1D_Y in "jlp_cover_mask.c" */
   bispec_1D_Y(image,mim,modsq,snrm,&nx,&ny,bispp,&ir,&nbeta,&ngamma);

/****************** End of loop with ifil (all input files processed) ****/
#ifndef ERIC_FORMAT
   free(image1);
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

   printf(" Calling prepare output bisp_1D \n");
   prepare_output_bisp_1D(bispp,modsq,snrm,nx,ny,xframes,xphot,
                          nbeta,ngamma,photon_correction,fp1);

/* Mean of the long integration: */
   for(i = 0; i < nxy; i++) long_int[i] /= xframes;

/*********************************************************/
/* Now output of the results : */

/* Comments: */
#ifdef FITS_CUBE
      sprintf(outcomments,"%s nfr: %d, from %d to %d, nfiles= %d",
           infile[0],nframes,iz_start,iz_end,nfiles);
#else
   if(argc == 5)
/* Case of a list: */
      sprintf(outcomments,"%s nfiles=%d",infile[0],nfiles);
   else
      sprintf(outcomments,"%s nfr: %d, from %d to %d, nfiles= %d",
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
   RECENT_FFT_1D_Y(modsq,modsq,&nx,&ny,&nx);
   sprintf(outfile,"modsq%s",file_ext);
   JLP_D_WRITEIMAG(modsq,&nx,&ny,&nx,outfile,outcomments);

/* SNR of squared modulus (actually 1/sigma if no photon correction): */
   RECENT_FFT_1D_Y(snrm,snrm,&nx,&ny,&nx);
   sprintf(outfile,"snrm%s",file_ext);
   JLP_D_WRITEIMAG(snrm,&nx,&ny,&nx,outfile,outcomments);

/* JLP99 */
   printf(" calling output bisp_1D \n");
/* Bispectrum : */
   nx_bisp = 3 * ngamma;
   isize = nx_bisp * nx * sizeof(float);
   float_array = (float *)malloc(isize);
   output_bisp_1D(bispp,float_array,ngamma,nx);
#ifdef DEBUG_2
  printf(" New sequence is: re[i],re[i+1]...im[i],im[i+1]...snr[i],snr[i+1]\n");
  for(i = 0; i < 5; i++)
      printf("bispp[%d]: (%e;%e)\n",i,bispp[i],bispp[i+ngamma]);
#endif
   for(i=0; i < 3; i++)
     printf("bispp[%d]=%f, float_array[%d]=%f\n",i,bispp[i],i,float_array[i]);
   free(bispp);
   sprintf(outfile,"bisp1%s",file_ext);
   JLP_WRITEIMAG(float_array,&nx_bisp,&nx,&nx_bisp,outfile,outcomments);
   free(float_array);

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
int ERIC_READIMAG1(tab, nx, ny, filename)
float *tab;
char *filename;
int *nx, *ny;
{
 int dim;
 register int i;
 FILE *fic,*fopen();

 *nx = 128; *ny = 128;
 dim = *nx * *ny;

 if ((fic = fopen(filename,"r")) != NULL)
 {
  for (i=0; i<dim; i++) *tab++ = (float) getc(fic);
  fclose(fic);
 }
 else 
 {
  printf("ERIC_READIMAG1/Fatal error opening file >%s< \n",filename);
  exit(-1);
 }

/* tab[i + j * dim] = 
Attention, inversion en Y par rapport aux fichiers TIFF
*/
}
/****************************************************************
Salut,

Bruno Lopez est en train de rediger un article general sur Mira  
regrouppant les observations de poussiere par optique adaptative  
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

int rotate_ima(im1, nx1, ny1, im2, nx, ny, istart, jstart, 
                c1, c2, c3, c4, c5, c6, rot_angle)
float *im1;
double *im2, rot_angle;
int nx, ny, nx1, ny1, istart, jstart;
float c1, c2, c3, c4, c5, c6;
{
register int i1, j1, i_rot, j_rot;

/* Case of invariance: */
if(istart == 0 && jstart == 0 && rot_angle == 0.)
 {
 for(j1 = 0; j1 < ny1; j1++)
  for(i1 = 0; i1 < nx1; i1++)
              im2[i1 + j1 * nx] = (double)im1[i1 + j1 * nx1];
 }
 else
 {
/* Add 0.5 to avoid steps... */
/* WARNING: should use new variables, since i is modified in first line... */
 for(j1 = 0; j1 < ny; j1++)
  {
  for(i1 = 0; i1 < nx; i1++)
    {
    i_rot = (int)(0.5 + c1 + c2 * (float)i1 + c3 * (float)j1);
    j_rot = (int)(0.5 + c4 + c5 * (float)i1 + c6 * (float)j1);
/* Transfer to image array, using working window: */
    i_rot -= istart;
    j_rot -= jstart;
    if(i_rot > 0 && i_rot < nx1 && j_rot > 0 && j_rot < ny1)
              im2[i1 + j1 * nx] = (double)im1[i_rot + j_rot * nx1];
    }
   }
  }
return(0);
}
