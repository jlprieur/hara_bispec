/******************************************************************
 autoco.c

 Routines to perform integration and autocorrelation
 from photon coordinates

Assumes:
  photon coord. on  9 + 9 bits.
  resolution reduced to 8 + 8 bits before computing the ac,
  The ac is in an 512xVERT array of short integers (i.e. 2 byte integers).

From LK (mars 1989)
JLP
Version 20-08-93
********************************************************************/

typedef photon long;
#define HORIZ 512
#define VERT 512

#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <jlp_ftoc.h>

/*
#define DEBUG
*/

#define MAXPHOTONS 12000 
#define IDIM 256
/* Note that DIM = IDIM * IDIM     */
#define DIM 65536 
/* Note that DIM2 = (2*IDIM + 1) * (2*IDIM + 1)     */
#define DIM2 263169

main()
{
/* Images IDIM*IDIM, irmax=25, ngmax=187566 : */

float out_auto[DIM], image[DIM], autocor[DIM2], long_int[DIM];
int xwork[MAXPHOTONS], ywork[MAXPHOTONS];
float sum, w1, w2;
FILE *fd1;
int ivalues, nvalues, nph, ixy, idim=IDIM;
unsigned char pacframe[MAXPHOTONS*2];
int i1, ix, iy, ixc, iyc, block_size, imax, kmax;
int ir, nbeta, ngamma, ngmax=NGMAX, in_descr, npack, nphotons;
int nshift, nframes, nx, ny, iframe, iph, istat;
auto int i, j, k;
unsigned length1;
char infile[41], nphfile[41], nphcomments[81], buffer[81];
char outfile[41], outcomments[81], cdescr[1025];

  printf(" Program decode_mam  Version 04-11-91 (%d x %d) \n",idim,idim);

/* Erasing the arrays : */
    for(i=0; i<DIM; i++) {
      long_int[i]=0.;
      out_auto[i]=0.;
      }
    for(i=0; i<DIM2; i++) autocor[i]=0.;

/* Reading the file with the number of photons */
  JLP_INQUIFMT();
  printf(" Total number of photons to be processed?\n");
  gets(buffer);sscanf(buffer,"%d",&nphotons);
  printf(" Number of photons per pack ?\n");
  gets(buffer);sscanf(buffer,"%d",&npack);
  npack = (npack > MAXPHOTONS - 10) ? MAXPHOTONS - 10 : npack;
  nx = 256; ny = 256;
  printf(" Will process %d photons in packs of %d photons \n",nphotons,npack);
  printf(" Output size is nx = %d, ny=%d \n",nx,ny);

/* Opening the data file */
  printf(" Input data file ?\n");
  gets(infile);
  printf(" Reading %s \n",infile);
  fd1=fopen(infile,"r");

/* Reading the first four-byte of the data file (integration time)*/
  nvalues = fread(pacframe,1,4,fd1);
  i = pacframe[0];
  j = pacframe[1];
  i = (i << 8) | j;
  j = pacframe[2];
  i = (i << 8) | j;
  j = pacframe[3];
  i = (i << 8) | j;
  printf(" Integration time is %d \n",i);

  printf(" Current pointer is %d \n",ftell(fd1));
  printf(" Number of bytes to shift?\n");
  gets(buffer);sscanf(buffer,"%d",&nshift);
  if(nshift == 1) nvalues = fread(pacframe,1,1,fd1);
  printf(" Current pointer is %d \n",ftell(fd1));

/* Reading the first 2*(npack-1) data values in the data file */
  nvalues = fread(pacframe,1,2*(npack-1),fd1);

#ifdef DEBUG
  printf(" pacframe[0] %d \n",pacframe[0]);
  printf(" Current pointer is %d \n",ftell(fd1));
#endif

  if(nvalues != 2*(npack-1)) 
   {printf(" Fatal error reading data file: nvalues=%d \n",nvalues);
    goto end1;}

/* Transfer to work array: */
   for(i=0; i<=npack-2; i++)
     {ywork[i] = pacframe[i*2+1];
     xwork[i] = pacframe[i*2];
     ixy = IDIM*ywork[i] + xwork[i];
#ifdef DEBUG
     printf(" pacframe[%d]= %d pacframe[%d]= %d\n",2*i,pacframe[2*i],2*i+1,
	      pacframe[2*i+1]);
     printf(" x, y %d %d \n",xwork[i],ywork[i]);
#endif
/* Long integration : */
#ifdef DEBUG
     if(ixy < 0 || ixy >= DIM)printf(" error: ixy = %d \n",ixy);
#endif
     long_int[ixy]++;}
  
/****************************************************************/
/* Start main loop: */
  iph = npack-1;
  while (iph < nphotons) {

/* Read by sets of (MAXPHOTONS-npack)*2 */
    block_size = (nphotons - iph > (MAXPHOTONS - npack)) ? 
		 (MAXPHOTONS - npack) : nphotons - iph;

    nvalues = fread(pacframe,1,2*block_size,fd1);
    if(nvalues == 0) 
     {printf(" End of data file: nphotons=%d \n",iph);
      break;}

    printf("iph=%d, block_size=%d nvalues=%d \n",iph,block_size,nvalues);

/* Main loop on all the values of packframe */
    kmax = nvalues/2;
    for(k=0; k < kmax; k++)
      {
      iph++;
/* Fills in last value of work arrays: */
      i1 = npack-1+k;
      xwork[i1] = pacframe[2*k];
      ywork[i1] = pacframe[2*k+1];
      ixy = IDIM*ywork[i1] + xwork[i1];

#ifdef DEBUG
    printf(" x, y %d %d \n",xwork[i1],ywork[i1]);
#endif

/* Long integration : */
      long_int[ixy]++;

/* Auto-correlation: */
      auto_cor(xwork,ywork,autocor,npack,k,nx,ny);

      }

/* Shift the array down by kmax pixels: */
      imax = npack+kmax-1;
      for(i = kmax; i < imax; i++) {
      i1=i-kmax;
      xwork[i1] = xwork[i]; 
      ywork[i1] = ywork[i]; 
      }

   }

   nphotons = iph;

/****************************************************************/
/* Now output the results : */

/* Autocorrelation : */
   symmetry_autocor(out_auto,autocor,nx,ny);
   strcpy(outfile,"autocor.bdf                               ");
   sprintf(outcomments," nph: %d; npack: %d ",nphotons,npack);
   JLP_WRITEIMAG(out_auto,&nx,&ny,&idim,outfile,outcomments);

/* Long integration : */

#ifdef DEBUG
   sum =0;
   for(i=0; i<DIM; i++) sum=sum+long_int[i];
   printf(" sum = %f \n",sum);
#endif

   strcpy(outfile,"long.bdf                                ");
   sprintf(outcomments," nph: %d; npack: %d ",nphotons,npack);
   JLP_WRITEIMAG(long_int,&nx,&ny,&idim,outfile,outcomments);

  printf(" Successful end \n");

/* End : */
  end1: fclose(fd1);
  JLP_END();
}
/*******************************************************
*  Output format of auto-correlation: 
*  Sum (x,y) and (-x,-y), and then symmetry (x,y) (-x,-y) 
*******************************************************/
int symmetry_autocor(out_auto,autocor,nx,ny)
float out_auto[],autocor[];
int nx, ny;
{
int idim2, ixc, iyc, ioff, joff;
auto int i, j;
/* Symmetry of autocorrelation 
 Pixel (i,j) is linked with (idim2-i,idim2-j)
*/
   idim2 = 2*IDIM + 1;
   ixc = nx;
   iyc = ny;
   for(j=1; j<=iyc; j++)
   {
    for(i=1; i<idim2; i++)
    {
    autocor[(idim2 - i) + (idim2 - j)*idim2] = autocor[i + j*idim2] 
	      + autocor[(idim2 - i) + (idim2 - j)*idim2];
    autocor[i + j*idim2] = autocor[(idim2 - i) + (idim2 - j)*idim2];
    }
   }
/* Troncation now: */
   ioff = nx/2;
   joff = ny/2;
   for(j=0; j<ny; j++)
   {
    for(i=0; i<nx; i++)
    {
     out_auto[i+j*IDIM]=autocor[(i+ioff)+(j+joff)*idim2];
    }
   }
return(0);
}
/**********************************************************************/
* Long integration 
* phot_ptr:     address of current photon
* ph_per_integ: number of photons per integration
* integ:        address of tab but (sum)
***********************************************************************/
void integ_from_buffer9 (phot_ptr, ph_per_integ, integ)
register photon	*phot_ptr;
register long	ph_per_integ;
register unsigned short *integ;
{
register short	*im_ptr;
register unsigned long	y, x, yx;

if (ph_per_integ < 1L) return;

while ( --ph_per_integ >= 0L) 
  {
/* Compute coordinates of pixel to increment: */ 
    im_ptr = (short *)( (long)integ + ((*phot_ptr) << 1) );
    (*im_ptr) ++;
    phot_ptr++;
  }
}

/* -------- integration d'autocorrelations en flux continu -------- */
void correlate9 (phot_ptr, ph_per_correl, ph_per_integ, ac_org)
photon	*phot_ptr;	/* pointer to current photon */
long	ph_per_correl;	/* photons per correlation time */
long	ph_per_integ;	/* photons per integration time */
unsigned short *ac_org; /* pointer to origin af target array (storing AC) */

/* integre l'ac  2-dim a partir de coordonnees de photons en memoire.
 renvoie le pointeur de photons mis a jour.
 L'ac est calculee par soustraction de coordonnee d'un photon a n autres photons.
On tire avantage du fait que :
1)  x & y sont cote a cote en memoire : on les traite en une seule operation.
Pour que la retenue eventuelle d'une soustraction des X n'affecte pas
les Y, il faut qu'il y ait un bit "vide" separant X de Y.
Ceci est fait au prix d'une baisse de resolution, en mettant a 0
le lsb de X et de Y, par la routine half_res9.
Maintenir la resolution initiale serait possible en doublant la
taille en X et en Y de la fenetre d'ac calculee.
L'interet de cette facon de faire est qu'il n'y a qu'une soustraction
et aucun test dans la boucle de calcul.
2) Le resultat d'une soustraction de coordonnees correspond directement 
a l'adresse d'un pixel dans le tableau contenant l'autocorrelation, a condition
que le nombre de pixels par ligne soit un puissance de 2.
L. Koechlin, Center for Astrophycs, Cambridge Mass, 1988
*/
{
register short	*ac_ptr;	/* points current pixel in ac */
long		integ_count;	/* photons left to integrate */
register long	correl_count;	/* photons left to correlate with current photon */
register unsigned long	xiyi;	/* coordonnee couplee d'un photon : xi & yi */
unsigned long		xi;	/* coordonnee seule */
unsigned long		yi;
register photon	*phot_ptrj;	/* will point ph.to correlate with current photon */
register long	ac = (long)ac_org;

integ_count = ph_per_integ;
if (integ_count < 1L)
	return;

while (--integ_count >= 0L)
{
  xiyi = *phot_ptr++;
/* ----- attention, ici on fait l'hypothese que xi est code sur 9 bits ----- */
  xi = xiyi & 0x1FF;	/* xi va de 0 a 510 par 2  (0 a HORIZ-2) */
  xi |= 0x200;	/* maintenant xi va de 512 a 1022 par 2 (HORIZ a 2*HORIZ-2) */
/* HORIZ doit imperativement etre une puissance de 2 */

  yi = xiyi >>9; /* yi va de 0 a VERT par 2 */
  yi += VERT;	 /* maintenant yi va de VERT a 2*VERT-2 par 2 */

  xiyi = xi | (yi<<9);	/* on remet ensemble xi et yi modifies */
/* ----- fin de la partie dependante du nb de bits de codage des x ----- */
		
  phot_ptrj = phot_ptr;	/* pointe sur xj & yj */
  correl_count = ph_per_correl;
  while (--correl_count >= 0L) {
  ac_ptr = (short *)(ac + xiyi - *phot_ptrj);
/* xj va de 0 a HORIZ-2 par pas de 2
yj va de 0 a VERT-2 par pas de 2.
	
xi-xj va de 0 a 2*HORIZ-2 par pas de 2.
yi-yj va de 0 a 2*VERT-2 par pas de 2.

Si HORIZ est une puissance de 2,
yi-yj correspond aux bits de poids forts
de l'adresse du pixel a incrementer dans l'ac,
et xi-xj aux bits de poids faibles.
	
Si l'ac est in tableau d'entiers sur 2 octets,
deux pixels consecutifs  ont une diff. d'adresse de 2,
comme les yi-yj et xi-xj sont par pas de 2, 
il n'y a pas de recalage a faire.

Ajouter a l'adresse de debut de tableau ac pour
trouver l'adresse du pixel a incrementer dans l'ac. */
  (*ac_ptr)++;			/* incrementer pixel but dans l'ac */
  phot_ptrj++;			/* xj & yj suivants */
  }
 }
}	

/* --------  integ d'autocorrelations image par image -------- */
void correlim9 (phot_ptri, ph_per_im, ac_org)

photon		*phot_ptri;	/* pointeur 1er photon de l'im*/
long		ph_per_im;	/* photons par image */
unsigned short	*ac_org;	/* pointeur d'origine du tableau but (AC) */
{
register short	*ac_ptr;	/* points current pixel in ac */
long		integ_count;	/* photons left to integrate */
register long	correl_count;	/* photons left to correlate with current photon */
register unsigned long	xiyi;	/* coordonnee couplee d'un photon : xi & yi */
unsigned long	xi;		/* coordonnee seule */
unsigned long	yi;
register photon	*phot_ptrj;	/* will point ph.to correlate with current photon */
register long	ac = (long)ac_org;

integ_count = ph_per_im;
if (integ_count < 1L)
    return;

while (--integ_count >= 0L)
{
   xiyi = *phot_ptri++;
/* ----- attention, ici on fait l'hypothese que xi est code sur 9 bits ----- */
   xi = xiyi & 0x1FF;	/* xi va de 0 a 510 par 2  (0  HORIZ-2) */
   xi |= 0x200;		/* maintenant xi va de 512 a 1022 par 2 (HORIZ  2*HORIZ-2) */
/* HORIZ doit imperativement etre une puissance de 2 */

  yi = xiyi >>9;		/* yi va de 0 a VERT par 2 */
  yi += VERT;		/* maintenant yi va de VERT a 2*VERT-2 par 2 */

  xiyi = xi | (yi<<9);	/* on remet ensemble xi et yi modifies */
/* ----- fin de la partie dependante du nb de bits de codage des x ----- */

  phot_ptrj = phot_ptri;	/* pointe sur xj & yj */
  correl_count = integ_count;
while (--correl_count >= 0L) {
  ac_ptr = (short *)(ac + xiyi - *phot_ptrj);
  (*ac_ptr)++;	/* increment target pixel in ac */
  phot_ptrj++;	/* next xj & yj */
  }
 }
}	

/* -------- inter-correlation Ft en flux continu -------- */

void de_correlate9 (phot_ptr, ph_per_correl, ph_per_integ, ph_per_decorr, de_ac_org)
photon	*phot_ptr;	/* pointer to current photon */
long	ph_per_correl;	/* equivalent to a correlation time */
long	ph_per_integ;	/* equivalent to an integration time */
long	ph_per_decorr;	/* equivalent to a de_correlation time */
unsigned short	*de_ac_org; /* pointer to origin af target array (ac storage) */
{
register short	*de_ac_ptr;	/* points current pixel in de_ac */
long		integ_count;	/* photons left to integrate */
register long	correl_count;	/* photons left to correlate with current photon */
register unsigned long	xiyi;	/* current photon coordinate : xi & yi */
unsigned long	xi;		/* separate coordinates */
unsigned long	yi;
register photon	*phot_ptrj;	/* will point ph.to correlate with current photon */

register long	de_ac = (long)de_ac_org;

integ_count = ph_per_integ;
if (integ_count < 1L)
	return;
	
while (--integ_count >= 0L)
{
  xiyi = *phot_ptr++;
/* ----- attention, ici on fait l'hypothese que xi est code sur 9 bits ----- */
  xi = xiyi & 0x1FF;	/* xi va de 0 a 510 par 2  (0  HORIZ-2) */
  xi |= 0x200;		/* maintenant xi va de 512 a 1022 par 2 (HORIZ a 2*HORIZ-2) */
/* HORIZ doit imperativement tre une puissance de 2 */

  yi = xiyi >>9;		/* yi va de 0 a VERT par 2 */
  yi += VERT;		/* maintenant yi va de VERT a 2*VERT-2 par 2 */

  xiyi = xi | (yi<<9);	/* on remet ensemble xi et yi modifies */
/* ----- fin de la partie dependante du nb de bits de codage des x ----- */

  phot_ptrj = phot_ptr + ph_per_decorr;	/* points to xj & yj */
  correl_count = ph_per_correl;
  while (--correl_count >= 0L)
  {
    de_ac_ptr = (short *)(de_ac + xiyi - *phot_ptrj);
    (*de_ac_ptr)++; 
    phot_ptrj++;
  }
 }
}

/* ----- inter-correlation between two photon data sequences (images) ----- */

void inter_correlim9 (phot_ptrA, phots_dans_A, phot_ptrB, phots_dans_B, de_ac_org)

photon	*phot_ptrA;	/* pointe vers le 1er ph. de la sequence A a intercor */
long	phots_dans_A;	/* nb de photons dans la seqence A */
photon	*phot_ptrB;	/* pointe vers le 1er ph. de la sequence B a intercor */
long	phots_dans_B;	/* nb de photons dans la seqence B */
unsigned short	*de_ac_org;	/* pointe vers origine tableau but (intercorrelation) */
{
register short	*de_ac_ptr;	/* pixel but dans intercorr */
long		ph_count_A;	/* photons restant a integ dans A */
register long	ph_count_B;	/* photons restant a integ dans B */
register unsigned long	xiyi;	/* couple de coordonnee : xi & yi */
unsigned long	xi;		/* coordonnee seule */
unsigned long	yi;
photon		*phot_ptri;
register photon	*phot_ptrj;

register long	de_ac = (long)de_ac_org;

phot_ptri = phot_ptrA;
ph_count_A = phots_dans_A;
	
if (ph_count_A < 1L)
	return;
	
while (--ph_count_A >= 0L)
 {
  xiyi = *phot_ptri++;
/* ----- attention, ici on fait l'hypothese que xi est code sur 9 bits ----- */
  xi = xiyi & 0x1FF;	/* xi va de 0 a 510 par 2  (0  HORIZ-2) */
  xi |= 0x200;		/* maintenant xi va de 512 a 1022 par 2 (HORIZ  2*HORIZ-2) */
/* HORIZ doit imperativement tre une puissance de 2 */
  yi = xiyi >>9;	/* yi va de 0 a VERT par 2 */
  yi += VERT;		/* maintenant yi va de VERT a 2*VERT-2 par 2 */
  xiyi = xi | (yi<<9);	/* on remet ensemble xi et yi modifies */
/* ----- fin de la partie dependante du nb de bits de codage des x ----- */
  phot_ptrj = phot_ptrB; /* points to xj & yj */
  ph_count_B = phots_dans_B;
  while (--ph_count_B >= 0L)
   {		
     de_ac_ptr = (short *)(de_ac + xiyi - *phot_ptrj);
     (*de_ac_ptr)++; 
     phot_ptrj++;
   }
 }
}

/* ---------------- intercorrelations Fx en flux continu ---------------- */
void inter_correlflux9 (phot_ptr, ph_per_correl, ph_per_integ, ac_org)
photon	*phot_ptr;	/* pointer to current photon */
long	ph_per_correl;	/* photons per correlation time */
long	ph_per_integ;	/* photons per integration time */
unsigned short	*ac_org;  /* pointer to origin af target array (storing AC) */

/* integre l'ic  2-dim a partir de coordonnees de photons en memoire.
renvoie le pointeur de photons mis a jour.
L'ic est calculee par soustraction de coordonnee d'un photon a n autres photons.
L'appartenance des coordonnees des photons a l'une ou l'autre des fenetres
d'intercorrelation a ete pre-marquee par la routine "marque et centre"
en mettant a 1 l'un ou/et l'autre des deux bits de poids fort (mot de 32 bits).
Ceci permet d'accelerer le calcul, en reduisant le temps necessaire aux tests
dans la boucle interne executee n*m fois : le marquage a necessite
au maximum n + m operations.
*/
{
register photon	*phot_ptrj;	/* pointe vers le ph.a correler avec ph. en cours */
register short	*ac_ptr;	/* pixel but dans l'ac */
long		integ_count;	/* photons restant a integrer */
register long	correl_count;	/* photons restant a correler avec photon en cours */
	
register unsigned long	xiyi;	/* coordonnee couplee d'un photon : xi & yi */
register long	xjyj;		/* pas "unsigned" pour que le test sur msb marche */
unsigned long	xi, yi;
unsigned long	bit_fenetre_i = 0X40000000; /* bit a 1 => appartenance a fen. i */
register unsigned long	ac = (long)ac_org;	/* origine a ajouter */
register unsigned long	mask = 0x3FFFFFFF;	/* pour couper les bits de poids forts (aoh!) */
	
integ_count = ph_per_integ;
if (integ_count < 1L)
  return;

while (--integ_count >= 0L) {
  xiyi = *phot_ptr++;
  if ((xiyi & bit_fenetre_i) != 0) { 	/* on est dans la fenetre i */
/* ----- attention, ici on fait l'hypothese que xi est code sur 9 bits ----- */
  xi = xiyi & 0x1FF;		/* xi va de 0  510 par 2  (0 a HORIZ-2) */
  yi = (xiyi >>9) & 0x1FF;/* yi va de 0 a VERT par 2 */
  xi |= 0x200;		/* maintenant xi va de 512 a 1022 par 2 (HORIZ a 2*HORIZ-2) */
  yi += VERT;		/* maintenant yi va de VERT a 2*VERT-2 par 2 */
  xiyi = xi | (yi<<9);  /* on remet ensemble xi et yi modifies */
/* ----- fin de la partie dependante du nb de bits de codage des x ----- */

phot_ptrj = phot_ptr;	/* pointe sur xj & yj */
correl_count = ph_per_correl;
while (--correl_count >= 0L) {
  xjyj = *phot_ptrj;
  if (xjyj < 0 ) {	/* le msb est marque : on est dans la fenetre j */
   xjyj &= mask;	/* ne garder que les coordonnees */
   ac_ptr = (short *)(ac + xiyi - xjyj);	/* faire l'intercorrelation */
   (*ac_ptr)++;		/* incrementer pixel dans ic */
   }
   phot_ptrj++;
  }
}
}
}

/* ------ dump ------ */
void dump_coor (phot_buf)
phot_buf_rec phot_buf;
{
long i;
long x, y, yx;
photon *ph_ptr, **im_ptr;

/* donner les etats */
printf ("adr buf %5lX  adr tab debuts_image %5lX\n",
         phot_buf.adr, phot_buf.debuts_image);
		
printf ("duree %5ld  ph. %5ld   format %5ld   nb ims %5ld\n",
         phot_buf.duree, phot_buf.photons, phot_buf.format, phot_buf.nb_images);
printf ("plein %5d  processed. %5d   reduced %5d   tries %5d\n",
         phot_buf.plein, phot_buf.processed, phot_buf.reduced, phot_buf.tries);	
	
/* Print some coordinates */
ph_ptr = phot_buf.adr;
for (i=0 ; i<lmin (10L,phot_buf.photons) ; i++) 
  {
   yx = *ph_ptr;
   x = yx & 0x1FF;
   y = (yx & 0x3FE00) >>9;
   printf ("ph %3ld   %9lX   y %3lX   x %3lX  ", i, yx, y, x);
   bindis (yx);
   printf ("\n");
   ph_ptr++;
  }
/* Print the addresses of the image beginning */
im_ptr = phot_buf.debuts_image;
for (i=0 ; i<lmin (10L, phot_buf.nb_images -1); i++) 
  {
   printf ("im %3ld   %9lX  ph/im %4ld\n",
            i+1, im_ptr [i], (im_ptr [i+1] - im_ptr [i]));
   }
}

/*----- Displays a number with binary code -----*/
bindis (a)
long a;
{
int i, j;
 for (j = 0 ; j <=3 ; j++) 
  {
    printf ("  ");	
    for (i = 0 ; i <=7 ; i++) 
     {
      printf ("%1ld", (a & 0x80000000) >>31);
      a <<=1;
     }
  }
}
