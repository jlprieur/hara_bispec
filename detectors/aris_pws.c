// Traitement de speckles d']toile double en comptage
// Calcul de densite spectrale - seuillage des images a 3 sigmas 
// Centrage des photons : la tache devient un delta pondere par son integr.
// necessite un jeu d'images double et un jeu de reference

#include <stdio.h>
#include <math.h>
#include "/Robert/aris/bib/mylib.c"

#define sqr(a) ((a)*(a))
#define DEBUG

// -------------------------------------------------------------------
// fonction de transfert d)un telescope circulaire
// -------------------------------------------------------------------
float transf(float f)
{
	if (abs(f)>1) return 0;
	else
	return 2/3.14*(acos(f)-f*sqrt(1-f*f));
}


seuillage(ima,dim)

int	dim;
float	*ima;
{
	int seuil=20;
	int i,j,imax,jmax;
	int ii,jj;
	float pmax;
	
	
	for (i=1; i<dim-2; i++) for (j=1; j<dim-2; j++) if (*(ima+i*dim+j) < seuil) *(ima+i*dim+j)=0;
	else
	/* Centrage */
	{
		pmax=0;
		for (ii=0; ii<3; ii++)
		for (jj=0; jj<3; jj++)
		if ( pmax < *(ima+(i+ii)*dim+j+jj) )
		{
			pmax = *(ima+(i+ii)*dim+j+jj);
			imax=i+ii;
			jmax=j+jj;
		};

		for (ii=-1; ii<2; ii++)	// on fait le vide
		for (jj=-1; jj<2; jj++)
		*(ima+(imax+ii)*dim+jmax+jj)=0;

		*(ima+imax*dim+jmax)=pmax;	// on pose un photon
		
	}

	/* gestion des bords */
	for (i=dim-2; i<dim; i++) for (j=0; j<dim; j++) *(ima+i*dim+j)=0;
	for (i=0; i<dim; i++) for (j=dim-2; j<dim; j++) *(ima+i*dim+j)=0;
	for (j=0; j<dim; j++)  *(ima+j)=0;		// i=0
	for (i=0; i<dim; i++)  *(ima+i*dim)=0;		// j=0
	
	return(0);
}


// -------------------------------------------------------------------
// normalisation au max et regularisation
// la regularisation se fait en mettant les termes
// au dela de la frequence de coupure a zero
// -------------------------------------------------------------------

norm(tabr,dim,fc)

float	*tabr;
int	dim;
float 	fc;

{
	float	moy=0;
	float	maxi=-1e30;
	float	r;		
	long	i,j;
	
	for (i=0; i<dim*dim; i++) {moy+=*(tabr+i)/dim/dim; maxi=MAX(maxi,*(tabr+i));}
	
	for (i=0; i<dim; i++) for (j=0; j<dim; j++) 
	{
		*(tabr+i*dim+j)/=maxi; 
		r=sqrt((i-dim/2)*(i-dim/2)+(j-dim/2)*(j-dim/2));
		if (r>=fc) *(tabr+i*dim+j)=0;
	}
}


// -------------------------------------------------------------------
// Calcul du spectre
// -------------------------------------------------------------------

isp(nomgen,nb,format,dim,pws,bias,normalisation)

char	*nomgen;
float	*pws;
int	dim,format,bias,normalisation;

{
	long	i,j;
	char	nomfic[80];
	float	*ima,*tmp;
	int	k,isign=1;
	float	val_norm=0,valeur_bis=0;
	
	ima=(float *)malloc(dim*dim*4);
	tmp=(float *)malloc(dim*dim*4);

	raz(pws,dim*dim);
	
	for (k=1; k<=nb; k++)
	{
		valeur_bis=0;
		sprintf(nomfic,"%s%d",nomgen,k);
		printf("Image no %d : %s\n",k,nomfic);
		switch(format)
		{
			case 0 : {lect(ima,nomfic,dim*dim); break;}
			case 1 : {ascopen(nomfic,dim,dim,ima);break;}
			case 2 : {ropen(nomfic,dim,dim,ima); break;}
		}
		raz(tmp,dim*dim);
		seuillage(ima,dim);
		
		isign=1; tf2d_mod(ima,tmp,dim,isign);
		
		if (bias)
			for (i=0; i<10; i++) for (j=0; j<4; j++) 
		   	valeur_bis+=*(ima+i*dim+j)/40.;
		printf("fond normalise=%e\n",valeur_bis/(*(ima+dim*dim/2+dim/2)));

		for (i=0; i<dim*dim;i++)
		{
				*(ima+i)-=valeur_bis;
				*(ima+i)*=*(ima+i); *(ima+i)/=nb*dim*dim;
				*(pws+i)+=*(ima+i);
		}

	}
	printf("------------------\n");
}

// -------------------------------------------------------------------
// Calcul de l'autocorrelation corrigee de la psf
// On utilise une apodisation par la fn de transert du tel.
// -------------------------------------------------------------------

acorr(pwsdbl,pwsref,acobj,dim,fc)

float	*pwsdbl,*pwsref,*acobj;
int	dim;
float	fc;

{
	float 	*tmp;
	float	den,pmin1=1e30,pmin2=1e30;
	int 	i,compteur=0,j;
	int	isign;

	tmp=(float *)malloc(dim*dim*sizeof(float));
	
//      Soustraction du min de pwsdbl et pwsref
//      ---------------------------------------
	for (i=0; i<dim*dim;i++)
	{
                pmin1=MIN(pmin1,*(pwsdbl+i));
                pmin2=MIN(pmin2,*(pwsref+i));
        }
	for (i=0; i<dim*dim;i++)
	{
                *(pwsdbl+i)-=pmin1;
                *(pwsref+i)-=pmin2;
        }

	for (i=0; i<dim;i++)
	for (j=0; j<dim;j++)
	{
		den=*(pwsref+i*dim+j);
		if (den==0) *(acobj+i*dim+j)= *(tmp+i*dim+j)=0; else
		{
			compteur++;
			*(acobj+i*dim+j)=*(pwsdbl+i*dim+j)/den*transf(hypot((float)i-dim/2,(float)j-dim/2)/fc);
			*(tmp+i*dim+j)=0;
		}
	}
	ascwrite("reel",dim,dim,acobj);
	ascwrite("imag",dim,dim,tmp);
	isign=1; tf2d_mod(acobj,tmp,dim,isign);
}


// ===============================================================================
// Programme principal
// ===============================================================================

main()
{
	char	nomgendbl[60],nomgenref[60];
	int	dim,format,nbdbl,nbref,bias,flag,normalisation;
	float	fc;		// rayon freq. de coupure (valeur demandee plus tard)
	float	*pwsref,*pwsdbl,*acobj;
	
//	Input des parametres
//	--------------------
	printf("[0] Calcul des spectres ou [1] skip : ");
	scanf("%d",&flag);
	
	if (flag==0)
	{
		printf("Nom generique images double =");
		scanf("%s",&nomgendbl);
	
		printf("Nom generique images reference =");
		scanf("%s",&nomgenref);

		printf("Nombre images double =");
		scanf("%d",&nbdbl);

		printf("Nombre images ref =");
		scanf("%d",&nbref);

		printf("Dimension des images =");
		scanf("%d",&dim);
	
		printf("format [0=bin, 1=asc, 2=float] :");
		scanf("%d",&format);

//		Le biais est noemalement soustrait en comptage de photons. 
//		il est calcule sur les points au dela de la fc

		printf("Traitement du biais  [0=non 1=oui] :");
		scanf("%d",&bias);

		normalisation=0;
		
	}
	else
	{
		printf("Densite spectrale double =");
		scanf("%s",&nomgendbl);
		
		printf("Densite spectrale simple =");
		scanf("%s",&nomgenref);
		
		printf("format [1=asc, 2=float 3=fits] :");
		scanf("%d",&format);

		if (format!=3)
		{
		  printf("Dimension des images =");
		  scanf("%d",&dim);
		}
	}	
	
	printf("Freq. de coupure (frequels) =");
	scanf("%f",&fc);
	

	if (format!=3)
	{
	  acobj=(float *)malloc(dim*dim*sizeof(float));
	  pwsref=(float *)malloc(dim*dim*sizeof(float));
	  pwsdbl=(float *)malloc(dim*dim*sizeof(float));
	}

	if (flag==0)
	{
// 		Calcul du spectre a partir des images
//		-------------------------------------

		printf("Calcul des spectres ------------\n");
		isp(nomgendbl,nbdbl,format,dim,pwsdbl,bias,normalisation);
		ascwrite("pwsdbl",dim,dim,pwsdbl);

		isp(nomgenref,nbref,format,dim,pwsref,bias,normalisation);
		ascwrite("pwsref",dim,dim,pwsref);
	}
	else
	{
		switch(format)
		{
			case 1 : {ascopen(nomgendbl,dim,dim,pwsdbl);
				  ascopen(nomgenref,dim,dim,pwsref);
				  break;}
			case 2 : {ropen(nomgendbl,dim,dim,pwsdbl);
				  ropen(nomgenref,dim,dim,pwsref);
				  break;}
			case 3 : {fitsopen(nomgendbl,&dim,&bias,&pwsdbl);
				  fitsopen(nomgenref,&dim,&bias,&pwsref);
				  acobj=(float *)malloc(dim*dim*sizeof(float));
				  break;}
		}
	
	}

	acorr(pwsdbl,pwsref,acobj,dim,fc);
	if (format==3) fitswrite("acorr.fit",dim,dim,acobj);
	ascwrite("acorr",dim,dim,acobj);

}

/* -----------------------------------*/
/* Lecture des fichiers non formattes */
/* ----------------------------------*/

lect(tab,nomfic,dim)

float *tab;
char nomfic[];
int dim;

{
 register int i;
 FILE *fic,*fopen();

 if ((fic=fopen(nomfic,"r"))!=NULL)
 {
  for (i=0;i<dim;i++) *tab++=(float) getc(fic);
  fclose(fic);
 }
 else puts("Pb ouverture fichier ");

/* tab[i + j * dim] = 
Attentiomn, inversion en Y opar rapport au foichiers TIFF
*/
}

