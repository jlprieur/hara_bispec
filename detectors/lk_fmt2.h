/************************************************************************* 
* "lk_fmt2.h" 
* Needed for LK formats
*
* From many include files from L. Koechlin (Version 04-07-95)
* JLP
* Version 04-03-01
***************************************************************************/

/* "Rect" structure: defined by Macintosh C compiler, and guessed by JLP: */
typedef struct{
int top;
int bottom;
int left; 
int right;
} Rect; 

/********************* "params.h" **********************************/

/* ------- Taille des fenetres et images (pixels affichage:512*512) ------ */
/* ------ fenetre AC et somme ------ */

#define HORIZ	(512)	/* nb de colones dans l'image (doit tre une puissance de 2) */
#define VERT	(512)	/* 3Mo -> 420 ; 2Mo -> 250 nb de l. dans l'ima.*/

/* ------ taille des buffers de photons en octets ------ */

#define PHOT_BUF_SIZE	100000L 
/* 3Mo->30000 2Mo->20000 taille buf en nb de photons*/
#define MAX_IMAGES	(PHOT_BUF_SIZE / 8)  /* minimum 8 photons par image en moyenne */

/* ---------- Cameras CP40 ---------- */
/* champ individuel d'un canal CP40 : 288 x 384 */
#define LGX_CP40 (512)		/* resolutions brutes (codage) */
#define LGY_CP40 (512)

#define DECALX_CP40 (2)		/* ordonnancement des donnees */
#define DECALY_CP40 (14)

#define	masc_c_pair	(0x800000)  /* position des bits indicateurs de canal */
#define	masc_c_impair	(0x800)
/********************* end of "params.h" **********************************/
/********************* "phottypes.h" **********************************/
/* ----- structures pour fichiers de photons en RAM et sur disque ----- */

/* JLP 2001 PC/Linux: int=4bytes, long=4bytes */
typedef long photon;
typedef short photon_papa;

/*
 cette structure contient les infos sur un fichier de coors. de photons sur disque.
 A l'ecriture sur disque, cette structure est recopiee en entete du fichier. 
 A la lecture, elle est relue.
 Les champs refNum, deja_lus, processed, ne sont pas lies aux photons contenus :
 ils sont modifies par l'ouverture du fichier ou par le traitement.
*/
typedef struct {
long	duree_fich;	/* duree correspondante au fichier en ms */
unsigned long	date;	/* peu utilise : redondant avec la date de creation ou de modif */
long	photons;	/* nombre de photons contenus */
long	clef_format;	/* contient 'FORM' si le format qui suit est valide */
/* JLP95: previously: long format; */
char	format[5];	/* nom du format : 'PAPA' 'MAMA' 'CP40' 'CP20' 'YX9' etc. */
long	nb_images;	/* si les donnees sont decoupees en images, 0 sinon */
short	refNum;		/* file ref number quand elle est ouverte */
long	deja_lus;	/* position du pointeur (photons deja lus) sur le fichier */
int	tout_lu;	/* vrai si on est arrive au bout du fichier */
} phot_file_rec;


/* cette structure contient les infos 
sur un buffer de coordonnees. de photons en RAM. */
			
typedef struct {
long	duree;		/* duree correspondante au contenu du buffer en ms */
long	photons;	/* nombre de photons contenus */
/* JLP95: previously: long format; */
char	format[5];	/* nom du format : 'PAPA' 'MAMA' 'CP40' 'CP20' 'YX9' etc. */
long	nb_images;	/* si les donnees sont decoupees en images, 0 sinon */
long	t_calcul;
photon	*adr;		/* adresse de base du buffer */
photon	**debuts_image;	/* pointe vers le tableau contenant les adresses
					du 1er photon de chaque image */
int	plein;		/* vrai si le buffer contient des coordonnees */
int	processed;	/* vrai si les coordonnees ont ete traitees */
int reduced;		/* vrai si les coors  a resolution reduite (pour AC) */
int zoom_y; 		/* vrai si les coors  zoomees en y */
int zoom_x; 		/* vrai si les coors  zoomees en x */
int zoom_yx; 		/* vrai si les coors  zoomees en yx */
float	echelle_zoomy;
int dedistor; 		/* vrai si les coors  corrigees de la distorsion CP40 */
int tries; 		/* vrai si les coors  triees en zones dans l'omage */
} phot_buf_rec;


/*
 coor de phots : cette structure va remplacer pour la TF les phot_buf_rec
 definis ci-dessus ; elle occupe plus de memoire mais les donnees
 sont plus rapidement accessibles.
*/
typedef struct {	/* un photon a coors separees */	
short	x;
short	y;
short	t;
} photon_xyt;


typedef struct {		/* liste de photons a coors separees */
photon_xyt	*ph_ptr;	/* pointe vers la liste des coors de phots */
long		nb_phots;	/* nombre de photons contenus */
long		lgx;		/* amplitude des valeurs possibles de x */
long		lgy;		/* amplitude des valeurs possibles de y */
long		lgt;		/* amplitude des valeurs possibles de t */
long		resol_t;
long	duree;			/* duree correspondante au contenu du buffer en ms */
} phots_xyt;



/* --------------- structure pour les differents types de camera CP40 -------------- */

typedef struct {
unsigned long debut_im;
unsigned long masc_debut_im;
/* indicateur de canal */
unsigned long masc_c2_4;		
unsigned long masc_c1_3;
int     h_c1;        /* recalage horizontal du canal 1 */ 
int     v_c1;        /* recalage vertical   du canal 1 */ 
int     h_c2;        /* recalage horizontal du canal 2 */ 
int     v_c2;        /* recalage vertical   du canal 2 */ 
int     h_c3;        /* recalage horizontal du canal 3 */ 
int     v_c3;        /* recalage vertical   du canal 3 */ 
int     h_c4;        /* recalage horizontal du canal 4 */ 
int     v_c4;        /* recalage vertical   du canal 4 */ 
long		*x_distor;
long		*y_distor;
} camera;


/* ----------- structures pour preparation des donnees ------------- */

typedef struct {
int	on;		/* vrai=> afficher les cadrages et clipper */
Rect	im;		/* diaphragme sur les images (champ A) */
Rect	ic;		/* champ B (a intercorreler avec le champ A) */
}
cadrage;
/********************* end of "phottypes.h" **********************************/

/* --------------- protos --------------- */
void set_CP4N (camera *);
void set_CP4A (camera *);
void set_CP4T (camera *);
void set_CQ40 (camera *);

void CP40_vers_YX9M (phot_buf_rec *obsA_ptr, camera  *quelle_CP40,
     double *image1, int nx1, int ny1, int ixstart, int iystart, int ireduc,
     double *re1, double *im1, double *re2, double *im2,
     double *modsq, double *snrm, double *yce1, double *long_int,
     float *ffield, int nx_ff, int ny_ff,
     int ir, int nbeta, int ngamma, long *iphot, long *ixy_1,
     long max_nphot, long distor_corr, int mode_test);

void CP40_distor_vers_YX9M (phot_buf_rec *obsA_ptr, camera *quelle_CP40, 
                            int mode_test);

void PAPA_vers_YX9 (phot_buf_rec *, int mode_test);
void YCAR_vers_YX9 (phot_buf_rec *, int mode_test);

void prepa9_ac	(phot_buf_rec *);
void zoom_Y9 	(phot_buf_rec *, Rect);
void zoom_X9	(phot_buf_rec *, Rect);
void zoom_YX9	(phot_buf_rec *, Rect);
void reduc_YX9	(phot_buf_rec *);
void filtre	(phot_buf_rec *source, phot_buf_rec *destination, Rect);

void dump_coor	(phot_buf_rec phot_buf);

void marque_et_centre (phot_buf_rec *, phot_buf_rec *, Rect, Rect);	
/* prepa a l'IC */
