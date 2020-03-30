/********************************************************************
* Set of routines for uv-coverage, spectral and bispectral lists:
*   
* Contains:
* int COVERA_MASK(mask,nx_mask,ny_mask,ir,max_nclosure,nbeta,ngamma)
* int COVERA(INT4 *ir, INT4 *max_nclosure, INT4 *nbeta, INT4 *ngamma)
* static int couv_mask(mask,nxm,nym,ir,nbeta,ngamma,nbc_dim,nbc_offset)
* static int affect_mask(mask,nxm,nym,isx,isy,irs,nb,ng,
*                       nbc_dim,nbc_offset,max_nclosure,dimension_only)
* int cover_mask_to_eric(bisp_list,spec_list,dimension_file,nbeta,ngamma)
* int modsq_to_eric(modsq,nx,ny,mod_ascii_file,nbeta)
* int bisp2D_to_2D_image(out_image,modsq,nx,ny,bisp,nx_bisp,nbeta,ngamma,iopt)
* int bisp1D_to_2D_image(out_image,nx_out,ny_out,modsq,nx_modsq,bisp,nx_bisp,
*                        nbeta,ngamma,line_to_output,iopt)
* int BISPEC3(re,im,modsq,snrm,nx,ny,bispp,ir,nbeta,ngamma)
* int BISPEC3_SINGLE(re,im,modsq,snrm,nx,ny,bispp,ir,nbeta,ngamma)
* int photon_corr(bispp,modsq,snrm,nx,ny,xphot,nbeta,ngamma,photon_correction)
* int photon_corr_mask(bispp,modsq,phot_modsq,nx,ny,bisp_dim,xphot,nbeta,ngamma,fraction)
* int calib_ref(bispp,modsq,bispp_ref,modsq_ref,nx,ny,bisp_dim,nbeta,ngamma)
* int rearrange_mask(bispp,ngamma,ydim1)
* int output_lists_coverage(nbeta,ngamma)
*
* *** 1D routines:
* int COVERA_MASK_1D(mask,nx_mask,ir,max_nclosure,nbeta,ngamma)
* static int couv_mask_1D(mask,nxm,ir,nbeta,ngamma,max_nclosure,dimension_only)
* int photon_corr_1D(bispp,modsq,nx,ny,bisp_dim,nbeta,ngamma)
* int bispec_1D_X(re,im,modsq,snrm,nx,ny,bispp,ir,nbeta,ngamma)
* int bispec_1D_Y(re,im,modsq,snrm,nx,ny,bispp,ir,nbeta,ngamma)
* static int affect_mask_1D(mask,nxm,isx,nb,ng,max_nclosure,dimension_only)
* int rearrange_mask_1D(bispp,ngamma,ydim1,ny)
* int output_lists_coverage_1D(nbeta,ngamma)
*
*
* Fortran interface to access to hidden variables:
* int COVER_NGT(ngt_val,index)
* int COVER_IXY(ixy1_val,ixy2_val,index)
* int COVER_IXY_1D(ixy1_val,index)
* int COVER_NBCOUV(nb_val,index1,index2,ir_max)
* int COVER_KLM(klm_val,klm_index,ng_index)
*
* C interface to access to hidden variables:
* int cover_ngt1(ngt_val,index)
* int cover_klm0(klm_val,klm_index,ng_index)
* int cover_klm1(klm_val,klm_index,ng_index)
*
* Changes:
* Dec 94: double precision arrays for modsq, re, im, and bisp1, ...)
*
* JLP
* Version 02-09-95
*******************************************************************/
#include <stdio.h>
#include <stdlib.h> /* malloc */
#include <math.h>
#include <jlp_ftoc.h>

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

#define IRMAX 50 
#define NBMAX 3922 

#include "jlp_cover_mask.h"
/*
* Defined in "jlp_cover_mask.h" :

int COVERA_MASK(float *mask, INT4 *nx_mask, INT4 *ny_mask, INT4 *ir,
                INT4 *max_nclosure, INT4 *nbeta, INT4 *ngamma);
int COVERA(INT4 *ir, INT4 *max_nclosure, INT4 *nbeta, INT4 *ngamma);
int COVER_NGT(INT4 *ngt_val, INT4 *index);
int COVER_IXY(INT4 *ixy1_val, INT4 *ixy2_val, INT4 *index);
int COVER_IXY_1D(INT4 *ixy1_val, INT4 *index);
int COVER_NBCOUV(INT4 *nb_val, INT4 *index1, INT4 *index2, INT4 *ir_max);
int COVER_KLM(INT4 *klm_val, INT4 *klm_index, INT4 *ng_index);
int cover_ngt1(INT4 *ng2, INT4 nb);
int cover_klm1(INT4 *k, INT4 kk, INT4 ng);
int cover_klm0(INT4 *k, INT4 kk, INT4 ng);
int corr_bisp_sky(double *bispp, double *modsq, INT4 *nx, INT4 *ny,
                  INT4 *bisp_dim, float *sky_level, INT4 *nbeta, 
                  INT4 *ngamma);
int photon_corr(double *bispp, double *modsq, double *snrm,
                INT4 *nx, INT4 *ny, float *xphot, INT4 *nbeta, INT4 *ngamma,
                INT4 *photon_correction);
int BISPEC3(double *re, double *im, double *modsq, double *snrm,
            INT4 *nx, INT4 *ny, double *bispp, INT4 *ir, 
            INT4 *nbeta, INT4 *ngamma);
int bispec_1D_Y(double *re, double *im, double *modsq, double *snrm,
                INT4 *nx, INT4 *ny, double *bispp, INT4 *ir, INT4 *nbeta, 
                INT4 *ngamma);
int cover_mask_to_eric(char *bisp_list, char *spec_list, char *dimension_file,
                       INT4 *nbeta, INT4 *ngamma);
int modsq_to_eric(double *modsq, INT4 *nx, INT4 *ny, char *mod_ascii_file,
                  INT4 *nbeta);
int bisp_to_eric(double *bisp, INT4 *idim, char *bisp_ascii_file, INT4 *ngamma);
int photon_corr_1D(float *bispp, float *modsq, INT4 *nx, INT4 *ny, 
                   INT4 *bisp_dim, INT4 *nbeta, INT4 *ngamma);
*/

static int couv_mask(float *mask, INT4 nxm, INT4 nym, INT4 *ir, INT4 *nbeta,
                     INT4 *ngamma, INT4 nbc_dim,
                     INT4 nbc_offset, INT4 max_nclosure, INT4 dimension_only);
static int couv_no_mask(INT4 *ir, INT4 *nbeta, INT4 *ngamma, INT4 nbc_dim,
                      INT4 nbc_offset, INT4 max_nclosure, INT4 dimension_only);
static int affect_mask(float *mask, INT4 nxm, INT4 nym, INT4 isx, INT4 isy, 
                       INT4 irs, INT4 *nb, INT4 *ng, INT4 nbc_dim, 
                       INT4 nbc_offset, INT4 max_nclosure, INT4 dimension_only);
static int affect_no_mask(INT4 isx, INT4 isy, INT4 irs, INT4 *nb, INT4 *ng, 
                          INT4 nbc_dim, INT4 nbc_offset, INT4 max_nclosure, 
                          INT4 dimension_only);
static int couv_mask_1D(float *mask, INT4 nxm, INT4 *ir, INT4 *nbeta, 
                        INT4 *ngamma, INT4 max_nclosure, INT4 dimension_only);
static int affect_mask_1D(float *mask, INT4 nxm, INT4 isx, INT4 *nb, INT4 *ng,
                          INT4 max_nclosure, INT4 dimension_only);


/* Keep the addresses of the arrays as hidden variables */
static INT4 *nbcouv, *ixy_1, *ixy_2, *ngt, *klm;

/********************************************************************
* Subroutine COVERA_MASK
* To compute the frequencies of the uv coverage, 
* and the elements of the A matrix.
* (To solve eventually the equation A*X=Y, 
*   i.e., to invert the bispectral relations)
*
* Compute u-v coverage with the mask and within a disk of radius ir.
*
* INPUT:
*  ir: maximum radius of the uv-coverage
*  max_nclosure: maximum number of closure relations
*  mask[nx * ny]: frequency mask ( >=1. if frequency accessible, 0. otherwise) 
*
* OUTPUT:
*  NBETA: Number of elements of the spectral list (Number of columns of A)
*  NGAMMA: Number of elements of the bispectral list (Number of rows of A)
*
* OUTPUT values that are hidden as static variables:
* INT4 *nbcouv, *ixy_1, *ixy_2, *ngt, *klm;
* 
*********************************************************************/
int COVERA_MASK(float *mask, INT4 *nx_mask, INT4 *ny_mask, INT4 *ir,
                INT4 *max_nclosure, INT4 *nbeta, INT4 *ngamma)
{
/* max_nclosure: max number of closure relations per spectral term */
/* Limitation of the maximum number of closure relations
* (to simulate Knox-Thompson and/or accelerate computations...) */ 
INT4 isize, nxm, nym, ir_max, nbc_dim, nbc_offset, iw1;
register INT4 i, k;

/* JLP 94
 ir_max = IRMAX;
*/
ir_max = *ir;
if(ir_max > IRMAX)
   {printf(" COVERA_MASK/Fatal error: ir max = %d \n",IRMAX);
    exit(-1);
   }
nbc_dim = 2 * ir_max + 1; 
nbc_offset = ir_max + 1; 

nxm = *nx_mask;
nym = *ny_mask;
/* Check first that mask is correct: */
if(nxm < 2 * (*ir) || nym < 2 * (*ir))
  {
   printf(" COVERA_MASK/Fatal error: mask is too small !!!\n");
   printf(" nxm = %d nym = %d ir = %d \n",nxm,nym,*ir);
   exit(-1);
  }
 
/* nbcouv( x from -IRMAX to IRMAX,  y from 0 to IRMAX ) */
/* JLP94 */
isize = (ir_max + 1) * nbc_dim * sizeof(INT4);
JLP_GVMI(&nbcouv,&isize);
/* ixy( 1 for x; 2 for y,  nb from 0 to NBMAX ) */
isize = (NBMAX + 1) * sizeof(INT4);
JLP_GVMI(&ixy_1,&isize);
JLP_GVMI(&ixy_2,&isize);
JLP_GVMI(&ngt,&isize);
/*JLP 94 */
for(i = 0; i < nbc_dim * (ir_max + 1); i++) nbcouv[i] = 0;
 
/* Computing the A matrix as generated by the uv-coverage 
   (defined only by IR for a full pupil) */

/* First call to compute bispectrum size NGMAX = ngamma: */
couv_mask(mask,nxm,nym,ir,nbeta,ngamma,nbc_dim,nbc_offset,*max_nclosure,1);

/* Allocation of memory for bispectrum: */
isize = 3 * (*ngamma) * sizeof(INT4);
JLP_GVMI(&klm,&isize);

/* Second call to compute spectral and bispectral lists: */
couv_mask(mask,nxm,nym,ir,nbeta,ngamma,nbc_dim,nbc_offset,*max_nclosure,0);
 
printf(" covera/uv-coverage, IR = %d \n",*ir);
printf("  NBETA (spectral list) = %d \n",*nbeta);
printf("  NGAMMA (bispec. list) = %d \n",*ngamma);
 
/* Computing the number of closure relations for some values of NB: */
for ( k = 2; k < 5 ; k++)
  {
   i = 3 + (*nbeta * k) /5;
   iw1 = sqrt((double)(ixy_1[i] * ixy_1[i] + ixy_2[i] * ixy_2[i]));
   printf(" NB = %d irad=%d Closure relations: %d \n",
             i,iw1,ngt[i] - ngt[i - 1]);
  }
 
return(0);
}
/********************************************************************
* Subroutine COVERA (same as COVERA_MASK, but without a mask)
* To compute the frequencies of the uv coverage, 
* and the elements of the A matrix.
* (To solve eventually the equation A*X=Y, 
*   i.e., to invert the bispectral relations)
*
* Compute u-v coverage within a disk of radius ir.
*
* INPUT:
*  ir: maximum radius of the uv-coverage
*  max_nclosure: maximum number of closure relations
*
* OUTPUT:
*  NBETA: Number of elements of the spectral list (Number of columns of A)
*  NGAMMA: Number of elements of the bispectral list (Number of rows of A)
*
* OUTPUT values that are hidden as static variables:
* INT4 *nbcouv, *ixy_1, *ixy_2, *ngt, *klm;
* 
*********************************************************************/
int COVERA(INT4 *ir, INT4 *max_nclosure, INT4 *nbeta, INT4 *ngamma)
{
/* max_nclosure: max number of closure relations per spectral term */
/* Limitation of the maximum number of closure relations
* (to simulate Knox-Thompson and/or accelerate computations...) */ 
INT4 isize, ir_max, nbc_dim, nbc_offset, iw1;
register INT4 i, k;

/* JLP 94
 ir_max = IRMAX;
*/
ir_max = *ir;
if(ir_max > IRMAX)
   {printf(" COVERA/Fatal error: ir max = %d \n",IRMAX);
    exit(-1);
   }
nbc_dim = 2 * ir_max + 1; 
nbc_offset = ir_max + 1; 

/* nbcouv( x from -IRMAX to IRMAX,  y from 0 to IRMAX ) */
/* JLP94 */
isize = (ir_max + 1) * nbc_dim * sizeof(INT4);
JLP_GVMI(&nbcouv,&isize);
/* ixy( 1 for x; 2 for y,  nb from 0 to NBMAX ) */
isize = (NBMAX + 1) * sizeof(INT4);
JLP_GVMI(&ixy_1,&isize);
JLP_GVMI(&ixy_2,&isize);
JLP_GVMI(&ngt,&isize);
/*JLP 94 */
for(i = 0; i < nbc_dim * (ir_max + 1); i++) nbcouv[i] = 0;
 
/* Computing the A matrix as generated by the uv-coverage 
   (defined only by IR for a full pupil) */

/* First call to compute bispectrum size NGMAX = ngamma: */
couv_no_mask(ir,nbeta,ngamma,nbc_dim,nbc_offset,*max_nclosure,1);

/* Allocation of memory for bispectrum: */
isize = 3 * (*ngamma) * sizeof(INT4);
JLP_GVMI(&klm,&isize);

/* Second call to compute spectral and bispectral lists: */
couv_no_mask(ir,nbeta,ngamma,nbc_dim,nbc_offset,*max_nclosure,0);
 
printf(" -- COVERA/uv-coverage, IR = %d \n",*ir);
printf(" -- NBETA (spectral list) = %d \n",*nbeta);
printf(" -- NGAMMA (bispec. list) = %d \n",*ngamma);
 
/* Computing the number of closure relations for some values of NB: */
for ( k = 2; k < 5 ; k++)
  {
   i = 3 + (*nbeta * k) /5;
   iw1 = sqrt((double)(ixy_1[i] * ixy_1[i] + ixy_2[i] * ixy_2[i]));
   printf(" -- NB = %d irad=%d Closure relations: %d \n",
             i,iw1,ngt[i] - ngt[i - 1]);
  }
 
return(0);
}
/*******************************************************************
* couv_mask defines the uv-coverage and the A matrix:
* Input:
* mask[nx * ny]: frequency mask ( >=1. if frequency accessible, 0. otherwise) 
* IR: maximum radius of the uv-coverage
* dimension_only: flag set to one if only bispectral list dimension is required
*
* Output:
* NBETA, NGAMMA: number of elements of the spectral and bispectral lists
*
* In common blocks (output): the uv-coverage is accessible from 2 sides:
*
* NBCOUV(I,J): uv-coverage (i.e. number of the spectral list for
*              the pixel(I,J) of the spectrum (I=0,J=0 for null frequency)
*
* IXY(1,I) and IXY(2,I) coordinates of the pixel number I in the spectral list
*              (this allows another entry for the uv-coverage)
*
*************************************************************************/
static int couv_mask(float *mask, INT4 nxm, INT4 nym, INT4 *ir, INT4 *nbeta,
                     INT4 *ngamma, INT4 nbc_dim,
                     INT4 nbc_offset, INT4 max_nclosure, INT4 dimension_only)
{
/* 
* nb: beta index, i.e. spectral list index
* ng: gamma index, i.e. bispectral list index
*/
INT4 ir2max, nb, ng;
INT4 i, i2, j2, ii2, irs, icent;
register INT4 j, ir2;
 
/* Coordinate of center of the mask: */
icent = nxm/2 + nxm * nym/2; 
 
/* Easy cases: */
/* First spectral value at (0,0) */
 nb = 0;
 i = 0; j = 0;
 if(mask[ i + icent + nxm * j] > 0.)
   {
   ixy_1[nb] = i; ixy_2[nb] = j;
   nbcouv[i + nbc_offset + nbc_dim * j] = nb;
   }
 
/* Second spectral value at (1,0), nb=1 */
 i = 1; j = 0;
 if(mask[ i + icent + nxm * j] > 0.)
   {
   nb++;
   ixy_1[nb] = i; ixy_2[nb] = j;
   nbcouv[i + nbc_offset + nbc_dim * j] = nb;
   }

/* Third spectral value at (0,1), nb=2 */
 i = 0; j = 1;
 if(mask[ i + icent + nxm * j] > 0.)
   {
   nb++;
   ixy_1[nb] = i; ixy_2[nb] = j;
   nbcouv[i + nbc_offset + nbc_dim * j] = nb;
   }
 
/* Reseting the total number of elements the  bispectral list */
 ng = 0;
 
/* Main loop: work with successive iterations on circles with
* increasing radii.
* Squared radius: IR2 = 2, 3, ... ,IR2MAX
*/
ir2max = *ir * *ir;
for (ir2 = 2; ir2 <= ir2max; ir2++) 
  {
/* Searching for the couples (I,J) such as: I**2 + J**2 = IR2 with I>=J */
    for ( j = 0; j <= *ir; j++) 
      {
       j2 = j*j; 
       i2 = ir2 - j2;
       if ( i2 < j2) break;
         i = (int)sqrt((double)i2);
         ii2 = i * i;
/* Selecting the points defined by each couple (i,j) 
   such that (i*i + j*j= ir2): */
         if(ii2 == i2)
            {
              irs = (int)sqrt((double)ir2);
               affect_mask(mask,nxm,nym,i,j,irs,&nb,&ng,
                           nbc_dim,nbc_offset,max_nclosure,dimension_only);
/* Now use the symmetry relations: */
              if( i != j) 
                affect_mask(mask,nxm,nym,j,i,irs,&nb,&ng,
                           nbc_dim,nbc_offset,max_nclosure,dimension_only);
              if( j != 0) 
                affect_mask(mask,nxm,nym,-j,i,irs,&nb,&ng,
                           nbc_dim,nbc_offset,max_nclosure,dimension_only);
              if( i != j && j != 0) 
                affect_mask(mask,nxm,nym,-i,j,irs,&nb,&ng,
                           nbc_dim,nbc_offset,max_nclosure,dimension_only);
            }
      }
  }
 
/* NBETA: Total number of the spectral list (Number of columns of the X matrix)
* NGAMMA: Total number of the bispectral list (Number of rows of the X matrix)
* (Remember, we have to solve    A*X = Y) */
  *nbeta = nb;
  *ngamma = ng;

  return(0);
}
/*******************************************************************
* couv_no_mask defines the uv-coverage and the A matrix:
* Same as couv_mask, but without a mask
*
* Input:
* IR: maximum radius of the uv-coverage
* dimension_only: flag set to one if only bispectral list dimension is required
*
* Output:
* NBETA, NGAMMA: number of elements of the spectral and bispectral lists
*
* In common blocks (output): the uv-coverage is accessible from 2 sides:
*
* NBCOUV(I,J): uv-coverage (i.e. number of the spectral list for
*              the pixel(I,J) of the spectrum (I=0,J=0 for null frequency)
*
* IXY(1,I) and IXY(2,I) coordinates of the pixel number I in the spectral list
*              (this allows another entry for the uv-coverage)
*
*************************************************************************/
static int couv_no_mask(INT4 *ir, INT4 *nbeta, INT4 *ngamma, INT4 nbc_dim,
                     INT4 nbc_offset, INT4 max_nclosure, INT4 dimension_only)
{
/* 
* nb: beta index, i.e. spectral list index
* ng: gamma index, i.e. bispectral list index
*/
INT4 ir2max, nb, ng;
INT4 i, i2, j2, ii2, irs;
register INT4 j, ir2;
 
/* Easy cases: */
/* First spectral value at (0,0) */
 nb = 0;
 i = 0; j = 0;
 ixy_1[nb] = i; ixy_2[nb] = j;
 nbcouv[i + nbc_offset + nbc_dim * j] = nb;
 
/* Second spectral value at (1,0), nb=1 */
 i = 1; j = 0;
 nb++;
 ixy_1[nb] = i; ixy_2[nb] = j;
 nbcouv[i + nbc_offset + nbc_dim * j] = nb;

/* Third spectral value at (0,1), nb=2 */
 i = 0; j = 1;
 nb++;
 ixy_1[nb] = i; ixy_2[nb] = j;
 nbcouv[i + nbc_offset + nbc_dim * j] = nb;
 
/* Reseting the total number of elements the  bispectral list */
 ng = 0;
 
/* Main loop: work with successive iterations on circles with
* increasing radii.
* Squared radius: IR2 = 2, 3, ... ,IR2MAX
*/
ir2max = *ir * *ir;
for ( ir2 = 2; ir2 <= ir2max; ir2++) 
  {
/* Searching for the couples (I,J) such as: I**2 + J**2 = IR2 with I>=J */
    for ( j = 0; j <= *ir; j++) 
      {
       j2 = j*j; 
       i2 = ir2 - j2;
       if ( i2 < j2) break;
         i = (int)sqrt((double)i2);
         ii2 = i * i;
/* Selecting the points defined by each couple (i,j) 
   such that (i*i + j*j= ir2): */
         if(ii2 == i2)
            {
              irs = (int)sqrt((double)ir2);
               affect_no_mask(i,j,irs,&nb,&ng,nbc_dim,nbc_offset,
                              max_nclosure,dimension_only);
/* Now use the symmetry relations: */
              if( i != j) 
                affect_no_mask(j,i,irs,&nb,&ng,nbc_dim,nbc_offset,
                               max_nclosure,dimension_only);
              if( j != 0) 
                affect_no_mask(-j,i,irs,&nb,&ng,nbc_dim,nbc_offset,
                               max_nclosure,dimension_only);
              if( i != j && j != 0) 
                affect_no_mask(-i,j,irs,&nb,&ng,nbc_dim,nbc_offset,
                               max_nclosure,dimension_only);
            }
      }
  }
 
/* NBETA: Total number of the spectral list (Number of columns of the X matrix)
* NGAMMA: Total number of the bispectral list (Number of rows of the X matrix)
* (Remember, we have to solve    A*X = Y) */
  *nbeta = nb;
  *ngamma = ng;

  return(0);
}
/*******************************************************************
* affect_mask 
*
* Compare with the mask to check if (isx,isy) is a valid uv vector 
* and can be added to the spectral list.
* Then increment nb by one, and look for all new closed triangles whose sides
* are constituted with previous vectors of the spectral list and the new 
* uv vector.
*
* INPUT:
*  mask: uv mask used to invalidate some uv vectors.
*  nxm, nym: X and Y size of mask array
*  isx, isy : coordinates of the uv vector to be validated
*  irs: current limiting radius for searching uv components
*  nb: previous value of the number of items in the spectral list
*  ng: previous value of the number of items in the bispectral list
*  nbc_dim: dimension of nbcouv array
*  nbc_offset: offset used for nbcouv (since there are negative indices)
*  max_nclosure: maximum number of closure relations allowed by the user
*  dimension_only: flag set to one if this routine is only called
*                  to estimate the value of NGAMMA and NBETA
*
* OUTPUT:
*  nb: new value of the number of items in the spectral list
*      (i.e. index of the new uv vector in the spectral list)
*  ng: new value of the number of items in the bispectral list
*      (incremented in this routine by the number of new triangles found
*      that involves the new uv vector) 
********************************************************************/
static int affect_mask(float *mask, INT4 nxm, INT4 nym, INT4 isx, INT4 isy, 
                       INT4 irs, INT4 *nb, INT4 *ng, INT4 nbc_dim, 
                       INT4 nbc_offset, INT4 max_nclosure, INT4 dimension_only)
{ 
INT4 nbs, nbk, nbr, itx, ity, icent;
register INT4 nbq, iklm;

/* Coordinate of center of the mask: */
icent = nxm/2 + nxm * nym/2; 
 
/* First condition: input point has to be accessible */
  if(mask[isx + icent + isy * nxm] <= 0.) return(-1);

/* If so, record the new point of the uv coverage (spectral list): */
 (*nb)++;
 ixy_1[*nb] = isx; ixy_2[*nb] = isy;
 nbcouv[isx + nbc_offset + nbc_dim * isy] = *nb;

/* Searching for the couples associated with the point NBS=NB
* Generating the bispectral list and building the rows of the A matrix:
*/
 nbs = *nb;
 
/* Loop on all the possible points (Q): */
 iklm = 3 * (*ng);
 for( nbq = 1; nbq < nbs; nbq++)
  { 
/* JLP 94: add an exit test when maximum number of closure relations
has been found (to simulate Knox-Thompson and/or accelerate computations...) */ 
   if((*ng - ngt[(*nb) -1]) == max_nclosure) break;

/* Coordinates of the vector T = S - Q */
   itx = isx - ixy_1[nbq];
   ity = isy - ixy_2[nbq];
 
/* Work within the circle of radius IRS, so we can a priori reject
* the points outside the window [-IRS,+IRS]:
*/
   if(itx >= -irs && itx <= irs && ity >= -irs && ity <= irs)
      {
       if(ity > 0 || ( ity == 0 && itx >= 0))
          {
/* Case number 1 (which could be : k=t, l=q, k+l=s) */

/* It should be accessible too: */
            if(mask[itx + icent + ity * nxm] > 0.)
            {
              nbk = nbcouv[ itx + nbc_offset + ity * nbc_dim];
 
/* We select this couple (U,V) if the vector NBK is in [0,NBQ] */
              if( nbk != 0 && nbk <= nbq)
                 {
/* Add new k,l,m coordinates if more than dimension is wanted */
                 if(!dimension_only)
                   {
                   klm[iklm] = nbk; 
                   iklm++;
                   klm[iklm] = nbq; 
                   iklm++;
                   klm[iklm] = nbs;
                   iklm++;
                   }
                 (*ng)++;
                 }
            }
          }
      else
          { 
/* Case number 2 (which could be : r=-t, s=s, r+s=m=q) */
/* r should be accessible too: */
            if(mask[-itx + icent - ity * nxm] > 0.)
            {
             nbr = nbcouv[ -itx + nbc_offset - ity * nbc_dim];
/* We select this couple (R,S) if the vector NBR is in [0,NBS] */
              if( nbr != 0 && nbr <= nbs)
                  {
/* Add new k,l,m coordinates if more than dimension is wanted */
                  if(!dimension_only)
                     {
                     klm[iklm] = nbr; 
                     iklm++;
                     klm[iklm] = nbs; 
                     iklm++;
                     klm[iklm] = nbq;
                     iklm++;
                     }
                  (*ng)++;
                  }
            }
/* Nota: we can only have L=NBS or M=NBS (never K=NBS) */
          }
      } 

/* End of loop on nbq */
}
 
/* NGT(NB) is the number of the last U,V couple of the group NB=NBS=S: */
   ngt[nbs] = *ng;

 return(0);
}
/*******************************************************************
* affect_no_mask 
* (same as affect_mask, but without taking a mask into account)
*
* Add the (isx,isy) uv vector to the spectral list.
* Then increment nb by one, and look for all new closed triangles whose sides
* are constituted with previous vectors of the spectral list and the new 
* uv vector.
*
* INPUT:
*  isx, isy : coordinates of the uv vector to be validated
*  irs: current limiting radius for searching uv components
*  nb: previous value of the number of items in the spectral list
*  ng: previous value of the number of items in the bispectral list
*  nbc_dim: dimension of nbcouv array
*  nbc_offset: offset used for nbcouv (since there are negative indices)
*  max_nclosure: maximum number of closure relations allowed by the user
*  dimension_only: flag set to one if this routine is only called
*                  to estimate the value of NGAMMA and NBETA
*
* OUTPUT:
*  nb: new value of the number of items in the spectral list
*      (i.e. index of the new uv vector in the spectral list)
*  ng: new value of the number of items in the bispectral list
*      (incremented in this routine by the number of new triangles found
*      that involves the new uv vector) 
********************************************************************/
static int affect_no_mask(INT4 isx, INT4 isy, INT4 irs, INT4 *nb, INT4 *ng, 
                          INT4 nbc_dim, INT4 nbc_offset, INT4 max_nclosure, 
                          INT4 dimension_only)
{ 
INT4 nbs, nbk, nbr, itx, ity;
register INT4 nbq, iklm;

/* Record the new point S of the uv coverage (spectral list): */
 (*nb)++;
 ixy_1[*nb] = isx; ixy_2[*nb] = isy;
 nbcouv[isx + nbc_offset + nbc_dim * isy] = *nb;

/* Searching for the couples associated with the point NBS=NB
* Generating the bispectral list and building the rows of the A matrix:
*/
 nbs = *nb;
 
/* Loop on all the possible points (Q): */
 iklm = 3 * (*ng);
 for( nbq = 1; nbq < nbs; nbq++)
  { 
/* JLP 94: add an exit test when maximum number of closure relations
has been found (to simulate Knox-Thompson and/or accelerate computations...) */ 
   if((*ng - ngt[(*nb) -1]) == max_nclosure) break;

/* Coordinates of the vector T = S - Q */
   itx = isx - ixy_1[nbq];
   ity = isy - ixy_2[nbq];
 
/* Work within the circle of radius IRS, so we can a priori reject
* the points outside the window [-IRS,+IRS]:
*/
   if(itx >= -irs && itx <= irs && ity >= -irs && ity <= irs)
      {
       if(ity > 0 || ( ity == 0 && itx >= 0))
          {
/* Case number 1 (which could be : k=t, l=q, k+l=s) */
/* (T=k) + (Q=l) = S */
            nbk = nbcouv[ itx + nbc_offset + ity * nbc_dim];
 
/* We select this couple (U,V) if the vector NBK is in [0,NBQ] */
              if( nbk != 0 && nbk <= nbq) {
/* Add new k,l,m coordinates if more than dimension is wanted */
                 if(!dimension_only) {
                   klm[iklm] = nbk; 
                   iklm++;
                   klm[iklm] = nbq; 
                   iklm++;
                   klm[iklm] = nbs;
                   iklm++;
                   }
                 (*ng)++;
                 }
          }
      else
          { 
/* Case number 2 (which could be : r=-t, s=s, r+s=m=q) */
/* (-T=r) + S = (Q=m) */
             nbr = nbcouv[ -itx + nbc_offset - ity * nbc_dim];
/* We select this couple (R,S) if the vector NBR is in [0,NBS] */
              if( nbr != 0 && nbr <= nbs)
                  {
/* Add new k,l,m coordinates if more than dimension is wanted */
                  if(!dimension_only)
                     {
                     klm[iklm] = nbr; 
                     iklm++;
                     klm[iklm] = nbs; 
                     iklm++;
                     klm[iklm] = nbq;
                     iklm++;
                     }
                  (*ng)++;
                  }
/* Nota: we can only have L=NBS or M=NBS (never K=NBS) */
          }
      } 

/* End of loop on nbq */
}
 
/* NGT(NB) is the number of the last U,V couple of the group NB=NBS=S: */
   ngt[nbs] = *ng;

 return(0);
}
/*******************************************************************
* Output spectral, bispectral list in Eric Anterrieu's format
*
* Input:
* bisp_list, spec_list, dimension_file: file names for bispectral list
*                                       spectral list, and dimension file
*
* Output:
*******************************************************************/
int cover_mask_to_eric(char *bisp_list, char *spec_list, char *dimension_file,
                       INT4 *nbeta, INT4 *ngamma)
{
FILE   *fp;
register INT4 ng, nb, nr;
int k1, k2, k3;
INT4 n_networks = 1;

/***** Output of bispectral list in a file ********************/
if ((fp = fopen(bisp_list,"w")) == NULL)
   {
   printf(" cover_mask_to_eric/Fatal error opening output file: %s \n",
            bisp_list);
   exit(-1);
   }

for (ng = 0; ng < *ngamma; ng++)
  {
  fprintf(fp,"%d\n",ng);

/* Bispectral redundancy (red_bispect) : */
/*
  for (j = 0; j < n_networks; j++)  fprintf(fp,"%d ",red_bispect[ng][j]);
*/
  for (nr = 0; nr < n_networks; nr++)  fprintf(fp,"%d ",1);

/* Spectral indices k, l, m, of the three components of bispectral term: 
  (from 0 to nb, to be compatible with Eric ?) 
 */
  cover_klm0(&k1,1,ng);
  cover_klm0(&k2,2,ng);
  cover_klm0(&k3,3,ng);
  fprintf(fp,"\n%d %d %d\n", k1, k2, k3);
  }
fclose(fp);

/***** Output of spectral list in a file ********************/
if ((fp = fopen(spec_list,"w")) == NULL)
   {
   printf(" cover_mask_to_eric/Fatal error opening output file: %s \n",
            spec_list);
   exit(-1);
   }

for (nb = 1; nb <= *nbeta; nb++)
  {
  fprintf(fp,"%d ",nb-1);
  fprintf(fp,"(%d,%d) ",ixy_1[nb],ixy_2[nb]);

/* Spectral redundancy: */
/*
  fprintf(fp,"%d ",red_spect[nb]);
*/
  fprintf(fp,"%d ",1);

  fprintf(fp,"%d\n",ngt[nb]);

  }

/* As the following is not understood (and not used either), I remove it: */
#ifdef UNDERSTAND
for (nr = 0; nr < n_networks; nr++)
  {
  fprintf(fp,"%d\n",nr);
  for (nb = 1; nb <= *nbeta; nb++)
    {
    fprintf(fp,"  %-d:\t",nb-1);
/* Redundancy, caracteristic function (for FTO) ... */
/*
    fprintf(fp,"(%d-%d)\t",u_v_res[nr][nb][2],u_v_res[nr][nb][1]);
    fprintf(fp,"%d\n",u_v_res[nr][nb][0]);
*/
    fprintf(fp,"(%d-%d)\t",1,1);
    fprintf(fp,"%d\n",1);
    }
  }
#endif

fclose(fp);

/***** Output of dimensions ****************************/
if ((fp = fopen(dimension_file,"w")) == NULL)
   {
   printf(" cover_mask_to_eric/Fatal error opening output file: %s\n",
            dimension_file);
   exit(-1);
   }

/* Number of networks */
fprintf(fp,"%d\n",n_networks);

/* Number of telescopes per network */
for (nr = 0; nr < n_networks; nr++) fprintf(fp,"%d ",1);

/* Number of spectral frequencies */
fprintf(fp,"\n%d\n",*nbeta);

/* Number of bispectral bi-frequencies */
fprintf(fp,"%d\n",*ngamma);

fclose(fp);

return(0);
}
/*******************************************************************
* Output modulus to Eric Anterrieu's format ASCII file: 
*
* Input:
* mod_ascii_file: file name for modulus 
*
* Output:
*******************************************************************/
int modsq_to_eric(double *modsq, INT4 *nx, INT4 *ny, char *mod_ascii_file,
                  INT4 *nbeta)
{
FILE   *fp;
double work;
INT4 icent, inull;
register INT4 nb;

/* Coordinate of the center: */
icent = (*nx)/2 + (*nx) * (*ny)/2; 
 
/***** Output of modulus to a file ********************/
if ((fp = fopen(mod_ascii_file,"w")) == NULL)
   {
   printf(" modsq_to_eric/Fatal error opening output file: %s \n",
            mod_ascii_file);
   exit(-1);
   }

inull = 0;
for (nb = 0; nb < *nbeta; nb++)
 {
  work = modsq[ ixy_1[nb] + icent + (*nx) * ixy_2[nb] ];
/* Compute square root: */
  if(work < 0)
    {
     work = 0.;
     inull++;
    }
  else
     work = sqrt(work);

  fprintf(fp,"%f\n",work);
 }

fclose(fp);

/* Display warning message to alert the user... */
if(inull)printf("mod_to_eric/Warning: %d negative values in modsq !!!\n",inull);

return(0);
}
/*******************************************************************
* Output bispectrum phasor to Eric Anterrieu's format ASCII file: 
*
* Input:
* bisp_ascii_file: file name for output bispectrum 
*
* Output:
*******************************************************************/
int bisp_to_eric(double *bisp, INT4 *idim, char *bisp_ascii_file, INT4 *ngamma)
{
FILE   *fp;
double r0, i0, work;
INT4 inull;
register INT4 ng;

/***** Output of bispectrum to a file ********************/
if ((fp = fopen(bisp_ascii_file,"w")) == NULL)
   {
   printf(" bisp_to_eric/Fatal error opening output file: %s \n",
            bisp_ascii_file);
   exit(-1);
   }

inull = 0;
for (ng = 0; ng < *ngamma; ng++)
 {
  r0 = bisp[ng]; i0 = bisp[ng + (*idim)];
  
/* Compute modulus: */
  work = r0*r0 + i0*i0;
  if(work == 0.)
    {
     inull++;
    }
  else
    {
     work = sqrt(work);
     r0 /= work; i0 /= work;
    }

  fprintf(fp,"%f %f\n",r0,i0);
 }

fclose(fp);

/* Display warning message to alert the user... */
if(inull)printf("bisp_to_eric/Warning: %d null values in bispectrum!\n",inull);

return(0);
}
/*******************************************************************
* BISPEC3 
* Same as bispec1, but newer version.
* Sum of full bispectrum terms (amplitude and phase)
* Reads a spectrum from an real image (from observations)
*
* Computes the bispectrum and spectrum lists from RE and IM
* already computed before calling this routine.
*
* Integrates the squared modulus of the spectrum and
* the phase factor of the bispectrum. 
* Does not correct from photon noise effects.
*
* Input:
* RE(NX,NY) (please note that idim has disappeared...)
* IM(NX,NY)
*
* Output:
* MODSQ: Sum of the modulus squared
* YCE1: phase factor of the bispectrum (real and imaginary)
* YCE1(.,2): sum of square real parts of the bispectrum
* YCE1(.,3): sum of square imaginary parts of the bispectrum
********************************************************************/
int BISPEC3(double *re, double *im, double *modsq, double *snrm,
            INT4 *nx, INT4 *ny, double *bispp, INT4 *ir, 
            INT4 *nbeta, INT4 *ngamma)
{
double w1;
INT4 k1, k2, k3;
double wr1, wr2, wr3, wi1, wi2, wi3;
double aar, aai, work1, work2; 
register INT4 ng, i, j, nb, yce_ng;
float xr[NBMAX], xi[NBMAX];

#ifdef DEBUG
/* Check the FFT (zero frequency is at (0,0)): */
w1 = re[0]*re[0] + im[0]*im[0];
w1 = sqrt(w1);
printf(" xphot1 = %f \n",w1);
#endif
 
/* IXY(1,NB) and IXY(2,NB) are the X,Y coordinates of the
* element NB of the spectral list. 
* As they (i.e. IXY) can be negative, and that zero frequency is at (0,0),
* we do the following transformation:
*/
  for(nb = 0; nb <= *nbeta; nb++)
     {
       i = ((ixy_1[nb] + *nx) % *nx); 
       j = ((ixy_2[nb] + *ny) % *ny); 
/* Debug:
       printf(" BISPEC3/debug: nb=%d ixy_1[nb]=%d ixy_2[nb]=%d i=%d j=%d\n", 
                nb, ixy_1[nb], ixy_2[nb], i, j);
*/
       xr[nb] = re[i + j * (*nx)];
       xi[nb] = im[i + j * (*nx)];
     }
 
/****************************************************************/
/* Phase factor of the bispectrum (with bispectral list):
* (No correction of photon noise)
*/
yce_ng = 0;
for(ng = 0; ng < *ngamma; ng++)
 {
     cover_klm0(&k1,1,ng);
     cover_klm0(&k2,2,ng);
     cover_klm0(&k3,3,ng);
/* debug:
     printf(" BISPEC3/debug: k1=%d k2=%d k3=%d \n",k1,k2,k3);
*/
     wr1 = xr[k1]; wr2 = xr[k2]; wr3 = xr[k3];
     wi1 = xi[k1]; wi2 = xi[k2]; wi3 = xi[k3];
/*  CC1=XC1*XC2*CONJG(XC3): AAR*WR3+AAI*WI3 */
     aar = wr1 * wr2 - wi1 * wi2;
     aai = wr2 * wi1 + wr1 * wi2;
     work1 = aar * wr3 + aai * wi3;
     work2 = aai * wr3 - aar * wi3;
     bispp[yce_ng] += work1;
     yce_ng++;
     bispp[yce_ng] += work2;

/* Estimation of the noise with the sum of squares: */
     yce_ng++;
     bispp[yce_ng] += work1*work1;
     yce_ng++;
     bispp[yce_ng] += work2*work2;
     yce_ng++;
  }
/* Debug:
     printf(" BISPEC3/debug: yce_ng max =%d\n",yce_ng-1);
*/
 
/****************************************************************/
/* Squared modulus for the output (not corrected for photon noise) */
    for( i = 0; i < (*nx) * (*ny); i++)
       {
        w1 = re[i]*re[i] + im[i]*im[i];
        modsq[i] += w1;
        snrm[i] += (w1 * w1);
       }
 
 return(0);
}
/*******************************************************************
* BISPEC3_SINGLE 
* Same as BISPEC3, but for single precision arrays 
********************************************************************/
int BISPEC3_SINGLE(float *re, float *im, float *modsq, float *snrm,
            INT4 *nx, INT4 *ny, float *bispp, INT4 *ir, 
            INT4 *nbeta, INT4 *ngamma)
{
float w1;
INT4 k1, k2, k3;
float wr1, wr2, wr3, wi1, wi2, wi3;
float aar, aai, work1, work2; 
register INT4 ng, i, j, nb, yce_ng;
float xr[NBMAX], xi[NBMAX];

/* IXY(1,NB) and IXY(2,NB) are the X,Y coordinates of the
* element NB of the spectral list. 
* As they (i.e. IXY) can be negative, and that zero frequency is at (0,0),
* we do the following transformation:
*/
  for(nb = 0; nb <= *nbeta; nb++)
     {
       i = ((ixy_1[nb] + *nx) % *nx); 
       j = ((ixy_2[nb] + *ny) % *ny); 
       xr[nb] = re[i + j * (*nx)];
       xi[nb] = im[i + j * (*nx)];
     }
 
/****************************************************************/
/* Phase factor of the bispectrum (with bispectral list):
* (No correction of photon noise)
*/
yce_ng = 0;
for(ng = 0; ng < *ngamma; ng++)
 {
     cover_klm0(&k1,1,ng);
     cover_klm0(&k2,2,ng);
     cover_klm0(&k3,3,ng);
/* DEBUG:
     if(ng <= 10) printf("ng=%d k1,k2,k3= %d %d %d\n", ng, k1,k2,k3);
*/
     wr1 = xr[k1]; wr2 = xr[k2]; wr3 = xr[k3];
     wi1 = xi[k1]; wi2 = xi[k2]; wi3 = xi[k3];
/*  CC1=XC1*XC2*CONJG(XC3): AAR*WR3+AAI*WI3 */
     aar = wr1 * wr2 - wi1 * wi2;
     aai = wr2 * wi1 + wr1 * wi2;
     work1 = aar * wr3 + aai * wi3;
     work2 = aai * wr3 - aar * wi3;
     bispp[yce_ng] += work1;
     yce_ng++;
     bispp[yce_ng] += work2;

/* Estimation of the noise with the sum of squares: */
     yce_ng++;
     bispp[yce_ng] += work1*work1;
     yce_ng++;
     bispp[yce_ng] += work2*work2;
     yce_ng++;
  }
 
/****************************************************************/
/* Squared modulus for the output (not corrected for photon noise) */
    for( i = 0; i < (*nx) * (*ny); i++)
       {
        w1 = re[i]*re[i] + im[i]*im[i];
        modsq[i] += w1;
        snrm[i] += (w1 * w1);
       }
 
 return(0);
}
/*******************************************************************
* bispec_1D_X
* (Fourier spectra ALONG THE LINES !).
*
* From BISPEC3, 1-D version for spectra
* Sum of full bispectrum terms (amplitude and phase)
* Reads a spectrum from an real image (from observations)
*
* Computes the bispectrum and spectrum lists from RE and IM
* already computed before calling this routine 
* (Fourier spectra ALONG THE LINES !).
*
* Integrates the squared modulus of the spectrum and
* the phase factor of the bispectrum. 
*
* Does not correct from photon noise effects. 
*
* Input:
* RE(NX,NY) (please note that idim has disappeared...)
* IM(NX,NY)
*
* Output:
* MODSQ: Sum of the modulus squared
* YCE1(.,0,iy) and YCE1(.,1,iy): phase factor of the bispectrum (real and imaginary)
* YCE1(.,2,iy): sum of square real parts of the bispectrum
* YCE1(.,3,iy): sum of square imaginary parts of the bispectrum
********************************************************************/
int bispec_1D_X(re,im,modsq,snrm,nx,ny,bispp,ir,nbeta,ngamma)
double re[], im[], modsq[], snrm[], bispp[];
INT4 *nx, *ny, *ir, *nbeta, *ngamma;
{
double w1;
INT4 k1, k2, k3;
double wr1, wr2, wr3, wi1, wi2, wi3;
double aar, aai, work1, work2; 
INT4 iix, ng, nb, iy_s, iy_b, bisp_dim;
register INT4 iy, i;
float xr[NBMAX], xi[NBMAX];

#ifdef DEBUG
/* Check the FFT in the middle line #(ny/2)
   (zero frequency is at (0,iy) for each line): */
i = 0 + (*ny / 2) * (*nx);
w1 = re[i]*re[i] + im[i]*im[i];
w1 = sqrt(w1);
printf(" xphot1(middle line) = %f \n",w1);
#endif

bisp_dim = *ngamma * 4;
for(iy = 0; iy < *ny; iy++)
  {
  iy_s = iy * (*nx);
  iy_b = iy * bisp_dim;
  for(nb = 0; nb <= *nbeta; nb++)
     {
/* Get coordinates (i,0) from nb index: */
       COVER_IXY_1D(&iix,&nb);
/* IXY(1,NB) and IXY(2,NB) are the X,Y coordinates of the
* element NB of the spectral list. 
* As they (i.e. IXY) can be negative, and that zero frequency is at (0,0),
* we do the following transformation:
*/
       iix = ((iix + *nx) % *nx); 
       xr[nb] = re[iix + iy_s];
       xi[nb] = im[iix + iy_s];
     }
 
/****************************************************************/
/* Phase factor of the bispectrum (with bispectral list):
* (No correction of photon noise)
*/
  for(ng = 0; ng < *ngamma; ng++)
    {
     cover_klm0(&k1,1,ng);
     cover_klm0(&k2,2,ng);
     cover_klm0(&k3,3,ng);
     wr1 = xr[k1]; wr2 = xr[k2]; wr3 = xr[k3];
     wi1 = xi[k1]; wi2 = xi[k2]; wi3 = xi[k3];
/*  CC1=XC1*XC2*CONJG(XC3): AAR*WR3+AAI*WI3 */
     aar = wr1 * wr2 - wi1 * wi2;
     aai = wr2 * wi1 + wr1 * wi2;
     work1 = aar * wr3 + aai * wi3;
     work2 = aai * wr3 - aar * wi3;
     bispp[ng + iy_b] += work1;
     bispp[ng + *ngamma + iy_b] += work2;

/* Estimation of the noise with the sum of squares: */
     bispp[ng + 2*(*ngamma) + iy_b] += work1*work1;
     bispp[ng + 3*(*ngamma) + iy_b] += work2*work2;
  }
 
/****************************************************************/
/* Squared modulus for the output (not corrected for photon noise) */
    for( i = 0; i < *nx; i++)
       {
        w1 = re[i + iy_s]*re[i + iy_s] + im[i + iy_s]*im[i + iy_s];
        modsq[i + iy_s] += w1;
        snrm[i + iy_s] += (w1 * w1);
       }
/* End of loop on iy: */
  }
 
 return(0);
}
/*******************************************************************
* bispec_1D_Y
* (Fourier spectra ALONG THE COLUMNS !).
*
* From BISPEC3, 1-D version for spectra
* Sum of full bispectrum terms (amplitude and phase)
* Reads a spectrum from an real image (from observations)
*
* Computes the bispectrum and spectrum lists from RE and IM
* already computed before calling this routine 
* (Fourier spectra ALONG THE COLUMNS !).
*
* Integrates the squared modulus of the spectrum and
* the phase factor of the bispectrum. 
*
* Does not correct from photon noise effects. 
*
* Input:
* RE(NX,NY) (please note that idim has disappeared...)
* IM(NX,NY)
*
* Output:
* MODSQ: Sum of the modulus squared
* YCE1(.,0,iy) and YCE1(.,1,iy): phase factor of the bispectrum (real and imaginary)
* YCE1(.,2,iy): sum of square real parts of the bispectrum
* YCE1(.,3,iy): sum of square imaginary parts of the bispectrum
********************************************************************/
int bispec_1D_Y(double *re, double *im, double *modsq, double *snrm,
                INT4 *nx, INT4 *ny, double *bispp, INT4 *ir, INT4 *nbeta, 
                INT4 *ngamma)
{
double w1;
INT4 k1, k2, k3;
double wr1, wr2, wr3, wi1, wi2, wi3;
double aar, aai, work1, work2; 
INT4 iiy, ng, nb, ix_b, bisp_dim;
register INT4 ix, i;
float xr[NBMAX], xi[NBMAX];

#ifdef DEBUG
/* Check the FFT in the middle column #(nx/2) 
   (zero frequency is at (ix,0) for each line): */
i = *nx /2 + 0 * (*nx);
w1 = re[i]*re[i] + im[i]*im[i];
w1 = sqrt(w1);
printf(" xphot1(middle line) = %f \n",w1);
#endif

bisp_dim = *ngamma * 4;
for(ix = 0; ix < *nx; ix++)
  {
  ix_b = ix * bisp_dim;
  for(nb = 0; nb <= *nbeta; nb++)
     {
/* Get coordinates (i,0) from nb index: */
       COVER_IXY_1D(&iiy,&nb);
/* IXY(1,NB) and IXY(2,NB) are the X,Y coordinates of the
* element NB of the spectral list. 
* As they (i.e. IXY) can be negative, and that zero frequency is at (0,0),
* we do the following transformation:
*/
       iiy = ((iiy + *ny) % *ny); 
       xr[nb] = re[ix + iiy*(*nx)];
       xi[nb] = im[ix + iiy*(*nx)];
     }
/* JLP999: same notation to be able to compare with "inv_bispec2_1D" */
if(ix == 0)
{
  for(i = 0; i < 3; i++)
  {
  w1 = sqrt(SQUARE(xr[i]) + SQUARE(xi[i]));
  printf("xc_re[%d]=%e xc_im[%d]=%e \n",i,xr[i]/w1,i,xi[i]/w1);
  }
}
 
/****************************************************************/
/* Phase factor of the bispectrum (with bispectral list):
* (No correction of photon noise)
*/
  for(ng = 0; ng < *ngamma; ng++)
    {
     cover_klm0(&k1,1,ng);
     cover_klm0(&k2,2,ng);
     cover_klm0(&k3,3,ng);
     wr1 = xr[k1]; wr2 = xr[k2]; wr3 = xr[k3];
     wi1 = xi[k1]; wi2 = xi[k2]; wi3 = xi[k3];
/*  CC1=XC1*XC2*CONJG(XC3): AAR*WR3+AAI*WI3 */
     aar = wr1 * wr2 - wi1 * wi2;
     aai = wr2 * wi1 + wr1 * wi2;
     work1 = aar * wr3 + aai * wi3;
     work2 = aai * wr3 - aar * wi3;
     bispp[ng + ix_b] += work1;
     bispp[ng + *ngamma + ix_b] += work2;
/* JLP999: same notation to be able to compare with "inv_bispec2_1D" */
    if(ix == 0 && ng < 3)
       {
        w1 = work1 * work1 + work2 * work2;
        w1 = sqrt(w1);
        if(w1 > 1.e-10)
        printf("line#%d: k=%d l=%d m=%d bisp[%d] = %e;%e or (%e;%e) modulus=%e\n",ix,k1,k2,k3,ng,work1,work2,work1/w1,work2/w1,w1);
        else
        printf("line#%d: k=%d l=%d m=%d bisp[%d] = %e;%e but modulus=%e\n",ix,k1,k2,k3,ng,work1,work2,w1);
       printf("bispp[1]= (%e,%e)\n",bispp[1],bispp[1 + *ngamma]);
       }

/* Estimation of the noise with the sum of squares: */
     bispp[ng + 2*(*ngamma) + ix_b] += work1*work1;
     bispp[ng + 3*(*ngamma) + ix_b] += work2*work2;
  }
 
/****************************************************************/
/* Squared modulus for the output (not corrected for photon noise) */
    for( i = 0; i < *ny; i++)
       {
        w1 = re[ix + i * (*nx)]*re[ix + i * (*nx)] + im[ix + i * (*nx)]*im[ix + i * (*nx)];
        modsq[ix + i * (*nx)] += w1;
        snrm[ix + i * (*nx)] += (w1 * w1);
       }
/* End of loop on ix: */
  }
 
/* JLP999 */
 printf("End: bispp[1]= (%e,%e)\n",bispp[1],bispp[1 + *ngamma]);
 return(0);
}
/*********************************************************************
*
* jlp_normalize_fft_1D
*********************************************************************/
int jlp_normalize_fft_1D(bispp,modsq,snrm,nx,ny,nbeta,ngamma)
double bispp[], modsq[], snrm[];
INT4 *nx, *ny, *nbeta, *ngamma;
{
double norm_fft, norm_fft2, norm_fft3, norm_fft4, norm_fft6;
INT4 bisp_dim, iy_s, iy_b;
register INT4 iy, i, ng;

norm_fft = sqrt((double)(*nx));
norm_fft2 = *nx;
norm_fft3 = norm_fft2 * norm_fft;
norm_fft4 = norm_fft2 * norm_fft2;
norm_fft6 = norm_fft3 * norm_fft3;

bisp_dim = *ngamma;

for(iy = 0; iy < *ny; iy++)
  {
    iy_s = iy * (*nx);
    iy_b = iy * bisp_dim * 4;
/* Normalizes FFT (for compatibility with FFT_2D instead of FFT_2D_FAST...*/
  for( i = 0; i < *nx; i++) 
     {
     modsq[i + iy_s] *= norm_fft2;
     snrm[i + iy_s] *= norm_fft4;
     }

   for(ng = 0; ng < *ngamma; ng++)
     {
     bispp[ng + iy_b] *= norm_fft3;
     bispp[ng + bisp_dim + iy_b] *= norm_fft3;
     bispp[ng + 2*bisp_dim + iy_b] *= norm_fft6;
     bispp[ng + 3*bisp_dim + iy_b] *= norm_fft6;
     }
  }

return(0);
}
/*********************************************************************
*
* jlp_normalize_fft
*********************************************************************/
int jlp_normalize_fft(bispp,modsq,snrm,nx,ny,nbeta,ngamma)
double bispp[], modsq[], snrm[];
INT4 *nx, *ny, *nbeta, *ngamma;
{
double norm_fft, norm_fft2, norm_fft3, norm_fft4, norm_fft6;
register INT4 i, ng, yce_ng;

norm_fft = sqrt((double)(*nx * *ny));
norm_fft2 = *nx * *ny;
norm_fft3 = norm_fft2 * norm_fft;
norm_fft4 = norm_fft2 * norm_fft2;
norm_fft6 = norm_fft3 * norm_fft3;

/* Normalizes FFT (for compatibility with FFT_2D instead of FFT_2D_FAST...*/
  for( i = 0; i < *ny * *nx; i++) 
     {
     modsq[i] *= norm_fft2;
     snrm[i] *= norm_fft4;
     }

   yce_ng = 0;
   for(ng = 0; ng < *ngamma; ng++)
     {
     bispp[yce_ng] *= norm_fft3;
     yce_ng++; 
     bispp[yce_ng] *= norm_fft3;
     yce_ng++; 
     bispp[yce_ng] *= norm_fft6;
     yce_ng++; 
     bispp[yce_ng] *= norm_fft6;
     yce_ng++; 
     }

return(0);
}
/**************************************************************** 
* Photon noise correction: (called by decode_car2 )
*
*  Input:
* modsq[]: mean normalized squared modulus
* xphot: mean photon flux per frame 
*
* Photon noise correction (cf JOSA 2, 14, Wirnitzer):
* (When spectrum normalized to one in the center)
* <i(u)>**2 = E(D...)/N**2 - 1/N
*
* <i(3)(u,v)> = E(D(3)(u,v)/N**3 - E(D(2)(u)/N**2)/N - E(D(2)(v)... +2/N**3)
*
*****************************************************************/
int photon_corr(double *bispp, double *modsq, double *snrm,
                INT4 *nx, INT4 *ny, float *xphot, INT4 *nbeta, INT4 *ngamma,
                INT4 *photon_correction)
{
double xphot2, xphot3, xphot4, xphot6;
double w2, work1;
float *xmodsq;
INT4 k1, k2, k3, iix, iiy;
register INT4 i, nb, ng, yce_ng;

/* Normalizes FFT (for compatibility with FFT_2D instead of FFT_2D_FAST...*/
jlp_normalize_fft(bispp,modsq,snrm,nx,ny,nbeta,ngamma);

/*******************************************************/
/* Return if no correction is needed: */
if(!*photon_correction) return(0);

/* Allocate memory for xmodsq */
xmodsq = (float *)malloc((*nbeta + 1) * sizeof(float));

/*******************************************************/
/* First the bispectrum */

  xphot2 = *xphot * *xphot;
/* IXY(1,NB) and IXY(2,NB) are the X,Y coordinates of the
* element NB of the spectral list.
* As they (i.e. IXY) can be negative, and that zero frequency is at (0,0),
* we do the following transformation:
*/
  for(nb = 0; nb <= *nbeta; nb++)
     {
       iix = ((ixy_1[nb] + *nx) % *nx);
       iiy = ((ixy_2[nb] + *ny) % *ny);
/* xmodsq is used to store the mean (not normalized) square modulus: */
       xmodsq[nb] = modsq[iix + iiy * *nx];
     }

/****************************************************************/
/* Phase factor of the bispectrum (with bispectral list):
* and correction from photon noise effects:
*/
 w2 = 2. *  *xphot;
 yce_ng = 0;
for(ng = 0; ng < *ngamma; ng++)
 {
     cover_klm0(&k1,1,ng);
     cover_klm0(&k2,2,ng);
     cover_klm0(&k3,3,ng);
/* Photon noise correction 
*/
     work1 = - (xmodsq[k1] + xmodsq[k2] + xmodsq[k3]) + w2;
     bispp[yce_ng] += work1;

/* Estimation of the noise with the sum of squares: */
     yce_ng++; yce_ng++;
     bispp[yce_ng] -= work1*work1;
     yce_ng++; yce_ng++;
  }

/*******************************************************/
/* Then correcting the squared modulus: */
  xphot4 = xphot2 * xphot2;
  for( i = 0; i < *ny * *nx; i++) 
     {
/* Photon noise correction: (biased_sq = xphot + xphot_sq * unbiased_sq)*/
      modsq[i] -= *xphot;
/* Normalization for compatibility with previous computations: */
      modsq[i] /= xphot2;
      snrm[i] /= xphot4;
    }

/*******************************************************/
xphot3 = *xphot * *xphot * *xphot;
xphot6 = xphot3 * xphot3;
yce_ng = 0;
for(ng = 0; ng < *ngamma; ng++)
 {
     bispp[yce_ng] /= xphot3;
     yce_ng++;
     bispp[yce_ng] /= xphot3;
     yce_ng++;
     bispp[yce_ng] /= xphot6;
     yce_ng++;
     bispp[yce_ng] /= xphot6;
     yce_ng++;
 }

free(xmodsq);
return(0);
}
/*******************************************************/
int rearrange_mask(bispp,ngamma,ydim1)
double bispp[];
INT4 *ngamma, *ydim1;
{
INT4 isize, ydim;
register INT4 ng, yce_ng;
double *work;

ydim = *ydim1;
isize = 4 * ydim * sizeof(double);
work = (double*) malloc(isize);

/* Erases work array: */
for(ng = 0; ng < 4 * ydim; ng++) work[ng] = 0.; 

/* Transfer to temporary array: */
yce_ng = 0;
for(ng = 0; ng < *ngamma; ng++)
 {
     work[ng] = bispp[yce_ng];
     yce_ng++;
     work[ng + ydim] = bispp[yce_ng];
     yce_ng++;
     work[ng + 2 * ydim] = bispp[yce_ng];
     yce_ng++;
     work[ng + 3 * ydim] = bispp[yce_ng];
     yce_ng++;
 }

/* Transfer back to bispp: */
for(ng = 0; ng < 4 * ydim; ng++) bispp[ng] = work[ng]; 

free(work);
return(0);
}
/**************************************************************** 
* Photon noise correction (Called by main program "phot_noise_mask.c"):
*
*  Input:
* modsq[]: mean normalized squared modulus
* phot_modsq[]: power spectrum of photon response  
* xphot: mean photon flux per frame 
* bispp[]: bispectrum (sorted in a "standard" way, and thus different
*         from what is needed by photon_corr...)
*
* Photon noise correction (cf JOSA 2, 14, Wirnitzer):
* (When spectrum normalized to one in the center)
* <i(u)>**2 = E(D...)/N**2 - 1/N
*
* <i(3)(u,v)> = E(D(3)(u,v)/N**3 - E(D(2)(u)/N**2)/N - E(D(2)(v)... +2/N**3)
* fraction: gain correction to be applied to phot_modsq 
*****************************************************************/
int photon_corr_mask(double *bispp, double *modsq, double *phot_modsq,
                     INT4 *nx, INT4 *ny, INT4 *bisp_dim, float *xphot,
                     INT4 *nbeta, INT4 *ngamma, float *fraction)
{
float xphot2, *phot_re, *xmodsq;
double w2, work, work1, epsilon = 1.e-8;
INT4 k1, k2, k3, iix, iiy, nb;
register INT4 i, ng;

if(fraction == 0) 
             {
             printf(" photon_corr_mask/fatal error: fraction=0!\n");
             exit(-1);
             }

 for( i = 0; i < (*ny) * (*nx); i++) phot_modsq[i] *= (*fraction);

/*******************************************************/
/* first division by the power spectrum phot_modsq: */
  for( i = 0; i < (*ny) * (*nx); i++) 
   if(phot_modsq[i] > epsilon)
         modsq[i] /= (phot_modsq[i] + epsilon);
   else
         modsq[i] = 0.;

if((phot_re = (float*)malloc((*nbeta + 1) * sizeof(float))) == NULL)
 {
 printf("photon_corr_mask/fatal error, allocating memory space \n");
 exit(-1);
 };

/* Allocate memory for xmodsq */
xmodsq = (float *)malloc((*nbeta + 1) * sizeof(float));

/*******************************************************/
/* First the bispectrum, then the power spectrum */

  xphot2 = *xphot * *xphot;
/* IXY(1,NB) and IXY(2,NB) are the X,Y coordinates of the
* element NB of the spectral list.
* As they (i.e. IXY) can be negative, and that zero frequency is at (0,0),
* we do the following transformation:
*/
  for(nb = 0; nb <= *nbeta; nb++)
     {
/* before 2008: 
       iix = ((ixy_1[nb] + *nx) % *nx);
       iiy = ((ixy_2[nb] + *ny) % *ny);
*/
/* Get coordinates (i,j) from nb index: */
       COVER_IXY(&iix, &iiy, &nb);
       iix = ((iix + *nx) % *nx);
       iiy = ((iiy + *ny) % *ny);
/* xmodsq is used to store the mean (not normalized) square modulus: */
       xmodsq[nb] = modsq[iix + iiy * *nx];
/* The photon response is symmetric and real, hence its FT is real: */
       phot_re[nb] = sqrt(phot_modsq[iix + iiy * *nx]);
#ifdef DEBUG
       if(nb < 4) printf(" nb = %d iix, iiy: %d %d \n",nb,iix,iiy);
       if(nb < 4) printf(" xmodsq[%d] = %f \n",nb,xmodsq[nb]);
#endif
     }

/* JLP96: */
/* Number of photons:  xphot**2 + xphot - modsq_0 = 0
* Delta = b**2 - 4 a*c = 1 + 4 * modsq_0
* Solutions: (- b +/- sqrt(Delta))/2a
*  i.e.,     (-1 + sqrt(1 + 4 * modsq_0) )/2
*/
  work = (-1. + sqrt((double)(1. + 4. * xmodsq[0])))/2.;
  printf(" My estimate of xphot (from modsq[0]) is: %f, whereas yours is: %f \n",
           work,*xphot);
/* End of JLP96. */
  printf("Correcting bispectrum with fraction=%f of square modulus \n",*fraction);

/****************************************************************/
/* Phase factor of the bispectrum (with bispectral list):
* and correction from photon noise effects:
*/
 w2 = 2. *  *xphot;
for(ng = 0; ng < *ngamma; ng++)
 {
     cover_klm0(&k1,1,ng);
     cover_klm0(&k2,2,ng);
     cover_klm0(&k3,3,ng);
/* JLP2000: division by the spectrum of the photon response: */
     work1 = phot_re[k1] * phot_re[k2] * phot_re[k3] + epsilon;
     bispp[ng] /= work1;
     bispp[ng + *bisp_dim] /= work1;

#ifdef DEBUG
     if(ng < 4) printf(" ng = %d k, l, m: %d %d %d \n",ng,k1,k2,k3);
     if(ng < 4) printf(" xmodsq1, xmodsq2, xmodsq3 %f %f %f \n",
                       xmodsq[k1],xmodsq[k2],xmodsq[k3]);
#endif
/* Photon noise correction 
*/
     work1 = - (xmodsq[k1] + xmodsq[k2] + xmodsq[k3]) * (*fraction) + w2;
#ifdef DEBUG
     if(ng < 4) printf(" ng= %d corr= %f \n",ng,work1);
#endif
     bispp[ng] += work1;
/* Phase factor: */
/* JLP96: Do not normalize it (to be able to visualize it with bisp_to_image)
     w1 = bispp[ng]*bispp[ng] + bispp[ng + *bisp_dim] * bispp[ng + *bisp_dim];
     if(w1 > 0) 
       {
       w1 = sqrt(w1);
       bispp[ng] /= w1;
       bispp[ng + *bisp_dim] /= w1;
       }
*/

  }
/*******************************************************/
/* Then correcting the squared modulus: */
/* (Normalization to obtain unity in the center (zero frequency is at 0,0 here))
  w2 = modsq[0] - *xphot;
 */
  for( i = 0; i < *ny * *nx; i++)
     {
/* Photon noise correction: (biased_sq = xphot + xphot_sq * unbiased_sq)*/
      modsq[i] -= *xphot;
/* Normalization for compatibility with previous computations: */
/* and get xphot in the central frequency: (instead of xphot**2): */
/* JLP96: Do not normalize it (to be able to visualize it with bisp_to_image)
      modsq[i] /= w2;
*/
    }

free(phot_re);
free(xmodsq);
return(0);
}
/*********************************************************************
* Correction of the sky background of long integration)
*
*********************************************************************/
int corr_bisp_sky(double *bispp, double *modsq, INT4 *nx, INT4 *ny,
                  INT4 *bisp_dim, float *sky_level, INT4 *nbeta, INT4 *ngamma)
{
float *xmodsq;
double w1, w2;
INT4 k1, k2, k3, iix, iiy;
register INT4 nb, ng;

/* Allocate memory for xmodsq */
xmodsq = (float *)malloc((*nbeta + 1) * sizeof(float));

/*******************************************************/
/* First the bispectrum */

/* IXY(1,NB) and IXY(2,NB) are the X,Y coordinates of the
* element NB of the spectral list.
* As they (i.e. IXY) can be negative, and that zero frequency is at (0,0),
* we do the following transformation:
*/
  for(nb = 0; nb <= *nbeta; nb++)
     {
       iix = ((ixy_1[nb] + *nx) % *nx);
       iiy = ((ixy_2[nb] + *ny) % *ny);
/* xmodsq is used to store the mean (not normalized) square modulus: */
       xmodsq[nb] = modsq[iix + iiy * *nx];
#define DEBUG
#ifdef DEBUG
       if(nb < 4) printf(" nb = %d iix, iiy: %d %d \n",nb,iix,iiy);
       if(nb < 4) printf(" xmodsq = %f \n",xmodsq[nb]);
#endif
     }

/****************************************************************/
/* Phase factor of the bispectrum (with bispectral list):
* and correction from photon noise effects:
*/
 w2 = (*sky_level) * sqrt((double)((*nx) * (*ny)));
for(ng = 0; ng < *ngamma; ng++)
 {
     cover_klm0(&k1,1,ng);
     cover_klm0(&k2,2,ng);
     cover_klm0(&k3,3,ng);
/*
     if(ng < 4) printf(" ng = %d k, l, m: %d %d %d \n",ng,k1,k2,k3);
     if(ng < 4) printf(" xmodsq1, xmodsq2, xmodsq3 %f %f %f \n",
                       xmodsq[k1],xmodsq[k2],xmodsq[k3]);
*/
     if(k1 == 0 || k2 == 0 || k3==0) printf(" ng = %d k, l, m: %d %d %d \n",ng,k1,k2,k3);
#ifdef TTR
/* Photon noise correction 
*/
     work1 = - (xmodsq[k1] + xmodsq[k2] + xmodsq[k3]) + w2;
#ifdef DEBUG
     if(ng < 4) printf(" ng= %d corr= %f \n",ng,work1);
#endif
     bispp[ng] += work1;
/* Phase factor: */
/* JLP96: Do not normalize it (to be able to visualize it with bisp_to_image)
     w1 = bispp[ng]*bispp[ng] + bispp[ng + *bisp_dim] * bispp[ng + *bisp_dim];
     if(w1 > 0) 
       {
       w1 = sqrt(w1);
       bispp[ng] /= w1;
       bispp[ng + *bisp_dim] /= w1;
       }
*/
#endif

  }
/*******************************************************/
/* Then correcting the squared modulus: */
/* At the center: |H(0)|^2 \times |O(0)|^2  + 2c \times H(0) \times O(0) + c^2
* Hence: y = x2 + 2cx + c2 = (x + c)^2
* x = - c + sqrt(y) 
*/
      w1 = sqrt((double)modsq[0]); 
/* Fourier transform normalized by sqrt(nx * ny), hence only sqrt(nx * ny)
* otherwise if not normalized, it would be nx * ny: */
      w2 = (*sky_level) * sqrt((double)((*nx) * (*ny)));
      printf(" w1 = %e w2 = %e modsq[0] = %e \n",w1,w2,modsq[0]);
      modsq[0] = (w1 - w2)*(w1-w2);

free(xmodsq);
return(0);
}
/**************************************************************** 
* Calibration by reference star
* (Called by main program "phot_noise_ma2.c"):
*
*  Input:
* modsq[]: mean normalized squared modulus
* xphot: mean photon flux per frame 
* bispp[]: bispectrum (sorted in a "standard" way, and thus different
*         from what is needed by photon_corr...)
*
*****************************************************************/
int calib_ref(bispp,modsq,bispp_ref,modsq_ref,nx,ny,bisp_dim,nbeta,ngamma)
double bispp[], modsq[], bispp_ref[], modsq_ref[]; 
INT4 *nx, *ny, *bisp_dim, *nbeta, *ngamma;
{
double w1, w2, work;
double cos1, sin1;
register INT4 ng;

/*******************************************************/
/* First the bispectrum */
/****************************************************************/
/* Phase factor of the bispectrum (with bispectral list):
*/
for(ng = 0; ng < *ngamma; ng++)
  {
  cos1 = bispp_ref[ng]; 
  sin1 = bispp_ref[ng + (*ngamma)]; 
  work = cos1 * cos1 + sin1 * sin1;
  work = sqrt(work);
  if(work == 0.)
    {
    bispp[ng] = 0.;
    }
  else
    {
    cos1 /= work;
    sin1 /= work;
    w1 = bispp[ng];
    w2 = bispp[ng + (*ngamma)];
    bispp[ng] = w1 * cos1 + w2 * sin1; 
    bispp[ng + (*ngamma)] = w2 * cos1 - w1 * sin1; 
    }
  }
/*******************************************************/
/* Squared modulus: */
/* Possibility of simple deconvolution:
for( i = 0; i < *ny * *nx; i++)
   {
   work = modsq_ref[i];
   if(work < 1.e-10)
      modsq[i] = 0.;
   else
      modsq[i] /= work;
   }
*/

return(0);
}
/**************************************************************** 
* Photon noise correction for slit spectra (Called by phot_noise_1D):
*
*  Input:
* modsq[]: mean (not normalized) squared modulus
* xphot: mean photon flux per frame 
* bispp[]: bispectrum (sorted in a "standard" way, and thus different
*         from what is needed by photon_corr...)
*
* Photon noise correction (cf JOSA 2, 14, Wirnitzer):
* (When spectrum normalized to one in the center)
* <i(u)>**2 = E(D...)/N**2 - 1/N
* <i(3)(u,v)> = E(D(3)(u,v)/N**3 - E(D(2)(u)/N**2)/N 
*               - E(D(2)(v)/N**2)/N -E(D(2)(-u-v)/N**2)/N +2/N**3)
* When it is not normalized:
* (cf Wirnitzer, JOSA A 1985, p16)
* Q3 = D3(u,v) - (D2(u) + D2(v) + D2(-u-v) - 2N)
* where D2= N + N**2 <|i(u)|**2>
*****************************************************************/
int photon_corr_1D(float *bispp, float *modsq, INT4 *nx, INT4 *ny, 
                   INT4 *bisp_dim, INT4 *nbeta, INT4 *ngamma)
/*
* WARNING: single precision !!!!! 
*/
{
float xphot, modsq_0, *xmodsq;
double w2, work1;
INT4 k1, k2, k3, iix, nb, ng, iy_b,iy_s;
register INT4 iy, i;

/* Allocate memory for xmodsq */
xmodsq = (float *)malloc((*nbeta + 1) * sizeof(float));

/*******************************************************/
/* First the bispectrum */

for(iy = 0; iy < *ny; iy++)
  { 
/* Index of bispectrum list: */
   iy_b = iy * (*bisp_dim);
/* Index of spectrum line: */
   iy_s = iy * (*nx); 
/* Number of photons:  N**2 + N at zero frequency 
*  (i.e., modsq_0=N**2+N) */
   modsq_0 = modsq[0+iy_s];
/* Case of empty line: */
   if(modsq_0 < 0) 
     {
     modsq[0+iy_s] = 0.;
/* Set modulus to zero and leave bispectrum unchanged */
     }
/* Case of non-empty line: */
   else
     {
/* Number of photons:  xphot**2 + xphot - modsq_0 = 0
* Delta = b**2 - 4 a*c = 1 + 4 * modsq_0
* Solutions: (- b +/- sqrt(Delta))/2a
*  i.e.,     (-1 + sqrt(1 + 4 * modsq_0) )/2
*/
     xphot = (-1. + sqrt((double)(1. + 4. * modsq_0)))/2.;
     for(nb = 0; nb <= *nbeta; nb++)
       {
/* Get coordinates (iix,iiy) from nb index: */
       COVER_IXY_1D(&iix,&nb);
/* IXY(1,NB) and IXY(2,NB) are the X,Y coordinates of the
* element NB of the spectral list. 
* As they (i.e. IXY) can be negative, and that zero frequency is at (0,0),
* we do the following transformation:
* (JLP96: transformation checked, OK)
*/
       iix = ((iix + *nx) % *nx); 
/* xmodsq is used to store the mean (not normalized) square modulus: */
       xmodsq[nb] = modsq[iix + iy_s];
       }

/****************************************************************/
/* Phase factor of the bispectrum (with bispectral list):
* and correction from photon noise effects:
*/
      w2 = 2. *  xphot;
      for(ng = 0; ng < *ngamma; ng++)
        {
        cover_klm0(&k1,1,ng);
        cover_klm0(&k2,2,ng);
        cover_klm0(&k3,3,ng);
/* Photon noise correction 
*/
        work1 = - (xmodsq[k1] + xmodsq[k2] + xmodsq[k3]) + w2;
        bispp[ng + iy_b] += work1;
/* To get the phase factor, need normalization : */
/* Commented out:
        w1 = bispp[ng + iy_b]*bispp[ng + iy_b] 
         + bispp[ng + *bisp_dim + iy_b] * bispp[ng + *bisp_dim + iy_b];
*/
/* Normalization: (no longer used now, JLP95)*/
/* Commented out:
        if(w1 > 0) 
         {
         w1 = sqrt(w1);
         if(w1 < 1.e-10)w1=1.e-10;
         bispp[ng + iy_b] /= w1;
         bispp[ng + *bisp_dim + iy_b] /= w1;
         }
*/
/* End of case: modsq > 0 */
    }

/* End of loop on ng: */
  }

/*******************************************************/
/* Then correcting the squared modulus: */
/* (Normalization to obtain unity in the center (zero frequency is at 0,y here)) */
     for( i = 0; i < *nx; i++) 
      {
/* Photon noise correction: (biased_sq = xphot + xphot_sq * unbiased_sq)*/
      modsq[i + iy_s] -= xphot;
      }
/* End of loop on iy */
 }

free(xmodsq);
return(0);
}
/***************************************************************
*
* Compute bispectrum slice:
* (Called by bisp_to_image.c)
*
* bisp(u,v) = bisp((u_1,u_2),(v_1,v_2)) 
* bisp(u,v) = spec(u) * spec(v) * spec(-u-v) 
* bisp(u,v) = spec(u_1,u_2) * spec(v_1,v_2) * spec(-u_1-v_1,-u_2-v_2) 
* Here we take u_2 = v_2 = 0
* bisp(u,v) = spec(u_1,0) * spec(v_1,0) * spec(-u_1-v_1,0) 
* Hence we display the image bisp(u_1,0,v_1,0), 
* which is a function of the two coordinates u_1 and v_1
*
* ShouldNormalize: flag set to 1 if modsq and bisp have to be normalized
* iopt: option set to 0 if real part of bispectrum has to be output,
*                     1 if imaginary part ...
*                     2 if SNR ...
*                     3 if phase
*                     4 if modulus i.e., sqrt(real^2 + imag^2) 
****************************************************************/
int bisp2D_to_2D_image(float *out_image, float *modsq, INT4 *nx, INT4 *ny, 
                       float *bisp, INT4 *nx_bisp, INT4 *nbeta, INT4 *ngamma, 
                       INT4 *ShouldNormalize, INT4 *iopt)
/* WARNING: Single precision for modsq and bisp arrays!!! */
{
double cos1, sin1, modulus;
float work, corr_factor, angle1;
register INT4 i, ng, ix, iy;
INT4 icent, kk, ll, mm, ix_cent, iy_cent;

/* If you want to check that all pixels are scanned in this routine: 
replace 0 by any arbitrary number:  */
/* JLP2008: some pixels are not set because the image has even pixels
* and symmetry cannot reach them... */
for(i = 0; i < (*nx)*(*ny); i++) out_image[i] = 0.;

  if(*nx != *ny){
  fprintf(stderr, 
  "bisp2D_to_2D_image/Fatal error, the image should be square: nx=%d ny=%d!\n", 
           *nx, *ny);
  exit(-1);
  }

ix_cent = *nx / 2;
iy_cent = *ny / 2;
/* Coordinate of the center: */
icent = (*nx)/2 + (*nx) * (*ny)/2; 

/* Fill main array: 
* iopt = 3, i.e., phase 
* iopt = 4, i.e., modulus
*/
if(*iopt == 3 || *iopt == 4)
 {
  for (ng = 0; ng < *ngamma; ng++)
    {
/* Spectral indices k, l, m, of the three components of bispectral term:
*/
     cover_klm0(&kk,1,ng);
     cover_klm0(&ll,2,ng);
     cover_klm0(&mm,3,ng);
/* Load Bisp[ u(kk,0), v(ll, 0) ]: */
     if(ixy_2[kk] == 0 && ixy_2[ll] == 0)
       {
         ix = ixy_1[kk];
         iy = ixy_1[ll];
         modulus = SQUARE(bisp[ng]) + SQUARE(bisp[ng+*ngamma]); 
         if(modulus == 0.)
            {
            angle1 = 0.; 
            }
         else
            {
            modulus = sqrt(modulus);
            cos1 = bisp[ng] / modulus;
            sin1 = bisp[ng + *ngamma] / modulus;
            angle1 = acos(cos1); 
            if(sin1 < 0.) angle1 = -angle1;
            }
         if(*iopt == 3) 
            out_image[ix + icent + (*nx) * iy] = angle1; 
         else
            out_image[ix + icent + (*nx) * iy] = modulus; 
       }
    }
 }
/* Fill main array, general case: iopt=0 (real), 1(imag), 2(snr) */
else
 {
  for (ng = 0; ng < *ngamma; ng++)
    {
/* Spectral indices k, l, m, of the three components of bispectral term:
*/
     cover_klm0(&kk,1,ng);
     cover_klm0(&ll,2,ng);
     cover_klm0(&mm,3,ng);
     if(ixy_2[kk] == 0 && ixy_2[ll] == 0)
       {
         ix = ixy_1[kk];
         iy = ixy_1[ll];
/* Load Bisp[ u(kk,0), v(0,ll) ]: */
         out_image[ix + icent + (*nx) * iy] = 
                              bisp[ng + (*ngamma)*(*iopt)];
       }
    }
/* end of iopt != 3 */
 }

/* Fill central line with 
* modsq, for since bisp(u,0)=modsq(u) for iopt=0 (real part), iopt=2 (snr) 
* or for iopt=4 (real^2 + imag^2)
* zero for iopt=1, (imaginary part) or iopt=3 (phase) */
  if(*iopt == 0 || *iopt == 2 || *iopt == 4)
   {

/* X and Y axes correspond to the normalized power spectrum: 
* (They are equal since B(u,v)=B(v,u) ...)
* X axis, i.e., v_2=0:
* bisp(u,v) = spec(u_1,0) * spec(0,0) * spec(-u_1,0) = |spec(u_1,0)|**2 
* Y axis, i.e., u_2=0:
* bisp(u,v) = spec(0,0) * spec(v_1,0) * spec(-v_1,0) = |spec(v_1,0)|**2 
*/

/* Upper part of Y axis:  (at the centre of out_image) */
/* JLP2008: I multiply with the value of the spectrum for zero frequency: */
  corr_factor = sqrt(modsq[icent]);
  for(ix = ix_cent; ix < (*nx); ix++) 
    out_image[ix_cent + (*nx) * ix] = 
                     corr_factor * modsq[ix + (*nx) * iy_cent];

  }
/* iopt = 1 (imaginary part) or iopt = 3 (phase) */
  else
  {
  for(ix = ix_cent; ix < *nx; ix++) out_image[ix_cent + (*nx) * ix] = 0.; 
  }

/* Now the upper left triangle in the upper right quadrant is filled */ 
/* Fill remaining half (lower right triangle) of the upper right quadrant with
*  axial symmetry: bisp(u,v)=bisp(v,u) */
  for(ix = 0; ix < ix_cent; ix++)
    for(iy = 0; iy < ix; iy++)
       {
       out_image[ix + ix_cent + (iy + iy_cent) * (*nx)] =
                  out_image[iy + ix_cent + (ix + iy_cent) * (*nx)];
       }

/* Fill the lower left quadrant with central symmetry:
* bisp(u,v)=conjugate of bisp(-u,-v) (since image values are real values)*/
  for(iy = 0; iy < iy_cent; iy++)
    for(ix = 0; ix < ix_cent; ix++)
       {
/* iopt=1 imag, iopt=3 phase: */
       if(*iopt == 1 || *iopt == 3)
       out_image[-ix + ix_cent + (-iy + iy_cent) * (*nx)] =
                  -out_image[ix + ix_cent + (iy + iy_cent) * (*nx)];
/* iopt=0 real, iopt=2 snrm, iopt=4 modulus: */
       else
       out_image[-ix + ix_cent + (-iy + iy_cent) * (*nx)] =
                  out_image[ix + ix_cent + (iy + iy_cent) * (*nx)];
       }

/* Fill half of the lower right quadrant with "special" symmetry:
* bisp(u,v)=s(u)s(v)s(-u-v)= s(u)s(-u-v)s(-u-(-u-v))=bisp(u,-u-v) */
  for(iy = 0; iy < iy_cent; iy++)
    for(ix = 0; ix < ix_cent - iy; ix++)
       {
       out_image[ix + ix_cent + (-ix-iy + iy_cent) * (*nx)] =
                  out_image[ix + ix_cent + (iy + iy_cent) * (*nx)];
       }

/* Again, use bisp(-u,v)=conjugate of bisp(u,-v) (since image values 
* are real values)
* to fill half of the upper left quadrant with central symmetry: */
  for(iy = 0; iy < iy_cent; iy++)
    for(ix = 0; ix < iy; ix++)
       {
/* iopt=1 imag, iopt=3 phase: */
       if(*iopt == 1 || *iopt == 3)
       out_image[-ix + ix_cent + (iy + iy_cent) * (*nx)] =
                  -out_image[ix + ix_cent + (-iy + iy_cent) * (*nx)];
/* iopt=0 real, iopt=2 snrm, iopt=4 modulus: */
       else
       out_image[-ix + ix_cent + (iy + iy_cent) * (*nx)] =
                  out_image[ix + ix_cent + (-iy + iy_cent) * (*nx)];
       }

/* Fill remaining half of the lower right quadrant with axial symmetry:
* bisp(v,-u)=bisp(-u,v) */
  for(iy = 0; iy < iy_cent; iy++)
    for(ix = 0; ix < iy; ix++)
       {
       out_image[iy + ix_cent + (-ix + iy_cent) * (*nx)] =
                  out_image[-ix + ix_cent + (iy + iy_cent) * (*nx)];
       }

/* Fill remaining half of the upper left quadrant with central symmetry:
* bisp(-v,u)=bisp(u,-v) */
  for(iy = 0; iy < iy_cent; iy++)
    for(ix = iy; ix < ix_cent; ix++)
       {
       if(*iopt == 1 || *iopt == 3)
       out_image[-ix + ix_cent + (iy + iy_cent) * (*nx)] =
                  -out_image[iy + ix_cent + (-ix + iy_cent) * (*nx)];
       else
       out_image[-ix + ix_cent + (iy + iy_cent) * (*nx)] =
                  out_image[iy + ix_cent + (-ix + iy_cent) * (*nx)];
       }

/* iopt=0 real, iopt=2 snrm, iopt=4 modulus: */
work = out_image[icent];
if(*ShouldNormalize && (*iopt == 0 || *iopt == 2 || *iopt == 4) && work != 0.){
      printf("Bisp_2D_to_2D_image/Normalizing output image with scale=%f\n", 
              work);
       for(i = 0; i < (*nx)*(*ny); i++) out_image[i] /= work;
      }

return(0);
}
/***********************************************************************
* Display bispectra for 1D-spectra (slit spectra)
* for a given line (line_to_output)
*
* iopt: option set to 0 if real part of bispectrum has to be output,
*                     1 if imaginary part ...
*                     2 if SNR ...
*                     3 if phase
***********************************************************************/
int bisp1D_to_2D_image(float *out_image, INT4 *nx_out, INT4 *ny_out, 
                       float *modsq, INT4 *nx_modsq, float *bisp, INT4 *nx_bisp,
                       INT4 *nbeta, INT4 *ngamma, INT4 *line_to_output, 
                       INT4 *normalize, INT4 *iopt)
/* WARNING: Single precision for modsq and bisp arrays!!! */
{
float corr_factor;
INT4 ix, iy, iy_b, ix_cent, iy_cent, iix, k1, k2, ng;

/* JLP96: */
 if(*iopt == 3) {printf("bisp1D_to_2Dimage/error, option not yet available\n");
                exit(-1);
               }

ix_cent = (*nx_out)/2;
iy_cent = (*ny_out)/2;
iy_b = (*line_to_output) * (*nx_bisp);

/* Display bispectrum of (u_1,0,v_1,0): */
  for(ng = 0; ng < *ngamma; ng++)
    {
/* Get nb spectral index (k1) for the first vector 
* of ng index bispectral list: */
     cover_klm0(&k1,1,ng);
/* Get coordinates (iix,0) from nb index: */
     COVER_IXY_1D(&iix,&k1);
     ix = ((iix + *nx_modsq) % *nx_modsq); 
/* Same for second vector : */
     cover_klm0(&k2,2,ng);
     COVER_IXY_1D(&iix,&k2);
     iy = ((iix + *nx_modsq) % *nx_modsq); 
/* I add (ix_cent,iy_cent) to get (0,0) bifrequency in the center: */
     out_image[ix + ix_cent + (iy + iy_cent) * (*nx_out)] = 
                  bisp[ng + (*iopt) * (*ngamma) + iy_b];
    }

/* Fill central line with modsq, since bisp(u,0)=modsq(u): 
* (except when iopt = 1, i.e., imaginary part is wanted) */
  if(*iopt != 1)
   {
/* For some unknown reason, modsq should be multiplied by 2 (symmetry ???): 
   and snrm by 4*/
   if(*iopt == 0)
     corr_factor = 2.;
   else
     corr_factor = 4.;
   for(iy = 0; iy < iy_cent; iy++) 
     out_image[ix_cent + (iy + iy_cent) * (*nx_out)] = 
      corr_factor * modsq[(*nx_modsq)/2 + iy + (*line_to_output) * (*nx_modsq)];
   }
/* iopt = 1: imaginary part is wanted */
  else
   {
   for(iy = 0; iy < iy_cent; iy++) 
       out_image[ix_cent + (iy + iy_cent) * (*nx_out)] = 0.; 
   }

/* Fill remaining half of the upper right quadrant with
*  axial symmetry: bisp(u,v)=bisp(v,u) */
  for(iy = 0; iy < iy_cent; iy++)
    for(ix = 0; ix < iy; ix++)
       {
       out_image[iy + ix_cent + (ix + iy_cent) * (*nx_out)] = 
                  out_image[ix + ix_cent + (iy + iy_cent) * (*nx_out)];
       }

/* Fill the lower left quadrant with central symmetry: 
* bisp(u,v)=conjugate of bisp(-u,-v) (since image values are real values)*/
  for(iy = 0; iy < iy_cent; iy++)
    for(ix = 0; ix < ix_cent; ix++)
       {
       if(*iopt == 1)
       out_image[-ix + ix_cent + (-iy + iy_cent) * (*nx_out)] = 
                  -out_image[ix + ix_cent + (iy + iy_cent) * (*nx_out)];
       else
       out_image[-ix + ix_cent + (-iy + iy_cent) * (*nx_out)] = 
                  out_image[ix + ix_cent + (iy + iy_cent) * (*nx_out)];
       }

/* Fill half of the lower right quadrant with "special" symmetry: 
* bisp(u,v)=s(u)s(v)s(-u-v)= s(u)s(-u-v)s(-u-(-u-v))=bisp(u,-u-v) */
  for(iy = 0; iy < iy_cent; iy++)
    for(ix = 0; ix < ix_cent - iy; ix++)
       {
       out_image[ix + ix_cent + (-ix-iy + iy_cent) * (*nx_out)] = 
                  out_image[ix + ix_cent + (iy + iy_cent) * (*nx_out)];
       }

/* Fill half of the upper left quadrant with central symmetry: 
* bisp(-u,v)=conjugate of bisp(u,-v) (since image values are real values)*/
  for(iy = 0; iy < iy_cent; iy++)
    for(ix = 0; ix < iy; ix++)
       {
       if(*iopt == 1)
       out_image[-ix + ix_cent + (iy + iy_cent) * (*nx_out)] = 
                  -out_image[ix + ix_cent + (-iy + iy_cent) * (*nx_out)];
       else
       out_image[-ix + ix_cent + (iy + iy_cent) * (*nx_out)] = 
                  out_image[ix + ix_cent + (-iy + iy_cent) * (*nx_out)];
       }

/* Fill remaining half of the lower right quadrant with axial symmetry: 
* bisp(v,-u)=bisp(-u,v) */
  for(iy = 0; iy < iy_cent; iy++)
    for(ix = 0; ix < iy; ix++)
       {
       out_image[iy + ix_cent + (-ix + iy_cent) * (*nx_out)] = 
                  out_image[-ix + ix_cent + (iy + iy_cent) * (*nx_out)];
       }

/* Fill remaining half of the upper left quadrant with central symmetry: 
* bisp(-v,u)=bisp(u,-v) */
  for(iy = 0; iy < iy_cent; iy++)
    for(ix = iy; ix < ix_cent; ix++)
       {
       if(*iopt == 1)
       out_image[-ix + ix_cent + (iy + iy_cent) * (*nx_out)] = 
                  -out_image[iy + ix_cent + (-ix + iy_cent) * (*nx_out)];
       else
       out_image[-ix + ix_cent + (iy + iy_cent) * (*nx_out)] = 
                  out_image[iy + ix_cent + (-ix + iy_cent) * (*nx_out)];
       }

return(0);
}
/***************************************************************
*
* Interface to Fortran programs
* Output internal static arrays
*
****************************************************************/

/********************************************
* Fortran array NGT(NBMAX)
* Index: 1 to NBMAX (but c compatible on Feb 2nd 1994)
********************************************/
int COVER_NGT(INT4 *ngt_val, INT4 *index)
{
*ngt_val = ngt[*index];
return(0);
}
/********************************************
* Fortran array IXY(2,NULL:NBMAX)
* Second index: like c programs
********************************************/
int COVER_IXY(INT4 *ixy1_val, INT4 *ixy2_val, INT4 *index)
{
*ixy1_val = ixy_1[*index];
*ixy2_val = ixy_2[*index];
return(0);
}
/********************************************
* Fortran array IXY(NULL:NBMAX)
* Second index: like c programs
********************************************/
int COVER_IXY_1D(INT4 *ixy1_val, INT4 *index)
{
*ixy1_val = ixy_1[*index];
/* JLP98: I set it to "identity", since array overflood somewhere...: */
/*  *ixy1_val = *index; */
/* end of JLP98 */
return(0);
}
/********************************************
* Fortran array NBCOUV(-IRMAX:IRMAX,NULL:IRMAX)
* Indices: (like c programs) 
********************************************/
int COVER_NBCOUV(INT4 *nb_val, INT4 *index1, INT4 *index2, INT4 *ir_max)
{
INT4 nbc_dim, nbc_offset;

nbc_dim = 2 * (*ir_max) + 1; 
nbc_offset = (*ir_max) + 1; 

*nb_val = nbcouv[(*index1) + nbc_offset + nbc_dim * (*index2)];
return(0);
}
/*******************************************
* Fortran array KLM(3,NGMAX)
* Input indices: 1 to 3 and 1 to NGMAX
* Output index: klm_val from 0 to NBMAX 
*******************************************/
int COVER_KLM(INT4 *klm_val, INT4 *klm_index, INT4 *ng_index)
{

*klm_val = klm[(*ng_index - 1) * 3 + (*klm_index - 1)];

return(0);
}
/********************************************************************
* Subroutine COVERA_MASK_1D
* To compute the frequencies of the uv coverage, 
* and the elements of the A matrix.
* (To solve eventually the equation A*X=Y, 
*   i.e., to invert the bispectral relations)
*
* Compute u-v coverage with the mask and within a segment of length 2*ir.
*
* Input:
*  ir: maximum radius of the uv-coverage
*  max_nclosure: maximum number of closure relations
*  mask[nx]: frequency mask ( >=1. if frequency accessible, 0. otherwise) 
*
* Output:
*  NBETA: Number of elements of the spectral list (Number of columns of A)
*  NGAMMA: Number of elements of the bispectral list (Number of rows of A)
*
* OUTPUT values that are hidden as static variables:
* INT4 *nbcouv, *ixy_1, *ixy_2, *ngt, *klm;
* 
*********************************************************************/
int COVERA_MASK_1D(float *mask, INT4 *nx_mask, INT4 *ir, INT4 *max_nclosure,
                   INT4 *nbeta, INT4 *ngamma)
{
/* max_nclosure: max number of closure relations per spectral term */
/* Limitation of the maximum number of closure relations
* (to simulate Knox-Thompson and/or accelerate computations...) */ 
INT4 isize, nxm, ir_max, iw1;
register INT4 i, k;

/* JLP 94
 ir_max = IRMAX;
*/
ir_max = *ir;
if(ir_max > IRMAX)
   {printf(" COVERA_MASK_1D/Fatal error: ir max = %d \n",IRMAX);
    exit(-1);
   }

nxm = *nx_mask;
/* Check first that mask is correct: */
if(nxm < 2 * (*ir) )
  {
   printf(" COVERA_MASK_1d/Fatal error: mask is too small !!!\n");
   printf(" nxm = %d ir = %d \n",nxm,*ir);
   exit(-1);
  }
 
/* nbcouv( x from 0 to IRMAX ) */
/* JLP94 */
isize = (ir_max + 1) * sizeof(INT4);
JLP_GVMI(&nbcouv,&isize);
/* ixy_1 for x;  nb from 0 to NBMAX ) */
isize = (NBMAX + 1) * sizeof(INT4);
JLP_GVMI(&ixy_1,&isize);
JLP_GVMI(&ngt,&isize);
/*JLP 94 */
for(i = 0; i < (ir_max + 1); i++) nbcouv[i] = 0;
 
/* Computing the A matrix as generated by the uv-coverage 
   (defined only by IR for a full pupil) */

/* First call to compute bispectrum size NGMAX = ngamma: */
couv_mask_1D(mask,nxm,ir,nbeta,ngamma,*max_nclosure,1);

/* Allocation of memory for bispectrum index: */
isize = 3 * (*ngamma) * sizeof(INT4);
JLP_GVMI(&klm,&isize);

/* Second call to compute spectral and bispectral lists: */
couv_mask_1D(mask,nxm,ir,nbeta,ngamma,*max_nclosure,0);
 
printf(" covera/uv-coverage, IR = %d \n",*ir);
printf("  NBETA (spectral list) = %d \n",*nbeta);
printf("  NGAMMA (bispec. list) = %d \n",*ngamma);
 
/* Computing the number of closure relations for some values of NB: */
for ( k = 2; k < 5 ; k ++)
  {
   i = 3 + *nbeta * k / 5;
   iw1 = sqrt((double)(ixy_1[i] * ixy_1[i]));
   printf(" NB = %d irad=%d Closure relations: %d \n",
             i,iw1,ngt[i] - ngt[i - 1]);
  }
 
return(0);
}
/*******************************************************************
* couv_mask_1D defines the uv-coverage and the A matrix:
* Input:
* mask[nx]: frequency mask ( >=1. if frequency accessible, 0. otherwise) 
* IR: maximum radius of the uv-coverage
* dimension_only: flag set to one if only bispectral list dimension is required
*
* Output:
* NBETA, NGAMMA: number of elements of the spectral and bispectral lists
*
* In common blocks (output): the uv-coverage is accessible from 2 sides:
*
* NBCOUV(I,J): uv-coverage (i.e. number of the spectral list for
*              the pixel(I,J) of the spectrum (I=0,J=0 for null frequency)
*
* IXY(1,I) X-coordinate of the pixel number I in the spectral list
*              (this allows another entry for the uv-coverage)
*
*************************************************************************/
static int couv_mask_1D(float *mask, INT4 nxm, INT4 *ir, INT4 *nbeta, 
                        INT4 *ngamma, INT4 max_nclosure, INT4 dimension_only)
{
/* 
* nb: beta index, i.e. spectral list index
* ng: gamma index, i.e. bispectral list index
*/
INT4 nb, ng, icent;
register INT4 i;
 
/* Coordinate of center of the mask: */
icent = nxm/2; 
 
/* Easy cases: */
/* First spectral value at (0) */
 nb = 0;
 i = 0; 
 if(mask[i + icent] > 0.)
   {
   ixy_1[nb] = i; 
   nbcouv[i] = nb;
   }
 
/* Second spectral value at (1), nb=1 */
 i = 1; 
 if(mask[i + icent] > 0.)
   {
   nb++;
   ixy_1[nb] = i;
   nbcouv[i] = nb;
   }

/* Reseting the total number of elements the  bispectral list */
 ng = 0;
 
/* Main loop: 
*/
for ( i = 2; i <= *ir; i++) 
  {
    affect_mask_1D(mask,nxm,i,&nb,&ng,max_nclosure,dimension_only);
  }
 
/* NBETA: Total number of the spectral list (Number of columns of the X matrix)
* NGAMMA: Total number of the bispectral list (Number of rows of the X matrix)
* (Remember, we have to solve    A*X = Y) */
  *nbeta = nb;
  *ngamma = ng;

  return(0);
}
/*******************************************************************
* affect_mask_1D 
*
* Check if (isx,isy) is a valid uv vector and can be added to the spectral list.
* Then increment nb by one, and look for all new closed triangles whose sides
* are constituted with previous vectors of the spectral list and the new 
* uv vector.
*
* INPUT:
*  mask: uv mask used to invalidate some uv vectors.
*  nxm: X size of mask array
*  isx : coordinate of the uv vector to be validated
*  nb: previous value of the number of items in the spectral list
*  ng: previous value of the number of items in the bispectral list
*  max_nclosure: maximum number of closure relations allowed by the user
*  dimension_only: flag set to one if this routine is only called
*                  to estimate the value of NGAMMA and NBETA
*
* OUTPUT:
*  nb: new value of the number of items in the spectral list
*      (i.e. index of the new uv vector in the spectral list)
*  ng: new value of the number of items in the bispectral list
*      (incremented in this routine by the number of new triangles found
*      that involves the new uv vector) 
********************************************************************/
static int affect_mask_1D(float *mask, INT4 nxm, INT4 isx, INT4 *nb, INT4 *ng,
                          INT4 max_nclosure, INT4 dimension_only)
{ 
INT4 nbs, nbk, nbr, itx, icent;
register INT4 nbq, iklm;

/* Coordinate of center of the mask: */
icent = nxm/2;
 
/* First condition: input point has to be accessible */
  if(mask[isx + icent] <= 0.) return(-1);

/* If so, record the new point of the uv coverage (spectral list): */
 (*nb)++;
 ixy_1[*nb] = isx;
 nbcouv[isx] = *nb;

/* Searching for the couples associated with the point NBS=NB
* Generating the bispectral list and building the rows of the A matrix:
*/
 nbs = *nb;
 
/* Loop on all the possible points (Q): */
 iklm = 3 * (*ng);
 for( nbq = 1; nbq < nbs; nbq++)
  { 
/* JLP 94: add an exit test when maximum number of closure relations
has been found (to simulate Knox-Thompson and/or accelerate computations...) */ 
   if((*ng - ngt[(*nb) -1]) == max_nclosure) break;

/* Coordinates of the vector T = S - Q */
   itx = isx - ixy_1[nbq];
 
/* Work within the segment of length IRS, so we can a priori reject
* the points outside the window [-IRS,+IRS]:
*/
   if(itx >= -isx && itx <= isx)
      {
/* Case number 1 only for 1D case (which could be : k=t, l=q, k+l=s) */

/* t should be accessible too: */
            if(mask[itx + icent] > 0.)
            {
              nbk = nbcouv[itx];
 
/* We select this couple (U,V) if the vector NBK is in [0,NBQ] */
              if( nbk != 0 && nbk <= nbq)
                 {
/* Add new k,l,m coordinates if more than dimension is wanted */
                 if(!dimension_only)
                   {
                   klm[iklm] = nbk; 
                   iklm++;
                   klm[iklm] = nbq; 
                   iklm++;
                   klm[iklm] = nbs;
                   iklm++;
                   }
                 (*ng)++;
                 }
          }
      else
          { 
/* Case number 2 (which could be : r=-t, s=s, r+s=m=q) */
/* r should be accessible too: */
            if(mask[-itx + icent] > 0.)
            {
             nbr = nbcouv[-itx];
/* We select this couple (R,S) if the vector NBR is in [0,NBS] */
              if( nbr != 0 && nbr <= nbs)
                  {
/* Add new k,l,m coordinates if more than dimension is wanted */
                  if(!dimension_only)
                     {
                     klm[iklm] = nbr; 
                     iklm++;
                     klm[iklm] = nbs; 
                     iklm++;
                     klm[iklm] = nbq;
                     iklm++;
                     }
                  (*ng)++;
                  }
            }
/* Nota: we can only have L=NBS or M=NBS (never K=NBS) */
          }
      } 

/* End of loop on nbq */
}
 
/* NGT(NB) is the number of the last U,V couple of the group NB=NBS=S: */
   ngt[nbs] = *ng;

 return(0);
}
/***************************************************************
* C interface with COVER_KLM
* It is designed for C calls with C-like arrays: 
* loops: for(ng=0; ng < ngamma; ng++)
* (COVER_KLM was designed for Fortran arrays)
*
****************************************************************/
int cover_klm0(INT4 *k, INT4 kk, INT4 ng)
{
INT4 kkk, nng;
kkk = kk;
nng = ng+1;
COVER_KLM(k,&kkk,&nng);
return(0);
}
/***************************************************************
* C interface with COVER_KLM
* It is designed for C calls with Fortran-like arrays: 
* loops: for(ng=1; ng <= ngamma; ng++)
* (COVER_KLM was designed for Fortran arrays)
*
****************************************************************/
int cover_klm1(INT4 *k, INT4 kk, INT4 ng)
{
INT4 kkk, nng;
kkk = kk;
nng = ng;
COVER_KLM(k,&kkk,&nng);
return(0);
}
/***************************************************************
* C interface with COVER_NGT
* It is designed for C calls with Fortran-like arrays: 
* (COVER_NGT was designed for Fortran arrays)
* loops: for(ng=1; ng <= ngamma; ng++)
*
****************************************************************/
int cover_ngt1(INT4 *ng2, INT4 nb)
{
INT4 nbb;
nbb = nb;
COVER_NGT(ng2,&nbb);
return(0);
}
/***************************************************************
* Routine to output spectral and bispectral lists (for debugging purpose) 
*
****************************************************************/
int output_lists_coverage(INT4 *nbeta, INT4 *ngamma)
{
register int nb, ng;
int k, l, m;
INT4 ng2;
char filename[40];
FILE *fp1;

strcpy(filename,"spectral_lists.dat");
if((fp1 = fopen(filename,"w")) == NULL)
{printf("Error opening list file %s \n",filename);
 return(-1);
}

fprintf(fp1,"#################################### \n");
fprintf(fp1,"# nbeta = %d  ngamma = %d \n", *nbeta, *ngamma);
fprintf(fp1,"#################################### \n");
fprintf(fp1,"# Spectral list: index, ix, iy \n");
fprintf(fp1,"#################################### \n");
 for(nb = 0; nb < *nbeta; nb++)
  {
    fprintf(fp1,"%d %d %d \n", nb, ixy_1[nb], ixy_2[nb]);
  }
fprintf(fp1,"#################################### \n");
fprintf(fp1,"# NGT: index ngt,... (to compute nb I use ng=ngt[nb-1]; ng<ngt[nb]; ng++)\n");
fprintf(fp1,"#################################### \n");
 for(nb = 0; nb <= *nbeta; nb++)
  {
    cover_ngt1(&ng2,(INT4)nb);
    fprintf(fp1,"%d %d \n", nb, (int)ng2);
  }
fprintf(fp1,"#################################### \n");
fprintf(fp1,"# Bispectral list: index, k, l, m \n");
fprintf(fp1,"#################################### \n");
 for(ng = 0; ng < *ngamma; ng++)
  {
    cover_klm0(&k,1,ng);
    cover_klm0(&l,2,ng);
    cover_klm0(&m,3,ng);
    fprintf(fp1,"%d %d %d %d \n", ng, k, l, m);
  }
fclose(fp1);

return(0);
}
/***************************************************************
* Routine to output spectral and bispectral lists (for debugging purpose) 
*
****************************************************************/
int output_lists_coverage_1D(INT4 *nbeta, INT4 *ngamma)
{
register int nb, ng;
INT4 k, l, m, ng2;
char filename[40];
FILE *fp1;

strcpy(filename,"spectral_lists_1D.dat");
if((fp1 = fopen(filename,"w")) == NULL)
{printf("Error opening list file %s \n",filename);
 return(-1);
}

fprintf(fp1,"#################################### \n");
fprintf(fp1,"# nbeta = %d  ngamma = %d \n", *nbeta, *ngamma);
fprintf(fp1,"#################################### \n");
fprintf(fp1,"# Spectral list: index, ix, iy \n");
fprintf(fp1,"#################################### \n");
 for(nb = 0; nb <= *nbeta; nb++)
  {
    fprintf(fp1,"%d %d \n", nb, ixy_1[nb]);
  }
fprintf(fp1,"#################################### \n");
fprintf(fp1,"# NGT: index ngt,... (to compute nb I use ng=ngt[nb-1]; ng<ngt[nb]; ng++)\n");
fprintf(fp1,"#################################### \n");
 for(nb = 0; nb <= *nbeta; nb++)
  {
    cover_ngt1(&ng2,(INT4)nb);
    fprintf(fp1,"%d %d \n", nb, (int)ng2);
  }
fprintf(fp1,"#################################### \n");
fprintf(fp1,"# Bispectral list: index, k, l, m \n");
fprintf(fp1,"#################################### \n");
 for(ng = 0; ng < *ngamma; ng++)
  {
    cover_klm0(&k,1,ng);
    cover_klm0(&l,2,ng);
    cover_klm0(&m,3,ng);
    fprintf(fp1,"%d %d %d %d \n", ng, k, l, m);
  }
fclose(fp1);

return(0);
}
