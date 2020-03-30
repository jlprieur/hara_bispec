/**************************************************************** 
*  correct_bisp_auto1
*  (derived from hege's method of fitting a model to the photon response...)
*
*  Contains: photon_corr_auto1()
*  called by correct_bisp.c
*
*  Fits the photon response to the power spectrum
*
* JLP
* Version 10-09-2008
*****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <jlp_ftoc.h>
#include "jlp_cover_mask.h"

/*
#define DEBUG
*/

/* Contained here and defined in jlp_cover_mask.h: 
int photon_corr_auto1(float *bispp, float *modsq, float *phot_modsq,
                          INT4 *nx, INT4 *ny, INT4 *bisp_dim, float *xphot,
                          INT4 *nbeta, INT4 *ngamma)
*/
static int compute_nphotons_test(float *modsq, float *bispp, int ngamma, 
                                 int nx, int ny, float a1, float *xphot);
static int compute_k_factor(float modsq_0, float phot_modsq_0, 
                            float a1, float k_range, float *k_factor);

/**************************************************************** 
* Photon noise correction found by JLP in May2008 
*
*
* Photon noise correction (cf JOSA 2, 14, Wirnitzer):
* (When spectrum normalized to one in the center)
* <i(u)>**2 = E(D...)/N**2 - 1/N
*
* <i(3)(u,v)> = E(D(3)(u,v)/N**3 - E(D(2)(u)/N**2)/N - E(D(2)(v)... +2/N**3)
*
* INPUT:
* modsq[]: mean normalized squared modulus (with zero frequency at nx/2,ny/2)
* phot_modsq[]: power spectrum of photon response (with 0 freq. at nx/2,ny/2) 
* bispp[]: bispectrum (sorted in a "standard" way, and thus different
*         from what is needed by photon_corr...)
*
* OUTPUT:
* xphot: mean photon flux per frame 
*****************************************************************/
int photon_corr_auto1(float *bispp, float *modsq, float *phot_modsq,
                      INT4 *nx, INT4 *ny, INT4 *bisp_dim, float *xphot,
                      INT4 *nbeta, INT4 *ngamma)
{
float *xmodsq, *phot_re, a1, epsilon=1.e-12; 
float k_range, k_factor, k_correction, rmin, rmax;
double work1;
int status;
INT4 k1, k2, k3, iix, iiy, nb, ixc, iyc;
INT4 ng, iklm;
register int i;

ixc = (*nx)/2;
iyc = (*ny)/2;

/* Old and DEBUG version: I used the value entered by the user...
k_factor = *modsq;
*/

/* JLP96: */
/* Number of photons:  xphot**2 + xphot - modsq_0 = 0
* (for photon-counting devices!)
* Delta = b**2 - 4 a*c = 1 + 4 * modsq_0
* Solutions: (- b +/- sqrt(Delta))/2a
*  i.e.,     (-1 + sqrt(1 + 4 * modsq_0) )/2
*/
  work1 = (-1. + sqrt((double)(1. 
       + 4. * modsq[ixc + (iyc * (*nx))]/phot_modsq[ixc + (iyc * (*nx))])))/2.;
  printf(" My rough estimate of xphot (from modsq[0 freq]) is: %f\n", work1);
/* End of JLP96. */

/* Fit phot_modsq = |g(v)|^2 to the outer parts of modsq: 
fit_photon_model_to_modsq(modsq, phot_modsq, &a0, &a1, *nx, *ny, *nx);
*/

/* New version with compute_nphotons (in "photon_model.c"): 
* Computes a1 = xphot / k_factor^2, were xphot is the number of photons per frame,
* by fitting a constant to the ratio modsq/phot_modsq 
* for the high frequencies
*
*/

/* JLP SEP08/DEBUG ONLY */
/* Attempt of fitting a Gaussian, but still not working (Sept 11th 2008) */
if(0) {
rmin=0.3;
rmax=100.;
printf(" compute_nphoton_gauss/ Fit a Gaussian\n");
compute_nphotons_gauss(modsq, phot_modsq, &a1, *nx, *ny, *nx, rmin, rmax);
k_range = 40.;
} else {
rmin = 0.1;
rmax = 0.3;
printf(" compute_nphoton/Fit the model contained in phot_gv7_128.fits\n");
compute_nphotons(modsq, phot_modsq, &a1, *nx, *ny, *nx, rmin, rmax);
k_range = 10.;
}

/* Compute k_factor from modsq_0 and a1: */
status = compute_k_factor(modsq[ixc + iyc * (*nx)], 
                          phot_modsq[ixc + iyc * (*nx)], a1, k_range, &k_factor);
if(status) {
   fprintf(stderr, "Fatal error in compute_k_factor: could not find a solution...\n");
   exit(-1);
  }

/* N = a1 * k^2 */
*xphot = a1 * SQUARE(k_factor);
printf("Number of photons/frame: %.3f\n", *xphot);

/* Correction of the experimental power spectrum: */
for(i = 0; i < (*nx) * (*ny); i++){ 
    modsq[i] = (modsq[i] / (a1 * phot_modsq[i] + epsilon)) - 1.;
/* Multiplication with N of the corrected power spectrum, to make it 
compatible with the bispectrum (for bisp_to_image...): */
    modsq[i] *= *xphot;
    }

/* Allocation of memory for xmodsq, modsq on the spectral list: */
if((xmodsq = (float *)malloc((*nbeta + 1) * sizeof(float))) == NULL)
 {
 printf("photon_corr_mask/fatal error, allocating memory space \n");
 exit(-1);
 }

if((phot_re = (float*)malloc((*nbeta + 1) * sizeof(float))) == NULL)
 {
 printf("photon_corr_mask/fatal error, allocating memory space \n");
 free(xmodsq);
 exit(-1);
 }

/* IXY(1,NB) and IXY(2,NB) are the X,Y coordinates of the
* element NB of the spectral list.
* As they (i.e. IXY) can be negative, and that zero frequency is at (0,0),
* we do the following transformation:
*/
  for(nb = 0; nb <= *nbeta; nb++)
     {
/* Get coordinates (i,j) from nb index: */
       COVER_IXY(&iix, &iiy, &nb);
/* OLD VERSION (assuming that zero freq. is at (0,0):
       iix = ((iix + *nx) % *nx);
       iiy = ((iiy + *ny) % *ny);
*/
/* NEW VERSION (since that zero freq. has been shifted to (nx/2ny/2):
*/
       iix += (*nx)/2;
       iiy += (*ny)/2;
/* The photon response is assumed to be symmetric and real, hence its FT is real: */
       phot_re[nb] = sqrt(phot_modsq[iix + iiy * (*nx)]);
/* xmodsq is used to store the corrected square modulus */
       xmodsq[nb] = modsq[iix + iiy * (*nx)];
#ifdef DEBUG
       if(nb < 4) {
         printf(" nb = %d iix, iiy: %d %d \n",nb,iix,iiy);
         printf(" xmodsq=%f phot_modsq=%f\n", xmodsq[nb], phot_modsq[nb]);
       }
#endif
     }


/****************************************************************/
/* Phase factor of the bispectrum (with bispectral list):
* and correction from photon noise effects:
*/

/* Division by the real and imaginary parts of the bispectrum 
* by the photon response: 
* (from 2000 onwards...)
*/
for(ng = 0; ng < *ngamma; ng++) {
     cover_klm0(&k1,1,ng);
     cover_klm0(&k2,2,ng);
     cover_klm0(&k3,3,ng);
     work1 = phot_re[k1] * phot_re[k2] * phot_re[k3] + epsilon;
     bispp[ng] /= work1;
     bispp[ng + *bisp_dim] /= work1;
  }

/* Compute xphot, the mean number of photons per frame: */
/*
 compute_nphotons_test(modsq, bispp, *ngamma, *nx, *ny, a1, xphot);
*/

 iklm = 0;
k_correction = 1. / (k_factor * k_factor * k_factor);
for(ng = 0; ng < *ngamma; ng++)
 {
     cover_klm0(&k1,1,ng);
     cover_klm0(&k2,2,ng);
     cover_klm0(&k3,3,ng);

/* Photon noise correction 
     work1 = (xmodsq[k1] + xmodsq[k2] + xmodsq[k3] - 2 * *xphot);
     if xmodsq is not corrected
     or
     work1 = (xmodsq[k1] + xmodsq[k2] + xmodsq[k3] + *xphot);
     if xmodsd is corrected (i.e. xmodsq := xmodsq - xphot)
*/
     work1 = (xmodsq[k1] + xmodsq[k2] + xmodsq[k3]  + *xphot) * k_correction; 
#ifdef DEBUG
     if(ng < 4) {
        printf(" ng = %d k, l, m: %d %d %d \n",ng,k1,k2,k3);
        printf(" xmodsq1, xmodsq2, xmodsq3 %f %f %f (k_corr=%f)\n",
                xmodsq[k1], xmodsq[k2], xmodsq[k3], k_correction);
        printf(" ng= %d bispp=%f correction= %f (xphot=%f)\n",
               ng, bispp[ng], work1, *xphot);
     }
#endif
     bispp[ng] -= work1;
  }

free(phot_re);
free(xmodsq);
return(0);
}
/********************************************************************** 
* Compute the mean number of photons per frame
* by comparing the corrected bispectrum and the corrected power spectrum
*
* NO LONGER USED!
***********************************************************************/
static int compute_nphotons_test(float *modsq, float *bispp, int ngamma, 
                            int nx, int ny, float a1, float *xphot)
{
double sum, sum_bisp, epsilon=1.e-8, norm1, norm2;
int npts, npts_bisp, ixc, iyc, ix1, iy1, ix2, iy2, k1, k2, ng;

ixc = nx/2;
iyc = ny/2;

sum = 0.;
sum_bisp = 0.;
npts = 0;
npts_bisp = 0;
for(ng = 0; ng < ngamma; ng++) {
     cover_klm0(&k1,1,ng);
     cover_klm0(&k2,2,ng);
/* Get coordinates (i,j) from nb index: */
     COVER_IXY(&ix1, &iy1, &k1);
     COVER_IXY(&ix2, &iy2, &k2);
     norm1 = SQUARE(ix1) + SQUARE(iy1);
     norm2 = SQUARE(ix2) + SQUARE(iy2);
/* OLD version:
     ix1 = ((ix1 + nx) % nx);
     iy1 = ((iy1 + ny) % ny);
     ix2 = ((ix2 + nx) % nx);
     iy2 = ((iy2 + ny) % ny);
*/
/* NEW VERSION (since that zero freq. has been shifted to (nx/2ny/2):
*/
     ix1 += nx/2;
     iy1 += ny/2;
     ix2 += nx/2;
     iy2 += ny/2;
     if(norm1 > 40 && norm2 > 40) {
/* If norm1 is small, I compare with modsq2 which corresponds to norm1=0 */
     if(norm1/(norm2 + epsilon) < 0.1) {
       sum += bispp[ng] - modsq[ix2 + iy2 * nx];
       npts++;
       if(npts < 20) printf(" ix1,iy1= %d %d ix2,iy2= %d %d norm1=%.1f norm2=%.1f bisp,modsq= %.1f %.1f \n", 
            ix1, iy1, ix2, iy2, norm1, norm2, bispp[ng], modsq[ix2 + iy2 * nx]);
/* If norm2 is small, I compare with modsq1 which corresponds to norm1=0 */
     } else if(norm2/(norm1 + epsilon) < 0.1) {
       sum += bispp[ng] - modsq[ix1 + iy1 * nx];
       npts++;
       if(npts < 20) printf(" ix1,iy1= %d %d ix2,iy2= %d %d norm1=%.1f norm2=%.1f bisp,modsq= %.1f %.1f \n", 
            ix1, iy1, ix2, iy2, norm1, norm2, bispp[ng], modsq[ix2 + iy2 * nx]);
    }
     if(norm1/(norm2 + epsilon) > 0.5 && norm2/(norm1 + epsilon) > 0.5) {
       sum_bisp += bispp[ng];
       npts_bisp++;
     }
    }
 }
printf("\n\n *************** npts=%d  npts_bisp=%d\n", npts, npts_bisp);
if(npts) {
   sum /= (float)npts;
   printf(" *************** mean=%f\n", sum);
} else {
sum = 1.;
}
if(npts_bisp) {
   sum_bisp /= (float)npts_bisp;
   printf(" *************** mean bisp =%f\n", sum_bisp);
}

/* TEST JLP20088 */
for(ng = 0; ng < ngamma; ng++) {
  bispp[ng] /= sum;
 }

/* DEBUG ONLY! */
 *xphot = a1;
return(0);
}
/**********************************************************************
*
* Solving equation a1^3 k^4 - k bisp_0 + 3 mod2_0 - 2 a1 = 0
*
* INPUT:
*  k_range: maximum value of the interval [0, k_range] where 
*           to look for the solution
**********************************************************************/
static int compute_k_factor(float modsq_0, float phot_modsq_0, 
                            float a1, float k_range, float *k_factor)
{
double mod2_0, bisp_0, kf, cte, work, old_value;
int status = -1, nsteps;
register int i;

*k_factor = 1.;

mod2_0 = modsq_0 / phot_modsq_0;
bisp_0 = pow(mod2_0, 1.5);

printf("modsq_0=%f phot_modsq_0=%f\n", modsq_0, phot_modsq_0);
printf(" mod2_0=%f bisp_0=%f a1=%f \n", mod2_0, bisp_0, a1);
printf(" bisp_0^{1/3}=%f mod_0^{1/2}=%f \n", pow(bisp_0, 0.333), sqrt(mod2_0));
cte = 3 * mod2_0 - 2. * a1;
work = bisp_0 - cte;
printf("cte = %f    bisp - cte = %f\n", cte, work);
work /= (a1 * a1 * a1);
printf(" bisp_prime = %f\n", work);

/* From 0 to 10, function is positive very close to zero
* then it is negative and positive again (this is the value we are looking for) */
old_value = +1;
nsteps = 1000.;
  for(i = 0; i < nsteps; i++){ 
    kf = i * k_range /(float)nsteps;
    work = a1 * a1 * a1 * kf * kf * kf *kf - kf * bisp_0 + cte; 
#ifdef DEBUG
    if(i < nsteps) printf(" %d work=%f k_f=%f \n", i, work, kf);
#endif
    if(old_value < 0. && work > 0.) {
     printf("Good value of k_f is k_factor = %f\n", kf);
     *k_factor = kf;
     status = 0;
     break;
     }
    old_value = work;
   }

if(status) fprintf(stderr,"compute_k_factor/Error, cannot find solution...\n");

return(status);
}
