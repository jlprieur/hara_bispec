/******************************************************************
* photon_model
*
* Fit a model of the background 
*
* JLP
* Version 03/09/2008
*******************************************************************/
#include "jlp_ftoc.h"

#define DEBUG

/* Contained here and defined in jlp_cover_mask.h: 
int fit_photon_model_to_modsq(float *modsq, float *phot_modsq, float *coeff0, 
                              float *coeff1, int nx, int ny, int idim);
int compute_nphotons(float *modsq, float *phot_modsq, float *a1,
                     int nx, int ny, int idim, float rmin, float rmax);
int compute_nphotons_gauss(float *modsq, float *phot_modsq, float *a1, 
                           int nx, int ny, int idim, float rmin, float rmax);
*/

/**************************************************************************
* Compute the coefficients to subtract the background model on the frame
* and subtract this model.
*
* Linear regression:
*
* Minimum of Sum ( z - a1 f - a0)^2
* is reached when gradient (d/da1, d/da2)is nul, i.e., when:
* sum_fz = a1 sum_ff + a0 sum_f
* sum_z  = a1 sum_f + a0 sum_1
*
* (Indeed: Sum of z2 + a1^2 f^2 + a0^2 - 2 a1 z f - 2 a0 z + 2 a1 a0 f
*  leads to d/da1 = 2 a1 f^2 - 2 z f + 2 a0 f = 0
*      and  d/da0 = 2 a0 - 2 z + 2 a1 f       = 0 )
*
* INPUT:
* modsq: power spectrum (with zero frequency at nx/2 ny/2) 
* phot_modsq: power spectrum of photon response
*
**************************************************************************/
int fit_photon_model_to_modsq(float *modsq, float *phot_modsq, float *coeff0, 
                              float *coeff1, int nx, int ny, int idim)
{
double wback, wdata; 
double sum_1, sum_f, sum_z, sum_fz, sum_ff, a0, a1, det;
float rad2, rad2_min, rmin;
int jj, ixc, iyc; 
register int i, j;

*coeff0 = 0.;
*coeff1 = 0.;

ixc = nx/2; iyc = ny/2;

/* Central pixels are "polluted" by the power spectrum of the object, 
* so I neutralize them for the fit: 
*/
rmin = 0.3;
rad2_min = SQUARE(rmin * (float)nx);

sum_1 = 0.;
sum_f = 0.;
sum_z = 0.;
sum_fz = 0.;
sum_ff = 0.;
for(j = 0; j < ny; j++) {
  jj = j * idim;
  for(i = 0; i < nx; i++) {
   wback = phot_modsq[i + jj];
   wdata = modsq[i + jj];
   rad2 = SQUARE(j-iyc) + SQUARE(i - ixc);
   if(rad2 > rad2_min){
     sum_1 += 1.;
     sum_f += wback;
     sum_z += wdata;
     sum_fz += wback * wdata;
     sum_ff += SQUARE(wback);
     }
  }
 }

/* Resolution with the determinant of the system:
* sum_fz = a1 sum_ff + a0 sum_f
* sum_z  = a1 sum_f + a0 sum_1
*/
det = sum_ff * sum_1 - sum_f * sum_f;
a1 = (sum_1 * sum_fz  - sum_f * sum_z) / det;
a0 = (sum_ff * sum_z  - sum_f * sum_fz) / det;

#ifdef DEBUG
printf("sum_ff=%f sum_1=%f sum_f=%f sum_fz=%f\n", sum_ff, sum_1, sum_f, sum_fz);
printf("subtract_photon_model/ a1=%f a0=%f (det=%f)\n", a1, a0, det);
#endif

*coeff0 = a0;
*coeff1 = a1;

return(0);
}
/**************************************************************************
* Contribute to determining the number of photons per frame,
* by fitting a constant a1 to the ratio modsq/phot_modsq 
* for the high frequencies
*
* INPUT:
* modsq: power spectrum (with zero frequency at nx/2 ny/2) 
* phot_modsq: power spectrum of photon response
* rmin, rmax: relative radii limits for fitting photon reponse
*             (as a fraction of nx)
*
**************************************************************************/
int compute_nphotons(float *modsq, float *phot_modsq, float *a1, 
                     int nx, int ny, int idim, float rmin, float rmax)
{
float work, epsilon = 1.e-10, sigma, *tmp, mean, median;
double sum, sumsq; 
float rad2_min, rad2_max, rad2;
int jj, ixc, iyc, npts; 
register int i, j;

/* Default value: */
*a1 = 1.;

ixc = nx/2; iyc = ny/2;

tmp = (float*)malloc(nx * ny * sizeof(float));

/* Central pixels are "polluted" by the power spectrum of the object, 
* so I neutralize them for the fit: 
*
* For ads11454 with 200 ph/frame (Paper VI)
rmin = 0.4;
rmax = 100;
rad2_min = SQUARE(0.4 * (float)nx);
rad2_max = 100;
* For ads11454 with 20 ph/frame
rmin = 0.1;
rmax = 0.3;
rad2_min = SQUARE(0.1 * (float)nx);
rad2_max = SQUARE(0.3 * (float)nx);
*/
rad2_min = SQUARE(rmin * (float)nx);
rad2_max = SQUARE(rmax * (float)nx);

sum = 0.;
sumsq = 0.;
npts = 0;
for(j = 0; j < ny; j++) {
  jj = j * idim;
  for(i = 0; i < nx; i++) {
   rad2 = SQUARE((j-iyc)/1.2) + SQUARE(i - ixc);
   if(rad2 > rad2_min && rad2 < rad2_max){
     work = modsq[i + jj] / (phot_modsq[i + jj] + epsilon);
     sum += work;
     sumsq += SQUARE(work);
     tmp[npts] = work;
     npts++;
     }
  }
 }

if(npts <= 5) {
  fprintf(stderr, "rad_min,max=%.3f,%.3f\n", rad2_min, rad2_max);
  fprintf(stderr, "compute_nphotons/Fatal error: two few points: npts=%d\n", npts);
  free(tmp);
  exit(-1);
  }

 mean = sum/(float)npts;
 sigma = sqrt(sumsq / (float)npts - SQUARE(mean));
 JLP_MEDIAN(tmp, npts, &median);

#ifdef DEBUG
printf("compute_nphotons/ mean=%f sigma=%f median=%f npts=%d\n", 
        mean, sigma, median, npts);
#endif

*a1 = median;

free(tmp);

/* JLP SEP08: */
/* SEP08: add a correction to the model (DEBUGGGG) */
#if 0
for(j = 0; j < ny; j++) 
  for(i = 0; i < nx; i++) 
    modsq[i + j * idim] = (modsq[i + j * idim] 
                          / ((*a1) * phot_modsq[i + j * idim] + epsilon)) - 1.;
for(j = 0; j < ny; j++) {
  jj = j * idim;
  for(i = 0; i < nx; i++) {
  rad2 = SQUARE((j - ny/2)/1.2) + SQUARE(i - nx/2);
  if(rad2 >= rad2_min)
        phot_modsq[i + jj] = (phot_modsq[i + jj] + modsq[i + jj] / (*a1))/2.;
  }
 }
#endif
return(0);
}
/*************************************************************************
*
**************************************************************************/
int subtract_photon_model(float *modsq, float *phot_modsq, float a0, 
                          float a1, int nx, int ny, int idim)
{
double sum_1, sum_f, sum_ff;
float variance, rad2, rad2_min, rmin;
int jj, ii, ixc, iyc;
register int i, j;

ixc = nx/2; iyc = ny/2;

/* Subtract scaled model: */
for(j = 0; j < ny; j++) {
  jj = j * idim;
  for(i = 0; i < nx; i++) {
   ii = i + jj;
     modsq[ii] -= (a1 * phot_modsq[ii] + a0);
  }
}

/* Compute standard deviation of the outer parts: */
sum_1 = 0;
sum_f = 0.;
sum_ff = 0.;
rmin = 0.3;
rad2_min = SQUARE(rmin * (float)nx);
for(j = 0; j < ny; j++) {
  jj = j * idim;
  for(i = 0; i < nx; i++) {
   ii = i + jj;
    rad2 = SQUARE(j-iyc) + SQUARE(i - ixc);
    if(rad2 > rad2_min){
     sum_1 += 1.;
     sum_f += modsq[ii];
     sum_ff += SQUARE(modsq[ii]);
    }
  }
}
if(sum_1 > 4.) {
  sum_f /= sum_1;
  sum_ff /= sum_1;
  variance = sum_ff - sum_f * sum_f;
  printf("photon_mode/variance=%f\n", variance);
  } else {
  fprintf(stderr,"Fatal error: sum_1=%d (reduce rad2_min in the program...)\n",
          (int)sum_1);
  exit(-1);
  }

return(0);
}
/******************************************************************
* Fit a 2D Gaussian to an image
*
* OUTPUT:
* phot_modsq: model of photon response
*****************************************************************/
int compute_nphotons_gauss(float *modsq, float *phot_modsq, float *a1, 
                           int nx, int ny, int idim, float rmin, float rmax)
{
INT4 npts, ifail;
float errors[4], sigx, sigy, xc, yc, rho;
float rad2, rad2_min, rad2_max;
float *xx, *yy, *f1;
register int i, j;

xx = (float *)malloc(nx * ny * sizeof(float));
yy = (float *)malloc(nx * ny * sizeof(float));
f1 = (float *)malloc(nx * ny * sizeof(float));

/* Transfer of data: */
rad2_min = SQUARE(rmin * (float)nx);
rad2_max = SQUARE(rmax * (float)nx);
npts = 0;
 for(j = 0; j < ny; j++)
   for(i = 0; i < nx; i++) {
    rad2 = SQUARE((j - ny/2)/1.2) + SQUARE(i - nx/2);
    if(rad2 > rad2_min && rad2 < rad2_max) {
     xx[npts] = (float)i;
     yy[npts] = (float)j;
     f1[npts] = modsq[i + j * idim];
     npts++;
     }
   }
/* fit_gauss:
* f1: array of values f1(xx[i],yy[i])
* xx, yy: arrays of coordinates x, y
* npts: number of points
*
* Parameters: 
* sigx, sigy, xc, yc, rho
* errors[4] : errors[sigx,sigy,xc,yc,rho]
*/
printf("rmin=%.3f rmax=%.3f npts=%d\n", rmin, rmax, npts);

/*
int jlp_fit_gauss_flt(float *xx0, float *yy0, float *ff0, INT4 *npts,
                      float *sigx0, float *sigy0, float *xc0, float *yc0,
                      float *rho0, float *errors0, INT4 *ifail)
*/
 jlp_fit_gauss_flt(xx, yy, f1, &npts, &sigx, &sigy, &xc, &yc, &rho, 
                   errors, &ifail);
/*
* ifail = 0 if correct
*         -3: all values are negative or null!
*/
if(ifail) {
  fprintf(stderr,"photon_gauss/Error in jlp_fit_gauss: ifail=%d\n", ifail);
  free(xx); free(yy); free(f1);
  return(-1);
  }
printf("fit_gauss: xc=%.3f yc=%.3f sigx=%.3f sigy=%.3f rho=%e\n",
        xc, yc, sigx, sigy, rho);

/* Build a photon model: */
 for(j = 0; j < ny; j++)
   for(i = 0; i < nx; i++) {
   rad2 = SQUARE(((float)i - xc) / sigx) + SQUARE (((float)j - yc) / sigy); 
   phot_modsq[i + j * idim] = rho * exp(-rad2); 
   }

*a1 = 1.;

free(xx); free(yy); free(f1);
return(0);
}
