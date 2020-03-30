/*******************************************************************
*
* Contains:
* bispec3
* correct_photon_noise
*
* Called by jlp_cover2.c and jlp_cover_mask.c
*
* JLP
* Version 24-08-93
*******************************************************************/
#include <stdio.h>
#include <math.h>
#include <jlp_ftoc.h>
/*********************************************************************
* BISPEC1
* Same as BISPEC2, but older version, and sum of full bispectrum terms.
* (amplitude and phase)
* Reads a spectrum from an real image (from observations)
*
* Computes the bispectrum and spectrum lists from RE and IM
* already computed before calling this routine.
*
* Integrates the squared modulus of the spectrum and
* the phase factor of the bispectrum. Corrects from photon noise effects.
*
* Photon noise correction (cf JOSA 2, 14, Wirnitzer):
* (When spectrum normalized to one in the center)
* <i(u)>**2 = E(D...)/N**2 - 1/N
*
* <i(3)(u,v)> = E(D(3)(u,v)/N**3 - E(D(2)(u)/N**2)/N - E(D(2)(v)... +2/N**3)
*
* Input:
* RE(NX,NY)
* IM(NX,NY)
*
* Output:
* MODSQ: Sum of the modulus squared
* YCE1: phase factor of the bispectrum (real and imaginary)
* YCE1(.,3): sum of square real parts of the bispectrum
* YCE1(.,4): sum of square imaginary parts of the bispectrum
********************************************************************/
int BISPEC3(image,im,modsq,snrm,nx,ny,idim,yce1,ir,nbeta,ngamma)

float image[], im[], modsq[], snrm[], yce1[];
long *nx, *ny, *idim, *ir, *nbeta, *ngamma;
{
}
