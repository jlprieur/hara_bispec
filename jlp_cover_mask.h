/********************************************************************
* jlp_cover_mask.h
*
* Set of routines for uv-coverage, spectral and bispectral lists
* contained in "jlp_cover_mask.c" :
* (and in correct_bisp_auto1.c")
*   
*
* JLP
* Version 16/04/2008 
*******************************************************************/
#ifndef _jlp_cover_mask   /* Sentry */
#define _jlp_cover_mask

#define BISPEC3_SINGLE bispec3_single__
#define BISPEC3 bispec3_

/* Contained in correct_bisp_auto1.c */
int photon_corr_auto1(float *bispp, float *modsq, float *phot_modsq,
                      INT4 *nx, INT4 *ny, INT4 *bisp_dim, float *xphot,
                      INT4 *nbeta, INT4 *ngamma);

/* Contained in photon_model.c */
int fit_photon_model_to_modsq(float *modsq, float *phot_modsq, float *coeff0, 
                              float *coeff1, int nx, int ny, int idim);
int compute_nphotons(float *modsq, float *phot_modsq, float *a1,
                     int nx, int ny, int idim, float rmin, float rmax);
int compute_nphotons_gauss(float *modsq, float *phot_modsq, float *a1, 
                           int nx, int ny, int idim, float rmin, float rmax);

/* Contained in jlp_cover_mask.c */
int photon_corr_mask(double *bispp, double *modsq, double *phot_modsq,
                     INT4 *nx, INT4 *ny, INT4 *bisp_dim, float *xphot,
                     INT4 *nbeta, INT4 *ngamma, float *fraction);
int COVERA_MASK(float *mask, INT4 *nx_mask, INT4 *ny_mask, INT4 *ir,
                INT4 *max_nclosure, INT4 *nbeta, INT4 *ngamma);
int COVERA_MASK_1D(float *mask, INT4 *nx_mask, INT4 *ir, INT4 *max_nclosure, 
                   INT4 *nbeta, INT4 *ngamma);
int COVERA(INT4 *ir, INT4 *max_nclosure, INT4 *nbeta, INT4 *ngamma);
int COVER_NGT(INT4 *ngt_val, INT4 *index);
int COVER_IXY(INT4 *ixy1_val, INT4 *ixy2_val, INT4 *index);
int COVER_IXY_1D(INT4 *ixy1_val, INT4 *index);
int COVER_NBCOUV(INT4 *nb_val, INT4 *index1, INT4 *index2, INT4 *ir_max);
int COVER_KLM(INT4 *klm_val, INT4 *klm_index, INT4 *ng_index);
int cover_ngt1(INT4 *ng2, INT4 nb);
int cover_klm1(INT4 *k, INT4 kk, INT4 ng);
int cover_klm0(INT4 *k, INT4 kk, INT4 ng);
int output_lists_coverage(INT4 *nbeta, INT4 *ngamma);
int corr_bisp_sky(double *bispp, double *modsq, INT4 *nx, INT4 *ny,
                  INT4 *bisp_dim, float *sky_level, INT4 *nbeta, 
                  INT4 *ngamma);
int photon_corr(double *bispp, double *modsq, double *snrm,
                INT4 *nx, INT4 *ny, float *xphot, INT4 *nbeta, INT4 *ngamma,
                INT4 *photon_correction);
int BISPEC3(double *re, double *im, double *modsq, double *snrm,
            INT4 *nx, INT4 *ny, double *bispp, INT4 *ir, 
            INT4 *nbeta, INT4 *ngamma);
int BISPEC3_SINGLE(float *re, float *im, float *modsq, float *snrm,
            INT4 *nx, INT4 *ny, float *bispp, INT4 *ir,
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
int bisp2D_to_2D_image(float *out_image, float *modsq, INT4 *nx, INT4 *ny, 
                       float *bisp, INT4 *nx_bisp, INT4 *nbeta, INT4 *ngamma, 
                       INT4 *normalize, INT4 *iopt);
int bisp1D_to_2D_image(float *out_image, INT4 *nx_out, INT4 *ny_out, 
                       float *modsq, INT4 *nx_modsq, float *bisp, INT4 *nx_bisp,
                       INT4 *nbeta, INT4 *ngamma, INT4 *line_to_output, 
                       INT4 *normalize, INT4 *iopt);
#endif
