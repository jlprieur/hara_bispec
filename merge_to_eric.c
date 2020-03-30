/********************************************************************
* Combination of sets of bispectrum and modulus in Eric Anterrieur's format 
*
* Calls: 
* (in jlp_cover_mask.c)
*  
* Warning: frequencies can be ordered slightly differently in input
* and output spectral lists (since convention is not always followed
* in input lists, whereas it is in output list)
*
* JLP
* Version 02-12-93
*******************************************************************/
#include <stdio.h>
#include <math.h>
#include <jlp_ftoc.h>

#define DEBUG

/* Maximum number of data sets: */
#define NSMAX 5 

main()
{
float *modulus, *bisp_real, *bisp_imag;
int *spec_red, *bisp_red, *from_set_to_merged;
int *u_freq, *v_freq, *k_bifreq, *l_bifreq, *m_bifreq, *ngt0, *ngt;
int ir, nbeta_max, ngamma_max, nbeta, ngamma, nbeta0[NSMAX], ngamma0[NSMAX];
int pntr1, isize, nsets;
register int i, iset;
char *pc, buffer[81];
char out_bisp_list[61], out_spec_list[61], out_dimension_file[61];
char out_bisp_data[61], out_mod_data[61], out_ext[41];
char in_bisp_list[61], in_spec_list[61], in_dimension_file[61];
char in_bisp_data[61], in_mod_data[61], in_ext[NSMAX][41];

printf(" merge_to_eric           -- Version 02-12-93 --\n");
printf(" To combine many sets of masked data\n");

printf(" Number of data sets (maxi=%d) : ",NSMAX); 
scanf("%d",&nsets);
 if(nsets > NSMAX)
   {
   printf("merge_to_eric/Fatal error too many sets (nsets=%d > NSMAX=%d)\n",
            nsets,NSMAX);
   exit(-1);
   }

/* Input of input file extensions */
for(iset = 0; iset < nsets; iset++)
  {
  printf(" Input file extension of data set #%d =< ",iset); 
  gets(in_ext[iset]); if(in_ext[iset][0] == '\0')gets(in_ext[iset]);
/* Replaces the first blank character encountered by zero (end of string): */
  pc = in_ext[iset];
  while(*pc != ' ' && *pc) pc++;
  *pc = '\0'; 
  }

  printf(" Output file extension =< "); 
  gets(out_ext);

/***************************************************************/
JLP_BEGIN();

  printf(" Maximum radius of uv-coverage (IR) in pixels =< "); 
  gets(buffer); sscanf(buffer,"%d",&ir);
#ifdef DEBUG
  printf(" OK, nsets=%d ir=%d \n",nsets,ir);
#endif

/****************************************************************
* Read dimensions (to fill arrays nbeta0 and ngamma0) 
*****************************************************************/
nbeta_max = 0;
ngamma_max = 0;
for(iset = 0; iset < nsets; iset++)
  {
   sprintf(in_dimension_file,"dimensions%s.txt",in_ext[iset]);
   read_dimen(in_dimension_file,&nbeta0[iset],&ngamma0[iset]);
   nbeta_max += nbeta0[iset];
   ngamma_max += ngamma0[iset];
  }
printf(" nbeta_max = %d, ngamma_max = %d\n",nbeta_max,ngamma_max);

/****************************************************************
* Merged spectral list
*****************************************************************/
/* Get memory space: */
   isize = (nbeta_max + 1) * sizeof(int);
   JLP_GVM(&ngt0,&isize);
   JLP_GVM(&ngt,&isize);
   JLP_GVM(&u_freq,&isize);
   JLP_GVM(&v_freq,&isize);

/* In a first step, add the contribution of the sets one by one
   to compute merged spectral list: u_freq, v_freq
 */
nbeta = 0;
for(iset = 0; iset < nsets; iset++)
  {
   printf(" Building spectral list with data set #%d (nbeta0=%d)\n",
            iset,nbeta0[iset]);

/* Spectral list: */
   sprintf(in_spec_list,"list_spectre%s.txt",in_ext[iset]);
   merge_spec_list(in_spec_list,u_freq,v_freq,nbeta0[iset],&nbeta);
  }

printf(" Merged spectral list has nbeta=%d frequencies \n",nbeta);

/****************************************************************
* Merged bispectral list
*****************************************************************/
/* Get memory space: */
   isize = ngamma_max * sizeof(int);
   JLP_GVM(&k_bifreq,&isize);
   JLP_GVM(&l_bifreq,&isize);
   JLP_GVM(&m_bifreq,&isize);
   isize = nbeta * sizeof(int) * nsets;
   JLP_GVM(&from_set_to_merged,&isize);

/* Initialization of ngt and from_set_to_merged: */
for (i = 0; i < nbeta; i++) ngt[i] = 0;
for (i = 0; i < nbeta * nsets; i++) from_set_to_merged[i] = -1;

/* In a second step, since the merged spectral list is correct
   determine the correspondance
   between set index and merged index of the spectral list,
   Then build the bispectral list.
   Build "merged" arrays: k_bifreq, l_bifreq, m_bifreq, ngt
 */
ngamma = 0; 
for(iset = 0; iset < nsets; iset++)
  {
   printf(" Building bispectral list with data set #%d (ngamma0=%d)\n",
            iset,ngamma0[iset]);

/* Spectral and bispectral lists: */
   sprintf(in_spec_list,"list_spectre%s.txt",in_ext[iset]);
   sprintf(in_bisp_list,"list_bispectre%s.txt",in_ext[iset]);
   merge_bispec_list(in_spec_list,in_bisp_list,
                u_freq,v_freq,k_bifreq,l_bifreq,m_bifreq,ngt0,ngt,
                &from_set_to_merged[iset*nbeta],
                nbeta0[iset],ngamma0[iset],nbeta,&ngamma);
  }


printf(" Merged bispectral list has ngamma=%d bifrequencies \n",ngamma);

/****************************************************************
* Redundancy arrays 
*****************************************************************/
/* Get memory space for redundancy arrays: */
   isize = nbeta * sizeof(int) * nsets;
   JLP_GVM(&spec_red,&isize);
   isize = ngamma * sizeof(int) * nsets;
   JLP_GVM(&bisp_red,&isize);

/* Initialization of arrays: */
for (i = 0; i < nbeta * nsets; i++) spec_red[i] = 0;
for (i = 0; i < ngamma * nsets; i++) bisp_red[i] = 0;

/* In a third step, since the merged spectral and bispectral lists are correct
   and since the correspondance between set index and merged index 
   of the spectral list is determined, we compute
   the spectral and bispectral redundancies
   of the individual sets relative to the merged spectral and bispectral lists
   Build "merged" arrays: spec_red and bisp_red 
 */
for(iset = 0; iset < nsets; iset++)
  {
/* Spectral and bispectral lists: */
   sprintf(in_spec_list,"list_spectre%s.txt",in_ext[iset]);
   redun_spec(in_spec_list,&from_set_to_merged[iset*nbeta],
              &spec_red[iset*nbeta],ngt0,nbeta0[iset],nbeta);
   sprintf(in_bisp_list,"list_bispectre%s.txt",in_ext[iset]);
   redun_bispec(in_bisp_list,&bisp_red[iset*ngamma],
                  k_bifreq,l_bifreq,m_bifreq,ngt0,ngt,
                &from_set_to_merged[iset*nbeta],
                nbeta0[iset],ngamma0[iset],nbeta,ngamma);
  }

/* Output spectral, bispectral list in Eric's format */
   sprintf(out_spec_list,"list_spectre%s.txt",out_ext);
   write_spec_list(out_spec_list,spec_red,u_freq,v_freq,ngt,
                   nsets,nbeta,nbeta,ngamma);

   sprintf(out_bisp_list,"list_bispectre%s.txt",out_ext);
   write_bispec_list(out_bisp_list,bisp_red,
                     k_bifreq,l_bifreq,m_bifreq,nsets,ngamma,ngamma);

   sprintf(out_dimension_file,"dimensions%s.txt",out_ext);
   write_dimen(out_dimension_file,nsets,nbeta,ngamma);

/****************************************************************
* Modulus data
*****************************************************************/
/* Get memory space (less than in previous loop since we know nbeta now): */
   isize = nbeta * sizeof(float);
   JLP_GVM(&modulus,&isize);
/* Initialization: */
   for (i = 0; i < nbeta; i++) modulus[i] = 0.;

/* Build the whole set of modulae using redundancy information 
*  and store it in "modulus": */ 
  for(iset = 0; iset < nsets; iset++)
    {
     sprintf(in_mod_data,"module_spectre_exp%s.txt",in_ext[iset]);
     read_modulus(in_mod_data,modulus,&from_set_to_merged[iset*nbeta],
                  nbeta0[iset],nbeta);
    }

/* Now output modulus to Eric's format ASCII file: */
   sprintf(out_mod_data,"module_spectre_exp%s.txt",out_ext);
   merge_modulus(out_mod_data,modulus,spec_red,nbeta,nbeta,nsets);

/* Free memory: */
   JLP_FVM(&modulus);
   JLP_FVM(&spec_red);

/****************************************************************
* Bispectrum data
*****************************************************************/

/* Get memory space: */
   isize = ngamma * sizeof(float) * nsets;
   JLP_GVM(&bisp_real,&isize);
   JLP_GVM(&bisp_imag,&isize);
/* Initialization: */
for (i = 0; i < ngamma * nsets; i++) 
   {bisp_real[i] = 0.; bisp_imag[i] = 0.;}

/* Build the whole set of bispectra using redundancy information 
*  and store it in "bisp_real" and "bisp_imag": */ 
  for(iset = 0; iset < nsets; iset++)
    {
     sprintf(in_spec_list,"list_spectre%s.txt",in_ext[iset]);
     sprintf(in_bisp_list,"list_bispectre%s.txt",in_ext[iset]);
     sprintf(in_bisp_data,"phaseur_bispectre_exp%s.txt",in_ext[iset]);
     read_bispect(in_bisp_data,in_spec_list,in_bisp_list,
                  u_freq,v_freq,k_bifreq,l_bifreq,m_bifreq,
                  &from_set_to_merged[iset*nbeta],ngt0,ngt,
                  bisp_real,bisp_imag,nbeta0[iset],nbeta,ngamma0[iset],ngamma);
    }


/* Now output bispectrum to Eric's format ASCII file: */
   sprintf(out_bisp_data,"phaseur_bispectre_exp%s.txt",out_ext);
   merge_bisp_data(out_bisp_data,bisp_real,bisp_imag,bisp_red,
               ngamma,ngamma,nsets);

/* Free memory: */
   JLP_FVM(&bisp_real);
   JLP_FVM(&bisp_imag);
   JLP_FVM(&bisp_red);

/* End of program: */
exit(0);
}
/*******************************************************************
* Merge a set of bispectra, according to bispectral redundancy 
* Output bispectrum phasor to Eric's format ASCII file: 
*
* Input:
* bisp_ascii_file: file name for output bispectrum 
*
* Output:
*******************************************************************/
int merge_bisp_data(out_bisp_data,bisp_real,bisp_imag,bisp_red,
               ngamma,ngamma_max,nsets)
float bisp_real[], bisp_imag[];
int bisp_red[];
int ngamma, ngamma_max, nsets;
char *out_bisp_data;
{
FILE   *fp;
float r0, i0, work;
int inull;
register int ng;

/***** Output of bispectrum to a file ********************/
if ((fp = fopen(out_bisp_data,"w")) == NULL)
   {
   printf(" merge_bisp/Fatal error opening output file: %s \n",
            out_bisp_data);
   exit(-1);
   }

inull = 0;
for (ng = 0; ng < ngamma; ng++)
 {
/* Normalizes bispectral terms: */
       r0 = bisp_real[ng]; 
       i0 = bisp_imag[ng]; 
       work = r0*r0 + i0*i0;
       if(work == 0.) 
          inull++;
       else
          {
           work = sqrt((double)work);
           r0 /= work; i0 /= work;
          }
  fprintf(fp,"%f %f\n",r0,i0);

/* End of loop on bispectral terms (ng) */
 }

fclose(fp);

/* Displays a warning message: */
   if(inull) printf("merge_bisp/Warning: %d null values! (ngamma=%d)\n",
                     inull,ngamma);

return(0);
}
/*************************************************************
* read_dimen
* To get nbeta0 and ngamma0 for a "dimension" file
*
* INPUT:
* in_dimension_file: file name
*
* OUTPUT:
* nbeta0, ngamma0
*
**************************************************************/
int read_dimen(in_dimension_file,nbeta0,ngamma0)
int *nbeta0, *ngamma0;
char *in_dimension_file;
{
FILE *fp;
int n_net, i;

/***** Opening file *********************/
if ((fp = fopen(in_dimension_file,"r")) == NULL)
  { printf(" read_dimen/Fatal error opening input file: %s \n",
            in_dimension_file);
    exit(-1);
  }

/***** number of networks **********************/
fscanf(fp,"%d\n",&n_net);
if(n_net > 1)
  {
   printf(" read_dimen/Fatal error, more than one network in input data set (file: %s)\n",
            in_dimension_file);
   exit(-1);
  }
/* Skip network information: */
for (i=0; i<n_net; i++)  fscanf(fp,"%*d ");

/***** Number of spectral frequencies ********/
fscanf(fp,"\n%d\n",nbeta0);

/***** Number of bispectral bi-frequencies ******/
fscanf(fp,"%d\n",ngamma0);

fclose(fp);

#ifdef DEBUG
  printf("read_dimen: %s successfully read nbeta0=%d ngamma0=%d\n",
          in_dimension_file,*nbeta0,*ngamma0);
#endif
}
/*******************************************************************
* Output spectral list in Eric's format
*
* INPUT:
* out_spec_list
* spec_red: spectral redundancy
*
*******************************************************************/
int write_spec_list(out_spec_list,spec_red,u_freq,v_freq,ngt,nsets,
                    nbeta_max,nbeta,ngamma)
char *out_spec_list;
int spec_red[], u_freq[], v_freq[], ngt[];
int nbeta_max,  nsets, nbeta;
{
FILE   *fp;
int redund;
register int nb, i;

/***** Output of spectral list to a file ********************/
if ((fp = fopen(out_spec_list,"w")) == NULL)
   {
   printf(" write_spec_list/Fatal error opening output file: %s\n",
            out_spec_list);
   exit(-1);
   }

 for (nb = 0; nb < nbeta; nb++)
   {
     fprintf(fp,"%d ",nb);
     fprintf(fp,"(%d,%d) ",u_freq[nb],v_freq[nb]);

/* Spectral redundancy: */
     redund = 0;
     for(i = 0; i < nsets; i++) redund += spec_red[nb + i * nbeta_max];
#ifdef DEBUG
     if(nb < 5) printf("write_spec_list/nb=%d, redund=%d\n",nb,redund);
#endif
     fprintf(fp,"%d ",redund);

/* Should write ngt[nb] index in bispectral list of last closure relation 
*  concerning spectral term #nb 
*  (all relations about nb are between 
     ngt[nb-1] (included) and ngt[nb] (excluded) )
* (Keep the following for consistency with Eric's programs...)
*/
        fprintf(fp,"%d\n",ngt[nb]);
   }

fclose(fp);
}
/*******************************************************************
* Output bispectral list in Eric's format
*
* INPUT:
* out_bisp_list
*
*******************************************************************/
int write_bispec_list(out_bisp_list,bisp_red,
                     k_bifreq,l_bifreq,m_bifreq,nsets,ngamma_max,ngamma)
char *out_bisp_list;
int bisp_red[], k_bifreq[], l_bifreq[], m_bifreq[];
int ngamma_max,  nsets, ngamma;
{
FILE   *fp;
register int ng, i;

/***** Output of bispectral list in a file ********************/
if ((fp = fopen(out_bisp_list,"w")) == NULL)
   {
   printf(" write_bispec_list/Fatal error opening output file: %s \n",
            out_bisp_list);
   exit(-1);
   }

for (ng = 0; ng < ngamma; ng++)
  {
  fprintf(fp,"%d\n",ng);

/* Bispectral redundancy (red_bispect) : */
  for (i = 0; i < nsets; i++)  
        fprintf(fp,"%d ",bisp_red[ng + ngamma_max * i]);

/* Bispectral indices k, l, m, of the three components of bispectral term: */
  fprintf(fp,"\n%d %d %d\n",k_bifreq[ng],l_bifreq[ng],m_bifreq[ng]);
  }

fclose(fp);
}
/*************************************************************
* write_dimen
* To write nsets, nbeta and ngamma in Eric's format
*
* INPUT:
* in_dimension_file: file name
* nstes, nbeta, ngamma
*
**************************************************************/
int write_dimen(out_dimension_file,nsets,nbeta,ngamma)
int nsets, nbeta, ngamma;
char *out_dimension_file;
{
FILE *fp;
register int i;

/***** Opening output file *********************/
if ((fp = fopen(out_dimension_file,"w")) == NULL)
  { printf(" read_dimen/Fatal error opening output file: %s \n",
            out_dimension_file);
    exit(-1);
  }

/* Number of networks */
fprintf(fp,"%d\n",nsets);

/* Number of telescopes per network */
for (i = 0; i < nsets; i++) fprintf(fp,"%d ",1);

/* Number of spectral frequencies */
fprintf(fp,"\n%d\n",nbeta);

/* Number of bispectral bi-frequencies */
fprintf(fp,"%d\n",ngamma);

fclose(fp);

return(0);
}
/*******************************************************************
* Output modulus to Eric's format ASCII file: 
*
* INPUT:
* mod_ascii_file: file name for modulus 
* modulus[nb]: sum of modulae for same frequency term #nb
*
*******************************************************************/
int merge_modulus(out_mod_data,modulus,spec_red,nbeta,nbeta_max,nsets)
char *out_mod_data;
float modulus[];
int spec_red[], nbeta, nbeta_max, nsets;
{
FILE   *fp;
float redund;
register int nb, i;

/***** Output of modulus to a file ********************/
if ((fp = fopen(out_mod_data,"w")) == NULL)
   {
   printf(" merge_modulus/Fatal error opening output file: %s \n",
            out_mod_data);
   exit(-1);
   }

/* Normalization of the modulus (with spectral redundancy) */
for (nb = 0; nb < nbeta; nb++)
 {
   redund = 0;
   for(i = 0; i < nsets; i++) redund += spec_red[nb + i * nbeta_max];
   if(redund > 0)
      modulus[nb] /= redund;
   else
/* Should never encounter this case, since we have rejected null redundancy
* in building the spectral list ....*/
      {
      printf("Merge_modulus/Fatal error: null redundancy for spectral term #%d\n",
             nb);
      exit(-1);
      }
/* Write modulus now: */
  fprintf(fp,"%f\n",modulus[nb]);
  }

fclose(fp);

return(0);
}

/*******************************************************************
* Read modulus of a data set (Eric's format ASCII file) 
* Load it by adding it at the right spot in merged spectral list
*
* Warning: frequencies can be ordered slightly differently in input
* and output spectral lists (since convention is not always followed
* in input lists, whereas it is in output list)
*
* INPUT:
* mod_ascii_file: file name for modulus 
* modulus[nb]: sum of modulae before adding present data set 
* nbeta0: number of frequencies of spectral list of current data set
* nbeta: number of frequencies of spectral list of merged data sets 
*
* OUTPUT:
* modulus[nb]: sum of modulae after addition of present data set 
*
*******************************************************************/
int read_modulus(in_mod_data,modulus,from_set_to_merged,nbeta0,nbeta)
float modulus[];
int from_set_to_merged[], nbeta0, nbeta;
char *in_mod_data;
{
FILE *fp;
float mod;
int nb_merged;
register int nb;

/***** Opening file *********************/
if ((fp = fopen(in_mod_data,"r")) == NULL)
  { printf(" read_modulus/Fatal error opening input file: %s \n",
            in_mod_data);
    exit(-1);
  }

for (nb = 0; nb < nbeta0; nb++)
  {
    nb_merged = from_set_to_merged[nb]; 
    if(nb_merged < 0)
     {
     printf("read_modulus/Fatal error: uncomplete spectral list\n");
     printf(" nb_set = %d\n",nb);
     exit(-1);
     }
    fscanf(fp,"%f\n",&mod);
    modulus[nb_merged] += mod;
  }

fclose(fp);
return(0);
}
/*******************************************************************
* Read modulus of a data set (Eric's format ASCII file)
* Load it by adding it at the right spot in merged spectral list
*
* Warning: frequencies can be ordered slightly differently in input
* and output spectral lists (since convention is not always followed
* in input lists, whereas it is in output list).
* Therefore, we read first the spectral list to get a correct ngt0
* and then the bispectral list to associate the bispectral data with
* the correct bispectral frequency in the merged bispectral list 
*
* INPUT:
* mod_ascii_file: file name for modulus
* bisp_real, bisp_imag: sum of bispectra before adding present data set
* ngamma0: number of bifrequencies of bispectral list of current data set
* ngamma: number of bifrequencies of bispectral list of merged data sets
*
* OUTPUT:
* bisp_real, bisp_imag: sum of bispectra after addition of present data set
*
*******************************************************************/
int read_bispect(in_bisp_data,in_spec_list,in_bisp_list,
                 u_freq,v_freq,k_bifreq,l_bifreq,m_bifreq,
                 from_set_to_merged,ngt0,ngt,
                 bisp_real,bisp_imag,nbeta0,nbeta,ngamma0,ngamma)
float bisp_real[], bisp_imag[];
int u_freq[], v_freq[], k_bifreq[], l_bifreq[], m_bifreq[];
int from_set_to_merged[], ngt0[], ngt[];
int nbeta0, nbeta, ngamma0, ngamma;
char *in_bisp_data, *in_spec_list, *in_bisp_list;
{
FILE *fp, *fp1;
float r0, i0;
int ifound, nb_merged, k_set, l_set, m_set, kk, ll, mm;
register int i, nb, ng;

/***** We first read the spectral list to get ngt0: */
if ((fp = fopen(in_spec_list,"r")) == NULL)
  { printf(" redun_spec/Fatal error opening input file: %s \n",
            in_spec_list); exit(-1);
  }

for (nb = 0; nb < nbeta0; nb++)
   {
    nb_merged = from_set_to_merged[nb];
/* Only interested in ngt0 */
       fscanf(fp,"%*d (%*d,%*d) %*d %d\n",&ngt0[nb]);
   }
fclose(fp);


/***** Opening bispectrum data file *********************/
if ((fp1 = fopen(in_bisp_data,"r")) == NULL)
  { printf(" read_bispect/Fatal error opening input file: %s \n",
            in_bisp_data);
    exit(-1);
  }

/***** Opening bispectral list *********************/
if ((fp = fopen(in_bisp_list,"r")) == NULL)
  { printf(" merge_bispec_list/Fatal error opening input file: %s \n",
            in_bisp_list);
    exit(-1);
  }

/****************************************************
* Main loop to associate the bispectral data with
* the correct bispectral frequency in the merged bispectral list 
*/

for (nb = 1; nb < nbeta0; nb++)
  {
   nb_merged = from_set_to_merged[nb];

   if(nb_merged > 0)
      {
      for(ng = ngt0[nb-1]; ng < ngt0[nb]; ng++)
        {
/* Read bispectral data: */
         fscanf(fp1,"%f %f\n",&r0,&i0);

/* Bifrequency number (in bispectral list): */
         fscanf(fp,"%*d\n");
/* Since only one network (input data set), I neglect the following loop:
*        for (i=0; i<n_networks; r++)  fscanf(fp,"%*d ");
* and put instead:
*/
         fscanf(fp,"%*d ");

/* Only interested in spectral indices of the three vectors */
         fscanf(fp,"\n%d %d %d\n",&k_set,&l_set,&m_set);

/* Check if bifrequency is already there
* and update k_bifreq, l_bifreq, m_bifreq, ngt and ngamma if needed */
         kk = from_set_to_merged[k_set];
         ll = from_set_to_merged[l_set];
         mm = from_set_to_merged[m_set];

/* Check that spectral frequencies are in "merged" spectral list: */
      if(kk < 0 || ll < 0 || mm < 0)
         {
         printf("read_bispec/Fatal error, uncomplete \"merged\" spectral list\n");
         printf(" ng = %d nb=%d \n",ng,nb);
         printf(" k_set = %d kk=%d, l_set=%d, ll=%d, m_set=%d mm=%d \n",
                  k_set,kk,l_set,ll,m_set,mm);
         exit(-1);
         }

      ifound = 0;
      for(i = ngt[nb_merged-1]; i < ngt[nb_merged]; i++)
         {
          if(k_bifreq[i] == kk && l_bifreq[i] == ll && m_bifreq[i] == mm)
            {ifound = 1; break;}
         }

/* Problem when kk,ll,mm not found in bispectral list: */
      if(!ifound)
         {
         printf("read_bispec/Fatal error, uncomplete \"merged\" bispectral list\n");
         printf(" kk=%d, ll=%d, mm=%d not found in bispectral list!\n",
                  kk,ll,mm);
         printf(" ngt[nb_merged-1] = %d ngt[nb_merged]=%d \n",
                  ngt[nb_merged-1],ngt[nb_merged]);
         printf(" nb_merged=%d nb_set = %d ng_set=%d\n",nb_merged,nb,ng);
         exit(-1);
         }

/* Now store the bispectral data at the right place: */
       bisp_real[i] += r0;
       bisp_imag[i] += i0;
     }

/* End of loop on ng */
    }
/* End of loop on nb */
  }

fclose(fp1);
fclose(fp);

/* Check that all bispectral bifrequencies have been read: */
 if(ng != ngamma0)
  {
  printf("read_bispec/Fatal error with %s, ng_max=%d, ngamma0=%d\n",
          in_bisp_list,ng,ngamma0);
  exit(-1);
  }

return(0);
}
/*******************************************************************
* merge_spec_list
* To build up the merged spectral list
*  
* Read spectral list of a data set (Eric's format ASCII file)
* Update the merged spectral list
*
* INPUT:
* in_spec_list, in_bisp_list: file names for spectral and bispectral lists 
* nbeta0, ngamma0: number of frequencies and bifrequencies of current data set
*
* OUTPUT:
* nbeta, ngamma: number of frequencies and bifrequencies of "merged" data set
*
*******************************************************************/
int  merge_spec_list(in_spec_list,u_freq,v_freq,nbeta0,nbeta)
int u_freq[], v_freq[];
int nbeta0, *nbeta;
char *in_spec_list;
{
int uu, vv, nb_merged, ifound;
FILE *fp;
register int nb, ng;

/***** Opening spectral list *********************/
if ((fp = fopen(in_spec_list,"r")) == NULL)
  { printf(" merge_spec_list/Fatal error opening input file: %s \n",
            in_spec_list);
    exit(-1);
  }

for (nb = 0; nb < nbeta0; nb++)
   {
/* Only interested in (u,v) components */
    fscanf(fp,"%*d (%d,%d) %*d %*d\n",&uu,&vv);

/* Check if frequency is already there 
* and update u_freq, v_freq and nbeta if needed 
* returns nb_merged, corresponding index in "merged list" */
     add_new_freq(uu,vv,nb,u_freq,v_freq,&nb_merged,nbeta,&ifound);
#ifdef DEBUG
      if(nb < 5)
        printf("nb=%d uu=%d vv=%d nb_merged=%d nbeta=%d ifound=%d\n",
               nb,uu,vv,nb_merged,*nbeta,ifound);
#endif
   }
fclose(fp);

#ifdef DEBUG
  printf("merge_spec_list: %s successfully read\n",in_spec_list);
#endif

return(0);
}
/*******************************************************************
* add_new_freq
* To build up the merged spectral list
* 
* Check that uu,vv is in the spectral list
* If not, update u_freq, v_freq and nbeta.
*
* INPUT:
* uu,vv: coordinates of spatial frequency to incorporate in "merged" list
* nb_set: index of spectral frequency (uu,vv) in the set list
* u_freq, v_freq: coordinates of spectral frequencies of current list 
* nbeta: current number of frequencies of the "merged" data set
*
* OUTPUT:
* nb_merged: index of spectral frequency (uu,vv) in the merged list
* nbeta: updated number of frequencies of the "merged" data set
* u_freq, v_freq: coordinates of spectral frequencies of updated list 
* ifound: =1 if already in current list, =0 otherwise
*
*******************************************************************/
int add_new_freq(uu,vv,nb_set,u_freq,v_freq,nb_merged,nbeta,ifound)
int u_freq[], v_freq[];
int uu, vv, nb_set, *nb_merged, *nbeta, *ifound;
{
int rad2, current_rad2;
register int nb, i;

/* Case of first call: */
if(*nbeta == 0)
  {
  u_freq[0] = uu; v_freq[0] = vv;
  *nb_merged = 0;
  *nbeta = 1;
  return(0);
  }

/* General case: */
*ifound = 0;
for (nb = 0; nb < *nbeta; nb++)
   {
   if(u_freq[nb] == uu && v_freq[nb] == vv)
     {
       *ifound = 1; *nb_merged = nb; break;
     }
   }

/* If found, exit: */
if(*ifound)return(0);

/* If not found, look for the right spot where to put the (uu,vv) vector: */
rad2 = uu * uu + vv * vv;
 
for (nb = 0; nb < *nbeta; nb++)
   {
     current_rad2 = u_freq[nb]*u_freq[nb] + v_freq[nb]*v_freq[nb];
/* Frequencies are sorted according to radius: */
     if(current_rad2 > rad2)
       {
        break;
       }
     else
       {
         if(rad2 == current_rad2)
         {
/* In case of same radius, sorted according to uu coordinate
* The first should be the one with the largest uu value,
*  but there may be some exceptions in the input lists 
* since we affect first the symmetrical in jlp_cover... 
* Anyway here we strictly adopt the rule:
*/
          while(u_freq[nb] > uu && nb < *nbeta) nb++; 
          break;
         }
       }
   }

 *nb_merged = nb;

/* Case before the end of the list: */
   if(nb < *nbeta)
    {
/* Now we have found the spot where to put the new (uu,vv) vector,
*  i.e. between nb and nb-1, i.e. should replace nb vector
* Therefore, we shift the remaining values of one step: */
     *nb_merged = nb + 1;
      for (i = *nbeta; i > nb; i--)
        {u_freq[i] = u_freq[i-1]; v_freq[i] = v_freq[i-1];}
    }


#ifdef DEBUG
  if(*nbeta > 1 && *nbeta < 5 && nb < *nbeta)
    {
      printf("nb=%d nbeta=%d, u_freq[nb] = %d v_freq[nb] = %d uu = %d, vv = %d\n",
               nb,*nbeta,u_freq[nb],v_freq[nb],uu,vv);
      if(nb < *nbeta - 1)
      printf("  u_freq[nb+1] = %d v_freq[nb+1] = %d \n",u_freq[nb+1],v_freq[nb+1]);
    }
#endif

/* Put new vector at the right spot and increment nbeta: */
u_freq[*nb_merged] = uu;
v_freq[*nb_merged] = vv;
(*nbeta)++;

return(0);
}
/*******************************************************************
* To build up the merged bispectral list
*
* Check that kk, ll, mm is in the bispectral list
* If not, update kk_bifreq, ll_bifreq, kk_bifreq and ngamma.
*
* INPUT:
* kk,ll,mm: coordinates of spatial bifrequency (in "merged spectral list")
            to incorporate in "merged" list
* (Bispectral indices k, l, m, of the three components of bispectral term)
* k_bifreq, l_bifreq, m_bifreq: coordinates of bispectral bifrequencies 
*                               of current list
* ngamma: current number of bifrequencies of the "merged" data set
*
* OUTPUT:
* ngamma: updated number of bifrequencies of the "merged" data set
* k_bifreq, l_bifreq, m_bifreq: coordinates of bispectral bifrequencies 
*                               of updated list
* ifound: =1 if already in current list, =0 otherwise
*
*******************************************************************/
int add_new_bifreq(kk,ll,mm,ng_set,nb_merged,k_bifreq,l_bifreq,m_bifreq,ngt,
                   nbeta,ngamma,ng_merged,ifound)
int kk, ll, mm, ng_set, nbeta, *ngamma, *ng_merged, *ifound;
int k_bifreq[], l_bifreq[], m_bifreq[], ngt[];
{
int rad2, current_rad2;
register int ng, i;

*ifound = 0;

/* Case of first call: */
if(*ngamma == 0)
  {
  k_bifreq[0] = kk; l_bifreq[0] = ll; m_bifreq[0] = mm;
  *ngamma = 1; 
  for (i = nb_merged; i < nbeta; i++) (ngt[i])++;
  return(0);
  }

/* General case: */
for (ng = ngt[nb_merged-1]; ng < ngt[nb_merged]; ng++)
   {
    if(k_bifreq[ng] == kk && l_bifreq[ng] == ll)
     {
       if(m_bifreq[ng] == mm)
          {
           *ifound = 1; *ng_merged = ng; break;
          }
       else
          {
          printf("add_new_bifreq/Fatal error k and l are the same, but not m!\n");
          printf("ng_set=%d, nb_merged=%d, k=%d l=%d m_bifreq=%d m=%d \n",
                  ng_set,nb_merged,kk,ll,m_bifreq[ng],mm);
          exit(-1);
          }
     }
}

/* If found, exit: */
if(*ifound)return(0);

/* If not found, simply put the (kk,ll,mm) vector at the end of present set 
* of relations concerning frequency #nb_merged:
*  i.e. between ngt[nb_merged]-1 and ngt[nb_merged]
*/

/* We shift the remaining values of one step: */
   for (i = *ngamma; i >= ngt[nb_merged]; i--)
     {
       k_bifreq[i] = k_bifreq[i-1];
       l_bifreq[i] = l_bifreq[i-1];
       m_bifreq[i] = m_bifreq[i-1];
     }
   *ng_merged = ngt[nb_merged];

/* Put new bifrequency coordinates at the right spot and increment ngamma: */
k_bifreq[*ng_merged] = kk;
l_bifreq[*ng_merged] = ll;
m_bifreq[*ng_merged] = mm;
(*ngamma)++;

/* We update "merged" ngt array: 
   -- ngt(nb_merged) increased since one new relation for frequency #nb_merged
   -- all subsequent nb indices increased because of shift
*/
  for (i = nb_merged; i < nbeta; i++) (ngt[i])++;

return(0);
}
/*******************************************************************
* merge_bispect_list
* To build up the merged bispectral list
*  
* Read spectral and bispectral lists of a data set (Eric's format ASCII file)
* Compute the new spectral and bispectral lists
*
* INPUT:
* in_spec_list, in_bisp_list: file names for spectral and bispectral lists 
* nbeta0, ngamma0: number of frequencies and bifrequencies of current data set
* from_set_to_merged: array of the indices of the spectral terms of current set
*                     in the merged spectral list
*                     (To ensure the corespondance for bispectral terms).
*
* OUTPUT:
* nbeta, ngamma: number of frequencies and bifrequencies of "merged" data set
*
*******************************************************************/
int  merge_bispec_list(in_spec_list,in_bisp_list,
                  u_freq,v_freq,k_bifreq,l_bifreq,m_bifreq,ngt0,ngt,
                  from_set_to_merged,nbeta0,ngamma0,nbeta,ngamma)
int u_freq[], v_freq[], from_set_to_merged[];
int k_bifreq[], l_bifreq[], m_bifreq[], ngt0[], ngt[]; 
int nbeta0, ngamma0, nbeta, *ngamma;
char *in_spec_list, *in_bisp_list;
{
int uu, vv, kk, ll, mm, k_set, l_set, m_set, nb_merged, ng_merged, ifound;
int nbeta1; 
FILE *fp;
register int nb, ng;

/***** Opening spectral list *********************/
if ((fp = fopen(in_spec_list,"r")) == NULL)
  { printf(" merge_bispec_list/Fatal error opening input file: %s \n",
            in_spec_list);
    exit(-1);
  }

nbeta1 = nbeta;
for (nb = 0; nb < nbeta0; nb++)
   {
/* Interested in (u,v) components and number of relations (ngt actually...)
 ngt[nb] = index of the last bispectral relation in bispectral list concerning 
           spectral frequency #nb
*/
  fscanf(fp,"%*d (%d,%d) %*d %d\n",&uu,&vv,&ngt0[nb]);

/* Check if frequency is already there 
* and update u_freq, v_freq and nbeta if needed 
* returns nb_merged, corresponding index in "merged list" */
     add_new_freq(uu,vv,nb,u_freq,v_freq,&nb_merged,&nbeta1,&ifound);
     if(!ifound)
         {
         printf("merge_bispec_list/Fatal error, frequency #%d not found in merged list!\n",
               nb);
         exit(-1);
         }
     from_set_to_merged[nb] = nb_merged;

   }
fclose(fp);

#ifdef DEBUG
  printf("\n merge_bispec_list: %s processed\n",in_spec_list);
#endif

/***** Opening bispectral list *********************/
if ((fp = fopen(in_bisp_list,"r")) == NULL)
  { printf(" merge_bispec_list/Fatal error opening input file: %s \n",
            in_bisp_list);
    exit(-1);
  }

/****************************************************
* Main loop to build bispectral list
*/

for (nb = 1; nb < nbeta0; nb++)
  {
   nb_merged = from_set_to_merged[nb];
   
   if(nb_merged > 0)
      {
      for(ng = ngt0[nb-1]; ng < ngt0[nb]; ng++)
        {
/* Bifrequency number: */
         fscanf(fp,"%*d\n");
/* Since only one network (input data set), I neglect the following loop:
*        for (i=0; i<n_networks; r++)  fscanf(fp,"%*d ");
* and put instead:
*/
         fscanf(fp,"%*d ");

/* Only interested in spectral indices of the three vectors */
         fscanf(fp,"\n%d %d %d\n",&k_set,&l_set,&m_set);

/* Check if bifrequency is already there 
* and update k_bifreq, l_bifreq, m_bifreq, ngt and ngamma if needed */
         kk = from_set_to_merged[k_set];
         ll = from_set_to_merged[l_set];
         mm = from_set_to_merged[m_set];

/* Check that spectral frequencies are in "merged" spectral list: */
      if(kk < 0 || ll < 0 || mm < 0)
         {
         printf("merge_bispec_list/Fatal error, uncomplete \"merged\" spectral list\n");
         printf(" ng = %d nb=%d \n",ng,nb);
         printf(" k_set = %d kk=%d, l_set=%d, ll=%d, m_set=%d mm=%d \n",
                  k_set,kk,l_set,ll,m_set,mm);
         exit(-1);
         }

          add_new_bifreq(kk,ll,mm,ng,nb_merged,k_bifreq,l_bifreq,m_bifreq,
                     ngt,nbeta,ngamma,&ng_merged,&ifound);
#ifdef DEBUG
         if(nb < 5)
          printf("nb=%d ng=%d kk=%d ll=%d mm=%d ngamma=%d ifound=%d\n",
                  nb,ng,kk,ll,mm,*ngamma,ifound);
#endif
/* End of loop on ng */
       } 
/* End of case: nb_merged > -1 */
      }
/* End of loop on nb */
  }

fclose(fp);

#ifdef DEBUG
  printf("merge_bispec_list: %s processed\n",in_bisp_list);
#endif

return(0);
}
/*******************************************************************
* redun_spec
* To determine spectral redundancy 
*  
* Read spectral list of a data set (Eric's format ASCII file)
* Update the redundancy spectral list
*
* INPUT:
* in_spec_list: file name for spectral list 
* nbeta0: number of frequencies of current data set
* from_set_to_merged: array of the indices of the spectral terms of current set
*                     in the merged spectral list
*                     (To ensure the corespondance for bispectral terms).
*
* OUTPUT:
* nbeta: number of frequencies of "merged" data set
*
*******************************************************************/
int redun_spec(in_spec_list,from_set_to_merged,spec_red,ngt0,nbeta0,nbeta)
int from_set_to_merged[], spec_red[], ngt0[];
int nbeta0, nbeta;
char *in_spec_list;
{
int nb_merged, ired;
FILE *fp;
register int nb;

/***** Opening spectral list *********************/
if ((fp = fopen(in_spec_list,"r")) == NULL)
  { printf(" redun_spec/Fatal error opening input file: %s \n",
            in_spec_list);
    exit(-1);
  }

for (nb = 0; nb < nbeta0; nb++)
   {
    nb_merged = from_set_to_merged[nb];
/* Check that this frequency has been incorporated in "merged" list */
    if(nb_merged > -1)
       {
/* Only interested in spectral redundancy and ngt0 */
       fscanf(fp,"%*d (%*d,%*d) %d %d\n",&ired,&ngt0[nb]);
       spec_red[nb_merged] += ired;
       }
    else
       {
         printf("redun_spec/Fatal error, uncomplete \"merged\" spectral list\n");
         printf(" nb_set = %d\n",nb);
         exit(-1);
       }
   }
fclose(fp);

/* Check that all spectral frequencies have been read: */
 if(nb != nbeta0)
  {
  printf("redun_spec/Fatal error with %s, nb_max=%d, nbeta0=%d\n",
          in_spec_list,nb,nbeta0);
  exit(-1);
  }

#ifdef DEBUG
  printf("redun_spec: %s processed, nb_max=%d, nbeta0=%d\n",
          in_spec_list,nb,nbeta0);
  for (nb = 0; nb < 5; nb++)
     printf("redun_spec: nb_merged=%d, spec_red_set=%d\n",nb,spec_red[nb]);
#endif
 
return(0);
}
/*******************************************************************
* redun_bispec
* To determine bispectral redundancy 
*  
* Read spectral and bispectral lists of a data set (Eric's format ASCII file)
* Compute the new spectral and bispectral lists
*
* INPUT:
* in_bisp_list: file name for bispectral lists 
* ngamma0: number of bifrequencies of current data set
* from_set_to_merged: array of the indices of the spectral terms of current set
*                     in the merged spectral list
*                     (To ensure the corespondance for bispectral terms).
*
* OUTPUT:
* nbeta, ngamma: number of frequencies and bifrequencies of "merged" data set
*
*******************************************************************/
int redun_bispec(in_bisp_list,bisp_red,
                  k_bifreq,l_bifreq,m_bifreq,ngt0,ngt,
                  from_set_to_merged,nbeta0,ngamma0,nbeta,ngamma)
int bisp_red[], from_set_to_merged[];
int k_bifreq[], l_bifreq[], m_bifreq[], ngt0[], ngt[]; 
int nbeta0, ngamma0, nbeta, *ngamma;
char *in_bisp_list;
{
int kk, ll, mm, k_set, l_set, m_set, nb_merged, ng_merged, ifound;
FILE *fp;
register int nb, ng, i;

/***** Opening bispectral list *********************/
if ((fp = fopen(in_bisp_list,"r")) == NULL)
  { printf(" redun_bispec/Fatal error opening input file: %s \n",
            in_bisp_list);
    exit(-1);
  }

/****************************************************
* Main loop to build bispectral list
*/

for (nb = 1; nb < nbeta0; nb++)
  {
   nb_merged = from_set_to_merged[nb];
   if(nb_merged < 0)printf("redun_bispec/Error: nb_merged=%d \n",nb_merged);
   
   for(ng = ngt0[nb-1]; ng < ngt0[nb]; ng++)
     {
/* Bifrequency number: */
      fscanf(fp,"%*d\n");
/* Since only one network (input data set), I neglect the following loop:
*     for (i=0; i<n_networks; r++)  fscanf(fp,"%*d ");
* and put instead:
*/
      fscanf(fp,"%*d ");

/* Only interested in spectral indices of the three vectors */
      fscanf(fp,"\n%d %d %d\n",&k_set,&l_set,&m_set);

/* Check if bifrequency is already there 
* and update k_bifreq, l_bifreq, m_bifreq, ngt and ngamma if needed */
      kk = from_set_to_merged[k_set];
      ll = from_set_to_merged[l_set];
      mm = from_set_to_merged[m_set];

/* Check that spectral frequencies are in "merged" spectral list: */
      if(kk < 0 || ll < 0 || mm < 0)
         {
         printf("redun_bispec/Fatal error, uncomplete \"merged\" spectral list\n");
         printf(" ng = %d nb=%d \n",ng,nb);
         printf(" k_set = %d kk=%d, l_set=%d, ll=%d, m_set=%d mm=%d \n",
                  k_set,kk,l_set,ll,m_set,mm);
         exit(-1);
         }

      ifound = 0;
      for(i = ngt[nb_merged-1]; i < ngt[nb_merged]; i++)
         {
          if(k_bifreq[i] == kk && l_bifreq[i] == ll && m_bifreq[i] == mm)
            {(bisp_red[i])++; ifound = 1; break;}
         }
/* Problem when kk,ll,mm not found in bispectral list: */
      if(!ifound)
         {
         printf("redun_bispec/Fatal error, uncomplete \"merged\" bispectral list\n");
         printf(" kk=%d, ll=%d, mm=%d not found in bispectral list!\n",
                  kk,ll,mm);
         printf(" ngt[nb_merged-1] = %d ngt[nb_merged]=%d \n",
                  ngt[nb_merged-1],ngt[nb_merged]);
         printf(" nb_merged=%d nb_set = %d ng_set=%d\n",nb_merged,nb,ng);
         exit(-1);
         }

/* End of loop on ng */
    } 
/* End of loop on nb */
  }
fclose(fp);

/* Check that all bispectral bifrequencies have been read: */
 if(ng != ngamma0)
  {
  printf("redun_bispec/Fatal error with %s, ng_max=%d, ngamma0=%d\n",
          in_bisp_list,ng,ngamma0);
  exit(-1);
  }

#ifdef DEBUG
  printf("redun_bispec: %s processed, ng_max=%d, ngamma0=%d\n",
          in_bisp_list,ng,ngamma0);
  for (ng = 0; ng < 5; ng++)
     printf("redun_bispec: ng_merged=%d, bisp_red_set=%d\n",ng,bisp_red[ng]);
#endif

return(0);
}
