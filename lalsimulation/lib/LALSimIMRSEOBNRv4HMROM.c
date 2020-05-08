
/*
 * Copyright (C) 2019 Michael Puerrer, John Veitch, Roberto Cotesta, Sylvain Marsat  
 * Copyright (C) 2014, 2015, 2016 Michael Puerrer, John Veitch
 *  Reduced Order Model for SEOBNR
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdbool.h>
#include <alloca.h>
#include <string.h>
#include <libgen.h>


#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_poly.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/Sequence.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALSimSphHarmMode.h>

#define NMODES 5
#define f_hyb_ini 0.003
#define f_hyb_end 0.004

#ifdef LAL_HDF5_ENABLED
#include <lal/H5FileIO.h>
static const char ROMDataHDF5[] = "SEOBNRv4HMROM.hdf5";
static const INT4 ROMDataHDF5_VERSION_MAJOR = 1;
static const INT4 ROMDataHDF5_VERSION_MINOR = 0;
static const INT4 ROMDataHDF5_VERSION_MICRO = 0;
#endif

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSEOBNRROMUtilities.c"

#include <lal/LALConfig.h>
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
#endif



#ifdef LAL_PTHREAD_LOCK
static pthread_once_t SEOBNRv4HMROM_is_initialized = PTHREAD_ONCE_INIT;
#endif


/*************** mode array definitions ******************/

const char mode_array[NMODES][3] = 
{"22",
"33",
"21",
"44",
"55"};

const UINT4 lmModes[NMODES][2] = {{2, 2}, {3, 3}, {2, 1}, {4, 4}, {5, 5}};

/*************** constant phase shifts ******************/

const REAL8 const_phaseshift_lm[NMODES] = {0.,-LAL_PI/2.,LAL_PI/2.,LAL_PI,LAL_PI/2.};
/* This constant phaseshift is coming from the relations between the phase
* of the lm mode phi_lm and the orbital phase phi. At leading order phi_lm = -m*phi + c_lm
* where c_lm are these constant phaseshift. They can be derived from the PN expression of the modes
* see Blanchet, L. Living Rev. Relativ. (2014) 17: 2. https://doi.org/10.12942/lrr-2014-2 Eqs(327-329)
*/

/*************** constant factors used to define fmax for each mode ******************/

const REAL8 const_fmax_lm[NMODES] = {1.7, 1.55, 1.7, 1.35, 1.25};

/*************** starting geometric frequency 22 mode non-hybridized ROM ***********************/

REAL8 Mf_low_22 = 0.0004925491025543576;

/*************** type definitions ******************/

typedef struct tagSEOBNRROMdataDS_coeff
{
  gsl_vector* c_real;
  gsl_vector* c_imag;
  gsl_vector* c_phase;
} SEOBNRROMdataDS_coeff;


struct tagSEOBNRROMdataDS_submodel
{
  gsl_vector* cvec_real;     // Flattened Re(c-mode) projection coefficients
  gsl_vector* cvec_imag;     // Flattened Im(c-mode) projection coefficients
  gsl_vector* cvec_phase;    // Flattened orbital phase projection coefficients
  gsl_matrix *Breal;         // Reduced SVD basis for Re(c-mode)
  gsl_matrix *Bimag;         // Reduced SVD basis for Im(c-mode)
  gsl_matrix *Bphase;        // Reduced SVD basis for orbital phase
  int nk_cmode;              // Number frequency points for c-mode
  int nk_phase;              // Number of frequency points for orbital phase
  gsl_vector *gCMode;        // Sparse frequency points for c-mode
  gsl_vector *gPhase;        // Sparse frequency points for orbital phase
  gsl_vector *qvec;          // B-spline knots in q
  gsl_vector *chi1vec;       // B-spline knots in chi1
  gsl_vector *chi2vec;       // B-spline knots in chi2
  int ncx, ncy, ncz;         // Number of points in q, chi1, chi2
  double q_bounds[2];        // [q_min, q_max]
  double chi1_bounds[2];     // [chi1_min, chi1_max]
  double chi2_bounds[2];     // [chi2_min, chi2_max]
};
typedef struct tagSEOBNRROMdataDS_submodel SEOBNRROMdataDS_submodel;

struct tagSEOBNRROMdataDS
{
  UINT4 setup;
  SEOBNRROMdataDS_submodel* hqhs;
  SEOBNRROMdataDS_submodel* hqls;
  SEOBNRROMdataDS_submodel* lqhs;
  SEOBNRROMdataDS_submodel* lqls;
  SEOBNRROMdataDS_submodel* lowf;
};
typedef struct tagSEOBNRROMdataDS SEOBNRROMdataDS;

static SEOBNRROMdataDS __lalsim_SEOBNRv4HMROMDS_data[NMODES];

typedef int (*load_dataPtr)(const char*, gsl_vector *, gsl_vector *, gsl_matrix *, gsl_matrix *, gsl_vector *);

typedef struct tagSplineData
{
  gsl_bspline_workspace *bwx;
  gsl_bspline_workspace *bwy;
  gsl_bspline_workspace *bwz;
} SplineData;

/**************** Internal functions **********************/

UNUSED static bool SEOBNRv4HMROM_IsSetup(UINT4);
UNUSED static void SEOBNRv4HMROM_Init_LALDATA(void);
UNUSED static int SEOBNRv4HMROM_Init(const char dir[],UINT4);
UNUSED static int SEOBNRROMdataDS_Init(SEOBNRROMdataDS *romdata, const char dir[], UINT4);
UNUSED static void SEOBNRROMdataDS_Cleanup(SEOBNRROMdataDS *romdata);

UNUSED static int SEOBNRROMdataDS_Init_submodel(
  UNUSED SEOBNRROMdataDS_submodel **submodel,
  UNUSED const char dir[],
  UNUSED const char grp_name[],
  UINT4 index_mode
);
UNUSED static void SEOBNRROMdataDS_Cleanup_submodel(SEOBNRROMdataDS_submodel *submodel);
UNUSED static void SplineData_Destroy(SplineData *splinedata);
UNUSED static void SplineData_Init(
  SplineData **splinedata,
  int ncx,                // Number of points in q  + 2
  int ncy,                // Number of points in chi1 + 2
  int ncz,                // Number of points in chi2 + 2
  const double *qvec,     // B-spline knots in q
  const double *chi1vec,  // B-spline knots in chi1
  const double *chi2vec   // B-spline knots in chi2
);

/**
 * Core function for computing the ROM waveform.
 * Interpolate projection coefficient data and evaluate coefficients at desired (q, chi).
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
UNUSED static int SEOBNRv4HMROMCoreModes(
  SphHarmFrequencySeries **hlm_list, /**< Spherical modes frequency series for the waveform */
  REAL8 phiRef, /**< Orbital phase at reference time */
  REAL8 fRef, /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance, /**< Distance of source (m) */
  REAL8 Mtot_sec, /**< Total mass in seconds **/
  REAL8 q, /**< Mass ratio **/
  REAL8 chi1, /**< Dimensionless aligned component spin 1 */
  REAL8 chi2, /**< Dimensionless aligned component spin 2 */
  const REAL8Sequence *freqs, /**< Frequency points at which to evaluate the waveform (Hz) */
  REAL8 deltaF,
  /**< If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
  INT4 nk_max, /**< truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1 */
  UINT4 nModes, /**< Number of modes to generate */
  REAL8 sign_odd_modes /**< Sign of the odd-m modes, used when swapping the two bodies */
);
UNUSED static void SEOBNRROMdataDS_coeff_Init(SEOBNRROMdataDS_coeff **romdatacoeff, int nk_cmode, int nk_phase);
UNUSED static void SEOBNRROMdataDS_coeff_Cleanup(SEOBNRROMdataDS_coeff *romdatacoeff);

static size_t NextPow2(const size_t n);
static UINT8 Setup_EOBROM__std_mode_array_structure(LALValue *ModeArray, UINT4 nModes);
static INT8 Check_EOBROM_mode_array_structure(LALValue *ModeArray,UINT4 nModes);
static int SEOBROMComputehplushcrossFromhlm(
    COMPLEX16FrequencySeries
        *hplusFS, /**<< Output: frequency series for hplus, already created */
    COMPLEX16FrequencySeries
        *hcrossFS,   /**<< Output: frequency series for hplus, already created */
    LALValue *ModeArray, /**<< Input: ModeArray structure with the modes to include */
    SphHarmFrequencySeries
        *hlm,  /**<< Input: list with frequency series for each mode hlm */
    REAL8 inc,  /**<< Input: inclination */
    REAL8 phi   /**<< Input: phase */
);

static int TP_Spline_interpolation_3d(
  REAL8 q,                // Input: eta-value for which projection coefficients should be evaluated
  REAL8 chi1,               // Input: chi1-value for which projection coefficients should be evaluated
  REAL8 chi2,               // Input: chi2-value for which projection coefficients should be evaluated
  gsl_vector *cvec,         // Input: data for spline coefficients
  int nk,                   // number of SVD-modes == number of basis functions
  int nk_max,               // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
  int ncx,                  // Number of points in eta  + 2
  int ncy,                  // Number of points in chi1 + 2
  int ncz,                  // Number of points in chi2 + 2
  const double *qvec,     // B-spline knots in eta
  const double *chi1vec,    // B-spline knots in chi1
  const double *chi2vec,    // B-spline knots in chi2
  gsl_vector *c_out        // Output: interpolated projection coefficients
);
UNUSED static UINT8 SEOBNRv4HMROM_Select_HF_patch(REAL8 q, REAL8 chi1);
UNUSED static UINT8 SEOBNRv4HMROM_phase_sparse_grid(gsl_vector* phase_f, REAL8 q, REAL8 chi1,REAL8 chi2, const char* freq_range, INT4 nk_max);
UNUSED static UINT8 SEOBNRv4HMROM_cmode_sparse_grid(gsl_vector* creal_f, gsl_vector* cimag_f, REAL8 q, REAL8 chi1,REAL8 chi2, const char* freq_range, UINT8 nMode,INT4 nk_max);
UNUSED static UINT8 SEOBNRv4HMROM_phase_sparse_grid_hybrid(gsl_vector* freq, gsl_vector* phase, gsl_vector* freq_lo, gsl_vector* freq_hi, REAL8 q, REAL8 chi1,REAL8 chi2,INT4 nk_max);
UNUSED static UINT8 SEOBNRv4HMROM_cmode_sparse_grid_hybrid(gsl_vector* freq, gsl_vector* cmode_real, gsl_vector* cmode_imag, gsl_vector* freq_lo, gsl_vector* freq_hi, REAL8 q, REAL8 chi1,REAL8 chi2, UINT8 nMode,INT8 nk_max);
UNUSED static UINT8 SEOBNRv4HMROM_freq_phase_sparse_grid_hybrid(UNUSED gsl_vector** freq, UNUSED gsl_vector** phase, REAL8 q, REAL8 chi1,REAL8 chi2,UINT8 flag_patch,INT4 nk_max);
UNUSED static UINT8 SEOBNRv4HMROM_freq_cmode_sparse_grid_hybrid(gsl_vector** freq_cmode, gsl_vector** cmode_real, gsl_vector** cmode_imag, REAL8 q, REAL8 chi1,REAL8 chi2, UINT8 nMode, UINT8 flag_patch,INT8 nk_max);
UNUSED static UINT8 SEOBNRv4HMROM_approx_phi_lm(gsl_vector* freq_mode_lm, gsl_vector* phase_approx_lm, gsl_vector* freq_carrier_hyb, gsl_vector* phase_carrier_hyb, UINT4 nMode);
/********************* Definitions begin here ********************/


/** Setup SEOBNRv4HMROM model using data files installed in $LAL_DATA_PATH
 */
UNUSED static void SEOBNRv4HMROM_Init_LALDATA()
{
  // Loop over all the available modes and check if they have already been initialized
  for(int i = 0; i < NMODES; i++) {
    if (SEOBNRv4HMROM_IsSetup(i)) return;
  }

  // Expect ROM datafile in a directory listed in LAL_DATA_PATH,
#ifdef LAL_HDF5_ENABLED
#define datafile ROMDataHDF5
  char *path = XLALFileResolvePathLong(datafile, PKG_DATA_DIR);
  if (path==NULL){
    XLAL_ERROR_VOID(XLAL_EIO, "Unable to resolve data file %s in $LAL_DATA_PATH\n", datafile);
  }  
  char *dir = dirname(path);
  
  UINT4 ret = XLAL_SUCCESS;
  for(int i = 0; i < NMODES; i++) {
    ret = SEOBNRv4HMROM_Init(dir, i);
    if(ret != XLAL_SUCCESS)
      XLAL_ERROR_VOID(XLAL_FAILURE, "Unable to find SEOBNRv4HMROM data \
      files in $LAL_DATA_PATH for the mode = %d\n", i);
  }
  XLALFree(path);
#else
  XLAL_ERROR_VOID(XLAL_EFAILED, "SEOBNRv4HMROM requires HDF5 support which is not enabled\n");
#endif
}

/** Helper function to check if the SEOBNRv4HMROM model has been initialised */
static bool SEOBNRv4HMROM_IsSetup(UINT4 index_mode) {
  if(__lalsim_SEOBNRv4HMROMDS_data[index_mode].setup){
    return true;
  }
  else{
    return false;
  }
}

// Setup B-spline basis functions for given points
static void SplineData_Init(
  SplineData **splinedata,
  int ncx,                // Number of points in q  + 2
  int ncy,                // Number of points in chi1 + 2
  int ncz,                // Number of points in chi2 + 2
  const double *qvec,   // B-spline knots in q
  const double *chi1vec,  // B-spline knots in chi1
  const double *chi2vec   // B-spline knots in chi2
)
{
  if(!splinedata) exit(1);
  if(*splinedata) SplineData_Destroy(*splinedata);

  (*splinedata)=XLALCalloc(1,sizeof(SplineData));

  // Set up B-spline basis for desired knots
  const size_t nbreak_x = ncx-2;  // must have nbreak = n-2 for cubic splines
  const size_t nbreak_y = ncy-2;  // must have nbreak = n-2 for cubic splines
  const size_t nbreak_z = ncz-2;  // must have nbreak = n-2 for cubic splines

  // Allocate a cubic bspline workspace (k = 4)
  gsl_bspline_workspace *bwx = gsl_bspline_alloc(4, nbreak_x);
  gsl_bspline_workspace *bwy = gsl_bspline_alloc(4, nbreak_y);
  gsl_bspline_workspace *bwz = gsl_bspline_alloc(4, nbreak_z);

  // Set breakpoints (and thus knots by hand)
  gsl_vector *breakpts_x = gsl_vector_alloc(nbreak_x);
  gsl_vector *breakpts_y = gsl_vector_alloc(nbreak_y);
  gsl_vector *breakpts_z = gsl_vector_alloc(nbreak_z);
  for (UINT4 i=0; i<nbreak_x; i++)
    gsl_vector_set(breakpts_x, i, qvec[i]);
  for (UINT4 j=0; j<nbreak_y; j++)
    gsl_vector_set(breakpts_y, j, chi1vec[j]);
  for (UINT4 k=0; k<nbreak_z; k++)
    gsl_vector_set(breakpts_z, k, chi2vec[k]);

  gsl_bspline_knots(breakpts_x, bwx);
  gsl_bspline_knots(breakpts_y, bwy);
  gsl_bspline_knots(breakpts_z, bwz);

  gsl_vector_free(breakpts_x);
  gsl_vector_free(breakpts_y);
  gsl_vector_free(breakpts_z);

  (*splinedata)->bwx=bwx;
  (*splinedata)->bwy=bwy;
  (*splinedata)->bwz=bwz;
}

/* Destroy B-spline basis functions for given points */
static void SplineData_Destroy(SplineData *splinedata)
{
  if(!splinedata) return;
  if(splinedata->bwx) gsl_bspline_free(splinedata->bwx);
  if(splinedata->bwy) gsl_bspline_free(splinedata->bwy);
  if(splinedata->bwz) gsl_bspline_free(splinedata->bwz);
  XLALFree(splinedata);
}

// Interpolate projection coefficients for either Re(c-mode), Im(c-mode) or orbital phase over the parameter space (q, chi).
// The multi-dimensional interpolation is carried out via a tensor product decomposition.
static int TP_Spline_interpolation_3d(
  REAL8 q,                  // Input: eta-value for which projection coefficients should be evaluated
  REAL8 chi1,               // Input: chi1-value for which projection coefficients should be evaluated
  REAL8 chi2,               // Input: chi2-value for which projection coefficients should be evaluated
  gsl_vector *cvec,         // Input: data for spline coefficients
  int nk,                   // number of SVD-modes == number of basis functions
  int nk_max,               // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
  int ncx,                  // Number of points in eta  + 2
  int ncy,                  // Number of points in chi1 + 2
  int ncz,                  // Number of points in chi2 + 2
  const double *qvec,       // B-spline knots in eta
  const double *chi1vec,    // B-spline knots in chi1
  const double *chi2vec,    // B-spline knots in chi2
  gsl_vector *c_out        // Output: interpolated projection coefficients
  ) {
  if (nk_max != -1) {
    if (nk_max > nk)
      XLAL_ERROR(XLAL_EDOM, "Truncation parameter nk_max %d must be smaller or equal to nk %d\n", nk_max, nk);
    else { // truncate SVD modes
      nk = nk_max;
    }
  }

  SplineData *splinedata=NULL;
  SplineData_Init(&splinedata, ncx, ncy, ncz, qvec, chi1vec, chi2vec);

  gsl_bspline_workspace *bwx=splinedata->bwx;
  gsl_bspline_workspace *bwy=splinedata->bwy;
  gsl_bspline_workspace *bwz=splinedata->bwz;

  int N = ncx*ncy*ncz;  // Size of the data matrix for one SVD-mode
  // Evaluate the TP spline for all SVD modes - amplitude
  for (int k=0; k<nk; k++) { // For each SVD mode
    gsl_vector v = gsl_vector_subvector(cvec, k*N, N).vector; // Pick out the coefficient matrix corresponding to the k-th SVD mode.
    REAL8 csum = Interpolate_Coefficent_Tensor(&v, q, chi1, chi2, ncy, ncz, bwx, bwy, bwz);
    gsl_vector_set(c_out, k, csum);
  }
  SplineData_Destroy(splinedata);

  return(0);
}

/** Setup SEOBNRv4HMROM mode using data files installed in dir
 */
static int SEOBNRv4HMROM_Init(const char dir[],UINT4 index_mode) {
  if(__lalsim_SEOBNRv4HMROMDS_data[index_mode].setup) {
    XLALPrintError("Error: SEOBNRv4HMROM data was already set up!");
    XLAL_ERROR(XLAL_EFAILED);
  }
  SEOBNRROMdataDS_Init(&__lalsim_SEOBNRv4HMROMDS_data[index_mode], dir, index_mode);

  if(__lalsim_SEOBNRv4HMROMDS_data[index_mode].setup) {
    return(XLAL_SUCCESS);
  }
  else {
    return(XLAL_EFAILED);
  }
  return(XLAL_SUCCESS);
}

/* Set up a new ROM mode, using data contained in dir */
int SEOBNRROMdataDS_Init(
  UNUSED SEOBNRROMdataDS *romdata,
  UNUSED const char dir[],
  UNUSED UINT4 index_mode)
{
  int ret = XLAL_FAILURE;

  /* Create storage for structures */
  if(romdata->setup) {
    XLALPrintError("WARNING: You tried to setup the SEOBNRv4HMROM model that was already initialised. Ignoring\n");
    return (XLAL_FAILURE);
  }

#ifdef LAL_HDF5_ENABLED
  // First, check we got the correct version number
  size_t size = strlen(dir) + strlen(ROMDataHDF5) + 2;
  char *path = XLALMalloc(size);
  snprintf(path, size, "%s/%s", dir, ROMDataHDF5);
  LALH5File *file = XLALH5FileOpen(path, "r");

  XLALPrintInfo("ROM metadata\n============\n");
  PrintInfoStringAttribute(file, "Email");
  PrintInfoStringAttribute(file, "Description");
  ret = ROM_check_version_number(file, ROMDataHDF5_VERSION_MAJOR,
                                 ROMDataHDF5_VERSION_MINOR,
                                 ROMDataHDF5_VERSION_MICRO);

  ret |= SEOBNRROMdataDS_Init_submodel(&(romdata)->hqhs, dir, "hqhs",index_mode);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : submodel high q high spins loaded sucessfully.\n", __func__);

  ret |= SEOBNRROMdataDS_Init_submodel(&(romdata)->hqls, dir, "hqls",index_mode);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : submodel high q low spins loaded sucessfully.\n", __func__);

  ret |= SEOBNRROMdataDS_Init_submodel(&(romdata)->lqhs, dir, "lqhs",index_mode);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : submodel low q high spins loaded sucessfully.\n", __func__);

  ret |= SEOBNRROMdataDS_Init_submodel(&(romdata)->lqls, dir, "lqls",index_mode);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : submodel low q low spins loaded sucessfully.\n", __func__);

  ret |= SEOBNRROMdataDS_Init_submodel(&(romdata)->lowf, dir, "lowf",index_mode);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : submodel low freqs loaded sucessfully.\n", __func__);


  if(XLAL_SUCCESS==ret){
    romdata->setup=1;
  }
   else
     SEOBNRROMdataDS_Cleanup(romdata);
#else
  XLAL_ERROR(XLAL_EFAILED, "HDF5 support not enabled");
#endif

  XLALFree(path);
  XLALH5FileClose(file);
  ret = XLAL_SUCCESS;

  return (ret);
}

/* Deallocate contents of the given SEOBNRROMdataDS structure */
static void SEOBNRROMdataDS_Cleanup(SEOBNRROMdataDS *romdata) {
  SEOBNRROMdataDS_Cleanup_submodel((romdata)->hqhs);
  XLALFree((romdata)->hqhs);
  (romdata)->hqhs = NULL;
  SEOBNRROMdataDS_Cleanup_submodel((romdata)->hqls);
  XLALFree((romdata)->hqls);
  (romdata)->hqls = NULL;
  SEOBNRROMdataDS_Cleanup_submodel((romdata)->lqhs);
  XLALFree((romdata)->lqhs);
  (romdata)->lqhs = NULL;
  SEOBNRROMdataDS_Cleanup_submodel((romdata)->lqls);
  XLALFree((romdata)->lqls);
  (romdata)->lqls = NULL;
  SEOBNRROMdataDS_Cleanup_submodel((romdata)->lowf);
  XLALFree((romdata)->lowf);
  (romdata)->lowf = NULL;
  romdata->setup=0;
}

/* Deallocate contents of the given SEOBNRROMdataDS_submodel structure */
static void SEOBNRROMdataDS_Cleanup_submodel(SEOBNRROMdataDS_submodel *submodel) {
  if(submodel->cvec_real) gsl_vector_free(submodel->cvec_real);
  if(submodel->cvec_imag) gsl_vector_free(submodel->cvec_imag);
  if(submodel->cvec_phase) gsl_vector_free(submodel->cvec_phase);
  if(submodel->Breal) gsl_matrix_free(submodel->Breal);
  if(submodel->Bimag) gsl_matrix_free(submodel->Bimag);
  if(submodel->Bphase) gsl_matrix_free(submodel->Bphase);
  if(submodel->gCMode)   gsl_vector_free(submodel->gCMode);
  if(submodel->gPhase) gsl_vector_free(submodel->gPhase);
  if(submodel->qvec)  gsl_vector_free(submodel->qvec);
  if(submodel->chi1vec) gsl_vector_free(submodel->chi1vec);
  if(submodel->chi2vec) gsl_vector_free(submodel->chi2vec);
}

/* Set up a new ROM submodel, using data contained in dir */
UNUSED static int SEOBNRROMdataDS_Init_submodel(
  SEOBNRROMdataDS_submodel **submodel,
  UNUSED const char dir[],
  UNUSED const char grp_name[],
  UINT4 index_mode
) {
  int ret = XLAL_FAILURE;
  //printf("Looking at mode %s = \n",mode_array[index_mode]);
  if(!submodel) exit(1);
  /* Create storage for submodel structures */
  if (!*submodel)
    *submodel = XLALCalloc(1,sizeof(SEOBNRROMdataDS_submodel));
  else
    SEOBNRROMdataDS_Cleanup_submodel(*submodel);

#ifdef LAL_HDF5_ENABLED
  size_t size = strlen(dir) + strlen(ROMDataHDF5) + 2;
  char *path = XLALMalloc(size);
  snprintf(path, size, "%s/%s", dir, ROMDataHDF5);

  LALH5File *file = XLALH5FileOpen(path, "r");
  LALH5File *sub = XLALH5GroupOpen(file, grp_name);

  // Read ROM coefficients

  //// c-modes coefficients
  char* path_to_dataset = concatenate_strings(3,"CF_modes/",mode_array[index_mode],"/coeff_re_flattened");
  ReadHDF5RealVectorDataset(sub, path_to_dataset, & (*submodel)->cvec_real);
  free(path_to_dataset);
  path_to_dataset = concatenate_strings(3,"CF_modes/",mode_array[index_mode],"/coeff_im_flattened");
  ReadHDF5RealVectorDataset(sub, path_to_dataset, & (*submodel)->cvec_imag);
  free(path_to_dataset);
  //// orbital phase coefficients
  //// They are used only in the 22 mode
  if(index_mode == 0){
    ReadHDF5RealVectorDataset(sub, "phase_carrier/coeff_flattened", & (*submodel)->cvec_phase);
  }
  

  // Read ROM basis functions

  //// c-modes basis
  path_to_dataset = concatenate_strings(3,"CF_modes/",mode_array[index_mode],"/basis_re");
  ReadHDF5RealMatrixDataset(sub, path_to_dataset, & (*submodel)->Breal);
  free(path_to_dataset);
  path_to_dataset = concatenate_strings(3,"CF_modes/",mode_array[index_mode],"/basis_im");
  ReadHDF5RealMatrixDataset(sub, path_to_dataset, & (*submodel)->Bimag);
  free(path_to_dataset);
  //// orbital phase basis
  //// Used only in the 22 mode
  if(index_mode == 0){
    ReadHDF5RealMatrixDataset(sub, "phase_carrier/basis", & (*submodel)->Bphase);
  }
  // Read sparse frequency points

  //// c-modes grid
  path_to_dataset = concatenate_strings(3,"CF_modes/",mode_array[index_mode],"/MF_grid");
  ReadHDF5RealVectorDataset(sub, path_to_dataset, & (*submodel)->gCMode);
  free(path_to_dataset);
  //// orbital phase grid
  //// Used only in the 22 mode
  if(index_mode == 0){
    ReadHDF5RealVectorDataset(sub, "phase_carrier/MF_grid", & (*submodel)->gPhase);
  }
  // Read parameter space nodes
  ReadHDF5RealVectorDataset(sub, "qvec", & (*submodel)->qvec);
  ReadHDF5RealVectorDataset(sub, "chi1vec", & (*submodel)->chi1vec);
  ReadHDF5RealVectorDataset(sub, "chi2vec", & (*submodel)->chi2vec);


  // Initialize other members
  (*submodel)->nk_cmode = (*submodel)->gCMode->size;
  //// Used only in the 22 mode
  if(index_mode == 0){
    (*submodel)->nk_phase = (*submodel)->gPhase->size; 
  }  
  (*submodel)->ncx = (*submodel)->qvec->size + 2;
  (*submodel)->ncy = (*submodel)->chi1vec->size + 2;
  (*submodel)->ncz = (*submodel)->chi2vec->size + 2;

  // Domain of definition of submodel
  (*submodel)->q_bounds[0] = gsl_vector_get((*submodel)->qvec, 0);
  (*submodel)->q_bounds[1] = gsl_vector_get((*submodel)->qvec, (*submodel)->qvec->size - 1);
  (*submodel)->chi1_bounds[0] = gsl_vector_get((*submodel)->chi1vec, 0);
  (*submodel)->chi1_bounds[1] = gsl_vector_get((*submodel)->chi1vec, (*submodel)->chi1vec->size - 1);
  (*submodel)->chi2_bounds[0] = gsl_vector_get((*submodel)->chi2vec, 0);
  (*submodel)->chi2_bounds[1] = gsl_vector_get((*submodel)->chi2vec, (*submodel)->chi2vec->size - 1);

  XLALFree(path);
  XLALH5FileClose(file);
  XLALH5FileClose(sub);
  ret = XLAL_SUCCESS;
#else
  XLAL_ERROR(XLAL_EFAILED, "HDF5 support not enabled");
#endif

  return ret;
}

/* Create structure for internal use to store coefficients for orbital phase and real/imaginary part of co-orbital modes */
static void SEOBNRROMdataDS_coeff_Init(SEOBNRROMdataDS_coeff **romdatacoeff, int nk_cmode, int nk_phase) {
  if(!romdatacoeff) exit(1);
  /* Create storage for structures */
  if(!*romdatacoeff)
    *romdatacoeff=XLALCalloc(1,sizeof(SEOBNRROMdataDS_coeff));
  else
    SEOBNRROMdataDS_coeff_Cleanup(*romdatacoeff);

  (*romdatacoeff)->c_real = gsl_vector_alloc(nk_cmode);
  (*romdatacoeff)->c_imag = gsl_vector_alloc(nk_cmode);
  (*romdatacoeff)->c_phase = gsl_vector_alloc(nk_phase);
}

/* Deallocate contents of the given SEOBNRROMdataDS_coeff structure */
static void SEOBNRROMdataDS_coeff_Cleanup(SEOBNRROMdataDS_coeff *romdatacoeff) {
  if(romdatacoeff->c_real) gsl_vector_free(romdatacoeff->c_real);
  if(romdatacoeff->c_imag) gsl_vector_free(romdatacoeff->c_imag);
  if(romdatacoeff->c_phase) gsl_vector_free(romdatacoeff->c_phase);
  XLALFree(romdatacoeff);
}
/* Return the closest higher power of 2  */
// Note: NextPow(2^k) = 2^k for integer values k.
static size_t NextPow2(const size_t n) {
  return 1 << (size_t) ceil(log2(n));
}
/* This function returns the orbital phase in the sparse grid hybridized between LF and HF ROM */
UNUSED static UINT8 SEOBNRv4HMROM_phase_sparse_grid_hybrid(gsl_vector* freq, gsl_vector* phase, gsl_vector* freq_lo, gsl_vector* freq_hi, REAL8 q, REAL8 chi1,REAL8 chi2,INT4 nk_max){
  
  gsl_vector *phase_f_lo = gsl_vector_alloc(freq_lo->size);
  gsl_vector *phase_f_hi = gsl_vector_alloc(freq_hi->size);

  /* Get LF and HF orbital phase */
  UNUSED UINT8 retcode = SEOBNRv4HMROM_phase_sparse_grid(phase_f_lo,q,chi1,chi2, "LF",nk_max);
  retcode = SEOBNRv4HMROM_phase_sparse_grid(phase_f_hi,q,chi1,chi2, "HF",nk_max);

  /* Align LF and HF orbital phase */
  retcode = align_wfs_window(freq_lo, freq_hi, phase_f_lo, phase_f_hi, f_hyb_ini, f_hyb_end);

  /* Blend LF and HF orbital phase together */
  retcode = blend_functions(freq,phase,freq_lo,phase_f_lo,freq_hi,phase_f_hi,f_hyb_ini, f_hyb_end);

  /* Cleanup */
  gsl_vector_free(phase_f_lo);
  gsl_vector_free(phase_f_hi);

  return XLAL_SUCCESS;

}

/* This function returns the frequency and orbital phase in a sparse grid. The frequency and orbital phase are hybridized between LF and HF ROM. */
UNUSED static UINT8 SEOBNRv4HMROM_freq_phase_sparse_grid_hybrid(UNUSED gsl_vector** freq, UNUSED gsl_vector** phase, REAL8 q, REAL8 chi1,REAL8 chi2,UINT8 flag_patch,INT4 nk_max){
  
  /* We point to the data for computing the orbital phase. */
  /* These data are stored in the structures for all modes */
  /* We use the 22 mode (index 0) structure because the 22 mode will always be present */
  SEOBNRROMdataDS *romdata=&__lalsim_SEOBNRv4HMROMDS_data[0];
  if (!SEOBNRv4HMROM_IsSetup(0)) {
    XLAL_ERROR(XLAL_EFAILED,
               "Error setting up SEOBNRv4ROM data - check your $LAL_DATA_PATH\n");
  }

  UNUSED UINT4 retcode=0;

  //We always need to glue two submodels together for this ROM
  SEOBNRROMdataDS_submodel *submodel_hi; // high frequency ROM
  SEOBNRROMdataDS_submodel *submodel_lo; // low frequency ROM

  /* Select ROM patches */
  //There is only one LF patch
  submodel_lo = romdata->lowf;

  //Select the HF patch
  if(flag_patch == 0){
    submodel_hi = romdata->hqls;
  }
  else if(flag_patch == 1){
    submodel_hi = romdata->hqhs;
  }
  else if(flag_patch == 2){
    submodel_hi = romdata->lqls;
  }
  else if(flag_patch == 3){
    submodel_hi = romdata->lqhs;
  }
  else{
    XLAL_ERROR( XLAL_EDOM );
    XLALPrintError("XLAL Error - %s: SEOBNRv4HM not currently available in this region!\n",
                  __func__);
  }

  /* Compute hybridized orbital phase */

  INT8 i_max_LF = 0;
  INT8 i_min_HF = 0;

  // Get sparse frequency array for the LF orbital phase
  gsl_vector* freq_lo = gsl_vector_alloc(submodel_lo->nk_phase);
  for(unsigned int i = 0; i < freq_lo->size; i++){
    freq_lo->data[i] = submodel_lo->gPhase->data[i];
  }
  gsl_vector* freq_hi = gsl_vector_alloc(submodel_hi->nk_phase);
  
  //Compute prefactor for rescaling the frequency. 
  //The HF ROM is build using freq/RD_freq, we need to compute RD_freq and undo the scaling

  REAL8 omegaQNM = Get_omegaQNM_SEOBNRv4(q,chi1,chi2,2,2);
  REAL8 inv_scaling = 1. / (2*LAL_PI/omegaQNM);
  for(unsigned int i = 0; i < freq_hi->size; i++){
    freq_hi->data[i] = inv_scaling*submodel_hi->gPhase->data[i];
  }

  //Initialize arrays for hybridized phase
  // Compute length of hybridized frequency array
  retcode = compute_i_max_LF_i_min_HF(&i_max_LF, &i_min_HF,freq_lo,freq_hi,f_hyb_ini);
  *freq = gsl_vector_alloc(i_max_LF + freq_hi->size - i_min_HF + 1);

  //It's important to initialize this function to 0
  *phase = gsl_vector_calloc((*freq)->size);

  // Fill the hybridized frequency array
  for(unsigned int i = 0; i <= i_max_LF;i++){
      (*freq)->data[i] = freq_lo->data[i];
  };
  for(unsigned int i = i_min_HF; i < freq_hi->size;i++){
      (*freq)->data[i_max_LF+i-i_min_HF+1] = freq_hi->data[i];
  };

  /* Get the hybridized orbital phase */

  retcode = SEOBNRv4HMROM_phase_sparse_grid_hybrid(*freq,*phase,freq_lo, freq_hi,q,chi1,chi2,nk_max);

  /* Cleanup */
  gsl_vector_free(freq_lo);
  gsl_vector_free(freq_hi);


  return XLAL_SUCCESS;

}

/* This function returns the frequency and approximate phase of nMode in a sparse grid hybridized between LF and HF ROM. */
/* Here we linearly extend the phase outside the regime where we have data */
/* The derivation of this approximate phase is in https://git.ligo.org/waveforms/reviews/SEOBNRv4HM_ROM/blob/master/documents/notes.pdf */
/* In this document this phase is the argument of the exponential in Eq.(22) */
UNUSED static UINT8 SEOBNRv4HMROM_approx_phi_lm(gsl_vector* freq_mode_lm, gsl_vector* phase_approx_lm, gsl_vector* freq_carrier_hyb, gsl_vector* phase_carrier_hyb, UINT4 nMode){

  UINT8 modeM = lmModes[nMode][1];

  REAL8 f_max_carrier = freq_carrier_hyb->data[freq_carrier_hyb->size -1];
  REAL8 phase_carrier_at_f_max = phase_carrier_hyb->data[phase_carrier_hyb->size -1];

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, phase_carrier_hyb->size);

  gsl_spline_init (spline, freq_carrier_hyb->data, phase_carrier_hyb->data, phase_carrier_hyb->size);
  REAL8 der_phase_carrier_at_f_max = gsl_spline_eval_deriv(spline, f_max_carrier, acc);

  /* Constant shift between modes */
  REAL8 const_phase_shift = const_phaseshift_lm[nMode] + (1.-(REAL8)modeM)*LAL_PI/4.;

  for(unsigned int i = 0; i < freq_mode_lm->size; i++){
    if(freq_mode_lm->data[i]/(REAL8)modeM < f_max_carrier){
      /* Compute the approximate phase of lm mode from the carrier phase */
      phase_approx_lm->data[i] = (REAL8)modeM*gsl_spline_eval (spline, freq_mode_lm->data[i]/(REAL8)modeM, acc) + const_phase_shift;
    }
    else{
      /* Linear extrapolation for frequency at which we don't have data for the orbital phase */
      phase_approx_lm->data[i] = modeM*(phase_carrier_at_f_max + der_phase_carrier_at_f_max*(freq_mode_lm->data[i]/modeM-f_max_carrier)) + const_phase_shift;
    }
  }

  /* Cleanup */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  return XLAL_SUCCESS;
}

/* This function returns the Re/Im part of c-modes in a sparse grid. Re/Im of c-modes are hybridized between LF and HF ROM. */
/* The derivation of the c-modes is in https://git.ligo.org/waveforms/reviews/SEOBNRv4HM_ROM/blob/master/documents/notes.pdf Eq. (22)*/
UNUSED static UINT8 SEOBNRv4HMROM_cmode_sparse_grid_hybrid(gsl_vector* freq, gsl_vector* cmode_real, gsl_vector* cmode_imag, gsl_vector* freq_lo, gsl_vector* freq_hi, REAL8 q, REAL8 chi1,REAL8 chi2, UINT8 nMode,INT8 nk_max){

  gsl_vector* real_f_lo = gsl_vector_alloc(freq_lo->size);
  gsl_vector* imag_f_lo = gsl_vector_alloc(freq_lo->size);
  gsl_vector* real_f_hi = gsl_vector_alloc(freq_hi->size);
  gsl_vector* imag_f_hi = gsl_vector_alloc(freq_hi->size);

  /* Generate LF and HF cmode */
  UNUSED UINT8 retcode = SEOBNRv4HMROM_cmode_sparse_grid(real_f_lo, imag_f_lo,q,chi1,chi2,"LF",nMode,nk_max);
  retcode = SEOBNRv4HMROM_cmode_sparse_grid(real_f_hi, imag_f_hi,q,chi1,chi2,"HF",nMode,nk_max);

  UINT8 modeM = lmModes[nMode][1];

  /* Blend LF and HF cmode together */
  retcode = blend_functions(freq,cmode_real,freq_lo,real_f_lo,freq_hi,real_f_hi,f_hyb_ini*(REAL8)modeM,f_hyb_end*(REAL8)modeM);
  retcode = blend_functions(freq,cmode_imag,freq_lo,imag_f_lo,freq_hi,imag_f_hi,f_hyb_ini*(REAL8)modeM,f_hyb_end*(REAL8)modeM);

  /* Cleanup */
  gsl_vector_free(real_f_lo);
  gsl_vector_free(imag_f_lo);
  gsl_vector_free(real_f_hi);
  gsl_vector_free(imag_f_hi);

  return XLAL_SUCCESS;

}

/* This function returns the frequency and Re/Im part of c-modes in a sparse grid. The frequency and Re/Im of c-modes are hybridized between LF and HF ROM. */
UNUSED static UINT8 SEOBNRv4HMROM_freq_cmode_sparse_grid_hybrid(gsl_vector** freq_cmode, gsl_vector** cmode_real, gsl_vector** cmode_imag, REAL8 q, REAL8 chi1,REAL8 chi2, UINT8 nMode, UINT8 flag_patch,INT8 nk_max){
  
  UINT8 modeL = lmModes[nMode][0];
  UINT8 modeM = lmModes[nMode][1];

  /* We point to the data for computing the cmode. */
  SEOBNRROMdataDS *romdata=&__lalsim_SEOBNRv4HMROMDS_data[nMode];
  if (!SEOBNRv4HMROM_IsSetup(nMode)) {
    XLAL_ERROR(XLAL_EFAILED,
               "Error setting up SEOBNRv4ROM data - check your $LAL_DATA_PATH\n");
  }

  UNUSED UINT4 retcode=0;

  /* We always need to glue two submodels together for this ROM */
  SEOBNRROMdataDS_submodel *submodel_hi; // high frequency ROM
  SEOBNRROMdataDS_submodel *submodel_lo; // low frequency ROM

    /* Select ROM patches */
  //There is only one LF patch
  submodel_lo = romdata->lowf;

  //Select the HF patch
  if(flag_patch == 0){
    submodel_hi = romdata->hqls;
  }
  else if(flag_patch == 1){
    submodel_hi = romdata->hqhs;
  }
  else if(flag_patch == 2){
    submodel_hi = romdata->lqls;
  }
  else if(flag_patch == 3){
    submodel_hi = romdata->lqhs;
  }
  else{
    XLAL_ERROR( XLAL_EDOM );
    XLALPrintError("XLAL Error - %s: SEOBNRv4HM not currently available in this region!\n",
                  __func__);
  }

  /* Compute hybridized Re/Im part c-mode */

  INT8 i_max_LF = 0;
  INT8 i_min_HF = 0;

  // Get sparse frequency array for the LF orbital phase
  gsl_vector* freq_lo = gsl_vector_alloc(submodel_lo->nk_cmode);
  gsl_vector* freq_hi = gsl_vector_alloc(submodel_hi->nk_cmode);

  for(unsigned int i = 0; i < freq_lo->size; i++){
    freq_lo->data[i] = submodel_lo->gCMode->data[i];
  }
  
  //Compute prefactor for rescaling the frequency. 
  //The HF ROM is build using freq/RD_freq, we need to compute RD_freq and undo the scaling
  REAL8 omegaQNM = Get_omegaQNM_SEOBNRv4(q,chi1,chi2,modeL,modeM);
  REAL8 inv_scaling = 1. / (2*LAL_PI/omegaQNM);
  for(unsigned int i = 0; i < freq_hi->size; i++){
    freq_hi->data[i] = inv_scaling*submodel_hi->gCMode->data[i];
  }

  //Initialize arrays for Re/Im part c-mode
  // Compute length of hybridized frequency array
  retcode = compute_i_max_LF_i_min_HF(&i_max_LF, &i_min_HF,freq_lo,freq_hi,f_hyb_ini*(REAL8)modeM);

  *freq_cmode = gsl_vector_alloc(i_max_LF + freq_hi->size - i_min_HF + 1);
  // It's important to initialize this function to 0

  for(unsigned int i = 0; i <= i_max_LF;i++){
      (*freq_cmode)->data[i] = freq_lo->data[i];
  };
  for(unsigned int i = i_min_HF; i < freq_hi->size;i++){
      (*freq_cmode)->data[i_max_LF+i-i_min_HF+1] = freq_hi->data[i];
  };

  *cmode_real = gsl_vector_calloc((*freq_cmode)->size);
  *cmode_imag = gsl_vector_calloc((*freq_cmode)->size);


  // Generate the hybridized Re/Im part of c-mode

  retcode = SEOBNRv4HMROM_cmode_sparse_grid_hybrid(*freq_cmode,*cmode_real,*cmode_imag,freq_lo, freq_hi,q,chi1,chi2,nMode,nk_max);


  /* Cleanup */
  gsl_vector_free(freq_lo);
  gsl_vector_free(freq_hi);

  return XLAL_SUCCESS;
}



/* This function returns the orbital phase in a sparse grid */
UNUSED static UINT8 SEOBNRv4HMROM_phase_sparse_grid(gsl_vector* phase_f, REAL8 q, REAL8 chi1,REAL8 chi2, const char* freq_range,INT4 nk_max){
  
  /* Loading the data for the 22 mode, the informations about orbital phase are stored here */
  SEOBNRROMdataDS *romdata=&__lalsim_SEOBNRv4HMROMDS_data[0];
  if (!SEOBNRv4HMROM_IsSetup(0)) {
    XLAL_ERROR(XLAL_EFAILED,
               "Error setting up SEOBNRv4ROM data - check your $LAL_DATA_PATH\n");
  }
  /* Check whether you are asking for the LF or the HF ROM */
  SEOBNRROMdataDS_submodel *submodel;
  if(strcmp("LF",freq_range) == 0){
    // In the case of LF there is only one patch
    submodel = romdata->lowf;
  }
  else{
    /* In the case of HF select the appropriate patch */
    UINT8 flag_patch = SEOBNRv4HMROM_Select_HF_patch(q,chi1);
    /* Based on the patch select the correct data */
    if(flag_patch == 0){
      submodel = romdata->hqls;
    }
    else if(flag_patch == 1){
      submodel = romdata->hqhs;
    }
    else if(flag_patch == 2){
      submodel = romdata->lqls;
    }
    else if(flag_patch == 3){
      submodel = romdata->lqhs;
    }
    else{
      XLAL_ERROR( XLAL_EDOM );
      XLALPrintError("XLAL Error - %s: SEOBNRv4HM not currently available in this region!\n",
                   __func__);
    }
  }
  /* Create and initialize the structure for interpolating the orbital phase */
  SEOBNRROMdataDS_coeff *romdata_coeff=NULL;
  SEOBNRROMdataDS_coeff_Init(&romdata_coeff, submodel->nk_cmode, submodel->nk_phase);

  /* Interpolating projection coefficients for orbital phase */
  UINT4 retcode=TP_Spline_interpolation_3d(
    q,                            // Input: q-value for which projection coefficients should be evaluated
    chi1,                         // Input: chi1-value for which projection coefficients should be evaluated
    chi2,                         // Input: chi2-value for which projection coefficients should be evaluated
    submodel->cvec_phase,         // Input: data for spline coefficients for amplitude
    submodel->nk_phase,           // number of SVD-modes == number of basis functions for amplitude
    nk_max,                       // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
    submodel->ncx,             // Number of points in q  + 2
    submodel->ncy,             // Number of points in chi1 + 2
    submodel->ncz,             // Number of points in chi2 + 2
    gsl_vector_const_ptr(submodel->qvec, 0),          // B-spline knots in q
    gsl_vector_const_ptr(submodel->chi1vec, 0),        // B-spline knots in chi1
    gsl_vector_const_ptr(submodel->chi2vec, 0),        // B-spline knots in chi2
    romdata_coeff->c_phase      // Output: interpolated projection coefficients for orbital phase
  );

  if(retcode!=0) {
    SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff);
    XLAL_ERROR(retcode);
  }

  /* Project interpolation coefficients onto the basis */
  /* phase_pts = B_A^T . c_A */
  gsl_blas_dgemv(CblasTrans, 1.0, submodel->Bphase, romdata_coeff->c_phase, 0.0, phase_f);

  /* Cleanup */
  SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff);
  
  return XLAL_SUCCESS;
}

/* This function returns the Re/Im part of c-modes in a sparse grid. */
UNUSED static UINT8 SEOBNRv4HMROM_cmode_sparse_grid(gsl_vector* creal_f, gsl_vector* cimag_f, REAL8 q, REAL8 chi1,REAL8 chi2, const char* freq_range,UINT8 nMode,INT4 nk_max){
  
  /* Loading the data for the lm mode */
  SEOBNRROMdataDS *romdata=&__lalsim_SEOBNRv4HMROMDS_data[nMode];
  if (!SEOBNRv4HMROM_IsSetup(nMode)) {
    XLAL_ERROR(XLAL_EFAILED,
               "Error setting up SEOBNRv4ROM data - check your $LAL_DATA_PATH\n");
  }
  /* Check whether you are asking for the LF or the HF ROM */
  SEOBNRROMdataDS_submodel *submodel;
  if(strcmp("LF",freq_range) == 0){
    /* In the case of LF there is only one patch */
    submodel = romdata->lowf;
  }
  else{
    /* In the case of HF select the appropriate patch */
    UINT8 flag_patch = SEOBNRv4HMROM_Select_HF_patch(q,chi1);
    /* Based on the patch select the correct data */
    if(flag_patch == 0){
      submodel = romdata->hqls;
    }
    else if(flag_patch == 1){
      submodel = romdata->hqhs;
    }
    else if(flag_patch == 2){
      submodel = romdata->lqls;
    }
    else if(flag_patch == 3){
      submodel = romdata->lqhs;
    }
    else{
      XLAL_ERROR( XLAL_EDOM );
      XLALPrintError("XLAL Error - %s: SEOBNRv4HM not currently available in this region!\n",
                   __func__);
    }
  }
  /* Create and initialize the structure for interpolating the coorbital mode */
  SEOBNRROMdataDS_coeff *romdata_coeff=NULL;
  SEOBNRROMdataDS_coeff_Init(&romdata_coeff, submodel->nk_cmode, submodel->nk_cmode);

  /* Interpolating projection coefficients for Re(cmode) */
  UINT4 retcode=TP_Spline_interpolation_3d(
    q,                            // Input: eta-value for which projection coefficients should be evaluated
    chi1,                         // Input: chi1-value for which projection coefficients should be evaluated
    chi2,                         // Input: chi2-value for which projection coefficients should be evaluated
    submodel->cvec_real,         // Input: data for spline coefficients for amplitude
    submodel->nk_cmode,           // number of SVD-modes == number of basis functions for amplitude
    nk_max,                       // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
    submodel->ncx,             // Number of points in q  + 2
    submodel->ncy,             // Number of points in chi1 + 2
    submodel->ncz,             // Number of points in chi2 + 2
    gsl_vector_const_ptr(submodel->qvec, 0),          // B-spline knots in eta
    gsl_vector_const_ptr(submodel->chi1vec, 0),        // B-spline knots in chi1
    gsl_vector_const_ptr(submodel->chi2vec, 0),        // B-spline knots in chi2
    romdata_coeff->c_real      // Output: interpolated projection coefficients for Re(c-mode)
  );

  /* Interpolating projection coefficients for Im(cmode) */
  retcode=TP_Spline_interpolation_3d(
    q,                            // Input: eta-value for which projection coefficients should be evaluated
    chi1,                         // Input: chi1-value for which projection coefficients should be evaluated
    chi2,                         // Input: chi2-value for which projection coefficients should be evaluated
    submodel->cvec_imag,         // Input: data for spline coefficients for amplitude
    submodel->nk_cmode,           // number of SVD-modes == number of basis functions for amplitude
    nk_max,                       // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
    submodel->ncx,             // Number of points in q  + 2
    submodel->ncy,             // Number of points in chi1 + 2
    submodel->ncz,             // Number of points in chi2 + 2
    gsl_vector_const_ptr(submodel->qvec, 0),          // B-spline knots in eta
    gsl_vector_const_ptr(submodel->chi1vec, 0),        // B-spline knots in chi1
    gsl_vector_const_ptr(submodel->chi2vec, 0),        // B-spline knots in chi2
    romdata_coeff->c_imag      // Output: interpolated projection coefficients for Im(c-mode)
  );

  if(retcode!=0) {
    SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff);
    XLAL_ERROR(retcode);
  }

  /* Project interpolation coefficients onto the basis */
  /* creal_pts = B_A^T . c_A */
  /* cimag_pts = B_A^T . c_A */
  gsl_blas_dgemv(CblasTrans, 1.0, submodel->Breal, romdata_coeff->c_real, 0.0, creal_f);
  gsl_blas_dgemv(CblasTrans, 1.0, submodel->Bimag, romdata_coeff->c_imag, 0.0, cimag_f);

  /* Cleanup */
  SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff);

  return XLAL_SUCCESS;
}


/* Select high frequency ROM submodel */
UNUSED static UINT8 SEOBNRv4HMROM_Select_HF_patch(REAL8 q, REAL8 chi1){
  if ((q> 3.) && (chi1<=0.8)){
    return 0;
  }
  else if ((q> 3.) && (chi1>0.8)){
    return 1;
  }
  else if ((q<= 3.) && (chi1<=0.8)){
    return 2;
  }    
  else if ((q<= 3.) && (chi1>0.8)){
    return 3;
  }
  else{
    XLAL_ERROR( XLAL_EDOM );
    XLALPrintError("XLAL Error - %s: SEOBNRv4HM not currently available in this region!\n",
                   __func__);
  }
}

/**
* Core function for computing the ROM modes.
* It rebuilds the modes starting from "orbital phase" and co-orbital modes.
* This function returns SphHarmFrequencySeries
*/
UNUSED static int SEOBNRv4HMROMCoreModes(
  UNUSED SphHarmFrequencySeries **hlm_list, /**<< Spherical modes frequency series for the waveform */
  UNUSED REAL8 phiRef, /**<< orbital reference phase */
  UNUSED REAL8 fRef, /**< Reference frequency (Hz); 0 defaults to fLow */
  UNUSED REAL8 distance, /**< Distance of source (m) */
  UNUSED REAL8 Mtot_sec, /**< Total mass in seconds **/
  UNUSED REAL8 q, /**< Mass ratio **/
  UNUSED REAL8 chi1, /**< Dimensionless aligned component spin 1 */
  UNUSED REAL8 chi2, /**< Dimensionless aligned component spin 2 */
  UNUSED const REAL8Sequence *freqs_in, /**< Frequency points at which to evaluate the waveform (Hz) */
  UNUSED REAL8 deltaF, /**< Frequency points at which to evaluate the waveform (Hz) */
  UNUSED INT4 nk_max, /**< truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1. We
  *nk_max == -1 is the default setting */
  UNUSED UINT4 nModes, /**<  Number of modes to generate */
  REAL8 sign_odd_modes /**<  Sign of the odd-m modes, used when swapping the two bodies */
  )
{

  /* Check output structure */
  if(!hlm_list)
    XLAL_ERROR(XLAL_EFAULT);

  /* 'Nudge' parameter values to allowed boundary values if close by */
  if (q < 1.0)     nudge(&q, 1.0, 1e-6);
  if (q > 50.0)     nudge(&q, 50.0, 1e-6);

  if ( chi1 < -1.0 || chi2 < -1.0 || chi1 > 1.0 || chi2 > 1.0) {
    XLALPrintError("XLAL Error - %s: chi1 or chi2 smaller than -1.0 or larger than 1.0!\n"
                   "SEOBNRv4HMROM is only available for spins in the range -1 <= a/M <= 1.0.\n",
                   __func__);
    XLAL_ERROR( XLAL_EDOM );
  }

  if (q<1.0 || q > 50.0) {
    XLALPrintError("XLAL Error - %s: q (%f) bigger than 50.0 or unphysical!\n"
                   "SEOBNRv4HMROM is only available for q in the range 1 <= q <= 50.\n",
                   __func__, q);
    XLAL_ERROR( XLAL_EDOM );
  } 
   

  /* Find frequency bounds */
  if (!freqs_in) XLAL_ERROR(XLAL_EFAULT);
  REAL8 fLow  = freqs_in->data[0];
  REAL8 fHigh = freqs_in->data[freqs_in->length - 1];

  UNUSED REAL8 fLow_geom = fLow * Mtot_sec;
  REAL8 fHigh_geom = fHigh * Mtot_sec;
  // UNUSED REAL8 fRef_geom = fRef * Mtot_sec;
  REAL8 deltaF_geom = deltaF * Mtot_sec;
  // The maximum available frequency is the one associated with the 55 mode
  REAL8 Mf_ROM_max = const_fmax_lm[4]*Get_omegaQNM_SEOBNRv4(q,chi1,chi2,5,5)/(2.*LAL_PI);

  // Enforce allowed geometric frequency range
  if (fHigh_geom == 0)
    fHigh_geom = Mf_ROM_max;
  if (fLow_geom < Mf_low_22)
    XLAL_ERROR(XLAL_EDOM, "Starting frequency Mflow=%g is smaller than lowest frequency in ROM Mf=%g.\n", fLow_geom, Mf_low_22);
  if (fHigh_geom < Mf_low_22)
    XLAL_ERROR(XLAL_EDOM, "End frequency %g is smaller than ROM starting frequency %g!\n", fHigh_geom, Mf_low_22);
  if (fHigh_geom <= fLow_geom)
    XLAL_ERROR(XLAL_EDOM, "End frequency %g is smaller than (or equal to) starting frequency %g!\n", fHigh_geom, fLow_geom);

  /* There are severals HF patches */
  /* Based on the patch select the correct data */
  UINT8 flag_patch = SEOBNRv4HMROM_Select_HF_patch(q,chi1);
  gsl_vector *freq_carrier_hyb = NULL;
  gsl_vector *phase_carrier_hyb = NULL;  
  
  /* Compute the orbital phase, this is the same for every mode */
  /* for this reason it's outside the for */

  UNUSED int retcode = SEOBNRv4HMROM_freq_phase_sparse_grid_hybrid(&freq_carrier_hyb, &phase_carrier_hyb, q, chi1,chi2,flag_patch,nk_max);
  if(retcode != XLAL_SUCCESS) XLAL_ERROR(retcode);

  /* Put the phase in the right convention */

  for(unsigned int i = 0; i < phase_carrier_hyb->size; i++){
    phase_carrier_hyb->data[i] = -phase_carrier_hyb->data[i]; 
  }

  // // With this variable we keep track of the phase shift to ensure phi_22(f_ref) = phiRef
  // REAL8 phase_change = 0.;

  /* Looping over cmodes */
  for(unsigned int nMode = 0; nMode < nModes; nMode++){

    gsl_vector *freq_cmode_hyb = NULL;
    gsl_vector *cmode_real_hyb = NULL;
    gsl_vector *cmode_imag_hyb = NULL;


    UNUSED UINT4 modeL = lmModes[nMode][0];
    UNUSED UINT4 modeM = lmModes[nMode][1];

    /* Compute Re and Im c-modes */
    retcode = SEOBNRv4HMROM_freq_cmode_sparse_grid_hybrid(&freq_cmode_hyb, &cmode_real_hyb, &cmode_imag_hyb,q,chi1,chi2,nMode,flag_patch,nk_max);

    gsl_vector *phase_approx_lm = gsl_vector_alloc(freq_cmode_hyb->size);
    /* Compute approximated phase from orbital phase */
    retcode = SEOBNRv4HMROM_approx_phi_lm(freq_cmode_hyb, phase_approx_lm,freq_carrier_hyb, phase_carrier_hyb,nMode);

    /* Compute the phase contribution coming from the cmode and unwrap it */
    COMPLEX16 complex_mode;
    gsl_vector *phase_cmode = gsl_vector_alloc(freq_cmode_hyb->size);
    for(unsigned int i = 0; i < freq_cmode_hyb->size;i++){
      complex_mode = cmode_real_hyb->data[i] + I*cmode_imag_hyb->data[i];
      phase_cmode->data[i] = carg(complex_mode);
    }

    /* Unwrap the phase contribution from the cmode */
    gsl_vector *unwrapped_phase_cmode = gsl_vector_alloc(freq_cmode_hyb->size);
    retcode = unwrap_phase(unwrapped_phase_cmode,phase_cmode);

    
    /* Reconstruct amplitude and phase */
    gsl_vector *reconstructed_phase = gsl_vector_alloc(freq_cmode_hyb->size);
    gsl_vector *reconstructed_amplitude = gsl_vector_alloc(freq_cmode_hyb->size);

    for(unsigned int i = 0; i < freq_cmode_hyb->size;i++){
        complex_mode = cmode_real_hyb->data[i] + I*cmode_imag_hyb->data[i];
        reconstructed_amplitude->data[i] = cabs(complex_mode);
        reconstructed_phase->data[i] = unwrapped_phase_cmode->data[i] - phase_approx_lm->data[i];
    }
  
    size_t npts = 0;
    LIGOTimeGPS tC = {0, 0};
    UINT4 offset = 0; // Index shift between freqs and the frequency series
    REAL8Sequence *freqs = NULL;
    // freqs contains uniform frequency grid with spacing deltaF; we start at frequency 0
    // I removed the if statement for the time being
    /* Set up output array with size closest power of 2 */
    npts = NextPow2(fHigh_geom / deltaF_geom) + 1;
    if (fHigh_geom < fHigh * Mtot_sec) /* Resize waveform if user wants f_max larger than cutoff frequency */
      npts = NextPow2(fHigh * Mtot_sec / deltaF_geom) + 1;

    XLALGPSAdd(&tC, -1. / deltaF);  /* coalesce at t=0 */
    COMPLEX16FrequencySeries *hlmtilde = XLALCreateCOMPLEX16FrequencySeries("hlmtilde: FD mode", &tC, 0.0, deltaF, &lalStrainUnit, npts);
    

    // Recreate freqs using only the lower and upper bounds
    // Use fLow, fHigh and deltaF rather than geometric frequencies for numerical accuracy
    double fHigh_temp = fHigh_geom / Mtot_sec;
    UINT4 iStart = (UINT4) ceil(fLow / deltaF);
    UINT4 iStop = (UINT4) ceil(fHigh_temp / deltaF);
    freqs = XLALCreateREAL8Sequence(iStop - iStart);
    if (!freqs) {
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    }
    for (UINT4 i=iStart; i<iStop; i++)
      freqs->data[i-iStart] = i*deltaF_geom;

    offset = iStart;
    

    if (!hlmtilde)	{
        XLALDestroyREAL8Sequence(freqs);
        gsl_vector_free(freq_cmode_hyb);
        gsl_vector_free(cmode_real_hyb);
        gsl_vector_free(cmode_imag_hyb);
        gsl_vector_free(phase_approx_lm);
        gsl_vector_free(phase_cmode);
        gsl_vector_free(unwrapped_phase_cmode);
        gsl_vector_free(freq_carrier_hyb);
        gsl_vector_free(phase_carrier_hyb);
        XLAL_ERROR(XLAL_EFUNC);
    }

    memset(hlmtilde->data->data, 0, npts * sizeof(COMPLEX16));
    COMPLEX16 *hlmdata=hlmtilde->data->data;

    REAL8 Mtot = Mtot_sec / LAL_MTSUN_SI;
    REAL8 amp0 = Mtot * Mtot_sec * LAL_MRSUN_SI / (distance); // Correct overall amplitude to undo mass-dependent scaling used in ROM

    // Assemble mode from amplitude and phase
    gsl_interp_accel *acc_amp = gsl_interp_accel_alloc ();
    gsl_spline *spline_amp = gsl_spline_alloc (gsl_interp_cspline, reconstructed_amplitude->size);
    gsl_spline_init (spline_amp, freq_cmode_hyb->data, reconstructed_amplitude->data, reconstructed_amplitude->size);

    // ROM was build with t = 0 at t_peak^22 -1000. To agree with LAL convention, we undo this shift to get t_peak^22 = 0
    REAL8 t_corr = 1000.;
    gsl_interp_accel *acc_phase = gsl_interp_accel_alloc ();
    gsl_spline *spline_phase = gsl_spline_alloc (gsl_interp_cspline, reconstructed_phase->size);
    gsl_spline_init (spline_phase, freq_cmode_hyb->data, reconstructed_phase->data, reconstructed_phase->size);    

    // REAL8 phiRef_internal = -phiRef;
    // // Since ROM internal function assume a different convention for the FFT, we have later to compute
    // // the conjugation of the modes. This means that the phiRef that has to be applied at this point
    // // correspond to -phiRef selected by the user
    // if((modeL==2) && (modeM==2)){
    //   // Evaluate reference phase for setting phiRef correctly
    //   phase_change = gsl_spline_eval(spline_phase, fRef_geom, acc_phase) - 2*phiRef_internal;
    // }

    // Maximum frequency at which we have data for the ROM
    REAL8 Mf_max_mode = const_fmax_lm[nMode]*Get_omegaQNM_SEOBNRv4(q,chi1,chi2,modeL,modeM)/(2.*LAL_PI);

    for (UINT4 i=0; i<freqs->length; i++) { // loop over frequency points in sequence
      REAL8 f = freqs->data[i];
      if (f > Mf_max_mode) continue; // We're beyond the highest allowed frequency; since freqs may not be ordered, we'll just skip the current frequency and leave zero in the buffer
      if (f <= Mf_low_22 * modeM/2.) continue; // We're above the lowest allowed frequency; since freqs may not be ordered, we'll just skip the current frequency and leave zero in the buffer
      int j = i + offset; // shift index for frequency series if needed
      REAL8 A = gsl_spline_eval(spline_amp, f, acc_amp);
      //REAL8 phase = gsl_spline_eval(spline_phase, f, acc_phase)-0.*phase_change*(REAL8)modeM/2.;
      REAL8 phase = gsl_spline_eval(spline_phase, f, acc_phase);
      hlmdata[j] = amp0*A * (cos(phase) + I*sin(phase));//cexp(I*phase);
      REAL8 phase_factor = -2.*LAL_PI*f*t_corr;
      COMPLEX16 t_factor = cos(phase_factor) + I*sin(phase_factor);
      hlmdata[j] *= t_factor;
      // We now return the (l,-m) mode that in the LAL convention has support for f > 0
      // We use the equation h(l,-m)(f) = (-1)^l h(l,m)*(-f) with f > 0
      hlmdata[j] = pow(-1.,modeL)*conj(hlmdata[j]);
      if(modeM%2 != 0){
        // This is changing the sign of the odd-m modes in the case m1 < m2, if m1>m2 sign_odd_modes = 1 and nothing changes
        hlmdata[j] = hlmdata[j]*sign_odd_modes;
      } 
    }
    /* Save the mode (l,-m) in the SphHarmFrequencySeries structure */
    *hlm_list = XLALSphHarmFrequencySeriesAddMode(*hlm_list, hlmtilde, modeL, -modeM);

    /* Cleanup inside of loop over modes */
    XLALDestroyREAL8Sequence(freqs);
    XLALDestroyCOMPLEX16FrequencySeries(hlmtilde);
    gsl_vector_free(freq_cmode_hyb);
    gsl_vector_free(cmode_real_hyb);
    gsl_vector_free(cmode_imag_hyb);
    gsl_vector_free(phase_approx_lm);
    gsl_vector_free(phase_cmode);
    gsl_vector_free(unwrapped_phase_cmode);
    gsl_vector_free(reconstructed_amplitude);
    gsl_vector_free(reconstructed_phase);
    gsl_spline_free(spline_amp);
    gsl_spline_free(spline_phase);
    gsl_interp_accel_free(acc_amp);
    gsl_interp_accel_free(acc_phase);
  }
  


  /* Cleanup outside of loop over modes */
  gsl_vector_free(freq_carrier_hyb);
  gsl_vector_free(phase_carrier_hyb);


  return(XLAL_SUCCESS);
}


/**
 * ModeArray is a structure which allows to select the modes to include
 * in the waveform.
 * This function will create a structure with the default modes for every model
 */
static UINT8 Setup_EOBROM__std_mode_array_structure(
    LALValue *ModeArray, UINT4 nModes)
{

  /* setup ModeArray */

      if (nModes == 5) {
        /* Adding all the modes of SEOBNRv4HM
        * i.e. [(2,2),(2,1),(3,3),(4,4),(5,5)]
        the relative -m modes are added automatically*/
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -1);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 3, -3);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 4, -4);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 5, -5);
      }
      else{
        /*All the other spin aligned model so far only have the 22 mode
        */
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
      }


    return XLAL_SUCCESS;
}

/**
 * ModeArray is a structure which allows to select the modes to include
 * in the waveform.
 * This function check if the selected modes are available
 */
static INT8 Check_EOBROM_mode_array_structure(
    LALValue *ModeArray, UINT4 nModes)
{
  INT4 flagTrue = 0;
  UINT4 modeL;
  INT4 modeM;
  /*Loop over all the possible modes
  *we only check -m modes, in the waveform both + and - m
  *modes will be included
  */
  for (UINT4 ELL = 2; ELL <= LAL_SIM_L_MAX_MODE_ARRAY; ELL++) {
    for (INT4 EMM = -ELL; EMM <= (INT4)ELL; EMM++) {
      if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ELL, EMM) == 1) {
            for (UINT4 k = 0; k < nModes; k++) {
              modeL  = lmModes[k][0];
              //lmModes includes the modes in EOB convention for which hlm(f) is non 0 for f>0
              //for this reason we need the - sign here
              modeM = -lmModes[k][1]; 
              if ((modeL == ELL)&&(modeM == EMM)) {
                flagTrue=1;
              }
              if ((modeL == ELL)&&(modeM == -EMM)) {
                flagTrue=2;
              }
            }
            /*For each active mode check if is available for the selected model
            */
            if (flagTrue == 0) {
              XLALPrintError ("Mode (%d,%d) is not available by the model SEOBNRv4HM_ROM\n", ELL,
                 EMM);
              return XLAL_FAILURE;
            }
            if (flagTrue == 2) {
              XLALPrintError ("Mode (%d,%d) is not available by the model SEOBNRv4HM_ROM.\n"
              "In this function you can only select (l,-|m|) modes that are directly modeled in the ROM.\n"
              "The (l,+|m|) mode will be automatically added to the waveform using symmetry arguments.\n", ELL,
                 EMM);
              return XLAL_FAILURE;
            }
            flagTrue = 0;
          }
      }
    }

    return XLAL_SUCCESS;
}

/**
 * This function combines the modes hlm frequencyseries with the sYlm to produce the
 * polarizations hplus, hcross.
 */
// NOTE: azimuthal angle of the observer entering the -2Ylm is pi/2-phi according to LAL conventions
static int SEOBROMComputehplushcrossFromhlm(
    COMPLEX16FrequencySeries
        *hplusFS, /**<< Output: frequency series for hplus, already created */
    COMPLEX16FrequencySeries
        *hcrossFS,   /**<< Output: frequency series for hplus, already created */
    LALValue *ModeArray, /**<< Input: ModeArray structure with the modes to include */
    SphHarmFrequencySeries
        *hlm,  /**<< Input: list with frequency series for each mode hlm */
    REAL8 inc,  /**<< Input: inclination */
    UNUSED REAL8 phi   /**<< Input: phase */
) {
  /* Loop over modes */
  SphHarmFrequencySeries *hlms_temp = hlm;
  while ( hlms_temp ) {
    if (XLALSimInspiralModeArrayIsModeActive(ModeArray, hlms_temp->l, hlms_temp->m) == 1){
      /*Here we check if the mode generated is in the ModeArray structure
      */
        XLALSimAddModeFD(hplusFS, hcrossFS, hlms_temp->mode, inc, LAL_PI/2. - phi, hlms_temp->l, hlms_temp->m, 1);
        //XLALSimAddModeFD(hplusFS, hcrossFS, hlms_temp->mode, inc, 0., hlms_temp->l, hlms_temp->m, 1);
      }
      hlms_temp = hlms_temp->next;

  }

  return XLAL_SUCCESS;
}

/**
 * Compute waveform in LAL format for the SEOBNRv4HM_ROM model.
 *
 * Returns the plus and cross polarizations as a complex frequency series with
 * equal spacing deltaF and contains zeros from zero frequency to the starting
 * frequency fLow and zeros beyond the cutoff frequency in the ringdown.
 */
int XLALSimIMRSEOBNRv4HMROM(
  UNUSED struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  UNUSED struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  UNUSED REAL8 phiRef,                                 /**< Phase at reference time */
  UNUSED REAL8 deltaF,                                 /**< Sampling frequency (Hz) */
  REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
  UNUSED REAL8 fHigh,                                  /**< End frequency */
  UNUSED REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  UNUSED REAL8 distance,                               /**< Distance of source (m) */
  UNUSED REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 chi1,                                   /**< Dimensionless aligned component spin 1 */
  REAL8 chi2,                                   /**< Dimensionless aligned component spin 2 */
  UNUSED INT4 nk_max,                                   /**< Truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1 */
  UNUSED UINT4 nModes,                                  /**< Number of modes to use. This should be 1 for SEOBNRv4_ROM and 5 for SEOBNRv4HM_ROM */
  LALDict *LALParams /**<< Dictionary of additional wf parameters, here is used to pass ModeArray */
)
{
  REAL8 sign_odd_modes = 1.;
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI) {
    // Swap m1 and m2
    REAL8 m1temp = m1SI;
    REAL8 chi1temp = chi1;
    m1SI = m2SI;
    chi1 = chi2;
    m2SI = m1temp;
    chi2 = chi1temp;
    sign_odd_modes = -1.; // When we swap m1 with m2 odd-m modes should get - sign because they depend on dm
  }

  /* Get masses in terms of solar mass */
  REAL8 mass1 = m1SI / LAL_MSUN_SI;
  REAL8 mass2 = m2SI / LAL_MSUN_SI;
  REAL8 Mtot = mass1+mass2;
  REAL8 q = mass1/mass2;
  UNUSED REAL8 Mtot_sec = Mtot * LAL_MTSUN_SI;       /* Total mass in seconds */

  /* Select modes to include in the waveform */
  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(LALParams);
  //ModeArray includes the modes chosen by the user
  if (ModeArray == NULL) {
    //If the user doesn't choose the modes, use the standard modes
    ModeArray = XLALSimInspiralCreateModeArray();
    Setup_EOBROM__std_mode_array_structure(ModeArray, nModes);
  }
  //Check that the modes chosen are available for the model
  if(Check_EOBROM_mode_array_structure(ModeArray, nModes) == XLAL_FAILURE){
    XLALPrintError ("Not available mode chosen.\n");
    XLAL_ERROR(XLAL_EFUNC);
  }



  // if(fRef==0.0)
  //   fRef=fLow;

  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequency we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = fLow;
  freqs->data[1] = fHigh;  

  /* Load ROM data if not loaded already */
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&SEOBNRv4HMROM_is_initialized, SEOBNRv4HMROM_Init_LALDATA);
#else
  SEOBNRv4HMROM_Init_LALDATA();
#endif
  SphHarmFrequencySeries *hlm = NULL;

  /* Generate modes */
  UNUSED UINT8 retcode = SEOBNRv4HMROMCoreModes(&hlm, phiRef, fRef, distance,
                           Mtot_sec, q, chi1, chi2, freqs,
                                deltaF, nk_max,nModes,sign_odd_modes);
  if(retcode != XLAL_SUCCESS) XLAL_ERROR(retcode);

  /* GPS time for output frequency series and modes */

  COMPLEX16FrequencySeries *hlm_for_hphc = XLALSphHarmFrequencySeriesGetMode(hlm, 2,-2);  
  LIGOTimeGPS tGPS = hlm_for_hphc->epoch;
  size_t npts = hlm_for_hphc->data->length;

  /* Create output frequencyseries for hplus, hcross */
  UNUSED COMPLEX16FrequencySeries *hPlustilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &tGPS, 0.0, deltaF,&lalStrainUnit, npts);
  UNUSED COMPLEX16FrequencySeries *hCrosstilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &tGPS, 0.0, deltaF,&lalStrainUnit, npts);
  memset(hPlustilde->data->data, 0, npts * sizeof(COMPLEX16));
  memset(hCrosstilde->data->data, 0, npts * sizeof(COMPLEX16));
  XLALUnitDivide(&hPlustilde->sampleUnits, &hPlustilde->sampleUnits, &lalSecondUnit);
  XLALUnitDivide(&hCrosstilde->sampleUnits, &hCrosstilde->sampleUnits, &lalSecondUnit);
  retcode = SEOBROMComputehplushcrossFromhlm(hPlustilde,hCrosstilde,ModeArray,hlm,inclination,phiRef);
  if(retcode != XLAL_SUCCESS) XLAL_ERROR(retcode);

  /* Point the output pointers to the relevant time series and return */
  (*hptilde) = hPlustilde;
  (*hctilde) = hCrosstilde;

  XLALDestroyREAL8Sequence(freqs);
  XLALDestroySphHarmFrequencySeries(hlm);
  XLALDestroyValue(ModeArray);

  // This is a debug enviromental variable, when activated is deleting the ROM dataset from memory
  const char* s = getenv("FREE_MEMORY_SEOBNRv4HMROM");
  if (s!=NULL){
    for(unsigned int index_mode = 0; index_mode < nModes;index_mode++){
      SEOBNRROMdataDS_Cleanup(&__lalsim_SEOBNRv4HMROMDS_data[index_mode]);
    }
  }  

  return(XLAL_SUCCESS);
}

/**
 * Compute modes in LAL format for the SEOBNRv4HM_ROM model.
 *
 * Returns hlm modes as a SphHarmFrequencySeries with
 * equal spacing deltaF and contains zeros from zero frequency to the starting
 * frequency fLow and zeros beyond the cutoff frequency in the ringdown.
 */
int XLALSimIMRSEOBNRv4HMROM_Modes(
  SphHarmFrequencySeries **hlm, /**< Output: Frequency-domain hlm */
  UNUSED REAL8 phiRef,                                 /**< Phase at reference time */
  UNUSED REAL8 deltaF,                                 /**< Sampling frequency (Hz) */
  REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
  UNUSED REAL8 fHigh,                                  /**< End frequency; 0 defaults to Mf=0.14 */
  REAL8 fRef,                                        /**< Reference frequency unused here */
  REAL8 distance,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 chi1,                                   /**< Dimensionless aligned component spin 1 */
  REAL8 chi2,                                   /**< Dimensionless aligned component spin 2 */
  UNUSED INT4 nk_max,                                   /**< Truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1 */
  UNUSED UINT4 nModes                                  /**< Number of modes to use. This should be 1 for SEOBNRv4_ROM and 5 for SEOBNRv4HM_ROM */
)
{
  REAL8 sign_odd_modes = 1.;
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI) {
    // Swap m1 and m2
    REAL8 m1temp = m1SI;
    REAL8 chi1temp = chi1;
    m1SI = m2SI;
    chi1 = chi2;
    m2SI = m1temp;
    chi2 = chi1temp;
    sign_odd_modes = -1.; // When we swap m1 with m2 odd-m modes should get - sign because they depend on dm
  }
  /* Check if the number of modes makes sense */

  if(nModes >5){
    XLAL_PRINT_ERROR("Requested number of modes not available. Set nModes = 0 to get all the available modes.\n");
    XLAL_ERROR(XLAL_EDOM);
  }

  /* Get masses in terms of solar mass */
  REAL8 mass1 = m1SI / LAL_MSUN_SI;
  REAL8 mass2 = m2SI / LAL_MSUN_SI;
  REAL8 Mtot = mass1+mass2;
  REAL8 q = mass1/mass2;
  UNUSED REAL8 Mtot_sec = Mtot * LAL_MTSUN_SI;       /* Total mass in seconds */

  if(fRef==0.0)
    fRef=fLow;

  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequency we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = fLow;
  freqs->data[1] = fHigh;  

  /* Load ROM data if not loaded already */
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&SEOBNRv4HMROM_is_initialized, SEOBNRv4HMROM_Init_LALDATA);
#else
  SEOBNRv4HMROM_Init_LALDATA();
#endif

  /* Generate modes */
  UNUSED UINT8 retcode; 
  
  if(nModes == 0){
    retcode = SEOBNRv4HMROMCoreModes(hlm, phiRef, fRef, distance,
                                Mtot_sec, q, chi1, chi2, freqs,
                                deltaF, nk_max,5,sign_odd_modes);
  }
  else{
    retcode = SEOBNRv4HMROMCoreModes(hlm, phiRef, fRef, distance,
                                Mtot_sec, q, chi1, chi2, freqs,
                                deltaF, nk_max,nModes,sign_odd_modes);
  }


  XLALDestroyREAL8Sequence(freqs);

  return(XLAL_SUCCESS);
}
