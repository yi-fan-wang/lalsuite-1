/******************** <lalVerbatim file="StochasticOptimalFilterTestCV">
Author: UTB Relativity Group; contact J. T. Whelan
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{StochasticOptimalFilterTest.c}}
\label{stochastic:ss:StochasticOptimalFilterTest.c}

Test suite for \texttt{LALStochasticOptimalFilter()}.

\subsubsection*{Usage}
\begin{verbatim}
./StochasticOptimalFilterTest [options]
Options:
  -h             print usage message
  -q             quiet: run silently
  -v             verbose: print extra information
  -d level       set lalDebugLevel to level
  -n length      frequency series contain length points
  -f fRef        set normalization reference frequency to fRef
  -w filename    read gravitational-wave spectrum from file filename
  -g filename    read overlap reduction function from file filename
  -i filename    read first inverse noise PSD from file filename
  -j filename    read second inverse noise PSD from file filename
  -s filename    read first half-whitened inverse noise PSD from file filename
  -t filename    read second half-whitened inverse noise PSD from file filename
  -o filename    print optimal filter to file filename
  -y             use normalization appropriate to heterodyned data
\end{verbatim}

\subsubsection*{Description}

This program tests the function \texttt{LALStochasticOptimalFilter()},
which generates a normalized optimal filter from a stochastic
gravitational-wave background spectrum
$h_{100}^2\Omega_{\scriptstyle{\rm GW}}(f)$, an overlap reduction
function $\gamma(f)$, and unwhitened and half-whitened noise power
spectral densities $\{P_i(f),P^{\scriptstyle{\rm HW}}_i(f)\}$ for a
pair of detectors.

First, it tests that the correct error codes 
(\textit{cf.}\ Sec.~\ref{stochastic:s:StochasticCrossCorrelation.h})
are generated for the following error conditions (tests in
\textit{italics} are not performed if \verb+LAL_NDEBUG+ is set, as
the corresponding checks in the code are made using the ASSERT macro):
\begin{itemize}
\item \textit{null pointer to input structure}
\item \textit{null pointer to output series}
\item \textit{null pointer to overlap reduction function}
\item \textit{null pointer to gravitational-wave spectrum}
\item \textit{null pointer to first inverse noise PSD}
\item \textit{null pointer to second inverse noise PSD}
\item \textit{null pointer to first half-whitened inverse noise PSD}
\item \textit{null pointer to second half-whitened inverse noise PSD}
\item \textit{null pointer to data member of output series}
\item \textit{null pointer to data member of overlap reduction function}
\item \textit{null pointer to data member of gravitational-wave spectrum}
\item \textit{null pointer to data member of first inverse noise PSD}
\item \textit{null pointer to data member of second inverse noise PSD}
\item \textit{null pointer to data member of first half-whitened inverse noise PSD}
\item \textit{null pointer to data member of second half-whitened inverse noise PSD}
\item \textit{null pointer to data member of data member of output series}
\item \textit{null pointer to data member of data member of overlap reduction function}
\item \textit{null pointer to data member of data member of gravitational-wave spectrum}
\item \textit{null pointer to data member of data member of first inverse noise PSD}
\item \textit{null pointer to data member of data member of second inverse noise PSD}
\item \textit{null pointer to data member of data member of first half-whitened inverse noise PSD}
\item \textit{null pointer to data member of data member of second half-whitened inverse noise PSD}
\item \textit{zero length}
\item \textit{negative frequency spacing}
\item \textit{zero frequency spacing}
\item negative start frequency
\item length mismatch between overlap reduction function and output series
\item length mismatch between overlap reduction function and gravitational-wave spectrum
\item length mismatch between overlap reduction function and first inverse noise PSD
\item length mismatch between overlap reduction function and second inverse noise PSD
\item length mismatch between overlap reduction function and first half-whitened inverse noise PSD
\item length mismatch between overlap reduction function and second half-whitened inverse noise PSD
\item frequency spacing mismatch between overlap reduction function and gravitational-wave spectrum
\item frequency spacing mismatch between overlap reduction function and first inverse noise PSD
\item frequency spacing mismatch between overlap reduction function and second inverse noise PSD
\item frequency spacing mismatch between overlap reduction function and first half-whitened inverse noise PSD
\item frequency spacing mismatch between overlap reduction function and second half-whitened inverse noise PSD
\item start frequency mismatch between overlap reduction function and gravitational-wave spectrum
\item start frequency mismatch between overlap reduction function and first inverse noise PSD
\item start frequency mismatch between overlap reduction function and second inverse noise PSD
\item start frequency mismatch between overlap reduction function and first half-whitened inverse noise PSD
\item start frequency mismatch between overlap reduction function and second half-whitened inverse noise PSD
\item reference frequency less than frequency spacing
\item reference frequency greater than maximum frequency
\end{itemize}

It then verifies that the correct optimal filter is generated
(checking the normalization by verifying that (\ref{stochastic:e:mu})
is satisfied) for each of the following simple test cases:
\begin{enumerate}
\item $\gamma(f) = h_{100}^2\Omega_{\scriptstyle{\rm GW}}(f) = P_1(f) 
  =P_2(f)=P^{\scriptstyle{\rm HW}}_1(f)=P^{\scriptstyle{\rm HW}}_2(f)=1$;   
  The expected optimal filter in this case is
  $\widetilde{Q}(f)\propto f^{-3}$.
\item $\gamma(f) = P_1(f) = P_2(f) = P^{\scriptstyle{\rm HW}}_1(f)
  = P^{\scriptstyle{\rm HW}}_2(f)=1$;
  $h_{100}^2\Omega_{\scriptstyle{\rm GW}}(f)=f^3$.
  The expected optimal filter in this case is
  $\widetilde{Q}(f)=\textrm{constant}$.
\end{enumerate}

\subsubsection*{Exit codes}
\input{StochasticOptimalFilterTestCE}

\subsubsection*{Uses}
\begin{verbatim}
LALStochasticOptimalFilter()
LALCheckMemoryLeaks()
LALCReadFrequencySeries()
LALSCreateVector()
LALSDestroyVector()
LALCCreateVector()
LALCDestroyVector()
LALCHARCreateVector()
LALCHARDestroyVector()
LALUnitAsString()
LALUnitCompare()
getopt()
printf()
fprintf()
freopen()
fabs()
\end{verbatim}

\subsubsection*{Notes}
\begin{itemize}
\item No specific error checking is done on user-specified data.  If
  \texttt{length} is missing, the resulting default will cause a bad
  data error.  If \texttt{fRef} is unspecified, a default value of
  1\,Hz is used.
\item The length of the user-provided series must be specified, even
  though it could in principle be deduced from the input file, because
  the data sequences must be allocated before the
  \texttt{LALCReadFrequencySeries()} function is called.
\item If some, but not all, of the \texttt{filename} arguments are
  present, the user-specified data will be silently ignored.
\end{itemize}

\vfill{\footnotesize\input{StochasticOptimalFilterTestCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <config.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <lal/StochasticCrossCorrelation.h>
#include <lal/AVFactories.h>
#include <lal/ReadFTSeries.h>
#include <lal/PrintFTSeries.h>
#include <lal/Units.h>

#include "CheckStatus.h"

NRCSID(STOCHASTICOPTIMALFILTERTESTC,
       "$Id$");

#define STOCHASTICOPTIMALFILTERTESTC_TRUE     1
#define STOCHASTICOPTIMALFILTERTESTC_FALSE    0
#define STOCHASTICOPTIMALFILTERTESTC_DELTAF   1.0
#define STOCHASTICOPTIMALFILTERTESTC_F0       0.0
#define STOCHASTICOPTIMALFILTERTESTC_FREF     1.0
#define STOCHASTICOPTIMALFILTERTESTC_LENGTH   8
#define STOCHASTICOPTIMALFILTERTESTC_TOL      1e-6

extern char *optarg;
extern int   optind;

int lalDebugLevel = LALNDEBUG;
/* int lalDebugLevel  = LALMSGLVL3; */
BOOLEAN optVerbose = STOCHASTICOPTIMALFILTERTESTC_FALSE;
BOOLEAN optHetero = STOCHASTICOPTIMALFILTERTESTC_FALSE;
REAL8 optFRef      = 1.0;
UINT4 optLength    = 0.0;
CHAR  optOverlapFile[LALNameLength]   = "";
CHAR  optOmegaFile[LALNameLength]     = "";
CHAR  optInvNoise1File[LALNameLength] = "";
CHAR  optInvNoise2File[LALNameLength] = "";
CHAR  optHwInvNoise1File[LALNameLength] = "";
CHAR  optHwInvNoise2File[LALNameLength] = "";
CHAR  optOptimalFile[LALNameLength]     = "";

INT4 code;

static void Usage (const char *program, int exitflag);
static void ParseOptions (int argc, char *argv[]);
static REAL8 mu(const REAL4FrequencySeries*, const REAL4FrequencySeries*,
                const COMPLEX8FrequencySeries*);

/************* <lalErrTable file="StochasticOptimalFilterTestCE"> */
#define STOCHASTICOPTIMALFILTERTESTC_ENOM 0
#define STOCHASTICOPTIMALFILTERTESTC_EARG 1
#define STOCHASTICOPTIMALFILTERTESTC_ECHK 2
#define STOCHASTICOPTIMALFILTERTESTC_EFLS 3
#define STOCHASTICOPTIMALFILTERTESTC_EUSE 4

#define STOCHASTICOPTIMALFILTERTESTC_MSGENOM "Nominal exit"
#define STOCHASTICOPTIMALFILTERTESTC_MSGEARG "Error parsing command-line arguments"
#define STOCHASTICOPTIMALFILTERTESTC_MSGECHK "Error checking failed to catch bad data"
#define STOCHASTICOPTIMALFILTERTESTC_MSGEFLS "Incorrect answer for valid data"
#define STOCHASTICOPTIMALFILTERTESTC_MSGEUSE "Bad user-entered data"
/***************************** </lalErrTable> */

int main(int argc, char *argv[])
{

  static LALStatus         status;

  StochasticOptimalFilterInput             input;

  REAL4FrequencySeries     realBadData;
  COMPLEX8FrequencySeries  complexBadData;
  REAL4*                   realTempPtr;
  COMPLEX8*                complexTempPtr;
  LIGOTimeGPS              epoch = {1,0};

  REAL4FrequencySeries     overlap;
  REAL4FrequencySeries     omegaGW;
  REAL4FrequencySeries     invNoise1;
  REAL4FrequencySeries     invNoise2;
  COMPLEX8FrequencySeries  hwInvNoise1;
  COMPLEX8FrequencySeries  hwInvNoise2;
  COMPLEX8FrequencySeries  optimal;

  REAL8                    omegaRef;

  INT4       i;
  REAL8      f;
  REAL8      muTest;  /*refers to the calculated mu to test normalization*/
  REAL8      testNum;     /*temporary value used to check optimal output */
 
  LALUnitPair              unitPair;
  BOOLEAN                  result;

  CHARVector               *unitString = NULL;

  StochasticOptimalFilterParameters  params;

  ParseOptions (argc, argv);

  /* define valid params */ 
  overlap.name[0] = '\0';
  overlap.f0     = STOCHASTICOPTIMALFILTERTESTC_F0;
  overlap.deltaF = STOCHASTICOPTIMALFILTERTESTC_DELTAF;
  overlap.epoch  = epoch;
  overlap.data   = NULL;
  overlap.sampleUnits = lalDimensionlessUnit;

  realBadData = omegaGW = invNoise1 = invNoise2 = overlap;

  invNoise1.sampleUnits.unitNumerator[LALUnitIndexStrain] = -2;
  invNoise1.sampleUnits.unitNumerator[LALUnitIndexSecond] = -1;
  invNoise1.sampleUnits.powerOfTen = 36;
  invNoise2.sampleUnits = invNoise1.sampleUnits;

  overlap.sampleUnits.unitNumerator[LALUnitIndexStrain] = 2;
 
  params.fRef = STOCHASTICOPTIMALFILTERTESTC_FREF;
  params.heterodyned = STOCHASTICOPTIMALFILTERTESTC_FALSE;

  strncpy(overlap.name, "", LALNameLength);
  hwInvNoise1.name[0] = '\0';
  hwInvNoise1.f0     = overlap.f0;
  hwInvNoise1.deltaF = overlap.deltaF;
  hwInvNoise1.epoch  = overlap.epoch;
  hwInvNoise1.data   = NULL;
  hwInvNoise1.sampleUnits = lalDimensionlessUnit;

  complexBadData = optimal = hwInvNoise2 = hwInvNoise1;

  hwInvNoise1.sampleUnits.unitNumerator[LALUnitIndexStrain] = -1;
  hwInvNoise1.sampleUnits.unitNumerator[LALUnitIndexADCCount] = -1;
  hwInvNoise1.sampleUnits.unitNumerator[LALUnitIndexSecond] = -1;
  hwInvNoise1.sampleUnits.powerOfTen = 18;
  hwInvNoise2.sampleUnits = hwInvNoise1.sampleUnits;

  /* allocate memory */
  LALSCreateVector(&status, &(overlap.data),
                   STOCHASTICOPTIMALFILTERTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "", 
			    STOCHASTICOPTIMALFILTERTESTC_EFLS,
			    STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALSCreateVector(&status, &(omegaGW.data),
                   STOCHASTICOPTIMALFILTERTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "", 
			    STOCHASTICOPTIMALFILTERTESTC_EFLS,
			    STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALSCreateVector(&status, &(invNoise1.data),
                   STOCHASTICOPTIMALFILTERTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "", 
			    STOCHASTICOPTIMALFILTERTESTC_EFLS,
			    STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALSCreateVector(&status, &(invNoise2.data),
                   STOCHASTICOPTIMALFILTERTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "", 
			    STOCHASTICOPTIMALFILTERTESTC_EFLS,
			    STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALCCreateVector(&status, &(hwInvNoise1.data),
                   STOCHASTICOPTIMALFILTERTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "", 
			    STOCHASTICOPTIMALFILTERTESTC_EFLS,
			    STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALCCreateVector(&status, &(hwInvNoise2.data),
                   STOCHASTICOPTIMALFILTERTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "", 
			    STOCHASTICOPTIMALFILTERTESTC_EFLS,
			    STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALCCreateVector(&status, &(optimal.data),
                   STOCHASTICOPTIMALFILTERTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "", 
			    STOCHASTICOPTIMALFILTERTESTC_EFLS,
			    STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }

  input.overlapReductionFunction = &overlap;
  input.omegaGW = &omegaGW;
  input.unWhitenedInverseNoisePSD1 = &invNoise1;
  input.unWhitenedInverseNoisePSD2 = &invNoise2;
  input.halfWhitenedInverseNoisePSD1 = &hwInvNoise1;
  input.halfWhitenedInverseNoisePSD2 = &hwInvNoise2;

  /* TEST INVALID DATA HERE -------------------------------------- */
#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    /* test behavior for null pointer to input structure */
    LALStochasticOptimalFilter(&status, &optimal, NULL, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to input structure results in error:\n \"%s\"\n",STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* test behavior for null pointer to output series */
    LALStochasticOptimalFilter(&status, NULL, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to output series results in error:\n \"%s\"\n",STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    
    /* test behavior for null pointer to overlap reduction function */ 
    input.overlapReductionFunction = NULL;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to overlap reduction function results in error:\n      \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.overlapReductionFunction = &overlap;
    
    /* test behavior for null pointer to gravitational-wave spectrum */ 
    input.omegaGW = NULL;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to gravitational-wave spectrum results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.omegaGW = &omegaGW;
    
    /* test behavior for null pointer to inverse noise 1 member */
    /* of input structure */
    input.unWhitenedInverseNoisePSD1 = NULL;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to first inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.unWhitenedInverseNoisePSD1 = &invNoise1;
    
    /* test behavior for null pointer to inverse noise 2 member */
    /* of input structure*/ 
    input.unWhitenedInverseNoisePSD2 = NULL;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to second inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.unWhitenedInverseNoisePSD2 = &invNoise2;
    
    /* test behavior for null pointer to half-whitened inverse noise 1 member */
    /* of input structure */
    input.halfWhitenedInverseNoisePSD1 = NULL;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to first half-whitened inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.halfWhitenedInverseNoisePSD1 = &hwInvNoise1;
    
    /* test behavior for null pointer to second half-whitened inverse noise PSD member */
    /* of input structure */
    input.halfWhitenedInverseNoisePSD2 = NULL;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to second half-whitened inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.halfWhitenedInverseNoisePSD2 = &hwInvNoise2;
    
    /* test behavior for null pointer to data member of overlap */
    input.overlapReductionFunction = &realBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of overlap reduction function results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.overlapReductionFunction = &overlap;
    
    /* test behavior for null pointer to data member of omega */
    input.omegaGW = &realBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of gravitational-wave spectrum results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.omegaGW = &omegaGW;
    
    /* test behavior for null pointer to data member of first inverse noise PSD */
    input.unWhitenedInverseNoisePSD1 = &realBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of first inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.unWhitenedInverseNoisePSD1 = &invNoise1;
    
    /* test behavior for null pointer to data member of second inverse noise PSD */
    input.unWhitenedInverseNoisePSD2 = &realBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of second inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.unWhitenedInverseNoisePSD2 = &invNoise2;
    
    /* test behavior for null pointer to data member of half-whitened */
    /* first inverse noise PSD */
    input.halfWhitenedInverseNoisePSD1 = &complexBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of first half-whitened noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.halfWhitenedInverseNoisePSD1 = &hwInvNoise1;
    
    /* test behavior for null pointer to data member of half-whitened */
    /* second inverse noise PSD */
    input.halfWhitenedInverseNoisePSD2 = &complexBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of second half-whitened noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.halfWhitenedInverseNoisePSD2 = &hwInvNoise2;
    
    /* test behavior for null pointer to data member of output */
    /* frequency series */ 
    LALStochasticOptimalFilter(&status, &complexBadData, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of output series results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    
    
    /* Create a vector for testing REAL4 null data-data pointers */
    LALSCreateVector(&status, &(realBadData.data), STOCHASTICOPTIMALFILTERTESTC_LENGTH);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOPTIMALFILTERTESTC_EFLS,
			      STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    realTempPtr = realBadData.data->data;
    realBadData.data->data = NULL;

    /* test behavior for null pointer to data member of data member of */
    /* output frequency series */
    
    LALStochasticOptimalFilter(&status, &complexBadData, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of output series results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    
    /* test behavior for null pointer to data member of data member of overlap */
    input.overlapReductionFunction = &realBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of overlap reduction function results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.overlapReductionFunction = &overlap;
    
    /* test behavior for null pointer to data member of data member of omega */
    input.omegaGW = &realBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of gravitational-wave spectrum results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.omegaGW = &omegaGW;
    
    /* test behavior for null pointer to data member of data member of */
    /* first inverse noise PSD */
    input.unWhitenedInverseNoisePSD1 = &realBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of first inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.unWhitenedInverseNoisePSD1 = &invNoise1;
    
    /* test behavior for null pointer to data member of data member of */
    /* second inverse noise PSD */
    input.unWhitenedInverseNoisePSD2 = &realBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of second inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.unWhitenedInverseNoisePSD2 = &invNoise2;
    
    /* Create a vector for testing COMPLEX8 null data-data pointers */ 
    LALCCreateVector(&status, &(complexBadData.data), STOCHASTICOPTIMALFILTERTESTC_LENGTH);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOPTIMALFILTERTESTC_EFLS,
			      STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    complexTempPtr = complexBadData.data->data;
    complexBadData.data->data = NULL;
    
    /* test behavior for null pointer to data member of data member of */
    /* first half-whitened inverse noise PSD */
    input.halfWhitenedInverseNoisePSD1 = &complexBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of first half-whitened noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.halfWhitenedInverseNoisePSD1 = &hwInvNoise1;
    
    /* test behavior for null pointer to data member of data member of */
    /* second half-whitened inverse noise PSD */
    input.halfWhitenedInverseNoisePSD2 = &complexBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of second half-whitened noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.halfWhitenedInverseNoisePSD2 = &hwInvNoise2;
    
    /** clean up **/
    realBadData.data->data = realTempPtr;
    LALSDestroyVector(&status, &(realBadData.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOPTIMALFILTERTESTC_EFLS,
			      STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    
    complexBadData.data->data = complexTempPtr;
    LALCDestroyVector(&status, &(complexBadData.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOPTIMALFILTERTESTC_EFLS,
			      STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    
    /* test behavior for zero length */
    overlap.data->length = 
      omegaGW.data->length = invNoise1.data->length = 
      invNoise2.data->length = 
      hwInvNoise1.data->length = 
      hwInvNoise2.data->length = 
      optimal.data->length = 0;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EZEROLEN,
			      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: zero length results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);
    /* reassign valid length */
    overlap.data->length = 
      omegaGW.data->length = invNoise1.data->length = 
      invNoise2.data->length = 
      hwInvNoise1.data->length = 
      hwInvNoise2.data->length =
      optimal.data->length =STOCHASTICOPTIMALFILTERTESTC_LENGTH;
    
    /* test behavior for negative frequency spacing */
    overlap.deltaF = omegaGW.deltaF =
      invNoise1.deltaF = invNoise2.deltaF = 
      hwInvNoise1.deltaF = 
      hwInvNoise2.deltaF = -3.5;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
			      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: negative frequency spacing results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);
    /* reassign valid frequency spacing */
    overlap.deltaF = omegaGW.deltaF = 
      invNoise1.deltaF = invNoise2.deltaF = 
      hwInvNoise1.deltaF = 
      hwInvNoise2.deltaF = STOCHASTICOPTIMALFILTERTESTC_DELTAF;
    
    /* test behavior for zero frequency spacing */
    overlap.deltaF = omegaGW.deltaF =
      invNoise1.deltaF = invNoise2.deltaF = 
      hwInvNoise1.deltaF = 
      hwInvNoise2.deltaF = 0;
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
			      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF,
			      STOCHASTICOPTIMALFILTERTESTC_ECHK,
			      STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: zero frequency spacing results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);
    /* reassign valid frequency spacing */
    overlap.deltaF = 
      omegaGW.deltaF = invNoise1.deltaF = 
      invNoise2.deltaF = 
      hwInvNoise1.deltaF = 
      hwInvNoise2.deltaF =  STOCHASTICOPTIMALFILTERTESTC_DELTAF;
  } /* if ( ! lalNoDebug ) */ 
#endif /* LAL_NDEBUG */

  /* test behavior for negative start frequency */
  overlap.f0 = omegaGW.f0 
    = invNoise1.f0 
    = invNoise2.f0 
    = hwInvNoise1.f0 
    = hwInvNoise2.f0 = -3.0;

  LALStochasticOptimalFilter(&status, &optimal, &input, &params);

  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN,
			    STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: negative start frequency results in error:\n      \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  overlap.f0 = omegaGW.f0 =
    invNoise1.f0 = invNoise2.f0 = 
    hwInvNoise1.f0 = 
    hwInvNoise2.f0 = STOCHASTICOPTIMALFILTERTESTC_F0;

  /* test behavior length mismatch between overlap reduction function and output */
  optimal.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH - 1;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: length mismatch between overlap reduction function and output series results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  optimal.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH;
  
  /* test behavior for length mismatch between overlap reduction function and omega */
  omegaGW.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH - 1;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: length mismatch between overlap reduction function and gravitational-wave spectrum results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  omegaGW.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH;
  
  /* test behavior length mismatch between overlap reduction function and first inverse noise PSD */
  invNoise1.data->length 
    = STOCHASTICOPTIMALFILTERTESTC_LENGTH - 1;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
   {
     return code;
   }
  printf("  PASS: length mismatch between overlap reduction function and first inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  invNoise1.data->length 
    = STOCHASTICOPTIMALFILTERTESTC_LENGTH;
  
  /* test behavior length mismatch between overlap reduction function and second inverse noise PSD */
  invNoise2.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH - 1;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: length mismatch between overlap reduction function and second inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  invNoise2.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH;
  
  /* test behavior length mismatch between overlap reduction function and half-whitened */
  /* first inverse noise PSD */
  hwInvNoise1.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH - 1;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: length mismatch between overlap reduction function and first half-whitened inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  hwInvNoise1.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH;
  
  /* test behavior length mismatch between overlap reduction function and half-whitened */
  /* first inverse noise PSD */
  hwInvNoise2.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH - 1;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: length mismatch between overlap reduction function and second half-whitened inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  hwInvNoise2.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH;
  
  /* test behavior for frequency spacing mismatch between overlap reduction function and omega */
  omegaGW.deltaF = 2.0 * STOCHASTICOPTIMALFILTERTESTC_DELTAF;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: frequency spacing mismatch between overlap reduction function and gravitational-wave spectrum results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  omegaGW.deltaF = STOCHASTICOPTIMALFILTERTESTC_DELTAF;
  
  /* test behavior frequency spacing mismatch between overlap reduction function and first inverse noise PSD */
  invNoise1.deltaF = 2.0 * STOCHASTICOPTIMALFILTERTESTC_DELTAF;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: frequency spacing mismatch between overlap reduction function and first inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  invNoise1.deltaF = STOCHASTICOPTIMALFILTERTESTC_DELTAF;
  
  /* test behavior frequency spacing mismatch between overlap reduction function and inverse */
  /* noise 2 */
  invNoise2.deltaF = 2.0 * STOCHASTICOPTIMALFILTERTESTC_DELTAF;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: frequency spacing mismatch between overlap reduction function and second inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  invNoise2.deltaF = STOCHASTICOPTIMALFILTERTESTC_DELTAF;
  
  /* test behavior frequency spacing mismatch between overlap reduction function and */
  /* first half-whitened inverse noise PSD */
  hwInvNoise1.deltaF = 2.0 * STOCHASTICOPTIMALFILTERTESTC_DELTAF;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: frequency spacing mismatch between overlap reduction function and first half-whitened inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  hwInvNoise1.deltaF = STOCHASTICOPTIMALFILTERTESTC_DELTAF;
  
  /* test behavior frequency spacing mismatch between overlap reduction function and */
  /* second half-whitened inverse noise PSD */
  hwInvNoise2.deltaF = 305;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: frequency spacing mismatch between overlap reduction function and second half-whitened inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  hwInvNoise2.deltaF = STOCHASTICOPTIMALFILTERTESTC_DELTAF;
  
  /* test behavior for start frequency mismatch between overlap reduction function and omega */
  omegaGW.f0 = 30;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: start frequency mismatch between overlap reduction function and gravitational-wave spectrum results in error:\n         \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  omegaGW.f0 = STOCHASTICOPTIMALFILTERTESTC_F0;
  
  /* test behavior start frequency mismatch between overlap reduction function and first inverse noise PSD */
  invNoise1.f0 = 30;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: start frequency mismatch between overlap reduction function and first inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  invNoise1.f0 = STOCHASTICOPTIMALFILTERTESTC_F0;
  
  /* test behavior start frequency mismatch between overlap reduction function and second inverse noise PSD */
  invNoise2.f0 = 30;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: start frequency mismatch between overlap reduction function and second inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  invNoise2.f0 = STOCHASTICOPTIMALFILTERTESTC_F0;
  
  /* test behavior start frequency mismatch between overlap reduction function and half-whitened */
  /* first inverse noise PSD */
  hwInvNoise1.f0 = 30;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: start frequency mismatch between overlap reduction function and first half-whitened inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  hwInvNoise1.f0 = STOCHASTICOPTIMALFILTERTESTC_F0;
  
  /* test behavior start frequency mismatch between overlap reduction function and half-whitened */
  /* second inverse noise PSD */
  hwInvNoise2.f0 = 305;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: start frequency mismatch between overlap reduction function and second half-whitened inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  hwInvNoise2.f0 = STOCHASTICOPTIMALFILTERTESTC_F0;
  
  /* test behavior for reference frequency to be less than frequency spacing */
  params.fRef = STOCHASTICOPTIMALFILTERTESTC_DELTAF/2;
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EOORFREF,
			    STOCHASTICCROSSCORRELATIONH_MSGEOORFREF,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: reference frequency less than frequency spacing results in error:\n       \"%s\"\n",STOCHASTICCROSSCORRELATIONH_MSGEOORFREF);
  params.fRef = STOCHASTICOPTIMALFILTERTESTC_FREF;
 
  /* test behavior for reference frequency to be greater than its maximum */
  params.fRef = (STOCHASTICOPTIMALFILTERTESTC_LENGTH*STOCHASTICOPTIMALFILTERTESTC_DELTAF);
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EOORFREF,
			    STOCHASTICCROSSCORRELATIONH_MSGEOORFREF,
			    STOCHASTICOPTIMALFILTERTESTC_ECHK,
			    STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: reference frequency greater than maximum frewquency results in error:\n       \"%s\"\n",STOCHASTICCROSSCORRELATIONH_MSGEOORFREF);
  params.fRef = STOCHASTICOPTIMALFILTERTESTC_FREF;
  

 /* VALID TEST DATA HERE ----------------------------------------- */ 
  
  /* create input to test */
  overlap.data->data[0] = 1;
  omegaGW.data->data[0] = 0;
  invNoise1.data->data[0] = 0;
  invNoise2.data->data[0] = 0;
  hwInvNoise1.data->data[0].re = 0;
  hwInvNoise1.data->data[0].im = 0;
  hwInvNoise2.data->data[0].re = 0;
  hwInvNoise2.data->data[0].im = 0;

  /** Test 1 **/
  params.fRef = STOCHASTICOPTIMALFILTERTESTC_FREF;
  for (i=1; i < STOCHASTICOPTIMALFILTERTESTC_LENGTH; i++)
  {
    f = i*STOCHASTICOPTIMALFILTERTESTC_DELTAF;
    
    overlap.data->data[i] = 1;
    omegaGW.data->data[i] = 1;
    invNoise1.data->data[i] = 1;
    invNoise2.data->data[i] = 1;
    hwInvNoise1.data->data[i].re = 1;
    hwInvNoise1.data->data[i].im = 0;
    hwInvNoise2.data->data[i].re = 1;
    hwInvNoise2.data->data[i].im = 0;
  }

  /* fill optimal input */
  input.overlapReductionFunction    = &(overlap);
  input.omegaGW                     = &(omegaGW);
  input.halfWhitenedInverseNoisePSD1 = &(hwInvNoise1);
  input.halfWhitenedInverseNoisePSD2 = &(hwInvNoise2);
  input.unWhitenedInverseNoisePSD1             = &(invNoise1);
  input.unWhitenedInverseNoisePSD2             = &(invNoise2);
 
  /* calculate optimal filter */
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status,0, "",
			    STOCHASTICOPTIMALFILTERTESTC_EFLS,
			    STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }

  /******   check output   ******/

  /* check output f0 */
  if (optVerbose)
  {
    printf("f0=%g, should be %g\n", optimal.f0,
	   STOCHASTICOPTIMALFILTERTESTC_F0);
  }
  if ( fabs(optimal.f0-STOCHASTICOPTIMALFILTERTESTC_F0)
       > STOCHASTICOPTIMALFILTERTESTC_TOL )
  {
    printf("  FAIL: Valid data test #1\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n",
	     STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }
  
  /* check output deltaF */
  if (optVerbose)
  {
    printf("deltaF=%g, should be %g\n", optimal.deltaF,
	   STOCHASTICOPTIMALFILTERTESTC_DELTAF);
  }
  if ( fabs(optimal.deltaF-STOCHASTICOPTIMALFILTERTESTC_DELTAF)
       / STOCHASTICOPTIMALFILTERTESTC_DELTAF 
	> STOCHASTICOPTIMALFILTERTESTC_TOL )
  {
    printf("  FAIL: Valid data test #1\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n", STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
     return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }
  
  /* check output units */
  
  unitPair.unitOne = lalDimensionlessUnit;
  unitPair.unitOne.unitNumerator[LALUnitIndexADCCount] = -2;
  unitPair.unitTwo = optimal.sampleUnits;
  LALUnitCompare(&status, &result, &unitPair);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICOPTIMALFILTERTESTC_EFLS,
			    STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }
  
  if (optVerbose) 
  {
    LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOPTIMALFILTERTESTC_EFLS,
			      STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    
    LALUnitAsString( &status, unitString, &(unitPair.unitTwo) );
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOPTIMALFILTERTESTC_EFLS,
			      STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "Units are \"%s\", ", unitString->data );
    
    LALUnitAsString( &status, unitString, &(unitPair.unitOne) );
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOPTIMALFILTERTESTC_EFLS,
			      STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "should be \"%s\"\n", unitString->data );
    
    LALCHARDestroyVector(&status, &unitString);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOPTIMALFILTERTESTC_EFLS,
			      STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
  }
  
  if (!result)
  {
    printf("  FAIL: Valid data test #1\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n", 
	     STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  /* check output values */
  if (optVerbose) 
  {
    printf("Q(0)=%g + %g i, should be 0\n",
	   optimal.data->data[0].re, optimal.data->data[0].im);
  }
  if ( fabs(optimal.data->data[0].re) > STOCHASTICOPTIMALFILTERTESTC_TOL
       || fabs(optimal.data->data[0].im) 
       > STOCHASTICOPTIMALFILTERTESTC_TOL )
  {
    printf("  FAIL: Valid data test #1\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n", STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }
  
  for (i=1; i < STOCHASTICOPTIMALFILTERTESTC_LENGTH; i++) 
  {
    f = i*STOCHASTICOPTIMALFILTERTESTC_DELTAF;
    testNum = (STOCHASTICOPTIMALFILTERTESTC_DELTAF
	       *STOCHASTICOPTIMALFILTERTESTC_DELTAF
	       *STOCHASTICOPTIMALFILTERTESTC_DELTAF)
      /(f*f*f);
    if (optVerbose) 
    {
      printf("Q(%g Hz)/Re(Q(%g Hz))=%g + %g i, should be %g\n",
	     f, STOCHASTICOPTIMALFILTERTESTC_DELTAF,
	     optimal.data->data[i].re/optimal.data->data[1].re,
	     optimal.data->data[i].im/optimal.data->data[1].re,
	     testNum);
    }
    if (fabs(optimal.data->data[i].re/optimal.data->data[1].re 
	     - testNum)/testNum
	> STOCHASTICOPTIMALFILTERTESTC_TOL
	|| fabs(optimal.data->data[i].im/optimal.data->data[1].re)
	> STOCHASTICOPTIMALFILTERTESTC_TOL)
    {
      printf("  FAIL: Valid data test #1\n");
      if (optVerbose)
      {
	printf("Exiting with error: %s\n",
	       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
      }
      return STOCHASTICOPTIMALFILTERTESTC_EFLS;
    }
  }

  /* normalization costant */
  omegaRef = 1.0;
  muTest = mu(&omegaGW, &overlap, &optimal);
  if (optVerbose) 
  {
      printf("mu=%g, should be %g\n", muTest, omegaRef);
  }
  if ( fabs(muTest - omegaRef)/omegaRef > STOCHASTICOPTIMALFILTERTESTC_TOL)
  {
     printf("  FAIL: Valid data test #1\n" );
     return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  printf("  PASS: Valid data test #1\n");

  /** Test 2 **/
  params.fRef = STOCHASTICOPTIMALFILTERTESTC_FREF;
  for (i=0; i < STOCHASTICOPTIMALFILTERTESTC_LENGTH; i++)
     {
       f = i*STOCHASTICOPTIMALFILTERTESTC_DELTAF;

       overlap.data->data[i] = 1;
       omegaGW.data->data[i] = pow(f,3);
       invNoise1.data->data[i] = 1;
       invNoise2.data->data[i] = 1;
       hwInvNoise1.data->data[i].re = 1;
       hwInvNoise1.data->data[i].im = 0;
       hwInvNoise2.data->data[i].re = 1;
       hwInvNoise2.data->data[i].im = 0;
    }

  /* fill optimal input */
  input.overlapReductionFunction    = &(overlap);
  input.omegaGW                     = &(omegaGW);
  input.halfWhitenedInverseNoisePSD1 = &(hwInvNoise1);
  input.halfWhitenedInverseNoisePSD2 = &(hwInvNoise2);
  input.unWhitenedInverseNoisePSD1             = &(invNoise1);
  input.unWhitenedInverseNoisePSD2             = &(invNoise2);
    
  /* calculate optimal filter */
  LALStochasticOptimalFilter(&status, &optimal, &input, &params);
  if ( ( code = CheckStatus(&status,0, "", 
			    STOCHASTICOPTIMALFILTERTESTC_EFLS,
			    STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }

  /******   check output   ******/

  /* check output f0 */
  if (optVerbose)
  {
    printf("f0=%g, should be %g\n", optimal.f0,
	   STOCHASTICOPTIMALFILTERTESTC_F0);
  }
  if ( fabs(optimal.f0-STOCHASTICOPTIMALFILTERTESTC_F0)
       > STOCHASTICOPTIMALFILTERTESTC_TOL )
  {
    printf("  FAIL: Valid data test #2\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n",
	     STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }
  
  /* check output deltaF */
  if (optVerbose)
  {
    printf("deltaF=%g, should be %g\n", optimal.deltaF,
	   STOCHASTICOPTIMALFILTERTESTC_DELTAF);
  }
  if ( fabs(optimal.deltaF-STOCHASTICOPTIMALFILTERTESTC_DELTAF)
       / STOCHASTICOPTIMALFILTERTESTC_DELTAF 
	> STOCHASTICOPTIMALFILTERTESTC_TOL )
  {
    printf("  FAIL: Valid data test #2\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n", STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
     return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }
  
  /* check output units */
  
  unitPair.unitOne = lalDimensionlessUnit;
  unitPair.unitOne.unitNumerator[LALUnitIndexADCCount] = -2;
  unitPair.unitTwo = optimal.sampleUnits;
  LALUnitCompare(&status, &result, &unitPair);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICOPTIMALFILTERTESTC_EFLS,
			    STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }
  
  if (optVerbose) 
  {
    LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOPTIMALFILTERTESTC_EFLS,
			      STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    
    LALUnitAsString( &status, unitString, &(unitPair.unitTwo) );
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOPTIMALFILTERTESTC_EFLS,
			      STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "Units are \"%s\", ", unitString->data );
    
    LALUnitAsString( &status, unitString, &(unitPair.unitOne) );
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOPTIMALFILTERTESTC_EFLS,
			      STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "should be \"%s\"\n", unitString->data );
    
    LALCHARDestroyVector(&status, &unitString);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOPTIMALFILTERTESTC_EFLS,
			      STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
  }
  
  if (!result)
  {
    printf("  FAIL: Valid data test #2\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n", 
	     STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  /* check output values */
  if (optVerbose) 
  {
    printf("Q(0)=%g + %g i, should be 0\n",
	   optimal.data->data[0].re, optimal.data->data[0].im);
  }
  if ( fabs(optimal.data->data[0].re) > STOCHASTICOPTIMALFILTERTESTC_TOL
       || fabs(optimal.data->data[0].im) 
       > STOCHASTICOPTIMALFILTERTESTC_TOL )
  {
    printf("  FAIL: Valid data test #2\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n", STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  testNum = 1.0;
  for (i=1; i < STOCHASTICOPTIMALFILTERTESTC_LENGTH; i++) 
  {
    f = i*STOCHASTICOPTIMALFILTERTESTC_DELTAF;
    if (optVerbose) 
    {
      printf("Q(%g Hz)/Re(Q(%g Hz))=%g + %g i, should be %g\n",
	     f, STOCHASTICOPTIMALFILTERTESTC_DELTAF,
	     optimal.data->data[i].re/optimal.data->data[1].re,
	     optimal.data->data[i].im/optimal.data->data[1].re,
	     testNum);
    }
    if (fabs(optimal.data->data[i].re/optimal.data->data[1].re 
	     - testNum)/testNum
	> STOCHASTICOPTIMALFILTERTESTC_TOL
	|| fabs(optimal.data->data[i].im/optimal.data->data[1].re)
	> STOCHASTICOPTIMALFILTERTESTC_TOL)
    {
      printf("  FAIL: Valid data test #2\n");
      if (optVerbose)
      {
	printf("Exiting with error: %s\n",
	       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
      }
      return STOCHASTICOPTIMALFILTERTESTC_EFLS;
    }
  }

  /* normalization costant */
  omegaRef = params.fRef * params.fRef * params.fRef;
  muTest = mu(&omegaGW, &overlap, &optimal);
  if (optVerbose) 
  {
      printf("mu=%g, should be %g\n", muTest, omegaRef);
  }
  if ( fabs(muTest - omegaRef)/omegaRef > STOCHASTICOPTIMALFILTERTESTC_TOL)
  {
     printf("  FAIL: Valid data test #2\n" );
     return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  printf("  PASS: Valid data test #2\n");

  /* clean up valid data */
  LALSDestroyVector(&status, &(overlap.data));
  if ( ( code = CheckStatus (&status, 0 , "", 
			     STOCHASTICOPTIMALFILTERTESTC_EFLS,
			     STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALSDestroyVector(&status, &(omegaGW.data));
  if ( ( code = CheckStatus (&status, 0 , "", 
			     STOCHASTICOPTIMALFILTERTESTC_EFLS,
			     STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALSDestroyVector(&status, &(invNoise1.data));
  if ( ( code = CheckStatus (&status, 0 , "", 
			     STOCHASTICOPTIMALFILTERTESTC_EFLS,
			     STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALSDestroyVector(&status, &(invNoise2.data));
  if ( ( code = CheckStatus (&status, 0 , "", 
			     STOCHASTICOPTIMALFILTERTESTC_EFLS,
			     STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALCDestroyVector(&status, &(hwInvNoise1.data));
  if ( ( code = CheckStatus (&status, 0 , "", 
			     STOCHASTICOPTIMALFILTERTESTC_EFLS,
			     STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALCDestroyVector(&status, &(hwInvNoise2.data));
  if ( ( code = CheckStatus (&status, 0 , "", 
			     STOCHASTICOPTIMALFILTERTESTC_EFLS,
			     STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALCDestroyVector(&status, &(optimal.data));
  if ( ( code = CheckStatus (&status, 0 , "", 
			     STOCHASTICOPTIMALFILTERTESTC_EFLS,
			     STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
       return code;
  }
  LALCheckMemoryLeaks();

 printf("PASS: all tests\n");


 /* VALID USER TEST DATA HERE ----------------------------------------- */

  if (optOverlapFile[0] &&  optOmegaFile[0] && optInvNoise1File[0] && optInvNoise2File[0]
      && optHwInvNoise1File[0] && optHwInvNoise2File[0] && optOptimalFile[0])
  { 
    
    /* allocate memory */ 
    overlap.data     = NULL;
    omegaGW.data     = NULL;
    invNoise1.data   = NULL;
    invNoise2.data   = NULL;
    hwInvNoise1.data = NULL;
    hwInvNoise2.data = NULL;
    optimal.data     = NULL;
    
    LALSCreateVector(&status, &(overlap.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "", 
			       STOCHASTICOPTIMALFILTERTESTC_EFLS,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALSCreateVector(&status, &(omegaGW.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "", 
			       STOCHASTICOPTIMALFILTERTESTC_EFLS,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSCreateVector(&status, &(invNoise1.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "", 
			       STOCHASTICOPTIMALFILTERTESTC_EFLS,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSCreateVector(&status, &(invNoise2.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "", 
			       STOCHASTICOPTIMALFILTERTESTC_EFLS,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCCreateVector(&status, &(hwInvNoise1.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "", 
			       STOCHASTICOPTIMALFILTERTESTC_EFLS,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCCreateVector(&status, &(hwInvNoise2.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "", 
			       STOCHASTICOPTIMALFILTERTESTC_EFLS,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCCreateVector(&status, &(optimal.data),optLength);
    if ( ( code = CheckStatus (&status, 0 , "", 
			       STOCHASTICOPTIMALFILTERTESTC_EFLS,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    
    /* Read input files */
    LALSReadFrequencySeries(&status, &overlap,   optOverlapFile);
    LALSReadFrequencySeries(&status, &omegaGW,   optOmegaFile);
    LALSReadFrequencySeries(&status, &invNoise1, optInvNoise1File);
    LALSReadFrequencySeries(&status, &invNoise2, optInvNoise2File);
    LALCReadFrequencySeries(&status, &hwInvNoise1, optHwInvNoise1File);
    LALCReadFrequencySeries(&status, &hwInvNoise2, optHwInvNoise2File);
    
    /* fill optimal input */
    input.overlapReductionFunction    = &(overlap);
    input.omegaGW                     = &(omegaGW);
    input.halfWhitenedInverseNoisePSD1 = &(hwInvNoise1);
    input.halfWhitenedInverseNoisePSD2 = &(hwInvNoise2);
    input.unWhitenedInverseNoisePSD1             = &(invNoise1);
    input.unWhitenedInverseNoisePSD2             = &(invNoise2);
    
    params.fRef = optFRef;
    params.heterodyned = optHetero;
    
    /* calculate optimal filter */
    LALStochasticOptimalFilter(&status, &optimal, &input, &params);
    if ( ( code = CheckStatus (&status, 0 , "", 
			       STOCHASTICOPTIMALFILTERTESTC_EUSE,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEUSE) ) )
    {
      return code;
    }
    
    LALCPrintFrequencySeries(&optimal, optOptimalFile);
    printf("=========== Optimal Filter Written to File %s ===========\n",
	   optOptimalFile);
    
    /* clean up */
    LALSDestroyVector(&status, &(overlap.data));
    if ( ( code = CheckStatus (&status, 0 , "",
			       STOCHASTICOPTIMALFILTERTESTC_EFLS,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSDestroyVector(&status, &(omegaGW.data));
    if ( ( code = CheckStatus (&status, 0 , "",
			       STOCHASTICOPTIMALFILTERTESTC_EFLS,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSDestroyVector(&status, &(invNoise1.data));
    if ( ( code = CheckStatus (&status, 0 , "", 
			       STOCHASTICOPTIMALFILTERTESTC_EFLS,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSDestroyVector(&status, &(invNoise2.data));
    if ( ( code = CheckStatus (&status, 0 , "", 
			       STOCHASTICOPTIMALFILTERTESTC_EFLS,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCDestroyVector(&status, &(hwInvNoise1.data));
    if ( ( code = CheckStatus (&status, 0 , "",
			       STOCHASTICOPTIMALFILTERTESTC_EFLS,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCDestroyVector(&status, &(hwInvNoise2.data));
    if ( ( code = CheckStatus (&status, 0 , "",
			       STOCHASTICOPTIMALFILTERTESTC_EFLS,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCDestroyVector(&status, &(optimal.data));
    if ( ( code = CheckStatus (&status, 0 , "",
			       STOCHASTICOPTIMALFILTERTESTC_EFLS,
			       STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
  }
  LALCheckMemoryLeaks();
  return 0;

};

/*----------------------------------------------------------------------*/
static REAL8 mu(const REAL4FrequencySeries* omegaGW, 
                const REAL4FrequencySeries* overlap,
                const COMPLEX8FrequencySeries* optimal)
{
  INT4       i;
  REAL8      f;
  REAL8      constant;
  REAL8      f3;
  REAL8      deltaF;
  INT4       length;
  REAL8      mu = 0.0;

  deltaF = omegaGW->deltaF;
  length = omegaGW->data->length;

  constant = 3.0L * ( (LAL_H0FAC_SI*1.0e+18) * (LAL_H0FAC_SI*1.0e+18) )
             / ( 20.0L * (LAL_PI*LAL_PI) );

  /* calculate mu */
  for(i=1; i < (length);i++)
    {
      f = i*deltaF;
      f3 = f*f*f;
      mu += (2.0  * constant * (omegaGW->data->data[i]) *
            (overlap->data->data[i]) * (optimal->data->data[i].re))/f3;  
    }
    return mu;
}

/* Usage () Message */
/* Prints a usage message for program program and exits with code exitcode.*/

static void Usage (const char *program, int exitcode)
{
  fprintf (stderr, "Usage: %s [options]\n", program);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "  -h             print this message\n");
  fprintf (stderr, "  -q             quiet: run silently\n");
  fprintf (stderr, "  -v             verbose: print extra information\n");
  fprintf (stderr, "  -d level       set lalDebugLevel to level\n");
  fprintf (stderr, "  -n length      frequency series contain length points\n");
  fprintf (stderr, "  -f fRef        set normalization reference frequency to fRef\n");
  fprintf (stderr, "  -w filename    read gravitational-wave spectrum from file filename\n");
  fprintf (stderr, "  -g filename    read overlap reduction function from file filename\n"); 
  fprintf (stderr, "  -i filename    read first inverse noise PSD from file filename\n");
  fprintf (stderr, "  -j filename    read second inverse noise PSD from file filename\n");
  fprintf (stderr, "  -s filename    read first half-whitened inverse noise PSD from file filename\n");
  fprintf (stderr, "  -t filename    read second half-whitened inverse noise PSD from file filename\n");
  fprintf (stderr, "  -o filename    print optimal filter to file filename\n");
  fprintf (stderr, "  -y             use normalization appropriate to heterodyned data\n");
  exit (exitcode);
}

/*
 * ParseOptions ()
 *
 * Parses the argc - 1 option strings in argv[].
 *
 */
static void
ParseOptions (int argc, char *argv[])
{
  while (1)
  {
    int c = -1;

    c = getopt (argc, argv, "hqvd:n:f:w:g:i:j:s:t:o:y");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'o': /* specify output file */
        strncpy (optOptimalFile, optarg, LALNameLength);
        break;
        
      case 'f': /* specify refernce frequency */
        optFRef = atoi (optarg);
        break;

      case 'n': /* specify number of points in frequency series */
        optLength = atoi (optarg);
        break;
        
      case 'w': /* specify omegaGW file */
        strncpy (optOmegaFile, optarg, LALNameLength);
        break;
        
      case 'g': /* specify overlap file */
        strncpy (optOverlapFile, optarg, LALNameLength);
        break;
        
      case 'i': /* specify InvNoise1 file */
        strncpy (optInvNoise1File, optarg, LALNameLength);
        break;

      case 'j': /* specify InvNoise2 file */
        strncpy (optInvNoise2File, optarg, LALNameLength);
        break;

      case 's': /* specify hwInvNoise1 file */
        strncpy (optHwInvNoise1File, optarg, LALNameLength);
        break;

      case 't': /* specify hwInvNoise2 file */
        strncpy (optHwInvNoise2File, optarg, LALNameLength);
        break;

      case 'd': /* set debug level */
        lalDebugLevel = atoi (optarg);
        break;

      case 'v': /* optVerbose */
        optVerbose = STOCHASTICOPTIMALFILTERTESTC_TRUE;
        break;

      case 'y': /* optHetero */
        optHetero = STOCHASTICOPTIMALFILTERTESTC_TRUE;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        freopen ("/dev/null", "w", stderr);
        freopen ("/dev/null", "w", stdout);
        break;

      case 'h':
        Usage (argv[0], 0);
        break;

      default:
        Usage (argv[0], 1);
    }

  }

  if (optind < argc)
  {
    Usage (argv[0], 1);
  }

  return;
}
