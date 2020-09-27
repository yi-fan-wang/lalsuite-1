#ifndef _DATETYPES_H
#define _DATETYPES_H

#ifdef newc
#include <ccomplex>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace std;
#else
#include <complex.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#endif
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>


/* Integer types */
typedef unsigned int UINT4;        /**< Four-byte unsigned integer. */

/* Real types */
typedef float REAL4;    /**< Single precision real floating-point number (4 bytes). */
typedef double REAL8;   /**< Double precision real floating-point number (8 bytes). */
typedef double complex COMPLEX16;
//typedef struct tagCOMPLEX16 { REAL8 re; REAL8 im; } COMPLEX16;

typedef int16_t  INT2;        /**< Two-byte signed integer */
typedef int32_t  INT4;        /**< Four-byte signed integer. */
typedef int64_t  INT8;        /**< Eight-byte signed integer; on some platforms this is equivalent to <tt>long int</tt> instead. */
typedef uint16_t UINT2;        /**< Two-byte unsigned integer. */
typedef uint64_t UINT8;        /**< Eight-byte unsigned integer; on some platforms this is equivalent to <tt>unsigned long int</tt> instead. */



//#define LAL_INT8_C INT64_C
#define LAL_INT8_C long int

//#define I (XLALCOMPLEX16Rect(0.0,1.0))
#define XLAL_BILLION_REAL8 1e9
#define XLAL_BILLION_INT8 LAL_INT8_C( 1000000000 )
#define EOB_RD_EFOLDS 10.0
#define LAL_MTSUN_SI  4.92549095e-6   /**< Geometrized solar mass, s */
#define LAL_PI_2      1.5707963267948966192313216916397514  /**< pi/2 */
#define LAL_1_PI      0.3183098861837906715377675267450287  /**< 1/pi */
#define LAL_MSUN_SI   1.98892e30      /**< Solar mass, kg */
#define LAL_GAMMA     0.5772156649015328606065120900824024  /**< gamma */
#define LAL_E         2.7182818284590452353602874713526625  /**< e */
#define LAL_MRSUN_SI  1.47662504e3    /**< Geometrized solar mass, m */
#define LAL_TWOPI     6.2831853071795864769252867665590058  /**< 2*pi */

#define LALFree(p) free(p)


typedef struct tagGSParams {
    int ampO;                 /**< twice PN order of the amplitude */
    REAL8 phiRef;             /**< phase at fRef */
    REAL8 deltaT;             /**< sampling interval */
    REAL8 m1;                 /**< mass of companion 1 */
    REAL8 m2;                 /**< mass of companion 2 */
    REAL8 f_min;              /**< start frequency */
    REAL8 e0;                 /**< eccentricity at start frequency */
    REAL8 f_max;              /**< end frequency */
    REAL8 distance;           /**< distance of source */
    REAL8 inclination;        /**< inclination of L relative to line of sight */
    REAL8 s1x;                /**< (x,y,z) components of spin of m1 body */
    REAL8 s1y;                /**< z-axis along line of sight, L in x-z plane */
    REAL8 s1z;                /**< dimensionless spin, Kerr bound: |s1| <= 1 */
    REAL8 s2x;                /**< (x,y,z) component ofs spin of m2 body */
    REAL8 s2y;                /**< z-axis along line of sight, L in x-z plane */
    REAL8 s2z;                /**< dimensionless spin, Kerr bound: |s2| <= 1 */
    int output;               /**< out put to file or not(0 if do not write to file) */
    char outname[256];        /**< file to which output should be written */
    int ampPhase;
    int verbose;
    char jobtag[256];
} GSParams;

typedef struct
tagLALUnit
{
    int  unitNumerator[7]; /**< Array of unit power numerators. */
    int  unitDenominatorMinusOne[7]; /**< Array of unit power denominators-minus-one. */
}
LALUnit;

/** Vector of type REAL8, see \ref ss_Vector for more details. */
typedef struct tagREAL8Vector
{
    UINT4  length; /**< Number of elements in array. */
    REAL8 *data; /**< Pointer to the data array. */
}
REAL8Vector;

typedef REAL8Vector     REAL8Sequence;    /**< See \ref ss_Sequence for documentation */

typedef struct
tagLIGOTimeGPS
{
    int gpsSeconds; /**< Seconds since 0h UTC 6 Jan 1980. */
    int gpsNanoSeconds; /**< Residual nanoseconds. */
}
LIGOTimeGPS;


/** Time series of REAL8 data, see \ref ss_TimeSeries for more details. */
typedef struct
tagREAL8TimeSeries
{
    char           name[64]; /**< The name of the time series. */
    LIGOTimeGPS    epoch; /**< The start time of the time series. */
    REAL8          deltaT; /**< The time step between samples of the time series in seconds. */
    REAL8          f0; /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit        sampleUnits; /**< The physical units of the quantity being sampled. */
    //REAL8Sequence *data; /**< The sequence of sampled data. */
    REAL8Vector *data; /**< The sequence of sampled data. */
}
REAL8TimeSeries;

typedef struct tagCOMPLEX16Vector
{
    UINT4      length; /**< Number of elements in array. */
    COMPLEX16 *data; /**< Pointer to the data array. */
}
COMPLEX16Vector;

typedef struct
tagEOBNonQCCoeffs
{
    REAL8 a1;
    REAL8 a2;
    REAL8 a3;
    REAL8 a3S;
    REAL8 a4;
    REAL8 a5;
    REAL8 b1;
    REAL8 b2;
    REAL8 b3;
    REAL8 b4;
} EOBNonQCCoeffs;

typedef struct
tagark4GSLIntegrator
{
    gsl_odeiv_step    *step;
    gsl_odeiv_control *control;
    gsl_odeiv_evolve  *evolve;
    
    gsl_odeiv_system  *sys;
    
    int (* dydt) (double t, const double y[], double dydt[], void * params);
    int (* stop) (double t, const double y[], double dydt[], void * params);
    
    int retries;        /* retries with smaller step when derivatives encounter singularity */
    int stopontestonly;    /* stop only on test, use tend to size buffers only */
    
    int returncode;
} ark4GSLIntegrator;

typedef struct
tagUINT4Vector
{
    UINT4  length; /**< Number of elements in array. */
    UINT4  *data; /**< Pointer to the data array. */
}
UINT4Vector;

/** Multidimentional array of REAL8, see \ref ss_Array for more details. */
typedef struct
tagREAL8Array
{
    UINT4Vector *dimLength; /**< Vector of array dimensions. */
    REAL8       *data; /**< Pointer to the data array. */
}
REAL8Array;

typedef struct
tagEOBACoefficients
{
    REAL8 n4;
    REAL8 n5;
    REAL8 d0;
    REAL8 d1;
    REAL8 d2;
    REAL8 d3;
    REAL8 d4;
    REAL8 d5;
}
EOBACoefficients;

typedef struct
tagFacWaveformCoeffs
{
    REAL8 delta22vh3;
    REAL8 delta22vh6;
    REAL8 delta22v8;
    REAL8 delta22vh9;
    REAL8 delta22v5;
    
    REAL8 rho22v2;
    REAL8 rho22v3;
    REAL8 rho22v4;
    REAL8 rho22v5;
    REAL8 rho22v6;
    REAL8 rho22v6l;
    REAL8 rho22v7;
    REAL8 rho22v8;
    REAL8 rho22v8l;
    REAL8 rho22v10;
    REAL8 rho22v10l;
    
    REAL8 delta21vh3;
    REAL8 delta21vh6;
    REAL8 delta21vh7;
    REAL8 delta21vh9;
    REAL8 delta21v5;
    REAL8 delta21v7;
    
    REAL8 rho21v1;
    REAL8 rho21v2;
    REAL8 rho21v3;
    REAL8 rho21v4;
    REAL8 rho21v5;
    REAL8 rho21v6;
    REAL8 rho21v6l;
    REAL8 rho21v7;
    REAL8 rho21v7l;
    REAL8 rho21v8;
    REAL8 rho21v8l;
    REAL8 rho21v10;
    REAL8 rho21v10l;
    
    REAL8 f21v1;
    
    REAL8 delta33vh3;
    REAL8 delta33vh6;
    REAL8 delta33vh9;
    REAL8 delta33v5;
    REAL8 delta33v7;
    
    REAL8 rho33v2;
    REAL8 rho33v3;
    REAL8 rho33v4;
    REAL8 rho33v5;
    REAL8 rho33v6;
    REAL8 rho33v6l;
    REAL8 rho33v7;
    REAL8 rho33v8;
    REAL8 rho33v8l;
    
    REAL8 f33v3;
    
    REAL8 delta32vh3;
    REAL8 delta32vh4;
    REAL8 delta32vh6;
    REAL8 delta32vh9;
    
    REAL8 rho32v;
    REAL8 rho32v2;
    REAL8 rho32v3;
    REAL8 rho32v4;
    REAL8 rho32v5;
    REAL8 rho32v6;
    REAL8 rho32v6l;
    REAL8 rho32v8;
    REAL8 rho32v8l;
    
    REAL8 delta31vh3;
    REAL8 delta31vh6;
    REAL8 delta31vh7;
    REAL8 delta31vh9;
    REAL8 delta31v5;
    
    REAL8 rho31v2;
    REAL8 rho31v3;
    REAL8 rho31v4;
    REAL8 rho31v5;
    REAL8 rho31v6;
    REAL8 rho31v6l;
    REAL8 rho31v7;
    REAL8 rho31v8;
    REAL8 rho31v8l;
    
    REAL8 f31v3;
    
    REAL8 delta44vh3;
    REAL8 delta44vh6;
    REAL8 delta44v5;
    
    REAL8 rho44v2;
    REAL8 rho44v3;
    REAL8 rho44v4;
    REAL8 rho44v5;
    REAL8 rho44v6;
    REAL8 rho44v6l;
    
    REAL8 delta43vh3;
    REAL8 delta43vh4;
    REAL8 delta43vh6;
    
    REAL8 rho43v;
    REAL8 rho43v2;
    REAL8 rho43v4;
    REAL8 rho43v5;
    REAL8 rho43v6;
    REAL8 rho43v6l;
    
    REAL8 f43v;
    
    REAL8 delta42vh3;
    REAL8 delta42vh6;
    
    REAL8 rho42v2;
    REAL8 rho42v3;
    REAL8 rho42v4;
    REAL8 rho42v5;
    REAL8 rho42v6;
    REAL8 rho42v6l;
    
    REAL8 delta41vh3;
    REAL8 delta41vh4;
    REAL8 delta41vh6;
    
    REAL8 rho41v;
    REAL8 rho41v2;
    REAL8 rho41v4;
    REAL8 rho41v5;
    REAL8 rho41v6;
    REAL8 rho41v6l;
    
    REAL8 f41v;
    
    REAL8 delta55vh3;
    REAL8 delta55v5;
    REAL8 rho55v2;
    REAL8 rho55v3;
    REAL8 rho55v4;
    REAL8 rho55v5;
    REAL8 rho55v6;
    
    REAL8 delta54vh3;
    REAL8 delta54vh4;
    REAL8 rho54v2;
    REAL8 rho54v3;
    REAL8 rho54v4;
    
    REAL8 delta53vh3;
    REAL8 rho53v2;
    REAL8 rho53v3;
    REAL8 rho53v4;
    REAL8 rho53v5;
    
    REAL8 delta52vh3;
    REAL8 delta52vh4;
    REAL8 rho52v2;
    REAL8 rho52v3;
    REAL8 rho52v4;
    
    REAL8 delta51vh3;
    REAL8 rho51v2;
    REAL8 rho51v3;
    REAL8 rho51v4;
    REAL8 rho51v5;
    
    REAL8 delta66vh3;
    REAL8 rho66v2;
    REAL8 rho66v3;
    REAL8 rho66v4;
    
    REAL8 delta65vh3;
    REAL8 rho65v2;
    REAL8 rho65v3;
    
    REAL8 delta64vh3;
    REAL8 rho64v2;
    REAL8 rho64v3;
    REAL8 rho64v4;
    
    REAL8 delta63vh3;
    REAL8 rho63v2;
    REAL8 rho63v3;
    
    REAL8 delta62vh3;
    REAL8 rho62v2;
    REAL8 rho62v3;
    REAL8 rho62v4;
    
    REAL8 delta61vh3;
    REAL8 rho61v2;
    REAL8 rho61v3;
    
    REAL8 delta77vh3;
    REAL8 rho77v2;
    REAL8 rho77v3;
    
    REAL8 rho76v2;
    
    REAL8 delta75vh3;
    REAL8 rho75v2;
    REAL8 rho75v3;
    
    REAL8 rho74v2;
    
    REAL8 delta73vh3;
    REAL8 rho73v2;
    REAL8 rho73v3;
    
    REAL8 rho72v2;
    
    REAL8 delta71vh3;
    REAL8 rho71v2;
    REAL8 rho71v3;
    
    REAL8 rho88v2;
    REAL8 rho87v2;
    REAL8 rho86v2;
    REAL8 rho85v2;
    REAL8 rho84v2;
    REAL8 rho83v2;
    REAL8 rho82v2;
    REAL8 rho81v2;
}
FacWaveformCoeffs;

typedef
struct tagNewtonMultipolePrefixes
{
    COMPLEX16 values[8+1][8+1];
}
NewtonMultipolePrefixes;

typedef
struct tagEOBParams
{
    REAL8 eta;
    REAL8 omega;
    REAL8 m1;
    REAL8 m2;
    EOBACoefficients        *aCoeffs;
    FacWaveformCoeffs       *hCoeffs;
    EOBNonQCCoeffs          *nqcCoeffs;
    NewtonMultipolePrefixes *prefixes;
}
EOBParams;

typedef struct
tagSpinEOBHCoeffs
{
    double KK;
    double k0;
    double k1;
    double k2;
    double k3;
    double k4;
    double b3;
    double bb3;
}
SpinEOBHCoeffs;

typedef struct
tagSpinEOBParams
{
    EOBParams               *eobParams;
    SpinEOBHCoeffs          *seobCoeffs;
    REAL8Vector             *sigmaStar;
    REAL8Vector             *sigmaKerr;
    REAL8                   a;
    int                     alignedSpins;
    int                     tortoise;
    
    // we assume m1>m2
    REAL8Vector             *Spin1;
    REAL8Vector             *Spin2;
}
SpinEOBParams;

typedef
struct tagSEOBRootParams
{
    REAL8          values[12]; /**<< Dynamical variables, x, y, z, px, py, pz, S1x, S1y, S1z, S2x, S2y and S2z */
    SpinEOBParams *params;     /**<< Spin EOB parameters -- physical, pre-computed, etc. */
    REAL8          omega;      /**<< Orbital frequency */
}
SEOBRootParams;


/* We need to encapsulate the data for the GSL derivative function */
typedef
struct tagHcapDerivParams
{
    const REAL8   *values;
    SpinEOBParams *params;
    UINT4         varyParam;
}
HcapDerivParams;

#define XLAL_REAL8_FAIL_NAN_INT LAL_INT8_C(0x7ff80000000001a1) /**< Hexadecimal representation of <tt>REAL8</tt> NaN failure bit pattern */
#define XLAL_IS_REAL8_FAIL_NAN(val) XLALIsREAL8FailNaN(val) /**< Tests if <tt>val</tt> is a XLAL <tt>REAL8</tt> failure NaN. */

/* We need to encapsulate the data for calculating spherical 2nd derivatives */
typedef
struct tagHcapSphDeriv2Params
{
    const REAL8     *sphValues;
    SpinEOBParams   *params;
    UINT4           varyParam1;
    UINT4           varyParam2;
}
HcapSphDeriv2Params;


/** XLAL error numbers and return values. */
enum XLALErrorValue {
    XLAL_SUCCESS =  0, /**< Success return value (not an error number) */
    XLAL_FAILURE = -1, /**< Failure return value (not an error number) */
    
    /* these are standard error numbers */
    XLAL_EIO     =  5,  /**< I/O error */
    XLAL_ENOMEM  = 12,  /**< Memory allocation error */
    XLAL_EFAULT  = 14,  /**< Invalid pointer */
    XLAL_EINVAL  = 22,  /**< Invalid argument */
    XLAL_EDOM    = 33,  /**< Input domain error */
    XLAL_ERANGE  = 34,  /**< Output range error */
    
    /* extended error numbers start at 128 ...
     * should be beyond normal errnos */
    
    /* these are common errors for XLAL functions */
    XLAL_EFAILED = 128, /**< Generic failure */
    XLAL_EBADLEN = 129, /**< Inconsistent or invalid length */
    XLAL_ESIZE   = 130, /**< Wrong size */
    XLAL_EDIMS   = 131, /**< Wrong dimensions */
    XLAL_ETYPE   = 132, /**< Wrong or unknown type */
    XLAL_ETIME   = 133, /**< Invalid time */
    XLAL_EFREQ   = 134, /**< Invalid freqency */
    XLAL_EUNIT   = 135, /**< Invalid units */
    XLAL_ENAME   = 136, /**< Wrong name */
    XLAL_EDATA   = 137, /**< Invalid data */
    
    /* user-defined errors */
    XLAL_EUSR0   = 200, /**< User-defined error 0 */
    XLAL_EUSR1   = 201, /**< User-defined error 1 */
    XLAL_EUSR2   = 202, /**< User-defined error 2 */
    XLAL_EUSR3   = 203, /**< User-defined error 3 */
    XLAL_EUSR4   = 204, /**< User-defined error 4 */
    XLAL_EUSR5   = 205, /**< User-defined error 5 */
    XLAL_EUSR6   = 206, /**< User-defined error 6 */
    XLAL_EUSR7   = 207, /**< User-defined error 7 */
    XLAL_EUSR8   = 208, /**< User-defined error 8 */
    XLAL_EUSR9   = 209, /**< User-defined error 9 */
    
    /* external or internal errors */
    XLAL_ESYS    = 254, /**< System error */
    XLAL_EERR    = 255, /**< Internal error */
    
    /* specific mathematical and numerical errors start at 256 */
    
    /* IEEE floating point errors */
    XLAL_EFPINVAL  = 256, /**< IEEE Invalid floating point operation, eg sqrt(-1), 0/0 */
    XLAL_EFPDIV0   = 257, /**< IEEE Division by zero floating point error */
    XLAL_EFPOVRFLW = 258, /**< IEEE Floating point overflow error */
    XLAL_EFPUNDFLW = 259, /**< IEEE Floating point underflow error */
    XLAL_EFPINEXCT = 260, /**< IEEE Floating point inexact error */
    
    /* numerical algorithm errors */
    XLAL_EMAXITER  = 261, /**< Exceeded maximum number of iterations */
    XLAL_EDIVERGE  = 262, /**< Series is diverging */
    XLAL_ESING     = 263, /**< Apparent singularity detected */
    XLAL_ETOL      = 264, /**< Failed to reach specified tolerance */
    XLAL_ELOSS     = 265, /**< Loss of accuracy */
    
    /* failure from within a function call: "or" error number with this */
    XLAL_EFUNC     = 1024 /**< Internal function call failed bit: "or" this with existing error number */
};

typedef struct
tagREAL8VectorSequence
{
    UINT4  length; /**< The number \a l of vectors. */
    UINT4  vectorLength; /**< The length \a n of each vector. */
    REAL8 *data; /**< Pointer to the data array.  Element \a i of vector \a j is \c data[ \a jn + \a i \c ]. */
}
REAL8VectorSequence;

enum
{
    LALUnitIndexMeter,     /**< The meter index. */
    LALUnitIndexKiloGram, /**< The kilogram index. */
    LALUnitIndexSecond,     /**< The second index. */
    LALUnitIndexAmpere,     /**< The ampere index. */
    LALUnitIndexKelvin,     /**< The kelvin index. */
    LALUnitIndexStrain,     /**< The strain index. */
    LALUnitIndexADCCount, /**< The ADC counts index. */
    LALNumUnits         /**< The number of units. */
};

#define LAL_PI        3.1415926535897932384626433832795029  /**< pi */

enum enumLALNameLength { LALNameLength = 64 };

#endif /* _DATE_H */




#ifndef _LALCONSTANTS_H
#define _LALCONSTANTS_H

#ifdef  __cplusplus
extern "C" {
#endif
    
    /** \name Floating-point constants
     * The following constants define the precision and range of
     * floating-point arithmetic in LAL.  They are taken from the IEEE
     * standard 754 for binary arithmetic.  All numbers are dimensionless. */
    /*@{*/
#define LAL_REAL4_MANT 24 /**< Bits of precision in the mantissa of a REAL4 */
#define LAL_REAL4_MAX 3.40282347e+38 /**< Largest REAL4 */
#define LAL_REAL4_MIN 1.17549435e-38 /**< Smallest nonzero REAL4 */
#define LAL_REAL4_EPS 1.19209290e-07 /**< 0.5^(LAL_REAL4_MANT-1), ie the difference between 1 and the next resolveable REAL4 */
#define LAL_REAL8_MANT 53 /**< Bits of precision in the mantissa of a REAL8 */
#define LAL_REAL8_MAX 1.7976931348623157e+308 /**< Largest REAL8 */
#define LAL_REAL8_MIN 2.2250738585072014e-308 /**< Smallest nonzero REAL8 */
#define LAL_REAL8_EPS 2.2204460492503131e-16  /**< 0.5^(LAL_REAL8_MANT-1), ie the difference between 1 and the next resolveable REAL8 */
    /*@}*/
    
    /** \name Mathematical constants
     * The following are fundamental mathematical constants.  They are mostly
     * taken from the GNU C <tt>math.h</tt> header (with the exception of
     * <tt>LAL_TWOPI</tt>, which was computed using Maple).  All numbers are
     * dimensionless. The value of exp(gamma) is taken from
     * http://www.research.att.com/~njas/sequences/A073004 */
    /*@{*/
#define LAL_E         2.7182818284590452353602874713526625  /**< e */
#define LAL_LOG2E     1.4426950408889634073599246810018922  /**< log_2 e */
#define LAL_LOG10E    0.4342944819032518276511289189166051  /**< log_10 e */
#define LAL_LN2       0.6931471805599453094172321214581766  /**< log_e 2 */
#define LAL_LN10      2.3025850929940456840179914546843642  /**< log_e 10 */
#define LAL_SQRT2     1.4142135623730950488016887242096981  /**< sqrt(2) */
#define LAL_SQRT1_2   0.7071067811865475244008443621048490  /**< 1/sqrt(2) */
#define LAL_GAMMA     0.5772156649015328606065120900824024  /**< gamma */
#define LAL_EXPGAMMA  1.7810724179901979852365041031071795  /**< exp(gamma) */
    /* Assuming we're not near a black hole or in Tennessee... */
#define LAL_PI        3.1415926535897932384626433832795029  /**< pi */
#define LAL_TWOPI     6.2831853071795864769252867665590058  /**< 2*pi */
#define LAL_PI_2      1.5707963267948966192313216916397514  /**< pi/2 */
#define LAL_PI_4      0.7853981633974483096156608458198757  /**< pi/4 */
#define LAL_1_PI      0.3183098861837906715377675267450287  /**< 1/pi */
#define LAL_2_PI      0.6366197723675813430755350534900574  /**< 2/pi */
#define LAL_2_SQRTPI  1.1283791670955125738961589031215452  /**< 2/sqrt(pi) */
#define LAL_PI_180    1.7453292519943295769236907684886127e-2 /**< pi/180 */
#define LAL_180_PI    57.295779513082320876798154814105170 /**< 180/pi */
    /*@}*/
    
    /** \name Exact physical constants
     * The following physical constants are defined to have exact values.
     * The values of \f$c\f$ and \f$g\f$ are taken from \ref Barnet_1996,
     * \f$p_\mathrm{atm}\f$ is from \ref Lang_1992, while \f$\epsilon_0\f$ and
     * \f$\mu_0\f$ are computed from \f$c\f$ using exact formulae.  The use
     * of a Julian year (365.25 days) as standard is specified by the IAU.
     * They are given in the SI units shown. */
    /*@{*/
#define LAL_C_SI      299792458 /**< Speed of light in vacuo, m s^-1 */
#define LAL_EPSILON0_SI  8.8541878176203898505365630317107503e-12 /**< Permittivity of free space, C^2 N^-1 m^-2 */
#define LAL_MU0_SI    1.2566370614359172953850573533118012e-6 /**< Permeability of free space, N A^-2 */
#define LAL_GEARTH_SI 9.80665 /**< Standard gravity, m s^-2 */
#define LAL_PATM_SI 101325 /**< Standard atmosphere, Pa */
#define LAL_YRJUL_SI 31557600 /**< Julian year, s */
#define LAL_LYR_SI 9.4607304725808e15 /**< (Julian) Lightyear, m */
    /*@}*/
    
    /** \name Physical constants
     * The following are measured fundamental physical constants, with values
     * given in \ref Barnet_1996.  When not dimensionless, they are given
     * in the SI units shown. */
    /*@{*/
#define LAL_G_SI      6.67259e-11    /**< Gravitational constant, N m^2 kg^-2 */
#define LAL_H_SI      6.6260755e-34  /**< Planck constant, J s */
#define LAL_HBAR_SI   1.05457266e-34 /**< Reduced Planck constant, J s */
#define LAL_MPL_SI    2.17671e-8     /**< Planck mass, kg */
#define LAL_LPL_SI    1.61605e-35    /**< Planck length, m */
#define LAL_TPL_SI    5.39056e-44    /**< Planck time, s */
#define LAL_K_SI      1.380658e-23   /**< Boltzmann constant, J K^-1 */
#define LAL_R_SI      8.314511       /**< Ideal gas constant, J K^-1 */
#define LAL_MOL       6.0221367e23   /**< Avogadro constant, dimensionless */
#define LAL_BWIEN_SI  2.897756e-3    /**< Wien displacement law constant, m K */
#define LAL_SIGMA_SI  5.67051e-8  /**< Stefan-Boltzmann constant, W m^-2 K^-4 */
#define LAL_AMU_SI    1.6605402e-27  /**< Atomic mass unit, kg */
#define LAL_MP_SI     1.6726231e-27  /**< Proton mass, kg */
#define LAL_ME_SI     9.1093897e-31  /**< Electron mass, kg */
#define LAL_QE_SI     1.60217733e-19 /**< Electron charge, C */
#define LAL_ALPHA  7.297354677e-3 /**< Fine structure constant, dimensionless */
#define LAL_RE_SI     2.81794092e-15 /**< Classical electron radius, m */
#define LAL_LAMBDAE_SI 3.86159323e-13 /**< Electron Compton wavelength, m */
#define LAL_AB_SI     5.29177249e-11 /**< Bohr radius, m */
#define LAL_MUB_SI    9.27401543e-24 /**< Bohr magneton, J T^-1 */
#define LAL_MUN_SI    5.05078658e-27 /**< Nuclear magneton, J T^-1 */
    /*@}*/
    
    /** \name Astrophysical parameters
     * The following parameters are derived from measured properties of the
     * Earth and Sun.  The values are taken from \ref Barnet_1996, except
     * for the obliquity of the ecliptic plane and the eccentricity of
     * Earth's orbit, which are taken from \ref Lang_1992.  All values are
     * given in the SI units shown.  Note that the ``year'' and
     * ``light-year'' have exactly defined values, and appear under
     * ``Exact physical constants''.
     */
    /*@{*/
#define LAL_REARTH_SI 6.378140e6      /**< Earth equatorial radius, m */
#define LAL_AWGS84_SI 6.378137e6      /**< Semimajor axis of WGS-84 Reference Ellipsoid, m */
#define LAL_BWGS84_SI 6.356752314e6   /**< Semiminor axis of WGS-84 Reference Ellipsoid, m */
#define LAL_MEARTH_SI 5.97370e24      /**< Earth mass, kg */
#define LAL_IEARTH    0.409092804     /**< Earth inclination (2000), radians */
#define LAL_EEARTH    0.0167          /**< Earth orbital eccentricity */
#define LAL_RSUN_SI   6.960e8         /**< Solar equatorial radius, m */
#define LAL_MSUN_SI   1.98892e30      /**< Solar mass, kg */
#define LAL_MRSUN_SI  1.47662504e3    /**< Geometrized solar mass, m */
#define LAL_MTSUN_SI  4.92549095e-6   /**< Geometrized solar mass, s */
#define LAL_LSUN_SI   3.846e26        /**< Solar luminosity, W */
#define LAL_AU_SI     1.4959787066e11 /**< Astronomical unit, m */
#define LAL_PC_SI     3.0856775807e16 /**< Parsec, m */
#define LAL_YRTROP_SI 31556925.2      /**< Tropical year (1994), s */
#define LAL_YRSID_SI  31558149.8      /**< Sidereal year (1994), s */
#define LAL_DAYSID_SI 86164.09053     /**< Mean sidereal day, s */
    /*@}*/
    
    /** \name Cosmological parameters
     * The following cosmological parameters are derived from measurements of
     * the Hubble expansion rate and of the cosmic background radiation
     * (CBR).  Data are taken from \ref Barnet_1996.  In what follows, the
     * normalized Hubble constant \f$h_0\f$ is equal to the actual Hubble
     * constant \f$H_0\f$ divided by \f$\langle H
     * \rangle=100\,\mathrm{km}\,\mathrm{s}^{-1}\mathrm{Mpc}^{-1}\f$.  Thus the
     * Hubble constant can be written as:
     * \f$H_0 = \langle H \rangle h_0\f$.
     * Similarly, the critical energy density \f$\rho_c\f$ required for spatial
     * flatness is given by: \f$\rho_c = \langle\rho\rangle h_0^2\f$.
     * Current estimates give \f$h_0\f$ a value of around 0.65, which is what is
     * assumed below.  All values are in the SI units shown. */
    /*@{*/
#define LAL_H0FAC_SI  3.2407792903e-18 /**< Hubble constant prefactor, s^-1 */
#define LAL_H0_SI     2e-18            /**< Approximate Hubble constant, s^-1 */
    /* Hubble constant H0 = h0*HOFAC, where h0 is around 0.65 */
#define LAL_RHOCFAC_SI 1.68860e-9   /**< Critical density prefactor, J m^-3 */
#define LAL_RHOC_SI   7e-10         /**< Approximate critical density, J m^-3 */
    /* Critical density RHOC = h0*h0*RHOCFAC, where h0 is around 0.65 */
#define LAL_TCBR_SI   2.726   /**< Cosmic background radiation temperature, K */
#define LAL_VCBR_SI   3.695e5 /**< Solar velocity with respect to CBR, m s^-1 */
#define LAL_RHOCBR_SI 4.177e-14 /**< Energy density of CBR, J m^-3 */
#define LAL_NCBR_SI   4.109e8   /**< Number density of CBR photons, m^-3 */
#define LAL_SCBR_SI   3.993e-14 /**< Entropy density of CBR, J K^-1 m^-3 */
    /*@}*/
    
    /*@}*/
#ifdef  __cplusplus
}
#endif

#endif /* _LALCONSTANTS_H */

#define PI M_PI


REAL8 XLALGPSGetREAL8( const LIGOTimeGPS *epoch );

COMPLEX16Vector* XLALCreateCOMPLEX16Vector(UINT4 length);
void XLALDestroyCOMPLEX16Vector(COMPLEX16Vector* v);

REAL8Vector* XLALCreateREAL8Vector(UINT4 length);
void XLALDestroyREAL8Vector(REAL8Vector* v);

UINT4Vector* XLALCreateUINT4Vector(UINT4 length);
void XLALDestroyUINT4Vector(UINT4Vector* v);

REAL8Array* XLALCreateREAL8ArrayL(UINT4 ndim,...);
void XLALDestroyREAL8Array(REAL8Array* v);

REAL8VectorSequence * XLALCreateREAL8VectorSequence ( UINT4 length, UINT4 veclen );
void XLALDestroyREAL8VectorSequence ( REAL8VectorSequence * vseq );

REAL8TimeSeries *XLALCreateREAL8TimeSeries ( const char *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaT, const LALUnit *sampleUnits, int length );
void XLALDestroyREAL8TimeSeries( REAL8TimeSeries * series );

int XLALSimIMREOBFinalMassSpin_re(
                               REAL8    *finalMass,    /**<< OUTPUT, the final mass (scaled by original total mass) */
                               REAL8    *finalSpin,    /**<< OUTPUT, the final spin (scaled by final mass) */
                               const REAL8     mass1,        /**<< The mass of the 1st component of the system */
                               const REAL8     mass2,        /**<< The mass of the 2nd component of the system */
                               const REAL8     spin1[3],    /**<< The spin of the 1st object; only needed for spin waveforms */
                               const REAL8     spin2[3]    /**<< The spin of the 2nd object; only needed for spin waveforms */
);

int XLALSimIMREOBGenerateQNMFreqV2_re(
                                   COMPLEX16Vector *modefreqs, /**<< OUTPUT, complex freqs of overtones in unit of Hz */
                                   const REAL8      mass1,     /**<< The mass of the 1st component (in Solar masses) */
                                   const REAL8      mass2,     /**<< The mass of the 2nd component (in Solar masses) */
                                   const REAL8      spin1[3],  /**<< The spin of the 1st object; only needed for spin waveforms */
                                   const REAL8      spin2[3],  /**<< The spin of the 2nd object; only needed for spin waveforms */
                                   UINT4            l,         /**<< The l value of the mode in question */
                                   UINT4            m,         /**<< The m value of the mode in question */
                                   UINT4            nmodes    /**<< The number of overtones that should be included (max 8) */
);

int XLALSimIMRSpinEOBCalculateSigmaStar(
                                        REAL8Vector *sigmaStar, /**<< OUTPUT, normalized (to total mass) spin of test particle */
                                        REAL8        mass1,     /**<< mass 1 */
                                        REAL8        mass2,     /**<< mass 2 */
                                        REAL8Vector *s1,        /**<< spin vector 1 */
                                        REAL8Vector *s2         /**<< spin vector 2 */);

int XLALSimIMRSpinEOBCalculateSigmaKerr(
                                        REAL8Vector *sigmaKerr, /**<< OUTPUT, normalized (to total mass) spin of deformed-Kerr */
                                        REAL8        mass1,     /**<< mass 1 */
                                        REAL8        mass2,     /**<< mass 2 */
                                        REAL8Vector *s1,        /**<< spin vector 1 */
                                        REAL8Vector *s2         /**<< spin vector 2 */);

int XLALSimIMREOBCalcSpinFacWaveformCoefficients(
                                                 FacWaveformCoeffs * const coeffs, /**< OUTPUT, pre-computed waveform coefficients */
                                                 const REAL8               m1,     /**< mass 1 */
                                                 const REAL8               m2,     /**< mass 2 */
                                                 const REAL8               eta,    /**< symmetric mass ratio */
                                                 const REAL8               a,      /**< Kerr spin parameter for test-particle terms */
                                                 const REAL8               chiS,   /**< (chi1+chi2)/2 */
                                                 const REAL8               chiA    /**< (chi1-chi2)/2 */
);

int XLALSimIMREOBComputeNewtonMultipolePrefixes(
                                                NewtonMultipolePrefixes *prefix, /**<< OUTPUT Structure containing the coeffs */
                                                const REAL8             m1,      /**<< Mass of first component */
                                                const REAL8             m2       /**<< Nass of second component */
);

int XLALSimIMRSpinEOBInitialConditions(
                                       REAL8Vector   *initConds, /**<< OUTPUT, Initial dynamical variables */
                                       const REAL8    mass1,     /**<< mass 1 */
                                       const REAL8    mass2,     /**<< mass 2 */
                                       const REAL8    fMin,      /**<< Initial frequency (given) */
                                       const REAL8    e0,        /**<< eccentricity at starting GW frequency (Hz) */
                                       const REAL8    inc,       /**<< Inclination */
                                       const REAL8    spin1[],   /**<< Initial spin vector 1 */
                                       const REAL8    spin2[],   /**<< Initial spin vector 2 */
                                       SpinEOBParams *params     /**<< Spin EOB parameters */
);

int
CalculateRotationMatrix(
                        gsl_matrix *rotMatrix,  /**< OUTPUT, rotation matrix */
                        gsl_matrix *rotInverse, /**< OUTPUT, rotation matrix inversed */
                        REAL8       r[],        /**< position vector */
                        REAL8       v[],        /**< velocity vector */
                        REAL8       L[]         /**< orbital angular momentum */
);
REAL8 CalculateCrossProduct( const int i, const REAL8 a[], const REAL8 b[] );

int NormalizeVector( REAL8 a[] );
int ApplyRotationMatrix(
                        gsl_matrix *rotMatrix, /**< rotation matrix */
                        REAL8      a[]         /**< OUTPUT, vector rotated */
);

int
XLALFindSphericalOrbit( const gsl_vector *x, /**<< Parameters requested by gsl root finder */
                       void *params,        /**<< Spin EOB parameters */
                       gsl_vector *f        /**<< Function values for the given parameters */
);

REAL8 XLALSpinHcapNumDerivWRTParam(
                                   const INT4 paramIdx,      /**<< Index of the parameters */
                                   const REAL8 values[],     /**<< Dynamical variables */
                                   SpinEOBParams *funcParams /**<< EOB Parameters */
);
double GSLSpinHamiltonianWrapper( double x, void *params );

int XLALSimIMRCalculateSpinEOBHCoeffs(
                                      SpinEOBHCoeffs *coeffs, /**<< OUTPUT, EOB parameters including pre-computed coefficients */
                                      const REAL8    eta,     /**<< symmetric mass ratio */
                                      const REAL8    a        /**<< Normalized deformed Kerr spin */
);


REAL8 XLALSimIMRSpinEOBHamiltonian(
                                   const REAL8    eta,                  /**<< Symmetric mass ratio */
                                   REAL8Vector    *x,         /**<< Position vector */
                                   REAL8Vector    *p,        /**<< Momentum vector (tortoise radial component pr*) */
                                   REAL8Vector    *sigmaKerr, /**<< Spin vector sigma_kerr */
                                   REAL8Vector    *sigmaStar, /**<< Spin vector sigma_star */
                                   INT4                      tortoise,  /**<< flag to state whether the momentum is the tortoise co-ord */
                                   SpinEOBHCoeffs *coeffs               /**<< Structure containing various coefficients */
);


int XLALSimIMRSpinEOBGetSpinFactorizedWaveform(
                                               COMPLEX16         * hlm,    /**< OUTPUT, hlm waveforms */
                                               REAL8Vector       * values, /**< dyanmical variables */
                                               const REAL8         v,               /**< velocity */
                                               const REAL8         Hreal,           /**< real Hamiltonian */
                                               const int          l,               /**< l mode index */
                                               const int          m,               /**< m mode index */
                                               SpinEOBParams     * params  /**< Spin EOB parameters */
);

REAL8 XLALCalculateSphHamiltonianDeriv2(
                                        const int      idx1,     /**<< Derivative w.r.t. index 1 */
                                        const int      idx2,     /**<< Derivative w.r.t. index 2 */
                                        const REAL8    values[], /**<< Dynamical variables in spherical coordinates */
                                        SpinEOBParams *params    /**<< Spin EOB Parameters */
);
int SphericalToCartesian(
                         REAL8 qCart[],      /**<< OUTPUT, position vector in Cartesean coordinates */
                         REAL8 pCart[],      /**<< OUTPUT, momentum vector in Cartesean coordinates */
                         const REAL8 qSph[], /**<< position vector in spherical coordinates */
                         const REAL8 pSph[]  /**<< momentum vector in spherical coordinates */
);
double GSLSpinHamiltonianDerivWrapper( double x,    /**<< Derivative at x */
                                      void  *params /**<< Function parameters */);
REAL8 XLALInspiralSpinFactorizedFlux(
                                     REAL8Vector           *values, /**< dynamical variables */
                                     const REAL8           omega,   /**< orbital frequency */
                                     SpinEOBParams         *ak,     /**< physical parameters */
                                     const REAL8            H,      /**< real Hamiltonian */
                                     const int             lMax    /**< upper limit of the summation over l */
);


int
XLALSimIMRSpinEOBCalculateNewtonianMultipole(
                                             COMPLEX16 *multipole, /**<< OUTPUT, Newtonian multipole */
                                             REAL8 x,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
                                             REAL8 r,       /**<< Orbital separation (units of total mass M */
                                             REAL8 phi,            /**<< Orbital phase (in radians) */
                                             UINT4  l,             /**<< Mode l */
                                             int  m,              /**<< Mode m */
                                             EOBParams *params     /**<< Pre-computed coefficients, parameters, etc. */
);
int
XLALScalarSphHarmThetaPiBy2(
                            COMPLEX16 *y, /**<< OUTPUT, Ylm(0,phi) */
                            int l,       /**<< Mode l */
                            int  m,      /**<< Mode m */
                            REAL8 phi     /**<< Orbital phase (in radians) */
);
REAL8
XLALAssociatedLegendreXIsZero( const int l,
                              const int m );

REAL8 XLALSimIMRSpinEOBHamiltonianDeltaR(
                                         SpinEOBHCoeffs *coeffs, /**<< Pre-computed coefficients which appear in the function */
                                         const REAL8    r,       /**<< Current orbital radius (in units of total mass) */
                                         const REAL8    eta,     /**<< Symmetric mass ratio */
                                         const REAL8    a        /**<< Normalized deformed Kerr spin */
);

REAL8 XLALSimIMRSpinEOBHamiltonianDeltaT(
                                         SpinEOBHCoeffs *coeffs, /**<< Pre-computed coefficients which appear in the function */
                                         const REAL8    r,       /**<< Current orbital radius (in units of total mass) */
                                         const REAL8    eta,     /**<< Symmetric mass ratio */
                                         const REAL8    a        /**<< Normalized deformed Kerr spin */
);
ark4GSLIntegrator *XLALAdaptiveRungeKutta4Init( int dim,
                                               int (* dydt) (double t, const double y[], double dydt[], void * params),  /* These are XLAL functions! */
                                               int (* stop) (double t, const double y[], double dydt[], void * params),
                                               double eps_abs, double eps_rel
                                               );
void XLALAdaptiveRungeKutta4Free( ark4GSLIntegrator *integrator );
int XLALSpinAlignedHcapDerivative(
                                  double  t,          /**< UNUSED */
                                  const REAL8   values[],   /**< dynamical varables */
                                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                                  void         *funcParams  /**< EOB parameters */
);
double GSLSpinAlignedHamiltonianWrapper( double x, void *params );

int
XLALEOBSpinAlignedStopCondition(double  t,  /**< UNUSED */
                                const double values[], /**< dynamical variable values */
                                double dvalues[],      /**< dynamical variable time derivative values */
                                void *funcParams       /**< physical parameters */
);

int XLALAdaptiveRungeKutta4( ark4GSLIntegrator *integrator,
                            void *params,
                            REAL8 *yinit,
                            REAL8 tinit, REAL8 tend, REAL8 deltat,
                            REAL8Array **yout
                            );
int XLALSpinAlignedHiSRStopCondition(double t,  /**< UNUSED */
                                     const double values[], /**< dynamical variable values */
                                     double dvalues[],      /**< dynamical variable time derivative values */
                                     void *funcParams       /**< physical parameters */
);
// scaling omega with (vy/(vy+vx))^gamma
REAL8
XLALSimIMRSpinAlignedEOBCalcOmega(
                                  const REAL8           values[],   /**<< Dynamical variables */
                                  SpinEOBParams         *funcParams /**<< EOB parameters */
);
REAL8
XMYSimIMRSpinAlignedEOBCalcOmega(
                                 const REAL8           values[],   /**<< Dynamical variables */
                                 SpinEOBParams         *funcParams /**<< EOB parameters */
);
LIGOTimeGPS * XLALGPSAdd( LIGOTimeGPS *epoch, REAL8 dt );
LIGOTimeGPS * XLALGPSSetREAL8( LIGOTimeGPS *epoch, REAL8 t );
LIGOTimeGPS * XLALGPSAddGPS( LIGOTimeGPS *epoch, const LIGOTimeGPS *dt );
LIGOTimeGPS * XLALGPSSet( LIGOTimeGPS *epoch, INT4 gpssec, INT8 gpsnan );
LIGOTimeGPS * XLALINT8NSToGPS( LIGOTimeGPS *epoch, INT8 ns );
INT8 XLALGPSToINT8NS( const LIGOTimeGPS *epoch );

int XLALSimIMRGetEOBCalibratedSpinNQC( EOBNonQCCoeffs *coeffs,
                                      INT4  l,
                                      INT4  m,
                                      REAL8 eta,
                                      REAL8 a );
int XLALSimIMRSpinEOBCalculateNQCCoefficients(
                                              REAL8Vector    *amplitude,   /**<< Waveform amplitude, func of time */
                                              REAL8Vector    *phase,       /**<< Waveform phase(rad), func of time */
                                              REAL8Vector    *rVec,        /**<< Position-vector, function of time */
                                              REAL8Vector    *prVec,       /**<< Momentum vector, function of time */
                                              REAL8Vector    *orbOmegaVec, /**<< Orbital frequency, func of time */
                                              INT4                      l,           /**<< Mode index l */
                                              INT4                      m,           /**<< Mode index m */
                                              REAL8                     timePeak,    /**<< Time of peak orbital frequency */
                                              REAL8                     deltaT,      /**<< Sampling interval */
                                              REAL8                     eta,         /**<< Symmetric mass ratio */
                                              REAL8                     a,           /**<< Normalized spin of deformed-Kerr */
                                              EOBNonQCCoeffs *coeffs       /**<< OUTPUT, NQC coefficients */);


REAL8 XLALSimIMREOBGetNRSpinPeakDeltaT(
                                       INT4 l,           /**<< Mode l */
                                       INT4 m,           /**<< Mode m */
                                       REAL8 eta, /**<< Symmetric mass ratio */
                                       REAL8 a           /**<< Dimensionless spin */
);
REAL8 GetNRSpinPeakOmega( INT4 l, INT4  m, REAL8  eta, REAL8 a );

REAL8 GetNRSpinPeakOmegaDot( INT4 l, INT4  m, REAL8 eta, REAL8 a );

int  XLALSimIMREOBNonQCCorrection(
                                  COMPLEX16      *nqc,    /**<< OUTPUT, The NQC correction */
                                  REAL8Vector    *values, /**<< Dynamics r, phi, pr, pphi */
                                  const REAL8               omega,  /**<< Angular frequency */
                                  EOBNonQCCoeffs *coeffs  /**<< NQC coefficients */
);

INT4 XLALGenerateHybridWaveDerivatives (
                                        REAL8Vector *rwave,      /**<< OUTPUT, values of the waveform at comb points */
                                        REAL8Vector *dwave,      /**<< OUTPUT, 1st deriv of the waveform at comb points */
                                        REAL8Vector *ddwave,     /**<< OUTPUT, 2nd deriv of the waveform at comb points */
                                        REAL8Vector *timeVec,    /**<< Vector containing the time */
                                        REAL8Vector *wave,       /**<< Last part of inspiral waveform */
                                        REAL8Vector *matchrange, /**<< Times which determine the size of the comb */
                                        REAL8           dt,          /**<< Sample time step */
                                        REAL8           mass1,       /**<< First component mass (in Solar masses) */
                                        REAL8           mass2        /**<< Second component mass (in Solar masses) */
);

COMPLEX16 XLALSpinWeightedSphericalHarmonic(
                                            REAL8 theta,  /**< polar angle (rad) */
                                            REAL8 phi,    /**< azimuthal angle (rad) */
                                            int s,        /**< spin weight */
                                            int l,        /**< mode number l */
                                            int m         /**< mode number m */
);

int
CalculateThisMultipolePrefix(
                             COMPLEX16 *prefix, /**<< OUTPUT, Prefix value */
                             const REAL8 m1,    /**<< mass 1 */
                             const REAL8 m2,    /**<< mass 2 */
                             int l,      /**<< Mode l */
                             int m       /**<< Mode m */
);
int XLALIsREAL8FailNaN(REAL8 val);
int CartesianToSpherical(
                         REAL8 qSph[],        /**<< OUTPUT, position vector in spherical coordinates */
                         REAL8 pSph[],        /**<< OUTPUT, momentum vector in Cartesean coordinates */
                         const REAL8 qCart[], /**<< position vector in spherical coordinates */
                         const REAL8 pCart[]  /**<< momentum vector in Cartesean coordinates */
);
REAL8
XLALSimIMRSpinAlignedEOBNonKeplerCoeff(
                                       const REAL8           values[],   /**<< Dynamical variables */
                                       SpinEOBParams         *funcParams /**<< EOB parameters */
);
REAL8 GetNRSpinPeakAmplitude( INT4 l, INT4 m, REAL8 eta, REAL8 a );
REAL8 GetNRSpinPeakADDot( INT4 l, INT4 m, REAL8 eta, REAL8 a );

COMPLEX16 XLALCOMPLEX16Rect (REAL8 x, REAL8 y);
double MYlog2(double x);
double MYcbrt(double x);
double carg (COMPLEX16 z);
COMPLEX16 cexp (COMPLEX16 a);
COMPLEX16 CX16polar(double r,double phi);
COMPLEX16 MYcpow(COMPLEX16 a,UINT4 n);
double cabs(COMPLEX16 z);


INT4 XLALSimIMREOBHybridRingdownWave(
                                     REAL8Vector          *rdwave1,   /**<< OUTPUT, Real part of ringdown waveform */
                                     REAL8Vector          *rdwave2,   /**<< OUTPUT, Imag part of ringdown waveform */
                                     const REAL8           dt,        /**<< Sampling interval */
                                     const REAL8           mass1,     /**<< First component mass (in Solar masses) */
                                     const REAL8           mass2,     /**<< Second component mass (in Solar masses) */
                                     REAL8VectorSequence  *inspwave1, /**<< Values and derivs of real part inspiral waveform */
                                     REAL8VectorSequence  *inspwave2, /**<< Values and derivs of imag part inspiral waveform */
                                     COMPLEX16Vector      *modefreqs, /**<< Complex freqs of ringdown (scaled by total mass) */
                                     REAL8Vector          *matchrange /**<< Times which determine the comb of ringdown attachment */
);
INT4 XLALSimIMREOBHybridAttachRingdown(
                                       REAL8Vector *signal1,    /**<< OUTPUT, Real of inspiral waveform to which we attach ringdown */
                                       REAL8Vector *signal2,    /**<< OUTPUT, Imag of inspiral waveform to which we attach ringdown */
                                       const INT4   l,          /**<< Current mode l */
                                       const INT4   m,          /**<< Current mode m */
                                       const REAL8  dt,         /**<< Sample time step (in seconds) */
                                       const REAL8  mass1,      /**<< First component mass (in Solar masses) */
                                       const REAL8  mass2,      /**<< Second component mass (in Solar masses) */
                                       const REAL8  spin1x,     /**<<The spin of the first object; only needed for spin waveforms */
                                       const REAL8  spin1y,     /**<<The spin of the first object; only needed for spin waveforms */
                                       const REAL8  spin1z,     /**<<The spin of the first object; only needed for spin waveforms */
                                       const REAL8  spin2x,     /**<<The spin of the second object; only needed for spin waveforms */
                                       const REAL8  spin2y,     /**<<The spin of the second object; only needed for spin waveforms */
                                       const REAL8  spin2z,     /**<<The spin of the second object; only needed for spin waveforms */
                                       REAL8Vector *timeVec,    /**<< Vector containing the time values */
                                       REAL8Vector *matchrange /**<< Time values chosen as points for performing comb matching */
);
int XLALSimSEOBNRE(
                   REAL8TimeSeries **hplus,     /**<< OUTPUT, +-polarization waveform */
                   REAL8TimeSeries **hcross,    /**<< OUTPUT, x-polarization waveform */
                   const REAL8     phiC,        /**<< coalescence orbital phase (rad) */
                   REAL8           deltaT,      /**<< sampling time step */
                   const REAL8     m1SI,        /**<< mass-1 in SI unit */
                   const REAL8     m2SI,        /**<< mass-2 in SI unit */
                   const REAL8     fMin,        /**<< starting frequency (Hz) */
                   const REAL8     e0,          /**<< eccentricity at starting GW frequency (Hz) */
                   const REAL8     distance,           /**<< distance in SI unit */
                   const REAL8     inc,         /**<< inclination angle */
                   const REAL8     spin1z,      /**<< z-component of spin-1, dimensionless */
                   const REAL8     spin2z       /**<< z-component of spin-2, dimensionless */
);
void PNwaveformPRD544813rdotc_22mode(double *hr,double *hi,
                                     const double x1,const double x2,const double x3, // test particle position
                                     const double v1,const double v2,const double v3, // test particle velocity
                                     const double eta); // symmetric mass ratio
REAL8 XLALInspiralSpinFactorizedFlux_elip(
                      REAL8Vector           *values, /**< dynamical variables */
                      const REAL8   rphivalues[],   /**< dynamical varables, r, phi, pr, pphi */
                      const REAL8   drphivalues[],   /**< dot dynamical varables, r, phi, pr, pphi */
                      const REAL8           omega,   /**< orbital frequency */
                      SpinEOBParams         *ak,     /**< physical parameters */
                      const REAL8            H,      /**< real Hamiltonian */
                      const int             lMax    /**< upper limit of the summation over l */
                     );
double EOBPNflux(const double m1,const double m2,         // we assume m1>m2
          const double s1x,const double s1y,const double s1z,
          const double s2x,const double s2y,const double s2z,
      const double x1,const double x2,const double x3, // relative position
      const double v1,const double v2,const double v3, // relative velocity
      double *EOBPNfluxwoNQC);