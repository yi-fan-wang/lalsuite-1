/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Teviet Creighton, John Whelan
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/*********************************************************************
 * PREAMBLE                                                          *
 *********************************************************************/

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/FindRoot.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>

/* Define some constants used in this module. */
#define MAXORDER 6        /* Maximum number of N and PN terms */
#define BUFFSIZE 1024     /* Number of timesteps buffered */
#define ACCURACY (1.0e-6) /* Accuracy of root finder */
#define TWOTHIRDS (0.6666666667) /* 2/3 */
#define ONEMINUSEPS (0.99999)    /* Something close to 1 */

/* A macro to computing the (normalized) frequency.  It appears in
   many places, including in the main loop, and I don't want the
   overhead of a function call.  The following variables are required
   to be defined and set outside of the macro:

   REAL4 c0, c1, c2, c3, c4, c5;   PN frequency coefficients
   BOOLEAN b0, b1, b2, b3, b4, b5; whether to include each PN term

   The following variables must be defined outside the macro, but are
   set inside it:

   REAL4 x2, x3, x4, x5;  the input x raised to power 2, 3, 4, and 5 */
#define FREQ( f, x )                                                 \
do {                                                                 \
  x2 = (x)*(x);                                                      \
  x3 = x2*(x);                                                       \
  x4 = x3*(x);                                                       \
  x5 = x4*(x);                                                       \
  (f) = 0;                                                           \
  if ( b0 )                                                          \
    (f) += c0;                                                       \
  if ( b1 )                                                          \
    (f) += c1*(x);                                                   \
  if ( b2 )                                                          \
    (f) += c2*x2;                                                    \
  if ( b3 )                                                          \
    (f) += c3*x3;                                                    \
  if ( b4 )                                                          \
    (f) += c4*x4;                                                    \
  if ( b5 )                                                          \
    (f) += c5*x5;                                                    \
  (f) *= x3;                                                         \
} while (0)

/* Definition of a data structure used by FreqDiff() below. */
typedef struct tagFreqDiffParamStruc {
  REAL4 *c;   /* PN coefficients of frequency series */
  BOOLEAN *b; /* whether to include each PN term */
  REAL4 y0;   /* normalized frequency being sought */
} FreqDiffParamStruc;

/* A function to compute the difference between the current and
   requested normalized frequency, used by the root bisector. */
static void
FreqDiff( LALStatus *stat, REAL4 *y, REAL4 x, void *p )
{
  FreqDiffParamStruc *par;        /* *p cast to its proper type */
  REAL4 c0, c1, c2, c3, c4, c5;   /* PN frequency coefficients */
  BOOLEAN b0, b1, b2, b3, b4, b5; /* whether each order is nonzero */
  REAL4 x2, x3, x4, x5;           /* x^2, x^3, x^4, and x^5 */

  INITSTATUS(stat);
  ASSERT( p, stat, 1, "Null pointer" );

  /* Set constants used by FREQ() macro. */
  par = (FreqDiffParamStruc *)( p );
  c0 = par->c[0];
  c1 = par->c[1];
  c2 = par->c[2];
  c3 = par->c[3];
  c4 = par->c[4];
  c5 = par->c[5];
  b0 = par->b[0];
  b1 = par->b[1];
  b2 = par->b[2];
  b3 = par->b[3];
  b4 = par->b[4];
  b5 = par->b[5];

  /* Evaluate frequency and compare with reference. */
  FREQ( *y, x );
  *y -= par->y0;
  RETURN( stat );
}

/* Definition of a data buffer list for storing the waveform. */
typedef struct tagPPNInspiralBuffer {
  REAL4 a[2*BUFFSIZE];               /* amplitude data */
  REAL8 phi[BUFFSIZE];               /* phase data */
  REAL4 f[BUFFSIZE];                 /* frequency data */
  struct tagPPNInspiralBuffer *next; /* next buffer in list */
} PPNInspiralBuffer;

/* Definition of a macro to free the tail of said list, from a given
   node onward. */
#define FREELIST( node )                                             \
do {                                                                 \
  PPNInspiralBuffer *herePtr = (node);                               \
  while ( herePtr ) {                                                \
    PPNInspiralBuffer *lastPtr = herePtr;                            \
    herePtr = herePtr->next;                                         \
    LALFree( lastPtr );                                              \
  }                                                                  \
} while (0)


/*********************************************************************
 * MAIN FUNCTION                                                     *
 *********************************************************************/

/**
 * \author Creighton, T. D.
 *
 * \brief Computes a parametrized post-Newtonian inspiral waveform.
 *
 * ### Description ###
 *
 * This function computes an inspiral waveform using the parameters in
 * <tt>*params</tt>, storing the result in <tt>*output</tt>.
 *
 * In the <tt>*params</tt> structure, the routine uses all the "input"
 * fields specified in \ref GeneratePPNInspiral.h, and sets all of the
 * "output" fields.  If <tt>params->ppn=NULL</tt>, a normal
 * post\f${}^2\f$-Newtonian waveform is generated; i.e.\ \f$p_0=1\f$, \f$p_1=0\f$,
 * \f$p_2=1\f$, \f$p_3=1\f$, \f$p_4=1\f$, \f$p_{5+}=0\f$.
 *
 * In the <tt>*output</tt> structure, the field <tt>output->h</tt> is
 * ignored, but all other pointer fields must be set to \c NULL.  The
 * function will create and allocate space for <tt>output->a</tt>,
 * <tt>output->f</tt>, and <tt>output->phi</tt> as necessary.  The
 * <tt>output->shift</tt> field will remain set to \c NULL, as it is
 * not required to describe a nonprecessing binary.
 *
 * ### Algorithm ###
 *
 * This function is a fairly straightforward calculation of
 * \eqref{eq_ppn_freq}--\eqref{eq_ppn_across} in
 * \ref GeneratePPNInspiral.h.  However, there are some nontrivial
 * issues involved, which are discussed in some depth in Secs.\ 6.4, 6.6,
 * and 6.9.2 of \cite GRASP2000 .  What follows is a brief
 * discussion of these issues and how this routine deals with them.
 *
 * <tt>Computing the start time:</tt>
 *
 * When building a waveform for data analysis, one would generally like
 * to start the waveform at some well-defined frequency where it first
 * enters the band of interest; one then defines the start time of the
 * integration by inverting \eqref{eq_ppn_freq} if
 * \ref GeneratePPNInspiral.h.  The current algorithm follows this
 * standard approach by requiring the calling routine to specify
 * <tt>params->fStartIn</tt>, which is then inverted to find
 * \f$\Theta_\mathrm{start}\f$.  This inversion is in fact the most
 * algorithmically complicated part of the routine, so we will discuss it
 * in depth.
 *
 * To help clarify the problem, let us rewrite the equation in
 * dimensionless parameters \f$y=8\pi fT_\odot m_\mathrm{tot}/M_\odot\f$ and
 * \f$x=\Theta^{-1/8}\f$:
 * \f{equation}{
 * \label{eq_ppn_fnorm}
 * y = \sum_{k=0}^{5} C_k x^{k+3} \; ,
 * \f}
 * where:
 * \f{eqnarray}{
 * C_0 & = & p_0 \;,\\
 * C_1 & = & p_1 \;,\\
 * C_2 & = & p_2\left(\frac{743}{2688}+\frac{11}{32}\eta\right) \;,\\
 * C_3 & = & p_3\frac{3\pi}{10} \;,\\
 * C_4 & = & p_4\left(\frac{1855099}{14450688}+\frac{56975}{258048}\eta+
 * \frac{371}{2048}\eta^2\right) \;,\\
 * C_5 & = & p_5\left(\frac{7729}{21504}+\frac{3}{256}\eta\right)\pi \;.
 * \f}
 * We note that \f$x\f$ is a time parameter mapping the range
 * \f$t=-\infty\rightarrow t_c\f$ to \f$x=0\rightarrow\infty\f$.
 *
 * In a normal post-Newtonian expansion it is possible to characterize
 * the general behaviour of this equation quite accurately, since the
 * values of \f$p_k\f$ are known and since \f$\eta\f$ varies only over the range
 * \f$[0,1/4]\f$.  In a parametrized post-Newtonian expansion, however, even
 * the relative orders of magnitude of the coefficients can vary
 * significantly, making a robust generic root finder impractical.
 * However, we are saved by the fact that we can restrict our search to
 * domains where the post-Newtonian expansion is a valid approximation.
 * We define the post-Newtonian expansion \e not to be valid if
 * \e any of the following conditions occur:
 * <ol>
 * <li> A higher-order term in the frequency expansion becomes larger in
 * magnitude than the leading (lowest-order nonzero) term.</li>
 * <li> The inferred orbital radius, approximated by
 * \f$r\sim4m_\mathrm{tot}\Theta^{1/4}\f$, drops below \f$2m_\mathrm{tot}\f$;
 * i.e.\ \f$\Theta<1/16\f$ or \f$x>\sqrt{2}\f$.</li>
 * <li> The frequency evolution becomes non-monotonic.</li>
 * </ol>
 * We can further require as a matter of convention that the lowest-order
 * nonzero coefficient in the frequency expansion be positive; this is
 * simply a sign convention stating that the frequency of a system be
 * positive at large radii.
 *
 * The first two conditions above allow us to set firm limits on the
 * range of the initial \f$x_\mathrm{start}\f$.  Let \f$C_j\f$ be the
 * lowest-order nonzero coefficient; then for every nonzero \f$C_{k>j}\f$
 * we can define a point \f$x_k=|C_j/C_k|^{1/(k-j)}\f$ where that term
 * exceeds the leading-order term in magnitude.  We can therefore limit
 * the range of \f$x\f$ to values less than \f$x_\mathrm{max}\f$, which is the
 * minimum of \f$\sqrt{2}\f$ and all \f$x_k\f$.  We note that even if we were
 * to extend the post-Newtonian expansion in \eqref{eq_ppn_fnorm} to an
 * infinite number of terms, this definition of \f$x_\mathrm{max}\f$ implies
 * that the frequency is guaranteed to be monotonic up to
 * \f$x_\mathrm{max}(5-\sqrt{7})/6\f$, and positive up to \f$x_\mathrm{max}/2\f$.
 * Thus we can confidently begin our search for \f$x_\mathrm{start}\f$ in the
 * domain \f$(0,0.39x_\mathrm{max})\f$, where the leading-order term
 * dominates, and end it if we ever exceed \f$x_\mathrm{max}\f$.
 *
 * We therefore bracket our value of \f$x_\mathrm{start}\f$ as follows: We
 * start with an initial guess
 * \f$x_\mathrm{guess}=(y_\mathrm{start}/C_j)^{1/(j+3)}\f$, or
 * \f$0.39x_\mathrm{max}\f$, whichever is less.  If
 * \f$y(x_\mathrm{guess})>y_\mathrm{start}\f$, we iteratively decrease \f$x\f$ by
 * factors of 0.95 until \f$y(x)<y_\mathrm{start}\f$; this is guaranteed to
 * occur within a few iterations, since we are moving into a regime where
 * the leading-order behaviour dominates more and more.  If
 * \f$y(x_\mathrm{guess})<y_\mathrm{start}\f$, we iteratively increase \f$x\f$ by
 * factors of 1.05 until \f$y(x)>y_\mathrm{start}\f$, or until
 * \f$x>x_\mathrm{max}\f$; this is also guaranteed to occur quickly because,
 * in the worst case, it only takes about 20 iterations to step from
 * \f$0.39x_\mathrm{max}\f$ to \f$x_\mathrm{max}\f$, and if \f$x_\mathrm{guess}\f$
 * were much lower than \f$0.39x_\mathrm{max}\f$ it would have been a pretty
 * good guess to begin with.  If at any point while increasing \f$x\f$ we
 * find that \f$y\f$ is decreasing, we determine that the starting frequency
 * is already in a regime where the post-Newtonian approximation is
 * invalid, and we return an error.  Otherwise, once we have bracketed
 * the value of \f$x_\mathrm{start}\f$, we use <tt>LALSBisectionFindRoot()</tt>
 * to pin down the value to an accuracy of a part in \f$10^6\f$.
 *
 * <tt>Computing the phase and amplitudes:</tt>
 *
 * Once we have \f$x_\mathrm{start}\f$, we can find \f$\Theta_\mathrm{start}\f$,
 * and begin incrementing it; at each timestep we compute \f$x\f$ and hence
 * \f$f\f$, \f$\phi\f$, \f$A_+\f$, and \f$A_\times\f$ according to
 * \eqref{eq_ppn_freq}, \eqref{eq_ppn_phi}, \eqref{eq_ppn_aplus},
 * and \eqref{eq_ppn_across}.  The routine progressively creates a list of
 * length-1024 arrays and fills them.  The process stops when any of the
 * following occurs:
 * <ol>
 * <li> The frequency exceeds the requested termination frequency.</li>
 * <li> The number of steps reaches the suggested maximum length in
 * <tt>*params</tt>.</li>
 * <li> The frequency is no longer increasing.</li>
 * <li> The parameter \f$x>x_\mathrm{max}\f$.</li>
 * <li> We run out of memory.</li>
 * </ol>
 * In the last case an error is returned; otherwise the waveform is
 * deemed "complete".  Output arrays are created of the appropriate
 * length and are filled with the data.
 *
 * Internally, the routine keeps a list of all coefficients, as well as a
 * list of booleans indicating which terms are nonzero.  The latter
 * allows the code to avoid lengthy floating-point operations (especially
 * the logarithm in the post\f${}^{5/2}\f$-Newtonian phase term) when these
 * are not required.
 *
 * When generating the waveform, we note that the sign of \f$\dot{f}\f$ is
 * the same as the sign of \f$y'/x^2 = \sum (k+3)C_k x^k\f$, and use this
 * series to test whether the frequency has stopped increasing.  (The
 * reason is that for waveforms far from coalescence \f$f\f$ is nearly
 * constant: numerical errors can cause positive <em>and negative</em>
 * fluctuations \f$\Delta f\f$ bewteen timesteps.  The analytic formulae for
 * \f$\dot{f}\f$ or \f$y'\f$ are less susceptible to this.)  The coefficients
 * \f$(k+3)C_k\f$ are also precomputed for added efficiency.
 *
 * <tt>Warnings and suggestions:</tt>
 *
 * If no post-Newtonian parameters are provided (i.e.\
 * <tt>params->ppn=NULL</tt>), we generate a post\f${}^2\f$-Newtonian waveform,
 * \e not a post\f${}^{5/2}\f$-Newtonian waveform.  This is done not only
 * for computationally efficiency, but also because the accuracy and
 * reliability of the post\f${}^{5/2}\f$-Newtonian waveform is actually
 * worse.  You can of course specify a post\f${}^{5/2}\f$-Newtonian waveform
 * with an appropriate assignment of <tt>params->ppn</tt>, but you do so at
 * your own risk!
 *
 * This routine also performs no sanity checking on the requested
 * sampling interval \f$\Delta t=\f$<tt>params->deltaT</tt>, because this
 * depends very much on how one intends to use the generated waveform.
 * If you plan to generate actual wave functions \f$h_{+,\times}(t)\f$ at the
 * same sample rate, then you will generally want a sampling interval
 * \f$\Delta t<1/2f_\mathrm{max}\f$; you can enforce this by specifying a
 * suitable <tt>params->fStopIn</tt>.
 *
 * However, many routines (such as those in \ref SimulateCoherentGW.h)
 * generate actual wave functions by linear interpolation of the
 * amplitude and phase data, which then need only be sampled on
 * timescales \f$\sim\dot{f}^{-1/2}\f$ rather than \f$\sim f^{-1}\f$.  More
 * precisely, we would like our interpolated phase to differ from the
 * actual phase by no more than some specified amount, say \f$\pi/2\f$
 * radians.  The largest deviation from linear phase evolution will
 * typically be on the order of \f$\Delta\phi\approx(1/2)\ddot{\phi}(\Delta
 * t/2)^2\approx(\pi/4)\Delta f\Delta t\f$, where \f$\Delta f\f$ is the
 * frequency shift over the timestep.  Thus in general we would like to
 * have
 * \f[
 * \Delta f \Delta t \lesssim 2
 * \f]
 * for our linear interpolation to be valid.  This routine helps out by
 * setting the output parameter field <tt>params->dfdt</tt> equal to the
 * maximum value of \f$\Delta f\Delta t\f$ encountered during the
 * integration.
 */
void
LALGeneratePPNInspiral( LALStatus     *stat,	/**< UNDOCUMENTED */
			CoherentGW    *output,	/**< UNDOCUMENTED */
		        PPNParamStruc *params 	/**< UNDOCUMENTED */
                        )
{

  /* System-derived constants. */
  BOOLEAN b0, b1, b2, b3, b4, b5; /* whether each order is nonzero */
  BOOLEAN b[MAXORDER];            /* vector of above coefficients */
  REAL4 c0, c1, c2, c3, c4, c5;   /* PN frequency coefficients */
  REAL4 c[MAXORDER];              /* vector of above coefficients */
  REAL4 d0, d1, d2, d3, d4, d5;   /* PN phase coefficients */
  REAL4 e0, e1, e2, e3, e4, e5;   /* PN dy/dx coefficients */
  REAL4 p[MAXORDER];              /* PN parameter values */
  REAL4 mTot, mu;      /* total mass and reduced mass */
  REAL4 eta, etaInv;   /* mass ratio and its inverse */
  REAL4 phiC;          /* phase at coalescence */
  REAL4 cosI;          /* cosine of system inclination */
  REAL4 fFac;          /* SI normalization for f and t */
  REAL4 f2aFac;        /* factor multiplying f in amplitude function */
  REAL4 apFac, acFac;  /* extra factor in plus and cross amplitudes */

  /* Integration parameters. */
  UINT4 i;         /* index over PN terms */
  UINT4 j;         /* index of leading nonzero PN term */
  UINT4 n, nMax;   /* index over timesteps, and its maximum + 1 */
  UINT4 nNext;     /* index where next buffer starts */
  REAL8 t, t0, dt; /* dimensionless time, start time, and increment */
  REAL4 tStop = 0.0625;  /* time when orbit reaches minimum radius */
  REAL4 x, xStart, xMax; /* x = t^(-1/8), and its maximum range */
  REAL4 y, yStart, yMax; /* normalized frequency and its range */
  REAL4 yOld, dyMax;     /* previous timestep y, and maximum y - yOld */
  REAL4 x2, x3, x4, x5;  /* x^2, x^3, x^4, and x^5 */
  REAL4 *a, *f; /* pointers to generated amplitude and frequency data */
  REAL8 *phi;   /* pointer to generated phase data */
  PPNInspiralBuffer *head, *here; /* pointers to buffered data */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /*******************************************************************
   * CHECK INPUT PARAMETERS                                          *
   *******************************************************************/

  /* Dumb initialization to shut gcc up. */
  head = here = NULL;
  b0 = b1 = b2 = b3 = b4 = b5 = 0;
  c0 = c1 = c2 = c3 = c4 = c5 = d0 = d1 = d2 = d3 = d4 = d5 = 0.0;

  /* Make sure parameter and output structures exist. */
  ASSERT( params, stat, GENERATEPPNINSPIRALH_ENUL,
	  GENERATEPPNINSPIRALH_MSGENUL );
  ASSERT( output, stat, GENERATEPPNINSPIRALH_ENUL,
	  GENERATEPPNINSPIRALH_MSGENUL );

  /* Make sure output fields don't exist. */
  ASSERT( !( output->a ), stat, GENERATEPPNINSPIRALH_EOUT,
	  GENERATEPPNINSPIRALH_MSGEOUT );
  ASSERT( !( output->f ), stat, GENERATEPPNINSPIRALH_EOUT,
	  GENERATEPPNINSPIRALH_MSGEOUT );
  ASSERT( !( output->phi ), stat, GENERATEPPNINSPIRALH_EOUT,
	  GENERATEPPNINSPIRALH_MSGEOUT );
  ASSERT( !( output->shift ), stat, GENERATEPPNINSPIRALH_EOUT,
	  GENERATEPPNINSPIRALH_MSGEOUT );

  /* Get PN parameters, if they are specified; otherwise use
     post2-Newtonian. */
  if ( params->ppn ) {
    ASSERT( params->ppn->data, stat, GENERATEPPNINSPIRALH_ENUL,
	    GENERATEPPNINSPIRALH_MSGENUL );
    j = params->ppn->length;
    if ( j > MAXORDER )
      j = MAXORDER;
    for ( i = 0; i < j; i++ )
      p[i] = params->ppn->data[i];
    for ( ; i < MAXORDER; i++ )
      p[i] = 0.0;
  } else {
    p[0] = 1.0;
    p[1] = 0.0;
    p[2] = 1.0;
    p[3] = 1.0;
    p[4] = 1.0;
    for ( i = 5; i < MAXORDER; i++ )
      p[i] = 0.0;
  }

  /*******************************************************************
   * COMPUTE SYSTEM PARAMETERS                                       *
   *******************************************************************/

  /* Compute parameters of the system. */
  mTot = params->mTot;
  ASSERT( mTot != 0.0, stat, GENERATEPPNINSPIRALH_EMBAD,
	  GENERATEPPNINSPIRALH_MSGEMBAD );
  eta = params->eta;
  ASSERT( eta != 0.0, stat, GENERATEPPNINSPIRALH_EMBAD,
	  GENERATEPPNINSPIRALH_MSGEMBAD );
  etaInv = 2.0 / eta;
  mu = eta*mTot;
  cosI = cos( params->inc );
  phiC = params->phi;

  /* Compute frequency, phase, and amplitude factors. */
  fFac = 1.0 / ( 4.0*LAL_TWOPI*LAL_MTSUN_SI*mTot );
  dt = -params->deltaT * eta / ( 5.0*LAL_MTSUN_SI*mTot );
  ASSERT( dt < 0.0, stat, GENERATEPPNINSPIRALH_ETBAD,
	  GENERATEPPNINSPIRALH_MSGETBAD );
  f2aFac = LAL_PI*LAL_MTSUN_SI*mTot*fFac;
  ASSERT( params->d != 0.0, stat, GENERATEPPNINSPIRALH_EDBAD,
	  GENERATEPPNINSPIRALH_MSGEDBAD );
  apFac = acFac = -2.0*mu*LAL_MRSUN_SI/params->d;
  apFac *= 1.0 + cosI*cosI;
  acFac *= 2.0*cosI;

  /* Compute PN expansion coefficients. */
  c0 = c[0] = p[0];
  c1 = c[1] = p[1];
  c2 = c[2] = p[2]*( 743.0/2688.0 + eta*11.0/32.0 );
  c3 = c[3] = -p[3]*( 3.0*LAL_PI/10.0 );
  c4 = c[4] = p[4]*( 1855099.0/14450688.0 + eta*56975.0/258048.0 +
		     eta*eta*371.0/2048.0 );
  c5 = c[5] = -p[5]*( 7729.0/21504.0 + eta*3.0/256.0 )*LAL_PI;

  /* Compute expansion coefficients for series in phi and dy/dx. */
  d0 = c0;
  d1 = c1*5.0/4.0;
  d2 = c2*5.0/3.0;
  d3 = c3*5.0/2.0;
  d4 = c4*5.0;
  d5 = c5*5.0/8.0;
  e0 = c0*3.0;
  e1 = c1*4.0;
  e2 = c2*5.0;
  e3 = c3*6.0;
  e4 = c4*7.0;
  e5 = c5*8.0;

  /* Use Boolean variables to exclude terms that are zero. */
  b0 = b[0] = ( c0 == 0.0 ? 0 : 1 );
  b1 = b[1] = ( c1 == 0.0 ? 0 : 1 );
  b2 = b[2] = ( c2 == 0.0 ? 0 : 1 );
  b3 = b[3] = ( c3 == 0.0 ? 0 : 1 );
  b4 = b[4] = ( c4 == 0.0 ? 0 : 1 );
  b5 = b[5] = ( c5 == 0.0 ? 0 : 1 );

  /* Find the leading-order frequency term. */
  for ( j = 0; ( j < MAXORDER ) && ( b[j] == 0 ); j++ )
    ;
  if ( j == MAXORDER ) {
    ABORT( stat, GENERATEPPNINSPIRALH_EPBAD,
	   GENERATEPPNINSPIRALH_MSGEPBAD );
  }

  /*******************************************************************
   * COMPUTE START TIME                                              *
   *******************************************************************/

  /* First, find the normalized start frequency, and the best guess as
     to the start time from the leading-order term.  We require the
     frequency to be increasing. */
  yStart = params->fStartIn / fFac;
  if ( params->fStopIn == 0.0 )
    yMax = LAL_REAL4_MAX;
  else {
    ASSERT( fabs( params->fStopIn ) > params->fStartIn, stat,
	    GENERATEPPNINSPIRALH_EFBAD, GENERATEPPNINSPIRALH_MSGEFBAD );
    yMax = fabs( params->fStopIn ) / fFac;
  }
  if ( ( c[j]*fFac < 0.0 ) || ( yStart < 0.0 ) || ( yMax < 0.0 ) ) {
    ABORT( stat, GENERATEPPNINSPIRALH_EPBAD,
	   GENERATEPPNINSPIRALH_MSGEPBAD );
  }
  xStart = pow( yStart/c[j], 1.0/( j + 3.0 ) );
  xMax = LAL_SQRT2;

  /* The above is exact if the leading-order term is the only one in
     the expansion.  Check to see if there are any other terms. */
  for ( i = j + 1; ( i < MAXORDER ) && ( b[i] == 0 ); i++ )
    ;
  if ( i < MAXORDER ) {
    /* There are other terms, so we have to use bisection to find the
       start time. */
    REAL4 xLow, xHigh; /* ultimately these will bracket xStart */
    REAL4 yLow, yHigh; /* the normalized frequency at these times */

    /* If necessary, revise the estimate of the cutoff where we know
       the PN approximation goes bad, and revise our initial guess to
       lie well within the valid regime. */
    for ( i = j + 1; i < MAXORDER; i++ )
      if ( b[i] != 0 ) {
	x = pow( fabs( c[j]/c[i] ), 1.0/(REAL4)( i - j ) );
	if ( x < xMax )
	  xMax = x;
      }
    if ( xStart > 0.39*xMax )
      xStart = 0.39*xMax;

    /* If we are ignoring PN breakdown, adjust xMax (so that it won't
       interfere with the start time search) and tStop. */
    if ( params->fStopIn < 0.0 ) {
      xMax = LAL_REAL4_MAX;
      tStop = 0.0;
    }

    /* If our frequency is too high, step backwards and/or forwards
       until we have bracketed the correct frequency. */
    xLow = xHigh = xStart;
    FREQ( yHigh, xStart );
    yLow = yHigh;
    while ( yLow > yStart ) {
      xHigh = xLow;
      yHigh = yLow;
      xLow *= 0.95;
      FREQ( yLow, xLow );
    }
    while ( yHigh < yStart ) {
      xLow = xHigh;
      yLow = yHigh;
      xHigh *= 1.05;
      FREQ( yHigh, xHigh );
      /* Check for PN breakdown. */
      if ( ( yHigh < yLow ) || ( xHigh > xMax ) ) {
	ABORT( stat, GENERATEPPNINSPIRALH_EFBAD,
	       GENERATEPPNINSPIRALH_MSGEFBAD );
      }
    }

    /* We may have gotten lucky and nailed the frequency right on.
       Otherwise, find xStart by root bisection. */
    if ( yLow == yStart )
      xStart = xLow;
    else if ( yHigh == yStart )
      xStart = xHigh;
    else {
      SFindRootIn in;
      FreqDiffParamStruc par;
      in.xmax = xHigh;
      in.xmin = xLow;
      in.xacc = ACCURACY;
      in.function = FreqDiff;
      par.c = c;
      par.b = b;
      par.y0 = yStart;
      TRY( LALSBisectionFindRoot( stat->statusPtr, &xStart, &in,
				  (void *)( &par ) ), stat );
    }
  }

  /* If we are ignoring PN breakdown, adjust xMax and tStop, if they
     haven't been adjusted already. */
  else if ( params->fStopIn < 0.0 ) {
    xMax = LAL_REAL4_MAX;
    tStop = 0.0;
  }

  /* Compute initial dimensionless time, record actual initial
     frequency (in case it is different), and record dimensional
     time-to-coalescence. */
  t0 = pow( xStart, -8.0 );
  FREQ( yStart, xStart );
  if ( yStart >= yMax ) {
    ABORT( stat, GENERATEPPNINSPIRALH_EFBAD,
	   GENERATEPPNINSPIRALH_MSGEFBAD );
  }
  params->fStart = yStart*fFac;
  params->tc = t0 * ( 5.0*LAL_MTSUN_SI*mTot ) / eta;

  /*******************************************************************
   * GENERATE WAVEFORM                                               *
   *******************************************************************/

  /* Set up data pointers and storage. */
  here = head = (PPNInspiralBuffer *)
    LALMalloc( sizeof(PPNInspiralBuffer) );
  if ( !here ) {
    ABORT( stat, GENERATEPPNINSPIRALH_EMEM,
	   GENERATEPPNINSPIRALH_MSGEMEM );
  }
  here->next = NULL;
  a = here->a;
  f = here->f;
  phi = here->phi;
  nMax = (UINT4)( -1 );
  if ( params->lengthIn > 0 )
    nMax = params->lengthIn;
  nNext = BUFFSIZE;
  if ( nNext > nMax )
    nNext = nMax;

  /* Start integrating!  Inner loop exits each time a new buffer is
     required.  Outer loop has no explicit test; when a termination
     condition is met, we jump directly from the inner loop using a
     goto statement.  All goto statements jump to the terminate: label
     at the end of the outer loop. */
  n = 0;
  t = t0;
  dyMax = 0.0;
  y = yOld = 0.0;
  x = xStart;
  while ( 1 ) {
    while ( n < nNext ) {
      REAL4 f2a; /* value inside 2/3 power in amplitude functions */
      REAL4 phase = 0.0; /* wave phase excluding overall constants */
      REAL4 dydx2 = 0.0; /* dy/dx divided by x^2 */

      /* Check if we're still in a valid PN regime. */
      if ( x > xMax ) {
	params->termCode = GENERATEPPNINSPIRALH_EPNFAIL;
	params->termDescription = GENERATEPPNINSPIRALH_MSGEPNFAIL;
	goto terminate;
      }

      /* Compute the normalized frequency.  This also computes the
         variables x2, x3, x4, and x5, which are used later. */
      FREQ( y, x );
      if ( y > yMax ) {
	params->termCode = GENERATEPPNINSPIRALH_EFSTOP;
	params->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;
	goto terminate;
      }

      /* Check that frequency is still increasing. */
      if ( b0 )
	dydx2 += e0;
      if ( b1 )
	dydx2 += e1*x;
      if ( b2 )
	dydx2 += e2*x2;
      if ( b3 )
	dydx2 += e3*x3;
      if ( b4 )
	dydx2 += e4*x4;
      if ( b5 )
	dydx2 += e5*x5;
      if ( dydx2 < 0.0 ) {
	params->termCode = GENERATEPPNINSPIRALH_EFNOTMON;
	params->termDescription = GENERATEPPNINSPIRALH_MSGEFNOTMON;
	goto terminate;
      }
      if ( y - yOld > dyMax )
	dyMax = y - yOld;
      *(f++) = fFac*y;

      /* Compute the amplitude from the frequency. */
      f2a = pow( f2aFac*y, TWOTHIRDS );
      *(a++) = apFac*f2a;
      *(a++) = acFac*f2a;

      /* Compute the phase. */
      if ( b0 )
	phase += d0;
      if ( b1 )
	phase += d1*x;
      if ( b2 )
	phase += d2*x2;
      if ( b3 )
	phase += d3*x3;
      if ( b4 )
	phase += d4*x4;
      if ( b5 )
	phase += d5*log(t)*x5;
      phase *= t*x3*etaInv;
      *(phi++) = phiC - phase;

      /* Increment the timestep. */
      n++;
      t = t0 + n*dt;
      yOld = y;
      if ( t <= tStop ) {
	params->termCode = GENERATEPPNINSPIRALH_ERTOOSMALL;
	params->termDescription = GENERATEPPNINSPIRALH_MSGERTOOSMALL;
	goto terminate;
      }
      x = pow( t, -0.125 );
    }

    /* We've either filled the buffer or we've exceeded the maximum
       length.  If the latter, we're done! */
    if ( n >= nMax ) {
      params->termCode = GENERATEPPNINSPIRALH_ELENGTH;
      params->termDescription = GENERATEPPNINSPIRALH_MSGELENGTH;
      goto terminate;
    }

    /* Otherwise, allocate the next buffer. */
    here->next =
      (PPNInspiralBuffer *)LALMalloc( sizeof(PPNInspiralBuffer) );
    here = here->next;
    if ( !here ) {
      FREELIST( head );
      ABORT( stat, GENERATEPPNINSPIRALH_EMEM,
	     GENERATEPPNINSPIRALH_MSGEMEM );
    }
    here->next = NULL;
    a = here->a;
    f = here->f;
    phi = here->phi;
    nNext += BUFFSIZE;
    if ( nNext > nMax )
      nNext = nMax;
  }

  /*******************************************************************
   * CLEANUP                                                         *
   *******************************************************************/

  /* The above loop only exits by triggering one of the termination
     conditions, which jumps to the following point for cleanup and
     return. */
 terminate:

  /* First, set remaining output parameter fields. */
  params->dfdt = dyMax*fFac*params->deltaT;
  params->fStop = yOld*fFac;
  params->length = n;

  /* Allocate the output structures. */
  if ( ( output->a = (REAL4TimeVectorSeries *)
	 LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
    FREELIST( head );
    ABORT( stat, GENERATEPPNINSPIRALH_EMEM,
	   GENERATEPPNINSPIRALH_MSGEMEM );
  }
  memset( output->a, 0, sizeof(REAL4TimeVectorSeries) );
  if ( ( output->f = (REAL4TimeSeries *)
	 LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
    FREELIST( head );
    LALFree( output->a ); output->a = NULL;
    ABORT( stat, GENERATEPPNINSPIRALH_EMEM,
	   GENERATEPPNINSPIRALH_MSGEMEM );
  }
  memset( output->f, 0, sizeof(REAL4TimeSeries) );
  if ( ( output->phi = (REAL8TimeSeries *)
	 LALMalloc( sizeof(REAL8TimeSeries) ) ) == NULL ) {
    FREELIST( head );
    LALFree( output->a ); output->a = NULL;
    LALFree( output->f ); output->f = NULL;
    ABORT( stat, GENERATEPPNINSPIRALH_EMEM,
	   GENERATEPPNINSPIRALH_MSGEMEM );
  }
  memset( output->phi, 0, sizeof(REAL8TimeSeries) );

  /* Allocate the output data fields. */
  {
    CreateVectorSequenceIn in;
    in.length = n;
    in.vectorLength = 2;
    LALSCreateVectorSequence( stat->statusPtr, &( output->a->data ), &in );
    BEGINFAIL( stat ) {
      FREELIST( head );
      LALFree( output->a );   output->a = NULL;
      LALFree( output->f );   output->f = NULL;
      LALFree( output->phi ); output->phi = NULL;
    } ENDFAIL( stat );
    LALSCreateVector( stat->statusPtr, &( output->f->data ), n );
    BEGINFAIL( stat ) {
      TRY( LALSDestroyVectorSequence( stat->statusPtr, &( output->a->data ) ),
	   stat );
      FREELIST( head );
      LALFree( output->a );   output->a = NULL;
      LALFree( output->f );   output->f = NULL;
      LALFree( output->phi ); output->phi = NULL;
    } ENDFAIL( stat );
    LALDCreateVector( stat->statusPtr, &( output->phi->data ), n );
    BEGINFAIL( stat ) {
      TRY( LALSDestroyVectorSequence( stat->statusPtr, &( output->a->data ) ),
	   stat );
      TRY( LALSDestroyVector( stat->statusPtr, &( output->f->data ) ),
	   stat );
      FREELIST( head );
      LALFree( output->a );   output->a = NULL;
      LALFree( output->f );   output->f = NULL;
      LALFree( output->phi ); output->phi = NULL;
    } ENDFAIL( stat );
  }

  /* Structures have been successfully allocated; now fill them.  We
     deallocate the list as we go along. */
  output->position = params->position;
  output->psi = params->psi;
  output->a->epoch = output->f->epoch = output->phi->epoch
    = params->epoch;
  output->a->deltaT = output->f->deltaT = output->phi->deltaT
    = params->deltaT;
  output->a->sampleUnits = lalStrainUnit;
  output->f->sampleUnits = lalHertzUnit;
  output->phi->sampleUnits = lalDimensionlessUnit;
  snprintf( output->a->name, LALNameLength, "PPN inspiral amplitudes" );
  snprintf( output->f->name, LALNameLength, "PPN inspiral frequency" );
  snprintf( output->phi->name, LALNameLength, "PPN inspiral phase" );
  a = output->a->data->data;
  f = output->f->data->data;
  phi = output->phi->data->data;
  here = head;
  while ( here && ( n > 0 ) ) {
    PPNInspiralBuffer *last = here;
    UINT4 nCopy = BUFFSIZE;
    if ( nCopy > n )
      nCopy = n;
    memcpy( a, here->a, 2*nCopy*sizeof(REAL4) );
    memcpy( f, here->f, nCopy*sizeof(REAL4) );
    memcpy( phi, here->phi, nCopy*sizeof(REAL8) );
    a += 2*nCopy;
    f += nCopy;
    phi += nCopy;
    n -= nCopy;
    here = here->next;
    LALFree( last );
  }

  /* This shouldn't happen, but free any extra buffers in the
     list. */
  FREELIST( here );

  /* Everything's been stored and cleaned up, so there's nothing left
     to do but quit! */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
