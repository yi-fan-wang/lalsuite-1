/************************************ <lalVerbatim file="UnitCompareCV">
Author: J. T. Whelan <jtwhelan@loyno.edu>
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{UnitCompare.c}}
\label{tools:ss:UnitCompare.c}

Compares two \texttt{LALUnit} structures, returning true if they are
equivalent, false otherwise.

\subsubsection*{Prototypes}
\input{UnitCompareCP}
\idx{LALUnitCompare()}

\subsubsection*{Description}

This function determines whether the units represented by
\texttt{*(input->unitOne)} and \texttt{*(input->unitTwo)} are the same (both
dimensionally and in the power-of-ten prefactor).  In this way,
programs and programmers can verify that quantities have the expected
units.
 
\subsubsection*{Algorithm}

The function first uses \texttt{LALUnitNormalize()} to bring both unit
structures into standard form, then compares the powers of ten and the
numerator and denominator of each exponent of a fundamental unit in
turn.

\subsubsection*{Uses}

\texttt{LALUnitNormalize()}

\subsubsection*{Notes}

\vfill{\footnotesize\input{UnitCompareCV}}

******************************************************* </lalLaTeX> */ 
#define TRUE 1
#define FALSE 0

#include <lal/LALStdlib.h>
#include <lal/Units.h>

NRCSID( UNITCOMPAREC, "$Id$" );


/* returns 1 if units are the same (after normalization) or 0 if not */
/* returns -1 if there is an error */
int XLALUnitCompare( const LALUnit *unit1, const LALUnit *unit2 )
{
  static const char *func = "XLALUnitCompare";
  LALUnit  unitOne, unitTwo;
  INT2 i;

  if ( ! unit1 || ! unit2 )
    XLAL_ERROR( func, XLAL_EFAULT );

  unitOne = *unit1;
  unitTwo = *unit2;

  /* normalize the units */
  if ( XLALUnitNormalize( &unitOne ) == XLAL_FAILURE )
    XLAL_ERROR( func, XLAL_EFUNC );
  if ( XLALUnitNormalize( &unitTwo ) == XLAL_FAILURE )
    XLAL_ERROR( func, XLAL_EFUNC );

  if (unitOne.powerOfTen != unitTwo.powerOfTen)
    return 0; /* false */

  for (i=0; i<LALNumUnits; ++i)
    if (
	unitOne.unitNumerator[i] != unitTwo.unitNumerator[i]
	|| unitOne.unitDenominatorMinusOne[i] 
	   != unitTwo.unitDenominatorMinusOne[i]
	)
      return 0; /* false */

  return 1; /* true */
}


/* <lalVerbatim file="UnitCompareCP"> */
void
LALUnitCompare (LALStatus *status, BOOLEAN *output, const LALUnitPair *input)
/* </lalVerbatim> */
     /* Compare two units structures, returning false if they differ */
{
  int code;

  INITSTATUS( status, "LALUnitCompare", UNITCOMPAREC );

  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );

  code = XLALUnitCompare( input->unitOne, input->unitTwo );
  if ( code == XLAL_FAILURE )
  {
    XLALClearErrno();
    ABORTXLAL( status );
  }

  *output = code ? TRUE : FALSE;

  RETURN(status);  
}
