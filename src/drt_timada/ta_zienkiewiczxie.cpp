/*======================================================================*/
/*!
\file ta_zienkiewczxie.cpp
\brief Generalized Alpha time integration for structural problems

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "ta_zienkiewiczxie.H"
#include "iostream"


/*======================================================================*/
/*!
\brief Constructor (public)
\author bborn
\date 10/07
*/
ZienkiewiczXie::ZienkiewiczXie() : TimeAdaptivity()
{
   return;
}

/*======================================================================*/
/*!
\brief Constructor with parameters (public)
\author bborn
\date 10/07
*/
ZienkiewiczXie::ZienkiewiczXie
(
   double maxstepsize,
   double minstepsize,
   double minsizeratio,
   double maxsizeratio,
   double saferatioscale,
   double errtol,
   int errpow
)
:  TimeAdaptivity
   (
      maxstepsize,
      minstepsize,
      minsizeratio,
      maxsizeratio,
      saferatioscale,
      errtol,
      errpow
   )
{
   return;
}

/*======================================================================*/
/*!
\brief Destructor
\author bborn
\date 10/07
*/
ZienkiewiczXie::~ZienkiewiczXie()
{
   return;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
