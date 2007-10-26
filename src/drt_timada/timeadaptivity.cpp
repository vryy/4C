/*======================================================================*/
/*!
\file timeadaptivity.cpp

\brief Time step adaptivity

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

#include "timeadaptivity.H"
#include "iostream"


/*======================================================================*/
/*!
\brief Constructor (public)
\author bborn
\date 10/07
*/
TimeAdaptivity::TimeAdaptivity()
:  maxstepsize_(0.1),
   minstepsize_(1.0e-6),
   minsizeratio_(0.1),
   maxsizeratio_(1.5),
   saferatioscale_(0.9),
   errtol_(1.0e-3),
   errpow_(1)
{
   return;
}  // TimeAdaptivity::TimeAdaptivity()

/*======================================================================*/
/*!
\brief Constructor (public)
\author bborn
\date 10/07
*/
TimeAdaptivity::TimeAdaptivity
(
   double maxstepsize,
   double minstepsize,
   double minsizeratio,
   double maxsizeratio,
   double saferatioscale,
   double errtol,
   int errpow
)
:  maxstepsize_(maxstepsize),
   minstepsize_(minstepsize),
   minsizeratio_(minsizeratio),
   maxsizeratio_(maxsizeratio),
   saferatioscale_(saferatioscale),
   errtol_(errtol),
   errpow_(errpow)
{ 
   return;
}  // TimeAdaptivity::TimeAdaptivity(...)


/*======================================================================*/
/*!
\brief Destructor
\author bborn
\date 10/07
*/
TimeAdaptivity::~TimeAdaptivity()
{
   return;
}  // TimeAdaptivity::~TimeAdaptivity


/*======================================================================*/
/*!
\brief Out stream
\author
\date 10/07
*/
std::ostream& operator<<
(
   ostream& str, 
   const TimeAdaptivity& ta
)
{
   return str
      << "TimeAdaptivity" << endl
      << "Max step size = " << ta.maxstepsize_ << endl
      << "Min step size = " << ta.minstepsize_ << endl
      << "Min size ratio = " << ta.minsizeratio_ << endl
      << "Max size ratio = " << ta.maxsizeratio_ << endl
      << "Safety ratio scale = " << ta.saferatioscale_ << endl
      << "Error tolerance = " << ta.errtol_ << endl;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
