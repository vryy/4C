/*----------------------------------------------------------------------*/
/*!
\file inpar_elch.H

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INPAR_ELCH_H
#define INPAR_ELCH_H

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for ELCH
  namespace ELCH
  {

    /// flag for the simulation of natural convection
    enum NatConv
    {
      natural_convection_no,
      natural_convection_substance,
      natural_convection_ion,
    };


  }  // namespace ELCH

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_ELCH_H
#endif  // #ifdef CCADISCRET
