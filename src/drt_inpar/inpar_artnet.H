/*----------------------------------------------------------------------*/
/*!
\file inpar_artnet.H

Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_ARTNET_H
#define INPAR_ARTNET_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/

// ToDo: move these enums to namespace INPAR::ARTNET etc.
//       is the typedef really needed?

/*!----------------------------------------------------------------------
\brief enum of arterial network dynamic types

<pre>                                                      ismail 11/08
This is the enumeration of all types of different integration schemes
</pre>

*-----------------------------------------------------------------------*/
typedef enum _ARTNET_DYNTYPE
{
  typ_tay_gal=0
} _ARTNET_DYNTYPE;


/*!----------------------------------------------------------------------
\brief enum of reduced dimensional airways dynamic types

<pre>                                                     ismail 01/10
This is the enumeration of all types of different integration schemes
</pre>

*-----------------------------------------------------------------------*/
typedef enum _RED_AIRWAYS_DYNTYPE
{
  typ_crank_nicolson,
  linear,
  nonlinear
} _RED_AIRWAYS_DYNTYPE;







namespace INPAR
{
  namespace ARTNET
  {
    /*!----------------------------------------------------------------------
    \brief enum of reduced dimensional relaxation type

    <pre>                                                     roth 10/13
    This is the enumeration of all types of different relaxation types
    </pre>

    *-----------------------------------------------------------------------*/
    enum RELAXTYPE_3D_0D
    {
      norelaxation,
      fixedrelaxation,
      Aitken,
      SD
    };

  }
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_ARTNET_H
