/*----------------------------------------------------------------------*/
/*!
\file matpar_parameter.cpp

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "matpar_parameter.H"
#include "matpar_material.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Parameter::Parameter(
  Teuchos::RCP<const MAT::PAR::Material> matdata
  )
: id_(matdata->Id()),
  type_(matdata->Type()),
  name_(matdata->Name())
{
}

/*----------------------------------------------------------------------*/
