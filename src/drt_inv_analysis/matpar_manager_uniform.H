/*----------------------------------------------------------------------*/
/*!
\file matpar_manager_unif.H
\brief Class to handle calls to material parameters from an optimization routine

<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef MATPAR_MANAGER_UNIFORM_H
#define MATPAR_MANAGER_UNIFORM_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "matpar_manager.H"
#include "../drt_inpar/inpar_material.H"
#include "Epetra_SerialDenseVector.h"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace STR
{

namespace INVANA

{
  class SMCParticle;



  class MatParManagerUniform  : public MatParManager
  {

  public:

  //! Constructor
  MatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret);

  //!destructor
  ~MatParManagerUniform(){};


  void ComputeParamsMultiVectorFromSMCParticlePosition(Teuchos::RCP<Epetra_MultiVector>& params, std::vector<double> my_global_params);

  };

}
}

/*----------------------------------------------------------------------*/
#endif
