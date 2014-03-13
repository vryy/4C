/*----------------------------------------------------------------------*/
/*!
\file matpar_manager_uniform.H
\brief Class to handle calls to material parameters from an optimization routine

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-10361
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef MATPAR_MANAGER_UNIFORM_H
#define MATPAR_MANAGER_UNIFORM_H

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace STR
{

namespace INVANA

{

  class MatParManagerUniform  : public MatParManager
  {

  protected:

    void FillParameters(Teuchos::RCP<Epetra_MultiVector> params);

    void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint, double val, int elepos, int parapos);

    void Consolidate(Teuchos::RCP<Epetra_MultiVector> dfint);

  public:

    //! Constructor
    MatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret);

    //!destructor
    ~MatParManagerUniform(){};
  };

  class MatParManagerPerElement  : public MatParManager
  {

  protected:

    void FillParameters(Teuchos::RCP<Epetra_MultiVector> params);

    void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint, double val, int elepos, int parapos);

    //! nothing to be done in case of elementwise parameters
    void Consolidate(Teuchos::RCP<Epetra_MultiVector> dfint){return;}

  public:

    //! Constructor
    MatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret);

    //!destructor
    ~MatParManagerPerElement(){};
  };

}
}

/*----------------------------------------------------------------------*/
#endif
