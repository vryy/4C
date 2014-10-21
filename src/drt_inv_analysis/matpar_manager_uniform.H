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
public:

  //! constructor
  MatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret);

  //! destructor
  virtual ~MatParManagerUniform(){}

  //! Setup specific parametrization layout
  virtual void Setup();

protected:

  virtual void FillParameters(Teuchos::RCP<Epetra_MultiVector> params);

  virtual void InitParameters(int parapos, double val);

  virtual void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint, double val, int elepos, int paraposglobal, int paraposlocal);

  virtual void Consolidate(Teuchos::RCP<Epetra_MultiVector> dfint);

};
}  // namespace INVANA
}  // namespace STR

/*----------------------------------------------------------------------*/
#endif
