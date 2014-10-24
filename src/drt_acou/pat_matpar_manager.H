/*----------------------------------------------------------------------*/
/*!
\file pat_matpar_manager.H
\brief Class to handle calls to material parameters from an optimization routine

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            089 - 289-15271
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef PAT_MATPAR_MANAGER_H
#define PAT_MATPAR_MANAGER_H

#include "../drt_inv_analysis/matpar_manager_uniform.H"
#include "../drt_inv_analysis/matpar_manager_elementwise.H"
#include "Epetra_MultiVector.h"

namespace ACOU
{

/*----------------------------------------------------------------------*/
class PatMatParManagerUniform : public INVANA::MatParManagerUniform
{
public:

  //! Constructor
  PatMatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret);

  //! Destructor
  ~PatMatParManagerUniform(){};

  //! evaluate derivative of the residual force vector w.r.t to the material parameters using duals
  void Evaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint, bool consolidate=false);

};

/*----------------------------------------------------------------------*/
class PatMatParManagerPerElement : public INVANA::MatParManagerPerElement
{
public:

  //! Constructor
  PatMatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret, bool scaleele);

  //! Destructor
  ~PatMatParManagerPerElement(){};

  //! evaluate derivative of the residual force vector w.r.t to the material parameters using duals
  void Evaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint, bool consolidate=false);

private:
  bool scalegradele_;
};

} // namespace ACOU

#endif /* PAT_MATPAR_MANAGER_H */
