/*----------------------------------------------------------------------*/
/*!
 * \file objective_funct_surfcurr.H
 *
<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef OBJECTIVE_FUNCT_SURFCURR_H
#define OBJECTIVE_FUNCT_SURFCURR_H

#include "objective_funct.H"

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "../drt_structure/strtimint_mstep.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace STR
{

namespace INVANA
{

class MatParManager;


  /*!
  \brief An objective function based on a surface current representation

  documentation will follow

  */
class ObjectiveFunctSurfCurr : public ObjectiveFunct
{
private:

  //! discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  //! dofrowmap of the discretization
  const Epetra_Map* dofrowmap_;

  //! dof mask on measurements
  Teuchos::RCP<Epetra_Vector> mask_;

  //! measured displacements
  Teuchos::RCP<Epetra_MultiVector> mdisp_;

  //! simulated timesteps
  Teuchos::RCP<std::vector<double> > timesteps_;

  //! number of vectors in MultiVectors mdisp_ and disp
  int msteps_;

public:

  //! constructor
  ObjectiveFunctSurfCurr(Teuchos::RCP<DRT::Discretization> discret, int steps, Teuchos::RCP<std::vector<double> > timesteps);

  //! standard destructor
  ~ObjectiveFunctSurfCurr(){};

  //! evaluate value of the objective function
  virtual void Evaluate(Teuchos::RCP<Epetra_MultiVector> disp, double& val);

  //! evaluate gradient of objective function w.r.t to displacements
  virtual void EvaluateGradient(Teuchos::RCP<Epetra_MultiVector> disp, Teuchos::RCP<Epetra_MultiVector> gradient);

}; // class ObjectiveFunctSurfCurr
}  // namespace INVANA
}  // namespace STR

#endif // OBJECTIVE_FUNCT_SURFCURR_H
