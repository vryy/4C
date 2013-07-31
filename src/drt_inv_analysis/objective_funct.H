/*----------------------------------------------------------------------*/
/*!
 * \file objective_funct.H
 *
<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef OBJECTIVE_FUNCT_H
#define OBJECTIVE_FUNCT_H

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

class ObjectiveFunct {

private:

  void ReadMonitor(std::string filename);

  // discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  // dofrowmap of the discretization
  const Epetra_Map* dofrowmap_;

  // dof mask on measurements
  Teuchos::RCP<Epetra_MultiVector> mask_;

  // measured displacements
  Teuchos::RCP<Epetra_MultiVector> mdisp_;

  // simulated timesteps
  Teuchos::RCP<std::vector<double> > timesteps_;

  // number of vectors in MultiVectors mdisp_ and disp
  int msteps_;


public:

  /// constructor
  ObjectiveFunct(Teuchos::RCP<DRT::Discretization> discret, int steps, Teuchos::RCP<std::vector<double> > timesteps);

  /// standard destructor
  ~ObjectiveFunct(){};

  //! evaluate difference measure
  double Evaluate(Teuchos::RCP<Epetra_MultiVector> disp, Teuchos::RCP<MatParManager> matman);

  //! evaluate gradient of objective function
  void EvaluateGradient(Teuchos::RCP<Epetra_MultiVector> disp, Teuchos::RCP<Epetra_MultiVector> gradient);

}; // class ObjectiveFunct
}  // namespace INVANA
}  // namespace STR

#endif /*INV_ANALYSIS_H_*/
