/*----------------------------------------------------------------------*/
/*!
 * \file objective_funct_disp.H
 *
<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef OBJECTIVE_FUNCT_DISP_H
#define OBJECTIVE_FUNCT_DISP_H

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


  /**
  \brief An objective function to be used by the inverse analysis

  This is the standard least square displacement based error measure.
  Dof's to be used for computation must be given in a monitor file using
  this specific format:

  --MONITOR FILE START (this line is not part of the monitor file)

  steps [numtimesteps to be monitored] nnodes [number of nodes being monitored]\n
  [baci internal nodeid 1] [numdofs of this nodes] [\<local dofids\>]\n
  ...\n
  [baci internal nodeid n] [numdofs of this nodes] [\<local dofids\>]\n
  #arbitrary number of comment lines beginning with [#]\n
  [time(1)] [[value of dof1 node(1)] ... [value of dofn node(1)]] ... [[value of dof1 node(n)] ... [value of dofn node(1)]]\n
  ...\n
  [time(n)] [[value of dof1 node(1)] ... [value of dofn node(1)]] ... [[value of dof1 node(n)] ... [value of dofn node(1)]]\n

  -- MONITOR FILE END (this line is not part of the monitor file)\n

  The number of monitored time steps must not exceed the number of timesteps simulated.

  This class provides access to the objective functions value as well as access to its derivative w.r.t
  the displacments.

  */
class ObjectiveFunctDisp : public ObjectiveFunct
{
private:

  void ReadMonitor(std::string monitorfilename);

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
  ObjectiveFunctDisp(Teuchos::RCP<DRT::Discretization> discret, int steps, Teuchos::RCP<std::vector<double> > timesteps);

  //! standard destructor
  ~ObjectiveFunctDisp(){};

  //! evaluate value of the objective function
  virtual void Evaluate(Teuchos::RCP<Epetra_MultiVector> disp, double& val);

  //! evaluate gradient of objective function w.r.t to displacements
  virtual void EvaluateGradient(Teuchos::RCP<Epetra_MultiVector> disp, Teuchos::RCP<Epetra_MultiVector> gradient);

}; // class ObjectiveFunctDisp
}  // namespace INVANA
}  // namespace STR

#endif // OBJECTIVE_FUNCT_DISP_H
