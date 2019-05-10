/*----------------------------------------------------------------------*/
/*!
 * \file objective_funct_disp.H
 * \brief Displacement based objective function
 *
 *
\level 3
\maintainer Sebastian Brandstaeter
            brandstaeter@lnm.mw.tum.de
            089 - 289-15276
*/
/*----------------------------------------------------------------------*/

#ifndef OBJECTIVE_FUNCT_DISP_H
#define OBJECTIVE_FUNCT_DISP_H

#include "objective_funct.H"
#include "../linalg/linalg_mapextractor.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace INVANA
{
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
  [time(1)] [[value of dof1 node(1)] ... [value of dofn node(1)]] ... [[value of dof1 node(n)] ...
  [value of dofn node(1)]]\n
  ...\n
  [time(n)] [[value of dof1 node(1)] ... [value of dofn node(1)]] ... [[value of dof1 node(n)] ...
  [value of dofn node(1)]]\n

  -- MONITOR FILE END (this line is not part of the monitor file)\n

  The number of monitored time steps must not exceed the number of timesteps simulated.

  This class provides access to the objective functions value as well as access to its derivative
  w.r.t the displacments.

  */
  class ObjectiveFunctDisp : public ObjectiveFunct
  {
   private:
    void ReadMonitor(std::string monitorfilename);

    //! discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! dofrowmap of the discretization
    const Epetra_Map* dofrowmap_;

    //! dofrowmap of measured dofs only
    Teuchos::RCP<Epetra_Map> mdofrowmap_;

    Teuchos::RCP<LINALG::MapExtractor> mextractor_;

    //! measured displacements/state
    Teuchos::RCP<Epetra_MultiVector> mstate_;

    //! simulated timesteps
    std::vector<double> timesteps_;

    //! scaling of the objective function?
    bool scaling_;

    //! scale factor
    double scalefac_;

    //! measurement variance estimation
    double var_estim_;


   public:
    //! constructor
    ObjectiveFunctDisp(Teuchos::RCP<DRT::Discretization> discret);

    //! standard destructor
    virtual ~ObjectiveFunctDisp(){};

    //! evaluate value of the objective function
    virtual void Evaluate(Teuchos::RCP<Epetra_Vector> state, double time, double& val);

    //! evaluate gradient of objective function w.r.t to displacements
    virtual void EvaluateGradient(
        Teuchos::RCP<Epetra_Vector> state, double time, Teuchos::RCP<Epetra_Vector> gradient);

    //! measured time steps
    virtual std::vector<double> MeasuredTime() { return timesteps_; }

    //! query for step of measurements of this time step
    int FindStep(double time);

    virtual double GetScaleFac() { return scalefac_; }

  };  // class ObjectiveFunctDisp
}  // namespace INVANA

#endif  // OBJECTIVE_FUNCT_DISP_H
