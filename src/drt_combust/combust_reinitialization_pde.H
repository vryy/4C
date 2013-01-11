/*!-----------------------------------------------------------------------------------------------*
\file combust_reinitialization_pde.H

\brief
<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef COMBUST_REINITIALIZATION_PDE_H
#define COMBUST_REINITIALIZATION_PDE_H


#include <Teuchos_RCP.hpp>


class Epetra_Vector;
class Epetra_Comm;

namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
}

namespace SCATRA
{
  class ScaTraTimIntImpl;
}

namespace Teuchos
{
  class ParameterList;
}

namespace COMBUST
{

//! Base class for PDE based reinitialization
/*!
  This class exists to provide the basic functions for PDE based reinitializations. This class holds the functions
  for the pseudo timeloop, the steady state checks, update for vectors and output.
  Within the constructor a new object of ScaTraBaseAlgorithm is created with a reinitswitch flag!

  see also "Level set transport and PDE-based reinitialization techniques" by Benedikt Schott
 */
class ReinitializationPDE
{
public:

  //! constructor
  ReinitializationPDE(const Epetra_Comm& comm);

  //! destructor
  virtual ~ReinitializationPDE();

  //! checks the steady state via numsteps or a gradient norm check
  void CheckSteadyState(double&    Gradient_Error_old,  //<! gradient based L1 error at old pseudo time step
                        bool&      STOP                 //<! boolian to abort
                        );

  //! calls the whole pde reinitialization after one global time step
  bool CallReinitialization(
      Teuchos::RCP<Epetra_Vector> phinp //<! pointer to phinp
      );

  //! routine to return scalar field phi at time step n+1
  Teuchos::RCP<Epetra_Vector> Phinp();

  //! output to file for reinitialized scatra field
  void OutputReinit(const int globalstep,   //<! time step of global problem
                    const double globaltime //<! time of global problem
                    );

private:

  /*========================================================================*/
  //! @name output functions
  /*========================================================================*/

  /// output reinitialization steps
  void OutputReinitializationSteps();

  /// write state vectors to Gmsh postprocessing files
  void OutputToGmshReinitializationSteps(
       const int reinit_step,           //<! reinitialization step
       const double reinit_pseudo_time  //<! reinitialization pseudo time
  );

  /// write state vectors (phinp) to Gmsh postprocessing files
  void OutputToGmshReinit(
       const int globalstep,    //<! time step of global problem
       const double globaltime, //<! time of global problem
       const int step,          //<! reinitialization step
       const double time        //<! reinitialization pseudo time
  );

  //! set the reinit vectors
  void SetPhiReinit(
      Teuchos::RCP<Epetra_Vector> phi //<! phi vector that has to be set
      );

  //! calls the inner timeloop for one reinitialization
  void TimeLoop_Reinit();

  //! print reinitialization info to screen
  void ReinitializeInfo();

  //! return the ScaTraTimIntImpl object of the ScaTraBaseAlgorithm
  SCATRA::ScaTraTimIntImpl& ScaTraReinitField();

  //! evaluate a gradient based error norm
  double EvaluateGradientNormError();


  //! ScaTraBaseAlgorithm object
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> reinit_;

  //! modified combustdyn parameter list with additional parameters used for reinitialization
  Teuchos::RCP<Teuchos::ParameterList> combustdynreinit_;

  //! is reinitialization accepted?
  bool reinitialization_accepted_;

  //! processor id
  int myrank_;

  //! maximal number of pseudo time steps
  int stepmax_;

}; // class ReinitializationPDE

} // namespace COMBUST

#endif //COMBUST_REINITIALIZATION_PDE_H
