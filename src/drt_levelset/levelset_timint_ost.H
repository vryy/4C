/*!----------------------------------------------------------------------
\file levelset_timint_ost.H
\brief

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifndef LEVELSET_TIMINT_OST_H
#define LEVELSET_TIMINT_OST_H

#include "../drt_scatra/scatra_timint_ost.H"
#include "levelset_algorithm.H"


namespace SCATRA
{

class LevelSetTimIntOneStepTheta : public LevelSetAlgorithm, public TimIntOneStepTheta
{

public:

  /// Standard Constructor
  LevelSetTimIntOneStepTheta(Teuchos::RCP<DRT::Discretization>       dis,
      Teuchos::RCP<LINALG::Solver>            solver,
      Teuchos::RCP<Teuchos::ParameterList>    params,
      Teuchos::RCP<Teuchos::ParameterList>    sctratimintparams,
      Teuchos::RCP<Teuchos::ParameterList>    extraparams,
      Teuchos::RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~LevelSetTimIntOneStepTheta();

  /// initialize time-integration scheme
  virtual void Init();

  /// read restart data
  virtual void ReadRestart(int step);

  /// redistribute the scatra discretization and vectors according to nodegraph
  void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph);

  /// Interpolate phi(alpha_F) values from OST calculations.
  Teuchos::RCP<Epetra_Vector> PhiafOst(const double alpha_F);

  /// Interpolate phi(alpha_F) values from OST calculations.
  Teuchos::RCP<Epetra_Vector> PhiamOst(const double alpha_M);

  /// Interpolate phi(alpha_F) values from OST calculations.
  Teuchos::RCP<Epetra_Vector> PhidtamOst(const double alpha_M);

  void SetAlphaF(const double alpha_F){alphaF_=alpha_F;};

  /// return phi at evaluation time
  virtual Teuchos::RCP<Epetra_Vector> EvaluationPhi(){
    if(alphaF_>=0.0){
      return PhiafOst(alphaF_);
    }
    else{
      return phinp_;
    }
  };

protected:

  /// Print information about current time step to screen (reimplementation for OST)
  void PrintTimeStepInfo();

  /// Initialization procedure before the first time step
  void PrepareFirstTimeStep();

  /// additional predictor not intended for level-set methods
  void ExplicitPredictor(){return;};

  /// Set the part of the righthandside belonging to the last timestep.
  void SetOldPartOfRighthandside();

  /// update state vectors
  /// current solution becomes old solution of next time step
  void UpdateState();

  /// update the solution after Solve()
  /// extended version for coupled level-set problems including reinitialization and particle correction if required
  void Update(const int num=0);

  /// update phi within the reinitialization loop
  void UpdateReinit();

private:

  //Member introduced to make the OST level set algorithm compatible with fluid gen alpha algorithm.
  double alphaF_;

}; // class LevelSetTimIntOneStepTheta

} // namespace SCATRA

#endif  // #ifndef LEVELSET_TIMINT_OST_H
