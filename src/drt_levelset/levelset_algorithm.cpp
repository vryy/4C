/*!-----------------------------------------------------------------------------------------------*
\file levelset_algorithm.cpp

\brief base level-set algorithm

    detailed description in header file levelset_algorithm.H

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 *------------------------------------------------------------------------------------------------*/

#include "levelset_algorithm.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_particle/scatra_particle_coupling.H"


SCATRA::LevelSetAlgorithm::LevelSetAlgorithm(
  const Epetra_Comm& comm,              ///< communicator
  const Teuchos::ParameterList& prbdyn, ///< problem-specific parameters
  bool isale,  ///< do we need an ALE formulation of the fields?
  const std::string disname, ///< scatra discretization name
  const Teuchos::ParameterList& solverparams)
{

  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(prbdyn,isale,disname,solverparams));
  scatra_ =  scatra->ScaTraFieldrcp();

  particle_ = Teuchos::rcp(new PARTICLE::ScatraParticleCoupling(scatra_, comm, prbdyn));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SCATRA::LevelSetAlgorithm::~LevelSetAlgorithm()
{
}
