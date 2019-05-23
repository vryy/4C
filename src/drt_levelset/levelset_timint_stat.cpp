/*!----------------------------------------------------------------------
\file levelset_timint_stat.cpp

\brief stationary time integration scheme for level-set problems (for coupled problems only)
       just a dummy

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236

*----------------------------------------------------------------------*/

#include "levelset_timint_stat.H"

//#include "../drt_scatra_ele/scatra_ele_action.H"
//#include <Teuchos_StandardParameterEntryValidators.hpp>
//#include <Teuchos_TimeMonitor.hpp>
//#include "../drt_io/io.H"
//#include "../drt_io/io_pstream.H"
//#include "../linalg/linalg_solver.H"
//
//#include "../drt_particle/scatra_particle_coupling.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                rasthofer 09/13 |
 *----------------------------------------------------------------------*/
SCATRA::LevelSetTimIntStationary::LevelSetTimIntStationary(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      LevelSetAlgorithm(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntStationary(actdis, solver, sctratimintparams, extraparams, output)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized
  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                              rasthofer 09/13 |
*-----------------------------------------------------------------------*/
SCATRA::LevelSetTimIntStationary::~LevelSetTimIntStationary() { return; }


/*----------------------------------------------------------------------*
 |  initialize time integration                             rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetTimIntStationary::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntStationary::Init();
  LevelSetAlgorithm::Init();

  if (myrank_ == 0)
  {
    std::cout << "\n*------------------------------------------------------------------------*"
              << std::endl;
    std::cout << "|            stationary level set for coupled problems only              |"
              << std::endl;
    std::cout << "*------------------------------------------------------------------------*\n"
              << std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*
 |  setup time integration                                  rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetTimIntStationary::Setup()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntStationary::Setup();
  LevelSetAlgorithm::Setup();

  if (myrank_ == 0)
  {
    std::cout << "\n*------------------------------------------------------------------------*"
              << std::endl;
    std::cout << "|            stationary level set for coupled problems only              |"
              << std::endl;
    std::cout << "*------------------------------------------------------------------------*\n"
              << std::endl;
  }

  return;
}
