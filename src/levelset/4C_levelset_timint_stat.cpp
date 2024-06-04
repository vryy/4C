/*----------------------------------------------------------------------*/
/*! \file

\brief stationary time integration scheme for level-set problems (for coupled problems only)
       just a dummy

\level 2


*----------------------------------------------------------------------*/

#include "4C_levelset_timint_stat.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Constructor (public)                                rasthofer 09/13 |
 *----------------------------------------------------------------------*/
SCATRA::LevelSetTimIntStationary::LevelSetTimIntStationary(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams,
    Teuchos::RCP<CORE::IO::DiscretizationWriter> output)
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

FOUR_C_NAMESPACE_CLOSE
