/*!------------------------------------------------------------------------------------------------*
 \file adapter_structure_timint_poroelast.cpp

 \brief Structure field adapter for poroelasticity

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
#include "adapter_structure_timint_poroelast.H"
//#include "../drt_io/io.H"
//#include "../drt_lib/drt_condition_utils.H"

//constructor
ADAPTER::StructureTimIntImplPoro::StructureTimIntImplPoro(Teuchos::RCP<
    STR::TimIntImpl> stii, Teuchos::RCP<Teuchos::ParameterList> ioparams, ///< I/O flags
    Teuchos::RCP<Teuchos::ParameterList> sdynparams, ///< input parameters
    Teuchos::RCP<Teuchos::ParameterList> xparams, ///< extra flags
    Teuchos::RCP<DRT::Discretization> discret, ///< current discretisation
    Teuchos::RCP<LINALG::Solver> solver, ///< the solver
    Teuchos::RCP<LINALG::Solver> contactsolver, ///< the contact solver
    Teuchos::RCP<IO::DiscretizationWriter> output ///< the output
) :
  StructureTimIntImpl(stii, ioparams, sdynparams, xparams, discret, solver,
      contactsolver, output), structure_(stii)
{
  //reset stiffness matrx and rhs
  structure_->PoroInitForceStiffResidual(sdynparams);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimIntImplPoro::Evaluate(Teuchos::RCP<
    const Epetra_Vector> disiterinc)
{
  structure_->UpdateIterIncrementally(disiterinc);

  // builds tangent, residual and applies DBC
  structure_->PoroEvaluateForceStiffResidual();
  structure_->PrepareSystemForNewtonSolve();
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
