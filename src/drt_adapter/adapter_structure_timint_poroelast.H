/*!------------------------------------------------------------------------------------------------*
 \file adapter_structure_timint_poroelast.H

 \brief Structure field adapter for poroelasticity

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

#ifndef ADAPTER_STRUCTURE_TIMINT_POROELAST_H_
#define ADAPTER_STRUCTURE_TIMINT_POROELAST_H_

/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_Vector.h>

#include "../drt_adapter/adapter_algorithmbase.H"
#include "adapter_structure_timint.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*/
namespace ADAPTER
{

  class StructureTimIntImplPoro: public StructureTimIntImpl
  {
    public:

      /// Constructor
      StructureTimIntImplPoro(Teuchos::RCP<STR::TimIntImpl> stii, ///< time integrator
          Teuchos::RCP<Teuchos::ParameterList> ioparams, ///< I/O flags
          Teuchos::RCP<Teuchos::ParameterList> sdynparams, ///< input parameters
          Teuchos::RCP<Teuchos::ParameterList> xparams, ///< extra flags
          Teuchos::RCP<DRT::Discretization> discret, ///< current discretisation
          Teuchos::RCP<LINALG::Solver> solver, ///< the solver
          Teuchos::RCP<LINALG::Solver> contactsolver, ///< the contact solver
          Teuchos::RCP<IO::DiscretizationWriter> output ///< the output
          );

      virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc);

    private:

      /// the actual structure algorithm
      Teuchos::RCP<STR::TimIntImpl> structure_;

  };//class StructureTimIntImplPoro

} // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif

#endif /* ADAPTER_STRUCTURE_TIMINT_POROELAST_H_ */
