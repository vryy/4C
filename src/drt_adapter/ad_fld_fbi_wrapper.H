/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field wrapper for fluid beam interactions

\level 2

\maintainer Nora Hagmeyer
*/
/*----------------------------------------------------------------------*/
#ifndef AD_FLD_FBI_WRAPPER_H
#define AD_FLD_FBI_WRAPPER_H

#include <Teuchos_RCP.hpp>

#include "ad_fld_fluid_fsi.H"

namespace LINALG
{
  class Solver;
}  // namespace LINALG

namespace IO
{
  class DiscretizationWriter;
}


namespace ADAPTER
{
  /*! \brief Fluid field adapter for fluid beam interaction
   *
   *
   *  Can only be used in conjunction with #FLD::FluidImplicitTimeInt
   */
  class FluidFBI : public FluidFSI
  {
   public:
    /// Constructor
    FluidFBI(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    /** \brief Pass in additional contributions from coupling terms for the system matrix
     *
     * To enforce weak dirichlet conditions as they arise from meshtying for example, such
     * contributions can be handed to the fluid, which will store the pointer on the coupling
     * contributions to assemble them into the system matrix in each Newton iteration.
     *
     * \param[in] matrix (size fluid_dof x fluid_dof) matrix containing weak dirichlet entries that
     * need to be assembled into the overall fluid system matrix
     */
    virtual void SetCouplingContributions(Teuchos::RCP<const LINALG::SparseMatrix> matrix);

    /**
     * \brief Resets the external forces acting on the fluid to zero
     */
    virtual void ResetExternalForces();
  };
}  // namespace ADAPTER

#endif
