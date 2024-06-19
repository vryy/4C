/*----------------------------------------------------------------------*/
/*! \file

 \brief Fluid field adapter for poroelasticity



\level 2

*----------------------------------------------------------------------*/


#ifndef FOUR_C_ADAPTER_FLD_PORO_HPP
#define FOUR_C_ADAPTER_FLD_PORO_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_fluid_fpsi.hpp"
#include "4C_fem_condition.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_poroelast_utils.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class FluidPoro : public FluidFPSI
  {
   public:
    //! Constructor
    FluidPoro(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    //! Evaluate no penetration constraint
    /*!
     \param Cond_RHS                  (o) condition part of rhs
     \param ConstraintMatrix          (o) static part of Fluid matrix associated with constraints
     \param struct_vel_constraint_matrix (o) transient part of Fluid matrix associated with
     constraints \param condIDs                   (o) vector containing constraint dofs \param
     coupltype                 (i) coupling type, determines which matrix is to be evaluated (0==
     fluid-fluid, 1== fluid -structure)
     */
    void evaluate_no_penetration_cond(Teuchos::RCP<Epetra_Vector> Cond_RHS,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> ConstraintMatrix,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> struct_vel_constraint_matrix,
        Teuchos::RCP<Epetra_Vector> condVector, Teuchos::RCP<std::set<int>> condIDs,
        PoroElast::Coupltype coupltype = PoroElast::fluidfluid);

    //! calls the VelPresSplitter on the time integrator
    virtual Teuchos::RCP<Core::LinAlg::MapExtractor> VelPresSplitter();

    /*!
      \brief Write extra output for specified step and time.
             Useful if you want to write output every iteration in partitioned schemes.
             If no step and time is provided, standard Output of fluid field is invoked.

      \param step (in) : Pseudo-step for which extra output is written
      \param time (in) : Pseudo-time for which extra output is written

      \note This is a pure DEBUG functionality. Originally used in immersed method development.

      \warning This method partly re-implements redundantly few lines of the common fluid output()
      routine. \return void
    */
    virtual void output(const int step = -1, const double time = -1);

   private:
    /// fluid field
    const Teuchos::RCP<Adapter::Fluid>& fluid_field() { return fluid_; }

    std::vector<Core::Conditions::Condition*>
        nopencond_;  ///< vector containing no penetration conditions
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
