/*----------------------------------------------------------------------*/
/*! \file

 \brief Fluid field adapter for poroelasticity



\level 2

*----------------------------------------------------------------------*/


#ifndef BACI_ADAPTER_FLD_PORO_H
#define BACI_ADAPTER_FLD_PORO_H

#include "baci_adapter_fld_fluid_fpsi.H"
#include "baci_linalg_mapextractor.H"
#include "baci_poroelast_utils.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Condition;
}

namespace ADAPTER
{
  class FluidPoro : public FluidFPSI
  {
   public:
    //! Constructor
    FluidPoro(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    //! Evaluate no penetration constraint
    /*!
     \param Cond_RHS                  (o) condition part of rhs
     \param ConstraintMatrix          (o) static part of Fluid matrix associated with constraints
     \param StructVelConstraintMatrix (o) transient part of Fluid matrix associated with constraints
     \param condIDs                   (o) vector containing constraint dofs
     \param coupltype                 (i) coupling type, determines which matrix is to be evaluated
                                         (0== fluid-fluid, 1== fluid -structure)
     */
    void EvaluateNoPenetrationCond(Teuchos::RCP<Epetra_Vector> Cond_RHS,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> ConstraintMatrix,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> StructVelConstraintMatrix,
        Teuchos::RCP<Epetra_Vector> condVector, Teuchos::RCP<std::set<int>> condIDs,
        POROELAST::coupltype coupltype = POROELAST::fluidfluid);

    //! calls the VelPresSplitter on the time integrator
    virtual Teuchos::RCP<CORE::LINALG::MapExtractor> VelPresSplitter();

    /*!
      \brief Write extra output for specified step and time.
             Useful if you want to write output every iteration in partitioned schemes.
             If no step and time is provided, standard Output of fluid field is invoked.

      \param step (in) : Pseudo-step for which extra output is written
      \param time (in) : Pseudo-time for which extra output is written

      \note This is a pure DEBUG functionality. Originally used in immersed method development.

      \warning This method partly re-implements redundantly few lines of the common fluid Output()
      routine. \return void
    */
    virtual void Output(const int step = -1, const double time = -1);

   private:
    /// fluid field
    const Teuchos::RCP<ADAPTER::Fluid>& FluidField() { return fluid_; }

    std::vector<DRT::Condition*> nopencond_;  ///< vector containing no penetration conditions
  };
}  // namespace ADAPTER

BACI_NAMESPACE_CLOSE

#endif
