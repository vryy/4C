/*!------------------------------------------------------------------------------------------------*
 \file ad_fld_poro.H

 \brief Fluid field adapter for poroelasticity

 \maintainer   Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251

*----------------------------------------------------------------------*/


#ifndef ADAPTER_FLUID_PORO_H
#define ADAPTER_FLUID_PORO_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "ad_fld_fluid_fpsi.H"
#include "../linalg/linalg_mapextractor.H"

#include "../drt_poroelast/poroelast_utils.H"

namespace DRT
{
  class Condition;
}

namespace ADAPTER {

  class FluidPoro : public FluidFPSI
  {
  public:

    //! Constructor
    FluidPoro(Teuchos::RCP<Fluid> fluid,
        Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output,
        bool isale,
        bool dirichletcond);

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
                           Teuchos::RCP<LINALG::SparseMatrix> ConstraintMatrix,
                           Teuchos::RCP<LINALG::SparseMatrix> StructVelConstraintMatrix,
                           Teuchos::RCP<Epetra_Vector> condVector,
                           Teuchos::RCP<std::set<int> > condIDs,
                           POROELAST::coupltype coupltype = POROELAST::fluidfluid);

    virtual void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

    virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

    virtual Teuchos::RCP<LINALG::MapExtractor> VelPresSplitter();

    /*!
      \brief Write extra output for specified step and time.
             Useful if you want to write output every iteration in partitioned schemes.
             If no step and time is provided, standard Output of fluid field is invoked.

      \param step (in) : Pseudo-step for which extra output is written
      \param time (in) : Pseudo-time for which extra output is written

      \note This is a pure DEBUG functionality. Originally used in immersed method development.

      \warning This method partly re-implements redundantly few lines of the common fluid Output() routine.
      \return void
    */
    virtual void Output(const int step=-1, const double time=-1);

  private:

    /// fluid field
    const Teuchos::RCP<ADAPTER::Fluid>& FluidField() { return fluid_; }

    std::vector<DRT::Condition*> nopencond_; ///< vector containing no penetration conditions

  };
}

#endif
