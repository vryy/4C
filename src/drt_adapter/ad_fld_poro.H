/*!------------------------------------------------------------------------------------------------*
 \file ad_fld_poro.H

 \brief Fluid field adapter for poroelasticity

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
*----------------------------------------------------------------------*/


#ifndef ADAPTER_FLUID_PORO_H
#define ADAPTER_FLUID_PORO_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include "ad_fld_fluid_fsi.H"
#include "../linalg/linalg_mapextractor.H"

namespace DRT
{
  class Condition;
}

namespace ADAPTER {

  class FluidPoro : public FluidFSI
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
                           Teuchos::RCP<std::set<int> > condIDs,
                           int coupltype=0);

    //! Evaluate no penetration constraint
    /*!
     \param FluidFluidMatrix          (o) fluid-fluid coupling matrix associated with boundary integral
     \param FluidStructMatrix         (o) fluid-structure coupling matrix associated with boundary integral
    */
    void EvaluateContiBoundaryCond(
                           Teuchos::RCP<LINALG::SparseMatrix> FluidFluidMatrix,
                           Teuchos::RCP<LINALG::SparseMatrix> FluidStructMatrix);

    virtual void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

  private:

    std::vector<DRT::Condition*> nopencond_; ///< vector containing no penetration conditions

  };
}

#endif
