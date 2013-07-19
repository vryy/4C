/*----------------------------------------------------------------------*/
/*!
\file ad_fld_base_algorithm.H

\brief Fluid Base Algorithm

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/



#ifndef ADAPTER_FLUID_BASE_ALGORITHM_H
#define ADAPTER_FLUID_BASE_ALGORITHM_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "ad_fld_fluid.H"


// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {

  /// fluid field solver
  class FluidBaseAlgorithm
  {
  public:

    /// constructor
    explicit FluidBaseAlgorithm(const Teuchos::ParameterList& prbdyn, bool isale);

    /// second constructor (special version for turbulent flows with separate inflow
    /// section for generation of turbulent inflow profiles)
    explicit FluidBaseAlgorithm(const Teuchos::ParameterList& prbdyn,
                                const Teuchos::RCP<DRT::Discretization> discret,
                                const Teuchos::RCP<std::map<int,std::vector<int> > > row_pbcmapmastertoslave,
                                const Teuchos::RCP<std::map<int,std::vector<int> > > col_pbcmapmastertoslave);

    /// virtual destructor to support polymorph destruction
    virtual ~FluidBaseAlgorithm();

    /// fluid field solver
    Fluid& FluidField() { return *fluid_; }
    /// return instance of fluid time integration class
    const Fluid& FluidField() const { return *fluid_; }
    /// return instance of fluid time integration class
    Teuchos::RCP<Fluid> FluidFieldrcp() { return fluid_; }

  private:

    //! setup fluid algorithm (overriding some fluid parameters with
    //! values specified in given problem-dependent ParameterList)
    /**
     * \note In this function the linear solver object is generated. For pure fluid problems or
     * fluid meshtying (no block matrix) the FLUID SOLVER block from the BACI dat file is used.
     * For fluid meshtying (block matrix) the MESHTYING SOLVER block is used as main solver object with a
     * block preconditioner (BGS or SIMPLE type). The block preconditioners use the information form the
     * FLUID SOLVER and the FLUID PRESSURE SOLVER block for the velocity and pressure dofs.
     */
    void SetupFluid(const Teuchos::ParameterList& prbdyn, bool& isale);

    /// setup second fluid algorithm (overriding some fluid parameters with
    /// values specified in given problem-dependent Turbulent Inflow ParameterList)
    /// separate discretization for inflow generation
    void SetupInflowFluid(const Teuchos::ParameterList& prbdyn,
                          const Teuchos::RCP<DRT::Discretization> discret,
                          const Teuchos::RCP<std::map<int,std::vector<int> > > row_pbcmapmastertoslave,
                          const Teuchos::RCP<std::map<int,std::vector<int> > > col_pbcmapmastertoslave);

    //! set parameters in list required for all schemes
    void SetGeneralParameters(const Teuchos::RCP<Teuchos::ParameterList> fluidtimeparams,
                              const Teuchos::ParameterList& prbdyn,
                              const Teuchos::ParameterList& fdyn);

    //! create a second solver
    void CreateSecondSolver(const Teuchos::RCP<LINALG::Solver> solver,
                            const Teuchos::ParameterList& fdyn);

    /// fluid field solver
    Teuchos::RCP<Fluid> fluid_;
  };

}

#endif
