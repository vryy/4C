/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for immersed-fluids

\level 3

*/
/*----------------------------------------------------------------------*/


#ifndef BACI_ADAPTER_FLD_FLUID_IMMERSED_H
#define BACI_ADAPTER_FLD_FLUID_IMMERSED_H

#include "baci_adapter_fld_base_algorithm.H"
#include "baci_adapter_fld_moving_boundary.H"

namespace ADAPTER
{
  // forward declarations
  class Coupling;

  /// fluid with moving interfaces
  class FluidImmersed : public FluidMovingBoundary
  {
   public:
    /// constructor
    explicit FluidImmersed(const Teuchos::ParameterList& prbdyn, std::string condname);

    /*========================================================================*/
    //! @name Misc
    /*========================================================================*/

    /// fluid field
    const Teuchos::RCP<ADAPTER::Fluid>& FluidField() override { return fluid_->FluidField(); }

    /// return the boundary discretization that matches the structure discretization
    Teuchos::RCP<DRT::Discretization> Discretization() override;

    /// communication object at the interface
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& Interface() const override;

    //@}

    /*========================================================================*/
    //! @name Time step helpers
    /*========================================================================*/

    /// start new time step
    void PrepareTimeStep() override;

    /// update at time step end
    void Update() override;

    /// output results
    void Output() override;

    /// read restart information for given time step
    double ReadRestart(int step) override;

    /*========================================================================*/
    //! @name Solver calls
    /*========================================================================*/

    /// nonlinear solve
    void NonlinearSolve(
        Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel) override;

    /// relaxation solve
    Teuchos::RCP<Epetra_Vector> RelaxationSolve(
        Teuchos::RCP<Epetra_Vector> idisp, double dt) override;
    //@}

    /*========================================================================*/
    //! @name Extract interface forces
    /*========================================================================*/

    /// After the fluid solve we need the forces at the FSI interface.
    Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces() override;
    //@}

    /*========================================================================*/
    //! @name extract helpers
    /*========================================================================*/

    /// extract the interface velocity at time t^(n+1)
    Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp() override;

    /// extract the interface velocity at time t^n
    Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln() override;
    //@}

    /*========================================================================*/
    //! @name Number of Newton iterations
    /*========================================================================*/

    //! For simplified FD MFNK solve we want to temporally limit the
    /// number of Newton steps inside the fluid solver
    /// Not used for IMMERSED yet !

    /// get the maximum number of iterations from the fluid field
    int Itemax() const override { return fluidadapter_->Itemax(); }

    /// set the maximum number of iterations for the fluid field
    void SetItemax(int itemax) override { FluidField()->SetItemax(itemax); }

    /// add dirichlet conditions to dirichlet condmap before next fluid solve
    virtual void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

    /// remove dirichlet conditions from dirichlet condmap after last fluid solve
    virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

    //@}

    /*========================================================================*/
    //! @name others
    /*========================================================================*/

    /// integrate the interface shape functions
    Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape() override;

    /// create the testing of fields
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest() override;


   private:
    /// fluid base algorithm object
    Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid_;
    Teuchos::RCP<ADAPTER::Fluid> fluidadapter_;
  };

}  // namespace ADAPTER

#endif
