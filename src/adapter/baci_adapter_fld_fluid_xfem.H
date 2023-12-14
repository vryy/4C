/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for xfem-fluids with moving boundaries

\level 1


*/
/*----------------------------------------------------------------------*/


#ifndef BACI_ADAPTER_FLD_FLUID_XFEM_H
#define BACI_ADAPTER_FLD_FLUID_XFEM_H

#include "baci_adapter_fld_base_algorithm.H"
#include "baci_adapter_fld_moving_boundary.H"

BACI_NAMESPACE_OPEN

namespace ADAPTER
{
  // forward declarations
  class Coupling;

  /// fluid with moving interfaces implemented by the XFEM
  class FluidXFEM : public FluidMovingBoundary
  {
   public:
    /// constructor
    explicit FluidXFEM(const Teuchos::ParameterList& prbdyn, std::string condname);

    /*========================================================================*/
    //! @name Misc
    /*========================================================================*/

    /// fluid field
    const Teuchos::RCP<ADAPTER::Fluid>& FluidField() override { return fluid_; }

    /// return the boundary discretization that matches the structure discretization
    Teuchos::RCP<DRT::Discretization> Discretization() override;

    /// return the boundary discretization that matches the structure discretization
    Teuchos::RCP<DRT::Discretization> BoundaryDiscretization();

    /// communication object at the interface
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& Interface() const override
    {
      return fluid_->Interface();
    }

    /// communication object at the struct interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor> const& StructInterface();

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

    /// get the maximum number of iterations from the fluid field
    int Itemax() const override { return fluid_->Itemax(); }

    /// set the maximum number of iterations for the fluid field
    void SetItemax(int itemax) override { fluid_->SetItemax(itemax); }

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
    Teuchos::RCP<Fluid> fluid_;
  };

}  // namespace ADAPTER

BACI_NAMESPACE_CLOSE

#endif
