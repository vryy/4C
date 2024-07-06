/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for xfem-fluids with moving boundaries

\level 1


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_ADAPTER_FLD_FLUID_XFEM_HPP
#define FOUR_C_ADAPTER_FLD_FLUID_XFEM_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_base_algorithm.hpp"
#include "4C_adapter_fld_moving_boundary.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
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
    const Teuchos::RCP<Adapter::Fluid>& fluid_field() override { return fluid_; }

    /// return the boundary discretization that matches the structure discretization
    Teuchos::RCP<Core::FE::Discretization> discretization() override;

    /// return the boundary discretization that matches the structure discretization
    Teuchos::RCP<Core::FE::Discretization> boundary_discretization();

    /// communication object at the interface
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& interface() const override
    {
      return fluid_->interface();
    }

    /// communication object at the struct interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor> const& struct_interface();

    //@}

    /*========================================================================*/
    //! @name Time step helpers
    /*========================================================================*/

    /// start new time step
    void prepare_time_step() override;

    /// update at time step end
    void update() override;

    /// output results
    void output() override;

    /// read restart information for given time step
    double read_restart(int step) override;

    /*========================================================================*/
    //! @name Solver calls
    /*========================================================================*/

    /// nonlinear solve
    void nonlinear_solve(
        Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel) override;

    /// relaxation solve
    Teuchos::RCP<Epetra_Vector> relaxation_solve(
        Teuchos::RCP<Epetra_Vector> idisp, double dt) override;
    //@}

    /*========================================================================*/
    //! @name Extract interface forces
    /*========================================================================*/

    /// After the fluid solve we need the forces at the FSI interface.
    Teuchos::RCP<Epetra_Vector> extract_interface_forces() override;
    //@}

    /*========================================================================*/
    //! @name extract helpers
    /*========================================================================*/

    /// extract the interface velocity at time t^(n+1)
    Teuchos::RCP<Epetra_Vector> extract_interface_velnp() override;

    /// extract the interface velocity at time t^n
    Teuchos::RCP<Epetra_Vector> extract_interface_veln() override;
    //@}

    /*========================================================================*/
    //! @name Number of Newton iterations
    /*========================================================================*/

    //! For simplified FD MFNK solve we want to temporally limit the
    /// number of Newton steps inside the fluid solver

    /// get the maximum number of iterations from the fluid field
    int itemax() const override { return fluid_->itemax(); }

    /// set the maximum number of iterations for the fluid field
    void set_itemax(int itemax) override { fluid_->set_itemax(itemax); }

    //@}

    /*========================================================================*/
    //! @name others
    /*========================================================================*/

    /// integrate the interface shape functions
    Teuchos::RCP<Epetra_Vector> integrate_interface_shape() override;

    /// create the testing of fields
    Teuchos::RCP<Core::UTILS::ResultTest> create_field_test() override;



   private:
    /// fluid base algorithm object
    Teuchos::RCP<Fluid> fluid_;
  };

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
