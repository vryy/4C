/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for immersed-fluids (beam)

\level 3

*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_ADAPTER_FLD_FBI_MOVINGBOUNDARY_HPP
#define FOUR_C_ADAPTER_FLD_FBI_MOVINGBOUNDARY_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_moving_boundary.hpp"
#include "4C_fluid_meshtying.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations

namespace Core::LinAlg
{
  class SparseOperator;
}

namespace Adapter
{
  class Coupling;
  class Fluid;

  /// fluid with moving interfaces
  class FBIFluidMB : public FluidMovingBoundary
  {
   public:
    /// constructor
    explicit FBIFluidMB(const Teuchos::ParameterList& prbdyn, std::string condname);

    /*========================================================================*/
    //! @name Misc
    /*========================================================================*/

    /// fluid field
    const Teuchos::RCP<Adapter::Fluid>& fluid_field() override { return fluidadapter_; }

    /// return the boundary discretization that matches the structure discretization
    Teuchos::RCP<Core::FE::Discretization> discretization() override;

    /// communication object at the interface
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& interface() const override;

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

    /// get the maximum number of iterations from the fluid field
    int itemax() const override;

    /// set the maximum number of iterations for the fluid field
    void set_itemax(int itemax) override;

    //@}

    /*========================================================================*/
    //! @name others
    /*========================================================================*/

    /// integrate the interface shape functions
    Teuchos::RCP<Epetra_Vector> integrate_interface_shape() override;

    /// create the testing of fields
    Teuchos::RCP<Core::UTILS::ResultTest> create_field_test() override;

    /// Get velocity at timestep n+1
    virtual Teuchos::RCP<const Epetra_Vector> velnp();

    virtual Teuchos::RCP<const FLD::Meshtying> get_meshtying();

    /** \brief Pass in additional contributions from coupling terms for the system matrix
     *
     * To enforce weak dirichlet conditions as they arise from meshtying for example, such
     * contributions can be handed to the fluid, which will store the pointer on the coupling
     * contributions to assemble them into the system matrix in each Newton iteration.
     *
     * \param[in] matrix (size fluid_dof x fluid_dof) matrix containing weak dirichlet entries that
     * need to be assembled into the overall fluid system matrix
     */
    virtual void set_coupling_contributions(
        Teuchos::RCP<const Core::LinAlg::SparseOperator> matrix);

    /**
     * \brief Pass additional contributions to the fluid residual to the fluid class
     *
     * \param[in] iforce contribution to the fluid residual
     *
     * \param[in] ivel unused in this implementation     *
     *
     */
    void apply_interface_values(Teuchos::RCP<Epetra_Vector> iforce,
        Teuchos::RCP<Epetra_Vector> ivel = Teuchos::null) override;

    /**
     * \brief Resets the external forces acting on the fluid to zero
     */
    virtual void reset_external_forces();

   private:
    /// fluid base algorithm object
    Teuchos::RCP<Adapter::Fluid> fluidadapter_;
  };

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
