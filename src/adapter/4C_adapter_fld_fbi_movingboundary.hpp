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

namespace CORE::LINALG
{
  class SparseOperator;
}

namespace ADAPTER
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
    const Teuchos::RCP<ADAPTER::Fluid>& FluidField() override { return fluidadapter_; }

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
    int Itemax() const override;

    /// set the maximum number of iterations for the fluid field
    void SetItemax(int itemax) override;

    //@}

    /*========================================================================*/
    //! @name others
    /*========================================================================*/

    /// integrate the interface shape functions
    Teuchos::RCP<Epetra_Vector> integrate_interface_shape() override;

    /// create the testing of fields
    Teuchos::RCP<CORE::UTILS::ResultTest> CreateFieldTest() override;

    /// Get velocity at timestep n+1
    virtual Teuchos::RCP<const Epetra_Vector> Velnp();

    virtual Teuchos::RCP<const FLD::Meshtying> GetMeshtying();

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
        Teuchos::RCP<const CORE::LINALG::SparseOperator> matrix);

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
    virtual void ResetExternalForces();

   private:
    /// fluid base algorithm object
    Teuchos::RCP<ADAPTER::Fluid> fluidadapter_;
  };

}  // namespace ADAPTER

FOUR_C_NAMESPACE_CLOSE

#endif
