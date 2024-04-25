/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for xfem-ale-fluids with moving boundaries

\level 2

*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_ADAPTER_FLD_FLUID_ALE_XFEM_HPP
#define FOUR_C_ADAPTER_FLD_FLUID_ALE_XFEM_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_fluid_ale.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ADAPTER
{
  // forward declarations
  class Coupling;

  /// fluid with moving interfaces implemented by the XFEM
  class FluidAleXFEM : public FluidAle
  {
   public:
    /// constructor
    explicit FluidAleXFEM(const Teuchos::ParameterList& prbdyn, std::string condname);

    /*========================================================================*/
    //! @name Misc
    /*========================================================================*/

    /// return the boundary discretization that matches the structure discretization
    Teuchos::RCP<DRT::Discretization> BoundaryDiscretization();

    /// communication object at the struct interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor> const& StructInterface();

    //@}

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
    //@}
  };

}  // namespace ADAPTER

FOUR_C_NAMESPACE_CLOSE

#endif
