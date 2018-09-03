/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_ale_xfem.H

\brief Fluid field adapter for xfem-ale-fluids with moving boundaries

<pre>
Maintainer: Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 - 15249
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ADAPTER_FLUID_ALE_XFEM_H
#define ADAPTER_FLUID_ALE_XFEM_H

#include "ad_fld_fluid_ale.H"

namespace ADAPTER
{
  /// forward declarations

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
    virtual void NonlinearSolve(
        Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel);

    /// relaxation solve
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(
        Teuchos::RCP<Epetra_Vector> idisp, double dt);
    //@}

    /*========================================================================*/
    //! @name Extract interface forces
    /*========================================================================*/

    /// After the fluid solve we need the forces at the FSI interface.
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();
    //@}

    /*========================================================================*/
    //! @name extract helpers
    /*========================================================================*/

    /// extract the interface velocity at time t^(n+1)
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp();

    /// extract the interface velocity at time t^n
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();
    //@}
    //@}
  };

}  // namespace ADAPTER

#endif
