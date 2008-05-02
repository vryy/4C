/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_xfem.H

\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_FLUID_XFEM_H
#define ADAPTER_FLUID_XFEM_H

#include <Epetra_Comm.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "adapter_fluid_moving_boundary.H"
#include "adapter_fluid_base_algorithm.H"
#include "adapter_coupling.H"

namespace ADAPTER
{

  /// fluid with moving interfaces implemented by the XFEM
  class FluidXFEM : public FluidMovingBoundary
  {
  public:

    FluidXFEM(const Teuchos::ParameterList& prbdyn,
                    std::string condname,
                    RCP<DRT::Discretization>);

    /// fluid field
    Fluid& FluidField() { return fluid_.FluidField(); }
    const Fluid& FluidField() const { return fluid_.FluidField(); }

    virtual Teuchos::RCP<DRT::Discretization> Discretization();
    virtual const LINALG::MapExtractor& Interface() const;

    virtual void PrepareTimeStep();
    virtual void Update();
    virtual void Output();
    virtual double ReadRestart(int step);

    virtual void NonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,Teuchos::RCP<Epetra_Vector> ivel);
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                        double dt);

    virtual void RobinNonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,
						   Teuchos::RCP<Epetra_Vector> ivel,
						   Teuchos::RCP<Epetra_Vector> iforce);

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceFluidVelocity();

    virtual int Itemax() const { return FluidField().Itemax(); }
    virtual void SetItemax(int itemax) { FluidField().SetItemax(itemax); }

    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape();

    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  protected:

    //! @name Transfer helpers

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> SolidToFluid(Teuchos::RCP<Epetra_Vector> iv) const;

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> SolidToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> FluidToSolid(Teuchos::RCP<Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> FluidToSolid(Teuchos::RCP<const Epetra_Vector> iv) const;

    //@}

    //Coupling& FluidAleFieldCoupling()             { return coupfa_; }
    //const Coupling& FluidAleFieldCoupling() const { return coupfa_; }

  private:

    FluidBaseAlgorithm fluid_;

    /// the interface map setup for interface <-> full translation
    LINALG::MapExtractor                   interface_;

    /// coupling of interface in to parallel distributions
    Coupling               icoupsf_;

    Teuchos::RCP<DRT::Discretization>      boundarydis_solidparalleldistrib_;
    //Teuchos::RCP<DRT::Discretization>      boundarydis_fluidparalleldistrib_;

  };

}

#endif
#endif
