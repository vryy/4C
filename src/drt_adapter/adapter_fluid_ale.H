/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_ale.H

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_FLUID_ALE_H
#define ADAPTER_FLUID_ALE_H

#include <Epetra_Comm.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "adapter_fluid_moving_boundary.H"
#include "adapter_fluid_base_algorithm.H"
#include "adapter_ale.H"
#include "adapter_coupling.H"

namespace ADAPTER
{

  /// fluid on ale
  class FluidAle : public FluidMovingBoundary
  {
  public:

    FluidAle(const Teuchos::ParameterList& prbdyn,
                    std::string condname);

    /// fluid field
    Fluid& FluidField() { return fluid_.FluidField(); }
    const Fluid& FluidField() const { return fluid_.FluidField(); }

    /// ale field
    Ale& AleField() { return ale_.AleField(); }
    const Ale& AleField() const { return ale_.AleField(); }

    virtual Teuchos::RCP<DRT::Discretization> Discretization();
    virtual const LINALG::MapExtractor& Interface() const;

    virtual void PrepareTimeStep();
    virtual void Update();
    virtual void Output();
    virtual double ReadRestart(int step);

    virtual void NonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,Teuchos::RCP<Epetra_Vector> ivel);

    virtual void RobinNonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<Epetra_Vector> iforce);

    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                        double dt);

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceFluidVelocity();

    virtual int Itemax() const { return FluidField().Itemax(); }
    virtual void SetItemax(int itemax) { FluidField().SetItemax(itemax); }

    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape();

    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  protected:

    //! @name Transfer helpers

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const;

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> AleToFluidField(Teuchos::RCP<const Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<const Epetra_Vector> iv) const;

    //@}

    //Coupling& FluidAleFieldCoupling()             { return coupfa_; }
    //const Coupling& FluidAleFieldCoupling() const { return coupfa_; }

  private:

    FluidBaseAlgorithm fluid_;
    AleBaseAlgorithm ale_;

    /// coupling of fluid and ale (whole field)
    Coupling               coupfa_;

    /// coupling of fluid and ale (interface only)
    Coupling               icoupfa_;

    /// coupling of fluid and ale at the free surface
    Coupling               fscoupfa_;

  };

}

#endif
#endif
