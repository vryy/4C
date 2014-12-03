/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_ale.H

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ADAPTER_FLUID_ALE_H
#define ADAPTER_FLUID_ALE_H

#include <Epetra_Comm.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "ad_fld_moving_boundary.H"
#include "ad_fld_base_algorithm.H"

namespace ADAPTER
{
  // forward declarations

  class Coupling;
  class AleFluidWrapper;


  /// fluid on ale
  class FluidAle : public FluidMovingBoundary
  {
  public:

    FluidAle(const Teuchos::ParameterList& prbdyn,
                    std::string condname);

    /// fluid field
    const Teuchos::RCP<ADAPTER::Fluid>& FluidField() { return fluid_; }

    /// ale field
    const Teuchos::RCP<ADAPTER::AleFluidWrapper>& AleField() const { return ale_; }

    virtual Teuchos::RCP<DRT::Discretization> Discretization();
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& Interface() const { return fluid_->Interface(); }

    virtual void PrepareTimeStep();
    virtual void Update();
    virtual void Output();
    virtual double ReadRestart(int step);

    virtual void NonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,Teuchos::RCP<Epetra_Vector> ivel);

    virtual void ApplyInterfaceValues(Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel);

    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                        double dt);

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();

    virtual int Itemax() const { return fluid_->Itemax(); }
    virtual void SetItemax(int itemax) { fluid_->SetItemax(itemax); }

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

    /// coupling of fluid and ale (whole field)
    Teuchos::RCP<Coupling>               coupfa_;

    /// coupling of fluid and ale (interface only)
    Teuchos::RCP<Coupling>               icoupfa_;

    /// coupling of fluid and ale at the free surface
    Teuchos::RCP<Coupling>               fscoupfa_;

    /// coupling of fluid and ale for the ale update condition
    Teuchos::RCP<Coupling>               aucoupfa_;

  private:

    /// problem-specific Fluid-wrapper
    Teuchos::RCP<ADAPTER::Fluid> fluid_;

    /// problem-specific ALE-wrapper
    Teuchos::RCP<ADAPTER::AleFluidWrapper>  ale_;

  };

}

#endif
