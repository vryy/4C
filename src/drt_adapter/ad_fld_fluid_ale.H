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
#include "../drt_ale/ale.H"

namespace ADAPTER
{
  // forward declarations

  class Coupling;


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
    ALE::Ale& AleField() { return ale_.AleField(); }
    const ALE::Ale& AleField() const { return ale_.AleField(); }

    virtual Teuchos::RCP<DRT::Discretization> Discretization();
    virtual const FLD::UTILS::MapExtractor& Interface() const;

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

  private:

    FluidBaseAlgorithm fluid_;
    ALE::AleBaseAlgorithm ale_;

    /// coupling of fluid and ale (whole field)
    Teuchos::RCP<Coupling>               coupfa_;

    /// coupling of fluid and ale (interface only)
    Teuchos::RCP<Coupling>               icoupfa_;

    /// coupling of fluid and ale at the free surface
    Teuchos::RCP<Coupling>               fscoupfa_;

    /// coupling of fluid and ale for the ale update condition
    Teuchos::RCP<Coupling>               aucoupfa_;

  };

}

#endif
