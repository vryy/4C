/*----------------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_ale.H

\brief Solver for fluid field on a moving ALE mesh

\level 1

\maintainer Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
*/
/*----------------------------------------------------------------------------*/

#ifndef ADAPTER_FLUID_ALE_H
#define ADAPTER_FLUID_ALE_H

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace ADAPTER
{
  class CouplingBase;
  class Coupling;
  class AleFluidWrapper;
}  // namespace ADAPTER

/*----------------------------------------------------------------------------*/
/* headers */
#include <Epetra_Comm.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "ad_fld_moving_boundary.H"
#include "ad_fld_base_algorithm.H"

//#include "../drt_fsi/fsi_dirichletneumann_volcoupl.H"

namespace FSI
{
  class InterfaceCorrector;
}

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace ADAPTER
{
  /// fluid on ale
  class FluidAle : public FluidMovingBoundary
  {
   public:
    FluidAle(const Teuchos::ParameterList& prbdyn, std::string condname);

    /// fluid field
    const Teuchos::RCP<ADAPTER::Fluid>& FluidField() { return fluid_; }

    /// ale field
    const Teuchos::RCP<ADAPTER::AleFluidWrapper>& AleField() const { return ale_; }

    /// discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization();

    /// fluid interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor> const& Interface() const
    {
      return fluid_->Interface();
    }

    /// Prepare a single time step
    virtual void PrepareTimeStep();

    /// Update to go from time step \f$t_n\f$ to \f$t_{n+1}\f$
    virtual void Update();

    /// Output current state of simulation
    virtual void Output();

    /// Read resatart data
    virtual double ReadRestart(int step  ///< step number to restart from
    );

    virtual void NonlinearSolve(
        Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel);

    virtual void NonlinearSolveVolCoupl(Teuchos::RCP<Epetra_Vector> idisp,
        Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<FSI::InterfaceCorrector> icorrector);

    virtual void ApplyInterfaceValues(
        Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel);

    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(
        Teuchos::RCP<Epetra_Vector> idisp, double dt);

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();

    virtual int Itemax() const { return fluid_->Itemax(); }
    virtual void SetItemax(int itemax) { fluid_->SetItemax(itemax); }

    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape();

    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

   protected:
    //! @name Transfer helpers
    //@{

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> AleToFluidField(
        Teuchos::RCP<Epetra_Vector> iv  ///< ALE vector (to be converted)
        ) const;

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> AleToFluidField(
        Teuchos::RCP<const Epetra_Vector> iv  ///< ALE vector (to be converted)
        ) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> FluidToAle(
        Teuchos::RCP<Epetra_Vector> iv  ///< Fluid vector (to be converted)
        ) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> FluidToAle(
        Teuchos::RCP<const Epetra_Vector> iv  ///< Fluid vector (to be converted)
        ) const;

    //@}

    /// coupling of fluid and ale (whole field)
    Teuchos::RCP<CouplingBase> coupfa_;

    /// coupling of fluid and ale (interface or volume...)
    Teuchos::RCP<CouplingBase> icoupfa_;

    /// coupling of fluid and ale at the free surface
    Teuchos::RCP<Coupling> fscoupfa_;

    /// coupling of fluid and ale for the ale update condition
    Teuchos::RCP<Coupling> aucoupfa_;

   private:
    /// problem-specific Fluid-wrapper
    Teuchos::RCP<ADAPTER::Fluid> fluid_;

    /// problem-specific ALE-wrapper
    Teuchos::RCP<ADAPTER::AleFluidWrapper> ale_;

    /// problem specific time parameter list
    const Teuchos::ParameterList& timeparams_;
  };

}  // namespace ADAPTER

#endif
