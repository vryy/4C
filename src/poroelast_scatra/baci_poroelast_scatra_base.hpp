/*----------------------------------------------------------------------*/
/*! \file

 \brief  base class for all poroelasticity scalar transport interaction algorithms

\level 2

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_POROELAST_SCATRA_BASE_HPP
#define FOUR_C_POROELAST_SCATRA_BASE_HPP

#include "baci_config.hpp"

#include "baci_adapter_algorithmbase.hpp"
#include "baci_coupling_adapter_volmortar.hpp"
#include "baci_poroelast_base.hpp"
#include "baci_scatra_algorithm.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/
namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
  class MortarVolCoupl;
  class FluidPoro;
  class FPSIStructureWrapper;
}  // namespace ADAPTER

// namespace POROELAST
//{
//  class PoroBase;
//}

/*----------------------------------------------------------------------*
 |                                                                       |
 *----------------------------------------------------------------------*/
namespace POROELASTSCATRA
{
  /// base class of algorithms for scalar transport in porous media
  class PoroScatraBase : public ADAPTER::AlgorithmBase
  {
   public:
    /// create using a Epetra_Comm
    explicit PoroScatraBase(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams);  // Problem builder

    //! Main time loop.
    virtual void Timeloop() = 0;

    //! prepare time step for single fields
    virtual void PrepareTimeStep(bool printheader = true)
    {
      dserror("not implemented in base class. override in subclass.");
    };

    //! perform iteration loop between fields
    virtual void Solve() { dserror("not implemented in base class. override in subclass."); };

    //! prepare output
    virtual void PrepareOutput()
    {
      dserror("not implemented in base class. override in subclass.");
    };

    //! update time step
    void Update() override { dserror("not implemented in base class. override in subclass."); };

    //! write output print to screen
    void Output() override { dserror("not implemented in base class. override in subclass."); };

    //! read and set fields needed for restart
    void ReadRestart(int restart) override = 0;

    //! setup for single fields
    virtual void SetupSystem();

    //! Build the combined dirichlet map of the monolithic poro problem
    virtual void BuildCombinedDBCMap() { poro_->BuildCombinedDBCMap(); };

    //! perform result test
    void TestResults(const Epetra_Comm& comm);

    //! apply solution of poro-problem to scatra
    void SetPoroSolution();

    //! apply solution of scatra to poro
    void SetScatraSolution();

    //! return pointer to porous medium problem
    const Teuchos::RCP<POROELAST::PoroBase>& PoroField() { return poro_; };

    //! return pointer to interstitial fluid
    const Teuchos::RCP<ADAPTER::FluidPoro>& FluidField() { return poro_->FluidField(); };

    //! return pointer to porous structure
    const Teuchos::RCP<ADAPTER::FPSIStructureWrapper>& StructureField()
    {
      return poro_->StructureField();
    };

    //! return pointer to scalar transport problem
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ScaTraField() { return scatra_->ScaTraField(); };

    //! return pointer to scalar problem adapter base class
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> ScaTraFieldBase() { return scatra_; };

    //! setup solver (for monolithic only)
    virtual bool SetupSolver() { return true; };

   protected:
    //! setup up of dofsets for two way coupling
    void ReplaceDofSets(Teuchos::RCP<DRT::Discretization> structdis,
        Teuchos::RCP<DRT::Discretization> fluiddis, Teuchos::RCP<DRT::Discretization> scatradis);

    //! setup up coupling objects if necessary
    void SetupCoupling(Teuchos::RCP<DRT::Discretization> structdis,
        Teuchos::RCP<DRT::Discretization> fluiddis, Teuchos::RCP<DRT::Discretization> scatradis);

    //! Pointer to the porous media problem. (poroelastic)
    Teuchos::RCP<POROELAST::PoroBase> poro_;
    //! Pointer to the ScaTra problem.     (scatra)
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;

    //! @name Volume Mortar stuff

    //! flag for matchinggrid
    const bool matchinggrid_;

    //! volume coupling (using mortar) adapter
    Teuchos::RCP<CORE::ADAPTER::MortarVolCoupl> volcoupl_structurescatra_;
    Teuchos::RCP<CORE::ADAPTER::MortarVolCoupl> volcoupl_fluidscatra_;
    //@}

   private:
    //! apply displacement fields to scatra
    void SetMeshDisp();
    //! apply velocity fields to scatra
    void SetVelocityFields();
  };
}  // namespace POROELASTSCATRA


BACI_NAMESPACE_CLOSE

#endif
